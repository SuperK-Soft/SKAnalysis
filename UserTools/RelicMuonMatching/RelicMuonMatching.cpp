#include "RelicMuonMatching.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "skheadC.h"
#include "ParticleCand.h"
#include <inttypes.h>
#include <algorithm>
#include <iomanip>
#include <bitset>

RelicMuonMatching::RelicMuonMatching():Tool(){}

extern "C" void event_tdiff_(int32_t* nevhwsk_cur, int32_t* it0xsk_cur, int32_t* nevhwsk_tar, int32_t* it0xsk_tar, int32_t* mode, double* timediff, int64_t* ticksdiff);
extern "C" void tdiff_muon_(int*, int*, int*, double*);

bool RelicMuonMatching::Initialise(std::string configfile, DataModel &data){
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("match_window", match_window);  // [s]
	match_window *= 1E9; // convert to [ns]
	match_window_ticks = match_window * COUNT_PER_NSEC;
	
	std::string rfmReaderName;
	m_variables.Get("rfmReaderName", rfmReaderName);
	if(m_data->Trees.count(rfmReaderName)==0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+rfmReaderName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// input reader
	rfmReader = m_data->Trees.at(rfmReaderName);
	
	// see if recording this as a cut, and if so make it
	get_ok = m_variables.Get("muSelectorName",muSelectorName);
	if(get_ok){
		m_data->AddCut(muSelectorName, "relic_mu_tdiff", "time to closest relic",true,"MU",-60,60);
		m_data->AddCut(muSelectorName, m_unique_name, "require at least one relic candidate within +-60s",true,1,1000);
	}
	// repeat for relics (n.b. just records the number of matched muons, no actual cut placed)
	get_ok = m_variables.Get("relicSelectorName",relicSelectorName);
	if(get_ok){
		m_data->AddCut(relicSelectorName, "relic_mu_tdiff", "time to closest relic",true,"MU",-60,60);
		m_data->AddCut(relicSelectorName, m_unique_name, "record number of muons within +-60s",true,1,1000);
	}
	
	// thought this might suppress ranlux printouts, but it seems not...?
	int zero=0;
	ran_verbosity_(&zero);
	// relic_sk4_ana has a modified version '$RELIC_WORK_DIR/lomufit/mufit/src/ranlux_nowrite.F'
	
	m_variables.Get("distrosFile",distros_file);
	if(!distros_file.empty()){
		hb.MakeFile(distros_file);
		hb.SaveHists(false);
	}
	
	return true;
}


bool RelicMuonMatching::Execute(){
	
	if(skhead_.nsubsk != currentSubRun){
		currentSubRun = skhead_.nsubsk;
		std::cout << "Subrun number:            " << skhead_.nsubsk << std::endl;
	}
	
	// if the toolchain made it here, the current rfm file entry is either a muon, a relic candidate,
	// or an AFT following one of the above.
	EventType lastEventType = eventType;
	m_data->vars.Get("eventType", eventType);
	
	/*
	logmessage<<m_unique_name<<" Matching event "<<skhead_.nevsk<<" as "<<eventType;
	if(skhead_.nevsk==(last_nevsk+1)){
		logmessage<<", lastEventType was "<<lastEventType<<std::endl;
	} else {
		logmessage<<", last event was unrelated"<<std::endl;
	}
	Log(logmessage.str(),v_debug,m_verbose);
	last_nevsk = skhead_.nevsk;
	*/
	
	// we want to save the AFT associated with any written out particles.
	// The AFT event should follow the corresponding prompt event in the output TTree,
	// and since we don't write out muons/relics until we've found all their matches,
	// we will need to postpone saving the AFT until that happens.
	// So for now just note when a previous event being checked had an AFT
	// N.B. since the preceding event may not have passed cuts, not every AFT
	// will need to be saved
	if(eventType==EventType::AFT){
		std::deque<ParticleCand>* thedeque=nullptr;
		if(lastEventType==EventType::Muon){
			Log(m_unique_name+" found AFT after muon",v_debug,m_verbose);
			thedeque = &m_data->muonCandDeque;
		} else {
			Log(m_unique_name+" found AFT after relic",v_debug,m_verbose);
			thedeque = &m_data->relicCandDeque;
			// the next entry in the output relic tree will be an AFT, so advance our relic entry counter
			++nextrelicentry;
			Log(m_unique_name+" Advancing relic entry to account for AFT after relic",v_debug,m_verbose);
		}
		if(thedeque->size() && thedeque->back().EventNumber==(skhead_.nevsk-1)){
			thedeque->back().hasAFT = true;
			Log(m_unique_name+" Setting AFT flag for "+(lastEventType==EventType::Muon ? "Muon " : "relic ")
			         +toString(thedeque->back().EventNumber),v_debug,m_verbose);
			if(lastEventType==EventType::Muon && thedeque->back().matchedParticleEvNum.size()>0){
				Log(m_unique_name+" Advancing muon entry to account for AFT after muon",v_debug,m_verbose);
				++nextmuentry;
			}
		}
		// and that's all we need to do for now
		return true;
	}
	
	if(eventType!=EventType::Muon && eventType!=EventType::LowE){
		Log(m_unique_name+" Error!!! Event is neither Mu nor LowE! Should not be here!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	if(eventType==EventType::LowE) ++reliccount;
	if(eventType==EventType::Muon) ++muoncount;
	
	// *************************** //
	// *** START SANITY CHECKS *** //
	// *************************** //
	
	// sanity check; calculate time since last readout
	//int64_t thiseventticks = ((skheadqb_.nevhwsk & ~0x1FFFF) << 15) + (int64_t(skheadqb_.it0sk) & 0xFFFFFFFF); // doesn't work?
	int64_t thiseventticks = (skheadqb_.nevhwsk & ~0x1FFFF);
	thiseventticks = thiseventticks << 15;
	int64_t iticks = *reinterpret_cast<uint32_t*>(&skheadqb_.it0sk);
	//std::cout<<"upticks:   "<<std::bitset<64>(thiseventticks)<<"\nlowticks:  "<<std::bitset<64>(iticks & 0xFFFFFFFF)<<std::endl;
	thiseventticks += iticks & 0xFFFFFFFF;
	//std::cout<<"thisticks: "<<std::bitset<64>(thiseventticks)<<std::endl;
	//std::cout<<"lastticks: "<<std::bitset<64>(lasteventticks)<<std::endl;
	int64_t ticksDiff = (thiseventticks - lasteventticks);
	//std::cout<<"diffticks: "<<std::bitset<64>(ticksDiff)<<std::endl;
	
	/*
	// insanity check -- doesn't trigger, we're good.
	double tdiff_ns;
	int pre=0;
	int64_t tdiff_ticks;
	event_tdiff_(&skheadqb_.nevhwsk, &skheadqb_.it0sk, &lastnevhwsk, &lastit0sk, &pre, &tdiff_ns, &tdiff_ticks);
	double tdiff_ns2;
	tdiff_muon_(&lastnevhwsk, &lastit0sk, &pre, &tdiff_ns2);
	if(tdiff_ns!=tdiff_ns2){
		std::cerr<<"*** SOMETHN GOT BROKE ***"<<std::endl;
	}
	*/
	
	//if(skheadqb_.nevhwsk < lastnevhwsk){
	//	std::cout<<"nevhwsk rollover at event "<<skhead_.nevsk<<std::endl;
	//}
	//if(skheadqb_.it0sk < lastit0sk){
	//	std::cout<<"it0sk rollover at event "<<skhead_.nevsk<<std::endl;
	//}
	
	if(ticksDiff<0){
		++num_rollovers;
		std::cerr<<"!!! ROLLOVER !!! at "<<skhead_.nevsk<<"\n\tthis evt nevhwsk: "<<skheadqb_.nevhwsk<<", it0sk: "<<skheadqb_.it0sk
			     <<", total ticks: "<<thiseventticks<<"\n\tlast evt nevhwsk: "<<lastnevhwsk<<", it0sk: "<<lastit0sk
			     <<", total ticks: "<<lasteventticks<<"\n\tdiff: "<<ticksDiff;
			ticksDiff += (int64_t(1) << 47);
		std::cerr<<", corrected ticks diff: "<<ticksDiff
			     <<"\n\tevts since last rollover: "<<(skhead_.nevsk-last_rollover_nevsk)<<std::endl;
		last_rollover_nevsk = skhead_.nevsk;
	}
	/*
	// as commented below, this sometimes triggers but the ticks difference is the same
	if((tdiff_ns-double(ticksDiff/COUNT_PER_NSEC)>10)){
		std::bitset<64> bitsdiff(ticksDiff);
		std::bitset<64> bitsdiff2(tdiff_ticks);
		std::cerr<<"***TDIFF ERROR: "<<bitsdiff<<" vs ref: "<<bitsdiff2<<std::endl;
	}
	*/
	double secs_since_last = double(ticksDiff/COUNT_PER_NSEC)/1.E9;
	//std::cout<<m_unique_name<<" secs to last event: "<<secs_since_last<<std::endl;
	
	if(!distros_file.empty()){
		if(eventType==EventType::Muon){
			ticksDiff = (thiseventticks - lastmuticks);
			if(ticksDiff<0) ticksDiff += (int64_t(1) << 47);
			secs_since_last = double(ticksDiff/COUNT_PER_NSEC)/1.E9;
			hb.Fill("mu_to_mu_secs", secs_since_last);
			lastmuticks = thiseventticks;
			
			ticksDiff = (thiseventticks - lastrelicticks);
			if(ticksDiff<0) ticksDiff += (int64_t(1) << 47);
			secs_since_last = double(ticksDiff/COUNT_PER_NSEC)/1.E9;
			hb.Fill("mu_to_relic_secs", secs_since_last);
			
		} else if(eventType==EventType::LowE){
			ticksDiff = (thiseventticks - lastrelicticks);
			if(ticksDiff<0) ticksDiff += (int64_t(1) << 47);
			secs_since_last = double(ticksDiff/COUNT_PER_NSEC)/1.E9;
			hb.Fill("relic_to_relic_secs", secs_since_last);
			lastrelicticks = thiseventticks;
			
			ticksDiff = (thiseventticks - lastmuticks);
			if(ticksDiff<0) ticksDiff += (int64_t(1) << 47);
			secs_since_last = double(ticksDiff/COUNT_PER_NSEC)/1.E9;
			hb.Fill("relic_to_mu_secs", secs_since_last);
			
		}
	}
	lastnevhwsk = skheadqb_.nevhwsk;
	lastit0sk = skheadqb_.it0sk;
	lasteventticks = thiseventticks;
	
	// ************************* //
	// *** END SANITY CHECKS *** //
	// ************************* //
	
	
	// for long-scale times we have a 47-bit clock that runs at 1.92 ticks per ns (#defined as COUNT_PER_NSEC)
	// the upper 32 bits of this clock are in skheadqb_.nevhwsk
	// the lower 15 bits are in the lower 15 bits of it0sk.
	// XXX for some reason in tdiff_muon the lower *17* bits of nevhwsk are masked, not just the lower 15 bits???
	// the 47-bit clock resets when a new run is manually started (automaton run changes do not reset it).
	// (we currently handle rollover manually, but don't do anything to handle matching
	//  across a manual run change that would have reset the 47-bit clock)
	int64_t currentTicks = (skheadqb_.nevhwsk & ~0x1FFFF);
	currentTicks = currentTicks << 15;
	/* equivalent to tdiff_muon's:
	currentTicks = skheadqb_.nevhwsk >> 17;
	currentTicks = currentTicks << 32;
	*/
	
	if(eventType==EventType::LowE){
		
		// combine with it0sk. The catch here is that it0sk is a 32-bit *signed* integer.
		// in fortran there does not appear to be unsigned integer types..!
		// In the 2's complement representation the upper bits of negative numbers
		// will be filled with 1's not 0's, so converting it to a 64-bit number
		// can result in the insertion of 1's in the upper 32 bits, which then screw up the result.
		// in c++ we can just interpret it as uint so it fills with 0's!
		currentTicks += *reinterpret_cast<uint32_t*>(&skheadqb_.it0sk);
		
		// match this relic candidate to any held muon candidates
		RelicMuonMatch(true, currentTicks, 0, 0);
		
	} else {
		// since we searched for muons using a manual subtrigger trigger scan, we may have
		// more than one muon in this event, each with different times.
		// to do things properly, we need to loop over all these muons and account for thir t0_sub offsets.
		
		// n.b. our matching window is +-60s, compared to a readout window of only ~2ms at the most
		// (T2K, AFT triggers). So it would probably be sufficient to just match against the primary trigger time...
		// but to save each muon as an independent 1.3us event we need properly shifted timestamps
		
		std::vector<int> muonTimes;
		m_data->CStore.Get("muonTimes", muonTimes);
		//std::cout<<"this muon event had "<<muonTimes.size()<<" muon times"<<std::endl;
		for(int i=0; i<muonTimes.size(); ++i){
			
			// muonTimes are 'swtrgt0ctr' value, which is t0_sub from get_sub_triggers.
			// this is a counts offset from it0sk
			int32_t it0xsk = skheadqb_.it0sk + muonTimes.at(i);
			currentTicks += *reinterpret_cast<uint32_t*>(&it0xsk);
			
			// match this muon candidate to any held relic candidates
			RelicMuonMatch(false, currentTicks, i, it0xsk);
		}
		
	}
	
	// prune any match candidates that have dropped off our window of interest
	Log(m_unique_name+" Relics to prune: "+toString(relicsToRemove.size())+
	                  ", muons to prune: "+toString(muonsToRemove.size()),v_debug,m_verbose);
	if(muonsToRemove.size() > 0){
		RemoveFromDeque(muonsToRemove, m_data->muonCandDeque);
	}
	if(relicsToRemove.size() > 0){
		RemoveFromDeque(relicsToRemove, m_data->relicCandDeque);
	}
	
	Log(m_unique_name+" Relics to Write out: "+toString(m_data->writeOutRelics.size())+
	                  ", muons to write out: "+toString(m_data->muonsToRec.size()),v_debug,m_verbose);
	
	return true;
}


bool RelicMuonMatching::Finalise(){
	
	// write out any remaining relics still being matched
	for(int i = 0; i < m_data->relicCandDeque.size(); i++){
		ParticleCand& targetCand = m_data->relicCandDeque.at(i);
		m_data->writeOutRelics.push_back(targetCand);
		if(!relicSelectorName.empty()){
			m_data->ApplyCut(relicSelectorName, m_unique_name,
			                 targetCand.matchedParticleEvNum.size());
		}
	}
	m_data->relicCandDeque.clear();
	
	// write out any remaining muons with a match
	for(int i = 0; i < m_data->muonCandDeque.size(); i++){
		ParticleCand& targetCand = m_data->muonCandDeque.at(i);
		if(targetCand.matchedParticleEvNum.size()) m_data->muonsToRec.push_back(targetCand);
		if(!muSelectorName.empty()){
			m_data->ApplyCut(muSelectorName, m_unique_name,
			                 targetCand.matchedParticleEvNum.size());
		}
	}
	m_data->muonCandDeque.clear();
	
	std::cout<<"checked "<<muoncount<<" muons and "<<reliccount<<" relics"<<std::endl;
	std::cout<<"compared "<<tdiffcount<<" muon-relic pairs and found "<<passing_tdiffcount<<" that were within 60s of each other"<<std::endl;
	
	/*
	// sanity check
	if(!muSelectorName.empty() && !relicSelectorName.empty()){
		// XXX ! currently returns 0 for both?
		std::cout<<"c.f. "<<m_data->Selectors.at(muSelectorName)->GetEntries("relic_mu_tdiff")
		                  <<" recorded muons and "
		                  <<m_data->Selectors.at(relicSelectorName)->GetEntries("relic_mu_tdiff")
		                  <<" recorded relics"<<std::endl;
	}
	*/
	
	if(!distros_file.empty()){
		hb.Save();
	}
	
	return true;
}

bool RelicMuonMatching::RemoveFromDeque(std::vector<int>& particlesToRemove, std::deque<ParticleCand>& particleDeque){
	for(int j=0; j<particlesToRemove.size(); ++j){
		for(int i=0; i<particleDeque.size(); ++i){
			if(particleDeque[i].EventNumber == particlesToRemove[j]){
				particleDeque.erase(particleDeque.begin() + i);
				break;
			}
		}
	}
	particlesToRemove.clear();
	return true;
}

bool RelicMuonMatching::RelicMuonMatch(bool loweEventFlag, int64_t currentTicks, int subtrg_num, int32_t it0xsk){
	
	// make a new ParticleCand to encapsulate the minimal info about this muon/relic candidate.
	ParticleCand currentParticle;
	currentParticle.EventNumber = skhead_.nevsk;
	currentParticle.SubTriggerNumber = subtrg_num;
	currentParticle.EventTicks = currentTicks;
	currentParticle.NumRollovers = num_rollovers;
	currentParticle.nevhwsk = skheadqb_.nevhwsk;
	currentParticle.it0xsk = it0xsk;
	currentParticle.InEntryNumber = rfmReader->GetEntryNumber();
	currentParticle.LowECommon = skroot_lowe_;
	currentParticle.hasAFT = false;
	
	// get the deque of in-memory targets to match this new event against
	// if this event is a muon then the targets are relic candidates, and vice versa
	std::deque<ParticleCand>* currentDeque = nullptr;
	std::deque<ParticleCand>* targetDeque = nullptr;
	if(loweEventFlag){
		currentParticle.PID = 1;
		currentDeque = &m_data->relicCandDeque;
		targetDeque = &m_data->muonCandDeque;
		// we save every relic, so can already assign its output ttree entry number
		currentParticle.OutEntryNumber = nextrelicentry;
		++nextrelicentry;
	} else {
		currentParticle.PID = 2;
		currentDeque = &m_data->muonCandDeque;
		targetDeque = &m_data->relicCandDeque;
		// we'll assign its output tree entry number if/when it gets matched to a relic
	}
	
	// scan over targets, oldest to newest
	
	Log(m_unique_name+" matching this "+(loweEventFlag ? "lowE" : "muon")+" candidate to "
	    +toString(targetDeque->size())+" targets",v_debug,m_verbose);
	
	bool firstmatch=true;
	for(int i = 0; i < targetDeque->size(); i++){
		ParticleCand& targetCand = targetDeque->at(i);
		
		//calculate time difference in ticks between this event and the target
		int64_t ticksDiff = (currentTicks - targetCand.EventTicks);
		
		// tdiff ought to be positive; if not counter must have rolled over
		if(ticksDiff<0) ticksDiff +=  (int64_t(1) << 47);
		
		if(subtrg_num==0 && i==0){
			Log(m_unique_name+" secs to oldest candidate "+toString(i)+": "
			   +toString(double(ticksDiff/COUNT_PER_NSEC)/1E9),v_error,m_verbose);
		}
		
		// validate: compare to tdiff_muon result (returns ns)
		/*
		double tdiffmu;
		int64_t tdiffmuticks;
		int pre_or_post = 0; //(targetCand.EventNumber > currentParticle.EventNumber);
		event_tdiff_(&currentParticle.nevhwsk, &currentParticle.it0xsk, &targetCand.nevhwsk, &targetCand.it0xsk, &pre_or_post, &tdiffmu, &tdiffmuticks);
		// for some reason there's actually a sizeable difference in result with fortran,
		// but the source is simply the conversion from idiff (ticks diff as integer)
		// to timediff (nanoseconds diff as dfloat), which doesn't give the same result as c++....
		if((tdiffmu - double(ticksDiff/COUNT_PER_NSEC)) > 100){
			std::cout<<"\ttdiff_muon: "<<tdiffmu<<", timeDiff: "<<double(ticksDiff/COUNT_PER_NSEC)<<std::endl;
		}
		*/
		
		// make a note of the time diff. The selector is just a recorder, so this won't
		// affect any actual selections, but we can use it to get the distribution of time diffs
		++tdiffcount;
		
		// get subtrigger number from target muon if the candidate is a relic
		if(loweEventFlag){
			subtrg_num = targetCand.SubTriggerNumber;
		}
		// fill relic matching tdiff histogram
		if(!relicSelectorName.empty()){
			m_data->ApplyCut(relicSelectorName, "relic_mu_tdiff", (ticksDiff/COUNT_PER_NSEC)/1.E9, subtrg_num);
		}
		// fill muon matching tdiff histogram
		if(!muSelectorName.empty()){
			m_data->ApplyCut(muSelectorName, "relic_mu_tdiff", (ticksDiff/COUNT_PER_NSEC)/1.E9, subtrg_num);
		}
		
		//If the time difference between the two events is less than 60 seconds then "match" the particles.
		//N.B. since events are time ordered, timediff is always positive
		if(ticksDiff < match_window_ticks){
			++passing_tdiffcount;
			
			/*
			std::cout<<"matching "<<((loweEventFlag) ? "relic" : "muon") << " event "
			         <<currentParticle.InEntryNumber<<" at "<<currentTicks<<" to target event "
			         <<targetCand.InEntryNumber<<" at "<<targetCand.EventTicks<<"; tdiff "
			         <<(ticksDiff/COUNT_PER_NSEC)<<"ns\n"
			         <<"\tnext muon entry: "<<nextmuentry<<"\n"
			         <<"\tnext relic entry: "<<nextrelicentry<<std::endl;
			*/
			
			// if this is the first match of this particle, set its event number in the output file
			// and increment the counter for the next event which will be written out
			if(firstmatch){
				Log(m_unique_name+" First match for current "+((loweEventFlag) ? "relic" : "muon"),v_debug,m_verbose);
				if(!loweEventFlag){
					currentParticle.OutEntryNumber = nextmuentry;
					++nextmuentry;
				}
				firstmatch=false;
			}
			if(targetCand.matchedParticleEvNum.size()==0){
				Log(m_unique_name+" First match for target "+((loweEventFlag) ? "muon" : "relic"),v_debug,m_verbose);
				// if this is the first match for a muon, we now know we'll be writing it out
				// so can set its output entry number and increment that for the next.
				if(loweEventFlag){
					targetCand.OutEntryNumber = nextmuentry;
					++nextmuentry;
					if(targetCand.hasAFT) ++nextmuentry;
				}
			}
			
			//add the event # of the current event to the target particle's "matched particle" list and add the
			//event # of the target particle to the current particle's "matched particle" list
			currentParticle.matchedParticleEvNum.push_back(targetCand.EventNumber);
			currentParticle.matchedParticleEntryNum.push_back(targetCand.OutEntryNumber);
			currentParticle.matchedParticleHasAFT.push_back(currentParticle.hasAFT);
			currentParticle.matchedParticleTimeDiff.push_back(ticksDiff / -COUNT_PER_NSEC);
			currentParticle.matchedParticleBSEnergy.push_back(targetCand.LowECommon.bsenergy);
			
			targetCand.matchedParticleEvNum.push_back(currentParticle.EventNumber);
			targetCand.matchedParticleEntryNum.push_back(currentParticle.OutEntryNumber);
			targetCand.matchedParticleHasAFT.push_back(currentParticle.hasAFT);
			targetCand.matchedParticleTimeDiff.push_back(ticksDiff / COUNT_PER_NSEC);
			targetCand.matchedParticleBSEnergy.push_back(currentParticle.LowECommon.bsenergy);
			
		//otherwise the current event came more than 60 seconds after the target event.
		} else {
			Log(m_unique_name+((loweEventFlag) ? "relic" : "muon")+" entry "
			    +toString(currentParticle.InEntryNumber)+" is >60s after target entry "
			    +toString(targetCand.InEntryNumber),v_debug,m_verbose);
			//any subsequent events will also be >60s after this target event;
			//which is to say we'll find no more matches for this target.
			if(loweEventFlag){
				Log(m_unique_name+" Muon "+toString(targetCand.InEntryNumber)+" matched to "
				    +toString(targetCand.matchedParticleEvNum.size())+" relics",v_debug,m_verbose);
				// we'll find a lot of muons, but we're only interested in ones matched to relic candidates.
				// only add it to the set of muons to record if it was matched to at least one relic.
				if(targetCand.matchedParticleEvNum.size()){
					m_data->muonsToRec.push_back(targetCand);
				}
				// remove it from the set of muons being matched
				muonsToRemove.push_back(targetCand.EventNumber);
				if(!muSelectorName.empty()){
					// make a note of this muon and its number of matches
					m_data->ApplyCut(muSelectorName, m_unique_name,
					                 targetCand.matchedParticleEvNum.size());
				}
			} else {
				Log(m_unique_name+" Relic "+toString(targetCand.InEntryNumber)+" matched to "
				    +toString(targetCand.matchedParticleEvNum.size())+" muons",v_debug,m_verbose);
				// add it to the set of relic candidates ready to write out
				m_data->writeOutRelics.push_back(targetCand);
				// remove it from the set of relic candidates being matched
				relicsToRemove.push_back(targetCand.EventNumber);
				if(!relicSelectorName.empty()){
					// make a note of this relic and its number of matches
					m_data->ApplyCut(relicSelectorName, m_unique_name,
					                 targetCand.matchedParticleEvNum.size());
				}
			}
		}
	}
	
	currentDeque->push_back(currentParticle);
	
	//There are ~2.5 cosmic ray muons interating in SK per second,
	//whereas relic candidates passing upstream cuts may be quite rare.
	//If we only prune muons when processing a relic,
	//we could end up accumulating an unreasonably large stack of muons.
	//We can safely prune any muons more than 60s older than the current event that have no matches.
	//only bother with this scan if we have >150 muons (~60s) of muons
	if(!loweEventFlag && currentDeque->size() > 150){
		Log(m_unique_name+" We have "+toString(currentDeque->size())
		    +" muons, dropping any older than the current one",v_debug,m_verbose);
		for(int i = 0; i < (int(currentDeque->size()) - 2); i++){
			ParticleCand& targetCand = currentDeque->at(i);
			
			/*
			std::cout<<"currentTicks:" <<currentTicks<<", target ticks: "<<targetCand.EventTicks<<std::endl;
			std::bitset<64> cticksbits(currentTicks);
			std::bitset<64> tticksbits(targetCand.EventTicks);
			std::cout<<"as bits: "<<cticksbits<<std::endl<<tticksbits<<std::endl;
			*/
			
			int64_t ticksDiff = (currentTicks - targetCand.EventTicks);
			
			/*
			std::cout<<"diff: "<<ticksDiff<<std::endl;
			std::bitset<64> diffbits(ticksDiff);
			std::cout<<"bitsdiff: "<<diffbits<<std::endl;
			*/
			
			if(ticksDiff<0) ticksDiff += int64_t(std::pow(2,47));  // fix rollover
			
			/*
			diffbits = std::bitset<64>(ticksDiff);
			double timdiff = double(ticksDiff)/double(COUNT_PER_NSEC);
			printf("to ns: %.7f\n",timdiff);
			
			double tdiffmu;
			int64_t tdiffmuticks;
			int pre_or_post = 0;
			event_tdiff_(&currentParticle.nevhwsk, &currentParticle.it0xsk, &targetCand.nevhwsk, &targetCand.it0xsk, &pre_or_post, &tdiffmu, &tdiffmuticks);
			int64_t ticksDiff1 = tdiffmu*COUNT_PER_NSEC;
			if(ticksDiff != ticksDiff1){
				std::cout<<"manual diff: "<<ticksDiff<<", event_tdiff: "<<ticksDiff1<<std::endl;
				exit(1);
			}
			*/
			
			// if this particle is more than 60s after our current one, we'll have found all its matches
			// so we can prune it now.
			if(ticksDiff > match_window_ticks){
				if(loweEventFlag){
					m_data->writeOutRelics.push_back(targetCand);
					relicsToRemove.push_back(targetCand.EventNumber);
					if(!relicSelectorName.empty()){
						m_data->ApplyCut(relicSelectorName, m_unique_name,
						                 targetCand.matchedParticleEvNum.size());
					}
				} else {
					if(targetCand.matchedParticleEvNum.size()) m_data->muonsToRec.push_back(targetCand);
					muonsToRemove.push_back(targetCand.EventNumber);
					if(!muSelectorName.empty()){
						m_data->ApplyCut(muSelectorName, m_unique_name,
						                 targetCand.matchedParticleEvNum.size());
					}
				}
			} else {
				break;
			}
		}
	}
	
	return true;
}

// DON'T LOOK AT THE CODE BELOW HERE SHHHHHHHHHHHHHHH

/*
float RelicMuonMatching::rollOver(unsigned long long int currentTicks, unsigned long long int targetTicks){
	unsigned long long int bitOne = 1;
	unsigned long long int tDiff;
	tDiff = currentTicks - targetTicks;
	tDiff = (bitOne << 47) + tDiff;
	
}

bool RelicMuonMatching::AddParticletoDeque(std::deque<ParticleCand>& addToThisDeque){
	unsigned long long int newTicks = bitshiftTime(skheadqb_.it0xsk, skheadqb_.nevhwsk);
	ParticleCand newParticle;
	newParticle.EventNumber = skhead_.nevsk;
	newParticle.EventTicks = newTicks;
	newParticle.InEntryNumber = rfmReader->GetEntryNumber();
	if(particleType == "LOWE"){
		newParticle.ReconEnergy = skroot_lowe_.bsenergy;
	} else{
		newParticle.ReconEnergy = 0.0;
	}
	addToThisDeque.push_back(newParticle);
	
	return true;
}

unsigned long long int RelicMuonMatching::bitshiftTime(unsigned long long int t0Time, unsigned long long int hardwareTime){
	
	unsigned long long int shiftedt0Time, shiftedhardwareTime, oneint;
	
	// equivalent to 00000000000000001111111111111111 in binary
	oneint = 65535;
	
	shiftedt0Time = t0Time >> 16;
	shiftedt0Time = shiftedt0Time << 16;
	
	shiftedt0Time = shiftedt0Time | (t0Time & oneint);
	
	shiftedhardwareTime = hardwareTime >> 17;
	shiftedhardwareTime = shiftedhardwareTime << 32;
	
	shiftedt0Time = shiftedt0Time + shiftedhardwareTime;
	
	t0Time = shiftedt0Time;
	
	
	return t0Time;
}
*/
