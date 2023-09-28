#include "RelicMuonMatching.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "skheadC.h"
#include "ParticleCand.h"
#include <inttypes.h>
#include <algorithm>

RelicMuonMatching::RelicMuonMatching():Tool(){}

extern "C" void tdiff_muon_(int* nevhwsk_tar, int* it0xsk_tar, int* mode, double* timediff);

bool RelicMuonMatching::Initialise(std::string configfile, DataModel &data){
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("match_window", match_window);  // [s]
	match_window *= 1E9; // convert to [ns]
	match_window_ticks = match_window * COUNT_PER_NSEC;
	
	std::string relicWriterName, muWriterName;
	m_variables.Get("rfmReaderName", rfmReaderName);
	if(m_data->Trees.count(rfmReaderName)==0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+rfmReaderName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// input reader
	rfmReaderLUN = m_data->GetLUN(rfmReaderName);
	rfmReader = m_data->Trees.at(rfmReaderName);
	
	// get LUNs needed to passing common block data from reco algorithms to TTrees
	m_variables.Get("muWriterName", muWriterName);
	m_variables.Get("relicWriterName", relicWriterName);
	muWriterLUN = m_data->GetLUN(muWriterName);
	if(muWriterLUN<0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+muWriterName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	relicWriterLUN = m_data->GetLUN(relicWriterName);
	if(relicWriterLUN<0){
		Log(m_unique_name+" Error! Failed to find "+relicWriterName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// Get output TTrees and add new branches to store matches
	TreeManager* muMgr = skroot_get_mgr(&muWriterLUN);
	TTree* muTree = muMgr->GetOTree();
	muTree->Branch("MatchedEvNums", &MatchedEvNums);
	muTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	muTree->Branch("MatchedParticleE", &MatchedParticleE);
	
	TreeManager* relicMgr = skroot_get_mgr(&relicWriterLUN);
	TTree* relicTree = relicMgr->GetOTree();
	relicTree->Branch("MatchedEvNums", &MatchedEvNums);
	relicTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	relicTree->Branch("MatchedParticleE", &MatchedParticleE);
	
	// see if recording this as a cut, and if so make it
	get_ok = m_variables.Get("muSelectorName",muSelectorName);
	if(get_ok){
		m_data->AddCut(muSelectorName, "relic_mu_tdiff", "time to closest relic",-60,60);
		m_data->AddCut(muSelectorName, m_unique_name, "require at least one relic candidate within +-60s",1,1000);
	}
	// repeat for relics (n.b. just records the number of matched muons, no actual cut placed)
	get_ok = m_variables.Get("relicSelectorName",relicSelectorName);
	if(get_ok){
		m_data->AddCut(relicSelectorName, "relic_mu_tdiff", "time to closest relic",-60,60);
		m_data->AddCut(relicSelectorName, m_unique_name, "record number of muons within +-60s",1,1000);
	}
	
	// thought this might suppress ranlux printouts, but it seems not...?
	int zero=0;
	ran_verbosity_(&zero);
	
	
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
	
	if(eventType==EventType::AFT){
		// we'll want to save the AFT associated with any written out particles
		// but the AFT event should follow the corresponding prompt event in the output TTree,
		// and since we don't write out muons/relics until we've found all their matches,
		// we will also need to postpone saving the AFT
		// So for now just note that the previous event had an AFT
		// But, note the preceding event may not have passed cuts!
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
	
	// FIXME tdiff_muon masks the middle 17 bits. What in the damn hell is going on. Is that a bug?
	
	// for long-scale times we have a 47-bit clock that runs at 1.92 ticks per ns (#defined as COUNT_PER_NSEC)
	// the upper 32 bits of this clock are in skheadqb_.nevhwsk
	// the lower 15 bits are in the lower 15 bits of it0sk.
	// this resets when a new run is manually started (automaton run changes do not reset it).
	// (we currently handle rollover manually, but don't do anything to handle matching
	//  across a manual run change that would have reset the 47-bit clock)
	//int64_t currentTime = (skheadqb_.nevhwsk << 15) + (skheadqb_.it0sk & int32_t(std::pow(2,15)));
	
	// i have no idea why, but let's match what happens in tdiff_muon:
	// mask the lower 17 bits of counter_32, then shift 15 bits up, and add to it0sk
	int64_t currentTime = ((skheadqb_.nevhwsk & 0x1FFFF) << 15) + skheadqb_.it0sk;
	
	if(eventType==EventType::LowE){
		// match this relic candidate to any held muon candidates
		RelicMuonMatch(true, currentTime, 0, 0);
	} else {
		// since we searched for muons using a manual subtrigger trigger scan, we may have
		// more than one muon in this event, each with different times.
		// to do things properly, we need to loop over all these muons and account for thir t0_sub offsets.
		
		// n.b. our matching window is +-60s, compared to a readout window of only ~2ms at the most
		// (T2K, AFT triggers). So it might be sufficient to just match against the primary trigger time...
		
		std::vector<int> muonTimes;
		m_data->CStore.Get("muonTimes", muonTimes);
		//std::cout<<"this muon event had "<<muonTimes.size()<<" muon times"<<std::endl;
		for(int i=0; i<muonTimes.size(); ++i){
			// muonTimes are 'swtrgt0ctr' value, which is t0_sub from get_sub_triggers.
			// this is a counts offset from it0sk, so:
			currentTime += muonTimes.at(i);
			//std::cout<<"\tmuon time: "<<currentTime<<std::endl;
			RelicMuonMatch(false, currentTime, i, muonTimes.at(i));
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
	
	// write finished candidates to file
	if(m_data->writeOutRelics.size()){
		WriteRelicInfo();
	}
	if(m_data->muonsToRec.size()){
		Log(m_unique_name+" We've got "+toString(m_data->muonsToRec.size()) +" muons to write!",v_debug,m_verbose);
		WriteMuonInfo();
	}
	
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
	
	// write candidates to file
	if(m_data->writeOutRelics.size()) WriteRelicInfo();
	if(m_data->muonsToRec.size()) WriteMuonInfo();
	
	
	if(!muSelectorName.empty() && !relicSelectorName.empty()){
		std::cout<<"compared "<<tdiffcount<<" muon-relic pairs and found "<<passing_tdiffcount<<" that were within 60s of each other"<<std::endl;
		std::cout<<"c.f. "<<m_data->Selectors.at(muSelectorName)->GetEntries("relic_mu_tdiff")<<" recorded muons and "
		         <<m_data->Selectors.at(relicSelectorName)->GetEntries("relic_mu_tdiff")<<" recorded relics"<<std::endl;
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

bool RelicMuonMatching::RelicMuonMatch(bool loweEventFlag, int64_t currentTime, int subtrg_num, int it0xsk){
	
	// make a new ParticleCand to encapsulate the minimal info about this muon/relic candidate.
	ParticleCand currentParticle;
	currentParticle.EventNumber = skhead_.nevsk;
	currentParticle.SubTriggerNumber = subtrg_num;
	currentParticle.EventTime = currentTime;
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
	bool firstmatch=true;
	for(int i = 0; i < targetDeque->size(); i++){
		ParticleCand& targetCand = targetDeque->at(i);
		
		//calculate time difference between this event and the target
		int64_t timeDiff = (currentTime - targetCand.EventTime);
		
		// to account for 48-bit clock rollover we can compare event numbers
		if(currentParticle.EventNumber > targetCand.EventNumber){
			// tdiff ought to be positive; correct for rollover if not
			if(timeDiff<0) timeDiff += int64_t(std::pow(2,47));
		} else {
			if(timeDiff>0) timeDiff -= int64_t(std::pow(2,47));
		}
		
		// compare to tdiff_muon result
		double tdiffmu;
		int pre_or_post = 0; //(targetCand.EventNumber > currentParticle.EventNumber);
		tdiff_muon_(&targetCand.nevhwsk, &targetCand.it0xsk, &pre_or_post, &tdiffmu);
		//std::cout<<"\ttdiff_muon: "<<tdiffmu<<", timeDiff: "<<(double(timeDiff/COUNT_PER_NSEC)/1.E9)<<std::endl;
		
		// ah sod it
		timeDiff = tdiffmu*COUNT_PER_NSEC;
		
		// make a note of the time diff. The selector is just a recorder, so this won't
		// affect any actual selections, but we can use it to get the distribution of time diffs
		++tdiffcount;
		if(!relicSelectorName.empty() && !muSelectorName.empty()){
			m_data->ApplyCut((loweEventFlag ? relicSelectorName : muSelectorName), "relic_mu_tdiff", timeDiff);
		}
		
		//If the time difference between the two events is less than 60 seconds then "match" the particles.
		//N.B. since events are time ordered, timediff is always positive
		if(timeDiff < match_window_ticks){
			++passing_tdiffcount;
			
			/*
			std::cout<<"matching "<<((loweEventFlag) ? "relic" : "muon") << " event "
			         <<currentParticle.InEntryNumber<<" at "<<currentTime<<" to target event "
			         <<targetCand.InEntryNumber<<" at "<<targetCand.EventTime<<"; tdiff "
			         <<timeDiff<<"\n"
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
			currentParticle.matchedParticleEvNum.push_back(targetCand.OutEntryNumber);
			currentParticle.matchedParticleTimeDiff.push_back(timeDiff * -COUNT_PER_NSEC);
			currentParticle.matchedParticleBSEnergy.push_back(targetCand.LowECommon.bsenergy);
			
			targetCand.matchedParticleEvNum.push_back(currentParticle.OutEntryNumber);
			targetCand.matchedParticleTimeDiff.push_back(timeDiff*COUNT_PER_NSEC);
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
	// FIXME re-enable, and only do for muons?
	//if(!loweEventFlag && currentDeque->size() > 150){
		for(int i = 0; i < currentDeque->size() - 2; i++){
			ParticleCand& targetCand = currentDeque->at(i);
			
			//int64_t timeDiff = (currentTime - targetCand.EventTime);
			//if(timeDiff<0) timeDiff += int64_t(std::pow(2,47));  // fix rollover
			
			double tdiffmu;
			int pre_or_post = 0;
			tdiff_muon_(&targetCand.nevhwsk, &targetCand.it0xsk, &pre_or_post, &tdiffmu);
			int64_t timeDiff = tdiffmu*COUNT_PER_NSEC;
			
			// if this particle is more than 60s after our current one, we'll have found all its matches
			// so we can prune it now.
			if(timeDiff > match_window_ticks){
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
	//}
	
	return true;
}


bool RelicMuonMatching::WriteRelicInfo(){
	
	int originalEntry = rfmReader->GetEntryNumber();
	
	std::vector<ParticleCand>& writeOutRelics = m_data->writeOutRelics;
	
	for(int writeEvent = 0; writeEvent < writeOutRelics.size(); writeEvent++){
		
		Log(m_unique_name+" Writing out next relic; entry "+
		    toString(writeOutRelics[writeEvent].OutEntryNumber),v_debug,m_verbose);
		
		// reload the TTree entry for this lowe event
		// (TODO not ideal as we're not reading linearly; can we buffer the required info?)
		int currentEntry = rfmReader->GetEntryNumber();
		if(writeOutRelics[writeEvent].InEntryNumber != currentEntry){
			m_data->getTreeEntry(rfmReaderName, writeOutRelics[writeEvent].InEntryNumber);
		}
		
		// for relics, keep the entire readout window. might be useful for spallation?
		//delete_outside_hits_();
		
		//pass header, tqreal and tqareal commons to output TTree branch variables
		skroot_set_tree_(&relicWriterLUN);
		
		// if we ran lfallfit before caching the relic, also put that into the output TTree
		skroot_lowe_ = writeOutRelics[writeEvent].LowECommon;
		skroot_set_lowe_(&relicWriterLUN,
		                 skroot_lowe_.bsvertex,
		                 skroot_lowe_.bsresult,
		                 skroot_lowe_.bsdir,
		                 skroot_lowe_.bsgood,
		                 &skroot_lowe_.bsdirks,
		                 skroot_lowe_.bseffhit,
		                 &skroot_lowe_.bsenergy,
		                 &skroot_lowe_.bsn50,
		                 &skroot_lowe_.bscossun,
		                 skroot_lowe_.clvertex,
		                 skroot_lowe_.clresult,
		                 skroot_lowe_.cldir,
		                 &skroot_lowe_.clgoodness,
		                 &skroot_lowe_.cldirks,
		                 skroot_lowe_.cleffhit,
		                 &skroot_lowe_.clenergy,
		                 &skroot_lowe_.cln50,
		                 &skroot_lowe_.clcossun,
		                 &skroot_lowe_.latmnum,
		                 &skroot_lowe_.latmh,
		                 &skroot_lowe_.lmx24,
		                 &skroot_lowe_.ltimediff,
		                 &skroot_lowe_.lnsratio,
		                 skroot_lowe_.lsdir,
		                 &skroot_lowe_.spaevnum,
		                 &skroot_lowe_.spaloglike,
		                 &skroot_lowe_.sparesq,
		                 &skroot_lowe_.spadt,
		                 &skroot_lowe_.spadll,
		                 &skroot_lowe_.spadlt,
		                 &skroot_lowe_.spamuyn,
		                 &skroot_lowe_.spamugdn,
		                 skroot_lowe_.posmc,
		                 skroot_lowe_.dirmc,
		                 skroot_lowe_.pabsmc,
		                 skroot_lowe_.energymc,
		                 &skroot_lowe_.darkmc,
		                 &skroot_lowe_.islekeep,
		                 &skroot_lowe_.bspatlik,
		                 &skroot_lowe_.clpatlik,
		                 &skroot_lowe_.lwatert,
		                 &skroot_lowe_.lninfo,
		                 skroot_lowe_.linfo);
		
		// update branch variables w/ info about matches
		MatchedEvNums = writeOutRelics[writeEvent].matchedParticleEvNum;
		MatchedTimeDiff = writeOutRelics[writeEvent].matchedParticleTimeDiff;
		MatchedParticleE = writeOutRelics[writeEvent].matchedParticleBSEnergy;
		
		//invoke TTree::Fill
		skroot_fill_tree_(&relicWriterLUN);
		
		// if the relic had an AFT trigger, write that out now
		if(writeOutRelics[writeEvent].hasAFT){
			Log(m_unique_name+" Appending AFT to relic tree",v_debug,m_verbose);
			m_data->getTreeEntry(rfmReaderName, writeOutRelics[writeEvent].InEntryNumber+1);
			skroot_set_tree_(&relicWriterLUN);
			skroot_fill_tree_(&relicWriterLUN);
		}
		
	}
	
	// reload last treeReader entry
	if(originalEntry != rfmReader->GetEntryNumber()){
		m_data->getTreeEntry(rfmReaderName, originalEntry);
	}
	
	writeOutRelics.clear();
	
	return true;
	
}

bool RelicMuonMatching::WriteMuonInfo(){
	
	std::vector<ParticleCand>& muonsToRec = m_data->muonsToRec;
	
	int currentEntry = rfmReader->GetEntryNumber();
	
	for(int i = 0; i < muonsToRec.size(); i++){
		
		Log(m_unique_name+" Writing out next muon; entry "+
		    toString(muonsToRec[i].OutEntryNumber),v_debug,m_verbose);
		
		// (re-)load the muon entry to save
		if(muonsToRec[i].InEntryNumber != currentEntry){
			// if we're going to call mufit we should disable bad channel masking.
			// Since we don't do muon fitting here any more, this may be redundant...
			// ...but will bad channel hits still be copied to the output file if they're masked?
			int current_badch_masking = combad_.imaskbadopt;
			int mufitbadopt = 0;
			skbadopt_(&mufitbadopt); // literally all this does is update the common block value.
			m_data->getTreeEntry(rfmReaderName, muonsToRec[i].InEntryNumber);
			// reset for the next read
			skbadopt_(&current_badch_masking);
			
			// if this was an untagged muon we should also shift the time window accordingly
			if(muonsToRec[i].SubTriggerNumber!=0){
				set_timing_gate_(&muonsToRec[i].it0xsk);
				int neglun = -std::abs(rfmReaderLUN);
				skcread_(&neglun, &get_ok);
				// get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
				if(get_ok!=0){
					Log(m_unique_name+" Error! skcread returned "+toString(get_ok)
					    +" when reloading muon subtrigger!",v_error,m_verbose);
					return false;
				}
			}
		}
		
		/* muon reconstruction now moved to later stage */
		
		// for muons, only keep hits around 1.3us trigger
		//delete_outside_hits_();
		// XXX FIXME XXX will this remove hits from any further muon subtriggers in this readout?
		
		// set header and tq info (epsecially updated hits)
		skroot_set_tree_(&muWriterLUN);
		
		// update branch variables w/ info about matches
		MatchedEvNums = muonsToRec[i].matchedParticleEvNum;
		MatchedTimeDiff = muonsToRec[i].matchedParticleTimeDiff;
		MatchedParticleE = muonsToRec[i].matchedParticleBSEnergy;
		
		// invoke TTree::Fill
		skroot_fill_tree_(&muWriterLUN);
		
		// sanity check
		TreeManager* mgr = skroot_get_mgr(&muWriterLUN);
		TTree* mutree = mgr->GetOTree();
		//std::cout<<"Muon tree now has "<<mutree->GetEntries()<<" entries"<<std::endl;
		
		// if the muon had an AFT trigger, write that out now
		if(muonsToRec[i].hasAFT){
			Log(m_unique_name+" Appending AFT to muon tree",v_debug,m_verbose);
			m_data->getTreeEntry(rfmReaderName, muonsToRec[i].InEntryNumber+1);
			skroot_set_tree_(&muWriterLUN);
			skroot_fill_tree_(&muWriterLUN);
			//std::cout<<"Added AFT: Muon tree now has "<<mutree->GetEntries()<<" entries"<<std::endl;
		} else {
			//std::cout<<"Muon "<<muonsToRec[i].EventNumber<<" has no AFT?!"<<std::endl;
		}
	}
	
	//return the previous entry so that no issues are caused with other tools
	if(currentEntry != rfmReader->GetEntryNumber()){
		m_data->getTreeEntry(rfmReaderName, currentEntry);
	}
	
	m_data->muonsToRec.clear();
	
	return true;
}


// DON'T LOOK AT THE CODE BELOW HERE SHHHHHHHHHHHHHHH

/*
float RelicMuonMatching::rollOver(unsigned long long int currentTime, unsigned long long int targetTime){
	unsigned long long int bitOne = 1;
	unsigned long long int tDiff;
	tDiff = currentTime - targetTime;
	tDiff = (bitOne << 47) + tDiff;
	
}

bool RelicMuonMatching::AddParticletoDeque(std::deque<ParticleCand>& addToThisDeque){
	unsigned long long int newTime = bitshiftTime(skheadqb_.it0xsk, skheadqb_.nevhwsk);
	ParticleCand newParticle;
	newParticle.EventNumber = skhead_.nevsk;
	newParticle.EventTime = newTime;
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
