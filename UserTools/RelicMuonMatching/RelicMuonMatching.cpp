#include "RelicMuonMatching.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "skheadC.h"
#include "ParticleCand.h"
#include <inttypes.h>
#include <algorithm>

RelicMuonMatching::RelicMuonMatching():Tool(){}


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
		m_data->AddCut(muSelectorName, m_unique_name, "require at least one relic candidate within +-60s");
	}
	// repeat for relics (n.b. just records the number of matched muons, no actual cut placed)
	get_ok = m_variables.Get("relicSelectorName",relicSelectorName);
	if(get_ok){
		m_data->AddCut(relicSelectorName, m_unique_name, "record number of muons within +-60s");
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
	
	// if the toolchain made it here, the current rfm file entry is either a muon or relic candidate!
	get_ok = m_data->vars.Get("newMuon", muonFlag);  // is it a muon?
	
	// FIXME tdiff_muon masks the middle 17 bits. What in the damn hell is going on. Is that a bug?
	
	// for long-scale times we have a 47-bit clock that runs at 1.92 ticks per ns (#defined as COUNT_PER_NSEC)
	// the upper 32 bits of this clock are in skheadqb_.nevhwsk
	// the lower 15 bits are in the lower 15 bits of it0sk.
	// this resets when a new run is manually started (automaton run changes do not reset it).
	// (we currently handle rollover manually, but don't do anything to handle matching
	//  across a manual run change that would have reset the 47-bit clock)
	int64_t currentTime = (skheadqb_.nevhwsk << 15) + (skheadqb_.it0sk & int32_t(std::pow(2,15)));
	
	if(muonFlag){
		// since we searched for muons using a manual subtrigger trigger scan, we may have
		// more than one muon in this event, each with different times.
		// to do things properly, we need to loop over all these muons and account for thir t0_sub offsets.
		
		// n.b. our matching window is +-60s, compared to a readout window of only ~2ms at the most
		// (T2K, AFT triggers). So it might be sufficient to just match against the primary trigger time...
		
		std::vector<int> muonTimes;
		m_data->CStore.Get("muonTimes", muonTimes);
		for(int i=0; i<muonTimes.size(); ++i){
			// muonTimes are 'swtrgt0ctr' value, which is t0_sub from get_sub_triggers.
			// this is a counts offset from it0sk, so:
			currentTime += muonTimes.at(i);
			RelicMuonMatch(muonFlag, currentTime, i, muonTimes.at(i));
		}
		
	} else {
		RelicMuonMatch(muonFlag, currentTime, 0, 0);
	}
	
	// prune any match candidates that have dropped off our window of interest
	if(muonsToRemove.size() > 0){
		RemoveFromDeque(muonsToRemove, m_data->muonCandDeque);
	}
	if(relicsToRemove.size() > 0){
		RemoveFromDeque(relicsToRemove, m_data->relicCandDeque);
	}
	
	// write finished candidates to file
	if(m_data->writeOutRelics.size()){
		WriteRelicInfo();
	}
	if(m_data->muonsToRec.size()){
		WriteMuonInfo();
	}
	
	return true;
}


bool RelicMuonMatching::Finalise(){
	
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

bool RelicMuonMatching::RelicMuonMatch(bool muonFlag, int64_t currentTime, int subtrg_num, int it0xsk){
	
	// make a new ParticleCand to encapsulate the minimal info about this muon/relic candidate.
	ParticleCand currentParticle;
	currentParticle.EventNumber = skhead_.nevsk;
	currentParticle.SubTriggerNumber = subtrg_num;
	currentParticle.EventTime = currentTime;
	currentParticle.it0xsk = it0xsk;
	currentParticle.InEntryNumber = rfmReader->GetEntryNumber();
	currentParticle.LowECommon = skroot_lowe_;
	
	// get the deque of in-memory targets to match this new event against
	// if this event is a muon then the targets are relic candidates, and vice versa
	std::deque<ParticleCand>* currentDeque = nullptr;
	std::deque<ParticleCand>* targetDeque = nullptr;
	if(muonFlag){
		currentParticle.PID = 2;
		currentDeque = &m_data->muonCandDeque;
		targetDeque = &m_data->relicCandDeque;
		// we'll assign its output tree entry number if/when it gets matched to a relic
	} else {
		currentParticle.PID = 1;
		currentDeque = &m_data->relicCandDeque;
		targetDeque = &m_data->muonCandDeque;
		// we save every relic, so can already assign its output ttree entry number
		currentParticle.OutEntryNumber = nextrelicentry;
		++nextrelicentry;
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
		//If the time difference between the two events is less than 60 seconds then "match" the particles.
		//N.B. since events are time ordered, timediff is always positive
		if(timeDiff < match_window_ticks ){
			
			std::cout<<"matching "<<((muonFlag) ? "muon" : "relic") << " event "
			         <<currentParticle.InEntryNumber<<" at "<<currentTime<<" to target event "
			         <<targetCand.InEntryNumber<<" at "<<targetCand.EventTime<<"; tdiff "
			         <<timeDiff<<"\n"
			         <<"\tnext muon entry: "<<nextmuentry<<"\n"
			         <<"\tnext relic entry: "<<nextrelicentry<<std::endl;
			
			// if this is the first match of this particle, set its event number in the output file
			// and increment the counter for the next event which will be written out
			if(firstmatch){
				std::cout<<"first match for current "<<((muonFlag) ? "muon" : "relic") <<std::endl;
				if(muonFlag){
					currentParticle.OutEntryNumber = nextmuentry;
					++nextmuentry;
				}
				firstmatch=false;
			}
			if(targetCand.matchedParticleEvNum.size()==0){
				std::cout<<"first match for target "<<((muonFlag) ? "relic" : "muon") <<std::endl;
				if(muonFlag){
					targetCand.OutEntryNumber = nextrelicentry;
					++nextrelicentry;
				} else {
					targetCand.OutEntryNumber = nextmuentry;
					++nextmuentry;
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
			std::cout<<((muonFlag) ? "muon" : "relic")<<" event "<<currentParticle.InEntryNumber
			         <<" is >60s after target event "<<targetCand.InEntryNumber<<std::endl;
			//any subsequent events will also be >60s after this target event;
			//which is to say we'll find no more matches for this target.
			if(muonFlag){
				// add it to the set of relic candidates ready to write out
				m_data->writeOutRelics.push_back(targetCand);
				// remove it from the set of relic candidates being matched
				relicsToRemove.push_back(targetCand.EventNumber);
				// make a note of this relic and its number of matches
				if(!relicSelectorName.empty()){
					m_data->ApplyCut(relicSelectorName, m_unique_name,
					                 targetCand.matchedParticleEvNum.size());
				}
			} else {
				//we'll find a lot of muons, but we're only interested in ones matched to relic candidates.
				// only add it to the set of muons to record if it was matched to at least one relic.
				if(targetCand.matchedParticleEvNum.size()){
					m_data->muonsToRec.push_back(targetCand);
				}
				// remove it from the set of muons being matched
				muonsToRemove.push_back(targetCand.EventNumber);
				// make a note of this muon and its number of matches
				if(!muSelectorName.empty()){
					m_data->ApplyCut(muSelectorName, m_unique_name,
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
	//if(currentDeque->size() > 150){
		for(int i = 0; i < currentDeque->size() - 1; i++){
			int64_t timeDiff = (currentTime - currentDeque->at(i).EventTime);
			if(currentParticle.EventNumber > currentDeque->at(i).EventNumber){
				if(timeDiff<0) timeDiff += int64_t(std::pow(2,47));
			} else {
				if(timeDiff>0) timeDiff -= int64_t(std::pow(2,47));
			}
			if(timeDiff > match_window_ticks && ! currentDeque->at(i).matchedParticleEvNum.size()){
				std::vector<int>* removedeq= (muonFlag) ? &muonsToRemove : &relicsToRemove;
				removedeq->push_back(currentDeque->at(i).EventNumber);
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
		
		std::cout<<"writing out next relic; entry "<<writeOutRelics[writeEvent].OutEntryNumber<<std::endl;
		
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
		// FIXME i believe we do not do this now, leaving lowe reconstruction for later
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
	}
	
	// reload last treeReader entry
	// FIXME maybe we don't need to do this if this is the end of the ToolChain?
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
		
		std::cout<<"writing out next muon; entry "<<muonsToRec[i].OutEntryNumber<<std::endl;
		
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
		delete_outside_hits_();
		// FIXME um.... will this remove hits from any further muon subtriggers in this readout?
		
		// set header and tq info (epsecially updated hits)
		skroot_set_tree_(&muWriterLUN);
		
		// update branch variables w/ info about matches
		MatchedEvNums = muonsToRec[i].matchedParticleEvNum;
		MatchedTimeDiff = muonsToRec[i].matchedParticleTimeDiff;
		MatchedParticleE = muonsToRec[i].matchedParticleBSEnergy;
		
		// invoke TTree::Fill
		skroot_fill_tree_(&muWriterLUN);
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
