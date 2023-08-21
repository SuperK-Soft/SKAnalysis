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
	
	std::string relicWriterName, muWriterName;
	m_variables.Get("rfmReaderName", rfmReaderName);
	if(m_data->Trees.count(rfmReaderName)==0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+rfmReaderName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// input reader
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
	
	TreeManager* relicMgr = skroot_get_mgr(&relicWriterLUN);
	TTree* relicTree = relicMgr->GetOTree();
	relicTree->Branch("MatchedEvNums", &MatchedEvNums);
	relicTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	
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
	// see which it is:
	get_ok = m_data->vars.Get("newMuon", muonFlag);
	
	if(muonFlag){
		// since we searched for muons using a software trigger scan, we may have found
		// more than one muon event in this readout window. Each needs to be matched independently.
		std::vector<int> muonTimes;
		m_data->CStore.Get("muonTimes", muonTimes);
		
		for(int i=0; i<muonTimes.size(); ++i){
			// FIXME muonTimes are 'swtrgt0ctr' value. Presume this is same units as counter_32?
			float muonTime = muonTimes.at(i) * 32768./1.92;
			RelicMuonMatch(muonFlag, muonTime, i);
		}
		
	} else {
		rfmReader->Get("HEADER", myHeader);
		float currentTime = myHeader->counter_32 * 32768./1.92;
		RelicMuonMatch(muonFlag, currentTime, 0);
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
			if(particleDeque[i].EntryNumber == particlesToRemove[j]){
				particleDeque.erase(particleDeque.begin() + i);
				break;
			}
		}
	}
	particlesToRemove.clear();
	return true;
}

bool RelicMuonMatching::RelicMuonMatch(bool muonFlag, float currentTime, int subtrg_num){
	
	// make a new ParticleCand to encapsulate the minimal info about this muon/relic candidate.
	ParticleCand currentParticle;
	currentParticle.EventNumber = myHeader->nevsk;
	currentParticle.SubTriggerNumber = subtrg_num;
	currentParticle.EventTime = currentTime;
	currentParticle.EntryNumber = rfmReader->GetEntryNumber();
	currentParticle.LowECommon = skroot_lowe_;
	
	// get the deque of in-memory targets to match this new event against
	// if this event is a muon then the targets are relic candidates, and vice versa
	std::deque<ParticleCand>* currentDeque = nullptr;
	std::deque<ParticleCand>* targetDeque = nullptr;
	if(muonFlag){
		currentParticle.PID = 2;
		currentDeque = &m_data->muonCandDeque;
		targetDeque = &m_data->relicCandDeque;
	} else {
		currentParticle.PID = 1;
		currentDeque = &m_data->relicCandDeque;
		targetDeque = &m_data->muonCandDeque;
	}
	
	// scan over targets, oldest to newest
	for(int i = 0; i < targetDeque->size(); i++){
		ParticleCand& targetCand = targetDeque->at(i);
		//calculate time difference between this event and the target
		float timeDiff = (currentTime - targetCand.EventTime);
		//If the time difference between the two events is less than 60 seconds then "match" the particles.
		//N.B. since events are time ordered, timediff is always positive
		if(timeDiff < 60. * pow(10,9)){
			//add the event # of the current event to the target particle's "matched particle" list and add the
			//event # of the target particle to the current particle's "matched particle" list
			currentParticle.matchedParticleEvNum.push_back(targetCand.EventNumber);
			currentParticle.matchedParticleTimeDiff.push_back(timeDiff * -1.);
			currentParticle.matchedParticleBSEnergy.push_back(targetCand.LowECommon.bsenergy);
			
			targetCand.matchedParticleEvNum.push_back(currentParticle.EventNumber);
			targetCand.matchedParticleTimeDiff.push_back(timeDiff);
			targetCand.matchedParticleBSEnergy.push_back(currentParticle.LowECommon.bsenergy);
			
		//otherwise the current event came more than 60 seconds after the target event.
		} else {
			//any subsequent events will also be >60s after this target event;
			//which is to say we'll find no more matches for this target.
			if(muonFlag){
				// add it to the set of relic candidates ready to write out
				m_data->writeOutRelics.push_back(targetCand);
				// remove it from the set of relic candidates being matched
				relicsToRemove.push_back(targetCand.EntryNumber);
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
				muonsToRemove.push_back(targetCand.EntryNumber);
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
	if(muonFlag && currentDeque->size() > 150){
		for(int i = 0; i < currentDeque->size() - 1; i++){
			float timeDiff = (currentTime - currentDeque->at(i).EventTime);
			if(timeDiff > 60.e+9 && ! currentDeque->at(i).matchedParticleEvNum.size()){
				muonsToRemove.push_back(currentDeque->at(i).EntryNumber);
			} else {
				break;
			}
		}
	}
	
	return true;
}


bool RelicMuonMatching::WriteRelicInfo(){
	
	int originalEntry = rfmReader->GetEntryNumber();
	
	std::vector<ParticleCand>& writeOutRelics = m_data->writeOutRelics;
	
	for(int writeEvent = 0; writeEvent < writeOutRelics.size(); writeEvent++){
		
		// reload the TTree entry for this lowe event
		// (TODO not ideal as we're not reading linearly; can we buffer the required info?)
		int currentEntry = rfmReader->GetEntryNumber();
		if(writeOutRelics[writeEvent].EntryNumber != currentEntry){
			m_data->getTreeEntry(rfmReaderName, writeOutRelics[writeEvent].EntryNumber);
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
		
		if(muonsToRec[i].EntryNumber != currentEntry){
			m_data->getTreeEntry(rfmReaderName, muonsToRec[i].EntryNumber);
		}
		
		int muyn_org, muynf;
		float mbentry [4];
		float mm_entry [36];
		int n_left;
		
		//Muon reconstruction developed by Tomoeda and Yamaguchi
		mfmuselect_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &skroot_mu_.mugoodness, &muyn_org);
		
		//muyn == 1 - good fit
		//muyn == 0 - bad fit
		
		if(muyn_org > 0){
			skroot_mu_.muyn = 1;
		} else if(muyn_org < 0) {
			skroot_mu_.muyn = 0;
		} else {
			Log("Muyn_org returning as == 0. Not supported yet", v_error, m_verbose);
			return false;
		}
		
		//Apply fast fit if mfmuselect has returned a bad fit
		if(skroot_mu_.muyn == 0){
			mffastfast_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &muynf);
			skroot_mu_.mufast_flag = 1;
		}else{
			skroot_mu_.mufast_flag = 0;
		}
		
		skroot_mu_.muyn = muyn_org;
		if(skroot_mu_.muyn == 0){
			skroot_mu_.muyn = muynf;
		}
		
		int muboy_class, muboy_numtracks, muboy_numleft;
		float muboy_entry [4], muboy_dir [3], muboy_otherentry [36];
		float muboy_tracklength, muboy_goodness;
		
		//Apply muboy
		muboy_zbs_(&skhead_.nevsk,
		           &skroot_mu_.muboy_status,
		           &mbentry,
		           &skroot_mu_.muboy_dir,
		           &skroot_mu_.muboy_length,
		           &skroot_mu_.muboy_goodness,
		           &skroot_mu_.muboy_ntrack,
		           &mm_entry,
		           &n_left);
		
		/*
		muboy_zbs_(&skhead_.nevsk, &muboy_class, &muboy_entry,
		           &muboy_dir, &muboy_tracklength, 
		           &muboy_goodness, &muboy_numtracks,
		           &muboy_otherentry, &muboy_numleft);
		
		std::cout << "muboy result:\n"
		          << "\tclass: " << muboy_class << "\n"
		          << "\tdir: (" << muboy_dir[0] << ", "
		                        << muboy_dir[1] << ", "
		                        << muboy_dir[2] << ")\n"
		          << "\tentry: (" << muboy_entry[0] << ", "
		                        << muboy_entry[1] << ","
		                        << muboy_entry[2] << ")\n"
		          << "\ttrack length: " << muboy_tracklength << std::endl;
		*/
		
		// get muon track entry position(s)
		for(int track = 0; track < skroot_mu_.muboy_ntrack; track++){
			if(track == 0){
				skroot_mu_.muboy_entpos[0][track] = 3;
				skroot_mu_.muboy_entpos[1][track] = mbentry[1];
				skroot_mu_.muboy_entpos[2][track] = mbentry[2];
				skroot_mu_.muboy_entpos[3][track] = mbentry[3];
			}else{
				skroot_mu_.muboy_entpos[0][track] = mm_entry[4*track - 8];
				skroot_mu_.muboy_entpos[1][track] = mm_entry[4*track - 7];
				skroot_mu_.muboy_entpos[2][track] = mm_entry[4*track - 6];
				skroot_mu_.muboy_entpos[3][track] = mm_entry[4*track - 5];
			}
		}
		
		// recalculate muon energy with run-wise water transparency
		// first get water transparency for the current run
		float watert;
		if(skhead_.nrunsk != lastRun){
			int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];
			lfwater_(&days_to_run_start, &watert);
			lastRun = skhead_.nrunsk;
		}
		
		float muentry[4];
		for(int j = 0; j < 4; j++){
			muentry[j] = skroot_mu_.muboy_entpos[j][0];
		}
		
		makededx_(&muentry,
		          &skroot_mu_.muboy_dir,
		          &skchnl_.ihcab,
		          &skq_.qisk,
		          &skt_.tisk,
		          &geopmt_.xyzpm,
		          &skq_.nqisk,
		          &skroot_mu_.muboy_dedx);
		
		makededx_intg_(&muentry,
		               &skroot_mu_.muboy_dir,
		               &skroot_mu_.muboy_length,
		               &skchnl_.ihcab,
		               &skq_.qisk,
		               &skt_.tisk,
		               &geopmt_.xyzpm,
		               &sktqz_.nqiskz,
		               &skhead_.nrunsk,
		               &skroot_mu_.muboy_dedx,
		               &sktqz_.ihtiflz,
		               &skhead_.nevsk);
		
		// try BFF if muboy failed and muon energy >12MeV
		if(skroot_mu_.muboy_status == 1 && skroot_mu_.muboy_goodness < 0.4 &&
		   *std::max_element(muonsToRec[i].matchedParticleBSEnergy.begin(),
		   muonsToRec[i].matchedParticleBSEnergy.end()) > 12.){
			
			std::cout << "Starting muon BFF" << std::endl;
			float bffpos[3];
			float hpos[3];
			float bffgood;
			float modd;
			newmufit_(&bffpos, &hpos, &bffgood);
			modd = sqrt( pow((hpos[0]-bffpos[0]),2) + pow((hpos[1]-bffpos[1]),2) + pow((hpos[2]-bffpos[2]),2) );
			for(int j = 0; j < 3; j++){
				skroot_mu_.mubff_entpos[j] = bffpos[j];
				skroot_mu_.mubff_dir[j] = (hpos[j] - bffpos[j])/modd;
			}
			std::cout << "Finished muon BFF" << std::endl;
			
			skroot_mu_.mubff_goodness = bffgood;
			if(skroot_mu_.mubff_goodness > 0.3){
				for(int j = 0; j < 3; j++){
					muentry[j] = skroot_mu_.mubff_entpos[j];
				}
			}
			
			makededx_(&muentry,
			          &skroot_mu_.muboy_dir,
			          &skchnl_.ihcab,
			          &skq_.qisk,
			          &skt_.tisk,
			          &geopmt_.xyzpm,
			          &skq_.nqisk,
			          &skroot_mu_.muboy_dedx);
			
			makededx_intg_(&muentry,
			               &skroot_mu_.muboy_dir,
			               &skroot_mu_.muboy_length,
			               &skchnl_.ihcab,
			               &skq_.qisk,
			               &skt_.tisk,
			               &geopmt_.xyzpm,
			               &sktqz_.nqiskz,
			               &skhead_.nrunsk,
			               &skroot_mu_.muboy_dedx,
			               &sktqz_.ihtiflz,
			               &skhead_.nevsk);
		}
		
		// pass reconstructed muon info to output Tree branch variables
		skroot_set_mu_(&muWriterLUN,
		               skroot_mu_.muentpoint,
		               skroot_mu_.mudir,
		               &skroot_mu_.mutimediff,
		               &skroot_mu_.mugoodness,
		               &skroot_mu_.muqismsk,
		               &skroot_mu_.muyn,
		               &skroot_mu_.mufast_flag,
		               &skroot_mu_.muboy_status,
		               &skroot_mu_.muboy_ntrack,
		               skroot_mu_.muboy_entpos,
		               skroot_mu_.muboy_dir,
		               &skroot_mu_.muboy_goodness,
		               &skroot_mu_.muboy_length,
		               skroot_mu_.muboy_dedx,
		               skroot_mu_.mubff_entpos,
		               skroot_mu_.mubff_dir,
		               &skroot_mu_.mubff_goodness,
		               &skroot_mu_.muninfo,
		               skroot_mu_.muinfo);
		
		// for muons, only keep hits around 1.3us trigger
		delete_outside_hits_();
		
		// set header and tq info (epsecially updated hits)
		skroot_set_tree_(&muWriterLUN);
		
		// update branch variables w/ info about matches
		MatchedEvNums = muonsToRec[i].matchedParticleEvNum;
		MatchedTimeDiff = muonsToRec[i].matchedParticleTimeDiff;
		
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

/*float RelicMuonMatching::rollOver(unsigned long long int currentTime, unsigned long long int targetTime){
	unsigned long long int bitOne = 1;
	unsigned long long int tDiff;
	tDiff = currentTime - targetTime;
	tDiff = (bitOne << 47) + tDiff;
	
}*/

/*bool RelicMuonMatching::AddParticletoDeque(std::deque<ParticleCand>& addToThisDeque){
	unsigned long long int newTime = bitshiftTime(skheadqb_.it0xsk, skheadqb_.nevhwsk);
	ParticleCand newParticle;
	newParticle.EventNumber = myHeader->nevsk;
	newParticle.EventTime = newTime;
	newParticle.EntryNumber = rfmReader->GetEntryNumber();
	if(particleType == "LOWE"){
		newParticle.ReconEnergy = skroot_lowe_.bsenergy;
	} else{
		newParticle.ReconEnergy = 0.0;
	}
	addToThisDeque.push_back(newParticle);
	
	return true;
	}*/
	
/*unsigned long long int RelicMuonMatching::bitshiftTime(unsigned long long int t0Time, unsigned long long int hardwareTime){
	
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
}*/
