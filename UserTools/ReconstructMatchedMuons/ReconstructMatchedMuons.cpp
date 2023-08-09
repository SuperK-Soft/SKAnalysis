#include "ReconstructMatchedMuons.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "Algorithms.h"

#include <algorithm>

ReconstructMatchedMuons::ReconstructMatchedMuons():Tool(){}


bool ReconstructMatchedMuons::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	m_variables.Get("treeReaderName", treeReaderName);
	m_variables.Get("treeWriterLUN", treeWriterLUN);

	if(m_data->Trees.count(treeReaderName)==0){
	Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
	return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
		
	TreeManager* mgr = skroot_get_mgr(&treeWriterLUN);
	WriteTree = mgr->GetOTree();
	
	MatchedEvNumsBranch = WriteTree->Branch("MatchedEvNums", &MatchedEvNums);
	MatchedTimeDiffBranch = WriteTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	PIDBranch = WriteTree->Branch("MuonTag", &PID);
	
	return true;
}


bool ReconstructMatchedMuons::Execute(){
	muonsToRec = m_data->muonsToRec;
	
	if(! muonsToRec.size()) return true;
	
	myTreeReader->Get("HEADER", myHeader);
	
	currentEntry = myTreeReader->GetEntryNumber();
	
	for(int i = 0; i < muonsToRec.size(); i++){
		if(muonsToRec[i].EntryNumber != currentEntry){
			m_data->getTreeEntry(treeReaderName, muonsToRec[i].EntryNumber);
		}
		
		WriteInfo(muonsToRec[i]);
		
		int muyn_org, muynf;
		
		float mbentry [4];
		float mm_entry [36];
		int n_left;
		
		int currentEventNum = myHeader->nevsk;
		
		//Muon reconstruction developed by Tomoeda and Yamaguchi
		mfmuselect_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &skroot_mu_.mugoodness, &muyn_org);
		
		//muyn == 1 - good fit
		//muyn == 0 - bad fit
		
		if(muyn_org > 0){
			skroot_mu_.muyn = 1;
		}else if(muyn_org < 0){
			skroot_mu_.muyn = 0;
		}else{
			Log("Muyn_org returning as == 0. Not supported yet", v_error, verbosity);
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
		
		//Apply muboy
/*		muboy_(&currentEventNum,
			  &skroot_mu_.muboy_status,
			  &mbentry,
			  &skroot_mu_.muboy_dir,
			  &skroot_mu_.muboy_length,
			  &skroot_mu_.muboy_goodness,
			  &skroot_mu_.muboy_ntrack,
			  &mm_entry,
			  &n_left);*/
		
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
		skroot_set_mu_(&treeWriterLUN, 
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
			delete_outside_hits_();
			skroot_set_tree_(&treeWriterLUN);
			skroot_fill_tree_(&treeWriterLUN);
	}
	

	
	//return the previous entry so that no issues are caused with other tools
	if(currentEntry != myTreeReader->GetEntryNumber()){
		m_data->getTreeEntry(treeReaderName, currentEntry);
	}
	
	m_data->muonsToRec.clear();
	return true;
}

bool ReconstructMatchedMuons::Finalise(){
	
	return true;
}

bool ReconstructMatchedMuons::WriteInfo(ParticleCand Event){
	
	MatchedEvNums = Event.matchedParticleEvNum;
	MatchedEvNumsBranch->Fill();
	
	MatchedTimeDiff = Event.matchedParticleTimeDiff;
	MatchedTimeDiffBranch->Fill();
	
	PID = Event.PID;

	PIDBranch->Fill();
	
return true;
}


/*		//now fit muons using mufit
		lfclear_all_();

		//tentative routine to adjust saturation charge
		//fix_maxqisk_();
		
		//apply mufit
		std::cout << __FILE__ << __LINE__ << std::endl;
		lfmufit_sk4_();
		
		float muEntry[4] = {};
		
		for(int j = 0; j < 4; j++){
			muEntry[j] = skroot_mu_.muboy_entpos[0][j];
		}
		
/*		//once per run water transparency update
		if(myHeader->nrunsk != lastRun){
			int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];
			lfwater_(&days_to_run_start, &watert);
		}






/*		

		//recalculate dE/dx using Kirk's makededx.F TODO This needs importing into ToolFramework somehow!
		makededx_(&muEntry, &skroot_mu_.muboy_dir, &skchnl_.ihcab, &skq_.qisk, &skt_.tisk, &geopmt_.xyzpm, &skq_.nqisk,
				 &skhead_.nrunsk, &watert, &skroot_mu_.muinfo);
		
		//recalculate dE/dx using Scott's makededx.G
		makededx_intg_(&muEntry, &skroot_mu_.muboy_dir, &skroot_mu_.muboy_length, &skchnl_.ihcab, &skq_.qisk,
					  &skt_.tisk, &geopmt_.xyzpm, &skq_.nqisk, &skhead_.nrunsk, &skroot_mu_.muboy_dedx,
					  &sktqz_.ihtiflz, &skhead_.nevsk);
		
		float BFFPos[3];
		float HPos[3];
		float BFFGood;
		float modd;
		
		//apply BBF (single mu with goodness < 0,4 && relic energy > 12 MeV)
		
		if(skroot_mu_.muboy_status == 1 && skroot_mu_.muboy_goodness < 0.4 && *std::max_element(muonsToRec[i].matchedParticleEne.begin(), muonsToRec[i].matchedParticleEne.end()) > 12){
			std::cout << "Starting BFF..." << std::endl;
			newmufit_(&BFFPos, &HPos, &BFFGood);
			modd = sqrt(pow((HPos[0]-BFFPos[0]), 2) + pow((HPos[1]-BFFPos[1]), 2) + pow((HPos[2] - BFFPos[2]), 2));
			for(int j = 0; j < 3; j++){
				skroot_mu_.mubff_entpos[j] = BFFPos[j];
				skroot_mu_.mubff_dir[j] = (HPos[j] - BFFPos[j])/modd;
			}
			skroot_mu_.mubff_goodness = BFFGood;
			std::cout << "Ending BFF" << std::endl;
			// new dE/dx based on the brute force fitter track
			if(BFFGood > 0.3){
				for(int j =0; j < 3; j++){
					muEntry[j] = skroot_mu_.mubff_entpos[j];
				}
				makededx_(&muEntry, &skroot_mu_.muboy_dir, &skchnl_.ihcab, &skq_.qisk, &skt_.tisk, &geopmt_.xyzpm, &skq_.nqisk,
					 &skhead_.nrunsk, &watert, &skroot_mu_.muinfo);
				
				makededx_intg_(&muEntry, &skroot_mu_.muboy_dir, &skroot_mu_.muboy_length, &skchnl_.ihcab, &skq_.qisk,
								&skt_.tisk, &geopmt_.xyzpm, &skq_.nqisk, &skhead_.nrunsk, &skroot_mu_.muboy_dedx,
								&sktqz_.ihtiflz, &skhead_.nevsk);
			}
		}
*/