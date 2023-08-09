#include "WriteSpallCand.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "Algorithms.h"

WriteSpallCand::WriteSpallCand():Tool(){}


bool WriteSpallCand::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	
	m_variables.Get("treeReaderName", treeReaderName);
	m_variables.Get("treeWriterName", treeWriterName);
	
	myTreeWriter = m_data->Trees.at(treeWriterName);
	
	if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
	return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	lun = 20;
	
	TreeManager* mgr = skroot_get_mgr(&lun);
	WriteTree = mgr->GetOTree();
	
	MatchedEvNumsBranch = WriteTree->Branch("MatchedEvNums", &MatchedEvNums);
	MatchedTimeDiffBranch = WriteTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	PIDBranch = WriteTree->Branch("MuonTag", &PID);
	
	return true;
}


bool WriteSpallCand::Execute(){
	
	if(! m_data->writeOutRelics.size()) return true;
	
	int originalEntry = myTreeReader->GetEntryNumber();
	
	for(int writeEvent = 0; writeEvent < m_data->writeOutRelics.size(); writeEvent++){
		int currentEntry = myTreeReader->GetEntryNumber();
		if(m_data->writeOutRelics[writeEvent].EntryNumber != currentEntry){
			// don't use this anymore. m_data->GetTreeEntry(treeReaderName, entryNum) has now been added and 
			// can be used instead
			//myTreeReader->GetEntry(m_data->writeOutRelics[writeEvent]);
			m_data->getTreeEntry(treeReaderName, m_data->writeOutRelics[writeEvent].EntryNumber);
		}
		int io = 1;
		
		for(int i; i < branchestoSkip.size(); i++){
			skroot_zero_branch_(&lun, &io, branchestoSkip[i].c_str(), branchestoSkip[i].size());
		}
		
		WriteInfo(m_data->writeOutRelics[writeEvent]);
		
		skroot_lowe_ = m_data->writeOutRelics[writeEvent].LowECommon;
		
		std::cout << "THE BSENERGY IS:                       " << skroot_lowe_.bsenergy << std::endl;
		
		skroot_set_lowe_(&lun,
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
		
		
		/*	skroot_set_mu_(&lun, 
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
						skroot_mu_.muinfo);*/
		
		delete_outside_hits_();
		skroot_set_tree_(&lun);
		skroot_fill_tree_(&lun);
	}
	
	
	if(originalEntry != myTreeReader->GetEntryNumber()){
		m_data->getTreeEntry(treeReaderName, originalEntry);
	}
	
	m_data->writeOutRelics.clear();
	
	return true;
}


bool WriteSpallCand::Finalise(){
	
	return true;
}

bool WriteSpallCand::WriteInfo(ParticleCand Event){
	MatchedEvNums = Event.matchedParticleEvNum;
	
	for(int eventnum: MatchedEvNums){
		std::cout << "MATCHED EV NUM FOR LOW: " << eventnum << std::endl;
	}
	
	MatchedEvNumsBranch->Fill();
	
	MatchedTimeDiff = Event.matchedParticleTimeDiff;
	MatchedTimeDiffBranch->Fill();
	
	PID = Event.PID;

	PIDBranch->Fill();
	
return true;
}