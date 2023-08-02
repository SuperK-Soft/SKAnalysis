#include "WriteMatchedInfo.h"

WriteMatchedInfo::WriteMatchedInfo():Tool(){}


bool WriteMatchedInfo::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	
	int lun = 20;
	
	TreeManager* mgr = skroot_get_mgr(&lun);
	WriteTree = mgr->GetOTree();
	
	MatchedEvNumsBranch = WriteTree->Branch("MatchedEvNums", &MatchedEvNums);
	MatchedTimeDiffBranch = WriteTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	PIDBranch = WriteTree->Branch("MuonTag", &PID);
	
return true;
}


bool WriteMatchedInfo::Execute(){
	
	std::vector<ParticleCand> LoweEvents = m_data->writeOutRelics;
	std::vector<ParticleCand> Muons = m_data->muonsToRec;
	
	if(! LoweEvents.size() && ! Muons.size()){
		return true;
	}else{
		if(LoweEvents.size()){
			WriteInfo(LoweEvents);
		}
		if(Muons.size()){
			WriteInfo(Muons);
		}
	}
	
	LoweEvents.clear();
	Muons.clear();
	
	m_data->writeOutRelics.clear();
	m_data->muonsToRec.clear();
	
return true;
}


bool WriteMatchedInfo::Finalise(){
	
return true;
}

bool WriteMatchedInfo::WriteInfo(std::vector<ParticleCand> Events){
	
	for(ParticleCand Event: Events){
		MatchedEvNums = Event.matchedParticleEvNum;
		MatchedEvNumsBranch->Fill();
		
		MatchedTimeDiff = Event.matchedParticleTimeDiff;
		MatchedTimeDiffBranch->Fill();
		
		PID = Event.PID;
		std::cout << "THE EVENT PID IS" << Event.PID << std::endl;
		PIDBranch->Fill();
	}
	
return true;
}