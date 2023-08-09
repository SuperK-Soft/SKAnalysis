#include "RelicCandidates.h"

RelicCandidates::RelicCandidates():Tool(){}


bool RelicCandidates::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("treeReaderName",treeReaderName);
	
	// if getting data from TTree, check the TreeReader
	if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	return true;
}


bool RelicCandidates::Execute(){
	//	std::cout << "Made it to the reliccandidates tool" << std::endl;
	bool muonEventFlag = false;
	
	m_data->vars.Get("newMuon", muonEventFlag);
	
	//if the chain has got to this point without either setting the newMuon flag or skipping the events then it
	//is considered to be a relic candidate, so set the newRelic flag
	if(!muonEventFlag){
		m_data->vars.Set("newRelic", true);
	}
	
	
	return true;
}


bool RelicCandidates::Finalise(){
	
	return true;
}
