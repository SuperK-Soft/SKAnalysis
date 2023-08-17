#include "IDChargeCut.h"

IDChargeCut::IDChargeCut():Tool(){}


bool IDChargeCut::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	m_variables.Get("hitLimit", hitLimit);
	
	get_ok = m_variables.Get("selectorName", selectorName);
	if(get_ok){
		std::string description = "cut events with nqisk > "+toString(hitLimit);
		m_data->AddCut(selectorName, m_unique_name, description);
	}
	
	return true;
}


bool IDChargeCut::Execute(){
	
	totalHits = skq_.nqisk;
	
	if(totalHits<1) Log(m_unique_name+": Warning! skq_ common block is empty!",v_error,m_verbose);
	
	bool muon = false;
	m_data->vars.Get("newMuon", muon);
	
	if(muon) return true;
	
	if(totalHits > hitLimit){
		m_data->vars.Set("Skip", true);
	}
	
	if(selectorName) m_data->AddPassingEvent(selectorName, m_unique_name);
	
	return true;
}


bool IDChargeCut::Finalise(){
	
	return true;
}
