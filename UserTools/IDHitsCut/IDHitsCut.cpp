#include "IDHitsCut.h"

IDHitsCut::IDHitsCut():Tool(){}


bool IDHitsCut::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	m_variables.Get("hitLimit", hitLimit);
	
	get_ok = m_variables.Get("selectorName", selectorName);
	if(get_ok){
		std::string description = "cut events with nqisk > "+toString(hitLimit);
		m_data->AddCut(selectorName, m_unique_name, description,true);
	}
	
	return true;
}


bool IDHitsCut::Execute(){
	
	totalHits = skq_.nqisk;
	
	if(totalHits<1) Log(m_unique_name+": Warning! skq_ common block is empty!",v_error,m_verbose);
	
	EventType eventType;
	m_data->vars.Get("eventType", eventType);
	if(eventType!=EventType::LowE) return true;
	
	if(totalHits > hitLimit){
		m_data->vars.Set("Skip", true);
	}
	
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, m_unique_name,totalHits);
	
	return true;
}


bool IDHitsCut::Finalise(){
	
	return true;
}
