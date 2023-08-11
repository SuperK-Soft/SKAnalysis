#include "HitsCut.h"

HitsCut::HitsCut():Tool(){
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}


bool HitsCut::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	m_variables.Get("hitLimit", hitLimit);
	
	return true;
}


bool HitsCut::Execute(){
	
	totalHits = skq_.nqisk;

	
	bool muon = false;
	m_data->vars.Get("newMuon", muon);
	
	if(muon){
		return true;
	}
	
	if(totalHits > hitLimit){
		m_data->vars.Set("Skip", true);
	}
	
	if(totalHits<=0){
		Log(m_unique_name+": Warning. No hits for this event were found in the skq_ common block",v_error,m_verbose);
	}
	
	totalHits = -1;
	
	return true;
}


bool HitsCut::Finalise(){
	
	return true;
}
