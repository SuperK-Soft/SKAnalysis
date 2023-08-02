#include "HitsCut.h"
#include "fortran_routines.h"

HitsCut::HitsCut():Tool(){}


bool HitsCut::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	
	return true;
}


bool HitsCut::Execute(){
	
	totalHits = skq_.nqisk;

	
	bool muon = false;
	m_data->vars.Get("newMuon", muon);
	
	if(muon){
		return true;
	}
	
	if(totalHits > 999){
		m_data->vars.Set("Skip", true);
	}
	
	if(totalHits == 0 || totalHits == -1){
		std::cout << "Warning: No hits for this event were found in the SKTQZ common block" << std::endl;
	}
	
	totalHits = -1;
	
	//THIS IS TEMPORARY PLEASE DON'T===
	//FORGET TO GET RID OF IT JACK=====
	//=================================
	//	m_data->vars.Set("newMuon", false);
	//=================================
	//=================================
	//=================================
	
	return true;
}


bool HitsCut::Finalise(){
	
	return true;
}
