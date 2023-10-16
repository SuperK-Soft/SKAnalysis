#include "ODCut.h"

#include "fortran_routines.h"
#include "SK_helper_functions.h"
#include "Constants.h"

#include <bitset>

ODCut::ODCut():Tool(){}


bool ODCut::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	m_variables.Get("verbosity",m_verbose);
	std::string treeReaderName, selectionsName;
	m_variables.Get("treeReaderName",treeReaderName);
	m_variables.Get("hitsThreshold", hitsThreshold);
	m_variables.Get("windowMinT",windowMinT);
	m_variables.Get("windowMaxT",windowMaxT);
	
	get_ok = m_variables.Get("selectorName",selectorName);
	if(get_ok){
		std::string description = "cut events with the OD trigger bit set in the primary trigger, or more than "
		                          +toString(hitsThreshold)+" hits within the time range "
		                          +toString(windowMinT)+" to "+toString(windowMaxT)
		                          +" ns from the primary trigger";
		m_data->AddCut(selectorName, m_unique_name, description, true);
	}
	
	return true;
}


bool ODCut::Execute(){
	
	// check if this is a relic candidate (low-e event), and if not bypass this cut
	EventType eventType;
	m_data->vars.Get("eventType", eventType);
	if(eventType!=EventType::LowE) return true;
	
	// otherwise check if the OD software trigger bit has been set, and if so throw the event out
	std::bitset<32> triggerID{skhead_.idtgsk};
	if(triggerID.test(4)){
		Nskipped++;
		m_data->vars.Set("Skip", true);
		return true;
	}
	
	// as a more stringent check we can also scan for OD activity
	int ODHits = sktqaz_.nhitaz;
	
	//check for number of OD PMT hits within a 500 - 1300 ns range
	// FIXME is this equivalent to `lfnhita`???
	float ODHitTime;
	int ODHitsInWindow = 0;
	for(int hitPMT; hitPMT < ODHits; hitPMT++){
		ODHitTime = sktqaz_.taskz[hitPMT];
		if(ODHitTime >= windowMinT && ODHitTime <= windowMaxT){
			ODHitsInWindow++;
		}
	}
	
	//if there are more than `hitsThreshold` hits in this range, throw out the event as a relic candidate
	if(ODHitsInWindow > hitsThreshold){
		Nskipped++;
		m_data->vars.Set("Skip", true);
	}
	
	// record whether this event passed the cut
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, m_unique_name, ODHitsInWindow);
	
	return true;
}


bool ODCut::Finalise(){
	
	std::cout << "Number of events skipped due to ODCut: " << Nskipped << std::endl;
	
	return true;
}
