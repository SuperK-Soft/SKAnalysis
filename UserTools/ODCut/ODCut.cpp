#include "ODCut.h"
#include "Constants.h"

#include "fortran_routines.h"
#include "softtrg_tbl.h"
#include "SK_helper_functions.h"

#include <bitset>

ODCut::ODCut():Tool(){}


bool ODCut::Initialise(std::string configfile, DataModel &data){
	
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


bool ODCut::Execute(){
	
	bool muonEventFlag = false;
	m_data->vars.Get("newMuon", muonEventFlag);
	
	//check if the muon flag has been set. if it has skip the event as a possible relic candidate
	//if not then check to see if the OD software trigger has been set. If it has then throw out the event as a
	//possible relic candidate and set the skip variable.
	if(muonEventFlag){
		return true;
	}
	else{
		myTreeReader->Get("HEADER", myHeader);
		triggerID = myHeader->idtgsk;
		if(triggerID.test(4)){
			Nskipped++;
			m_data->vars.Set("Skip", true);
			return true;
		}
	}
	
	//if not a muon and no OD trigger set then need to check the number of OD hits in a time window
	//get the number of PMT hits in the OD for this event
	int ODHits = sktqaz_.nhitaz;
	
	//check for number of OD PMT hits within a 500 - 1300 ns range
	float ODHitTime;
	int ODHitsInWindow = 0;
	for(int hitPMT; hitPMT < ODHits; hitPMT++){
		ODHitTime = sktqaz_.taskz[hitPMT];
		if(ODHitTime >= 500. && ODHitTime <= 1300.){
			ODHitsInWindow++;
		}
	}
	
	//if there are more than 20 hits in the 500 - 1300 ns range then throw out the event as a relic candidate
	if(ODHitsInWindow > 20){
		Nskipped++;
		m_data->vars.Set("Skip", true);
	}
	
	return true;
}


bool ODCut::Finalise(){
	
	std::cout << "Number of events skipped due to ODCut: " << Nskipped << std::endl;
	
	return true;
}
