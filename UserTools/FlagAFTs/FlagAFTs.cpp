#include "FlagAFTs.h"
#include <bitset>

FlagAFTs::FlagAFTs():Tool(){}


bool FlagAFTs::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	
	return true;
}


bool FlagAFTs::Execute(){
	
	std::bitset<32> triggerID(skhead_.idtgsk);
	
	/*
	// get the eventType from the end of the last loop execution
	EventType lastEventType;
	get_ok = m_data->vars.Get("eventType", lastEventType);
	if(!get_ok){
		std::cerr<<m_unique_name<<" Error getting eventType!"<<std::endl;
	}
	*/
	
	// set the EventType of the current event to AFT if applicable
	//std::cout<<"Trigger bits for event: "<<skhead_.nevsk<<": "<<GetTriggerNames(skhead_.idtgsk)<<std::endl;
	if(triggerID.test(TriggerType::AFT_or_Cal) || triggerID.test(TriggerType::AFT)){
		m_data->vars.Set("eventType", EventType::AFT);
		//std::cout<<m_unique_name<<" flagging event "<<skhead_.nevsk<<" as AFT"<<std::endl;
		//std::cout<<" -- last event type was "<<lastEventType<<std::endl;
	} else {
		m_data->vars.Set("eventType", EventType::Unknown);
		//std::cout<<m_unique_name<<" flagging event "<<skhead_.nevsk<<" as unknown"<<std::endl;
		//std::cout<<" -- last event type was "<<lastEventType<<std::endl;
	}
	
	
	return true;
	
	/*
	
	// set flag for WriteSkEvent to save AFT triggers
	// along with which file to save to based on the preceding event type
	if(triggerID.test(TriggerType::AFT)){
		
		bool dosave=true;
		m_data->vars.Set("saveEvent", dosave);
		
		// if AFT follows a muon, save it to the muon tree.
		// (we'll need to scan it for SLE triggers for the neutron cloud cut)
		
		// Otherwise save it to the relic tree.
		// (we'll need to scan it for neutron tagging with the BDT)
		
		std::string treeReaderName;
		if(lastEventType==EventType::Muon){
			treeReaderName = "muWriter";
		} else {
			treeReaderName = "relicWriter";
		}
		
		int LUN = m_data->GetLUN(treeReaderName);
		if(LUN<0){
			Log(m_unique_name+" Error! Failed to find TreeReader "+treeReaderName+
			    " in DataModel!",v_error,m_verbose);
			m_data->vars.Set("StopLoop",1); // fatal error
			return false;
		}
		
		m_data->vars.Set("WriteSkEventLUN", LUN);
		
		return true;
	}
	*/
	
	return true;
}


bool FlagAFTs::Finalise(){
	
	return true;
}
