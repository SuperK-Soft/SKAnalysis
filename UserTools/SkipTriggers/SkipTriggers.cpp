#include "SkipTriggers.h"

#include "fortran_routines.h" // for access to skhead_.idtgsk
#include "Constants.h"
#include "MTreeSelection.h"
#include <bitset>

SkipTriggers::SkipTriggers():Tool(){}


bool SkipTriggers::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	
	if(!ParseOptions()) return false;
	
	// if saving the passing events to an MTreeSelection...
	get_ok = m_variables.Get("selectorName", selectorName);
	if(get_ok){
		std::string description="";
		if(allowedTriggers.size()){
			description="allowed triggers: {";
			for(int i=0; i<allowedTriggers.size(); ++i){
				if(i>0) description+=", ";
				description+=std::to_string(allowedTriggers.at(i));
			}
		} else {
			description="rejected triggers: {";
			for(int i=0; i<skippedTriggers.size(); ++i){
				if(i>0) description+=", ";
				description+=std::to_string(skippedTriggers.at(i));
			}
		}
		m_data->AddCut(selectorName, m_unique_name, description);
	}
	
	return true;
}


bool SkipTriggers::Execute(){
	
	std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
	
	// debug prints
	if(m_verbose >= v_debug) PrintTriggerBits();
	
	bool skipit=false;
	// specify only the types we do want
	if(allowedTriggers.size()){
		skipit = true;
		for(size_t bit_i=0; bit_i<allowedTriggers.size(); ++bit_i){
			int test_bit = allowedTriggers.at(bit_i);
			if(trigger_bits.test(test_bit)){
				skipit=false;
			}
		}
	}
	
	// skip all events matching trigger types in skippedTriggers
	for(size_t bit_i=0; bit_i<skippedTriggers.size(); ++bit_i){
		int test_bit = skippedTriggers.at(bit_i);
		if(trigger_bits.test(test_bit)){
			skipit=true; // skip this event
		}
	}
	
	if(skipit) m_data->vars.Set("Skip", true); // skip this event
	
	if(!selectorName.empty() && !skipit) m_data->AddPassingEvent(selectorName, m_unique_name);
	
	return true;
}


bool SkipTriggers::Finalise(){
	
	return true;
}

bool SkipTriggers::ParseOptions(){
	
	std::string allowedTriggersString="";
	std::string skippedTriggersString="";
	m_variables.Get("skippedTriggers", skippedTriggersString);
	m_variables.Get("allowedTriggers", allowedTriggersString);
	if(allowedTriggersString!="" && skippedTriggerString!=""){
		// seems too much like user error
		Log(m_unique_name+" Error! Please specify either allowed or skipped triggers, not both",v_error,verbosity);
		return false;
	}
	std::stringstream allowedTriggersSS(allowedTriggersString);
	std::stringstream skippedTriggersSS(skippedTriggersString);
	// parse the allowedTriggersString for allowed triggers
	std::string next_trig;
	while(allowedTriggersSS >> next_trig){
		if(next_trig[0]=='#') break; // trailing comments
		if(TriggerNameToID(next_trig)>=0){
			allowedTriggers.push_back(TriggerNameToID(next_trig));
			continue;
		}
		try{
			int next_bit = stoi(next_trig);
			allowedTriggers.push_back(next_bit);
		} catch (...) {
			Log(m_unique_name+" error parsing allowed trigger '"+next_trig+"'",v_error,verbosity);
			return false;
		}
	}
	// same for skippedTriggerString
	while(skippedTriggersSS >> next_trig){
		if(next_trig[0]=='#') break; // trailing comments
		if(TriggerNameToID(next_trig)>=0){
			skippedTriggers.push_back(TriggerNameToID(next_trig));
			continue;
		}
		try{
			int next_bit = stoi(next_trig);
			skippedTriggers.push_back(next_bit);
		} catch (...) {
			Log(m_unique_name+" error parsing skipped trigger '"+next_trig+"'",v_error,verbosity);
			return false;
		}
	}
	
	if(skippedTriggers.empty() && allowedTriggers.empty()){
		Log(m_unique_name+" Error! No 'skippedTriggers' nor 'allowedTriggers' in config file!",v_error,m_verbose);
		return false;
	}
	
	return true;
}


void SkipTriggers::PrintTriggers(){
	std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
	
	Log(m_unique_name+" Trigger word for the active entry is: "
		+trigger_bits.to_string(),v_debug,m_verbose);
	if(m_verbose>(v_debug+1)){
		for(int i=0; i<(sizeof(int)*8); ++i){
			if(trigger_bits.test(i)) std::cout<<"Trigger "<<TriggerIDToName(i)<<" ("<<i<<") set"<<std::endl;
		}
	}
}
