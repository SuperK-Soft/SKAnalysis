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
	
	if(!ParseOptions(configfile)) return false;
	
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
		m_data->AddCut(selectorName, m_unique_name, description,false);
	}
	
	return true;
}


bool SkipTriggers::Execute(){
	
	std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
	
	// debug prints
	if(m_verbose >= v_debug) PrintTriggers();
	
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

bool SkipTriggers::ParseOptions(std::string configfile){
	
	std::string allowedTriggersString="";
	std::string skippedTriggersString="";
	
	std::ifstream ifile(configfile);
	if(!ifile.is_open()){
		Log(m_unique_name+" Error! Failed to open config file "+configfile,v_error,m_verbose);
		return false;
	}
	std::string line;
	std::string key;
	std::stringstream ss;
	while(getline(ifile, line)){
		if(line.empty()) continue;
		if(line[0]=='#') continue;
		ss.clear();
		ss.str(line);
		if(!(ss >> key)) continue;
		if(key[0]=='#') continue;
		if(key=="skippedTriggers"){
			skippedTriggersString = line.substr(line.find("skippedTriggers")+key.length(),std::string::npos);
		} else if(key=="allowedTriggers"){
			allowedTriggersString = line.substr(line.find("allowedTriggers")+key.length(),std::string::npos);
		}
	}
	ifile.close();
	
	if(allowedTriggersString!="" && skippedTriggersString!=""){
		// seems too much like user error
		Log(m_unique_name+" Error! Please specify either allowed or skipped triggers, not both",v_error,m_verbose);
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
			Log(m_unique_name+" error parsing allowed trigger '"+next_trig+"'",v_error,m_verbose);
			return false;
		}
	}
	// same for skippedTriggersString
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
			Log(m_unique_name+" error parsing skipped trigger '"+next_trig+"'",v_error,m_verbose);
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
	Log(m_unique_name+" Trigger word for this event: "+trigger_bits.to_string()
	   +" = "+GetTriggerNames(skhead_.idtgsk),v_debug,m_verbose);
}
