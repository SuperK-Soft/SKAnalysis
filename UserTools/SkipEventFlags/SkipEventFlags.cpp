#include "SkipEventFlags.h"

#include "fortran_routines.h" // for access to skhead_.ifevsk
#include "Constants.h"
#include "MTreeSelection.h"
#include <bitset>

SkipEventFlags::SkipEventFlags():Tool(){}


bool SkipEventFlags::Initialise(std::string configfile, DataModel &data){
	
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
		if(allowedFlags.size()){
			description="allowed flags: {";
			for(int i=0; i<allowedFlags.size(); ++i){
				if(i>0) description+=", ";
				description+=std::to_string(allowedFlags.at(i));
			}
		} else {
			description="rejected flags: {";
			for(int i=0; i<skippedFlags.size(); ++i){
				if(i>0) description+=", ";
				description+=std::to_string(skippedFlags.at(i));
			}
		}
		m_data->AddCut(selectorName, m_unique_name, description,false);
	}
	
	return true;
}


bool SkipEventFlags::Execute(){
	
	std::bitset<sizeof(int)*8> flags_bits = skhead_.ifevsk;
	
	// debug prints
	if(m_verbose >= v_debug) PrintFlags();
	
	bool skipit=false;
	// specify only the types we do want
	if(allowedFlags.size()){
		skipit = true;
		for(size_t bit_i=0; bit_i<allowedFlags.size(); ++bit_i){
			int test_bit = allowedFlags.at(bit_i);
			if(flags_bits.test(test_bit)){
				skipit=false;
			}
		}
	}
	
	// skip all events matching flags types in skippedFlags
	for(size_t bit_i=0; bit_i<skippedFlags.size(); ++bit_i){
		int test_bit = skippedFlags.at(bit_i);
		if(flags_bits.test(test_bit)){
			skipit=true; // skip this event
		}
	}
	
	if(skipit) m_data->vars.Set("Skip", true); // skip this event
	
	if(!selectorName.empty() && !skipit) m_data->AddPassingEvent(selectorName, m_unique_name);
	
	return true;
}


bool SkipEventFlags::Finalise(){
	
	return true;
}

bool SkipEventFlags::ParseOptions(std::string configfile){
	
	std::string allowedFlagsString="";
	std::string skippedFlagsString="";
	
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
		if(key=="skippedFlags"){
			skippedFlagsString = line.substr(line.find("skippedFlags")+key.length(),std::string::npos);
		} else if(key=="allowedFlags"){
			allowedFlagsString = line.substr(line.find("allowedFlags")+key.length(),std::string::npos);
		}
	}
	ifile.close();
	
	std::stringstream allowedFlagsSS(allowedFlagsString);
	std::stringstream skippedFlagsSS(skippedFlagsString);
	// parse the allowedFlagsString for allowed flags
	std::string next_flag;
	while(allowedFlagsSS >> next_flag){
		if(next_flag[0]=='#') break; // trailing comments
		if(constants::string_to_flag_SKIV.count(next_flag)){
			allowedFlags.push_back(constants::string_to_flag_SKIV.at(next_flag));
			continue;
		}
		if(constants::string_to_flag_SKI_III.count(next_flag)){
			allowedFlags.push_back(constants::string_to_flag_SKI_III.at(next_flag));
			continue;
		}
		try{
			int next_bit = stoi(next_flag);
			allowedFlags.push_back(next_bit);
		} catch (...) {
			Log(m_unique_name+" error parsing allowed flag '"+next_flag+"'",v_error,m_verbose);
			return false;
		}
	}
	// same for skippedFlagString
	while(skippedFlagsSS >> next_flag){
		if(next_flag[0]=='#') break; // trailing comments
		if(constants::string_to_flag_SKIV.count(next_flag)){
			skippedFlags.push_back(constants::string_to_flag_SKIV.at(next_flag));
			continue;
		}
		if(constants::string_to_flag_SKI_III.count(next_flag)){
			skippedFlags.push_back(constants::string_to_flag_SKI_III.at(next_flag));
			continue;
		}
		try{
			int next_bit = stoi(next_flag);
			skippedFlags.push_back(next_bit);
		} catch (...) {
			Log(m_unique_name+" error parsing skipped flag '"+next_flag+"'",v_error,m_verbose);
			return false;
		}
	}
	
	if(skippedFlags.empty() && allowedFlags.empty()){
		Log(m_unique_name+" Error! No 'skippedFlags' nor 'allowedFlags' in config file!",v_error,m_verbose);
		return false;
	}
	
	return true;
}


void SkipEventFlags::PrintFlags(){
	std::bitset<sizeof(int)*8> flags_bits = skhead_.ifevsk;
	Log(m_unique_name+" Event flags for this event: "+flags_bits.to_string()
	   +" = "+GetEventFlagNames(skhead_.ifevsk),v_debug,m_verbose);
}
