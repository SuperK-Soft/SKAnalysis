/* vim:set noexpandtab tabstop=4 wrap */
#include "MyToolM.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

MyToolM::MyToolM():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool MyToolM::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	
	return true;
}


bool MyToolM::Execute(){
	
	return true;
}


bool MyToolM::Finalise(){
	
	return true;
}

