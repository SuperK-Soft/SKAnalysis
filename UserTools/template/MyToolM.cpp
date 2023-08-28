/* vim:set noexpandtab tabstop=4 wrap */
#include "MyToolM.h"

#include "Algorithms.h"
#include "Constants.h"

MyToolM::MyToolM():Tool(){}

bool MyToolM::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);            // how verbose to be
	
	return true;
}


bool MyToolM::Execute(){
	
	return true;
}


bool MyToolM::Finalise(){
	
	return true;
}

