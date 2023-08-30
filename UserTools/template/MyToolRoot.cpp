/* vim:set noexpandtab tabstop=4 wrap */
#include "MyToolRoot.h"

#include "Algorithms.h"
#include "Constants.h"

MyToolRoot::MyToolRoot():Tool(){}

bool MyToolRoot::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);            // how verbose to be
	std::string treeReaderName="";
	m_variables.Get("treeReaderName",treeReaderName);  // the name of the TTree to process
	
	// open the input TFile and TTree
	// ------------------------------
	if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+": Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	return true;
}


bool MyToolRoot::Execute(){
	
	Log(m_unique_name+" executing...",v_debug,m_verbose);
	
	// retrieve desired branches
	get_ok = GetBranchValues();
	
	return true;
}


bool MyToolRoot::Finalise(){
	
	return true;
}

int MyToolRoot::GetBranchValues(){
	int success = (
//		(myTreeReader->GetBranchValue("filename",filename))                         &&
//		(myTreeReader->GetBranchValue("gamma_time",gamma_time))
		true
	);
	
	return success;
}

