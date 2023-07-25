/* vim:set noexpandtab tabstop=4 wrap */
#include "MyToolRoot.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

MyToolRoot::MyToolRoot():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool MyToolRoot::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	std::string treeReaderName="";
	m_variables.Get("treeReaderName",treeReaderName);  // the name of the TTree to process
	
	// open the input TFile and TTree
	// ------------------------------
	if(m_data->Trees.count(treeReaderName)==0){
		Log(toolName+": Failed to find TreeReader "+treeReaderName+" in DataModel!",0,0);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	return true;
}


bool MyToolRoot::Execute(){
	
	Log(toolName+" executing...",v_debug,verbosity);
	
	// retrieve desired branches
	get_ok = GetBranches();
	
	// process the data
	try{
		Analyse();
	}
	catch(std::exception& e){
		// catch any exceptions to ensure we always increment the event number
		// and load the next entry. This prevents us getting stuck in a loop
		// forever processing the same broken entry!
		Log(toolName+" encountered error "+e.what()+" during Analyse()",v_error,verbosity);
	}
	
	return true;
}


bool MyToolRoot::Finalise(){
	
	return true;
}

bool MyToolRoot::Analyse(){
	
	return true;
}

int MyToolRoot::GetBranches(){
	int success = (
//		(myTreeReader->GetBranchValue("filename",filename))                         &&
//		(myTreeReader->GetBranchValue("gamma_time",gamma_time))
		true
	);
	
	return success;
}

