/* vim:set noexpandtab tabstop=4 wrap */
#include "TreeReaderDemo.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <time.h>  // for 'struct tm'
#include <numeric> // for std::accumulate

TreeReaderDemo::TreeReaderDemo():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool TreeReaderDemo::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("treeReaderName",treeReaderName);
	
	if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	return true;
}


bool TreeReaderDemo::Execute(){
	
	// retrieve variables from TreeReader
	GetBranchValues();
	
	// print event info, just for demonstration
	// ----------------------------------------
	
	// run and subrun number
	std::cout<<std::endl<<"Event "<<myTreeReader->GetEntryNumber()
			 <<" was from run "<<header->nrunsk<<", subrun "<<header->nsubsk<<std::endl;
	
	// event time and date
	tm rundate = {0};
	rundate.tm_year = header->ndaysk[0];
	rundate.tm_mon = header->ndaysk[1] - 1;
	rundate.tm_mday = header->ndaysk[2];
	rundate.tm_hour = header->ntimsk[0];
	rundate.tm_min = header->ntimsk[1];
	rundate.tm_sec = header->ntimsk[2];
	time_t runtime = mktime(&rundate);         // need to use mktime to derive day of week
	std::string timestring = ctime(&runtime);  // format timestamp into a string
	timestring.pop_back();                     // drop trailing newline
	std::cout<<"The event occurred at "<<timestring<<std::endl;
	
	// total hits and charge
	double total_charge = std::accumulate(tqreal->Q.begin(),tqreal->Q.end(),0);
	std::cout<<"The event saw "<<tqreal->nhits
			 <<" PMTs hit with a total charge of "<<total_charge<<" nC"<<std::endl;
	
	bool has_secondary_info = myTreeReader->Get("SECONDARY",secondaries);
	if(has_secondary_info){
		std::cout<<"calling secondary print"<<std::endl;
		secondaries->Print();
	} else {
		Log(toolName+" No SECONDARY branch - cannot print secondary information at this time.",v_debug,verbosity);
	}
	
	// reconstructed energy and vertex
	std::cout<<"Low-E reconstruction indicated an energy of "<<loweinfo->bsenergy<<" MeV"<<std::endl;
	basic_array<float> bsvertex(loweinfo->bsvertex);
	std::cout<<"with a primary vertex at {X,Y,Z,T} = {";
	for(auto&& vtxcomp : bsvertex){ std::cout<<vtxcomp<<", "; };
	std::cout<<"\b\b}"<<std::endl;
	
	return true;
}


bool TreeReaderDemo::Finalise(){
	return true;
}

bool TreeReaderDemo::GetBranchValues(){
	// retrieve variables from branches
	// include a list here of all branches you wish to retrieve.
	bool success = 
	(myTreeReader->Get("HEADER", header)) &&      // TreeReader::Get("branchname", variable);
	(myTreeReader->Get("LOWE", loweinfo)) &&
	(myTreeReader->Get("TQREAL", tqreal));
	
	return success;
}
