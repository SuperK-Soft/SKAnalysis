/* vim:set noexpandtab tabstop=4 wrap */
#include "RootReadTest.h"

#include "type_name_as_string.h"

#include <iostream>
#include <vector>
#include <sstream>

RootReadTest::RootReadTest():Tool(){}


bool RootReadTest::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	m_variables.Get("inputFile",inputFile);
	m_variables.Get("treeName",treeName);
	m_variables.Get("testFileType",testFileType);
	m_variables.Get("maxEvents",maxEvents);
	
	myTreeReader.Load(inputFile, treeName);  // will pre-load entry 0
	myTreeReader.SetVerbosity(1);
	entrynum=0;
	
	return true;
}

bool RootReadTest::Execute(){
	std::cout<<"ReadRootTest getting entry "<<entrynum<<std::endl;
	
	// ttree entry is already pre-loaded,
	// so just retrieve the desired branches
	if(testFileType=="official_ntuple"){
		get_ok = GetBranchesNtuple();
	}
	else if(testFileType=="SKROOT"){
		get_ok = GetBranchesSKROOT();
	}
	
	// process the data
	if(testFileType=="official_ntuple"){
		CheckEntryNtuple();
	}
	else if(testFileType=="SKROOT"){
		CheckEntrySKROOT();
	}
	
	// move to next entry
	entrynum++;
	// check if we've hit the user-requested entry limit
	if((maxEvents>0)&&(entrynum==maxEvents)){
		std::cout<<"hit max events, setting StopLoop"<<std::endl;
		m_data->vars.Set("StopLoop",1);
		return 1;
	}
	
	// otherwise, pre-load the next ttree entry.
	// we do this now so that we can prevent the next iteration (by setting StopLoop) should we
	// find that we've hit the end of the TChain, or there is an error reading the next entry.
	// This ensures downstream tools always have valid data to process.
	get_ok = ReadEntry(entrynum);
	if(get_ok==0){
		return 1; // end of file
	} else if (get_ok<0){
		return 0; // read error
	}
	
	return get_ok;
}


bool RootReadTest::Finalise(){
	
	return true;
}

int RootReadTest::ReadEntry(long entry_number){
	// load next entry data from TTree
	int bytesread = myTreeReader.GetEntry(entrynum);
	
	// stop loop if we ran off the end of the tree
	if(bytesread<1&&bytesread>-3){
		std::cout<<"ReadRootTest Hit end of input file, stopping loop"<<std::endl;
		m_data->vars.Set("StopLoop",1);
	}
	// stop loop if we had an error of some kind
	else if(bytesread<0){
		 if(bytesread==-1) std::cerr<<"ReadRootTest IO error loading next input entry!"<<std::endl;
		 if(bytesread==-10) std::cerr<<"ReadRootTest AutoClear error loading next input entry!"<<std::endl;
		 if(bytesread <-2) std::cerr<<"ReadRootTest Unknown error loading next input entry!"<<std::endl;
		 m_data->vars.Set("StopLoop",1);
	}
	
	return bytesread;
}

int RootReadTest::GetBranchesNtuple(){
	int success = 
	(myTreeReader.GetBranchValue("nscndprt", n_secondaries_2)) &&
	(myTreeReader.GetBranchValue("iprtscnd", secondary_PDG_code_2)) &&
	(myTreeReader.GetBranchValue("vtxscnd", secondary_start_vertex_2));
	
	return success;
}

int RootReadTest::GetBranchesSKROOT(){
	int success = 
	(myTreeReader.GetBranchValue("MC", mc_info)) &&
	(myTreeReader.GetBranchValue("HEADER", file_header));
	
	return success;
}


// Official Ntuple
// ---------------
int RootReadTest::CheckEntryNtuple(){
	std::cout<<"we had "<<n_secondaries_2<<" secondaries"<<std::endl;
	
	std::cout<<"we had "<<secondary_PDG_code_2.size()<<" secondary PDG codes: {";
	for(auto&& asecondary : secondary_PDG_code_2){
		std::cout<<asecondary<<", ";
	}
	std::cout<<"\b\b}"<<std::endl;
	
	std::cout<<"we had "<<secondary_start_vertex_2.size()<<" secondary start vertices: {";
	for(int i=0; i<secondary_start_vertex_2.size(); ++i){
		auto&& avertex = secondary_start_vertex_2.at(i);
		std::cout<<"[";
		for(auto&& aval : avertex){
			std::cout<<aval<<", ";
		}
		std::cout<<"\b\b], ";
	}
	std::cout<<"\b\b}"<<std::endl;
	return 1;
}

// SKROOT
// ------
int RootReadTest::CheckEntrySKROOT(){
	std::cout<<"MCInfo is at "<<mc_info<<std::endl;
	std::cout<<"we had "<<mc_info->nvc<<" primaries"<<std::endl;
	
	for(int primary_i=0; primary_i<mc_info->nvc; ++primary_i){
		std::cout<<"primary "<<primary_i<<" had pdg code "<<mc_info->ipvc[primary_i]
				 <<" and initial momentum ("<<mc_info->pvc[primary_i][0]
				 <<", "<<mc_info->pvc[primary_i][1]<<", "<<mc_info->pvc[primary_i][2]<<"); ";
	}
	
	std::cout<<"Header is at "<<file_header<<", and indicates nevsk "
			 <<file_header->nevsk<<" and swtrig_id "<<file_header->swtrg_id<<std::endl;
	
	return 1;
}

