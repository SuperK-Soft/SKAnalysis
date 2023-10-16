#include "CutRecorder.h"

CutRecorder::CutRecorder():Tool(){}


bool CutRecorder::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	std::string selectorName, treeReaderName, selectionsFile, distributionsFile;
	m_variables.Get("selectorName",selectorName);            // a name for downstream Tools
	m_variables.Get("selectionsFile",selectionsFile);        // output file to generate
	m_variables.Get("distributionsFile",distributionsFile);  // output file to generate
	m_variables.Get("treeReaderName",treeReaderName);        // TreeReader for input data
	
	MTreeReader* thereader = m_data->Trees.at(treeReaderName);
	myTreeSelections.SetTreeReader(thereader);
	myTreeSelections.MakeOutputFile(selectionsFile, distributionsFile);
	
	m_data->Selectors.emplace(selectorName, &myTreeSelections);
	
	myTreeSelections.AddCut("all", "all events",false);
	
	return true;
}


bool CutRecorder::Execute(){
	
	myTreeSelections.AddPassingEvent("all");
	++nentries;
	
	return true;
}


bool CutRecorder::Finalise(){
	
	Log(m_unique_name+" event counts trace: ", v_warning, m_verbose);
	//myTreeSelections.SetEntries(nentries);
	myTreeSelections.PrintCuts();
	
	return true;
}

CutRecorder::~CutRecorder(){
	// write out the event numbers that passed each cut
	myTreeSelections.Write();
}
