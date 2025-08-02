/* vim:set noexpandtab tabstop=4 wrap */
#include "SK2p2MeV_ntag.h"
#include "SK2p2MeV.h"
#include "SK2p2MeV_ambe.h"
#include "SK2p2MeV_mc.h"
#include "SK2p2MeV_relic.h"
#include "SK2p2MeV_t2k.h"
#include "SK2p2MeV_merged.h"

#include "type_name_as_string.h"

SK2p2MeV_ntag::SK2p2MeV_ntag():Tool(){}

// -----------------------------

// Tool to apply the SK2p2MeV neutron tagging algorithm.
// There are a few derived classes that perform additional initialization of the data handed
// to the actual tagging algorithm (SK2p2MeV::NeutronSearch())
// Results are written to an output TTree.
// TODO: move results to the DataModel, where they can be accessed by downstream tools
// TODO: break up the SK2p2MeV::NeutronSearch() algorithm into its component steps,
// such that they can be swapped in with alternatives and intermediate values can be inspected.

bool SK2p2MeV_ntag::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	std::string treeReaderName;
	std::string derived_classname;
	std::string outFile = "SK2p2MeV.root";
	m_variables.Get("verbosity",m_verbose);                   // how verbose to be
	m_variables.Get("treeReaderName",treeReaderName);         // reader for input tree
	m_variables.Get("derived_classname", derived_classname);  // derived class of SK2p2MeV
	m_variables.Get("outFile",outFile);                       // output file for TTree
	m_variables.Get("isWIT",isWIT);                           // is the input file a WIT file?
	
	
	if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	// make output file
	if(outFile == nullptr) outFile = "SK2p2MeV_output.root";
	fout = new TFile(outFile.c_str(), "RECREATE");
	
	// make the SK2p2MeV class instance
	// XXX base class constructor initializes bonsai - will this interfere with other Tools using bonsai?
	if(derived_classname == "relic"){
		get_ok = InitRelic();
	} else if(derived_classname == "mc"){
		get_ok = InitMC();
	} else if(derived_classname == "t2k"){
		get_ok = InitT2k();
	} else if(derived_classname == "ambe"){
		get_ok = InitAmBe();
	} else if(derived_classname == "merged"){
		get_ok = InitMerged();
	} else {
		Log(m_unique_name+" Unrecognised derived class "+derived_classname,v_error,m_verbose);
		return false;
	}
	
	if(not get_ok){
		Log(m_unique_name+" Failed to Initialise SK2p2MeV class instance!",v_error,m_verbose);
		return false;
	}
	
	// so that downstream Tools can process the outputs without having to go via disk
	// in a separate ToolChain run, we'll also create a TreeReader that allows downstream Tools
	// to access this tree as if it were being read by a TreeReader Tool
	outTreeReader.Load(theOTree);
	m_data->RegisterReader("sk2p2_OutTree", &outTreeReader);
	// tell our TreeReader that it doesn't own the file/tree and shouldn't close them in destructor
	outTreeReader.SetOwnsFile(false);
	
	// set base class options
	// TODO explain these in the README for this tool
	// TODO move all base configuration options here
	//      and have defaults provided via default config files.
	//      for now, only the options used by all current callers are set here.
	
	// initial defaults
	int N10Threshold = 6;
	int N10CutThreshold = 6;
	
	// update with any values given in config file
	m_variables.Get("N10Threshold",N10Threshold);
	m_variables.Get("N10CutThreshold",N10CutThreshold);
	
	// pass to the SK2p2MeV class instance
	ntagger->SetVerbosity(std::min(m_verbose,2));
	ntagger->SetN10Threshold(N10Threshold);
	ntagger->SetN10CutThreshold(N10CutThreshold);
	
	return true;
}

// -----------------------------

bool SK2p2MeV_ntag::Execute(){
	
	get_ok = ntagger->GetBranchValues();
	// important; intialise results to 0
	ntagger->Clear();
	/*
	// ignore this for now - SK2p2MeV currently gets more branches than may be needed,
	// depending on the input source (data / mc). If branches are not present, it just
	// ignores them and continues anyway. TODO better handling.
	if(not get_ok){
		Log(m_unique_name+" Error getting input variables from TreeReader!",v_error,m_verbose);
		return false;
	}
	*/
	
	// get current entry num, and whether it is the last entry in the TChain
	long entry = myTreeReader->GetEntryNumber();
	bool last_entry;
	m_data->vars.Get("StopLoop", last_entry);
	
	// call SK2p2MeV_XXX::Analyze to find the neutrons!
	// results are written to the output TChain
	// TODO write them into the DataModel, so that downstream Tools can access them!
	ntagger->Analyze(entry, last_entry);
	// TODO return something that indicates error - for now SK2p2MeV::Analyze returns void.
	
	return true;
}

// -----------------------------

bool SK2p2MeV_ntag::Finalise(){
	
	// Save results
	fout->Write();
	delete theOTree;
	fout->Close();
	delete fout;
	theOTree = 0;
	fout = 0;
	
	return true;
}

// =============================
// -----------------------------
// =============================

bool SK2p2MeV_ntag::InitRelic(){
	
	// make the tagger
	std::cout<<"Constructing SK2p2MeV_relic"<<std::endl;
	SK2p2MeV_relic* the_ntagger = new SK2p2MeV_relic(geopmt_.xyzpm);
	std::cout<<"Constructed!"<<std::endl;
	
	// make output TTree
	theOTree = the_ntagger->MakeOutputTree(isWIT);
	
	// get additional config variables
	bool fakeData = false;
	int seed = 333;
	int timebin = -1;
	float AFTGate = 1000000.;
	std::string T2KData = "";
	std::string T2KDir = "";
	
	m_variables.Get("AFTGate",  AFTGate);
	m_variables.Get("fakeData", fakeData);   // whether this is fake data(?)
	m_variables.Get("seed",     seed);       // 
	m_variables.Get("timebin",  timebin);    // 
	m_variables.Get("T2KData",  T2KData);    // path to a text file containing T2K run #s
	m_variables.Get("T2KDir",   T2KDir);     // path to the directory containing T2K run data files
	
	if(fakeData && (T2KData=="" || T2KDir=="")){
		Log(m_unique_name+" Error! Using SK2p2MeV_relic with fake data flag requires T2KData and T2KDir to be set!",
		    v_error,m_verbose);
		return false;
	}
	
	the_ntagger->SetAFTGate(AFTGate);
	
	the_ntagger->Initialise(myTreeReader, fakeData, seed, T2KData.c_str(), T2KDir.c_str(), timebin);
	
	ntagger = the_ntagger;
	
	return true;
}

// =============================
// -----------------------------
// =============================

bool SK2p2MeV_ntag::InitMC(){
	
	// make the tagger
	SK2p2MeV_mc* the_ntagger = new SK2p2MeV_mc(geopmt_.xyzpm);
	
	// make output TTree
	theOTree = the_ntagger->MakeOutputTree(isWIT);
	
	// get additional config variables
	int seed = 333;
	bool useTrueVertex=true;                                 // is this a sensible default...?
	bool smearTrueVertex=true;
	float vtxResX=16.f, vtxResY=16.f, vtxResZ=30.f;
	float kshift = -1;
	
	m_variables.Get("seed",       seed);                     // random seed
	m_variables.Get("useTrueVertex",     useTrueVertex);     // whether to use bonsai vtx or MC truth vtx
	m_variables.Get("smearTrueVertex",   smearTrueVertex);   // whether to smear the true MC vertex
	m_variables.Get("vtxResX",           vtxResX);           // if smearing the vertex, how much
	m_variables.Get("vtxResY",           vtxResY);           // to smear it by. stddev of a gaussian
	m_variables.Get("vtxResZ",           vtxResZ);           // random vector added, units cm.
	m_variables.Get("kshift",            kshift);            // a fixed propagation length from MC vtx
	// if kshift>0, a vector of given length in a random direction is added
	// to the MC true vertex to define the vertex for T-ToF calculation.
	
	// There are 17 missing channel in SK-IV
	//mask also some high rates PMT
	const int NMIS    = 17;
	const int MISCH[] = {
		7667,
		7686,
		8937,
		8980,
		9283,
		9312,
		9339,
		9362,
		9415,
		9434,
		9817,
		10685,
		10728,
		11031,
		11060,
		11087,
		11110
	};
	
	the_ntagger->SetMisch (NMIS, MISCH);
	the_ntagger->SetBadch (combad_.nbad, combad_.isqbad);
	the_ntagger->SetVertexResolution(vtxResX, vtxResY, vtxResZ);
	the_ntagger->SetSmearFlag (smearTrueVertex);
	the_ntagger->SetMCFlag (useTrueVertex);
	the_ntagger->SetSeed (seed);
	//the_ntagger->SetDarkCut(false);
	//Linyan
	if (kshift > 0){
		the_ntagger->SetShiftFlag (kTRUE);
		the_ntagger->SetShiftDistance (kshift);
	} else {
		the_ntagger->SetShiftFlag (kFALSE);
	}
	
	get_ok = the_ntagger->Initialise(myTreeReader);
	
	ntagger = the_ntagger;
	
	return get_ok;
}

// =============================
// -----------------------------
// =============================

bool SK2p2MeV_ntag::InitT2k(){
	
	// make the tagger
	SK2p2MeV_t2k* the_ntagger = new SK2p2MeV_t2k(geopmt_.xyzpm);
	
	// make output TTree
	theOTree = the_ntagger->MakeOutputTree(isWIT);
	
	// get additional config variables
	bool random = true;
	Float_t VX = 0.;
	Float_t VY = 0.;
	Float_t VZ = 0.;
	
	m_variables.Get("random", random);   // whether to use random distributed vertices...
	m_variables.Get("VX", VX);           // fixed vertex X
	m_variables.Get("VY", VY);           // fixed vertex Y
	m_variables.Get("VZ", VZ);           // fixed vertex Z
	
	if(!random){
		the_ntagger->SetVertex (VX, VY, VZ);
	}
	
	get_ok = the_ntagger->Initialise(myTreeReader, random);
	
	ntagger = the_ntagger;
	
	return get_ok;
}

// =============================
// -----------------------------
// =============================

bool SK2p2MeV_ntag::InitAmBe(){
	
	// make the tagger
	SK2p2MeV_ambe* the_ntagger = new SK2p2MeV_ambe(geopmt_.xyzpm);
	
	// make output TTree
	theOTree = the_ntagger->MakeOutputTree(isWIT);
	
	// get additional config variables
	float AFTGate = 850000.;
	std::string position="z";   // FIXME what's the third option conventionally called?
	
	m_variables.Get("AFTGate",  AFTGate);
	m_variables.Get("position", position);
	
	// Center
	const Float_t VX_CENTER = 35.3;
	const Float_t VY_CENTER = -70.7;
	const Float_t VZ_CENTER = 0.;
	
	// Z=15m
	const Float_t VX_Z15M = 35.3;
	const Float_t VY_Z15M = -70.7;
	const Float_t VZ_Z15M = 1500.;
	
	// Y=12m
	const Float_t VX_Y12M = 35.3;
	const Float_t VY_Y12M = -1201.9;
	const Float_t VZ_Y12M = 0.;
	
	if      (position == "c") the_ntagger->SetVertex (VX_CENTER, VY_CENTER, VZ_CENTER);
	else if (position == "y") the_ntagger->SetVertex (VX_Y12M, VY_Y12M, VZ_Y12M);
	else                      the_ntagger->SetVertex (VX_Z15M, VY_Z15M, VZ_Z15M);
	
	the_ntagger->SetAFTGate(AFTGate);
	
	get_ok = the_ntagger->Initialise(myTreeReader);
	
	ntagger = the_ntagger;
	
	return get_ok;
}

// =============================
// -----------------------------
// =============================

bool SK2p2MeV_ntag::InitMerged(){
	
	// make the tagger
	std::cout<<"Constructing SK2p2MeV_merged"<<std::endl;
	SK2p2MeV_merged* the_ntagger = new SK2p2MeV_merged(geopmt_.xyzpm);
	std::cout<<"Constructed!"<<std::endl;
	
	// make output TTree
	theOTree = the_ntagger->MakeOutputTree(isWIT);
	
	// get additional config variables
	float AFTGate = 1000000.;
	m_variables.Get("AFTGate",  AFTGate);
	the_ntagger->SetAFTGate(AFTGate);
	
	the_ntagger->Initialise(myTreeReader);
	
	ntagger = the_ntagger;
	
	return true;
}
