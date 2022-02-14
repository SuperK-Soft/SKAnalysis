/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef TruthNeutronCaptures_v2_H
#define TruthNeutronCaptures_v2_H

#include <string>
#include <iostream>
#include <vector>
#include <array>

#include "Tool.h"

#include "MTreeReader.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

class TApplication;
class TFile;
class TTree;
class TVector3;
class TLorentzVector;

/**
* \class TruthNeutronCaptures_v2
*
* A simple tool to plot basic characteristics of neutron captures, reading outputs from skdetsim.
* 1. checkout and compile modified SKOFL with secondaries added to MCInfo class.
* 2. compile modified skdetsim with NEUTRON environmental variable defined to enable saving secondaries.
* 3. run skdetsim with option 'SKCNTL-OUTPUTTYPE 1' and 'SKCNTL-FILEFORMAT 1' to generate SKROOT output.
* This file processes the output root files from step 3.
*
* $Author: M.O'Flahery $
* $Date: 2020/08/12 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class TruthNeutronCaptures_v2: public Tool {
	
	public:
	
	TruthNeutronCaptures_v2(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool purpose.
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	std::string toolName;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// file stuff
	std::string inputFile;                      // if just passing a single filename directly to this tool
	std::vector<std::string> input_file_names;  // if using upstream LoadFileList tool
	std::string outputFile;                     // name of output file to write
	int MAX_EVENTS=-1;                          // max n events to process
	int WRITE_FREQUENCY=10;                     // update output file every N fills
	
	// handy constants
	double neutron_mass;
	
	// output file
	TFile* outfile=nullptr;
	TTree* outtree=nullptr;
	
	// Functions
	// =========
	int CalculateVariables();
	int GenerateHistograms();
	int ReadEntry(long entry_number);
	int CreateOutputFile(std::string outputFile);
	void ClearOutputTreeBranches();
	void PrintBranches();
	int WriteTree();
	void CloseFile();
	int DisableUnusedBranches();
	
	// Member variables
	// ================
	int entry_number=0; // input TTree entry, potentially different to event_number, which is a branch variable
	
	// variables to read in
	// ====================
	MTreeReader myTreeReader;                                   // the TTree reader
	const Header   *run_header    = nullptr;
	const MCInfo   *mc_info       = nullptr;
	const SecondaryInfo *sec_info = nullptr;
	
	// variables to write out
	// ======================
	// Each TTree entry will correspond to one event, which may have multiple neutrons
	// Each neutron capture may have multiple gammas, and each gamma may have multiple PMT hits
	
	// file-wise
	std::string out_filename;
	int out_skdetsim_version;
	int out_tba_table_version;
	float out_water_transparency;
	
	// event-wise
	int out_run_number;
	int out_entry_number;   // TTree entry number to be able to identify the source event
	int out_subevent_number;
	
	// primary particle - an event can have many primary particles
	std::vector<int> out_primary_pdg;
	std::vector<double> out_primary_energy;                          // [MeV]
	std::vector<TVector3> out_primary_start_mom;                     // [MeV/c] redundant with above?
	std::vector<TLorentzVector> out_primary_start_pos;               // [cm, ns]
	std::vector<TLorentzVector> out_primary_end_pos;                 // [cm, ns]
	
	// parent nuclide - potentially many per event
	std::vector<int> out_nuclide_parent_pdg;
	std::vector<int> out_nuclide_daughter_pdg;
	
	// neutrons  - potentially many per nuclide? or just one?
	std::vector<TLorentzVector> out_neutron_start_pos; // [cm, ns]
	std::vector<TLorentzVector> out_neutron_end_pos;   // [cm, ns]
	std::vector<double> out_neutron_start_energy;      // [MeV]
	std::vector<double> out_neutron_end_energy;        // [MeV] - on capture/decay
	std::vector<int> out_neutron_end_process;          // geant3 process code
	std::vector<int> out_neutron_ndaughters;
	
	// gammas - potentially many per neutron
	std::vector<std::vector<double> > out_gamma_energy; // [MeV]
	std::vector<std::vector<double> > out_gamma_time;   // [ns] since?
	
	// internal conversion electrons - potentially many per neutron
	std::vector<std::vector<double> > out_electron_energy; // [MeV]
	std::vector<std::vector<double> > out_electron_time;   // [ns] since?
	
	// detector information
	// total charge? time distribution of hits?
	// build a timestamp and calculate time since last event? using PrevT0?
	
};


#endif
