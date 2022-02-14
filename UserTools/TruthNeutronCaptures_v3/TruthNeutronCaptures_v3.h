/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef TruthNeutronCaptures_v3_H
#define TruthNeutronCaptures_v3_H

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
* \class TruthNeutronCaptures_v3
*
* A simple tool to plot basic characteristics of neutron captures in SKG4, reading outputs from Harada-san's
* SecondaryInfo arrays. These represent largely the same variable set as from the atmpd information
* available via fillnt_simple and h2root. There are, however, minor differences accounted for in this tool
* compared to the TruthNeutronCaptures tool used for skdetsim->fillnt_simple->h2root files.
*
* $Author: M.O'Flahery $
* $Date: 2020/08/12 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class TruthNeutronCaptures_v3: public Tool {
	
	public:
	
	TruthNeutronCaptures_v3(); ///< Simple constructor
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
	
	TFile* outfile=nullptr;
	TTree* outtree=nullptr;
	
	// Functions
	// =========
	void CopyVariables();
	int CalculateVariables();
	int GenerateHistograms();
	int ReadEntryNtuple(long entry_number);
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
	
	const MCInfo* mc_info = nullptr;
	const SecondaryInfo* sec_info = nullptr;
	
	// run meta info
	//int simulation version??
	double water_transparency;                                  // [cm]
	
	// event meta info
	int run_number;
	int subrun_number;
	int event_number;
	int subevent_number;
	
	// event level detector info
	int N_hit_ID_PMTs;                                          // "nqisk"
	int total_ID_pes;                                           // "qismsk"
	int max_ID_PMT_pes;                                         // "qimxsk", max # PEs from a single ID PMT?
	
	// primary event
	basic_array<float> primary_event_vertex;                    // [cm]
	float primary_event_dist_from_wall;                         // [cm]
	int n_outgoing_primaries;                                   // should be (tot_n_primaries - 2)...
	
	// following are arrays of size n_outgoing_primaries
	basic_array<int*> primary_PDG_code;                         // MCInfo stores int pdg codes in ipvc
	basic_array<float(*)[3]> primary_start_mom;                 // [units?] this and ipv are arrays of size npar
	
	// secondaries - first secondaries arrays...
	int n_secondaries_1;
	// the following are arrays of size npar2
	basic_array<int*> secondary_G3_code_1;                      // 
	basic_array<float(*)[3]> secondary_start_vertex_1;          // array of 3, [cm?] what about time?
	basic_array<float*> secondary_start_dist_from_wall_1;       // [cm?]
	basic_array<float(*)[3]> secondary_start_mom_1;             // [units?]
	basic_array<int*> secondary_origin_1;                       // what is "origin"?
	
	// secondaries - second secondaries array...
	int n_secondaries_2;
	// the following are arrays of size nscndprt
	basic_array<int*> secondary_PDG_code_2;                     //
	basic_array<float(*)[3]> secondary_start_vertex_2;          // [units?]
	basic_array<float*> secondary_start_time_2;                 // [ns]? relative to event start?
	basic_array<float(*)[3]> secondary_start_mom_2;             // [units?]
	basic_array<int*> secondary_gen_process;                    // use constants::G3_process_code_to_string
	basic_array<int*> secondary_n_daughters;                    // 
	basic_array<int*> secondary_first_daugher_index;            // if >0, 1-based index in this array
	basic_array<int*> parent_index;                             // if >0, 1-based index in this array
	
	// further parentage information - Useful?
//	basic_array<int*> parent_G3_code;                           // or is it a PDG code?
	basic_array<float(*)[3]> parent_mom_at_sec_creation;        // use w/daughter γ to see n energy @ capture
	basic_array<float(*)[3]> parent_init_pos;                   // [cm?] position of parent @ birth
	basic_array<float(*)[3]> parent_init_mom;                   // [MeV?] momentum of parent @ birth
	basic_array<int*> parent_trackid;                           // maybe primary parent index???
	
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
//	int out_subrun_number;
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
