/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef MakeNCaptTree_H
#define MakeNCaptTree_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"

class TFile;
class TTree;
#include "TVector3.h"
class TLorentzVector;

/**
* \class MakeNCaptTree
*
* This tool reads files generated by the TruthNeutronCaptures tool and makes an auxiliary file with a Tree
* that adds some derived variables... could probably be made part of TruthNeutronCaptures?
*
* $Author: M.O'Flahery $
* $Date: 2020/08/12 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class MakeNCaptTree: public Tool {
	
	public:
	MakeNCaptTree(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();        ///< Executre function used to perform Tool perpose. 
	bool Finalise();       ///< Finalise funciton used to clean up resorces.
	
	private:
	// functions
	// =========
	int GetBranches();
	int MakeHistos();
	int FillFriend();
	void ClearOutputTreeBranches();
	int WriteTree();
	
	// tool variables
	// ==============
	std::string toolName;
	std::string outputFile; // or just add to the input file?
	int maxEvents=-1;
	int WRITE_FREQUENCY;
	
	int entrynum=0;
	std::map<std::string,int> capture_nuclide_vs_count;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// variables to write out
	// ======================
	TFile* outfile=nullptr;
	TTree* friendTree=nullptr;
	TVector3 neutrino_momentum;
	TVector3 muon_momentum;
	std::vector<double> neutron_longitudinal_travel;     // relative to neutrino dir
	std::vector<double> neutron_perpendicular_travel;    // relative to neutrino dir
	std::vector<double> total_gamma_energy;              // decay gammas
	std::vector<double> total_electron_energy;           // decay electrons from internal conversion
	std::vector<double> total_daughter_energy;           // sum of above
	std::vector<int> neutron_n_gammas;                   // count, per capture
	std::vector<int> neutron_n_electrons;                // count, per capture
	std::vector<int> neutron_n_daughters;                // count, per capture
	
	TVector3* neutrino_momentump=&neutrino_momentum;     // as far as i can tell we ought not to need these
	TVector3* muon_momentump=&muon_momentum;             // but TTree::Branch("name",TVector3* obj) segfaults???
	
	// variables to read in
	// ====================
	MTreeReader* myTreeReader=nullptr; // the TTree reader
	TTree* intree=nullptr;
	
	// file-wise
	std::string filename;
	float water_transparency;
	int skdetsim_version;
	int tba_table_version;
	
	// event-wise
	int entry_number;     // TTree entry number to be able to identify the source event
	int subevent_number;  // 
	
	// primary particle - an event can have many primary particles
	const std::vector<int>* primary_pdg=nullptr;                                 // 
	const std::vector<double>* primary_energy=nullptr;                           // [MeV]
	const std::vector<TLorentzVector>* primary_start_pos=nullptr;                // [cm, ns]
	const std::vector<TLorentzVector>* primary_end_pos=nullptr;                  // [cm, ns]
	const std::vector<TVector3>* primary_start_mom=nullptr;                      // [MeV/c]
	
	// parent nuclide - potentially many per event
	const std::vector<int>* nuclide_pdg=nullptr;
	const std::vector<TLorentzVector>* nuclide_creation_pos=nullptr;              // [cm, ns]
	const std::vector<TLorentzVector>* nuclide_decay_pos=nullptr;                 // [cm, ns]
	const std::vector<int>* nuclide_daughter_pdg=nullptr;                         // 
	
	// neutrons  - potentially many per nuclide? or just one?
	const std::vector<TLorentzVector>* neutron_start_pos=nullptr;   // [cm, ns]
	const std::vector<TLorentzVector>* neutron_end_pos=nullptr;     // [cm, ns]
	const std::vector<double>* neutron_start_energy=nullptr;        // [MeV]
	const std::vector<double>* neutron_end_energy=nullptr;          // [MeV] - on capture/decay
	const std::vector<int>* neutron_end_process=nullptr;            // geant3 process code
	const std::vector<int>* neutron_ndaughters=nullptr;              // num gammas
	
	// gammas - potentially many per neutron
	const std::vector<std::vector<double> >* gamma_energy=nullptr;  // [MeV]
	const std::vector<std::vector<double> >* gamma_time=nullptr;    // [ns]
	
	// electrons - potentially many per neutron
	const std::vector<std::vector<double> >* electron_energy=nullptr;  // [MeV]
	const std::vector<std::vector<double> >* electron_time=nullptr;    // [ns]
	
	// detector information
	// total charge? time distribution of hits?
	// build a timestamp and calculate time since last event? using PrevT0?
	
};


#endif