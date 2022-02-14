/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef LoadBetaSpectraFluka_H
#define LoadBetaSpectraFluka_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"

#include "TH1.h"

/**
* \class LoadBetaSpectraFluka
*
* Read output from FLUKA simulation describing beta energy spectra of spallation isotopes
* Also includes reconstructed energy following propagation in skdetsim & reconstruction with bonsai
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class LoadBetaSpectraFluka: public Tool {
	
	public:
	LoadBetaSpectraFluka();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	bool Analyse();   ///< Perform analysis of this entry
	
	private:
	// functions
	// =========
	int ReadEntry(long entry_number);
	int GetBranches();
	int DisableUnusedBranches();
	
	// tool variables
	// ==============
	std::string toolName;
	std::string inputFile;
	std::string histosFile="beta_spectra.root";
	std::string mapsFile="beta_spectra.bs";
	int maxEvents;
	int entrynum=0;
	double bonsai_goodness_cut=0.4;
	double bonsai_maxE_cut=1000;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// variables to read in
	// ====================
	MTreeReader myTreeReader; // the TTree reader
	// input branch variables
//	int muonID;
//	int event_num;
	int Z;
	int A;
//	double decay_time;
//	int npart;
//	double true_decay_pos[3];
	double true_beta_E;
	double true_photon_E;
//	double bonsai_pos[3];
	double bonsai_energy;
	double bonsai_goodness;
	
	std::map<std::string,int> unknown_isotopes;
	
	// variables to write out
	// ======================
	std::map<std::string, TH1D> true_spectra;
	std::map<std::string, TH1D> reco_spectra;
	std::map<std::string, TH1D> reco_over_true_spectra;
	
	// how many events had a true energy above/below thresholds
	// to determine detection efficiency
	// use a map so we get a clean cut at threshold - could use a histo counts provided
	// we had a bin edge at the threshold but this is easier....probably
	std::map<std::string,int> true_events_below_8MeV;
	std::map<std::string,int> true_events_above_8MeV;
	std::map<std::string,int> true_events_below_6MeV;
	std::map<std::string,int> true_events_above_6MeV;
	
	// same but using energy from bonsai
	std::map<std::string,int> reco_events_below_8MeV;
	std::map<std::string,int> reco_events_above_8MeV;
	std::map<std::string,int> reco_events_below_6MeV;
	std::map<std::string,int> reco_events_above_6MeV;
	
};


#endif
