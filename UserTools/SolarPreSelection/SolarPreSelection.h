#ifndef SolarPreSelection_H
#define SolarPreSelection_H

#include <string>
#include <iostream>
#include <set>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "SolarRelic.h"

/**
 * \class SolarPreSelection
*/

class SolarPreSelection: public Tool {
	
	public:
	
	SolarPreSelection(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Execute function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	void FillRelic(SolarRelic& a_relic);
	
	private:
	
	MTreeReader* rfmReader = nullptr;
	std::string solarSelectorName;
	TTree* solarTree=nullptr;
	
	std::deque<SolarRelic> relics_this_run;
	TFile* fout=nullptr;
	TTree* relic_tree=nullptr;
	SolarRelic out_relic;
	
	Long64_t thiseventticks=0;
	double match_window = 60; // [seconds]
	int64_t match_window_ticks;
	int solar_nqisk_precut_thresh; // nhits
	
	std::vector<int32_t> matched_ev_nums;
	std::vector<uint64_t> matched_out_entry_nums;
	std::vector<size_t> matched_indices; // for forwarding to WriteSolarMatches tool
	std::vector<double> matched_tdiffs;  // seconds
	int solarcount=0;
	int matchedsolarcount=0;
	int erased_count=0;
	
};


#endif
