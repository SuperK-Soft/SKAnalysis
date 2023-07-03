#ifndef NCaptInfo_H
#define NCaptInfo_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "NCaptCandidate.h"
#include "MTreeReader.h"

/**
 * \class NCaptInfo
 *
 * This is an abstract base class for a Tool that converts neutron tagging information
 * into a generalised form. 
 * Add support for your own ntagger by making a class that derives from this, and overrides
 * the pure virtual function GetCandidates. See existing cases for examples.
*
* $Author: M.O $
* $Date: 2023/05/09 $
*/

class TH1D;
class TH2D;

class NCaptInfo: public Tool {
	
	public:
	
	NCaptInfo(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	virtual bool GetCandidates(std::vector<NCaptCandidate>& candidates)=0;
	bool MatchToTrueCaptures();
	bool MatchMistags();
	bool PrintCandidates();
	bool MakePlots(int step);
	
	protected:
	MTreeReader* myTreeReader;
	std::string m_unique_name=""; // TODO remove once we upgrade ToolFrameworkCore
	bool match_mistags=false;
	
	double time_match_tolerance;
	double dist_match_tolerance;
	double likelihood_threshold;
	
	// histograms
	TH1D* h_likelihood = nullptr;
	TH1D* h_tdiff = nullptr;
	TH1D* h_xdiff = nullptr;
	TH1D* h_ydiff = nullptr;
	TH1D* h_zdiff = nullptr;
	TH1D* h_rdiff = nullptr;
	TH2D* h_tdiff_vs_metric = nullptr;
	TH2D* h_rdiff_vs_metric = nullptr;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
};


#endif