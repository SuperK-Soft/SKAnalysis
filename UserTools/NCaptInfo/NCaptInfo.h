#ifndef NCaptInfo_H
#define NCaptInfo_H

#include <string>
#include <iostream>
#include <unordered_map>

#include "Tool.h"
#include "NCaptCandidate.h"
#include "MTreeReader.h"

/**
 * \class NCaptInfo
 *
 * This is an abstract base class for a Tool that converts neutron tagging information
 * into a generalised form. 
 * Add support for your own ntagger by making a class that derives from this, and overrides
 * the pure virtual functions InitCandidateReader and GetCandidates. See existing cases for examples.
*
* $Author: M.O $
* $Date: 2023/05/09 $
*/

class TFile;
class TH1D;
class TH2D;

class NCaptInfo: public Tool {
	
	public:
	
	NCaptInfo(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool purpose.
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	virtual bool InitCandidateReader()=0;
	virtual bool GetCandidates(std::vector<NCaptCandidate>& candidates)=0;
	bool MatchToTrueCaptures();
	bool MatchMistags();
	bool PrintCandidates();
	bool MakePlots(int step);
	
	protected:
	MTreeReader* myTreeReader;
	bool match_mistags=false;
	std::string candidates_file="";
	std::string mctruth_file="";
	
	double time_match_tolerance;
	double dist_match_tolerance;
	double likelihood_threshold=-999;
	
	// histograms
	std::string outfilename="";
	TFile* out_file = nullptr;
	TH1D* h_likelihood = nullptr;
	TH1D* h_tdiff = nullptr;
	TH1D* h_xdiff = nullptr;
	TH1D* h_ydiff = nullptr;
	TH1D* h_zdiff = nullptr;
	TH1D* h_rdiff = nullptr;
	TH2D* h_tdiff_vs_metric = nullptr;
	TH2D* h_rdiff_vs_metric = nullptr;
	
	TTree* candtree = nullptr;
	TTree* matchtree = nullptr;
	std::unordered_map<std::string, int> cibranchvars;
	std::unordered_map<std::string, double> cdbranchvars;
	std::unordered_map<std::string, int> mibranchvars;
	std::unordered_map<std::string, double> mdbranchvars;
	
};


#endif
