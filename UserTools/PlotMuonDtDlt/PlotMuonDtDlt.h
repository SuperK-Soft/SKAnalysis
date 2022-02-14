/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef PlotMuonDtDlt_H
#define PlotMuonDtDlt_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.
#include "basic_array.h"

class MTreeReader;
class MTreeSelection;
class THStack;

/**
* \class PlotMuonDtDlt
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class PlotMuonDtDlt: public Tool {
	
	public:
	PlotMuonDtDlt();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	bool GetBranchValues();
	bool PlotMuonDt();
	bool PlotMuonDlt();
	bool PlotPaperDt(THStack& ourplots);
	bool PlotPaperDlt(THStack& ourplots);
	//void MeasureDltSystematic(); // TODO
	
	// tool variables
	// ==============
	std::string toolName;
	std::string outputFile="";
	std::string treeReaderName;
	MTreeReader* myTreeReader=nullptr;
	MTreeSelection* myTreeSelections=nullptr;
	
	// maps of muboy class vs a histogram of mu->lowe time and transverse distance
	std::vector<std::vector<float>> dlt_vals_pre{6};
	std::vector<std::vector<float>> dt_vals_pre{6};
	std::vector<std::vector<float>> dlt_vals_post{6};
	std::vector<std::vector<float>> dt_vals_post{6};
	
	// varying muon->lowe dt thresholds for assessing systematic of dlt cut
	int num_dt_cuts=5; // this is the size of the spall_lifetimes vector in PurewaterLi9Rate tool
	// moreover it's the number of cuts named "pre/post_mu_dt_cut_%d" we have.
	// TODO retrieve the list of cuts from the MTreeSelection, count how many we have of this type?
	std::vector<std::vector<float>> dlt_systematic_dt_cuts_pre{5};
	std::vector<std::vector<float>> dlt_systematic_dt_cuts_post{5};
	// each entry is a different dt cut, inner vector is the dlts of passing events
	// for a given dt cut, the difference between values gives the distribution of *spallation* dlt.
	// across the various dt cuts, the difference between spallation dlt distributions gives
	// the systematic error in dlt cut efficiency.
	// We're particularly interested in the variation in the bin corresponding to cut value of dlt=200cm,
	// so make sure we have a bin edge at that value
	
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
	basic_array<int*>   mu_class;                    // muboy muon classification (see enum class at top)
	basic_array<float*> dt_mu_lowe;                  // time between muon and lowe event [seconds]
	basic_array<float*> dlt_mu_lowe;                 // transverse distance between muon and lowe event [cm]
	
	// variables to write out
	// ======================
	
};


#endif
