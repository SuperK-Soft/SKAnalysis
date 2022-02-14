/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef FitLi9Lifetime_H
#define FitLi9Lifetime_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "basic_array.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

class MTreeReader;
class MTreeSelection;
class TH1F;

/**
* \class FitLi9Lifetime
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/

class FitLi9Lifetime: public Tool {
	
	public:
	FitLi9Lifetime();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// config variables
	// ================
	float li9_lifetime_dtmin;         // range of dt_mu_lowe values to accept
	float li9_lifetime_dtmax;         // for Li9 candidates, seconds
	std::string outputFile="";
	std::string treeReaderName;
	MTreeReader* myTreeReader=nullptr;
	MTreeSelection* myTreeSelections=nullptr;
	
	// functions
	// =========
	bool GetBranchValues();
	bool PlotLi9BetaEnergy();
	bool PlotLi9LifetimeDt();
	double BinnedLi9DtChi2Fit(TH1F* li9_muon_dt_hist);
	
	// tool variables
	// ==============
	std::string toolName;
	std::vector<float> li9_e_vals;
	std::vector<float> li9_muon_dt_vals;
	
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
	const LoweInfo *LOWE  = new LoweInfo;
	basic_array<float*> dt_mu_lowe;                  // time between muon and lowe event [seconds]
	
	// variables to write out
	// ======================
	
};


#endif
