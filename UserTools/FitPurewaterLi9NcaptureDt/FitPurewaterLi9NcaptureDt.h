/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef FitPurewaterLi9NcaptureDt_H
#define FitPurewaterLi9NcaptureDt_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.
#include "basic_array.h"

class TH1F;
class MTreeReader;
class MTreeSelection;

/**
* \class FitPurewaterLi9NcaptureDt
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class FitPurewaterLi9NcaptureDt: public Tool {
	
	public:
	FitPurewaterLi9NcaptureDt();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	bool GetBranchValues();
	bool PlotNcaptureDt();
	double BinnedNcapDtChi2Fit(TH1F* li9_ncap_dt_hist);
	bool UnbinnedNcapDtLogLikeFit(TH1F* li9_ncap_dt_hist, double num_li9_events);
	double ncap_lifetime_loglike(double* x, double* par);
	
	// tool variables
	// ==============
	std::string toolName;
	std::vector<float> li9_ntag_dt_vals;
	std::string outputFile="";
	float ncap_dtmin;                 // range of dt_mu_ncap values to accept
	float ncap_dtmax;                 // for Li9 abundance extraction, **microseconds**
	std::string treeReaderName;
	MTreeReader* myTreeReader=nullptr;
	MTreeSelection* myTreeSelections=nullptr;
	
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
	int num_neutron_candidates;                      // num pulses: i.e. num neutrons in AFT window
	basic_array<float*> dt_lowe_n;                   // dt positron->neutron
	
	// variables to write out
	// ======================
	
};


#endif
