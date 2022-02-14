/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef FitSpallationDt_H
#define FitSpallationDt_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "basic_array.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

#include "ColourWheel.h"

class TH1;
class TF1;
class MTreeReader;
class MTreeSelection;

/**
* \class FitSpallationDt
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class FitSpallationDt: public Tool {
	
	public:
	FitSpallationDt();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	bool Analyse();
	bool GetBranchValues();
	bool Analyse_Laura();
	bool GetBranchValuesLaura();
	bool GetEnergyCutEfficiencies();
	bool PlotSpallationDt();
	bool FitDtDistribution(TH1& dt_mu_lowe_hist_short, int rangenum);
	// helper functions used in FitSpallationDt
	std::vector<double> MakeLogBins(double xmin, double xmax, int nbins);
	void FixLifetime(TF1& func, std::string isotope);
	void PushFitAmp(TF1& func, std::string isotope);
	void PushFitAmp(double amp, std::string isotope);
	void PullFitAmp(TF1& func, std::string isotope, bool fix=true);
	void PullPaperAmp(TF1& func, std::string isotope, bool threshold_scaling=true, double fixed_scaling=1.);
	double GetPaperAmp(std::string isotope, bool threshold_scaling, double fixed_scaling);
	void BuildPaperPlot();
	TF1 BuildFunction(std::vector<std::string> isotopes, double func_min=0, double func_max=30);
	TF1 BuildFunctionHack(std::vector<std::string> isotopes, double func_min=0, double func_max=30);
	TF1 BuildFunctionNoHack(std::vector<std::string> isotopes, double func_min=0, double func_max=30);
	
	// tool variables
	// ==============
	std::string toolName;
	std::string outputFile="";
	std::string valuesFileMode="";
	std::string valuesFile="";
	double paper_scaling = 1.;
	bool useHack=false;
	std::string efficienciesFile="FlukaBetaEfficiencies.bs";    // BoostStore of efficiencies of energy thresholds
	
	std::string treeReaderName;
	MTreeReader* myTreeReader=nullptr; // the TTree reader
	MTreeSelection* myTreeSelections=nullptr;
	
	int run_min=1;
	int run_max=9999999;
	
	std::vector<float> dt_mu_lowe_vals;  // data to fit
	double livetime=0;
	double binwidth;
	std::string hist_to_fit="log";
	std::string laurasfile="";
	bool fix_const=false;
	bool use_par_limits=false;
	bool split_iso_pairs=false;
	int n_dt_bins=5000;
	int binning_type=0;
	bool random_subtract=false;
	
	// energy threshold comparison
	// ===========================
	std::map<std::string, double> true_effs_6mev;
	std::map<std::string, double> true_effs_8mev;
	std::map<std::string, double> true_effs_scaling;
	std::map<std::string, double> reco_effs_6mev;
	std::map<std::string, double> reco_effs_8mev;
	std::map<std::string, double> reco_effs_scaling;
	
	// results used in fitting of the number of isotope events
	std::map<std::string,double> fit_amps;
	
	ColourWheel colourwheel;
	
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
	// for nominal processing all we need is the array of times between muon and lowe event [seconds]
	basic_array<float*> dt_mu_lowe;
	// when fitting data from laura's file, each mu-lowe pair is a single TTree entry (not an array)
	// and we may apply additional cuts to ensure the same dataset is being fit.
	int muboy_status;
	float dt;
	float lt;
	float energy_relic;
	float bsgood_relic;
	int nrunsk;
	
	// variables to write out
	// ======================
	
};


#endif
