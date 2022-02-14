#ifndef PurewaterSpallAbundanceCuts_H
#define PurewaterSpallAbundanceCuts_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "Algorithms.h"
#include "MTreeReader.h"
#include "MTreeSelection.h"
#include "SkrootHeaders.h"    // MCInfo, Header etc.
#include "thirdredvars.h"     // ThirdRed class

/**
* \class PurewaterSpallAbundanceCuts
*
* Measure rate of Li9 production by extracting the number of Li9 events following muons
* Methods based on 2015 paper by Yang Zhang
*
* $Author: M.O'Flaherty $
* $Date: 2020/12/11 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/

class PurewaterSpallAbundanceCuts: public Tool {
	
	public:
	
	PurewaterSpallAbundanceCuts(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose.
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	
	// file IO stuff
	// =============
	std::string treeReaderName="";                              // name of MTreeReader for input
	MTreeReader* myTreeReader;                                  // the TTree reader
	MTreeSelection myTreeSelections;                            // record what passes what cuts
	int entry_number=0;                                         // input TTree entry
	std::string outputFile="li9_cuts.root";                     // output file to write
	
	// cut configurations
	// ==================
	int run_min=0;
	int run_max=0;
	float max_closest_muon_dt=0.001;  // 1ms
	float max_closest_lowe_dx=490;    // 490cm
	float ntag_FOM_threshold=0.95;    // BDT FOM threshold for Li9 ntagging
	
	// functions
	// =========
	bool GetBranchValues();           // retrieve tree branches
	bool Analyse();                   // main body
	bool apply_third_reduction(const ThirdRed *th, const LoweInfo *LOWE);
	
	// variables to read in
	// ====================
	const Header  *HEADER = new Header;
	const LoweInfo *LOWE  = new LoweInfo;
	const ThirdRed *thirdredvars = new ThirdRed;
	int num_neutron_candidates;                      // num pulses: i.e. num neutrons in AFT window
	int max_hits_200ns_AFT;                          // max num hits in 200ns sliding window within AFT trigger
	basic_array<float*> ntag_FOM;                    // ntag BDT output figure of merit
	basic_array<float*> dt_lowe_n;                   // dt positron->neutron
	int num_pre_muons;                               // num muons in 30s preceding lowe event
	int num_post_muons;                              // num muons in 30s following lowe event
//	basic_array<float*> spaloglike;                  // log likelihood of being spallation based on...?
//	basic_array<float*> spaloglike_shfld;            // version based on...?
//	basic_array<float*> spaloglike_kirk;             // version based on...?
//	basic_array<float*> spaloglike_kirk_shfld;       // version based on...?
//	basic_array<float*> spaloglike_li9;              // version based on...?
//	basic_array<float*> spaloglike_kirk_li9;         // version based on...?
	basic_array<int*>   mu_class;                    // muboy muon classification (see enum class at top)
//	basic_array<int*>   mubntrack;                   // num muons found by muboy
	basic_array<int*>   mu_index;                    // index of this muon, of those found by muboy
	basic_array<float*> mu_fit_goodness;             // muboy goodness of fit: >0.4 is ok for single thru muons
//	basic_array<float*> mubffgood;                   // brute force fitter goodness of fit
	basic_array<float*> dt_mu_lowe;                  // time between muon and lowe event [seconds]
//	basic_array<float*> spadt_li9;                   // version based on...?
	basic_array<float*> dlt_mu_lowe;                 // transverse distance between muon and lowe event [cm]
//	basic_array<float*> spadll;                      // longitudinal distance between muon and lowe event [cm]
//	basic_array<float*> spadll_kirk;                 // version based on...?
//	basic_array<float*> sparesq;                     // surplus charge over a MIP with reco track length
//	basic_array<float*> spaqpeak;                    // ? 
//	basic_array<float*> spaqpeak_kirk;               // version based on...?
//	basic_array<float*> spamuqismsk;                 // ?
//	basic_array<float*> spamuqismsk_pertrack;        // ?
//	basic_array<float*> spadts;                      // ?
//	basic_array<int*>   candidates;                  // ?
//	basic_array<int*>   muindex;                     // ?
//	basic_array<int*>   neut_flag;                   // ?
//	basic_array<int*>   mult_flag;                   // ?
//	basic_array<float*[3]> neut_shift;               // ?
//	basic_array<float*> neutdiff;                    // ?
	basic_array<float*> closest_lowe_60s;            // closest distance to another lowe event within 60s??
	
	// crude livetime tracker
	// ======================
	struct tm runstart = {0};
	struct tm runend = {0};
	int current_run=0;
	double livetime=0;
	void IncrementLivetime();
	void AddLastRunTime();
	
	// standard tool stuff
	// ===================
	std::string toolName;
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
