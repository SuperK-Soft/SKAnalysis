#ifndef CombinedFitter_H
#define CombinedFitter_H

#include <string>
#include <iostream>
#include <vector>

#include "Tool.h"

#include "pmt_geometry.h"
#include "bonsaifit.h"
#include "pairlikelihood.h"

#include <TH1.h>

/**
 * \class CombinedFitter
 *
 * This is a vertex fitter class to run the bonsai fit either locally
 * or using bonsai in skofl
 *
 * $Author: L.Kneale $
 * Contact: e.kneale@sheffield.ac.uk
 */
class CombinedFitter: public Tool {


	public:

	CombinedFitter(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Execute function used to perform Tool purpose. 
	bool Finalise(); ///< Finalise function used to clean up resources.



	private:

	// tool functions
	// Functions to get the prompt and delayed hits
	void SetPromptHits(std::vector<float>& charges, std::vector<float>& times, std::vector<int>& cableIDs, int &nhit);
	int SetAftHits(float prompt_trigger_t0, int SLE_threshold, int AddNoise, int numPMTs, float darkmc, float last_prompt, std::vector<float> chargesRaw, std::vector<float> timesRaw, std::vector<int> cableIDsRaw, int nhitsRaw, std::vector<float>& charges_AFT, std::vector<float>& times_AFT, std::vector<int>& cableIDs_AFT);

	// Fitting functions
	void lbfset0(void *dat, int *numdat); // Resets the fitting variables
	int CalculateNX(int timewindow, float* vertex, std::vector<int> cableIDs, std::vector<float> times, std::vector<int>& cableIDsNX); // Calculates the energy analogue (e.g. bn50)
	void SingleEventFit(std::vector<float> charges, std::vector<float> times, std::vector<int> cableIDs, int nhit,float *bsvertex,float *bsresult,float *bsgood);
	
	// Functions to fill the trees
	void FillLoweTree(); // Fills the LOWE tree in the input skroot file
	void FillCoreTree(); // Fills the ntuples in the output file

	// tool variables
	// ==============
	//	std::string treeReaderName;
	//	MTreeReader* myTreeReader=nullptr; 
	// skroot tree
	TTree* t;
	TFile* fout = nullptr;
	// file and tree for ntuple output
	TTree *outputTree =nullptr;
	// variables to fill and write out
	// fitting
	float x, y, z = 0;
	float x_prev, y_prev, z_prev = 0;
	float x_combined, y_combined, z_combined = 0;
	float tgood, tgood_prev = 0;
	float tgood_combined, tgood_combined_prev = 0;
	// mc
	float mcx, mcy, mcz = 0;
	float mcx_prev, mcy_prev, mcz_prev = 0;
	int pdg, pdg_prev = 0;
	float mcx_ncapture, mcy_ncapture, mcz_ncapture = 0;
	float mc_energy, mc_energy_prev = 0;	
	float mct_ncapture=0;
	int nhitsAFT_raw=0;
	int bsn50 = 0;		int bsn10 = 0;
	int bsn50AFT = 0;	int bsn10AFT = 0;
	int n50_prev = 0;	int n10_prev = 0;
	int n50_aft = 0;	int n10_aft = 0;

	// variables to use in the tool
	struct HitInfo
	{
		float charge;
		float time;
		int cableID;
	};
	int SLE_threshold = 0;
	int nread=0;          // just track num loops for printing
	int nrunsk_last=0;    // to know when to read in new transparency data at start of each new run
	int nsubsk_last=0;    // same for new badch list, loaded on new run and subrun
	float darkmc=0;
	float watert;         // water transparency
	int numPMTs;          // total number of PMTs


	// variables to read in
	// ====================
	std::string readerName="";
	std::string outputFile = "";
	int lun=0;
	int nhitcut;
	int ev=0;
	bool MC=false;
	bool addNoise = true;
	int dataSrc=0;     // 0=sktqz_ common block, 1=TQReal branch
	int bonsaiSrc = 0; // 0= built-in bonsai calls; 1 = direct bonsai functions

	// bonsai
	// ======
	pmt_geometry* bsgeom;
	pairlikelihood* bspairlike;
	bonsaifit *bspairfit;


};


#endif
