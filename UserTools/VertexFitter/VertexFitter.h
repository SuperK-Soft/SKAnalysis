#ifndef VertexFitter_H
#define VertexFitter_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "pmt_geometry.h"
#include "bonsaifit.h"
#include "likelihood.h"
#include "goodness.h"
#include "fourhitgrid.h"

#include <TH1.h>
/**
 * \class VertexFitter
 *
 * This is a vertex fitter class to run the bonsai fit either locally
 * or using bonsai in skofl
 *
 * $Author: L.Kneale $
 * Contact: e.kneale@sheffield.ac.uk
 */
class VertexFitter: public Tool {


 public:

  VertexFitter(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Executre function used to perform Tool perpose. 
  bool Finalise(); ///< Finalise funciton used to clean up resorces.



 private:

	// tool functions
  	void lbfset0(void *dat, int *numdat);
	int CalculateNX(int timewindow, float* vertex, int cableIDs[], float times[], int (&cableIDs_twindow)[500]);
    
	// tool variables
    // ==============
//	std::string treeReaderName;
//	MTreeReader* myTreeReader=nullptr; 
	TTree* t;
    int nread=0;          // just track num loops for printing
    int nrunsk_last=0;    // to know when to read in new transparency data at start of each new run
    int nsubsk_last=0;    // same for new badch list, loaded on new run and subrun
    float watert;         // water transparency
	int numPMTs;          // total number of PMTs

	float prev_t_nsec=0;
	TH1D *ht = new TH1D();
    // variables to read in
    // ====================
    std::string readerName="";
    int lun=0;
	int ev=0;
    bool MC=false;
	int dataSrc=0;     // 0=sktqz_ common block, 1=TQReal branch
	int bonsaiSrc = 0; // 0= built-in bonsai calls; 1 = direct bonsai functions

	// bonsai
	// ======
    pmt_geometry* bsgeom;
	likelihood* bslike;
	bonsaifit *bsfit;

};


#endif
