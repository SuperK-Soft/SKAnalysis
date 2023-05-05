#ifndef ntag_BDT_H
#define ntag_BDT_H

#ifndef PYTHON
// dummy class
#include "DummyTool.h"
class ntag_BDT : public DummyTool{ public: ntag_BDT(){
#error "ntag_BDT tool requires python and bind!"
}; };
#else

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.
#include "DataModel.h"
#include "Algorithms.h"

#include "TFile.h"
#include "TTree.h"

/*#include <cstdlib>*/
/*#include <vector>*/
/*#include <map>*/

/*#if not defined(__CINT__) || defined(__MAKECINT__)*/
/*#include "TMVA/Tools.h"*/
/*#include "TMVA/Reader.h"*/
/*#include "TMVA/MethodCuts.h"*/
/*#endif*/

namespace py = pybind11;
using namespace pybind11::literals;

/**
 * \class ntag_BDT
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class ntag_BDT : public Tool {
	
	public:
	
	ntag_BDT(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	bool GetBranchValues();
	Int_t GetNlowIndex(Float_t rsqred, Float_t z, const Int_t init);
	
	MTreeReader* myTreeReader = nullptr;
	TFile* outfile = nullptr;
	TTree* treeout = nullptr;
	int WRITE_FREQUENCY = 500;
	
	// BDT configuration variables
	int N10TH;
	int NLOWINDEX;
	
	// BDT model
	py::object predict_proba5;
	
	// variables read from input file
	const Header *HEADER = nullptr;
	const LoweInfo *LOWE = nullptr;
	basic_array<float*> smearedvertex;
	basic_array<float*> trmsold;
	basic_array<float*> phi;
	basic_array<float*> theta;
	basic_array<float*> dthetarms;
	basic_array<float*> dqrms;
	basic_array<float*> dqmean;
	basic_array<float*> trmsdiff;
	basic_array<float*> mintrms6;
	basic_array<float*> mintrms3;
	basic_array<float*> bswall;
	basic_array<float*> bse;
	basic_array<float*> fpdist;
	basic_array<float*> bfdist;
	basic_array<float*> fwall;
	basic_array<float*> dt;
	basic_array<int*> n10;
	int np;
	int type;
	int nnhits;
	int N200M;
	int T200M;
	basic_array<int*> nc;
	basic_array<int*> nnback;
	basic_array<int*> n300;
	basic_array<int*> nnhighq;
	basic_array<int*> nnlowtheta;
	basic_array<int*> n10d;
	basic_array<float*> vx;
	basic_array<float*> vy;
	basic_array<float*> vz;
	// vector for branches 'Nlow1','Nlow2','Nlow3'... unknown number of such branches
	std::vector<basic_array<int*>> Nlow;
	
	// output variables - arrays of size defined by branch 'np'
	int MAX_EVENTS=500;
	float* neutron5 = nullptr;
	int* nlow = nullptr;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	std::string toolName;
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
};


#endif // no python
#endif
