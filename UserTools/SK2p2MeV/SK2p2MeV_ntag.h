/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef SK2p2MeV_ntag_H
#define SK2p2MeV_ntag_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

/**
* \class SK2p2MeV_ntag
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/

class SK2p2MeV;

class SK2p2MeV_ntag: public Tool {
	
	public:
	SK2p2MeV_ntag();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	bool InitRelic();
	bool InitMC();
	bool InitT2k();
	bool InitAmBe();
	
	// tool variables
	// ==============
	std::string toolName;
	MTreeReader* myTreeReader = nullptr;
	SK2p2MeV* ntagger = nullptr;
	TFile* fout = nullptr;
	TTree* theOTree = nullptr;
	bool isWIT = false;
	
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
	
	// variables to write out
	// ======================
	
};


#endif
