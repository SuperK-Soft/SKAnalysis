/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef SK2p2MeV_ntag_H
#define SK2p2MeV_ntag_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.
#include "MTreeReader.h"

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
	MTreeReader* myTreeReader = nullptr;
	MTreeReader outTreeReader{"SK2p2MeV_OutTree"};
	SK2p2MeV* ntagger = nullptr;
	TFile* fout = nullptr;
	TTree* theOTree = nullptr;
	bool isWIT = false;
	
};


#endif
