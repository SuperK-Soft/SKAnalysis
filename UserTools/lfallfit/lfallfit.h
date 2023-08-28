/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef lfallfit_H
#define lfallfit_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

#include "basic_array.h"

/**
* \class lfallfit
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class lfallfit: public Tool {
	
	public:
	lfallfit();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	
	// tool variables
	// ==============
	
	int nread=0;                     // just track num loops for printing
	int nrunsk_last=0;               // to know when to read in new transparency data at start of each new run
	int nsubsk_last=0;               // same for new badch list, loaded on new run and subrun
	float watert;                    // water transparency
	int reference_watert_run=85609;  // reference run for water properties if MC
	
	// variables to read in
	// ====================
	std::string readerName="";
	int lun=0;
	bool MC=false;
	int writeout=0;
	bool delete_outside_hits=false;
	
	// variables to write out
	// ======================
	
};


#endif
