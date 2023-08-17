/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef lfallfit_simple_H
#define lfallfit_simple_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

/**
* \class lfallfit_simple
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class lfallfit_simple: public Tool {
	
	public:
	lfallfit_simple();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	
	// tool variables
	// ==============
	std::string fname_in;
	std::string fname_out;
	
	int lun = 10;         // a unique ID to give the TreeManager for accessing this file.
	int nread=0;          // just track num loops for printing
	int nrunsk_last=0;    // to know when to read in new transparency data at start of each new run
	float watert;         // water transparency
	
};


#endif
