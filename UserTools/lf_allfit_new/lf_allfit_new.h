/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef lf_allfit_new_H
#define lf_allfit_new_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

#include "basic_array.h"

/**
* \class lf_allfit_new
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class lf_allfit_new: public Tool {
	
	public:
	lf_allfit_new();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	
	// tool variables
	// ==============
	std::string toolName;
	
	int nread=0;          // just track num loops for printing
	int nrunsk_last=0;    // to know when to read in new transparency data at start of each new run
	int nsubsk_last=0;    // same for new badch list, loaded on new run and subrun
	float watert;         // water transparency
	
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
	std::string readerName="";
	int lun=0;
	bool MC=false;
	
	// variables to write out
	// ======================
	
};


#endif
