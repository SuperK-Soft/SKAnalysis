/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef lf_allfit_H
#define lf_allfit_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

#include "basic_array.h"

class MTreeReader;
class MTreeSelection;

/**
* \class lf_allfit
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class lf_allfit: public Tool {
	
	public:
	lf_allfit();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	
	// tool variables
	// ==============
	std::string toolName;
	// for reading input data
	//std::string treeReaderName;
	//MTreeReader* myTreeReader=nullptr;
	//MTreeSelection* myTreeSelections=nullptr;
	//bool GetBranchValues();
	std::string fname_in;
	std::string fname_out;
	
	int lun = 10;         // a unique ID to give the TreeManager for accessing this file.
	int nread=0;          // just track num loops for printing
	int nrunsk_last=0;    // to know when to read in new transparency data at start of each new run
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
	
	// variables to write out
	// ======================
	
};


#endif
