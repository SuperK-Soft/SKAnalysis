/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef MYTOOLRoot_H
#define MYTOOLRoot_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"

/**
* \class MyToolRoot
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class MyToolRoot: public Tool {
	
	public:
	MyToolRoot();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	int GetBranchValues();
	
	// tool variables
	// ==============
	
	// variables to read in
	// ====================
	MTreeReader* myTreeReader=nullptr; // the TTree reader
	
	// variables to write out
	// ======================
	
};


#endif
