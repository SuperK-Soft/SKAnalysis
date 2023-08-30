/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef TreeReaderDemo_H
#define TreeReaderDemo_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

/**
* \class TreeReaderDemo
* A simple tool to demonstrate processing events in a ROOT file with the TreeReader tool.
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class TreeReaderDemo: public Tool {
	
	public:
	TreeReaderDemo();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	bool GetBranchValues();
	
	// tool variables
	// ==============
	std::string treeReaderName;
	MTreeReader* myTreeReader=nullptr;
	
	int event_num=0;
	
	// variables to read in
	// ====================
	const Header* header=nullptr;
	const LoweInfo* loweinfo=nullptr;
	const TQReal* tqreal=nullptr;
	SecondaryInfo* secondaries=nullptr;
	
	// variables to write out
	// ======================
	
};


#endif
