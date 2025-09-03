/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef lfallfit_H
#define lfallfit_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

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
	int nread=0;  // track num loops for printing
	std::string readerName="";
	int lun=0;
	bool MC=false;
	bool update_lowe_branch=true;
	int max_nqisk_for_clusfit;
	int flag_skip=0;   // what reconstruction steps to do (or what to skip). Default (0): everything.
	bool checkEventType=false; // whether to only apply when EventType is LowE
};


#endif
