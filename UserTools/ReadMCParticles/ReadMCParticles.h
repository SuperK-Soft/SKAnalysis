/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef ReadMCParticles_H
#define ReadMCParticles_H

#include <string>
#include <iostream>

#include "Tool.h"

/**
* \class ReadMCParticles
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class ReadMCParticles: public Tool {
	
	public:
	ReadMCParticles();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	
	bool GetSecondaryInfo();
	bool GetSecondaryVectors();
	bool GetAtmpdInfo();
	bool PrintEvent();
	bool PrintSecondaryInfo();
	bool PrintSecondaryVectors(bool checkconsistency=false);
	MTreeReader* myTreeReader=nullptr;
	int dataSrc=0; // 0= SecondaryInfo atmpd arrays, 1=SecondaryInfo vectors
	int debugEntryNum=-1; // crank up verbosity for this event
	
	const SecondaryInfo * sec_info = nullptr;
	const MCInfo* mc_info = nullptr;
	
};


#endif
