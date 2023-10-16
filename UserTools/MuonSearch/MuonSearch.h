#ifndef MuonSearch_H
#define MuonSearch_H

#include <string>
#include <iostream>
#include <bitset>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "ConnectionTable.h"


/**
* \class MuonSearch
*
* Applies the software trigger to search for coincident HE+OD triggers.
* If found, marks the event as containing one or muons.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class MuonSearch: public Tool {
	
	
	public:
	
	MuonSearch(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param 	data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	
	double coincidence_threshold=100;
	std::string selectorName;
	EventType eventType;
	
};


#endif
