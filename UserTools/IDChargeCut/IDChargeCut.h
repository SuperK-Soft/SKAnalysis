#ifndef IDChargeCut_H
#define IDChargeCut_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
 * \class IDChargeCut
 *
 * Skips the remainder of the toolchain if the 'newMuon' flag is not present, and nqisk is above a given limit.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class IDChargeCut: public Tool {
	
	
	public:
	
	IDChargeCut(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Execute function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	
	private:
	
	int totalHits = 0;
	int hitLimit = 999;
	
	std::string selectorName;
	
};


#endif
