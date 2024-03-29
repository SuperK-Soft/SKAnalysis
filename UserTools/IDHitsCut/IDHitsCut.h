#ifndef IDHitsCut_H
#define IDHitsCut_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
 * \class IDHitsCut
 *
 * Skips the remainder of the toolchain if the EventType is LowE, and nqisk is above a given limit.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class IDHitsCut: public Tool {
	
	
	public:
	
	IDHitsCut(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Execute function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	
	private:
	
	int totalHits = 0;
	int hitLimit = 999;
	
	std::string selectorName;
	
};


#endif
