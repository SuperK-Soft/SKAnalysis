#ifndef ODCut_H
#define ODCut_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "ConnectionTable.h"


/**
 * \class ODCut
 *
 * Remove events without the 'newMuon' flag set, which have either the OD trigger bit set within the primary trigger window, or more than hitsThreshold OD hits within 500-1300us following the primary trigger.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class ODCut: public Tool {
	
	public:
	
	ODCut(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Execute function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	int hitsThreshold=20;
	float windowMinT = 500.;
	float windowMaxT = 1300.;
	int Nskipped = 0;
	std::string selectorName;
	
};


#endif
