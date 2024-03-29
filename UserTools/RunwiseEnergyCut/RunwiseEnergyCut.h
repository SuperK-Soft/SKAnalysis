#ifndef RunwiseEnergyCut_H
#define RunwiseEnergyCut_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include "Tool.h"

/**
* \class RunwiseEnergyCut
*
* Apply a cut on reconstructed (bonsai) energy for events marked as EventType::LowE, where the
* cut threshold may vary depending on the run number.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class RunwiseEnergyCut: public Tool {
	
	public:
	
	RunwiseEnergyCut(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	bool ParseOptions(std::string& configfile);
	
	private:
	std::string selectorName;
	int Nskipped = 0;
	std::vector<std::pair<std::pair<int,int>,std::pair<float,float>>> cuts;
};


#endif
