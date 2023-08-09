#ifndef WriteMatchedInfo_H
#define WriteMatchedInfo_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "ParticleCand.h"
#include "fortran_routines.h"


/**
 * \class WriteMatchedInfo
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class WriteMatchedInfo: public Tool {


public:
	
	WriteMatchedInfo(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	
private:
	bool WriteInfo(std::vector<ParticleCand>);
	
	int m_verbose;
	
	TTree* WriteTree=nullptr;
	TBranch* MatchedEvNumsBranch=nullptr;
	TBranch* MatchedTimeDiffBranch=nullptr;
	TBranch* PIDBranch=nullptr;
	
	std::vector<int> MatchedEvNums;
	std::vector<float> MatchedTimeDiff;
	
	int PID;
	
};


#endif
