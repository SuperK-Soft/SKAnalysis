#ifndef RelicMuonMatching_H
#define RelicMuonMatching_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "ParticleCand.h"

/**
 * \class RelicMuonMatching
*/

class RelicMuonMatching: public Tool {
	
	
	public:
	
	RelicMuonMatching(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	
	private:
	
	MTreeReader* rfmReader = nullptr;
	const Header* myHeader=nullptr;
	
	int muWriterLUN;
	int relicWriterLUN;
	
	bool RemoveFromDeque(std::vector<int>& particlesToRemove, std::deque<ParticleCand>& particleDeque);
	bool RelicMuonMatch(bool muonFlag, float currentTime, int subtrg_num=0);
	
	bool WriteRelicInfo();
	bool WriteMuonInfo();
	
	bool muonFlag;
	int currentSubRun;
	int lastRun = 0;
	
	std::vector<int> relicsToRemove;
	std::vector<int> muonsToRemove;
	
	std::vector<int> MatchedEvNums;
	std::vector<float> MatchedTimeDiff;
	
	std::string muSelectorName;
	std::string relicSelectorName;
	
	// defunct....
	//unsigned long long int bitshiftTime(unsigned long long int t0Time, unsigned long long int hardwareTime);
};


#endif
