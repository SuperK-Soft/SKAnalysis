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
	
	int rfmReaderLUN;
	int muWriterLUN;
	int relicWriterLUN;
	
	bool RemoveFromDeque(std::vector<int>& particlesToRemove, std::deque<ParticleCand>& particleDeque);
	bool RelicMuonMatch(bool loweEventFlag, int64_t currentTime, int subtrg_num=0, int32_t it0xsk=0);
	
	bool WriteRelicInfo();
	bool WriteMuonInfo();
	
	EventType eventType;
	int currentSubRun;
	int lastRun = 0;
	int last_nevsk=0;
	double match_window = 60; // [seconds]
	int64_t match_window_ticks;
	
	std::vector<int> relicsToRemove;
	std::vector<int> muonsToRemove;
	
	std::vector<int> MatchedEvNums;
	std::vector<float> MatchedTimeDiff;
	std::vector<float> MatchedParticleE;
	
	int nextmuentry=0;
	int nextrelicentry=0;
	
	std::string muSelectorName;
	std::string relicSelectorName;
	std::string rfmReaderName;
	
	// track how many muon/relic pairs we check
	uint64_t tdiffcount=0;
	// and how many pass our +-60s tdiff cut.
	uint64_t passing_tdiffcount=0;
	
	// defunct....
	//unsigned long long int bitshiftTime(unsigned long long int t0Time, unsigned long long int hardwareTime);
};


#endif
