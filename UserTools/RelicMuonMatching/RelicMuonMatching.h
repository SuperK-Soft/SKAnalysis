#ifndef RelicMuonMatching_H
#define RelicMuonMatching_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "ParticleCand.h"
#include "HistogramBuilder.h"

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
	
	std::string muSelectorName;
	std::string relicSelectorName;
	MTreeReader* rfmReader = nullptr;
	
	bool RemoveFromDeque(std::vector<int>& particlesToRemove, std::deque<ParticleCand>& particleDeque);
	bool RelicMuonMatch(bool loweEventFlag, int64_t currentTicks, int subtrg_num=0, int32_t it0xsk=0);
	
	EventType eventType;
	int currentSubRun;
	int lastRun = 0;
	int last_nevsk=0;
	int num_rollovers=0;
	double match_window = 60; // [seconds]
	int64_t match_window_ticks;
	
	std::vector<int> relicsToRemove;
	std::vector<int> muonsToRemove;
	
	int32_t lastnevhwsk, lastit0sk, last_rollover_nevsk;
	int64_t lasteventticks, lastmuticks, lastrelicticks;
	
	int nextmuentry=0;
	int nextrelicentry=0;
	
	uint64_t muoncount=0;
	uint64_t reliccount=0;
	// track how many muon/relic pairs we check
	uint64_t tdiffcount=0;
	// and how many pass our +-60s tdiff cut.
	uint64_t passing_tdiffcount=0;
	
	HistogramBuilder hb;
	std::string distros_file;
	
	// defunct....
	//unsigned long long int bitshiftTime(unsigned long long int t0Time, unsigned long long int hardwareTime);
};


#endif
