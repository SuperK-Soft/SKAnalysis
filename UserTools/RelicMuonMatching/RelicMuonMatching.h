#ifndef RelicMuonMatching_H
#define RelicMuonMatching_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "ConnectionTable.h"
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
	
	std::string treeReaderName;
	
	MTreeReader* myTreeReader = nullptr;
	
	const Header* myHeader=nullptr;
	const HardwareTrigger* myHardwareTrigger=nullptr;
	const LoweInfo* myLowe=nullptr;
	
	
	int verbosity = 1;
	int m_verbose;
	int v_error = 0;
	
	int lun;
	
	bool RemoveFromDeque(std::vector<int>& particlesToRemove, std::deque<ParticleCand>& particleDeque);
	unsigned long long int bitshiftTime(unsigned long long int t0Time, unsigned long long int hardwareTime);
	bool RelicMuonMatch(std::string particleType, std::deque<ParticleCand>& currentDeque, std::deque<ParticleCand>& targetDeque);
	
	bool WriteRelicInfo();
	bool WriteMuonInfo();
	bool WriteInfo(ParticleCand Event);
	bool Makededx(float (&muentry)[4], float (&mdir)[3], int (&ihcab)[11146], float (&qisk)[11146], float (&tisk)[11146], float (&xyzpm)[11146][3], int &nqiskm, float (&dedx)[200]);
	int mujecttime(float (&v)[4], float (&d)[3], float (&p)[4], double (&dist)[4], float (&cosang)[2]);
	
	bool muonFlag;
	bool relicFlag;
	int currentSubRun;
	
	std::vector<std::string> branchestoSkip;
	
	unsigned long long int muonTime;
	unsigned long long int relicTime;
	float timeDiff;
	unsigned long long int bitOne = 1;
	
	std::vector<int> relicsToRemove;
	std::vector<int> muonsToRemove;
	
	TTree* WriteTree=nullptr;
	TBranch* MatchedEvNumsBranch=nullptr;
	TBranch* MatchedTimeDiffBranch=nullptr;
	TBranch* PIDBranch=nullptr;
	
	std::vector<int> MatchedEvNums;
	std::vector<float> MatchedTimeDiff;
	
	int PID;
	
	int lastRun = 0;
	
	float bffpos[3];
	float hpos[3];
	float bffgood;
	float modd;
	
};


#endif
