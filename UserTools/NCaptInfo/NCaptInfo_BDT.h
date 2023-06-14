#ifndef NCaptInfo_BDT_H
#define NCaptInfo_BDT_H

#include "NCaptInfo.h"
#include "Algorithms.h"

class NCaptInfo_BDT : public NCaptInfo {
	public:
	NCaptInfo_BDT() : NCaptInfo(){};
	private:
	bool GetCandidates(std::vector<NCaptCandidate>& candidates);
};

#endif
