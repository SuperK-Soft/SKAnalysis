#ifndef NCaptInfo_NTAG_H
#define NCaptInfo_NTAG_H

#include "NCaptInfo.h"
#include "Algorithms.h"

class NCaptInfo_NTag : public NCaptInfo {
	public:
	NCaptInfo_NTag() : NCaptInfo(){};
	private:
	bool GetCandidates(std::vector<NCaptCandidate>& candidates);
};

#endif
