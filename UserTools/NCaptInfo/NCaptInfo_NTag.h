#ifndef NCaptInfo_NTAG_H
#define NCaptInfo_NTAG_H

#include "NCaptInfo.h"
#include "Algorithms.h"

class NCaptInfo_NTag : public NCaptInfo {
	public:
	NCaptInfo_NTag() : NCaptInfo(){};
	private:
	bool InitCandidateReader();
	bool GetCandidates(std::vector<NCaptCandidate>& candidates);
	MTreeReader* candidatesTreeReader=nullptr;
	MTreeReader* variablesTreeReader=nullptr;
};

#endif
