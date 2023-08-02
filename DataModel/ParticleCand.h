#ifndef PARTICLECAND_H
#define PARTICLECAND_H

#include "skroot_loweC.h"

struct ParticleCand {
	int EventNumber;
	float EventTime;
	int EntryNumber;
	int PID = 0; //0 = muon 1 = LowE
	std::vector<int> matchedParticleEvNum;
	std::vector<float> matchedParticleTimeDiff;
	std::vector<float> matchedParticleBSEnergy;
	skroot_lowe_common LowECommon;
};

#endif
