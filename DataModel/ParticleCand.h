#ifndef PARTICLECAND_H
#define PARTICLECAND_H

#include "skroot_loweC.h"

struct ParticleCand {
	int EventNumber;
	int SubTriggerNumber;
	int nevhwsk;
	int it0xsk;
	int64_t EventTicks;
	int NumRollovers;
	int InEntryNumber;  // entry number in input file TTree
	int OutEntryNumber; // entry number in output file TTree
	int PID = 0; //0 = muon 1 = LowE
	std::vector<int> matchedParticleEvNum;         // nevsk of matches
	std::vector<int> matchedParticleInEntryNum;    // TTree entry number of matches in input tree
	std::vector<int> matchedParticleOutEntryNum;   // TTree entry number of matches in output tree
	std::vector<bool> matchedParticleHasAFT;
	std::vector<float> matchedParticleTimeDiff;
	std::vector<float> matchedParticleBSEnergy;
	skroot_lowe_common LowECommon;
	bool hasAFT;
};

#endif
