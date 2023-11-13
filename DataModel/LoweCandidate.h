#ifndef LOWECANDIDATE_H
#define LOWECANDIDATE_H
#include "TVector3.h"
#include "BStore.h"
#include "Constants.h"
#include "MParticle.h"

class DataModel;

class LoweCandidate {
	public:
	LoweCandidate();
	~LoweCandidate();
	LoweCandidate(LoweCandidate&) = delete;
	LoweCandidate(LoweCandidate&&);
	// XXX if modifying matchType be sure to modify the matchTypes map accordingly!
	enum class matchType {kNotSet=-1,kMistag=0,kIBDPositron=1,kSpallation=2,kDecaye=3,kOther=4};
	static const std::map<matchType, std::string> matchTypes;
	friend std::ostream& operator << (std::ostream& os, const LoweCandidate::matchType& obj);
	
	// generated from lowe reco algorithm (bonsai, clusfit...)
	std::string algo;
	double goodness_metric = -1;  // usually bsgood[1], for example
	double event_time = 9999; // ns
	TVector3 event_pos{9999,9999,9999};  // cm
	double event_energy = 9999; // MeV
	
	matchType matchtype=matchType::kNotSet;
	BStore* featureMap=nullptr; //{true,constants::BSTORE_BINARY_FORMAT};  // features used by tagging algorithm
	BStore* recoVars=nullptr; //{true,constants::BSTORE_BINARY_FORMAT};    // additional reconstruction variables
	
	std::string GetMatchType();
	bool SetTrueParticleIdx(int idx);
	MParticle* GetTrueParticle();
	double* GetTerr();
	double* GetPosErr();
	double* GetEnergyErr();
	void Print(bool printRecoVarsMap=false, bool printRecoVarsMapValues=false, bool printFeatureMap=false, bool printFeatureMapValues=false);
	
	private:
	int trueparticle_idx=-1; // index in MParticles if matched to a true lowe particle
	bool got_pos_err=false;
	bool got_t_err=false;
	bool got_e_err=false;
	double t_err;
	double pos_err;
	double e_err;
	DataModel* m_data=nullptr;
	
};

#endif
