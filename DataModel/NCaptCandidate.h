#ifndef NCAPTCANDIDATE_H
#define NCAPTCANDIDATE_H
#include "TVector3.h"
#include "BStore.h"
#include "Constants.h"
#include "NCapture.h"

class DataModel;
class LoweCandidate;

class NCaptCandidate {
	public:
	NCaptCandidate();
	~NCaptCandidate();
	NCaptCandidate(NCaptCandidate&) = delete; // no copy construction until BStore supports it
	NCaptCandidate(NCaptCandidate&&);         // we need to handle move construction manually
	
	// XXX if modifying matchType be sure to modify the matchTypes map accordingly!
	enum class matchType {kNotSet=-1,kMistag=0,kHCapture=1,kGdCapture=2,kUnknownCapture=3,kDecaye=4};
	static const std::map<matchType, std::string> matchTypes;
	friend std::ostream& operator << (std::ostream& os, const NCaptCandidate::matchType& obj);
	
	// generated from n-tagger
	std::string algo;
	double capture_likelihood_metric = -1;  // 0->1, 1 is high confidence, 0 is no confidence.
	double capture_time = -1; // ns
	TVector3 capture_pos;  // cm
	double capture_E = -1; // MeV
	
	matchType matchtype=matchType::kNotSet;
	BStore* featureMap=nullptr; //{true,constants::BSTORE_BINARY_FORMAT};  // features used by tagging algorithm
	BStore* recoVars=nullptr; //{true,constants::BSTORE_BINARY_FORMAT};    // additional reconstruction variables
	
	std::string GetMatchType();
	bool SetTrueCaptureIdx(int idx);
	bool SetPromptEvent(int idx);
	NCapture* GetTrueCapture();
	LoweCandidate* GetPromptEvent();
	double* GetCaptTerr();
	double* GetCaptPosErr();
	double* GetCaptEnergyErr();
	void Print(bool printPrompt=true, bool printRecoVarsMap=false, bool printRecoVarsMapValues=false, bool printFeatureMap=false, bool printFeatureMapValues=false);
	
	private:
	int truecap_idx=-1; // index in NCapturesTrue if matched to a true capture
	int prompt_idx=-1;  // index in LoweCandidates of associated prompt event
	bool got_cap_pos_err=false;
	bool got_cap_t_err=false;
	bool got_cap_e_err=false;
	double cap_t_err;
	double cap_pos_err;
	double cap_e_err;
	DataModel* m_data=nullptr;
	
};

#endif
