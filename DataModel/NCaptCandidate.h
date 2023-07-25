#ifndef NCAPTCANDIDATE_H
#define NCAPTCANDIDATE_H
#include "TVector3.h"
#include "BStore.h"
#include "Constants.h"
#include "NCapture.h"

class DataModel;

class NCaptCandidate {
	public:
	NCaptCandidate();
	// XXX if modifying matchType be sure to modify the matchTypes map accordingly!
	enum class matchType {kNotSet=-1,kMistag=0,kHCapture=1,kGdCapture=2,kUnknownCapture=3,kDecaye=4};
	static const std::map<matchType, std::string> matchTypes;
	friend std::ostream& operator << (std::ostream& os, const NCaptCandidate::matchType& obj);
	
	// generated from n-tagger
	std::string algo;
	double likelihood_metric = -1;  // 0->1, 1 is high confidence, 0 is no confidence.
	double capture_time = -1; // ns
	TVector3 capture_pos;  // cm
	double prompt_time = -1; // ns
	TVector3 prompt_pos{999.,999.,999}; // cm
	
	matchType matchtype=matchType::kNotSet;
	BStore featureMap{true,constants::BSTORE_BINARY_FORMAT};  // features used by tagging algorithm
	
	std::string GetMatchType();
	bool SetTrueCaptureIdx(int idx);
	NCapture* GetTrueCapture();
	double* GetCaptTerr();
	double* GetCaptPosErr();
	double* GetPromptTerr();
	double* GetPromptPosErr();
	void Print(bool printFeatureMap=false, bool printFeatureMapValues=false);
	
	private:
	int truecap_idx=-1; // index in NCapturesTrue if matched to a true capture
	bool got_cap_pos_err=false;
	bool got_cap_t_err=false;
	bool got_prompt_pos_err=false;
	bool got_prompt_t_err=false;
	double cap_t_err;
	double cap_pos_err;
	double* cap_t_err_p=nullptr;
	double* cap_pos_err_p=nullptr;
	double prompt_t_err;
	double prompt_pos_err;
	double* prompt_t_err_p=nullptr;
	double* prompt_pos_err_p=nullptr;
	DataModel* m_data=nullptr;
	
};

#endif
