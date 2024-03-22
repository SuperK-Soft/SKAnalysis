#include "NCaptCandidate.h"
#include "DataModel.h"
#include "LoweCandidate.h"

NCaptCandidate::NCaptCandidate(){
	m_data = DataModel::GetInstance();
	featureMap = new BStore{true,constants::BSTORE_BINARY_FORMAT};
	recoVars = new BStore{true,constants::BSTORE_BINARY_FORMAT};
}

NCaptCandidate::~NCaptCandidate(){
	if(featureMap) delete featureMap;
	if(recoVars) delete recoVars;
	featureMap = nullptr;
	recoVars = nullptr;
}

NCaptCandidate::NCaptCandidate(NCaptCandidate&& rhs){
	
	m_data= rhs.m_data;
	
	algo = rhs.algo;
	capture_likelihood_metric = rhs.capture_likelihood_metric;
	capture_time = rhs.capture_time;
	capture_pos = rhs.capture_pos;
	capture_E = rhs.capture_E;
	matchtype = rhs.matchtype;
	
	truecap_idx = rhs.truecap_idx;
	prompt_idx = rhs.prompt_idx;
	
	cap_t_err = rhs.cap_t_err;
	cap_pos_err = rhs.cap_pos_err;
	cap_e_err = rhs.cap_e_err;
	got_cap_pos_err = rhs.got_cap_pos_err;
	got_cap_t_err = rhs.got_cap_t_err;
	got_cap_e_err = rhs.got_cap_e_err;
	
	// here's the reason we need to manually specify the move constructor:
	// we need to nullify the pointers to any on-heap data members
	// from the moved-from object after the move.
	featureMap = rhs.featureMap;
	rhs.featureMap = nullptr;
	recoVars = rhs.recoVars;
	rhs.recoVars = nullptr;
	
}

// XXX if modifying this, be sure to also modify the enum class in the header file!
//enum class matchType {kNotSet=-1,kMistag=0,kHCapture=1,kGdCapture=2,kUnknownCapture=3,kDecaye=4};
const std::map<NCaptCandidate::matchType, std::string> NCaptCandidate::matchTypes {
        {NCaptCandidate::matchType::kNotSet,"Unmatched"},
        {NCaptCandidate::matchType::kMistag,"Mistag"},
        {NCaptCandidate::matchType::kHCapture,"H-capture"},
        {NCaptCandidate::matchType::kGdCapture,"Gd-capture"},
        {NCaptCandidate::matchType::kUnknownCapture,"?-capture"},  // matched to true capture on unknown isotope
        {NCaptCandidate::matchType::kDecaye,"Decay-e"}
};

std::string NCaptCandidate::GetMatchType(){
	return (matchTypes.count(matchtype)) ? matchTypes.at(matchtype) : std::to_string((int)(matchtype));
}

bool NCaptCandidate::SetTrueCaptureIdx(int idx){
	if(idx<0 || idx>= m_data->NCapturesTrue.size()){
		return false;
	}
	truecap_idx = idx;
	// reset these when true match changes
	got_cap_pos_err=false;
	got_cap_t_err=false;
	got_cap_e_err=false;
	return true;
}

bool NCaptCandidate::SetPromptEvent(int idx){
	if(idx<0 || idx >= m_data->LoweCandidates.size()){
		return false;
	}
	prompt_idx = idx;
	return true;
}

NCapture* NCaptCandidate::GetTrueCapture(){
	if(truecap_idx<0 || truecap_idx >= m_data->NCapturesTrue.size()){
		return nullptr;
	}
	return &m_data->NCapturesTrue.at(truecap_idx);
}

LoweCandidate* NCaptCandidate::GetPromptEvent(){
	if(prompt_idx<0 || prompt_idx >= m_data->LoweCandidates.size()){
		return nullptr;
	}
	return &m_data->LoweCandidates.at(prompt_idx);
}

double* NCaptCandidate::GetCaptTerr(){
	if(got_cap_t_err) return &cap_t_err;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetTime()==nullptr) return nullptr;
	cap_t_err = (capture_time - *truecap->GetTime());
	got_cap_t_err=true;
	return &cap_t_err;
}

double* NCaptCandidate::GetCaptPosErr(){
	if(got_cap_pos_err) return &cap_pos_err;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetPos()==nullptr) return nullptr;
	cap_pos_err = (capture_pos - *truecap->GetPos()).Mag();
	got_cap_pos_err=true;
	return &cap_pos_err;
}

double* NCaptCandidate::GetCaptEnergyErr(){
	if(got_cap_e_err) return &cap_e_err;
	NCapture* truecap = GetTrueCapture();
	double truecap_E;
	if(truecap==nullptr || truecap->SumGammaE(truecap_E)==false) return nullptr;
	cap_e_err = (capture_E - truecap_E);
	return &cap_e_err;
}

void NCaptCandidate::Print(bool printPrompt, bool printRecoVarsMap, bool printRecoVarsMapValues, bool printFeatureMap, bool printFeatureMapValues){
	std::string cap_t_err_s="?";
	std::string cap_pos_err_s="?";
	if(GetCaptTerr()!=nullptr) cap_t_err_s = toString(cap_t_err/1000.);
	if(GetCaptPosErr()!=nullptr) cap_pos_err_s = toString(cap_pos_err);
	
	std::cout<<"\tcapture likelihood: "<<capture_likelihood_metric<<"\n"
	         <<"\tcapture time [us]: "<<capture_time/1000.<<"\n"
	         <<"\tcapture position [cm]: "<<toString(capture_pos)<<"\n"
	         <<"\tmatched true capture index: "<<truecap_idx<<"\n"
	         <<"\tmatch type: "<<GetMatchType()<<"\n"
	         <<"\treconstructed capture time error [us]: "<<cap_t_err_s<<"\n"
	         <<"\treconstructed capture position error [cm]: "<<cap_pos_err_s<<"\n";
	
	if(printRecoVarsMap){
		std::cout<<"\tlist of reconstruction variables: \n";
		if(printRecoVarsMapValues && recoVars){
			recoVars->Print();
		} else if(recoVars){
			std::vector<std::string> recovariables = recoVars->GetKeys();
			for(const std::string& recovar : recovariables){
				std::cout<<"\t\t"<<recovar<<"\n";
			}
		}
	}
	
	if(printFeatureMap){
		std::cout<<"\tlist of feature variables: \n";
		if(printFeatureMapValues && featureMap){
			featureMap->Print();
		} else if(featureMap){
			std::vector<std::string> features = featureMap->GetKeys();
			for(const std::string& featurename : features){
				std::cout<<"\t\t"<<featurename<<"\n";
			}
		}
	}
	
	if(printPrompt && GetPromptEvent()!=nullptr) GetPromptEvent()->Print();
	
}

std::ostream& operator << (std::ostream& os, const NCaptCandidate::matchType& obj){
	///os << static_cast<std::underlying_type<A>::type>(obj);
	if(NCaptCandidate::matchTypes.count(obj)) os << NCaptCandidate::matchTypes.at(obj);
	else os << "?";
	return os;
}
