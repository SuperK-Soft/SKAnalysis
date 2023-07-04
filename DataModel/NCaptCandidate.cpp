#include "NCaptCandidate.h"
#include "DataModel.h"

NCaptCandidate::NCaptCandidate(){
	m_data = DataModel::GetInstance();
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
	got_prompt_pos_err=false;
	got_prompt_t_err=false;
	return true;
}

NCapture* NCaptCandidate::GetTrueCapture(){
	if(truecap_idx<0 || truecap_idx >= m_data->NCapturesTrue.size()){
		return nullptr;
	}
	return &m_data->NCapturesTrue.at(truecap_idx);
}

double* NCaptCandidate::GetCaptTerr(){
	if(got_cap_t_err) return cap_t_err_p;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetTime()==nullptr) return nullptr;
	cap_t_err = (capture_time - *truecap->GetTime());
	cap_t_err_p = &cap_t_err;
	return cap_t_err_p;
}

double* NCaptCandidate::GetCaptPosErr(){
	if(got_cap_pos_err) return cap_pos_err_p;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetPos()==nullptr) return nullptr;
	cap_pos_err = (capture_pos - *truecap->GetPos()).Mag();
	cap_pos_err_p = &cap_pos_err;
	return cap_pos_err_p;
}

double* NCaptCandidate::GetPromptTerr(){
	if(got_prompt_t_err) return prompt_t_err_p;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetNeutron()==nullptr ||
	   truecap->GetNeutron()->GetStartTime()==nullptr) return nullptr;
	prompt_t_err = (prompt_time - *truecap->GetNeutron()->GetStartTime());
	prompt_t_err_p = &prompt_t_err;
	return prompt_t_err_p;
}

double* NCaptCandidate::GetPromptPosErr(){
	if(got_prompt_pos_err) return prompt_pos_err_p;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetNeutron()==nullptr ||
	   truecap->GetNeutron()->GetStartPos()==nullptr) return nullptr;
	prompt_pos_err = (prompt_pos - *truecap->GetNeutron()->GetStartPos()).Mag();
	prompt_pos_err_p = &prompt_pos_err;
	return prompt_pos_err_p;
}

void NCaptCandidate::Print(bool printFeatureMap, bool printFeatureMapValues){
	std::string cap_t_err_s="?";
	std::string cap_pos_err_s="?";
	if(GetCaptTerr()!=nullptr) cap_t_err_s = toString(*cap_t_err_p/1000.);
	if(GetCaptPosErr()!=nullptr) cap_pos_err_s = toString(*cap_pos_err_p);
	
	std::string prompt_t_err_s="?";
	std::string prompt_pos_err_s="?";
	if(GetPromptTerr()!=nullptr) prompt_t_err_s = toString(*prompt_t_err_p/1000.);
	if(GetPromptPosErr()!=nullptr) prompt_pos_err_s = toString(*prompt_pos_err_p);
	
	std::cout<<"\tlikelihood: "<<likelihood_metric<<"\n"
	         <<"\tprompt time [us]: "<<prompt_time/1000.<<"\n"
	         <<"\tprompt position [cm]: "<<toString(prompt_pos)<<"\n"
	         <<"\tcapture time [us]: "<<capture_time/1000.<<"\n"
	         <<"\tcapture position [cm]: "<<toString(capture_pos)<<"\n"
	         <<"\tmatched true capture index: "<<truecap_idx<<"\n"
	         <<"\treconstructed prompt time error [us]: "<<prompt_t_err_s<<"\n"
	         <<"\treconstructed prompt position error [cm]: "<<prompt_pos_err_s<<"\n"
	         <<"\treconstructed capture time error [us]: "<<cap_t_err_s<<"\n"
	         <<"\treconstructed capture position error [cm]: "<<cap_pos_err_s<<"\n"
	         <<"\tmatch type: "<<GetMatchType()<<"\n";
	if(printFeatureMap){
		std::cout<<"\tlist of feature variables: \n";
		if(printFeatureMapValues){
			featureMap.Print();
		} else {
			std::vector<std::string> features = featureMap.GetKeys();
			for(const std::string& featurename : features){
				std::cout<<"\t\t"<<featurename<<"\n";
			}
		}
	}
}
