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
	got_pos_err=false;
	got_t_err=false;
	return true;
}

NCapture* NCaptCandidate::GetTrueCapture(){
	if(truecap_idx<0 || truecap_idx >= m_data->NCapturesTrue.size()){
		return nullptr;
	}
	return &m_data->NCapturesTrue.at(truecap_idx);
}

double* NCaptCandidate::GetTerr(){
	if(got_t_err) return t_err_p;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetTime()==nullptr) return nullptr;
	t_err = (capture_time - *truecap->GetTime());
	t_err_p = &t_err;
	return t_err_p;
}

double* NCaptCandidate::GetPosErr(){
	if(got_pos_err) return pos_err_p;
	NCapture* truecap = GetTrueCapture();
	if(truecap==nullptr || truecap->GetPos()==nullptr) return nullptr;
	pos_err = (*truecap->GetPos() - capture_pos).Mag();
	pos_err_p = &pos_err;
	return pos_err_p;
}

void NCaptCandidate::Print(bool printFeatureMap, bool printFeatureMapValues){
	std::string t_err_s="?";
	std::string pos_err_s="?";
	if(GetTerr()!=nullptr) t_err_s = toString(*t_err_p/1000.);
	if(GetPosErr()!=nullptr) pos_err_s = toString(*pos_err_p);
	
	std::cout<<"\tlikelihood: "<<likelihood_metric<<"\n"
	         <<"\ttime [us]: "<<capture_time/1000.<<"\n"
	         <<"\tposition [cm]: "<<toString(capture_pos)<<"\n"
	         <<"\tmatched true capture index: "<<truecap_idx<<"\n"
	         <<"\treconstructed time error [us]: "<<t_err_s<<"\n"
	         <<"\treconstructed position error [cm]: "<<pos_err_s<<"\n"
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
