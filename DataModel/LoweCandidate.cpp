#include "LoweCandidate.h"
#include "DataModel.h"

LoweCandidate::LoweCandidate(){
	m_data = DataModel::GetInstance();
	featureMap=new BStore{true,constants::BSTORE_BINARY_FORMAT};
	recoVars=new BStore{true,constants::BSTORE_BINARY_FORMAT};
}

LoweCandidate::~LoweCandidate(){
	if(featureMap) delete featureMap;
	if(recoVars) delete recoVars;
	featureMap=nullptr;
	recoVars=nullptr;
}

LoweCandidate::LoweCandidate(LoweCandidate&& rhs){
	
	m_data = rhs.m_data;
	algo = rhs.algo;
	goodness_metric = rhs.goodness_metric;
	event_time = rhs.event_time;
	event_pos = rhs.event_pos;
	event_energy = rhs.event_energy;
	matchtype = rhs.matchtype;
	
	trueparticle_idx = rhs.trueparticle_idx;
	got_pos_err = rhs.got_pos_err;
	got_t_err = rhs.got_t_err;
	got_e_err = rhs.got_e_err;
	t_err = rhs.t_err;
	pos_err = rhs.pos_err;
	e_err = rhs.e_err;
	
	featureMap = rhs.featureMap;
	rhs.featureMap = nullptr;
	recoVars = rhs.recoVars;
	rhs.recoVars = nullptr;
	
}

// XXX if modifying this, be sure to also modify the enum class in the header file!
//enum class matchType {kNotSet=-1,kMistag=0,kHCapture=1,kGdCapture=2,kUnknownCapture=3,kDecaye=4};
const std::map<LoweCandidate::matchType, std::string> LoweCandidate::matchTypes {
        {LoweCandidate::matchType::kNotSet,"Unmatched"},
        {LoweCandidate::matchType::kMistag,"Mistag"},  // noise, flasher, afterpulse...
        {LoweCandidate::matchType::kIBDPositron,"IBD Positron"},
        {LoweCandidate::matchType::kSpallation,"Spallation"},
        {LoweCandidate::matchType::kDecaye,"Decay-e"},
        {LoweCandidate::matchType::kOther,"Other"} // solar (e-)? reactor (IBD)? NCQE (gamma?)?
};

std::string LoweCandidate::GetMatchType(){
	return (matchTypes.count(matchtype)) ? matchTypes.at(matchtype) : std::to_string((int)(matchtype));
}

bool LoweCandidate::SetTrueParticleIdx(int idx){
	if(idx<0 || idx>= m_data->eventParticles.size()){
		return false;
	}
	trueparticle_idx = idx;
	// reset these when true match changes
	got_pos_err=false;
	got_t_err=false;
	got_e_err=false;
	return true;
}

MParticle* LoweCandidate::GetTrueParticle(){
	if(trueparticle_idx<0 || trueparticle_idx >= m_data->eventParticles.size()){
		return nullptr;
	}
	return &m_data->eventParticles.at(trueparticle_idx);
}

double* LoweCandidate::GetTerr(){
	if(got_t_err) return &t_err;
	MParticle* trueparticle = GetTrueParticle();
	if(trueparticle==nullptr || trueparticle->GetStartTime()==nullptr) return nullptr;
	t_err = (event_time - *trueparticle->GetStartTime());
	got_t_err=true;
	return &t_err;
}

double* LoweCandidate::GetPosErr(){
	if(got_pos_err) return &pos_err;
	MParticle* trueparticle = GetTrueParticle();
	if(trueparticle==nullptr || trueparticle->GetStartPos()==nullptr) return nullptr;
	pos_err = (event_pos - *trueparticle->GetStartPos()).Mag();
	got_pos_err=true;
	return &pos_err;
}

double* LoweCandidate::GetEnergyErr(){
	if(got_e_err) return &e_err;
	MParticle* trueparticle = GetTrueParticle();
	if(trueparticle==nullptr || trueparticle->GetStartE()==nullptr) return nullptr;
	e_err = (event_energy - *trueparticle->GetStartE());
	got_e_err=true;
	return &e_err;
}

void LoweCandidate::Print(bool printRecoVarsMap, bool printRecoVarsMapValues, bool printFeatureMap, bool printFeatureMapValues){
	std::string t_err_s="?";
	std::string pos_err_s="?";
	if(GetTerr()!=nullptr) t_err_s = toString(t_err);
	if(GetPosErr()!=nullptr) pos_err_s = toString(pos_err);
	
	std::cout<<"\tgoodness: "<<goodness_metric<<"\n"
	         <<"\time [ns]: "<<event_time<<"\n"
	         <<"\tposition [cm]: "<<toString(event_pos)<<"\n"
	         <<"\tmatched true particle index: "<<trueparticle_idx<<"\n"
	         <<"\tmatch type: "<<GetMatchType()<<"\n"
	         <<"\treconstructed time error [us]: "<<t_err_s<<"\n"
	         <<"\treconstructed position error [cm]: "<<pos_err_s<<"\n";
	
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
}

std::ostream& operator << (std::ostream& os, const LoweCandidate::matchType& obj){
	///os << static_cast<std::underlying_type<A>::type>(obj);
	if(LoweCandidate::matchTypes.count(obj)) os << LoweCandidate::matchTypes.at(obj);
	else os << "?";
	return os;
}
