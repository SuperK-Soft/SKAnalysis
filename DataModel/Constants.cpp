/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#include "Constants.h"
#include "DataModel.h"
#include <bitset>

std::set<std::string> fundamental_types{
	"bool",
	"char",
	"char16_t",
	"char32_t",
	"wchar_t",
	"signed char",
	"short int",
	"int",
	"long int",
	"long long int",
	"unsigned char",
	"unsigned short int",
	"unsigned int",
	"unsigned long int",
	"unsigned long long int",
	"float",
	"double",
	"long double",
//	"std::nullptr_t",
//	"decltype(nullptr)",
	"void*"
};

std::set<std::string> container_types{
	"pair"        // not strictly a container
	"array",
	"vector",
	"deque",
	"forward_list",
	"list",
	"map",
	"multimap",
	"queue",
	"priority_queue",
	"set",
	"multiset",
	"stack",
	"unordered_map",
	"unordered_set",
	"unordered_multimap",
	"unordered_multiset",
	"valarray",   // not strictly a container
	"bitset"      // not strictly a container
};

std::map<TInterpreter::EErrorCode, std::string> TInterpreterErrors{
	{TInterpreter::kNoError,"kNoError"},
	{TInterpreter::kRecoverable,"kRecoverable"},
	{TInterpreter::kDangerous,"kDangerous"},
	{TInterpreter::kFatal,"kFatal"},
	{TInterpreter::kProcessing,"kProcessing"}
};

std::string G3_process_code_to_string(int process_code){
	if(constants::G3_process_code_to_string.count(process_code)){
		return constants::G3_process_code_to_string.at(process_code);
	}
	return "?";
}

std::string G4_process_code_to_string(int process_code){
	if(constants::G4_process_code_to_string.count(process_code)){
		return constants::G4_process_code_to_string.at(process_code);
	}
	return "?";
}

std::string numnu_code_to_string(int numnu_code){
	if(constants::numnu_code_to_string.count(numnu_code)){
		return constants::numnu_code_to_string.at(numnu_code);
	} else if(numnu_code>5){
		return "fs_other";
	}
	return "?";
}

std::string neut_mode_to_string(int neut_code){
	if(constants::neut_mode_to_string.count(neut_code)){
		return constants::neut_mode_to_string.at(neut_code);
	}
	return "?";
}

std::string PdgToString(int code){
	const TDatabasePDG* pdgdb = TDatabasePDG::Instance();
	if(pdgdb!=nullptr){
		const TParticlePDG* particle = pdgdb->GetParticle(code);
		if(particle!=nullptr) return particle->GetName();
	}
	if(constants::pdg_to_string->count(code)!=0){
		return constants::pdg_to_string->at(code);
	}
	return std::to_string(code);
}

int StringToPdg(std::string name){
	const TDatabasePDG* pdgdb = TDatabasePDG::Instance();
	if(pdgdb!=nullptr){
		const TParticlePDG* particle = pdgdb->GetParticle(name.c_str());
		if(particle!=nullptr) return particle->PdgCode();
	}
	// we have some additional custom ones
	if(constants::string_to_pdg.count(name)!=0){
		return constants::string_to_pdg.at(name);
	}
	return -1;
}

std::string G3ParticleCodeToString(int code){
	if(constants::g3_particle_code_to_string.count(code)){
		return constants::g3_particle_code_to_string.at(code);
	}
	return std::to_string(code);
}

int StringToG3ParticleCode(std::string name){
	if(constants::string_to_g3_particle_code.count(name)){
		return constants::string_to_g3_particle_code.at(name);
	}
	return -1;
}

int G3ParticleCodeToPdg(int code){
	if(constants::g3_particle_code_to_pdg.count(code)){
		return constants::g3_particle_code_to_pdg.at(code);
	}
	return -1;
}

int PdgToG3ParticleCode(int code){
	if(constants::pdg_to_g3_particle_code.count(code)){
		return constants::pdg_to_g3_particle_code.at(code);
	}
	return -1;
}

double PdgToMass(int code){
	const TDatabasePDG* pdgdb = TDatabasePDG::Instance();
	if(pdgdb!=nullptr){
		const TParticlePDG* particle = pdgdb->GetParticle(code);
		if(particle!=nullptr) return particle->Mass()*1000.;      // converted to MeV
	}
	// nuclei aren't in the particle database. To first order we can assume
	// the mass of the nucleus is the sum of masses of nucleons
	// pdg codes for nucleons are 10-digit numbers ±10LZZZAAAI, giving us Z and A
	std::string pdgasstring = std::to_string(code);
	if(pdgasstring.length()==10){
		int nprotons = std::stoi(pdgasstring.substr(3,3));
		int nnucleons = std::stoi(pdgasstring.substr(6,3));
		int nneutrons = nnucleons-nprotons;
		static const double protonmass = pdgdb->GetParticle(2212)->Mass()*1000.;
		static const double neutronmass = pdgdb->GetParticle(2112)->Mass()*1000.;
		return ((nprotons*protonmass)+(nneutrons*neutronmass));
	}
	// else doesn't appear to be a nucleus.
	std::cerr<<"PdgToMass Warning! Particle with pdg "<<code<<" not found!"<<std::endl;
	return -1;
}

std::string RunModeToName(int mdrnsk){
	if(constants::mdrnsk_to_runtype.count(mdrnsk)){
		return constants::mdrnsk_to_runtype.at(mdrnsk);
	}
	std::cerr<<"RunModeToName: unkown mdrnsk "<<mdrnsk<<std::endl;
	return std::to_string(mdrnsk);
}

std::string TriggerIDToName(int code){
	if(constants::Trigger_ID_To_Trigger.count(code)){
		return constants::Trigger_ID_To_Trigger.at(code);
	}
	return "?";
}

int TriggerNameToID(std::string trigname){
	if(constants::Trigger_To_Trigger_ID.count(trigname)){
		return constants::Trigger_To_Trigger_ID.at(trigname);
	}
	return -1;
}

std::string NEUTInteractionModeToString(const int& mode){
  if (constants::Interaction_Mode_To_String.count(mode)){
    return constants::Interaction_Mode_To_String.at(mode);
  }
  throw std::invalid_argument("BAD NEUT INTERACTION MODE");  
}

std::string GetNEUTModeCharge(const int& mode){
  const int abs_mode = abs(mode);
  if (abs_mode == 0){return "mode = 0";}
  if (abs_mode >= 1 && abs_mode <= 27){
    return "CC";
  }
  else if (abs_mode >= 27 && abs_mode <= 52){
    return "NC";
  }
  throw std::invalid_argument("BAD NEUT INTERACTION MODE");
}

std::string GetNEUTModeProcess(const int& mode){
  const int abs_mode = abs(mode);
  if (abs_mode == 0){return "mode = 0";}
  if (abs_mode == 1){return "CCQE";}
  if (abs_mode >= 11 && abs_mode <= 13){return "CC Single Pi From Delta Resonance";}
  if (abs_mode == 16){return "CC Coherent Pion Production";}
  if (abs_mode == 21){return "CC Multi-Pion Production";}
  if (abs_mode == 22){return "CC Single Eta From Delta Resonance";}
  if (abs_mode ==23){return "CC Single K From Delta Resonance";}
  if (abs_mode == 26){return "CC Deep Inelastic";}
  if (abs_mode == 27){return "CC Diffractive Pion Production";}
  if (abs_mode >= 31 && abs_mode <= 34){return "NC Single Pi From Delta Resonance";}
  if (abs_mode == 36){return "NC Coherent Pion Production";}
  if (abs_mode == 41){return "NC Multi-Pion Production";}
  if (abs_mode == 42 || abs_mode == 43){return "NC Single Eta From Delta Resonance";}
  if (abs_mode == 44 || abs_mode == 45){return "NC Single K From Delta Resonance";}
  if (abs_mode == 46){return "NC Deep Inelastic";}
  if (abs_mode == 47){return "NC Diffractive Pion Production";}
  if (abs_mode == 51 || abs_mode == 52){return "NCQE";}
  if (abs_mode <= 52 && abs_mode >= 0){return "MYSTERY MODE";} //yeah... not ideal - need to figure out what these are.
  throw std::invalid_argument("BAD NEUT INTERACTION MODE");
}

std::string GetTriggerNames(int32_t trigid){
	std::bitset<32> triggerID(trigid);
	std::string Trigs="";
	for(int i=0; i<=31; i++){
		if(triggerID.test(i)){
			if(Trigs!="") Trigs += ", ";
			Trigs += TriggerIDToName(i);
		}
	}
	Trigs = "[" + Trigs + "]";
	return Trigs;
}

std::string EventFlagToString(int ifevsk, int sk_geometry){
	if(sk_geometry<0){
		// autodetect based on Header
		sk_geometry = skheadg_.sk_geometry;
	}
	const std::map<int,std::string>* flag_to_string = nullptr;
	if(sk_geometry>0 && sk_geometry<4){
		flag_to_string = &constants::flag_to_string_SKI_III;
	} else if(sk_geometry>3){
		flag_to_string = &constants::flag_to_string_SKIV;
	} else {
		std::cerr<<"EventFlagToString: unknown sk_geometry "<<sk_geometry<<std::endl;
	}
	if(flag_to_string){
		if(flag_to_string->count(ifevsk)) return flag_to_string->at(ifevsk);
		else std::cerr<<"EventFlagToString: unknown ifevsk "<<ifevsk<<std::endl;
	}
	return std::to_string(ifevsk);
}

std::string GetEventFlagNames(int32_t flagid){
	std::bitset<32> flagID(flagid);
	std::string Flags="";
	for(int i=0; i<=31; i++){
		if(flagID.test(i)){
			if(Flags!="") Flags += ", ";
			Flags += EventFlagToString(i);
		}
	}
	Flags = "[" + Flags + "]";
	return Flags;
}

std::unordered_map<std::string,int> GetHitFlagNames(int32_t ihtiflz, std::string* list){
	// from $SKOFL_ROOT/inc/sktq.h, meaning of bits of hitflags:
	// n.b. ihtiflz = ID, ihtflz = OD
	std::unordered_map<std::string,int> labels;
	labels.emplace("in 1.3us",(ihtiflz & 1));
	labels.emplace("in gate",(ihtiflz & 2));
	// bits 2->3 define trigger ID; only vals 0,1,2 used.
	int tmp = (ihtiflz >> 2) & 3;
	labels.emplace("narrow trigger",(tmp==0));
	labels.emplace("wide trigger",(tmp==1));
	labels.emplace("pedestal trigger",(tmp==2));
	// bits 4->5 define charge range
	tmp = (ihtiflz >> 4) & 3;
	labels.emplace("small charge range",(tmp==0));
	labels.emplace("medium charge range",(tmp==1));
	labels.emplace("large charge range",(tmp==2));
	// bits 6->11 are "(# of TRG EVENT COUNTER - 1) * 64 (0-63)" ...
	tmp = (ihtiflz >> 6);
	labels.emplace("trigger event counter",tmp);
	
	// TODO ehh, this mixes up all the ordering if that matters
	if(list!=nullptr){
		// combine flags
		for(auto&& label : labels){
			if(label.second==0) continue;
			if((*list)!="") *list+=", ";
			*list+=label.first;
		}
		*list="["+*list+"]";
	}
	return labels;
}

std::unordered_map<std::string,int> GetHitChargeAndFlags(int32_t iqiskz, int& adc_counts, int32_t& flags, std::string* list){
	// from $SKOFL_ROOT/inc/sktq.h
	// iqiskz combines raw ADC charge in bits 0-10...
	adc_counts = iqiskz & 0x7FF;
	// ...and hit flags in bits 11-15
	flags = iqiskz >> 11;
	// these flags are a subset of the flags in ihtiflz
	// (no 1.3us flag in bit 0, no TRG event counter in upper bits)
	
	std::unordered_map<std::string,int> labels;
	labels.emplace("in gate",(flags & 1));
	// bits 1->2 define trigger ID; only vals 0,1,2 used.
	int tmp = (flags >> 1) & 3;
	labels.emplace("narrow trigger",(tmp==0));
	labels.emplace("wide trigger",(tmp==1));
	labels.emplace("pedestal trigger",(tmp==2));
	// bits 4->5 define charge range
	tmp = (flags >> 4) & 3;
	labels.emplace("small charge range",(tmp==0));
	labels.emplace("medium charge range",(tmp==1));
	labels.emplace("large charge range",(tmp==2));
	
	// TODO ehh, this mixes up all the ordering if that matters
	if(list!=nullptr){
		// combine flags
		for(auto&& label : labels){
			if(label.second==0) continue;
			if((*list)!="") *list+=", ";
			*list+=label.first;
		}
		*list="["+*list+"]";
	}
	return labels;
}

int GetTriggerThreshold(int trigbit){
	if(constants::default_trig_thresholds.count(trigbit)){
		return constants::default_trig_thresholds.at(trigbit);
	}
	return 0;
}
