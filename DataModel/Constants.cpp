/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#include "Constants.h"
#include "DataModel.h"

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
	} else {
		return "unknown";
	}
}

std::string numnu_code_to_string(int numnu_code){
	if(constants::numnu_code_to_string.count(numnu_code)){
		return constants::numnu_code_to_string.at(numnu_code);
	} else if(numnu_code>5){
		return "fs_other";
	} else {
		return "unknown";
	}
}

std::string neut_mode_to_string(int neut_code){
	if(constants::neut_mode_to_string.count(neut_code)){
		return constants::neut_mode_to_string.at(neut_code);
	} else {
		return "unknown";
	}
}

std::string PdgToString(int code){
	if(constants::pdg_to_string.count(code)!=0){
		return constants::pdg_to_string.at(code);
	} else {
		return std::to_string(code);
	}
}

int StringToPdg(std::string name){
	if(constants::string_to_pdg.count(name)!=0){
		return constants::string_to_pdg.at(name);
	} else {
		return -1;
	}
}

std::string G3ParticleCodeToString(int code){
	if(constants::g3_particle_code_to_string.count(code)){
		return constants::g3_particle_code_to_string.at(code);
	} else {
		return std::to_string(code);
	}
}

int StringToG3ParticleCode(std::string name){
	if(constants::string_to_g3_particle_code.count(name)){
		return constants::string_to_g3_particle_code.at(name);
	} else {
		return -1;
	}
}

int G3ParticleCodeToPdg(int code){
	if(constants::g3_particle_code_to_pdg.count(code)){
		return constants::g3_particle_code_to_pdg.at(code);
	} else {
		return -1;
	}
}

int PdgToG3ParticleCode(int code){
	if(constants::pdg_to_g3_particle_code.count(code)){
		return constants::pdg_to_g3_particle_code.at(code);
	} else {
		return -1;
	}
}

double PdgToMass(int code){
	DataModel* m_data = DataModel::GetInstance();
	auto particle = m_data->pdgdb->GetParticle(code);
	if(particle==nullptr) return -1;
	return particle->Mass()*1000.;      // converted to MeV
}

std::string TriggerIDToTrigger(int code){
	if(constants::Trigger_ID_To_Trigger.count(code)){
		return constants::Trigger_ID_To_Trigger.at(code);
	}
	return "unknown";
}
