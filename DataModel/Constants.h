/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <set>
#include <string>
#include <map>
#include <sstream>
#include <unordered_map>
#include <unordered_map>
#include "TInterpreter.h" // TInterpreter::EErrorCode
#include "TDatabasePDG.h"
//#include <regex>    // std::regex doesn't work for older g++ versions

#include "fortran_routines.h"

static double SOL_IN_CM_PER_NS_IN_WATER = 22.484996; //speed of light in cm/ns

extern std::set<std::string> fundamental_types;
extern std::set<std::string> container_types;
extern std::map<TInterpreter::EErrorCode, std::string> TInterpreterErrors;

// functions
std::string G3_process_code_to_string(int process_code);
std::string G4_process_code_to_string(int process_code);
std::string numnu_code_to_string(int numnu_code);
std::string neut_mode_to_string(int neut_code);
std::string RunModeToName(int mdrnsk);
std::string EventFlagToString(int ifevsk, int sk_geometry=-1);
std::string GetEventFlagNames(int32_t flagid);
std::unordered_map<std::string,int> GetHitFlagNames(int32_t ihtiflz, std::string* list=nullptr);
std::unordered_map<std::string,int> GetHitChargeAndFlags(int32_t iqiskz, int& adc_counts, int32_t& flags, std::string* list=nullptr);
std::string PdgToString(int code);
int StringToPdg(std::string name);
int PdgToG3ParticleCode(int code);
int G3ParticleCodeToPdg(int code);
int StringToG3ParticleCode(std::string name);
std::string G3ParticleCodeToString(int code);
double PdgToMass(int code);
std::string TriggerIDToName(int code);
std::string GetTriggerNames(int32_t trigid);
int TriggerNameToID(std::string trigname);
const std::unordered_map<int,std::string>* const GetParticleNameMap();

enum class SKROOTMODE : int { NONE = 4, ZEBRA = 3, READ = 2, WRITE = 1, COPY = 0 };

// we often need to perform actions only on lowe events, muon events or AFTs
// when cuts on lowe/muon reconstruction, or including untagged muons, we can't just use the trigger
enum class EventType : int { Unknown=0, LowE=1, Muon=2, AFT=3 };

inline std::ostream& operator<< (std::ostream& out, EventType eventType){
	switch(eventType){
		case EventType::Unknown: out << "Unknown"; break;
		case EventType::LowE: out << "LowE"; break;
		case EventType::Muon: out << "Muon"; break;
		case EventType::AFT: out << "AFT"; break;
	}
	return out;
}

inline std::istream& operator>>(std::istream& str, EventType& eventType){
	std::string eventType_str;
	if(!(str >> eventType_str)){
		std::cout<<"error reading into string"<<std::endl;
	} else {
		if(eventType_str=="Unknown") eventType=EventType::Unknown;
		if(eventType_str=="LowE") eventType=EventType::LowE;
		if(eventType_str=="Muon") eventType=EventType::Muon;
		if(eventType_str=="AFT") eventType=EventType::AFT;
	}
	return str;
}


class TriggerType {
	public:
	enum TriggerTypeEnum {
		LE=0,
		HE=1,
		SLE=2,
		OD_or_Fission=3,
		Periodic=4,
		AFT_or_Cal=5,
		Veto_Start=6,
		Veto_Stop=7,
		unknown_8=8,
		unknown_9=9,
		unknown_10=10,
		Random_Wide=11,
		ID_Laser=12,
		LED=13,
		Ni=14,
		OD_Laser=15,
		LE_hitsum=16,
		HE_hitsum=17,
		SLE_hitsum=18,
		OD_hitsum=19,
		unknown_20=20,
		unknown_21=21,
		SN_Burst=22,
		mue_Decay=23,
		LINAC=24,
		LINAC_RF=25,
		unknown_26=26,
		Periodic_simple=27,
		SHE=28,
		AFT=29,
		Pedestal=30,
		T2K=31
	};
	TriggerTypeEnum trigtype;
	TriggerType() = delete;
	
	public:
	constexpr TriggerType(const TriggerTypeEnum& trig) : trigtype(TriggerTypeEnum(trig)) {};
	constexpr TriggerType(const int& trig) : trigtype(TriggerTypeEnum(trig)) {};
	int operator()() const { return static_cast<int>(trigtype); };
	TriggerType& operator=(const TriggerTypeEnum& trig) { trigtype = trig; return *this;}
	bool operator==(const TriggerTypeEnum trig) const { return trigtype == trig; }
	bool operator!=(const TriggerTypeEnum trig) const { return trigtype != trig; }
	operator int() const { return static_cast<int>(trigtype); }
	operator std::string() const { std::stringstream ss; ss << *this; return ss.str(); }
	
};

class EventFlagSKIV{
	public:
	enum EventFlagSKIVEnum{
		QBEE_TQ=0,
		HARD_TRG=1,
		QBEE_STAT=2,
		DB_STAT_BLOCK=3,
		CORRUPTED_CHECKSUM=4,
		MISSING_SPACER=5,
		PED_HIST_BLOCK=6,
		unknown_7=7,
		unknown_8=8,
		PEDESTAL_ON=9,
		RAW_AMT_BLOCK=10,
		GPS_DATA=11,
		PEDESTAL_CHECK=12,
		SEND_BLOCK=13,
		INNER_SLOW_DATA=14,
		RUN_INFORMATION=15,
		PREV_T0_BLOCK=16,
		unknown_17=17,
		FE_TRL_BLOCK=18,
		SPACER_BLOCK=19,
		INCOMPLETE_TQ=20,
		CORRUPT_TQ_BLOCK=21,
		TRG_MISMATCH_TQ=22,
		QBEE_ERROR=23,
		SORT_BLOCK=24,
		CORRUPTED_BLOCK=25,
		LED_BURST_ON=26,
		EVNT_TRAILER=27,
		INNER_DETECTOR_OFF=28,
		ANTI_DETECTOR_OFF=29,
		T2K_GPS=30,
		EVNT_HDR_AND_SOFTWARE_TRG=31
	};
	EventFlagSKIVEnum trigtype;
	EventFlagSKIV() = delete;
	
	public:
	constexpr EventFlagSKIV(const EventFlagSKIVEnum& trig) : trigtype(EventFlagSKIVEnum(trig)) {};
	constexpr EventFlagSKIV(const int& trig) : trigtype(EventFlagSKIVEnum(trig)) {};
	int operator()() const { return static_cast<int>(trigtype); };
	EventFlagSKIV& operator=(const EventFlagSKIVEnum& trig) { trigtype = trig; return *this;}
	bool operator==(const EventFlagSKIVEnum trig) const { return trigtype == trig; }
	bool operator!=(const EventFlagSKIVEnum trig) const { return trigtype != trig; }
	operator int() const { return static_cast<int>(trigtype); }
	operator std::string() const { std::stringstream ss; ss << *this; return ss.str(); }
};

enum class muboy_class{
	misfit=0,            // too few valid hits (<10) after cuts
	single_thru_going=1, // good confidence in single throughgoing
	single_stopping=2,   // 
	multiple_mu_1=3,     // Scott says: "80% of multiple muons are this type..."
	multiple_mu_2=4,     // "...and 20% are of this type". Maybe 3 is lower SNR?
	corner_clipper=5     // see lowe school slides.
};

namespace constants{
	
	/// @brief The format code for a multievent binary BStore
	constexpr int BSTORE_BINARY_FORMAT = 0;
	constexpr int BSTORE_ASCII_FORMAT = 1;
	constexpr int BSTORE_MULTIEVENT_FORMAT = 2;
	
	// https://root.cern.ch/root/html532/src/TDatabasePDG.cxx.html
	// https://root.cern/doc/v608/classTDatabasePDG.html
	
	static const std::map<int,std::string> G3_process_code_to_string{
		// full list from Geant3 manual page 445
		// see p420+ for a full list of Geant3 common blocks and their parameters
		{1, "Volume Boundary"},
		{2, "Multiple Scattering"},
		{3, "Continuous Energy Loss"},
		{4, "Magnetic Field Bending"},
		{5, "Decay"},
		{6, "Pair production"},
		{7, "Compton Scatter"},
		{8, "Photo-electric"},
		{9, "Bremsstrahlung"},
		{10, "Delta-ray Production"},
		{11, "Positron Annihilation"},
		{12, "Hadronic Interaction"},
		{13, "Hadronic Elastic Coherent"},
		{14, "Nuclear Evaporation"},
		{15, "Nuclear Fission"},
		{16, "Nuclear Absorption"},
		{17, "Anti-Proton Annihilation"},
		{18, "Neutron Capture"},
		{19, "Hadronid Elastic Incoherent"},
		{20, "Hadronic Inelastic"},
		{21, "Muon-Nuclear Interaction"},
		{22, "Exceeded Time of Fight Cut"},
		{23, "Nuclear Photo-Fission"},
		{24, "Bending in a Magnetic Field"},
		{25, "Rayleigh Effect"},
		{26, "Parametrisation Activated"},
		{27, "Error matrix computed (GEANE tracking)"},
		{28, "Not Used"},
		{29, "No Mechanism Active (usually at entrance of a new volume)"},
		{30, "Below tracking threshold"},
		{101, "Cerenkov Absorption"},
		{102, "Cerenkov Refection/Refraction"},
		{103, "Step Limited by STEMAX"},
		{104, "Correction against Loss of Precision in Boundary Crossing"},
		{105, "Cerenkov Generation"},
		{106, "Cerenkov Reflection"},
		{107, "Cerenkov Refraction"},
		{108, "Synchrotron Generation"},
		{109, "PAI or ASHO model used for energy loss fuctuations"}
	};
	
	// from skheadC.h
	static const std::map<int,std::string> mdrnsk_to_runtype{
		{0,"Monte Carlo"},
		{1,"Normal Data"},
		{2,"Laser Calibration"},
		{3,"Pedestal"},
		{4,"Xenon Lamp"},
		{5,"Nickel Calibration"},
		{6,"Random Trigger"},
		{7,"Linac Calibration"}
	};
	
	static const std::map<int,std::string> G4_process_code_to_string{
		{1,"CoulombScat"},
		{2,"Ionisation"},
		{3,"Brems"},
		{4,"PairProdCharged"},
		{5,"Annih"},
		{6,"AnnihToMuMu"},
		{7,"AnnihToHad"},
		{8,"NuclearStopp"},
		{9,"ElectronSuper"},
		{10,"Msc"},
		{11,"Rayleigh"},
		{12,"PhotoElectric"},
		{13,"Compton"},
		{14,"Conv"},
		{15,"ConvToMuMu"},
		{16,"GammaSuper"},
		{21,"Cerenkov"},
		{22,"Scintillation"},
		{23,"SynchRad"},
		{24,"TransRad"},
		{31,"OpAbsorb"},
		{32,"OpBoundary"},
		{33,"OpRayleigh"},
		{34,"OpWLS"},
		{35,"OpMieHG"},
		{51,"DNAElastic"},
		{52,"DNAExcit"},
		{53,"DNAIonisation"},
		{54,"DNAVibExcit"},
		{55,"DNAAttachment"},
		{56,"DNAChargeDec"},
		{57,"DNAChargeInc"},
		{58,"DNAElectronSolvatation"},
		{59,"DNAMolecularDecay"},
		{60,"ITTransportation"},
		{61,"DNABrownianTransportation"},
		{62,"DNADoubleIonisation"},
		{63,"DNADoubleCapture"},
		{64,"DNAIonisingTransfer"},
		{91,"Transportation"},
		{92,"CoupleTrans"},
		{111,"HadElastic"},
		{121,"HadInElastic"},
		{131,"HadCapture"},
		{132,"MuAtomicCapture"},
		{141,"HadFission"},
		{151,"HadAtRest"},
		{161,"HadCEX"},
		{201,"Decay"},
		{202,"DecayWSpin"},
		{203,"DecayPiSpin"},
		{210,"DecayRadio"},
		{211,"DecayUnKnown"},
		{221,"DecayMuAtom"},
		{231,"DecayExt"},
		{401,"StepLimiter"},
		{402,"UsrSepcCuts"},
		{403,"NeutronKiller"},
		{491,"ParallelWorld"}
	};
	
	// from skdetsim source file 'gt2pd.h'
	static const std::map<int,std::string> g3_particle_code_to_string{
		{1,"Gamma"},
		{2,"Positron"},
		{3,"Electron"},
		{4,"Neutrino"},
		{5,"Muon +"},
		{6,"Muon -"},
		{7,"Pion 0"},
		{8,"Pion +"},
		{9,"Pion -"},
		{10,"Kaon 0 Long"},
		{11,"Kaon +"},
		{12,"Kaon -"},
		{13,"Neutron"},
		{14,"Proton"},
		{15,"Antiproton"},
		{16,"Kaon 0 Short"},
		{17,"Eta"},
		{18,"Lambda"},
		{19,"Sigma +"},
		{20,"Sigma 0"},
		{21,"Sigma -"},
		{22,"Xi 0"},
		{23,"Xi -"},
		{24,"Omega -"},
		{25,"Antineutron"},
		{26,"Antilambda"},
		{27,"Antisigma -"},
		{28,"Antisigma 0"},
		{29,"Antisigma +"},
		{30,"Antixi 0"},
		{31,"Antixi +"},
		{32,"Antiomega +"},
		{33,"Tau +"},
		{34,"Tau -"},
		{35,"D +"},
		{36,"D -"},
		{37,"D 0"},
		{38,"AntiD 0"},
		{39,"Ds +"},
		{40,"Anti Ds -"},
		{41,"Lambda c +"},
		{42,"W +"},
		{43,"W -"},
		{44,"Z 0"},
		{45,"Deuteron"},
		{46,"Tritium"},
		{47,"Alpha"},
		{48,"Geantino"},
		{49,"He3"},
		{50,"Cerenkov"},
		{69,"O16"}
	};
	
	static const std::map<std::string,int> string_to_g3_particle_code{
		{"Gamma",1},
		{"Positron",2},
		{"Electron",3},
		{"Neutrino",4},
		{"Muon +",5},
		{"Muon -",6},
		{"Pion 0",7},
		{"Pion +",8},
		{"Pion -",9},
		{"Kaon 0 Long",10},
		{"Kaon +",11},
		{"Kaon -",12},
		{"Neutron",13},
		{"Proton",14},
		{"Antiproton",15},
		{"Kaon 0 Short",16},
		{"Eta",17},
		{"Lambda",18},
		{"Sigma +",19},
		{"Sigma 0",20},
		{"Sigma -",21},
		{"Xi 0",22},
		{"Xi -",23},
		{"Omega -",24},
		{"Antineutron",25},
		{"Antilambda",26},
		{"Antisigma -",27},
		{"Antisigma 0",28},
		{"Antisigma +",29},
		{"Antixi 0",30},
		{"Antixi +",31},
		{"Antiomega +",32},
		{"Tau +",33},
		{"Tau -",34},
		{"D +",35},
		{"D -",36},
		{"D 0",37},
		{"AntiD 0",38},
		{"Ds +",39},
		{"Anti Ds -",40},
		{"Lambda c +",41},
		{"W +",42},
		{"W -",43},
		{"Z 0",44},
		{"Deuteron",45},
		{"Tritium",46},
		{"Alpha",47},
		{"Geantino",48},
		{"He3",49},
		{"Cerenkov",50},
		{"O16",69}
	};
	
	static const std::map<int,int> g3_particle_code_to_pdg{
		// from https://root.cern.ch/root/html532/src/TDatabasePDG.cxx.html#228
		// or use Int_t TDatabasePDG::ConvertGeant3ToPdg(Int_t Geant3number)
		{1,22},        // photon
		{25,-2112},    // anti-neutron
		{2,-11},       // e+
		{26,-3122},    // anti-Lambda
		{3,11},        // e-
		{27,-3222},    // Sigma-
		{4,12},        // e-neutrino : Geant3 just has "neutrino"... which to map it to?
		{28,-3212},    // Sigma0
		{5,-13},       // mu+
		{29,-3112},    // Sigma+ (PB)*/
		{6,13},        // mu-
		{30,-3322},    // Xi0
		{7,111},       // pi0
		{31,-3312},    // Xi+
		{8,211},       // pi+
		{32,-3334},    // Omega+ (PB)
		{9,-211},      // pi-
		{10,130},      // K long
		{33,-15},      // tau+     << these are switched in gt2pd.h
		{34,15},       // tau-     << ROOT says 33=15, 34=-15
		{11,321},      // K+
		{35,411},      // D+
		{12,-321},     // K-
		{36,-411},     // D-
		{13,2112},     // n
		{37,421},      // D0
		{14,2212},     // p
		{38,-421},     // D0
		{15,-2212},    // anti-proton
		{39,431},      // Ds+
		{16,310},      // K short
		{40,-431},     // anti Ds-
		{17,221},      // eta
		{41,4122},     // Lamba_c+
		{18,3122},     // Lambda
		{42,24},       // W+
		{19,3222},     // Sigma+
		{43,-24},      // W-
		{20,3212},     // Sigma0
		{44,23},       // Z
		{21,3112},     // Sigma-
		{22,3322},     // Xi0
		{23,3312},     // Xi-
		{24,3334},     // Omega- (PB)
		// XXX The following are from skdetsim gt2pd.h, but are NOT legitimate PDG codes! XXX
		/*
		{45,100045},   // deuteron
		{46,100046},   // triton
		{47,100047},   // alpha
		{48,100048}    // Geantino
		{49,100049},   // He3
		{69,100069},   // 16O nucleus
		*/
		// the legitimate codes for the above
		{45,1000010020},  // deuteron
		{46,1000010030},  // triton
		{47,1000020040},  // alpha
		{48,0},            // Geantino
		{49,1000020030},  // He3
		{69,1000080160}   // 16O
		// "what are the g3 codes for gadolinium nuclei?" - there are none! no Gd nucleus is saved.
	};
	
	static const std::map<int,int> pdg_to_g3_particle_code{
		// from https://root.cern.ch/root/html532/src/TDatabasePDG.cxx.html#228
		// or use Int_t TDatabasePDG::ConvertPdgToGeant3(Int_t pdgNumber)
		{22,1},        // photon
		{-2112,25},    // anti-neutron
		{-11,2},       // e+
		{-3122,26},    // anti-Lambda
		{11,3},        // e-
		{-3222,27},    // Sigma-
		{12,4},        // e-neutrino : Geant3 just has "neutrino"
		{14,4},        // mu-neutrino
		{16,4},        // tau-neutrino
		{-12,4},       // if we add these mappings we can't distinguish
		{-14,4},       // neutrinos from antineutrinos, but perhaps
		{-16,4},       // that's better than simply not knowing anything?
		{-3212,28},    // Sigma0
		{-13,5},       // mu+
		{-3112,29},    // Sigma+ (PB)*/
		{13,6},        // mu-
		{-3322,30},    // Xi0
		{111,7},       // pi0
		{-3312,31},    // Xi+
		{211,8},       // pi+
		{-3334,32},    // Omega+ (PB)
		{-211,9},      // pi-
		{-15,33},      // tau+
		{130,10},      // K long
		{15,34},       // tau-
		{321,11},      // K+
		{411,35},      // D+
		{-321,12},     // K-
		{-411,36},     // D-
		{2112,13},     // n
		{421,37},      // D0
		{2212,14},     // p
		{-421,38},     // D0
		{-2212,15},    // anti-proton
		{431,39},      // Ds+
		{310,16},      // K short
		{-431,40},     // anti Ds-
		{221,17},      // eta
		{4122,41},     // Lamba_c+
		{3122,18},     // Lambda
		{24,42},       // W+
		{3222,19},     // Sigma+
		{-24,43},      // W-
		{3212,20},     // Sigma0
		{23,44},       // Z
		{3112,21},     // Sigma-
		{3322,22},     // Xi0
		{3312,23},     // Xi-
		{3334,24},     // Omega- (PB)
		// XXX The following are from skdetsim gt2pd.h, but are NOT legitimate PDG codes! XXX
		{100045,45},   // deuteron - real PDG code is '1000010020'
		{100046,46},   // triton
		{100047,47},   // alpha
		{100048,48},   // Geantino
		{100049,49},   // He3
		{100069,69},   // 16O
		// the legitimate codes for the above
		{1000010020,45},  // deuteron
		{1000010030,46},  // triton
		{1000020040,47},  // alpha
		{0,48},           // Geantino
		{1000020030,49},  // He3
		{1000080160,69}   // 16O
	};
	
	static const std::map<int,std::string> numnu_code_to_string{
		{1,"is_neutrino"},   // initial state (incident) neutrino
		{2,"is_target"},     // initial state (struck) target
		{3,"fs_lepton"},     // final state (outgoing) target
		{4,"fs_target"},     // final state (outgoing) lepton
		{5,"fs_other"}       // final state (outgoing) other particle (codes 5 and up)
	};
	
	static const std::map<int,std::string> neut_mode_to_string{
		{1, "CC quasi-elastic"},
		{11, "CC single pi from delta resonance"},
		{12, "CC single pi from delta resonance"},
		{13, "CC single pi from delta resonance"},
		{16, "CC coherent pi production"},
		{21, "CC multi pi production"},
		{27, "CC diffractive pion production"},
		{31, "NC single pi from delta resonance"},
		{32, "NC single pi from delta resonance"},
		{33, "NC single pi from delta resonance"},
		{34, "NC single pi from delta resonance"},
		{36, "NC coherent pi"},
		{41, "NC multi pi production"},
		{47, "NC diffractive pion production"},
		{51, "NC elastic"},
		{52, "NC elastic"}
	};
	
	static const std::map<muboy_class,std::string> muboy_class_to_name{
		{muboy_class::misfit,"misfit"},
		{muboy_class::single_thru_going,"single_thru_going"},
		{muboy_class::single_stopping,"single_stopping"},
		{muboy_class::multiple_mu_1,"multiple_mu_1"},
		{muboy_class::multiple_mu_2,"multiple_mu_2"},
		{muboy_class::corner_clipper,"corner_clipper"}
	};
	
	static const std::map<std::string,muboy_class> muboy_name_to_class{
		{"misfit",muboy_class::misfit},
		{"single_thru_going",muboy_class::single_thru_going},
		{"single_stopping",muboy_class::single_stopping},
		{"multiple_mu_1",muboy_class::multiple_mu_1},
		{"multiple_mu_2",muboy_class::multiple_mu_2},
		{"corner_clipper",muboy_class::corner_clipper}
	};
	
	static const std::unordered_map<int,std::string>* const pdg_to_string = GetParticleNameMap();
	/* populating this is now done in pdg_to_name_nuclei.cc, as its real big
	static const std::map<int,std::string> pdg_to_string{
		// FIXME use TParticlePDG for greater coverage, but need to add nuclei
		{2212,"Proton"},
		{-2212,"Anti Proton"},
		{11,"Electron"},
		{-11,"Positron"},
		{12,"Electron Neutrino"},
		{-12,"Anti Electron Neutrino"},
		{22,"Gamma"},
		{2112,"Neutron"},
		{-2112,"Anti Neutron"},
		{-13,"Muon+"},
		{13,"Muon-"},
		{130,"Kaonlong"},
		{211,"Pion+"},
		{-211,"Pion-"},
		{321,"Kaon+"},
		{-321,"Kaon-"},
		{3122,"Lambda"},
		{-3122,"Antilambda"},
		{310,"Kaonshort"},
		{3112,"Sigma-"},
		{3222,"Sigma+"},
		{3212,"Sigma0"},
		{111,"Pion0"},
		{311,"Kaon0"},
		{-311,"Antikaon0"},
		{14,"Muon Neutrino"},
		{-14,"Anti Muon Neutrino"},
		{-3222,"Anti Sigma-"},
		{-3212,"Anti Sigma0"},
		{-3112,"Anti Sigma+"},
		{3322,"Xsi0"},
		{-3322,"Anti Xsi0"},
		{3312,"Xsi-"},
		{-3312,"Xsi+"},
		{3334,"Omega-"},
		{-3334,"Omega+"},
		{-15,"Tau+"},
		{15,"Tau-"},
		{100,"OpticalPhoton"},
		// nuclei
		{1000641530,"Gd153"},
		{1000641550,"Gd155"},     // ground state
		{1000641550,"Gd155*"},    // excited 121.050 state
		{1000641560,"Gd156"},     // ground state
		{1000641561,"Gd156*"},    // excited 2137.600 state
		{1000641570,"Gd157"},     // ground state
		{1000641571,"Gd157*"},    // excited 426.600 state
		{1000641580,"Gd158"},     // ground state
		{1000080170,"O17"},
		// support skdetsim "custom" pdg codes
		{100045,"Deuterium"},
		{100046,"Tritium"},
		{100047,"Alpha"},
		{100048,"Geantino"},
		{100049,"He3"},
		{100069,"O16"}
	};
	*/
	
	static const std::map<std::string,int> string_to_pdg{
		// FIXME use TParticlePDG for greater coverage, but need to add nuclei
		{"Proton",2212},
		{"Anti Proton",-2212},
		{"Electron",11},
		{"Positron",-11},
		{"Electron Neutrino",12},
		{"Anti Electron Neutrino",-12},
		{"Gamma",22},
		{"Neutron",2112},
		{"Anti Neutron",-2112},
		{"Muon+",-13},
		{"Muon-",13},
		{"Kaonlong",130},
		{"Pion+",211},
		{"Pion-",-211},
		{"Kaon+",321},
		{"Kaon-",-321},
		{"Lambda",3122},
		{"Antilambda",-3122},
		{"Kaonshort",310},
		{"Sigma-",3112},
		{"Sigma+",3222},
		{"Sigma0",3212},
		{"Pion0",111},
		{"Kaon0",311},
		{"Antikaon0",-311},
		{"Muon Neutrino",14},
		{"Anti Muon Neutrino",-14},
		{"Anti Sigma-",-3222},
		{"Anti Sigma0",-3212},
		{"Anti Sigma+",-3112},
		{"Xsi0",3322},
		{"Anti Xsi0",-3322},
		{"Xsi-",3312},
		{"Xsi+",-3312},
		{"Omega-",3334},
		{"Omega+",-3334},
		{"Tau+",-15},
		{"Tau-",15},
		{"OpticalPhoton",100},
		// nuclei
		{"Deuterium",1000010020},
		{"Tritium",1000010030},
		{"Alpha",1000020040},
		{"He3",1000020030},
		{"O16",1000080160},
		{"Gd155",1000641550},     // ground state
		{"Gd155*",1000641550},    // excited 121.050 state
		{"Gd156",1000641560},     // ground state
		{"Gd156*",1000641561},    // excited 2137.600 state
		{"Gd157",1000641570},     // ground state
		{"Gd157*",1000641571},    // excited 426.600 state
		{"Gd158",1000641580}      // ground state
	};
	
	// TODO TParticleTable allows a 'DecayList' of daughter nuclides
	// could we use this to connect daughter nuclei to their parent?
	// we will probably need to build this decay list ourselves though.
	
	static const std::map<int, std::string> Trigger_ID_To_Trigger{
		// from skheadC.h
		// "# (Unknown)" entries are empty ID numbers in skheadC.h
		{-1, "?"},
		{0, "LE"},
		{1, "HE"},
		{2, "SLE"},
		{3, "OD_or_Fission"},  // in normal run/SK-IV+: OD trigger. in SKI-III Nickel run: fission trigger
		{4, "Periodic"},       // one of {nothing (null) trigger, TQ map laser, water atten. laser, Xenon ball}
		{5, "AFT_or_Cal"},     // in normal run: AFT trigger, in calib run, one of: {laser, Xenon, Nickel, Linac}
		{6, "Veto_Start"},
		{7, "Veto_Stop"},
		{8, "unknown_8"},   //  (cable 15005)
		{9, "unknown_9"},   //  (cable 15006)
		{10, "unknown_10"}, //  (cable 15007)
		{11, "Random_Wide"},
		{12, "ID_Laser"},   // (Usho Laser)
		{13, "LED"},
		{14, "Ni"},
		{15, "OD_Laser"},   // (AutoTQlaser)
		{16, "LE_hitsum"},
		{17, "HE_hitsum"},
		{18, "SLE_hitsum"},
		{19, "OD_hitsum"},
		{20, "unknown_20"},
		{21, "unknown_21"},
		{22, "SN_Burst"},
		{23, "mue_Decay"},
		{24, "LINAC"},
		{25, "LINAC_RF"},  // (LINAC Microwave)
		{26, "unknown_26"}, // (cable 15023)
		{27, "Periodic_simple"},
		{28, "SHE"},      // (SW)
		{29, "AFT"},      // (SW)
		{30, "Pedestal"}, // (SW)
		{31, "T2K"}       // (SW)
	};
	
	static const std::map<std::string, int> Trigger_To_Trigger_ID{
		{"?", -1},
		{"LE", 0},
		{"HE", 1},
		{"SLE", 2},
		{"OD_or_Fission", 3},
		{"Periodic", 4},
		{"AFT_or_Cal", 5},
		{"Veto_Start", 6},
		{"Veto_Stop", 7},
		{"unknown_8", 8},   // (cable # 15005)
		{"unknown_9", 9},   // (cable # 15006)
		{"unknown_10", 10}, // (cable # 15007)
		{"Random_Wide", 11},
		{"ID_Laser", 12},   //  (Usho Laser)
		{"LED", 13},
		{"Ni", 14},
		{"OD_Laser", 15}, // (AutoTQlaser)
		{"LE_hitsum", 16},
		{"HE_hitsum", 17},
		{"SLE_hitsum", 18},
		{"OD_hitsum", 19},
		{"unknown_20", 20},
		{"unknown_21", 21},
		{"SN_Burst", 22},
		{"mue_Decay", 23},
		{"LINAC", 24},
		{"LINAC_RF", 25},   // (LINAC Microwave)
		{"unknown_26", 26}, // (cable # 15023)
		{"Periodic_simple", 27},
		{"SHE", 28}, // (SW)
		{"AFT", 29}, // (SW)
		{"Pedestal", 30}, // (SW)
		{"T2K", 31} // (SW)
	};
	
	static const std::map<int, std::string> flag_to_string_SKI_III{
		{0,"ATM"},
		{1,"TRG"},
		{2,"SMP REGISTER"},
		{3,"SCALER"},
		{4,"PEDESTAL START"},
		{5,"PEDESTAL DATA(ATM)"},
		{6,"PEDESTAL HISTOGRAM"},
		{7,"PEDESTAL END"},
		{8,"END OF RUN"},
		{9,"PEDESTAL(ON)"},
		{10,"10 (unknown)"},
		{11,"GPS DATA"},
		{12,"CAMAC ADC"},
		{13,"ANTI DATA"},
		{14,"INNER SLOW DATA"},
		{15,"RUN INFORMATION"},
		{16,"ERROR (TKO-PS)"},
		{17,"ERROR (HV-PS)"},
		{18,"ERROR (TEMPERARTURE)"},
		{19,"19 (unknown)"},
		{20,"UNCOMPLETED ATM DATA"},
		{21,"INVALID ATM DATA"},
		{22,"22 (unknown)"},
		{23,"23 (unknown)"},
		{24,"ERROR (DATA)"},
		{25,"UNREASONABLE DATA"},
		{26,"LED BURST ON"},
		{27,"27 (unknown)"},
		{28,"INNER DETECTOR OFF"},
		{29,"ANTI  DETECTOR OFF"},
		{30,"30 (unknown)"},
		{31,"TRG IS AVAILABLE"}
	};
	
	static const std::map<int, std::string> flag_to_string_SKIV{
		{0,"QBEE TQ"},
		{1,"HARD TRG"},
		{2,"QBEE STAT"},
		{3,"DB_STAT_BLOCK"},
		{4,"CORRUPTED_CHECKSUM"},
		{5,"MISSING SPACER"},
		{6,"PED_HIST_BLOCK"},
		{7,"7 (unknown)"},
		{8,"8 (unknown)"},
		{9,"PEDESTAL ON"},
		{10,"RAW_AMT_BLOCK"},
		{11,"GPS DATA"},
		{12,"PEDESTAL_CHECK"},
		{13,"SEND_BLOCK"},
		{14,"INNER SLOW DATA"},
		{15,"RUN INFORMATION"},
		{16,"PREV T0 BLOCK"},
		{17,"17 (unknown)"},
		{18,"FE_TRL_BLOCK"},
		{19,"SPACER_BLOCK"},
		{20,"INCOMPLETE TQ"},
		{21,"CORRUPT TQ BLOCK"},
		{22,"TRG MISMATCH TQ"},
		{23,"QBEE ERROR"},
		{24,"SORT_BLOCK"},
		{25,"CORRUPTED_BLOCK"},
		{26,"LED BURST ON"},
		{27,"EVNT TRAILER"},
		{28,"INNER DETECTOR OFF"},
		{29,"ANTI  DETECTOR OFF"},
		{30,"T2K GPS"},
		{31,"(EVNT HDR)&(SOFTWARE TRG)"}
	};
	
	static const std::map<std::string, int> string_to_flag_SKI_III{
		{"ATM", 0},
		{"TRG", 1},
		{"SMP_REGISTER", 2},
		{"SCALER", 3},
		{"PEDESTAL_START", 4},
		{"PEDESTAL_DATA(ATM)", 5},
		{"PEDESTAL_HISTOGRAM", 6},
		{"PEDESTAL_END", 7},
		{"END_OF_RUN", 8},
		{"PEDESTAL_ON", 9},
		{"unknown_10", 10},
		{"GPS_DATA", 11},
		{"CAMAC_ADC", 12},
		{"ANTI_DATA", 13},
		{"INNER_SLOW_DATA", 14},
		{"RUN_INFORMATION", 15},
		{"ERROR_TKO-PS", 16},
		{"ERROR_HV-PS", 17},
		{"ERROR_TEMPERATURE", 18},
		{"unknown_19", 19},
		{"UNCOMPLETED_ATM_DATA", 20},
		{"INVALID_ATM_DATA", 21},
		{"unknown_22", 22},
		{"unknown_23", 23},
		{"ERROR_DATA", 24},
		{"UNREASONABLE_DATA", 25},
		{"LED_BURST_ON", 26},
		{"unknown_27", 27},
		{"INNER_DETECTOR_OFF", 28},
		{"ANTI_DETECTOR_OFF", 29},
		{"unknown_30", 30},
		{"TRG_IS_AVAILABLE", 31}
	};
	
	static const std::map<std::string, int> string_to_flag_SKIV{
		{"QBEE_TQ", 0},
		{"HARD_TRG", 1},
		{"QBEE_STAT", 2},
		{"DB_STAT_BLOCK", 3},
		{"CORRUPTED_CHECKSUM", 4},
		{"MISSING_SPACER", 5},
		{"PED_HIST_BLOCK", 6},
		{"unknown_7", 7},
		{"unknown_8", 8},
		{"PEDESTAL_ON", 9},
		{"RAW_AMT_BLOCK", 10},
		{"GPS_DATA", 11},
		{"PEDESTAL_CHECK", 12},
		{"SEND_BLOCK", 13},
		{"INNER_SLOW_DATA", 14},
		{"RUN_INFORMATION", 15},
		{"PREV_T0_BLOCK", 16},
		{"unknown_17", 17},
		{"FE_TRL_BLOCK", 18},
		{"SPACER_BLOCK", 19},
		{"INCOMPLETE_TQ", 20},
		{"CORRUPT_TQ_BLOCK", 21},
		{"TRG_MISMATCH_TQ", 22},
		{"QBEE_ERROR", 23},
		{"SORT_BLOCK", 24},
		{"CORRUPTED_BLOCK", 25},
		{"LED_BURST_ON", 26},
		{"EVNT_TRAILER", 27},
		{"INNER_DETECTOR_OFF", 28},
		{"ANTI_DETECTOR_OFF", 29},
		{"T2K_GPS", 30},
		{"EVNT_HDR_&_SOFTWARE_TRG", 31}
	};

  
	// geometry information from geotnkC.h: #define'd constants
	// ah... but we can't #include geotnkC.h because it conflicts with SK2p2MeV.h.... -__-
	// someone reused the #define'd macro names as variables
	// (e.g. `const Float_t RINTK`, when geotnkC.h does `#define RINTK`)
	/*
	double SK_TANK_DIAMETER = DITKTK;
	double SK_TANK_HEIGHT = HITKTK;  // height of water volume
	double SK_ID_DIAMETER = DIINTK;
	double SK_ID_HEIGHT = HIINTK;
	double SK_OD_THICKNESS_BOTTOM = TBATTK;
	double SK_OD_THICKNESS_BARREL = TWATTK;
	double SK_OD_THICKNESS_TOP = TTATTK;
	*/
	
	// TODO: mapping of calibration / veto PMT numbers in skvetoC.h
	
	
	// maps of regular expression error codes to descriptive strings
	// -------------------------------------------------------------
	/*
	const std::map<std::regex_constants::error_type,std::string> regex_err_strings{
		{std::regex_constants::error_collate, "The expression contained an invalid collating element name."},
		{std::regex_constants::error_ctype, "The expression contained an invalid character class name."},
		{std::regex_constants::error_escape, "The expression contained an invalid escaped character, or a trailing escape."},
		{std::regex_constants::error_backref, "The expression contained an invalid back reference."},
		{std::regex_constants::error_brack, "The expression contained mismatched brackets ([ and ])."},
		{std::regex_constants::error_paren, "The expression contained mismatched parentheses (( and ))."},
		{std::regex_constants::error_brace, "The expression contained mismatched braces ({ and })."},
		{std::regex_constants::error_badbrace, "The expression contained an invalid range between braces ({ and })."},
		{std::regex_constants::error_range, "The expression contained an invalid character range."},
		{std::regex_constants::error_space, "There was insufficient memory to convert the expression into a finite state machine."},
		{std::regex_constants::error_badrepeat, "The expression contained a repeat specifier (one of *?+{) that was not preceded by a valid regular expression."},
		{std::regex_constants::error_complexity, "The complexity of an attempted match against a regular expression exceeded a pre-set level."},
		{std::regex_constants::error_stack, "There was insufficient memory to determine whether the regular expression could match the specified character sequence."}
	};
	
	// boost version
	const std::map<boost::regex_constants::error_type,std::string> bregex_err_strings{
		{boost::regex_constants::error_collate, "The expression contained an invalid collating element name."},
		{boost::regex_constants::error_ctype, "The expression contained an invalid character class name."},
		{boost::regex_constants::error_escape, "The expression contained an invalid escaped character, or a trailing escape."},
		{boost::regex_constants::error_backref, "The expression contained an invalid back reference."},
		{boost::regex_constants::error_brack, "The expression contained mismatched brackets ([ and ])."},
		{boost::regex_constants::error_paren, "The expression contained mismatched parentheses (( and ))."},
		{boost::regex_constants::error_brace, "The expression contained mismatched braces ({ and })."},
		{boost::regex_constants::error_badbrace, "The expression contained an invalid range between braces ({ and })."},
		{boost::regex_constants::error_range, "The expression contained an invalid character range."},
		{boost::regex_constants::error_space, "There was insufficient memory to convert the expression into a finite state machine."},
		{boost::regex_constants::error_badrepeat, "The expression contained a repeat specifier (one of *?+{) that was not preceded by a valid regular expression."},
		{boost::regex_constants::error_complexity, "The complexity of an attempted match against a regular expression exceeded a pre-set level."},
		{boost::regex_constants::error_stack, "There was insufficient memory to determine whether the regular expression could match the specified character sequence."},
		{boost::regex_constants::error_bad_pattern, "Invalid regex pattern."}
	};
	*/
}

// from the PDG on Monte Carlo Codes:
/*
14. Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.
For a nucleus consisting of np protons, nn neutrons and
nΛ Λ’s, A = np+nn+nΛ gives the total baryon number, Z = np
the total charge and L = nΛ the total number of strange quarks.
I gives the isomer level, with I = 0 corresponding to the ground
state and I > 0 to excitations, see[8], where states denoted
m,n,p,q translate to I = 1 − 4. As examples, the deuteron
is 1000010020 and 235U is 1000922350.
*/

#endif // define CONSTANTS_H
