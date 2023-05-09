/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <set>
#include <string>
#include <map>
#include "TInterpreter.h" // TInterpreter::EErrorCode
#include "TDatabasePDG.h"
//#include <regex>    // std::regex doesn't work for older g++ versions

extern std::set<std::string> fundamental_types;
extern std::set<std::string> container_types;
extern std::map<TInterpreter::EErrorCode, std::string> TInterpreterErrors;

// functions
std::string G3_process_code_to_string(int process_code);
std::string numnu_code_to_string(int numnu_code);
std::string neut_mode_to_string(int neut_code);
std::string PdgToString(int code);
int StringToPdg(std::string name);
int PdgToG3ParticleCode(int code);
int G3ParticleCodeToPdg(int code);
int StringToG3ParticleCode(std::string name);
std::string G3ParticleCodeToString(int code);
double PdgToMass(int code);
std::string TriggerIDToTrigger(int code);

enum class SKROOTMODE : int { NONE = 4, ZEBRA = 3, READ = 2, WRITE = 1, COPY = 0 };

namespace constants{
	
	/// @brief The format code for a multievent binary BoostStore
	constexpr int BOOST_STORE_BINARY_FORMAT = 0;
	constexpr int BOOST_STORE_ASCII_FORMAT = 1;
	constexpr int BOOST_STORE_MULTIEVENT_FORMAT = 2;
	
	static const TDatabasePDG* particleDb = TDatabasePDG::Instance();
	// https://root.cern.ch/root/html532/src/TDatabasePDG.cxx.html
	// https://root.cern/doc/v608/classTDatabasePDG.html
	
	static const std::map<int,std::string> G3_process_code_to_string{
		// full list from Geant3 manual page 445
		// see p420+ for a full list of Geant3 common blocks and their parameters
		{5, "Decay"},
		{6, "Pair production"},
		{7, "Compton Scatter"},
		{8, "Photo-electric"},
		{9, "Bremsstrahlung"},
		{12, "Hadronic Interaction"},
		{13, "Hadronic Elastic Coherent Scattering"},
		{18, "Neutron Capture"},
		{20, "Hadronic Inelastic"},
		{21, "Muon-Nuclear Interaction"},
		{23, "Photonuclear"},
		{30, "Below tracking threshold"}
	};
		/*
		NEXT, 1, particle has reached the boundary of current volume
		MULS, 2, multiple scattering
		LOSS, 3, continuous energy loss
		FIEL, 4, bending in magnetic feld
		DCAY, 5, particle decay
		PAIR, 6, photon pair-production or muon direct pair production
		COMP, 7, Compton scattering
		PHOT, 8, photoelectric effect
		BREM, 9, bremsstrahlung
		DRAY, 10, δ-ray production
		ANNI, 11, positron annihilation
		HADR, 12, hadronic interaction
		ECOH, 13, hadronic elastic coherent scattering
		EVAP, 14, nuclear evaporation
		FISS, 15, nuclear fssion
		ABSO, 16, nuclear absorption
		ANNH, 17, anti-proton annihilation
		CAPT, 18, neutron capture
		EINC, 19, hadronic elastic incoherent scattering
		INHE, 20, hadronic inelastic scattering
		MUNU, 21, muon-nuclear interaction
		TOFM, 22, exceeded time of fight cut
		PFIS, 23, nuclear photo-fssion
		SCUT, 24, the particle was unexpectedly crossing volume boundaries due to bending in a magnetic feld and the step has been halved to avoid this
		RAYL, 25, Rayleigh effect
		PARA, 26, parametrisation activated
		PRED, 27, error matrix computed (GEANE tracking)
		LOOP, 28, not used
		NULL, 29, no mechanism is active, usually at the entrance of a new volume
		STOP, 30, particle has fallen below energy threshold and tracking stops
		LABS, 101, Cerenkov photon absorption
		LREF, 102, Cerenkov photon refection/refraction
		SMAX, 103, step limited by STEMAX
		SCOR, 104, correction against loss of precision in boundary crossing
		CKOV, 105, Cerenkov photon generation
		REFL, 106, Cerenkov photon refection
		REFR, 107, Cerenkov photon refraction
		SYNC, 108, synchrotron radiation generation
		STRA, 109, PAI or ASHO model used for energy loss fuctuations.
		*/
	
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
		{45,100045},   // deuteron
		{46,100046},   // triton
		{47,100047},   // alpha
		{48,100048},   // Geantino
		{49,100049},   // He3
		{69,100069}    // 16O nucleus
//		// the legitimate codes for the above  - since this is a map we cannot store both...!
//		{45,1000010020},  // deuteron
//		{46,1000010030},  // triton
//		{47,1000020040},  // alpha
//		{48,0},           // Geantino - no PDG code
//		{49,1000020030},  // He3
//		{69,1000080160}   // 16O
//		// XXX what are the g3 codes for gadolinium nuclei???
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
		{12,4},        // e-neutrino : Geant3 just has "neutrino"... which to map it to?
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
		//{0,48},         // Geantino - no PDG code
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
		{1000010020,"Deuterium"},
		{1000010030,"Tritium"},
		{1000020040,"Alpha"},
		{1000020030,"He3"},
		{1000080160,"O16"},
		{1000641550,"Gd155"},     // ground state
		{1000641550,"Gd155*"},    // excited 121.050 state
		{1000641560,"Gd156"},     // ground state
		{1000641561,"Gd156*"},    // excited 2137.600 state
		{1000641570,"Gd157"},     // ground state
		{1000641571,"Gd157*"},    // excited 426.600 state
		{1000641580,"Gd158"},     // ground state
		// support skdetsim "custom" pdg codes
		{100045,"Deuterium"},
		{100046,"Tritium"},
		{100047,"Alpha"},
		{100048,"Geantino"},
		{100049,"He3"},
		{100069,"O16"}
	};
	
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
	// could we use this to connect skdetsim daughter nuclei to their parent?
	// we will probably need to build this decay list ourselves though.
	
	enum muboy_classes{ misfit=0, single_thru_going=1, single_stopping=2, multiple_mu=3, also_multiple_mu=4, corner_clipper=5};
	const std::map<muboy_classes,std::string> muboy_class_to_name{
		{muboy_classes::misfit,"misfit"},
		{muboy_classes::single_thru_going,"single_thru_going"},
		{muboy_classes::single_stopping,"single_stopping"},
		{muboy_classes::multiple_mu,"multiple_mu"},
		{muboy_classes::also_multiple_mu,"also_multiple_mu"}
	};
	const std::map<std::string,muboy_classes> muboy_name_to_class{
		{"misfit",muboy_classes::misfit},
		{"single_thru_going",muboy_classes::single_thru_going},
		{"single_stopping",muboy_classes::single_stopping},
		{"multiple_mu",muboy_classes::multiple_mu},
		{"also_multiple_mu",muboy_classes::also_multiple_mu}
	};
	
	static const std::map<int, std::string> Trigger_ID_To_Trigger{
		// from skhead.h
		// "BLANK #" entries are empty ID numbers in skhead.h
		{0, "LE (sftw. trig.)"},
		{1, "HE (sftw. trig.)"},
		{2, "SLE (sftw. trig.)"},
		{3, "OD (sftw. trig.)"},
		{4, "Periodic (SKI-III)"},
		{5, "After/CAL (SKI-III)"},
		{6, "VETO START"},
		{7, "VETO STOP"},
		{8, "BLANK 8"},
		{9, "BLANK 9"},
		{10, "BLANK 10"},
		{11, "Random Wide Trigger"},
		{12, "Laser (ID, Usho Laser)"},
		{13, "LED"},
		{14, "Ni"},
		{15, "Laser (OD, AutoTQlaser)"},
		{16, "LE (hitsum)"},
		{17, "HE (histum)"},
		{18, "SLE (hitsum)"},
		{19, "OD (hitsum)"},
		{20, "BLANK 20"},
		{21, "BLANK 21"},
		{22, "SN Burst"},
		{23, "mu->e Decay"},
		{24, "LINAC"},
		{25, "LINAC Microwave"},
		{26, "BLANK 26"},
		{27, "Periodic (simple)"},
		{28, "SHE"},
		{29, "AFT"},
		{30, "Pedestal"},
		{31, "T2K"}
	};
	
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
