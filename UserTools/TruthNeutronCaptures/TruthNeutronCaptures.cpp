/* vim:set noexpandtab tabstop=4 wrap */
#include "TruthNeutronCaptures.h"
#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <algorithm> // std::find

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"

TruthNeutronCaptures::TruthNeutronCaptures():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool TruthNeutronCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEvents",MAX_EVENTS);           // terminate after processing at most this many events
	m_variables.Get("writeFrequency",WRITE_FREQUENCY); // how many events to TTree::Fill between TTree::Writes
	
	// get the list of input files from the CStore
	// -------------------------------------------
	// filled if using LoadFileList tool
	// TODO yet to implement support for this in MTreeReader
	if(inputFile==""){
		get_ok = m_data->CStore.Get("InputFileList", input_file_names);
		if(not get_ok){
			Log(toolName+" Error: No inputFile given and no InputFileList in CStore!",v_error,verbosity);
			return false;
		}
	}
	
	// open the input TFile and TTree
	// ------------------------------
	get_ok = myTreeReader.Load(inputFile, "h1"); // official ntuple TTree is descriptively known as 'h1'
	DisableUnusedBranches();
	if(get_ok) ReadEntryNtuple(0);
	
	// create the output TFile and TTree
	// ---------------------------------
	CreateOutputFile(outputFile);
	
	return true;
}

bool TruthNeutronCaptures::Execute(){
	
	Log(toolName+" processing entry "+toString(entry_number),v_debug,verbosity);
	
	// clear output vectors so we don't carry anything over
	Log(toolName+" clearing output vectors",v_debug,verbosity);
	ClearOutputTreeBranches();
	
	// Copy over directly transferred variables
	Log(toolName+" copying output variables",v_debug,verbosity);
	CopyVariables();
	
	// Calculate derived variables
	Log(toolName+" calculating output variables",v_debug,verbosity);
	CalculateVariables();
	
	// print the current event
	if(verbosity>1) PrintBranches();
	
	// Fill the output tree
	Log(toolName+" filling output TTree entry",v_debug,verbosity);
	outtree->Fill();
	
	// update the output file so we don't lose everything if we crash
	if((entry_number%WRITE_FREQUENCY)==0) WriteTree();
	
	// stop at user-defined limit to the number of events to process
	++entry_number;
	if((MAX_EVENTS>0)&&(entry_number>=MAX_EVENTS)){
		Log(toolName+" reached MAX_EVENTS, setting StopLoop",v_error,verbosity);
		m_data->vars.Set("StopLoop",1);
	} else {
		// Pre-Load next input entry so we can stop the toolchain
		// if we're about to run off the end of the tree or encounter a read error
		get_ok = ReadEntryNtuple(entry_number);
		if(get_ok<1&&get_ok>-3){
			m_data->vars.Set("StopLoop",1);
			Log(toolName+" Hit end of input file, stopping loop",v_warning,verbosity);
		}
		else if(get_ok==-10){
			Log(toolName+" Error during AutoClear while loading next input ntuple entry!",v_error,verbosity);
			return false;
		}
		else if(get_ok<0){
			Log(toolName+" IO error loading next input ntuple entry!",v_error,verbosity);
			return false;
		}
	}
	
	return true;
}


bool TruthNeutronCaptures::Finalise(){
	
	// ensure everything is written to the output file
	// -----------------------------------------------
	get_ok = WriteTree();
	if(not get_ok){
		Log(toolName+" Error writing output TTree!",v_error,verbosity);
	}
	
	// Close and delete the file handle
	// --------------------------------
	CloseFile();
	
	return true;
}

void TruthNeutronCaptures::CopyVariables(){
	// copy over variables from the input tree
	// ---------------------------------------
	// those we want to keep in the output tree without modification
	
	out_filename = myTreeReader.GetTree()->GetCurrentFile()->GetName();
	out_skdetsim_version =  -1;
	out_tba_table_version = -1;
	out_water_transparency = water_transparency;
	
	// event level
	out_run_number = run_number;
//	out_subrun_number = subrun_number;
//	out_event_number = event_number;
	out_entry_number = entry_number;
	out_subevent_number = subevent_number;
	
}

int TruthNeutronCaptures::CalculateVariables(){
	// calculate remaining variables
	// -----------------------------
	// those we want to save to the output tree but need to derive
	
	int neutron_pdg = 2112; //StringToPdg("Neutron");
	int neutron_g3 = 13; //StringToG3ParticleCode("Neutron");
	int gamma_pdg = 22; //StringToPdg("Gamma");
	double neutron_mass = PdgToMass(neutron_pdg);
	
	// note the primary event vertex. Primary particles don't have individual start vertices,
	// but they *should* all start from the primary event vertex ... unless the geant3 event
	// had primaries at different locations. Which is entirely reasonable. In that case,
	// we may need to use the daughters' "parent vertex at birth" branch to find our primary
	// start vertices.... FIXME for events that should have all particles originating from a
	// single primary event vertex, these are not consistent... why not??
	TLorentzVector primary_vertex_tvector(primary_event_vertex.at(0),
										  primary_event_vertex.at(1),
										  primary_event_vertex.at(2),
										  0.f); // is this correct? does this define T=0?
	
	// loop over event particles and scan for neutrons
	// we want to nest them in their parent nuclides, but we need to find the neutrons first
	// (no way to tell if a nuclide has a daughter neutron)
	// so scan for the neutrons first, then group them by parent after.
	// This means storing the neutron info in temporay container for now
	
	// some helper vectors
	std::map<int,int> primary_n_ind_to_loc;
	std::map<int,int> secondary_n_ind_to_loc;
	std::vector<int> neutron_parent_indices;         // index of neutron parent
	std::vector<bool> neutron_terminfo_unknown;      // whether we've yet to set the neutron stopping info
	
	// ==========================
	// SCAN FOR PRIMARY NEUTRONS
	// ==========================
	// loop over the primary particles first, because for IBD events we have a primary neutron
	Log(toolName+" event had "+toString(n_outgoing_primaries)+" primary particles",v_debug,verbosity);
	for(int primary_i=0; primary_i<n_outgoing_primaries; ++primary_i){
		Log(toolName+" primary "+toString(primary_i)+" had G3 code "+toString((int)primary_G3_code.at(primary_i))
				+" ("+G3ParticleCodeToString(primary_G3_code.at(primary_i))+")",v_debug,verbosity);
		if(primary_G3_code.at(primary_i)==neutron_g3){
			// we found a neutron!
			Log(toolName+" NEUTRON!", v_debug,verbosity);
			primary_n_ind_to_loc.emplace(primary_i,out_neutron_start_energy.size());
			// note neutron info
			out_neutron_start_pos.push_back(primary_vertex_tvector);
			double startE = primary_start_mom.at(primary_i);
			// FIXME primary momenta appear to be in different units to secondaries...? Need to convert to MeV
			out_neutron_start_energy.push_back(startE);
			out_neutron_ndaughters.push_back(0);  // XXX which branch???
			neutron_parent_indices.push_back(-1); // n/a for primary
			
			// termination information isn't saved, so we'll have to use the creation
			// info of a matched secondary gamma to infer neutron termination information
			// for now, to keep synchronization, fill with placeholders
			out_neutron_end_pos.push_back(TLorentzVector(0,0,0,0));
			out_neutron_end_energy.push_back(0);
			out_neutron_end_process.push_back(0); // do we have any way to determine this?
			neutron_terminfo_unknown.push_back(true);
			
			// a neutron may have many daughters, so we'll have to find those in a subsequent scan
			out_gamma_energy.push_back(std::vector<double>{});
			out_gamma_time.push_back(std::vector<double>{});
			// same for conversion electrons
			out_electron_energy.push_back(std::vector<double>{});
			out_electron_time.push_back(std::vector<double>{});
			// likewise we'll note the daughter nuclide from the capture in the subsequent scan
			out_nuclide_daughter_pdg.push_back(-1);
			out_nuclide_parent_pdg.push_back(-1);
		}
	}
	Log(toolName+" found "+toString(out_neutron_start_energy.size())+" primary neutrons", v_debug,verbosity);
	
	// ===========================
	// SCAN FOR SECONDARY NEUTRONS
	// ===========================
	Log(toolName+" event had "+toString(n_secondaries_2)+" secondary particles",v_debug,verbosity);
	for(int secondary_i=0; secondary_i<n_secondaries_2; ++secondary_i){
		Log(toolName+" secondary "+toString(secondary_i)+" had pdg "+toString(secondary_PDG_code_2.at(secondary_i))
				+" ("+PdgToString(secondary_PDG_code_2.at(secondary_i))+")",v_debug,verbosity);
		if(secondary_PDG_code_2.at(secondary_i)==neutron_pdg){
			// we found a neutron!
			Log(toolName+" NEUTRON!", v_debug,verbosity);
			secondary_n_ind_to_loc.emplace(secondary_i,out_neutron_start_energy.size());
			// note neutron info
			TLorentzVector startpos(secondary_start_vertex_2.at(secondary_i).at(0),
									secondary_start_vertex_2.at(secondary_i).at(1),
									secondary_start_vertex_2.at(secondary_i).at(2),
									secondary_start_time_2.at(secondary_i));
			out_neutron_start_pos.push_back(startpos);
			TVector3 startmom(secondary_start_mom_2.at(secondary_i).at(0),
							  secondary_start_mom_2.at(secondary_i).at(1),
							  secondary_start_mom_2.at(secondary_i).at(2));
			double startE = startmom.Mag2() / (2.*neutron_mass);
			out_neutron_start_energy.push_back(startE);
			//out_neutron_ndaughters.push_back(secondary_n_daughters.at(secondary_i));
			out_neutron_ndaughters.push_back(0); // manually count since can't do this for primaries
			
			// note its parent information
			// XXX we may have secondary neutrons from inelastic scattering...
			// these may well be the "same" neutron (one in, one out)
			// to avoid double-counting, we should try to identify and skip these...
			// should we skip neutrons whose parent is also a neutron?
			// can we use 'secondary_first_daugher_index' to quickly check for/grab
			// the daughter and check its creation process to determine this neutron's
			// termination process?
			neutron_parent_indices.push_back(parent_index.at(secondary_i)-1);
			
			// termination information isn't saved, so we'll have to use the creation
			// info of a matched secondary gamma to infer neutron termination information
			// for now, to keep synchronization, fill with placeholders
			out_neutron_end_pos.push_back(TLorentzVector(0,0,0,0));
			out_neutron_end_energy.push_back(0);
			out_neutron_end_process.push_back(0); // do we have any way to determine this?
			neutron_terminfo_unknown.push_back(true);
			
			// a neutron may have many daughters, so we'll have to find those in a subsequent scan
			out_gamma_energy.push_back(std::vector<double>{});
			out_gamma_time.push_back(std::vector<double>{});
			// same for conversion electrons
			out_electron_energy.push_back(std::vector<double>{});
			out_electron_time.push_back(std::vector<double>{});
			out_nuclide_parent_pdg.push_back(-1);
		}
	}
	Log(toolName+" found "+toString(out_neutron_start_energy.size())+" primary+secondary neutrons", v_debug,verbosity);
	
	// ==================================
	// SCAN FOR GAMMAS AND CAPTURE NUCLEI
	// ==================================
	// scan again, this time looking for gammas
	int n_gammas=0;  // for debug
	int n_electrons=0;
	for(int secondary_i=0; secondary_i<n_secondaries_2; ++secondary_i){
		// ---------------
		// Scan for Gammas
		// ---------------
		if(secondary_PDG_code_2.at(secondary_i)==gamma_pdg){
			Log(toolName+" GAMMA!", v_debug,verbosity);
			n_gammas++;
			// we found a gamma! See if it came from neutron capture (G3 process 18)
			bool from_ncapture = (secondary_gen_process.at(secondary_i)==18);
			Log(toolName+" from ncapture="+toString(from_ncapture),v_debug,verbosity);
			//if(not from_ncapture) continue;  // if not, not interested, skip it
			// for now, as a sanity check, checks if its parent was an identified neutron first
			// parent may either be a primary particle or secondary particle
			// sanity check: it should have one or the other, but not both
			int neutron_parent_loc = -1;
			int primary_parent_index = parent_trackid.at(secondary_i);
			int secondary_parent_index = parent_index.at(secondary_i);
			Log(toolName+" primary parent index "+toString(primary_parent_index)
						+" secondary parent index "+toString(secondary_parent_index),v_debug,verbosity);
			// first check if it has a valid SECONDARY parent
			if(secondary_parent_index>0){
				// its parent was a secondary: check if it's in our list of secondary neutrons
				if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){ // -1 as indices are 1-based
					neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
				}
				// else its parent was a secondary, but not one we know
				else if(from_ncapture){
					// if it came from ncapture of a secondary neutron, why don't we know about that neutron?
					Log(toolName+" WARNING, GAMMA FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT INDEX "
							+toString(secondary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			// only fall-back to getting parent PRIMARY if parent secondary index = 0
			// this is because parent primary index is carried over, so daughters of secondaries
			// will have the same primary parent index
			else if(primary_parent_index>0){
				// its parent was a primary: check if it's in our list of primary neutrons
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){  // -1 as indices are 1-based
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				}
				// else its parent was a primary, but not one we know
				else if(from_ncapture){
					// if it came from ncapture of a primary neutron, why don't we know about that neutron?
					Log(toolName+"WARNING, GAMMA FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT INDEX "
							 +toString(primary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			if((neutron_parent_loc<0)&&(from_ncapture)){
				// if we got here, it suggests this gamma is from ncapture
				// but had neither a valid primary or secondary parent index
				Log(toolName+" WARNING, GAMMA FROM NCAPTURE WITH NO PRIMARY OR SECONDARY PARENT INDEX ",
						v_warning,verbosity);
				continue;
			}
			if((neutron_parent_loc>=0)&&(not from_ncapture)){
				// gamma with a parent matched to one of our known neutrons, but not from ncapture?
				// maybe from fast neutron scattering? Leave these for later
				Log(toolName+" WARNING, GAMMA WITH NEUTRON PARENT BUT NOT FROM NCAPTURE",v_warning,verbosity);
				continue;
			}
			// one final check extracts the ones we're after
			if(neutron_parent_loc>=0){
				// neutron from ncapture with known neutron parent! hurray!
				if(neutron_terminfo_unknown.at(neutron_parent_loc)){
					// this is the first daughter gamma we've found for this neutron,
					// so use it's creation info to update the termination info for the neutron
					TLorentzVector gamma_start_pos(secondary_start_vertex_2.at(secondary_i).at(0),
												   secondary_start_vertex_2.at(secondary_i).at(1),
												   secondary_start_vertex_2.at(secondary_i).at(2),
												   secondary_start_time_2.at(secondary_i));
					out_neutron_end_pos.at(neutron_parent_loc) = gamma_start_pos;
					// get momentum of parent at time of creation of secondary (neutron E at capture)
					TVector3 neutron_end_mom(parent_mom_at_sec_creation.at(secondary_i).at(0),
											 parent_mom_at_sec_creation.at(secondary_i).at(1),
											 parent_mom_at_sec_creation.at(secondary_i).at(2));
					double endE = neutron_end_mom.Mag();
					out_neutron_end_energy.at(neutron_parent_loc) = endE;
					out_neutron_end_process.at(neutron_parent_loc) = secondary_gen_process.at(secondary_i);
					neutron_terminfo_unknown.at(neutron_parent_loc) = false;
					
					// double check - for primary neutrons we assume the neutron start pos is
					// the primary event vertex. Compare with the "parent position at creation"
					TVector3 parent_neutron_start_pos(parent_init_pos.at(secondary_i).at(0),
													  parent_init_pos.at(secondary_i).at(1),
													  parent_init_pos.at(secondary_i).at(2));
					if(parent_neutron_start_pos!=out_neutron_start_pos.at(neutron_parent_loc).Vect()){
						std::cout<<"WARNING, NEUTRON START LOC FROM PARENT POSITION AT BIRTH ("
								 <<parent_neutron_start_pos.X()<<", "<<parent_neutron_start_pos.Y()
								 <<", "<<parent_neutron_start_pos.Z()<<") "
								 <<" DIFFERS FROM NEUTRON START LOC FROM NEUTRON ITSELF ("
								 <<out_neutron_start_pos.at(neutron_parent_loc).X()<<", "
								 <<out_neutron_start_pos.at(neutron_parent_loc).Y()<<", "
								 <<out_neutron_start_pos.at(neutron_parent_loc).Z()<<") "
								 <<std::endl;
					}
					
					// do the same with initial momentum via parent_init_mom
					TVector3 parent_neutron_start_mom(parent_init_mom.at(secondary_i).at(0),
													  parent_init_mom.at(secondary_i).at(1),
													  parent_init_mom.at(secondary_i).at(2));
					double nStartE = parent_neutron_start_mom.Mag2() / (2.*neutron_mass);
					if(nStartE!=out_neutron_start_energy.at(neutron_parent_loc)){
						std::cout<<"WARNING, NEUTRON START ENERGY FROM PARENT AT BIRTH ("<<nStartE
								 <<") DIFFERS FROM NEUTRON START ENERGY FROM NEUTRON ITSELF ("
								 <<out_neutron_start_energy.at(neutron_parent_loc)<<") "<<std::endl;
					}
				}
				
				// ok, now add the gamma info
				TVector3 startmom(secondary_start_mom_2.at(secondary_i).at(0),
								  secondary_start_mom_2.at(secondary_i).at(1),
								  secondary_start_mom_2.at(secondary_i).at(2));
				double startE = startmom.Mag();
				double startT = secondary_start_time_2.at(secondary_i);
				out_gamma_energy.at(neutron_parent_loc).push_back(startE);
				out_gamma_time.at(neutron_parent_loc).push_back(startT);
				out_neutron_ndaughters.at(neutron_parent_loc)++;
			}
			// otherwise this is just a plain old gamma. Valid parent, primary or secondary but not both,
			// not from ncapture, and whose parent is not a neutron. Not interested.
		}
		// -----------------------------
		// Scan for Conversion Electrons
		// -----------------------------
		else if(secondary_gen_process.at(secondary_i)==18 && 
			    secondary_PDG_code_2.at(secondary_i)==11){
			Log(toolName+" Conversion Electron!", v_debug,verbosity);
			n_electrons++;
			// for now, as a sanity check, checks if its parent was an identified neutron first
			// parent may either be a primary particle or secondary particle
			// sanity check: it should have one or the other, but not both
			int neutron_parent_loc = -1;
			int primary_parent_index = parent_trackid.at(secondary_i);
			int secondary_parent_index = parent_index.at(secondary_i);
			Log(toolName+" primary parent index "+toString(primary_parent_index)
						+" secondary parent index "+toString(secondary_parent_index),v_debug,verbosity);
			// first check if it has a valid SECONDARY parent
			if(secondary_parent_index>0){
				// its parent was a secondary: check if it's in our list of secondary neutrons
				if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){ // -1 as indices are 1-based
					neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
				}
				// else its parent was a secondary, but not one we know
				else {
					// if it came from ncapture of a secondary neutron, why don't we know about that neutron?
					Log(toolName+" WARNING, CONVERSION ELECTRON FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT"
							+toString(secondary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			// only fall-back to getting parent PRIMARY if parent secondary index = 0
			// this is because parent primary index is carried over, so daughters of secondaries
			// will have the same primary parent index
			else if(primary_parent_index>0){
				// its parent was a primary: check if it's in our list of primary neutrons
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){  // -1 as indices are 1-based
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				}
				// else its parent was a primary, but not one we know
				else {
					// if it came from ncapture of a primary neutron, why don't we know about that neutron?
					Log(toolName+"WARNING, CONVERSION ELECTRON FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT"
							 +toString(primary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			if(neutron_parent_loc<0){
				// if we got here, it suggests this gamma is from ncapture
				// but had neither a valid primary or secondary parent index
				Log(toolName+" WARNING, CONVERSION ELECTRON FROM NCAPTURE WITH NO PRIMARY OR SECONDARY PARENT",
						v_warning,verbosity);
				continue;
			}
			// one final check extracts the ones we're after
			if(neutron_parent_loc>=0){
				// neutron from ncapture with known neutron parent! hurray!
				if(neutron_terminfo_unknown.at(neutron_parent_loc)){
					// this is the first daughter we've found for this neutron,
					// so use it's creation info to update the termination info for the neutron
					TLorentzVector elec_start_pos(secondary_start_vertex_2.at(secondary_i).at(0),
												  secondary_start_vertex_2.at(secondary_i).at(1),
												  secondary_start_vertex_2.at(secondary_i).at(2),
												  secondary_start_time_2.at(secondary_i));
					out_neutron_end_pos.at(neutron_parent_loc) = elec_start_pos;
					// get momentum of parent at time of creation of secondary (neutron E at capture)
					TVector3 neutron_end_mom(parent_mom_at_sec_creation.at(secondary_i).at(0),
											 parent_mom_at_sec_creation.at(secondary_i).at(1),
											 parent_mom_at_sec_creation.at(secondary_i).at(2));
					double endE = neutron_end_mom.Mag();
					out_neutron_end_energy.at(neutron_parent_loc) = endE;
					out_neutron_end_process.at(neutron_parent_loc) = secondary_gen_process.at(secondary_i);
					neutron_terminfo_unknown.at(neutron_parent_loc) = false;
					
					// double check - for primary neutrons we assume the neutron start pos is
					// the primary event vertex. Compare with the "parent position at creation"
					TVector3 parent_neutron_start_pos(parent_init_pos.at(secondary_i).at(0),
													  parent_init_pos.at(secondary_i).at(1),
													  parent_init_pos.at(secondary_i).at(2));
					if(parent_neutron_start_pos!=out_neutron_start_pos.at(neutron_parent_loc).Vect()){
						std::cout<<"WARNING, NEUTRON START LOC FROM PARENT POSITION AT BIRTH ("
								 <<parent_neutron_start_pos.X()<<", "<<parent_neutron_start_pos.Y()
								 <<", "<<parent_neutron_start_pos.Z()<<") "
								 <<" DIFFERS FROM NEUTRON START LOC FROM NEUTRON ITSELF ("
								 <<out_neutron_start_pos.at(neutron_parent_loc).X()<<", "
								 <<out_neutron_start_pos.at(neutron_parent_loc).Y()<<", "
								 <<out_neutron_start_pos.at(neutron_parent_loc).Z()<<") "
								 <<std::endl;
					}
					
					// do the same with initial momentum via parent_init_mom
					TVector3 parent_neutron_start_mom(parent_init_mom.at(secondary_i).at(0),
													  parent_init_mom.at(secondary_i).at(1),
													  parent_init_mom.at(secondary_i).at(2));
					double nStartE = parent_neutron_start_mom.Mag2() / (2.*neutron_mass);
					if(nStartE!=out_neutron_start_energy.at(neutron_parent_loc)){
						std::cout<<"WARNING, NEUTRON START ENERGY FROM PARENT AT BIRTH ("<<nStartE
								 <<") DIFFERS FROM NEUTRON START ENERGY FROM NEUTRON ITSELF ("
								 <<out_neutron_start_energy.at(neutron_parent_loc)<<") "<<std::endl;
					}
				}
				
				// ok, now add the gamma info
				TVector3 startmom(secondary_start_mom_2.at(secondary_i).at(0),
								  secondary_start_mom_2.at(secondary_i).at(1),
								  secondary_start_mom_2.at(secondary_i).at(2));
				double startE = startmom.Mag();
				double startT = secondary_start_time_2.at(secondary_i);
				out_electron_energy.at(neutron_parent_loc).push_back(startE);
				out_electron_time.at(neutron_parent_loc).push_back(startT);
				out_neutron_ndaughters.at(neutron_parent_loc)++;
			}
		}
		// ------------------------
		// Scan for Daughter Nuclei
		// ------------------------
		else if(secondary_gen_process.at(secondary_i)==18 &&
			    secondary_PDG_code_2.at(secondary_i)!=11){
			// if it's not a gamma or internal conversion electron,
			// but came from neutron capture it's the daughter nuclide
			// its parent will also be the captured neutron, same as the decay gamma.
			int neutron_parent_loc=-1;
			int primary_parent_index = parent_trackid.at(secondary_i);
			int secondary_parent_index = parent_index.at(secondary_i);
			if(secondary_parent_index>0){
				if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){ // -1 as indices are 1-based
					neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
				} else {
					// came from capture of a neutron we don't know?
					Log(toolName+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT (NEUTRON) INDEX "
							+toString(secondary_parent_index),v_warning,verbosity);
					continue;
				}
			} else if(primary_parent_index>0){
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){  // -1 as indices are 1-based
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				} else {
					// came from capture of a neutron we don't know?
					Log(toolName+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT (NEUTRON) INDEX "
							+toString(primary_parent_index),v_warning,verbosity);
					continue;
				}
			}
			if(neutron_parent_loc<0){
				Log(toolName+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH NO PARENT (NEUTRON) INDEX ",v_warning,verbosity);
				continue;
			} else {
				// nuclide from ncapture with known neutron parent! hurray!
				// Convert from secondary index (i.e. position in array of all secondaries)
				// into neutron index (i.e. position in our array of neutrons)
				if(out_nuclide_daughter_pdg.at(neutron_parent_loc)>0){
					Log(toolName+" ERROR, FOUND SECOND DAUGHTER NUCLIDE FROM NEUTRON CAPTURE."
						+" FIRST DAUGHTER PDG: "+toString(out_nuclide_daughter_pdg.at(neutron_parent_loc))
						+" ("+PdgToString(out_nuclide_daughter_pdg.at(neutron_parent_loc))+") "
						+" SECOND DAUGHTER PDG: "+toString(secondary_PDG_code_2.at(secondary_i))
						+" ("+PdgToString(secondary_PDG_code_2.at(secondary_i))+")",v_error,verbosity);
					assert(false);
					continue;
				}
				out_nuclide_daughter_pdg.at(neutron_parent_loc) = secondary_PDG_code_2.at(secondary_i);
			}
		} // end if from ncapture
	} // end scan for gammas/daughter nuclides
	Log(toolName+" found "+toString(n_gammas)+" gammas and "+toString(n_electrons)
				+" conversion electrons", v_debug,verbosity);
	
	// record all primaries...? do we need this info?
	for(int primary_i=0; primary_i<n_outgoing_primaries; ++primary_i){
		int primary_pdg_code = G3ParticleCodeToPdg(primary_G3_code.at(primary_i));
		out_primary_pdg.push_back(primary_pdg_code);
		double mom_sq = pow(primary_start_mom.at(primary_i),2.);
		double mass = PdgToMass(primary_pdg_code);
		mass = (mass==0) ? 1 : mass;
		out_primary_energy.push_back(primary_start_mom.at(primary_i));
		TVector3 start_mom_dir(primary_start_mom_dir.at(primary_i).at(0),
							   primary_start_mom_dir.at(primary_i).at(1),
							   primary_start_mom_dir.at(primary_i).at(2));
		TVector3 start_mom = primary_start_mom.at(primary_i)*start_mom_dir.Unit();
		out_primary_start_mom.push_back(start_mom);
		out_primary_start_pos.push_back(primary_vertex_tvector); // not sure about the validity of this
		out_primary_end_pos.push_back(TLorentzVector(0,0,0,0)); // need to get this from a daughter
	}
	
	return 1;
}

int TruthNeutronCaptures::ReadEntryNtuple(long entry_number){
	int bytesread = myTreeReader.GetEntry(entry_number);
	if(bytesread<=0) return bytesread;
	
	int success = 
	// file level
	// simulation version?
	(myTreeReader.GetBranchValue("wlen",water_transparency)) &&  // [cm]
	
	// event meta info
	(myTreeReader.GetBranchValue("nrun",run_number))         &&
	(myTreeReader.GetBranchValue("nsub",subrun_number))      &&
	(myTreeReader.GetBranchValue("nev",event_number))        &&
	(myTreeReader.GetBranchValue("nsube",subevent_number))   &&  // how does this relate to after trigger?
//	(myTreeReader.GetBranchValue("date",event_date))         &&  // [year,month,day]
//	(myTreeReader.GetBranchValue("time",time))               &&  // [hour,minute,second,?]
	
	// event level detector info
	(myTreeReader.GetBranchValue("nhit",N_hit_ID_PMTs))      &&  // "nqisk"
	(myTreeReader.GetBranchValue("potot",total_ID_pes))      &&  // "qismsk"
	(myTreeReader.GetBranchValue("pomax",max_ID_PMT_pes))    &&  // "qimxsk", presumably max # PEs from an ID PMT?
	
	// numnu is 0 even when npar is >3...
//	// neutrino interaction info - first primaries array includes neutrino and target (index 0 and 1)
//	(myTreeReader.GetBranchValue("mode",nu_intx_mode))       &&  // see neut_mode_to_string(mode)
//	(myTreeReader.GetBranchValue("numnu",tot_n_primaries))   &&  // both ingoing and outgoing
//	
//	// following are arrays of size numnu
//	(myTreeReader.GetBranchValue("ipnu",primary_pdg))        &&  // see constants::numnu_code_to_string
//	(myTreeReader.GetBranchValue("pnu",primary_momentum))    &&  // [GeV/c]
	
	// primary event - second primaries array includes more info
	(myTreeReader.GetBranchValue("posv",primary_event_vertex))          &&  // [cm]
//	(myTreeReader.GetBranchValue("wallv",primary_event_dist_from_wall)) &&  // [cm]
	(myTreeReader.GetBranchValue("npar",n_outgoing_primaries))          &&  // should be (tot_n_primaries - 2)?
	
	// following are arrays of size npar
	(myTreeReader.GetBranchValue("ipv",primary_G3_code))                &&  // see constants::g3_to_pdg
	(myTreeReader.GetBranchValue("dirv",primary_start_mom_dir))         &&  // 
	(myTreeReader.GetBranchValue("pmomv",primary_start_mom))            &&  // [units?]
	
//	// secondaries - first secondaries arrays...
//	(myTreeReader.GetBranchValue("npar2",n_secondaries_1))              &&
	// npar2 is 0 even when nscndprt is not???
//	
//	// following are arrays of size npar2
//	(myTreeReader.GetBranchValue("ipv2",secondary_G3_code_1))                &&  // 
//	(myTreeReader.GetBranchValue("posv2",secondary_start_vertex_1))          &&  // [cm?] what about time?
//	(myTreeReader.GetBranchValue("wallv2",secondary_start_dist_from_wall_1)) &&  // [cm?]
//	(myTreeReader.GetBranchValue("pmomv2",secondary_start_mom_1))            &&  // [units?]
//	(myTreeReader.GetBranchValue("iorg",secondary_origin_1))                 &&  // what is "origin"?
	
	// secondaries - second secondaries array...
	(myTreeReader.GetBranchValue("nscndprt",n_secondaries_2))          &&
	
	// following are arrays of size nscndprt
	(myTreeReader.GetBranchValue("iprtscnd",secondary_PDG_code_2))     &&  //
	(myTreeReader.GetBranchValue("vtxscnd",secondary_start_vertex_2))  &&  // [units?]
	(myTreeReader.GetBranchValue("tscnd",secondary_start_time_2))      &&  // [ns]? relative to event start?
	(myTreeReader.GetBranchValue("pscnd",secondary_start_mom_2))       &&  // [units?]
	(myTreeReader.GetBranchValue("lmecscnd",secondary_gen_process))    &&  // constants::G3_process_code_to_string
	(myTreeReader.GetBranchValue("nchilds",secondary_n_daughters))     &&  // 
	(myTreeReader.GetBranchValue("iprntidx",parent_index))             &&  // if >0, 1-based index in this array
//	(myTreeReader.GetBranchValue("ichildidx",secondary_first_daugher_index)) &&  // if >0, 1-based index in this
	
	// further parentage information - still arrays of size nscndprt. Useful?
//	(myTreeReader.GetBranchValue("iprntprt",parent_G3_code))           &&  // or is it a PDG code?
	(myTreeReader.GetBranchValue("pprnt",parent_mom_at_sec_creation))  &&  // use w/γ to get n energy @ capture
	(myTreeReader.GetBranchValue("vtxprnt",parent_init_pos))           &&  // [cm?] parent pos @ birth
	(myTreeReader.GetBranchValue("pprntinit",parent_init_mom))         &&  // [MeV?] parent mom @ birth
//	(myTreeReader.GetBranchValue("itrkscnd",parent_G3_trackid))        &&  // how do we use this?
//	(myTreeReader.GetBranchValue("istakscnd",parent_G3_stack_trackid)) &&  // how do we use this?
	(myTreeReader.GetBranchValue("iprnttrk",parent_trackid));        //&&  // relates secondaries to primaries
	// NOTE this is carried over to daughters of secondaries, so only use as parent if iprntidx==0
//	(myTreeReader.GetBranchValue("iorgprt",parent_track_pid_code))     &&  // i'm so confused
	
	// XXX for efficiency, add all used branches to DisableUnusedBranches XXX
	
	return success;
}

int TruthNeutronCaptures::DisableUnusedBranches(){
	std::vector<std::string> used_branches{
		"wlen",
		"nrun",
		"nsub",
		"nev",
		"nsube",
//		"date",
//		"time",
		"nhit",
		"potot",
		"pomax",
//		"mode",
//		"numnu",
//		"ipnu",
//		"pnu",
		"posv",
//		"wallv",
		"npar",
		"ipv",
		"dirv",
		"pmomv",
//		"npar2",
//		"ipv2",
//		"posv2",
//		"wallv2",
//		"pmomv2",
//		"iorg",
		"nscndprt",
		"iprtscnd",
		"vtxscnd",
		"tscnd",
		"pscnd",
		"lmecscnd",
		"nchilds",
		"iprntidx",
//		"ichildidx",
//		"iprntprt",
		"pprnt",
		"vtxprnt",
		"pprntinit",
//		"itrkscnd",
		"istakscnd",
		"iprnttrk"  //,
//		"iorgprt",
	};
	
	return myTreeReader.OnlyEnableBranches(used_branches);
}

int TruthNeutronCaptures::CreateOutputFile(std::string filename){
	// create the output ROOT file and TTree for writing
	// =================================================
	outfile = new TFile(filename.c_str(), "RECREATE");
	outtree = new TTree("eventtree", "Events with Neutron Captures");
	
	// create branches
	// ---------------
	// file level
	outtree->Branch("filename",&out_filename);
	outtree->Branch("skdetsim_version",&out_skdetsim_version);    // where?
	outtree->Branch("tba_table_version",&out_tba_table_version);  // where?
	outtree->Branch("water_transparency",&out_water_transparency);
	
	// event level
	outtree->Branch("run_number",&out_run_number);
//	outtree->Branch("subrun_number",&out_subrun_number);
	outtree->Branch("entry_number",&out_entry_number);
	outtree->Branch("subevent_num",&out_subevent_number);
	
	// primary particle
	outtree->Branch("primary_pdg",&out_primary_pdg,32000,0);
	outtree->Branch("primary_energy",&out_primary_energy,32000,0);
	outtree->Branch("primary_start_mom",&out_primary_start_mom,32000,0);
	outtree->Branch("primary_start_pos",&out_primary_start_pos,32000,0);
	outtree->Branch("primary_end_pos",&out_primary_end_pos,32000,0);
	
	// parent nuclide - one for each neutron
	outtree->Branch("nuclide_parent_pdg",&out_nuclide_parent_pdg,32000,0);
	outtree->Branch("nuclide_daughter_pdg",&out_nuclide_daughter_pdg,32000,0);
	
	// neutron
	outtree->Branch("neutron_start_pos",&out_neutron_start_pos,32000,0);
	outtree->Branch("neutron_end_pos",&out_neutron_end_pos,32000,0);
	outtree->Branch("neutron_start_energy",&out_neutron_start_energy,32000,0);
	outtree->Branch("neutron_end_energy",&out_neutron_end_energy,32000,0);
	outtree->Branch("neutron_end_process",&out_neutron_end_process,32000,0);
//	outtree->Branch("neutron_n_gammas",&out_neutron_n_gammas,32000,0);
//	outtree->Branch("neutron_n_electrons",&out_neutron_n_electrons,32000,0);
//	outtree->Branch("neutron_n_daughters",&out_neutron_ndaughters,32000,0);
	
	
	// gamma
	outtree->Branch("gamma_energy",&out_gamma_energy,32000,0);
	outtree->Branch("gamma_time",&out_gamma_time,32000,0);
	// use outtree->Draw("Length$(gamma_energy[])"); to draw gamma multiplicity (equivalent to neutron_n_daughters)
	// use outtree->Draw("Sum$(gamma_energy[])");    to draw total gamma energy from a capture
	
	// electron
	outtree->Branch("electron_energy",&out_electron_energy,32000,0);
	outtree->Branch("electron_time",&out_electron_time,32000,0);
	
	return 1;
}

void TruthNeutronCaptures::ClearOutputTreeBranches(){
	// clear any vector branches
	
	out_primary_pdg.clear();
	out_primary_energy.clear();
	out_primary_start_mom.clear();
	out_primary_start_pos.clear();
	out_primary_end_pos.clear();
	
	// parent nuclide
	out_nuclide_parent_pdg.clear();
	out_nuclide_daughter_pdg.clear();
	
	// neutrons
	out_neutron_start_pos.clear();
	out_neutron_end_pos.clear();
	out_neutron_start_energy.clear();
	out_neutron_end_energy.clear();
	out_neutron_end_process.clear();
	out_neutron_ndaughters.clear();
	
	// gammas
	out_gamma_energy.clear();
	out_gamma_time.clear();
	
	// electrons
	out_electron_energy.clear();
	out_electron_time.clear();
	
	return;
}

void TruthNeutronCaptures::PrintBranches(){
	std::cout<<"==========================================================="<<std::endl;
	std::cout<<"PRINTING EVENT"<<std::endl;
	std::cout<<"filename: "<<out_filename<<std::endl;
	std::cout<<"water transparency: "<<out_water_transparency<<std::endl;
	
//	std::cout<<"skdetsim_version: "<<out_skdetsim_version<<std::endl
//			 <<"tba_table_version: "<<out_tba_table_version<<std::endl
//			 <<"neutron_process_map: {";
//	for(std::map<int,std::string>::const_iterator aprocess=neutron_process_map.begin();
//	    aprocess!=neutron_process_map.end(); ++aprocess){
//			std::cout<<"["<<aprocess->first<<"="<<aprocess->second<<"], ";
//	}
//	std::cout<<"\b\b}"<<std::endl;
	std::cout<<"entry_number: "<<out_entry_number<<std::endl
			 <<"subevent_number:" <<out_subevent_number<<std::endl
			 <<"num primaries: "<<out_primary_pdg.size()<<std::endl;
	if(out_primary_pdg.size()){
		std::cout<<"primary vertex:"
				 <<" ("<<out_primary_start_pos.at(0).X()
				 <<", "<<out_primary_start_pos.at(0).Y()
				 <<", "<<out_primary_start_pos.at(0).Z()<<")"<<std::endl;
	}
	
	// print primaries
	for(int primary_i=0; primary_i<out_primary_pdg.size(); ++primary_i){
		std::cout<<"\tprimary ("<<primary_i<<"): "<<std::endl
				 <<"\t\tprimary pdg: "<<out_primary_pdg.at(primary_i)<<std::endl
				 <<"\t\tprimary energy: "<<out_primary_energy.at(primary_i)<<std::endl
				 <<"\t\tprimary momentum:"
				 <<" ("<<out_primary_start_mom.at(primary_i).X()
				 <<", "<<out_primary_start_mom.at(primary_i).Y()
				 <<", "<<out_primary_start_mom.at(primary_i).Z()<<")"<<std::endl;
//		std::cout<<"\t\tprimary start pos:"
//				 <<" ("<<out_primary_start_pos.at(primary_i).X()
//				 <<", "<<out_primary_start_pos.at(primary_i).Y()
//				 <<", "<<out_primary_start_pos.at(primary_i).Z()<<")"<<std::endl
//				 <<"\t\tprimary end pos:"
//				 <<" ("<<out_primary_end_pos.at(primary_i).X()
//				 <<", "<<out_primary_end_pos.at(primary_i).Y()
//				 <<", "<<out_primary_end_pos.at(primary_i).Z()<<")"<<std::endl;
	}
	
	int total_neutrons=0;
	int total_gammas=0;
	int total_electrons=0;
	// print neutron captures
	std::cout<<"num neutrons: "<<out_neutron_start_energy.size()<<std::endl;
	for(int neutron_i=0; neutron_i<out_neutron_start_energy.size(); ++neutron_i){
		std::cout<<"\tneutron "<<neutron_i<<": "<<std::endl
				 <<"\t\tneutron start energy: "
				 <<out_neutron_start_energy.at(neutron_i)<<std::endl
				 <<"\t\tneutron start pos:"
				 <<" ("<<out_neutron_start_pos.at(neutron_i).X()
				 <<", "<<out_neutron_start_pos.at(neutron_i).Y()
				 <<", "<<out_neutron_start_pos.at(neutron_i).Z()<<")"<<std::endl
				 <<"\t\tneutron end pos:"
				 <<" ("<<out_neutron_end_pos.at(neutron_i).X()
				 <<", "<<out_neutron_end_pos.at(neutron_i).Y()
				 <<", "<<out_neutron_end_pos.at(neutron_i).Z()<<")"<<std::endl
				 <<"\t\tneutron end energy: "
				 <<out_neutron_end_energy.at(neutron_i)<<std::endl
				 <<"\t\tneutron end process: "<<out_neutron_end_process.at(neutron_i)
				 <<std::endl
				 
				 // each neutron has one nuclide
				 <<"\t\tparent nuclide pdg: "<<out_nuclide_parent_pdg.at(neutron_i)<<std::endl
				 <<"\t\tdaughter nuclide pdg: "<<out_nuclide_daughter_pdg.at(neutron_i)<<std::endl
				 
				 // each neutron may have multiple decay gammas
				 <<"\t\tnum gammas from this neutron: "
				 <<out_gamma_energy.at(neutron_i).size()<<std::endl;
		for(int gamma_i=0; gamma_i<out_gamma_energy.at(neutron_i).size(); ++gamma_i){
			std::cout<<"\t\tgamma "<<gamma_i<<": "<<std::endl
					 <<"\t\t\tgamma energy: "<<out_gamma_energy.at(neutron_i).at(gamma_i)
					 <<std::endl
					 <<"\t\t\tgamma time: "<<out_gamma_time.at(neutron_i).at(gamma_i)
					 <<std::endl;
			total_gammas++;
		}
		std::cout<<"\t\tnum internal conversion electrons from this neutron: "
				 <<out_electron_energy.at(neutron_i).size()<<std::endl;
		for(int electron_i=0; electron_i<out_electron_energy.at(neutron_i).size(); ++electron_i){
			std::cout<<"\t\tgamma "<<electron_i<<": "<<std::endl
					 <<"\t\t\tgamma energy [MeV]: "<<out_electron_energy.at(neutron_i).at(electron_i)
					 <<std::endl
					 <<"\t\t\tgamma time [ns]: "<<out_electron_time.at(neutron_i).at(electron_i)
					 <<std::endl;
			total_electrons++;
		}
	}
	std::cout<<"total gammas in event: "<<total_gammas<<std::endl;
	std::cout<<"total conversion electrons in event: "<<total_electrons<<std::endl;
	std::cout<<"==========================================================="<<std::endl;
}

int TruthNeutronCaptures::WriteTree(){
	Log(toolName+" writing TTree",v_debug,verbosity);
	outfile->cd();
	// TObject::Write returns the total number of bytes written to the file.
	// It returns 0 if the object cannot be written.
	int bytes = outtree->Write("",TObject::kOverwrite);
	if(bytes<=0){
		Log(toolName+" Error writing TTree!",v_error,verbosity);
	} else if(verbosity>2){
		Log(toolName+ " Wrote "+toString(get_ok)+" bytes",v_debug,verbosity);
	}
	return bytes;
};

void TruthNeutronCaptures::CloseFile(){
	outtree->ResetBranchAddresses();
	// XXX note that while these do persist in the tree:
	// 1. they can't be used in a loop, since you can't use SetBranchAddress on an Alias.
	//    this is probably less of a concern as they can easily be calculated within a loop.
	// 2. they will not show up in TTree::Print(), so you need to know they're there...!
	// Would it be better to store the redundant information? Is there a better way?
	outtree->SetAlias("neutron_travel_dist","sqrt(pow(neutron_start_pos.X()-neutron_end_pos.X(),2)+pow(neutron_start_pos.Y()-neutron_end_pos.Y(),2)+pow(neutron_start_pos.Z()-neutron_end_pos.Z(),2))");
	outtree->SetAlias("neutron_travel_time","neutron_end_pos.T()-neutron_start_pos.T()");
	outtree->SetAlias("neutron_n_daughters","Length$(gamma_energy[])");
	outtree->SetAlias("neutron_tot_gammaE","Sum$(gamma_energy[])");
	outfile->Write("*",TObject::kOverwrite);
	outfile->Close();
	delete outfile;
	outfile = nullptr;
};


