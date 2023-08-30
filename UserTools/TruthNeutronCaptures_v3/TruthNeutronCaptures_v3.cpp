/* vim:set noexpandtab tabstop=4 wrap */
#include "TruthNeutronCaptures_v3.h"
#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <algorithm> // std::find

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"

TruthNeutronCaptures_v3::TruthNeutronCaptures_v3():Tool(){}

bool TruthNeutronCaptures_v3::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);            // how verbose to be
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEvents",MAX_EVENTS);           // terminate after processing at most this many events
	m_variables.Get("writeFrequency",WRITE_FREQUENCY); // how many events to TTree::Fill between TTree::Writes
	std::string treeReaderName;
	m_variables.Get("treeReaderName",treeReaderName);  // reader name for input data
	
	// Get input TreeReader
	// --------------------
	if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",0,0);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	// create the output TFile and TTree
	// ---------------------------------
	CreateOutputFile(outputFile);
	
	return true;
}

bool TruthNeutronCaptures_v3::Execute(){
	
	// enclose processing in a try-catch loop: we must catch any thrown errors
	// to prevent the tool crashing out, so that we always pre-load the next entry.
	// otherwise we'll just try to re-process it next Execute() and crash again!
	try {
		Log(m_unique_name+" processing entry "+toString(entry_number),v_debug,m_verbose);
		
		// clear output vectors so we don't carry anything over
		Log(m_unique_name+" clearing output vectors",v_debug,m_verbose);
		ClearOutputTreeBranches();
		
		// Copy over directly transferred variables
		Log(m_unique_name+" copying output variables",v_debug,m_verbose);
		CopyVariables();
		
		// Calculate derived variables
		Log(m_unique_name+" calculating output variables",v_debug,m_verbose);
		CalculateVariables();
		
		// print the current event
		if(m_verbose>1) PrintBranches();
		
	} catch (...){
		Log(m_unique_name+" ERROR PROCESSING ENTRY "+toString(entry_number),v_error,m_verbose);
	}
	
	// Fill the output tree
	Log(m_unique_name+" filling output TTree entry",v_debug,m_verbose);
	outtree->Fill();
	
	// update the output file so we don't lose everything if we crash
	if((entry_number%WRITE_FREQUENCY)==0) WriteTree();
	
	// stop at user-defined limit to the number of events to process
	++entry_number;
	if((MAX_EVENTS>0)&&(entry_number>=MAX_EVENTS)){
		Log(m_unique_name+" reached MAX_EVENTS, setting StopLoop",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
	} else {
		// Pre-Load next input entry so we can stop the toolchain
		// if we're about to run off the end of the tree or encounter a read error
		get_ok = ReadEntryNtuple(entry_number);
		if(get_ok<1&&get_ok>-3){
			m_data->vars.Set("StopLoop",1);
			Log(m_unique_name+" Hit end of input file, stopping loop",v_warning,m_verbose);
		}
		else if(get_ok==-10){
			Log(m_unique_name+" Error during AutoClear while loading next input ntuple entry!",v_error,m_verbose);
			return false;
		}
		else if(get_ok<0){
			Log(m_unique_name+" IO error loading next input ntuple entry!",v_error,m_verbose);
			return false;
		}
	}
	
	return true;
}


bool TruthNeutronCaptures_v3::Finalise(){
	
	// ensure everything is written to the output file
	// -----------------------------------------------
	get_ok = WriteTree();
	if(not get_ok){
		Log(m_unique_name+" Error writing output TTree!",v_error,m_verbose);
	}
	
	// Close and delete the file handle
	// --------------------------------
	CloseFile();
	
	return true;
}

void TruthNeutronCaptures_v3::CopyVariables(){
	// copy over variables from the input tree
	// ---------------------------------------
	// those we want to keep in the output tree without modification
	
	out_filename = myTreeReader->GetTree()->GetCurrentFile()->GetName();
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

int TruthNeutronCaptures_v3::CalculateVariables(){
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
	Log(m_unique_name+" event had "+toString(n_outgoing_primaries)+" primary particles",v_debug,m_verbose);
	for(int primary_i=0; primary_i<n_outgoing_primaries; ++primary_i){
		Log(m_unique_name+" primary "+toString(primary_i)+" had PDG code "+toString((int)primary_PDG_code.at(primary_i))
				+" ("+PdgToString(primary_PDG_code.at(primary_i))+")",v_debug,m_verbose);
		if(primary_PDG_code.at(primary_i)==neutron_pdg){
			// we found a neutron!
			Log(m_unique_name+" NEUTRON!", v_debug,m_verbose);
			primary_n_ind_to_loc.emplace(primary_i,out_neutron_start_energy.size());
			// note neutron info
			out_neutron_start_pos.push_back(primary_vertex_tvector);
			TVector3 start_mom( primary_start_mom.at(primary_i).at(0),
								primary_start_mom.at(primary_i).at(1),
								primary_start_mom.at(primary_i).at(2));
			double startE = sqrt(start_mom.Mag2() + pow(neutron_mass,2.)) - neutron_mass;
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
	Log(m_unique_name+" found "+toString(out_neutron_start_energy.size())+" primary neutrons", v_debug,m_verbose);
	
	// ===========================
	// SCAN FOR SECONDARY NEUTRONS
	// ===========================
	Log(m_unique_name+" event had "+toString(n_secondaries_2)+" secondary particles",v_debug,m_verbose);
	for(int secondary_i=0; secondary_i<n_secondaries_2; ++secondary_i){
		Log(m_unique_name+" secondary "+toString(secondary_i)+" had pdg "+toString(secondary_PDG_code_2.at(secondary_i))
				+" ("+PdgToString(secondary_PDG_code_2.at(secondary_i))+")",v_debug,m_verbose);
		if(secondary_PDG_code_2.at(secondary_i)==neutron_pdg){
			// we found a neutron!
			Log(m_unique_name+" NEUTRON! at "+toString(secondary_i), v_debug,m_verbose);
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
			double startE = sqrt(startmom.Mag2() + pow(neutron_mass,2.)) - neutron_mass;
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
			out_nuclide_daughter_pdg.push_back(-1);
			out_nuclide_parent_pdg.push_back(-1);
		}

		else if(secondary_gen_process.at(secondary_i)==18 &&
			    secondary_PDG_code_2.at(secondary_i)!=11  &&
			    secondary_PDG_code_2.at(secondary_i)!=22 ){
			std::cout<<"Found non-gamma non-electron product from ncapture!"<<std::endl;
			std::cout<<"PDG: "+toString(secondary_PDG_code_2.at(secondary_i))
					 <<" ("+PdgToString(secondary_PDG_code_2.at(secondary_i))+") "
					 <<" created at ("<<secondary_start_vertex_2.at(secondary_i).at(0)
					 <<", "<<secondary_start_vertex_2.at(secondary_i).at(1)
					 <<", "<<secondary_start_vertex_2.at(secondary_i).at(2)
					 <<", "<<secondary_start_time_2.at(secondary_i)<<")"
					 <<", terminated at (";
			// ah, but to get termination info we need to find a daughter from its termination process aughh
			auto it = std::find(parent_index.begin(),parent_index.end(),secondary_i);
			if(it!=parent_index.end()){
				int daughter_index = std::distance(parent_index.begin(),it);
				std::cout<<secondary_start_vertex_2.at(daughter_index).at(0)
				         <<", "<<secondary_start_vertex_2.at(daughter_index).at(1)
				         <<", "<<secondary_start_vertex_2.at(daughter_index).at(2)
				         <<", "<<secondary_start_time_2.at(daughter_index)<<")";
			} else {
				std::cout<<"?,?,?,?)";
			}
			std::cout<<std::endl;
		}

	}
	Log(m_unique_name+" found "+toString(out_neutron_start_energy.size())+" primary+secondary neutrons", v_debug,m_verbose);
	
	DumpMCInfo();
	
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
			Log(m_unique_name+" GAMMA!", v_debug,m_verbose);
			n_gammas++;
			// we found a gamma! See if it came from neutron capture (G3 process 18)
			bool from_ncapture = (secondary_gen_process.at(secondary_i)==18);
			Log(m_unique_name+" from ncapture="+toString(from_ncapture),v_debug,m_verbose);
			//if(not from_ncapture) continue;  // if not, not interested, skip it
			// for now, as a sanity check, checks if its parent was an identified neutron first
			// parent may either be a primary particle or secondary particle
			// sanity check: it should have one or the other, but not both
			int neutron_parent_loc = -1;
			// secondary parent indexes should be fixed when PR is merged
			int primary_parent_index = parent_trackid.at(secondary_i); // not yet implemented
			int secondary_parent_index = parent_index.at(secondary_i);
			Log(m_unique_name+" primary parent index "+toString(primary_parent_index)
						+" secondary parent index "+toString(secondary_parent_index),v_debug,m_verbose);
			// first check if it has a valid SECONDARY parent
			if(std::abs(secondary_parent_index)>0){
				if(secondary_parent_index>0){
					// its parent was a secondary: check if it's in our list of secondary neutrons
					if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){
						neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
					}
					// else its parent was a secondary, but not one we know
					else if(from_ncapture){
						// if it came from ncapture of a secondary neutron, why don't we know about that neutron?
						Log(m_unique_name+" WARNING, GAMMA FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT INDEX "
								+toString(secondary_parent_index),v_warning,m_verbose);
						continue;
					}
				} else {
					// its parent was a primary: check if it's in our list of primary neutrons
					secondary_parent_index = std::abs(secondary_parent_index);
					if(primary_n_ind_to_loc.count(secondary_parent_index-1)){
						neutron_parent_loc = primary_n_ind_to_loc.at(secondary_parent_index-1);
					}
					// else its parent was a primary, but not one we know
					else if(from_ncapture){
						// if it came from ncapture of a primary neutron, why don't we know about that neutron?
						Log(m_unique_name+" WARNING, GAMMA FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT INDEX "
								+toString(secondary_parent_index),v_warning,m_verbose);
						continue;
					}
				}
			}
			// no known secondary parent: fall-back to getting parent PRIMARY
			// this is because parent primary index is carried over, so daughters of secondaries
			// will have the same primary parent index
			else if(primary_parent_index>0){
				// its parent was a primary: check if it's in our list of primary neutrons
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				}
				// else its parent was a primary, but not one we know
				else if(from_ncapture){
					// if it came from ncapture of a primary neutron, why don't we know about that neutron?
					Log(m_unique_name+" WARNING, GAMMA FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT INDEX "
							 +toString(primary_parent_index),v_warning,m_verbose);
					continue;
				}
			}
			if((neutron_parent_loc<0)&&(from_ncapture)){
				// if we got here, it suggests this gamma is from ncapture
				// but had neither a valid primary or secondary parent index
				Log(m_unique_name+" WARNING, GAMMA FROM NCAPTURE WITH NO PRIMARY OR SECONDARY PARENT INDEX ",
						v_warning,m_verbose);
				continue;
			}
			if((neutron_parent_loc>=0)&&(not from_ncapture)){
				// gamma with a parent matched to one of our known neutrons, but not from ncapture?
				// maybe from fast neutron scattering? Leave these for later
				Log(m_unique_name+" WARNING, GAMMA WITH NEUTRON PARENT BUT NOT FROM NCAPTURE",v_warning,m_verbose);
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
					double endE = sqrt(neutron_end_mom.Mag2() + pow(neutron_mass,2.)) - neutron_mass;
					out_neutron_end_energy.at(neutron_parent_loc) = endE;
					out_neutron_end_process.at(neutron_parent_loc) = secondary_gen_process.at(secondary_i);
					neutron_terminfo_unknown.at(neutron_parent_loc) = false;
					
					// double check - for primary neutrons we assume the neutron start pos is
					// the primary event vertex. Compare with the "parent position at creation"
					// XXX HACK: /10 for agreeable units... but which is correct?
					TVector3 parent_neutron_start_pos(parent_init_pos.at(secondary_i).at(0)/10.,
													  parent_init_pos.at(secondary_i).at(1)/10.,
													  parent_init_pos.at(secondary_i).at(2)/10.);
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
					double nStartE = sqrt(parent_neutron_start_mom.Mag2() + pow(neutron_mass,2.)) - neutron_mass;
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
			Log(m_unique_name+" Conversion Electron!", v_debug,m_verbose);
			n_electrons++;
			// for now, as a sanity check, checks if its parent was an identified neutron first
			// parent may either be a primary particle or secondary particle
			// sanity check: it should have one or the other, but not both
			int neutron_parent_loc = -1;
			// secondary parent indexes should be fixed when SKG4 PR is merged
			int primary_parent_index = parent_trackid.at(secondary_i); // not yet implemented
			int secondary_parent_index = parent_index.at(secondary_i);
			Log(m_unique_name+" primary parent index "+toString(primary_parent_index)
						+" secondary parent index "+toString(secondary_parent_index),v_debug,m_verbose);
			// first check if it has a valid SECONDARY parent
			if(std::abs(secondary_parent_index)>0){
				if(secondary_parent_index>0){
					// its parent was a secondary: check if it's in our list of secondary neutrons
					if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){
						neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
					}
					// else its parent was a secondary, but not one we know
					else {
						// if it came from ncapture of a secondary neutron, why don't we know about that neutron?
						Log(m_unique_name+" WARNING, CONVERSION ELECTRON FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT"
								+toString(secondary_parent_index),v_warning,m_verbose);
						continue;
					}
				} else {
					// its parent was a primary: check if it's in our list of primary neutrons
					secondary_parent_index = std::abs(secondary_parent_index);
					if(primary_n_ind_to_loc.count(secondary_parent_index-1)){
						neutron_parent_loc = primary_n_ind_to_loc.at(secondary_parent_index-1);
					}
					// else its parent was a secondary, but not one we know
					else {
						// if it came from ncapture of a secondary neutron, why don't we know about that neutron?
						Log(m_unique_name+" WARNING, CONVERSION ELECTRON FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT"
								+toString(secondary_parent_index),v_warning,m_verbose);
						continue;
					}
				}
			}
			// only fall-back to getting parent PRIMARY if parent secondary index = 0
			// this is because parent primary index is carried over, so daughters of secondaries
			// will have the same primary parent index
			else if(primary_parent_index>0){
				// its parent was a primary: check if it's in our list of primary neutrons
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				}
				// else its parent was a primary, but not one we know
				else {
					// if it came from ncapture of a primary neutron, why don't we know about that neutron?
					Log(m_unique_name+"WARNING, CONVERSION ELECTRON FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT"
							 +toString(primary_parent_index),v_warning,m_verbose);
					continue;
				}
			}
			if(neutron_parent_loc<0){
				// if we got here, it suggests this gamma is from ncapture
				// but had neither a valid primary or secondary parent index
				Log(m_unique_name+" WARNING, CONVERSION ELECTRON FROM NCAPTURE WITH NO PRIMARY OR SECONDARY PARENT",
						v_warning,m_verbose);
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
					double endE = sqrt(neutron_end_mom.Mag2() + pow(neutron_mass,2.)) - neutron_mass;
					out_neutron_end_energy.at(neutron_parent_loc) = endE;
					out_neutron_end_process.at(neutron_parent_loc) = secondary_gen_process.at(secondary_i);
					neutron_terminfo_unknown.at(neutron_parent_loc) = false;
					
					// double check - for primary neutrons we assume the neutron start pos is
					// the primary event vertex. Compare with the "parent position at creation"
					// XXX HACK: /10 for agreeable units... but which is correct?
					TVector3 parent_neutron_start_pos(parent_init_pos.at(secondary_i).at(0)/10.,
													  parent_init_pos.at(secondary_i).at(1)/10.,
													  parent_init_pos.at(secondary_i).at(2)/10.);
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
					double nStartE = sqrt(parent_neutron_start_mom.Mag2() + pow(neutron_mass,2.)) - neutron_mass;
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
			// secondary parent indexes should be fixed when SKG4 PR is merged
			int primary_parent_index = parent_trackid.at(secondary_i); // not yet implemented
			int secondary_parent_index = parent_index.at(secondary_i);
			if(std::abs(secondary_parent_index)>0){
				if(secondary_parent_index>0){
					// secondary parent
					if(secondary_n_ind_to_loc.count(secondary_parent_index-1)){
						neutron_parent_loc = secondary_n_ind_to_loc.at(secondary_parent_index-1);
					} else {
						// came from capture of a neutron we don't know?
						Log(m_unique_name+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
							+" FROM NCAPTURE WITH UNKNOWN SECONDARY PARENT (NEUTRON) INDEX "
								+toString(secondary_parent_index),v_warning,m_verbose);
						continue;
					}
				} else {
					// primary parent
					secondary_parent_index = std::abs(secondary_parent_index);
					if(primary_n_ind_to_loc.count(secondary_parent_index-1)){
						neutron_parent_loc = primary_n_ind_to_loc.at(secondary_parent_index-1);
					} else {
						// came from capture of a neutron we don't know?
						Log(m_unique_name+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
							+" FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT (NEUTRON) INDEX "
								+toString(secondary_parent_index),v_warning,m_verbose);
						continue;
					}
				}
			} else if(primary_parent_index>0){
				if(primary_n_ind_to_loc.count(primary_parent_index-1)){
					neutron_parent_loc = primary_n_ind_to_loc.at(primary_parent_index-1);
				} else {
					// came from capture of a neutron we don't know?
					Log(m_unique_name+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH UNKNOWN PRIMARY PARENT (NEUTRON) INDEX "
							+toString(primary_parent_index),v_warning,m_verbose);
					continue;
				}
			}
			if(neutron_parent_loc<0){
				Log(m_unique_name+" WARNING, "+PdgToString(secondary_PDG_code_2.at(secondary_i))
						+" FROM NCAPTURE WITH NO PARENT (NEUTRON) INDEX ",v_warning,m_verbose);
				continue;
			} else {
				// nuclide from ncapture with known neutron parent! hurray!
				// Convert from secondary index (i.e. position in array of all secondaries)
				// into neutron index (i.e. position in our array of neutrons)
				if(out_nuclide_daughter_pdg.at(neutron_parent_loc)>0){
					Log(m_unique_name+" ERROR, FOUND SECOND DAUGHTER NUCLIDE FROM NEUTRON CAPTURE."
						+" FIRST DAUGHTER PDG: "+toString(out_nuclide_daughter_pdg.at(neutron_parent_loc))
						+" ("+PdgToString(out_nuclide_daughter_pdg.at(neutron_parent_loc))+") "
						+" SECOND DAUGHTER PDG: "+toString(secondary_PDG_code_2.at(secondary_i))
						+" ("+PdgToString(secondary_PDG_code_2.at(secondary_i))+")",v_error,m_verbose);
					assert(false);
					continue;
				}
				out_nuclide_daughter_pdg.at(neutron_parent_loc) = secondary_PDG_code_2.at(secondary_i);
			}
		} // end if from ncapture
	} // end scan for gammas/daughter nuclides
	Log(m_unique_name+" found "+toString(n_gammas)+" gammas and "+toString(n_electrons)
				+" conversion electrons", v_debug,m_verbose);
	
	// record all primaries...? do we need this info?
	for(int primary_i=0; primary_i<n_outgoing_primaries; ++primary_i){
		out_primary_pdg.push_back(primary_PDG_code.at(primary_i));
		TVector3 start_mom( primary_start_mom.at(primary_i).at(0),
							primary_start_mom.at(primary_i).at(1),
							primary_start_mom.at(primary_i).at(2));
		double primary_mass = PdgToMass(primary_PDG_code.at(primary_i));
		double startE = sqrt(start_mom.Mag2() + pow(primary_mass,2.)) - primary_mass;
		out_primary_energy.push_back(startE);
		out_primary_start_mom.push_back(start_mom);
		out_primary_start_pos.push_back(primary_vertex_tvector); // not sure about the validity of this
		out_primary_end_pos.push_back(TLorentzVector(0,0,0,0)); // need to get this from a daughter
	}
	
	return 1;
}

// for comparison of our results with those extracted by the ReadMCInfo Tool.
int TruthNeutronCaptures_v3::DumpMCInfo(){
	m_data->eventPrimaries.DumpAllElements();
	m_data->eventSecondaries.DumpAllElements();
	return 0;
}

int TruthNeutronCaptures_v3::ReadEntryNtuple(long entry_number){
	int bytesread = myTreeReader->GetEntry(entry_number);
	if(bytesread<=0) return bytesread;
	
	int success = (myTreeReader->GetBranchValue("SECONDARY",sec_info)) &&
				  (myTreeReader->GetBranchValue("MC",mc_info));
	
	// Print method should have been const-qualified. Hack around it.
	//MCInfo* mci = const_cast<MCInfo*>(mc_info);
	//mci->Print();
	
	// file level - get these from MCInfo / other SKROOT classes
	// simulation version?
	// water transparency?
	
	// event meta info
	// is more of this info in the mcninfo[] array? Or those other acronym members of mcinfo?
	run_number = mc_info->mcrun;
	// subrun_number
	// event_number
	// subevent_number
	
	// event level detector info
	// N_hit_ID_PMTs
	// total_ID_pes
	// max_ID_PMT_pes
	
	// primary event - Harada-san does not store the primary arrays from the ATMPD-format files
	// (npar, posv, ipv, pmomv), but this information should be available from the MCInfo class.
	// in fact, we have more information here (i think) as we have all the primary vertices,
	// whereas the ATMPD format only seems to store one primary vertex (posv).
	n_outgoing_primaries = mc_info->nvc;    //    (npar)
	// since the tool only supports one primary vertex, take the first for now
	primary_event_vertex = basic_array<float>(intptr_t(mc_info->pvtxvc[0]),3); // (posv)
	
	// following are arrays of size npar
	primary_PDG_code = basic_array<int*>(intptr_t(mc_info->ipvc),n_outgoing_primaries);  // (ipv)
	// MCInfo stores the primary momentum vector, not a separate magnitude and unit direction.
	// the tool has been modified to account for it, hence differences compared to TruthNeutronCaptures.cc
	primary_start_mom = basic_array<float(*)[3]>(intptr_t(mc_info->pvc),n_outgoing_primaries);    // (pmomv)
	
	// secondaries - second secondaries array...
	n_secondaries_2 = sec_info->nscndprt;
	
	// following are arrays of size nscndprt
	secondary_PDG_code_2 = basic_array<int*>(intptr_t(sec_info->iprtscnd),n_secondaries_2);
	secondary_start_vertex_2 = basic_array<float(*)[3]>(intptr_t(sec_info->vtxscnd),n_secondaries_2);
	secondary_start_time_2 = basic_array<float*>(intptr_t(sec_info->tscnd),n_secondaries_2);
	secondary_start_mom_2 = basic_array<float(*)[3]>(intptr_t(sec_info->pscnd),n_secondaries_2);
	secondary_gen_process = basic_array<int*>(intptr_t(sec_info->lmecscnd),n_secondaries_2);
	secondary_n_daughters = basic_array<int*>(intptr_t(sec_info->nchilds),n_secondaries_2);
	parent_index = basic_array<int*>(intptr_t(sec_info->iprntidx),n_secondaries_2);
	
	// further parentage information - still arrays of size nscndprt. Useful?
	parent_mom_at_sec_creation = basic_array<float(*)[3]>(intptr_t(sec_info->pprnt),n_secondaries_2);
	parent_init_pos = basic_array<float(*)[3]>(intptr_t(sec_info->vtxprnt),n_secondaries_2);
	parent_init_mom = basic_array<float(*)[3]>(intptr_t(sec_info->pprntinit),n_secondaries_2);
	parent_trackid = basic_array<int*>(intptr_t(sec_info->iprnttrk),n_secondaries_2); // not implemented in SKG4
	
	return success;
}

int TruthNeutronCaptures_v3::DisableUnusedBranches(){
	std::vector<std::string> used_branches{
	"MC",
	"SECONDARY"
	};
	
	return myTreeReader->OnlyEnableBranches(used_branches);
}

int TruthNeutronCaptures_v3::CreateOutputFile(std::string filename){
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

void TruthNeutronCaptures_v3::ClearOutputTreeBranches(){
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

void TruthNeutronCaptures_v3::PrintBranches(){
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

int TruthNeutronCaptures_v3::WriteTree(){
	Log(m_unique_name+" writing TTree",v_debug,m_verbose);
	outfile->cd();
	// TObject::Write returns the total number of bytes written to the file.
	// It returns 0 if the object cannot be written.
	int bytes = outtree->Write("",TObject::kOverwrite);
	if(bytes<=0){
		Log(m_unique_name+" Error writing TTree!",v_error,m_verbose);
	} else if(m_verbose>2){
		Log(m_unique_name+ " Wrote "+toString(get_ok)+" bytes",v_debug,m_verbose);
	}
	return bytes;
};

void TruthNeutronCaptures_v3::CloseFile(){
	outtree->ResetBranchAddresses();
	// XXX note that while these do persist in the tree:
	// 1. they can't be used in a loop, since you can't use SetBranchAddress on an Alias.
	//    this is probably less of a concern as they can easily be calculated within a loop.
	// 2. they will not show up in TTree::Print(), so you need to know they're there...!
	// Would it be better to store the redundant information? Is there a better way?
	outtree->SetAlias("neutron_travel_dist","sqrt(pow(neutron_start_pos.X()-neutron_end_pos.X(),2)+pow(neutron_start_pos.Y()-neutron_end_pos.Y(),2)+pow(neutron_start_pos.Z()-neutron_end_pos.Z(),2))");
	outtree->SetAlias("neutron_travel_time","(neutron_end_pos.T()-neutron_start_pos.T())/1000."); // [us]
	outtree->SetAlias("neutron_n_daughters","Length$(gamma_energy[])");
	outtree->SetAlias("neutron_tot_gammaE","Sum$(gamma_energy[])");
	outfile->Write("*",TObject::kOverwrite);
	outfile->Close();
	delete outfile;
	outfile = nullptr;
};


