/* vim:set noexpandtab tabstop=4 wrap */
#include "TruthNeutronCaptures_v2.h"
#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <algorithm> // std::find

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"

TruthNeutronCaptures_v2::TruthNeutronCaptures_v2():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool TruthNeutronCaptures_v2::Initialise(std::string configfile, DataModel &data){
	
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
	get_ok = myTreeReader.Load(inputFile, "data"); // SKROOT TTree is named descriptively 'data'
	DisableUnusedBranches();
	if(get_ok) ReadEntry(0);
	
	// note a couple of file-wise constants
	out_skdetsim_version = mc_info->ivmcp;    // version of monte carlo program (SK_GEOMETRY+1000)
	out_tba_table_version = mc_info->ivabl;   // version of absorbtion length
	out_water_transparency = -1.f;
	out_subevent_number = -1;
	
	// create the output TFile and TTree
	// ---------------------------------
	CreateOutputFile(outputFile);
	
	return true;
}

bool TruthNeutronCaptures_v2::Execute(){
	
	Log(toolName+" processing entry "+toString(entry_number),v_debug,verbosity);
	
	// clear output vectors so we don't carry anything over
	Log(toolName+" clearing output vectors",v_debug,verbosity);
	ClearOutputTreeBranches();
	
	// Set output variables
	Log(toolName+" setting output variables",v_debug,verbosity);
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
		Log(toolName+" reading entry "+toString(entry_number),v_debug,verbosity);
		get_ok = ReadEntry(entry_number);
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
	
	// handy constants
	neutron_mass = PdgToMass(StringToPdg("Neutron"));
	
	return true;
}


bool TruthNeutronCaptures_v2::Finalise(){
	
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

int TruthNeutronCaptures_v2::CalculateVariables(){
	
	out_filename = myTreeReader.GetTree()->GetCurrentFile()->GetName();
	out_run_number = mc_info->mcrun;
	out_entry_number = entry_number;    // TTree entry number to be able to identify the source event
	
	std::map<int,int> neutrons_map;
	
	Log(toolName+" looping over "+toString(sec_info->track_g3_code.size())
	            +" particles in event "+toString(entry_number),v_debug,verbosity);
	for(unsigned int particle_i=0; particle_i<sec_info->track_g3_code.size(); ++particle_i){
		// get indices for vertex info
		int creation_vtx_index = sec_info->track_creation_vtx.at(particle_i);
		int termination_vtx_index = sec_info->track_termination_vtx.at(particle_i);
		Log(toolName+" particle "+toString(particle_i)+" creation vtx index "+toString(creation_vtx_index)
				 +", termination_vtx_index "+toString(termination_vtx_index),v_debug,verbosity);
		
		// store primary particles
		if(sec_info->track_parent.at(particle_i)<0){
			Log(toolName+" storing primary",v_debug,verbosity);
			out_primary_pdg.push_back(G3ParticleCodeToPdg(sec_info->track_g3_code.at(particle_i)));
			// store initial momentum
			out_primary_start_mom.emplace_back(sec_info->track_ini_momentum.at(particle_i)[0],
			                                   sec_info->track_ini_momentum.at(particle_i)[1],
			                                   sec_info->track_ini_momentum.at(particle_i)[2]);
			// calculate initial energy - is this relativistic?
			float particle_e = sqrt(
				pow(sec_info->track_ini_momentum.at(particle_i)[0],2.f)+
				pow(sec_info->track_ini_momentum.at(particle_i)[1],2.f)+
				pow(sec_info->track_ini_momentum.at(particle_i)[2],2.f));
			out_primary_energy.push_back(particle_e*1000.0f); // convert GeV to MeV
			// store initial and final positions
			if((creation_vtx_index>=0)&&(creation_vtx_index<int(sec_info->vertex_time.size()))){
				out_primary_start_pos.emplace_back(sec_info->vertex_pos.at(creation_vtx_index)[0],
					                               sec_info->vertex_pos.at(creation_vtx_index)[1],
					                               sec_info->vertex_pos.at(creation_vtx_index)[2],
					                               sec_info->vertex_time.at(creation_vtx_index)
					                               +sec_info->track_creation_toffset.at(particle_i));
			} else {
				out_primary_start_pos.emplace_back(0,0,0,0);
			}
			if((termination_vtx_index>=0)&&(termination_vtx_index<int(sec_info->vertex_time.size()))){
				out_primary_end_pos.emplace_back(sec_info->vertex_pos.at(termination_vtx_index)[0],
					                             sec_info->vertex_pos.at(termination_vtx_index)[1],
					                             sec_info->vertex_pos.at(termination_vtx_index)[2],
					                             sec_info->vertex_time.at(termination_vtx_index));
			} else {
				out_primary_end_pos.emplace_back(0,0,0,0);
			}
		}
		
		// store neutrons and parent nuclides
		if(sec_info->track_g3_code.at(particle_i)==13){
			Log(toolName+" storing neutron",v_debug,verbosity);
			// TODO could store initial momentum, now that we have it
			// but do we need it? re-think output file given info
			// we now have available...
			
			// calculate initial energy - is this relativistic?
			float particle_e = sqrt(
				pow(sec_info->track_ini_momentum.at(particle_i)[0],2.f)+
				pow(sec_info->track_ini_momentum.at(particle_i)[1],2.f)+
				pow(sec_info->track_ini_momentum.at(particle_i)[2],2.f));
			out_neutron_start_energy.push_back(particle_e*1000.0f); // convert GeV to MeV
			
			// energy at capture
			if((termination_vtx_index>=0)&&(termination_vtx_index<int(sec_info->vertex_time.size()))){
				particle_e = sqrt(
					pow(sec_info->vertex_incident_particle_momentum.at(termination_vtx_index)[0],2.f)+
					pow(sec_info->vertex_incident_particle_momentum.at(termination_vtx_index)[1],2.f)+
					pow(sec_info->vertex_incident_particle_momentum.at(termination_vtx_index)[2],2.f));
				out_neutron_end_energy.push_back(particle_e*1000.0f); // convert GeV to MeV
				// XXX always 1GeV exactly... same as initial energy... bug?
			} else {
				out_neutron_end_energy.push_back(-1.);
			}
			// initial position
			if((creation_vtx_index>=0)&&(creation_vtx_index<int(sec_info->vertex_time.size()))){
				out_neutron_start_pos.emplace_back(sec_info->vertex_pos.at(creation_vtx_index)[0],
					                               sec_info->vertex_pos.at(creation_vtx_index)[1],
					                               sec_info->vertex_pos.at(creation_vtx_index)[2],
					                               sec_info->vertex_time.at(creation_vtx_index)
					                               +sec_info->track_creation_toffset.at(particle_i));
			} else {
				out_neutron_start_pos.emplace_back(0,0,0,0);
			}
			// final position and termination process
			if((termination_vtx_index>=0)&&(termination_vtx_index<int(sec_info->vertex_time.size()))){
				out_neutron_end_pos.emplace_back(sec_info->vertex_pos.at(termination_vtx_index)[0],
					                             sec_info->vertex_pos.at(termination_vtx_index)[1],
					                             sec_info->vertex_pos.at(termination_vtx_index)[2],
					                             sec_info->vertex_time.at(termination_vtx_index));
				out_neutron_end_process.push_back(sec_info->vertex_g3_process_codes.at(termination_vtx_index).at(0));
			} else {
				out_neutron_end_pos.emplace_back(0,0,0,0);
				out_neutron_end_process.push_back(-1);
			}
			// number of daughters
			neutrons_map.emplace(particle_i,out_neutron_ndaughters.size());
			out_neutron_ndaughters.push_back(0); // placeholder
			// capture nuclide pdg
			if((creation_vtx_index>=0) && (creation_vtx_index<int(sec_info->vertex_time.size()))){
				out_nuclide_parent_pdg.push_back(sec_info->vertex_target_g3_code.at(creation_vtx_index));
			}
			// daughter nuclide pdg
			out_nuclide_daughter_pdg.push_back(0); // placeholder
		}
	} // end first loop over particles
	
	// pre-allocate a vector of gammas for each neutron
	out_gamma_energy.resize(out_neutron_start_pos.size());
	out_gamma_time.resize(out_neutron_start_pos.size());
	out_electron_energy.resize(out_neutron_start_pos.size());
	out_electron_time.resize(out_neutron_start_pos.size());
	Log(toolName+" pre-allocating "+toString(out_neutron_start_pos.size())
	            +" products vectors for neutron captures",v_debug,verbosity);
	
	// scan again, this time filling gamma and daughter nuclide information
	// the order of recording may be such that this re-scan isn't necessary
	// (if daughters always follow the neutrons they came from)
	// but for now, i can't remember if that's necessarily the case
	for(unsigned int particle_i=0; particle_i<sec_info->track_g3_code.size(); ++particle_i){
		int creation_vtx_index = sec_info->track_creation_vtx.at(particle_i);
		int parent_index = sec_info->track_parent.at(particle_i);
		bool from_ncap=false;  // check all process codes for the creation vertex for any that are ncapture
		for(auto&& aproc: sec_info->vertex_g3_process_codes.at(creation_vtx_index)){
			from_ncap = from_ncap || (aproc==18);
		}
		// TODO add error checking on these
		// 1. valid creation vertex
		// 2. from neutron capture
		// 3. with valid parent neutron
		// 4. with valid neutron capture event
		if( (creation_vtx_index>=0) && (creation_vtx_index<int(sec_info->vertex_time.size())) &&
			from_ncap &&
			(parent_index>=0) && (parent_index<int(sec_info->track_g3_code.size())) &&
			neutrons_map.count(parent_index)
			){
			
			// get index of parent neutron in the array of neutrons
			int parent_neutron_index = neutrons_map.at(parent_index);
			
			// save gammas
			if(sec_info->track_g3_code.at(particle_i)==1){
				Log(toolName+" storing decay gamma",v_debug,verbosity);
				// add this gamma to the appropriate position
//				out_gamma_time.at(parent_neutron_index).push_back(sec_info->vertex_time.at(creation_vtx_index)
//				                                                 +sec_info->track_creation_toffset.at(particle_i));
				out_gamma_time.at(parent_neutron_index).push_back(sec_info->track_creation_toffset.at(particle_i));
				// calculate energy
				float particle_e = sqrt(
					pow(sec_info->track_ini_momentum.at(particle_i)[0],2.f)+
					pow(sec_info->track_ini_momentum.at(particle_i)[1],2.f)+
					pow(sec_info->track_ini_momentum.at(particle_i)[2],2.f));
				out_gamma_energy.at(parent_neutron_index).push_back(particle_e*1000.0f); // convert GeV to MeV
				out_neutron_ndaughters.at(parent_neutron_index)++;
				if(out_neutron_ndaughters.at(parent_neutron_index)>1){
					// for "captures on hydrogen" we sometimes see more than one gamma from normal skdetsim
					// these come from steps where there are multiple lmec entries, typically 12 and 18
					// i.e. in one step the neutron underwent hadronic scattering, generating a gamma,
					// and neutron capture, generating another gamma. Which is which? Who knows.
					// Moreover, in skdetsim-gd captures on hydrogen frequently have many gammas,
					// although their sum consistently adds up to 2.2MeV
					Log("WARNING! event "+toString(entry_number)+" had "
							 +toString(out_neutron_ndaughters.at(parent_neutron_index))
							 +" daughters from neutron capture!",v_message,verbosity);
					//assert(false);
				}
			}
			// save conversion electrons: these will also contribute to total deexcitation energy
			else if(sec_info->track_g3_code.at(particle_i)==3){
				Log(toolName+" storing electron",v_debug,verbosity);
				// add this gamma to the appropriate position
//				out_electron_time.at(parent_neutron_index).push_back(sec_info->vertex_time.at(creation_vtx_index)
//				                                                 +sec_info->track_creation_toffset.at(particle_i));
				out_electron_time.at(parent_neutron_index).push_back(sec_info->track_creation_toffset.at(particle_i));
				// calculate energy
				float particle_e = sqrt(
					pow(sec_info->track_ini_momentum.at(particle_i)[0],2.f)+
					pow(sec_info->track_ini_momentum.at(particle_i)[1],2.f)+
					pow(sec_info->track_ini_momentum.at(particle_i)[2],2.f));
				out_electron_energy.at(parent_neutron_index).push_back(particle_e*1000.0f); // convert GeV to MeV
				out_neutron_ndaughters.at(parent_neutron_index)++;
				if(out_neutron_ndaughters.at(parent_neutron_index)>1){
					Log("WARNING! event "+toString(entry_number)+" had "
							 +toString(out_neutron_ndaughters.at(parent_neutron_index))
							 +" daughters from neutron capture!",v_message,verbosity);
					//assert(false);
				}
			}
			// save daughter nuclide
			else {
				Log(toolName+" storing daughter nuclide",v_debug,verbosity);
				// sanity check that there really is only one daughter nuclide
				if(out_nuclide_daughter_pdg.at(parent_neutron_index)!=0){
					Log(toolName+" ERROR: more than one daughter nuclide for neutron "+toString(parent_neutron_index)
						+"; current is "+toString(out_nuclide_daughter_pdg.at(parent_neutron_index))+", new is "
						+toString(G3ParticleCodeToPdg(sec_info->track_g3_code.at(particle_i))),v_error,verbosity);
				}
				out_nuclide_daughter_pdg.at(parent_neutron_index) = G3ParticleCodeToPdg(sec_info->track_g3_code.at(particle_i));
			}
		} // if daughter from neutron capture
	} // end second loop over particles
	
	return 1;
}

int TruthNeutronCaptures_v2::ReadEntry(long entry_number){
	int bytesread = myTreeReader.GetEntry(entry_number);
	if(bytesread<=0) return bytesread;
	
	int success = 
	(myTreeReader.GetBranchValue("HEADER",run_header)) &&
	(myTreeReader.GetBranchValue("MC",mc_info)) &&
	(myTreeReader.GetBranchValue("SECONDARY",sec_info));
	
	
	// XXX for efficiency, add all used branches to DisableUnusedBranches XXX
	
	return success;
}

int TruthNeutronCaptures_v2::DisableUnusedBranches(){
	std::vector<std::string> used_branches{
		// branchname     || classname
		// ---------------------------
		"HEADER",         // Header
//		"TQREAL",         // TQReal
//		"TQAREAL",        // TQReal
//		"LOWE",           // LoweInfo
//		"ATMPD",          // AtmpdInfo
//		"UPMU",           // UpmuInfo
//		"MU",             // MuInfo
//		"SLE",            // SLEInfo
//		"SWTRGLIST",      // SoftTrgList
//		"IDODXTLK",       // IDOD_Xtlk
		"MC",             // MCInfo
		"SECONDARY"       // SecondaryInfo
		// ---------------------------
		// the following branches are added during analysis
//		"TQLIST",         // TClonesArray
//		"ODTQLIST",       // TClonesArray
//		"HWTRGLIST",      // TClonesArray
//		"PEDESTALS",      // Pedestal
//		"EVENTHEADER",    // EventHeader
//		"EVENTTRAILER",   // EventTrailer
//		"SOFTWARETRG",    // SoftwareTrigger
//		"QBEESTATUS",     // QBeeStatus
//		"DBSTATUS",       // QBeeStatus
//		"SPACERS",        // Spacer
//		"PREVT0",         // PrevT0
//		"MISMATCHEDHITS", // MismatchedHit
//		"GPSLIST",        // TClonesArray
//		"T2KGPSLIST"      // TClonesArray
		// N.B. RunInfo (RINFO) from 'RareList' is extracted
		// from TTree UserInfo during TreeManager::Initialize;
		// probably easiest to use TreeManager for this.
		// N.B. SlowControl (SLWCTRL) TreeManager Getter
		// returns a nullptr so is not included here.
	};
	
	return myTreeReader.OnlyEnableBranches(used_branches);
}

int TruthNeutronCaptures_v2::CreateOutputFile(std::string filename){
	// create the output ROOT file and TTree for writing
	// =================================================
	outfile = new TFile(filename.c_str(), "RECREATE");
	outtree = new TTree("eventtree", "Events with Neutron Captures");
	
	// create branches
	// ---------------
	// file level
	outtree->Branch("filename",&out_filename);
	outtree->Branch("skdetsim_version",&out_skdetsim_version);
	outtree->Branch("tba_table_version",&out_tba_table_version);
	outtree->Branch("water_transparency",&out_water_transparency);  // where?
	
	// event level
	outtree->Branch("run_number",&out_run_number);
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
	outtree->Branch("neutron_n_daughters",&out_neutron_ndaughters,32000,0);
	
	
	// gamma
	outtree->Branch("gamma_energy",&out_gamma_energy,32000,0);
	outtree->Branch("gamma_time",&out_gamma_time,32000,0);
	// use outtree->Draw("Length$(gamma_energy[])"); to draw gamma multiplicity FOR THE EVENT (combines captures!)
	// use outtree->Draw("Sum$(gamma_energy[])");    to draw total gamma energy FOR THE EVENT (combines captures!)
	
	// electron
	outtree->Branch("electron_energy",&out_electron_energy,32000,0);
	outtree->Branch("electron_time",&out_electron_time,32000,0);
	
	return 1;
}

void TruthNeutronCaptures_v2::ClearOutputTreeBranches(){
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

void TruthNeutronCaptures_v2::PrintBranches(){
	std::cout<<"==========================================================="<<std::endl;
	std::cout<<"PRINTING EVENT"<<std::endl;
	std::cout<<"filename: "<<out_filename<<std::endl
			 <<"skdetsim_version: "<<out_skdetsim_version<<std::endl
			 <<"tba_table_version: "<<out_tba_table_version<<std::endl;
	
	std::cout<<"entry_number: "<<out_entry_number<<std::endl
			 <<"subevent_number:" <<out_subevent_number<<std::endl
			 <<"num primaries: "<<out_primary_pdg.size()<<std::endl;
	if(out_primary_pdg.size()){
		std::cout<<"primary vertex [cm]:"
				 <<" ("<<out_primary_start_pos.at(0).X()
				 <<", "<<out_primary_start_pos.at(0).Y()
				 <<", "<<out_primary_start_pos.at(0).Z()<<")"<<std::endl;
	}
	
	// print primaries
	for(int primary_i=0; primary_i<out_primary_pdg.size(); ++primary_i){
		std::cout<<"\tprimary ("<<primary_i<<"): "<<std::endl
				 <<"\t\tprimary pdg: "<<out_primary_pdg.at(primary_i)<<std::endl
				 <<"\t\tprimary energy [MeV]: "<<out_primary_energy.at(primary_i)<<std::endl
				 <<"\t\tprimary momentum [MeV/c]:"
				 <<" ("<<out_primary_start_mom.at(primary_i).X()
				 <<", "<<out_primary_start_mom.at(primary_i).Y()
				 <<", "<<out_primary_start_mom.at(primary_i).Z()<<")"<<std::endl;
		std::cout<<"\t\tprimary start pos [cm]:"
				 <<" ("<<out_primary_start_pos.at(primary_i).X()
				 <<", "<<out_primary_start_pos.at(primary_i).Y()
				 <<", "<<out_primary_start_pos.at(primary_i).Z()<<")"<<std::endl
				 <<"\t\tprimary end pos [cm]:"
				 <<" ("<<out_primary_end_pos.at(primary_i).X()
				 <<", "<<out_primary_end_pos.at(primary_i).Y()
				 <<", "<<out_primary_end_pos.at(primary_i).Z()<<")"<<std::endl;
	}
	
	int total_neutrons=0;
	int total_gammas=0;
	int total_electrons=0;
	// print neutron captures
	std::cout<<"num neutrons: "<<out_neutron_start_energy.size()<<std::endl;
	for(int neutron_i=0; neutron_i<out_neutron_start_energy.size(); ++neutron_i){
		std::cout<<"\tneutron "<<neutron_i<<": "<<std::endl
				 <<"\t\tneutron start energy [MeV]: "
				 <<out_neutron_start_energy.at(neutron_i)<<std::endl
				 <<"\t\tneutron start pos [cm]:"
				 <<" ("<<out_neutron_start_pos.at(neutron_i).X()
				 <<", "<<out_neutron_start_pos.at(neutron_i).Y()
				 <<", "<<out_neutron_start_pos.at(neutron_i).Z()<<")"<<std::endl
				 <<"\t\tneutron start time [ns]:"
				 <<out_neutron_start_pos.at(neutron_i).T()<<std::endl
				 <<"\t\tneutron end pos [cm]:"
				 <<" ("<<out_neutron_end_pos.at(neutron_i).X()
				 <<", "<<out_neutron_end_pos.at(neutron_i).Y()
				 <<", "<<out_neutron_end_pos.at(neutron_i).Z()<<")"<<std::endl
				 <<"\t\tneutron end time [ns]:"
				 <<out_neutron_end_pos.at(neutron_i).T()<<std::endl
				 <<"\t\tneutron end energy [MeV]: "
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
					 <<"\t\t\tgamma energy [MeV]: "<<out_gamma_energy.at(neutron_i).at(gamma_i)
					 <<std::endl
					 <<"\t\t\tgamma time [ns]: "<<out_gamma_time.at(neutron_i).at(gamma_i)
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

int TruthNeutronCaptures_v2::WriteTree(){
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

void TruthNeutronCaptures_v2::CloseFile(){
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
