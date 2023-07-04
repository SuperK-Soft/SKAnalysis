/* vim:set noexpandtab tabstop=4 wrap */
#include "MakeNCaptTree.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"


MakeNCaptTree::MakeNCaptTree():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool MakeNCaptTree::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	std::string treeReaderName;
	m_variables.Get("treeReaderName",treeReaderName);  // name of input TreeReader
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEvents",maxEvents);            // user limit to number of events to process
	m_variables.Get("writeFrequency",WRITE_FREQUENCY); // how many events to TTree::Fill between TTree::Writes
	
	// get the TreeReader(s)
	// ---------------------
	if(m_data->Trees.count(treeReaderName)==0){
		Log(toolName+" failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
		intree= myTreeReader->GetTree();
	}
	
	// open the output TFile and TTree
	// -------------------------------
	outfile = new TFile(outputFile.c_str(),"RECREATE");
	friendTree = new TTree("ntree","Process Variables");
	friendTree->Branch("neutrino_momentum",&neutrino_momentump,32000,0);
	friendTree->Branch("muon_momentum",&muon_momentump,32000,0);
	// travel distance components relative to neutrino dir
	friendTree->Branch("neutron_longitudinal_travel",&neutron_longitudinal_travel,32000,0);
	friendTree->Branch("neutron_perpendicular_travel",&neutron_perpendicular_travel,32000,0);
	
	// to break these down on a per-capture (not per-event) basis we need to store them in branches
	friendTree->Branch("neutron_n_gammas",&neutron_n_gammas,32000,0);
	friendTree->Branch("neutron_n_electrons",&neutron_n_electrons,32000,0);
	friendTree->Branch("neutron_n_daughters",&neutron_n_daughters,32000,0);
	friendTree->Branch("neutron_tot_gammaE",&total_gamma_energy,32000,0);
	friendTree->Branch("neutron_tot_electronE",&total_electron_energy,32000,0);
	friendTree->Branch("neutron_tot_daughterE",&total_daughter_energy,32000,0);
//	friendTree->Branch("neutron_travel_dist",&placeholder,32000,0);
//	friendTree->Branch("neutron_travel_time",&placeholder,32000,0);
//	friendTree->Branch("neutron_n_daughters",&placeholder,32000,0);
	
	return true;
}


bool MakeNCaptTree::Execute(){
	Log(toolName+" processing entry "+toString(entrynum),v_debug,verbosity);
	entrynum++;
	
	get_ok = GetBranches();
	if(!get_ok) return false;
	
	get_ok = FillFriend();
	return get_ok;
}

int MakeNCaptTree::FillFriend(){
	// don't carry anything over
	Log(toolName+" clearing output ttree variables",v_debug,verbosity);
	ClearOutputTreeBranches();
	
	// loop over primaries and extract the neutrino and primary muon,
	// since we want their momenta for later derived values
	Log(toolName+" getting primary nu/mu momenta",v_debug,verbosity);
	int neutrino_pdg = 12;    // skdetsim has just one neutrino type, which gets saved as ν-e FIXME for SKG4?
	int muon_pdg = 13;
	for(size_t primary_i=0; primary_i<primary_pdg->size(); ++primary_i){
		int next_primary_pdg = primary_pdg->at(primary_i);
		if(next_primary_pdg==neutrino_pdg){
			neutrino_momentum = primary_start_mom->at(primary_i);
		} else if(next_primary_pdg==muon_pdg){
			muon_momentum = primary_start_mom->at(primary_i);
		}
	}
	
	// loop over neutrons in this entry and build the auxilliary info for the friend tree
	Log(toolName+" calculating neutron travel components",v_debug,verbosity);
	for(size_t neutron_i=0; neutron_i<neutron_start_pos->size(); ++neutron_i){
		// longitudinal distance = (neutron_travel_vector).(neutrino_direction_vector)
		TVector3 neutron_travel_vector = 
			neutron_end_pos->at(neutron_i).Vect() - neutron_start_pos->at(neutron_i).Vect();
		double next_neutron_longitudinal_travel = neutron_travel_vector.Dot(neutrino_momentum.Unit());
		double next_neutron_perpendicular_travel = 
			sqrt(neutron_travel_vector.Mag2()-pow(next_neutron_longitudinal_travel,2.)); // TODO fix sign?
		neutron_longitudinal_travel.push_back(next_neutron_longitudinal_travel);
		neutron_perpendicular_travel.push_back(next_neutron_perpendicular_travel);
		
		double total_gamma_E=0;
		for(auto&& agamma : gamma_energy->at(neutron_i)){
			total_gamma_E+=agamma;
		}
		total_gamma_energy.push_back(total_gamma_E);
		neutron_n_gammas.push_back(gamma_energy->at(neutron_i).size());
		
		double total_electron_E=0;
		for(auto&& aelectron : electron_energy->at(neutron_i)){
			total_electron_E+=aelectron;
		}
		total_electron_energy.push_back(total_electron_E);
		neutron_n_electrons.push_back(electron_energy->at(neutron_i).size());
		neutron_n_daughters.push_back(neutron_n_gammas.back()+neutron_n_electrons.back());
		total_daughter_energy.push_back(total_electron_E+total_gamma_E);
		
	}
	
	// XXX any further event-wise info we want to add to the friend tree?
	Log(toolName+" filling friendTree",v_debug,verbosity);
	friendTree->Fill();
	if((entrynum%WRITE_FREQUENCY)==0) WriteTree();
	
//	// angle between vector 'd' (from nu intx vertex and neutron capture)
//	// and "inferred neutron momentum (direction)", 'p', calculated somehow from CCQE assumption...?
//	// expect neutron and muon to have sum of transverse momentum = 0,
//	// -> neutron should be emitted in plane of muon/neutrino tracks,
//	// if we know neutrino direction and muon momentum, we can know muon pt,
//	// which should be balanced by neutron pt (CCQE assumption)
//	// .... but then what? How do we know pl to know expected neutron angle?
//	// is it just angle from neutrino direction??
//	theta_inferred = acos(d.p)/|d|
	
	// TODO:
	// check neutron travel distance as a function of concentration?
	// check neutron capture fraction on nuclei as a function of concentration?
	// compare travel distance to expected distribution for given energy distribution+concentration?
	
	// TODO stretch goals:
	// calculate Evis
	// calculate detection efficiency (efficiency of successful identification? reconstruction? ) (expect 90% on Gd, 25% on H)
	// capture vertex resolution
	// travel distance resolution
	// neutron multiplicity? vs muon pt? - need to revert TruthNeutronCaptures tool to version that keeps nuclides
	
	return 1;
}


bool MakeNCaptTree::Finalise(){
	
	// write out the friend tree
	Log(toolName+" writing output TTree",v_debug,verbosity);
	outfile->cd();
	friendTree->Write("",TObject::kOverwrite);
	
	// this doesn't persist to the output file, so a bit pointless really...
	// in fact, it also relies on branches from the main tree.
	friendTree->SetAlias("neutron_travel_dist",
	                     "sqrt(pow(neutron_longitudinal_travel,2.)+pow(neutron_perpendicular_travel,2.))");
	
	// apparently the list of Friends do in fact get written to file...
	friendTree->AddFriend(intree->GetName(), intree->GetCurrentFile()->GetName());
	// of course the corresponding file will probably not remain at that location...
	
	Log(toolName+" cleanup",v_debug,verbosity);
	friendTree->ResetBranchAddresses();
	if(outfile){ outfile->Close(); delete outfile; outfile=nullptr; }
	
	return true;
}

int MakeNCaptTree::GetBranches(){
	int success = (
//	(myTreeReader->GetBranchValue("filename",filename))                         &&
//	(myTreeReader->GetBranchValue("water_transparency",water_transparency))     &&
//	(myTreeReader->GetBranchValue("skdetsim_version",skdetsim_version))         &&
//	(myTreeReader->GetBranchValue("tba_table_version",tba_table_version))       &&
//	(myTreeReader->GetBranchValue("entry_number",entry_number))                 &&
//	(myTreeReader->GetBranchValue("subevent_num",subevent_number))              &&
	(myTreeReader->GetBranchValue("primary_pdg",primary_pdg))                   &&
//	(myTreeReader->GetBranchValue("primary_energy",primary_energy))             &&
	(myTreeReader->GetBranchValue("primary_start_mom",primary_start_mom))       &&
//	(myTreeReader->GetBranchValue("primary_start_pos",primary_start_pos))       &&
//	(myTreeReader->GetBranchValue("primary_end_pos",primary_end_pos))           &&
//	(myTreeReader->GetBranchValue("nuclide_parent_pdg",nuclide_parent_pdg))     &&
//	(myTreeReader->GetBranchValue("nuclide_creation_pos",nuclide_creation_pos)) &&
//	(myTreeReader->GetBranchValue("nuclide_decay_pos",nuclide_decay_pos))       &&
	(myTreeReader->GetBranchValue("nuclide_daughter_pdg",nuclide_daughter_pdg)) &&
	(myTreeReader->GetBranchValue("neutron_start_pos",neutron_start_pos))       &&
	(myTreeReader->GetBranchValue("neutron_end_pos",neutron_end_pos))           &&
//	(myTreeReader->GetBranchValue("neutron_start_energy",neutron_start_energy)) &&
//	(myTreeReader->GetBranchValue("neutron_end_energy",neutron_end_energy))     &&
//	(myTreeReader->GetBranchValue("neutron_end_process",neutron_end_process))   &&
	(myTreeReader->GetBranchValue("gamma_energy",gamma_energy))                 &&
//	(myTreeReader->GetBranchValue("gamma_time",gamma_time))                     &&
	(myTreeReader->GetBranchValue("electron_energy",electron_energy))         //&&
//	(myTreeReader->GetBranchValue("electron_time",electron_time))
	);
	if(!success){
		Log(toolName+" Error! failed to get all required branch values!",v_error,verbosity);
	}
	return success;
}

void MakeNCaptTree::ClearOutputTreeBranches(){
	neutrino_momentum.SetXYZ(0,0,0);
	muon_momentum.SetXYZ(0,0,0);
	neutron_longitudinal_travel.clear();
	neutron_perpendicular_travel.clear();
	total_gamma_energy.clear();
	total_electron_energy.clear();
	total_daughter_energy.clear();
	neutron_n_gammas.clear();
	neutron_n_electrons.clear();
	neutron_n_daughters.clear();
}

int MakeNCaptTree::WriteTree(){
	Log(toolName+" writing TTree",v_debug,verbosity);
	outfile->cd();
	// TObject::Write returns the total number of bytes written to the file.
	// It returns 0 if the object cannot be written.
	int bytes = friendTree->Write("",TObject::kOverwrite);
	if(bytes<=0){
		Log(toolName+" Error writing TTree!",v_error,verbosity);
	} else if(verbosity>2){
		Log(toolName+ " Wrote "+toString(get_ok)+" bytes",v_debug,verbosity);
	}
	return bytes;
}

