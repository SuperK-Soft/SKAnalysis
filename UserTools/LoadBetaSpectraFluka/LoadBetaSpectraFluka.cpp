/* vim:set noexpandtab tabstop=4 wrap */
#include "LoadBetaSpectraFluka.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include "TFile.h"

LoadBetaSpectraFluka::LoadBetaSpectraFluka():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

const std::map<std::pair<int,int>,std::string> isotope_AZ{
	{{11,4},"11Be"},
	{{16,7},"16N"},
	{{15,6},"15C"},
	{{8,3},"8Li"},
	{{8,5},"8B"},
	{{9,3},"9Li"},
	{{9,6},"9C"},
	{{8,2},"8He"},
	{{12,5},"12B"},
	{{12,7},"12N"}
};

// we only need this
const std::map<int,std::string> isotope_Z{
	{1,"H"},
	{2,"He"},
	{3,"Li"},
	{4,"Be"},
	{5,"B"},
	{6,"C"},
	{7,"N"},
	{8,"O"}
};

bool LoadBetaSpectraFluka::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("histosFile",histosFile);          // root file to write
	m_variables.Get("mapsFile",mapsFile);              // booststore file to write
	m_variables.Get("maxEvents",maxEvents);            // user limit to number of events to process
	m_variables.Get("bonsai_goodness_cut",bonsai_goodness_cut); // skip poorly reconstructed events
	m_variables.Get("bonsai_maxE_cut",bonsai_maxE_cut); // skip events with very high reconstructed E
	
	// open the input TFile and TTree
	// ------------------------------
	std::cout<<"loading file"<<std::endl;
	get_ok = myTreeReader.Load(inputFile, "decTree");
	std::cout<<"load file: "<<get_ok<<std::endl;
	
	// initialize our histos and maps
	for(auto&& an_isotope : isotope_AZ){
		std::string isotope = an_isotope.second; // name
		true_spectra[isotope] = TH1D{(isotope+"_true").c_str(),(isotope+"_true").c_str(),100,0.,16.};
		reco_spectra[isotope] = TH1D{(isotope+"_reco").c_str(),(isotope+"_reco").c_str(),100,0.,16.};
		reco_over_true_spectra[isotope] = TH1D{(isotope+"_ratio").c_str(),(isotope+"_ratio").c_str(),100,0.,16.};
		
		true_events_below_8MeV.emplace(isotope,0);
		true_events_above_8MeV.emplace(isotope,0);
		true_events_below_6MeV.emplace(isotope,0);
		true_events_above_6MeV.emplace(isotope,0);
		
		reco_events_below_8MeV.emplace(isotope,0);
		reco_events_above_8MeV.emplace(isotope,0);
		reco_events_below_6MeV.emplace(isotope,0);
		reco_events_above_6MeV.emplace(isotope,0);
	}
	
	return true;
}


bool LoadBetaSpectraFluka::Execute(){
	
	if((entrynum%1000)==0) Log(toolName+" getting entry "+toString(entrynum),v_debug,verbosity);
	
	// retrieve desired branches
	get_ok = GetBranches();
	
	// process the data
	Analyse();
	
	// move to next entry
	entrynum++;
	// check if we've hit the user-requested entry limit
	if((maxEvents>0)&&(entrynum==maxEvents)){
		Log(toolName+" hit max events, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
		return 1;
	}
	
	// pre-load the next ttree entry
	Log(toolName+" Preloading next entry",v_debug,verbosity);
	get_ok = ReadEntry(entrynum);
	if(get_ok==0){
		return 1; // end of file
	} else if (get_ok<0){
		return 0; // read error
	}
	
	Log(toolName+" Done",v_debug,verbosity);
	return true;
}

bool LoadBetaSpectraFluka::Analyse(){
	
	std::string isotope = std::to_string(A)+isotope_Z.at(Z);
	//Log(toolName+" checking isotope "+isotope,v_debug,verbosity);
	if(true_spectra.count(isotope)==0){
		if(unknown_isotopes.count(isotope)==0){
			std::cout<<"unknown isotope "<<isotope<<std::endl;
			unknown_isotopes.emplace(isotope,0);
		}
		return true;
	}
	
	// skip poorly reconstructed events to prevent distortion of the reconstructed spectra
	// from affecting the selected efficiency
	Log(toolName+" Bonsai goodness is "+toString(bonsai_goodness),v_debug,verbosity);
	if(bonsai_goodness<bonsai_goodness_cut) return true;
	
	// we also seem to have events where the bonsai energy is up at like 10k
	// even with a goodness > 0.4, so skip those
	if(bonsai_energy>bonsai_maxE_cut) return true;
	
	// when comparing true and reconstructed energy
	// for decays with both gammas and betas, both contribute
	// so we should use true energy that is the sum of both
	double true_total_energy = true_photon_E + true_beta_E;
	
	Log(toolName+" Filling histos",v_debug,verbosity);
	true_spectra.at(isotope).Fill(true_total_energy);
	reco_spectra.at(isotope).Fill(bonsai_energy);
	reco_over_true_spectra.at(isotope).Fill(bonsai_energy/true_total_energy);
	
	Log(toolName+" incrementing counters",v_debug,verbosity);
	if(true_total_energy < 8.0) ++true_events_below_8MeV.at(isotope);
	else ++true_events_above_8MeV.at(isotope);
	
	if(true_total_energy < 6.0) ++true_events_below_6MeV.at(isotope);
	else ++true_events_above_6MeV.at(isotope);
	
	if(bonsai_energy < 8.0) ++reco_events_below_8MeV.at(isotope);
	else ++reco_events_above_8MeV.at(isotope);
	
	if(bonsai_energy < 6.0) ++reco_events_below_6MeV.at(isotope);
	else ++reco_events_above_6MeV.at(isotope);
	
	return true;
}


bool LoadBetaSpectraFluka::Finalise(){
	
	Log(toolName+" Making output ROOT file "+histosFile,v_debug,verbosity);
	TFile* fout = new TFile(histosFile.c_str(),"RECREATE");
	fout->cd();
	
	// write out our histos and maps
	Log(toolName+" Writing histos",v_debug,verbosity);
	for(auto&& an_isotope : isotope_AZ){
		std::string isotope = an_isotope.second; // name
		
		true_spectra.at(isotope).Write();
		reco_spectra.at(isotope).Write();
		reco_over_true_spectra.at(isotope).Write();
		
		// h'mm, how about these....
//		true_events_below_8MeV.emplace(isotope,0);
//		true_events_above_8MeV.emplace(isotope,0);
//		true_events_below_6MeV.emplace(isotope,0);
//		true_events_above_6MeV.emplace(isotope,0);
//		
//		reco_events_below_8MeV.emplace(isotope,0);
//		reco_events_above_8MeV.emplace(isotope,0);
//		reco_events_below_6MeV.emplace(isotope,0);
//		reco_events_above_6MeV.emplace(isotope,0);
	}
	
	Log(toolName+" Closing output ROOT file",v_debug,verbosity);
	fout->Write("*",TObject::kOverwrite);
	fout->Close();
	delete fout;
	fout = nullptr;
	
	// let's put the maps in a store, since that's what we'll need anyway
	Log(toolName+" Making output BoostStore",v_debug,verbosity);
	BStore outStore(true);
	outStore.Initnew(mapsFile, uncompressed, true);
	outStore.Set("true_events_below_8MeV",true_events_below_8MeV);
	outStore.Set("true_events_above_8MeV",true_events_above_8MeV);
	outStore.Set("true_events_below_6MeV",true_events_below_6MeV);
	outStore.Set("true_events_above_6MeV",true_events_above_6MeV);
	
	outStore.Set("reco_events_below_8MeV",reco_events_below_8MeV);
	outStore.Set("reco_events_above_8MeV",reco_events_above_8MeV);
	outStore.Set("reco_events_below_6MeV",reco_events_below_6MeV);
	outStore.Set("reco_events_above_6MeV",reco_events_above_6MeV);
	
	// save BoostStore
	Log(toolName+" Saving BoostStore",v_debug,verbosity);
	outStore.Save();
	outStore.Close(); // necessary to complete the file write!
	
	Log(toolName+" Done",v_debug,verbosity);
	return true;
}

int LoadBetaSpectraFluka::ReadEntry(long entry_number){
	// load next entry data from TTree
	int bytesread = myTreeReader.GetEntry(entry_number);
	
	// stop loop if we ran off the end of the tree
	if(bytesread<1&&bytesread>-3){
		Log(toolName+" hit end of input file, stopping loop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
	}
	// stop loop if we had an error of some kind
	else if(bytesread<0){
		 if(bytesread==-1) Log(toolName+" IO error loading next input entry!",v_error,verbosity);
		 if(bytesread==-10) Log(toolName+" AutoClear error loading next input entry!",v_error,verbosity);
		 if(bytesread <-2) Log(toolName+" Unknown error "+toString(bytesread)
		                       +" loading next input entry!",v_error,verbosity);
		 m_data->vars.Set("StopLoop",1);
	}
	
	return bytesread;
}

int LoadBetaSpectraFluka::GetBranches(){
	int success = (
//		(myTreeReader.GetBranchValue("muonID",muonID))               && // primary muon number
//		(myTreeReader.GetBranchValue("nevents",event_num))           && // seem to be a duplicate of entry num
		(myTreeReader.GetBranchValue("Z",Z))                         && // isotope Z
		(myTreeReader.GetBranchValue("A",A))                         && // isotope A
//		(myTreeReader.GetBranchValue("dectime",decay_time))          && // relative to muon time
//		(myTreeReader.GetBranchValue("npart",npart))                 && // ==1 if no photon, ==2 if photons
//		(myTreeReader.GetBranchValue("true_pos",true_decay_pos))     && // 
		(myTreeReader.GetBranchValue("true_betaene",true_beta_E))    && // 
		(myTreeReader.GetBranchValue("true_phene",true_photon_E))    && // 0=none. Seems at most 1 given by Fluka.
//		(myTreeReader.GetBranchValue("bons_pos",bonsai_pos))         && // position reconstructed by Bonsai
		(myTreeReader.GetBranchValue("bons_ene",bonsai_energy))      && // energy reconstructed by Bonsai
		(myTreeReader.GetBranchValue("bons_good",bonsai_goodness))      // goodness of Bonsai reconstruction
	);
	
	return success;
}

int LoadBetaSpectraFluka::DisableUnusedBranches(){
	std::vector<std::string> used_branches{
		// list actively used branches here
//		"muonID",
//		"nevents",
		"Z",
		"A",
//		"dectime",
//		"npart",
//		"true_pos",
		"true_betaene",
//		"true_phene",
//		"bons_pos",
		"bons_ene"
		"bons_good"
	};
	
	return myTreeReader.OnlyEnableBranches(used_branches);
}

/* about the file:
File courtesy of Alice Coffani, generated via FLUKA simulations of primary muons.
Muons were generated with a XXX energy spectra and injected into water? XXX
All isotopes in the corresponding shower were then recorded, along with the energies of
photons and beta particles emitted by their subsequent decay.
These secondary products were then (independently, one decay at a time? XXX)
propagated in skdetsim and reconstructed with bonsai
The resulting file has one entry per spallation decay.
File originally located at:
/home/acoffani/skdetsim/simulations/decayprod_rootfiles/output/myrootfiles/mergedall_myroot_newsampling.root
*/
