/* vim:set noexpandtab tabstop=4 wrap */
#include "PlotNeutronCaptures.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <memory>  // unique_ptr
#include <map>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPie.h"
#include "TStyle.h"

// TODO move to DataModel RootAlgorithms or something
std::unique_ptr<TPie> GeneratePieFromHisto(TH1F* histo, int verbose=0);
std::unique_ptr<TPie> GeneratePieFromHisto(std::string histoname, int verbose=0);

PlotNeutronCaptures::PlotNeutronCaptures():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool PlotNeutronCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("inputFile",inputFile);            // a single specific input file
	m_variables.Get("outputFile",outputFile);          // output file to write
	m_variables.Get("maxEvents",maxEvents);            // user limit to number of events to process
	m_variables.Get("writeFrequency",WRITE_FREQUENCY); // how many events to TTree::Fill between TTree::Writes
	
	// open the input TFile and TTree
	// ------------------------------
	get_ok = myTreeReader.Load(inputFile, "eventtree");
	intree = myTreeReader.GetTree();
	
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


bool PlotNeutronCaptures::Execute(){
	
	Log(toolName+" processing entry "+toString(entrynum),v_debug,verbosity);
	
	// ttree entry is already loaded so just retrieve the desired branches
	Log(toolName+" getting data from input branches",v_debug,verbosity);
	get_ok = GetBranches();
	
	// process the data
	Log(toolName+" calculating processed data for output tree",v_debug,verbosity);
	get_ok = FillFriend();
	
	// move to next entry
	Log(toolName+" checking if stopping the toolchain",v_debug,verbosity);
	entrynum++;
	// check if we've hit the user-requested entry limit
	if((maxEvents>0)&&(entrynum==maxEvents)){
		Log(toolName+" hit max events, setting StopLoop",v_message,verbosity);
		m_data->vars.Set("StopLoop",1);
		return 1;
	}
	
	// pre-load the next ttree entry
	Log(toolName+" pre-loading entry "+toString(entrynum),v_debug,verbosity);
	get_ok = ReadEntry(entrynum);
	if(get_ok==0){
		return 1; // end of file
	} else if (get_ok<0){
		return 0; // read error
	}
	
	return get_ok;
}

int PlotNeutronCaptures::FillFriend(){
	// don't carry anything over
	Log(toolName+" clearing output ttree variables",v_debug,verbosity);
	ClearOutputTreeBranches();
	
	// loop over primaries and extract the neutrino and primary muon,
	// since we want their momenta for later derived values
	Log(toolName+" getting primary nu/mu momenta",v_debug,verbosity);
	int neutrino_pdg = 12;    // recall that skdetsim has just one neutrino type, which gets saved as ν-e
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
			neutron_travel_vector.Mag()-abs(next_neutron_longitudinal_travel); // TODO fix sign?
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
		
		// keep a map with capture nuclides to num capture events
		int capture_nuclide_pdg = nuclide_daughter_pdg->at(neutron_i);
		std::string capture_nuclide_name = PdgToString(nuclide_daughter_pdg->at(neutron_i));
		if(capture_nuclide_vs_count.count(PdgToString(nuclide_daughter_pdg->at(neutron_i)))){
			capture_nuclide_vs_count.at(PdgToString(nuclide_daughter_pdg->at(neutron_i)))++;
		} else {
			capture_nuclide_vs_count.emplace(PdgToString(nuclide_daughter_pdg->at(neutron_i)),1);
		}
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


bool PlotNeutronCaptures::Finalise(){
	
	// write out the friend tree
	Log(toolName+" writing output TTree",v_debug,verbosity);
	outfile->cd();
	friendTree->Write("",TObject::kOverwrite);
	
	// make and write out histograms
	Log(toolName+" making histograms",v_debug,verbosity);
	MakeHistos();
	
	Log(toolName+" cleanup",v_debug,verbosity);
	if(friendTree) friendTree->ResetBranchAddresses();
	if(outfile){ outfile->Close(); delete outfile; outfile=nullptr; }
	
	return true;
}

int PlotNeutronCaptures::ReadEntry(long entry_number){
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

int PlotNeutronCaptures::MakeHistos(){
	// the lazy way (proper would be to call TH1::Fill during FillFriend)
	// ============
	outfile->cd();
	
	// ======================
	// cumulative plots
	// ======================
	Log(toolName+" making aggregate plots",v_debug,verbosity);
	// neutron energy
	TH1D hNeutronE("hNeutronE","Neutron Energy;Neutron Energy [MeV];Num Events",100,0,250);
	intree->Draw("neutron_start_energy>>hNeutronE");
	hNeutronE.Write();
	
	// neutron travel distance
	TH1D hNeutronTravelDist("hNeutronTravelDist","Neutron Travel Distance;Distance [cm];Num Events",100,0,125);
	intree->Draw("neutron_travel_dist>>hNeutronTravelDist");
	hNeutronTravelDist.Write();
	
	// neutron travel time
	TH1D hNeutronTravelTime("hNeutronTravelTime","Neutron Travel Time;Time [ns];Num Events",100,0,1800E3);
	intree->Draw("neutron_travel_time>>hNeutronTravelTime");
	hNeutronTravelTime.Write();
	
	// gamma mulitiplicity
	TH1D hNumGammas("hNumGammas","Gamma Multiplicity (All Nuclides);Num Gammas Emitted;Num Events",100,0,30);
	friendTree->Draw("neutron_n_gammas>>hNumGammas");
	hNumGammas.Write();
	
	// gamma energy
	TH1D hGammaE("hGammaE", "Gamma Energy (All Nuclides);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE");
	hGammaE.Write();
	
	// total gamma energy from the neutron capture
	TH1D hSumGammaE("hSumGammaE","Total Emitted Gamma Energy (All Nuclides);Sum of Gamma Energy [MeV];Num Events",100,0,10);
	friendTree->Draw("neutron_tot_gammaE>>hSumGammaE");
	hSumGammaE.Write();
	
	// gamma emission time (parent lifetime)
	TH1D hGammaT("hGammaT", "Gamma Emission Time (All Nuclides);Gamma Emission Time [ns];Num Events",100,0,1800E3);
	intree->Draw("gamma_time>>hGammaT");
	hGammaT.Write();
	
	// electron mulitiplicity
	TH1D hNumElectrons("hNumElectrons","Electron Multiplicity (All Nuclides);Num Electrons Emitted;Num Events",100,0,30);
	friendTree->Draw("neutron_n_electrons>>hNumElectrons");
	hNumElectrons.Write();
	
	// electron energy
	TH1D hElectronE("hElectronE", "Electron Energy (All Nuclides);Electron Energy [MeV];Num Events",100,0,10);
	intree->Draw("electron_energy>>hElectronE");
	hElectronE.Write();
	
	// total electron energy from the neutron capture
	TH1D hSumElectronE("hSumElectronE","Total Emitted Electron Energy (All Nuclides);Sum of Electron Energy [MeV];Num Events",100,0,10);
	friendTree->Draw("neutron_tot_electronE>>hSumElectronE");
	hSumElectronE.Write();
	
	// electron emission time (parent lifetime)
	TH1D hElectronT("hElectronT", "Electron Emission Time (All Nuclides);Electron Emission Time [ns];Num Events",100,0,1800E3);
	intree->Draw("electron_time>>hElectronT");
	hElectronT.Write();
	
	// total daughter multiplicity
	TH1D hNumDaughters("hNumDaughters", "Daughter Multiplicity;Num Daughters;Num Events",100,0,30);
	friendTree->Draw("neutron_n_daughters>>hNumDaughters");
	hNumDaughters.Write();
	
	// total daughter energy
	TH1D hSumDaughterE("hSumDaughterE", "Total Daughter Energy;Total Energy [MeV];Num Events",100,0,10);
	friendTree->Draw("neutron_tot_daughterE>>hSumDaughterE");
	hSumDaughterE.Write();
	
	// pie chart of capture nuclei
	Log(toolName+" making pie chart",v_debug,verbosity);
	TPie pCaptureNuclidePdg = TPie("pCaptureNuclidePdg", "Captures by Nuclide",capture_nuclide_vs_count.size());
	int nuclide_i=0;
	for(auto&& anuclide : capture_nuclide_vs_count){
		pCaptureNuclidePdg.SetEntryLabel(nuclide_i,anuclide.first.c_str());
		pCaptureNuclidePdg.SetEntryVal(nuclide_i,anuclide.second);
		++nuclide_i;
	}
	// making it look nice
//	pCaptureNuclidePdg.SetAngularOffset(333);
	pCaptureNuclidePdg.SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
	pCaptureNuclidePdg.SetValueFormat("%4.0f");
	pCaptureNuclidePdg.SetPercentFormat("%3.0f");
	pCaptureNuclidePdg.SetCircle(0.5, 0.4702026, 0.3302274);
	pCaptureNuclidePdg.SetTextSize(0.03455766);
	// saving to file
	pCaptureNuclidePdg.Draw();
	pCaptureNuclidePdg.Write();
	
	// ==============================
	// broken down by capture nucleus
	// ==============================
	Log(toolName+" making total gamma energy stack",v_debug,verbosity);
	auto statsboxdefault = gStyle->GetOptStat();
	gStyle->SetOptStat(0); // turn off stats box; overlaps with legends
	// Total Gamma Energy
	// ------------------
	// capture on H
	TH1D hTotGammaE_H("hTotGammaE_H", "Total Gamma Energy (Capture on H);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("neutron_tot_gammaE>>hTotGammaE_H","nuclide_daughter_pdg==100045");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hTotGammaE_Gd_155("hTotGammaE_Gd_155", "Total Gamma Energy (Capture on Gd-155);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("neutron_tot_gammaE>>hTotGammaE_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hTotGammaE_Gd_157("hTotGammaE_Gd_157", "Total Gamma Energy (Capture on Gd-157);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("neutron_tot_gammaE>>hTotGammaE_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hTotGammaE_Stack("hTotGammaE_Stack","Gamma Spectrum by Capture Nucleus;Gamma Energy [MeV];Num Events");
	hTotGammaE_H.SetLineColor(kRed);
	hTotGammaE_Gd_155.SetLineColor(kBlue);
	hTotGammaE_Gd_157.SetLineColor(kMagenta);
	hTotGammaE_Stack.Add(&hTotGammaE_H);
	hTotGammaE_Stack.Add(&hTotGammaE_Gd_155);
	hTotGammaE_Stack.Add(&hTotGammaE_Gd_157);
	TLegend StackLegend(0.65,0.7,0.88,0.88,NULL);
	StackLegend.SetFillStyle(0);
	StackLegend.SetLineStyle(0);
	StackLegend.AddEntry(&hTotGammaE_H,"Hydrogen","l");
	StackLegend.AddEntry(&hTotGammaE_Gd_155,"Gd-155","l");
	StackLegend.AddEntry(&hTotGammaE_Gd_157,"Gd-157","l");
	hTotGammaE_Stack.Draw();
	StackLegend.Draw();
	// add the legend to the list of functions so that it gets saved on Write call
	// a THStack doesn't have a list of functions, so we have to add it to a component histo
	hTotGammaE_H.GetListOfFunctions()->Add(&StackLegend);
	hTotGammaE_Stack.Write();
	// we need to remove it afterwards, though, otherwise the histograms thinks it owns it now,
	// and tries to delete it when the function returns, causing a segfault.
	hTotGammaE_H.GetListOfFunctions()->Clear();
	
	// Gamma Multiplicity
	// -------------------
	Log(toolName+" making gamma multiplicity stack",v_debug,verbosity);
	// capture on H
	TH1D hNumGammas_H("hNumGammas_H", "Gamma Multiplicity (Capture on H);Num Gammas;Num Events",100,0,10);
	intree->Draw("Length$(gamma_energy[])>>hNumGammas_H","nuclide_daughter_pdg==100045");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hNumGammas_Gd_155("hNumGammas_Gd_155", "Gamma Multiplicity (Capture on Gd-155);Num Gammas;Num Events",100,0,10);
	intree->Draw("Length$(gamma_energy[])>>hNumGammas_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hNumGammas_Gd_157("hNumGammas_Gd_157", "Gamma Multiplicity (Capture on Gd-157);Num Gammas;Num Events",100,0,10);
	intree->Draw("Length$(gamma_energy[])>>hNumGammas_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hNumGammas_Stack("hNumGammas_Stack","Gamma Multiplicity by Capture Nucleus;Num Gammas;Num Events");
	hNumGammas_H.SetLineColor(kRed);
	hNumGammas_Gd_155.SetLineColor(kBlue);
	hNumGammas_Gd_157.SetLineColor(kMagenta);
	hNumGammas_Stack.Add(&hNumGammas_H);
	hNumGammas_Stack.Add(&hNumGammas_Gd_155);
	hNumGammas_Stack.Add(&hNumGammas_Gd_157);
	hNumGammas_Stack.Draw();
	StackLegend.Draw();
	hNumGammas_H.GetListOfFunctions()->Add(&StackLegend);
	hNumGammas_Stack.Write();
	hNumGammas_H.GetListOfFunctions()->Clear();
	
	// Gamma Energy Spectrum
	// ---------------------
	Log(toolName+" making gamma spectrum stack",v_debug,verbosity);
	// capture on H
	TH1D hGammaE_H("hGammaE_H", "Gamma Energy (Capture on H);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE_H","nuclide_daughter_pdg==100045");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hGammaE_Gd_155("hGammaE_Gd_155", "Gamma Energy (Capture on Gd-155);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hGammaE_Gd_157("hGammaE_Gd_157", "Gamma Energy (Capture on Gd-157);Gamma Energy [MeV];Num Events",100,0,10);
	intree->Draw("gamma_energy>>hGammaE_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hGammaE_Stack("hGammaE_Stack","Gamma Spectrum by Capture Nucleus;Gamma Energy [MeV];Num Events");
	hGammaE_H.SetLineColor(kRed);
	hGammaE_Gd_155.SetLineColor(kBlue);
	hGammaE_Gd_157.SetLineColor(kMagenta);
	hGammaE_Stack.Add(&hGammaE_H);
	hGammaE_Stack.Add(&hGammaE_Gd_155);
	hGammaE_Stack.Add(&hGammaE_Gd_157);
	hGammaE_Stack.Draw();
	StackLegend.Draw();
	hGammaE_H.GetListOfFunctions()->Add(&StackLegend);
	hGammaE_Stack.Write();
	hGammaE_H.GetListOfFunctions()->Clear();
	
	// restore stats box behaviour
	gStyle->SetOptStat(statsboxdefault);
	
	// Misc
	// ====
	friendTree->SetAlias("neutron_travel_dist",
	                     "sqrt(pow(neutron_longitudinal_travel,2.)+pow(neutron_perpendicular_travel,2.))");
	
	return 1;
}


int PlotNeutronCaptures::GetBranches(){
	int success = (
//	(myTreeReader.GetBranchValue("filename",filename))                         &&
//	(myTreeReader.GetBranchValue("water_transparency",water_transparency))     &&
//	(myTreeReader.GetBranchValue("skdetsim_version",skdetsim_version))         &&
//	(myTreeReader.GetBranchValue("tba_table_version",tba_table_version))       &&
//	(myTreeReader.GetBranchValue("entry_number",entry_number))                 &&
//	(myTreeReader.GetBranchValue("subevent_num",subevent_number))              &&
	(myTreeReader.GetBranchValue("primary_pdg",primary_pdg))                   &&
//	(myTreeReader.GetBranchValue("primary_energy",primary_energy))             &&
	(myTreeReader.GetBranchValue("primary_start_mom",primary_start_mom))       &&
//	(myTreeReader.GetBranchValue("primary_start_pos",primary_start_pos))       &&
//	(myTreeReader.GetBranchValue("primary_end_pos",primary_end_pos))           &&
//	(myTreeReader.GetBranchValue("nuclide_parent_pdg",nuclide_parent_pdg))     &&
//	(myTreeReader.GetBranchValue("nuclide_creation_pos",nuclide_creation_pos)) &&
//	(myTreeReader.GetBranchValue("nuclide_decay_pos",nuclide_decay_pos))       &&
	(myTreeReader.GetBranchValue("nuclide_daughter_pdg",nuclide_daughter_pdg)) &&
	(myTreeReader.GetBranchValue("neutron_start_pos",neutron_start_pos))       &&
	(myTreeReader.GetBranchValue("neutron_end_pos",neutron_end_pos))           &&
//	(myTreeReader.GetBranchValue("neutron_start_energy",neutron_start_energy)) &&
//	(myTreeReader.GetBranchValue("neutron_end_energy",neutron_end_energy))     &&
//	(myTreeReader.GetBranchValue("neutron_end_process",neutron_end_process))   &&
	(myTreeReader.GetBranchValue("gamma_energy",gamma_energy))                 &&
//	(myTreeReader.GetBranchValue("gamma_time",gamma_time))                     &&
	(myTreeReader.GetBranchValue("electron_energy",electron_energy))         //&&
//	(myTreeReader.GetBranchValue("electron_time",electron_time))
	);
	
	return success;
}

void PlotNeutronCaptures::ClearOutputTreeBranches(){
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

int PlotNeutronCaptures::WriteTree(){
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

// Produce pie chart of nuclei that captured neutrons
// ==================================================
std::unique_ptr<TPie> GeneratePieFromHisto(std::string histoname, int verbose){
	TH1F* histo = (TH1F*)gROOT->FindObject(histoname.c_str());
	if(histo==nullptr){
		std::cerr<<"GeneratePieFromHisto could not find histo "<<histoname<<std::endl;
		return nullptr;
	}
	return GeneratePieFromHisto(histo, verbose);
}

std::unique_ptr<TPie> GeneratePieFromHisto(TH1F* histo, int verbose){
	std::string histoname = std::string(histo->GetName());
	if(verbose) std::cout<<"creating pie chart from histo "<<histoname<<", which has "
						 <<histo->GetNbinsX()<<" bins with contents: "<<std::endl;
	std::vector< std::pair<std::string,float> > histbins;
	for(int bini=0; bini<histo->GetNbinsX(); bini++){
		TString binlabel = histo->GetXaxis()->GetBinLabel(bini+1);
		double binconts = histo->GetBinContent(bini+1);
		if(binconts<0.01) binconts = 0.0f;  // round floats. useful if the histo has been scaled.
		if(verbose && binconts!=0.0f) std::cout<<binlabel.Data()<<" : "<<binconts<<std::endl;
		if(binconts<0) std::cerr<<"error converting "<<histoname<<" to pie chart: bin "<<binlabel.Data()
			<<" has "<<binconts<<" entries!"<<std::endl;
		if(binconts!=0) histbins.emplace_back(binlabel.Data(),binconts);
	}
	
	auto thepie = std::unique_ptr<TPie>(new TPie(TString::Format("%sPie",histoname.c_str()), TString::Format("%s",histoname.c_str()), histbins.size()));
	
	for(size_t bini=0; bini<histbins.size(); bini++){
		std::pair<std::string,float> abin = histbins.at(bini);
		std::string thebinlabel = abin.first;
		float thebincontents = abin.second;
		if(thebincontents<0) std::cerr<<"error converting "<<histoname<<" to pie chart: bin "<<thebinlabel
			<<" has "<<thebincontents<<" entries!"<<std::endl;
		thepie->SetEntryVal(bini,thebincontents);  // NO +1 - TPie's have no underflow bin!
		thepie->SetEntryLabel(bini,thebinlabel.c_str());
	}
	return thepie;
}
