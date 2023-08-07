/* vim:set noexpandtab tabstop=4 wrap */
#include "PlotNeutronCaptures.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <map>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPie.h"
#include "TStyle.h"
#include "TLeaf.h"

PlotNeutronCaptures::PlotNeutronCaptures():Tool(){
	// get the name of the tool from its class name
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}

bool PlotNeutronCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);                        // how verbose to be
	std::string inputFile, friendFile;
	m_variables.Get("inputFile",inputFile);                        // name of a file directly
	m_variables.Get("friendFile",friendFile);                      // name of a friend file if rqd
	m_variables.Get("outputFile",outputFile);                      // output file to write
	
	// get the Trees
	// -------------
	TFile* infile = m_data->OpenFileForReading(inputFile);
	if(infile==nullptr || infile->IsZombie()) return false;
	
	inTree = (TTree*)infile->Get("eventtree");
	if(!inTree){
		Log(m_unique_name+" Couldn't find 'eventtree' in file "+inputFile,v_error,verbosity);
		return false;
	}
	
	// depending on source (NCaptInfo or TruthNeutronCaptures) we may need another TTree
	// with additional branches
	if(friendFile!=""){
		TFile* friendfile = m_data->OpenFileForReading(friendFile);
		if(friendfile==nullptr || friendfile->IsZombie()) return false;
		friendTree = (TTree*)friendfile->Get("ntree");
		if(!friendTree){
			Log(m_unique_name+" Couldn't find 'ntree' in file "+friendFile,v_error,verbosity);
			return false;
		}
	}
	
	outfile = m_data->OpenFileForWriting(outputFile);
	if(!outfile) return false;
	
	return true;
}

bool PlotNeutronCaptures::Execute(){
	return true;
}

bool PlotNeutronCaptures::Finalise(){
	
	// make and write out histograms
	Log(m_unique_name+" making histograms",v_debug,verbosity);
	MakeHistos();
	
	Log(m_unique_name+" cleanup",v_debug,verbosity);
	m_data->CloseFile(outfile);
	
	return true;
}

int PlotNeutronCaptures::MakeHistos(){
	
	outfile->cd();
	
	if(friendTree) inTree->AddFriend(friendTree);
	
	// ======================
	// cumulative plots
	// ======================
	Log(m_unique_name+" making aggregate plots",v_debug,verbosity);
	// neutron energy
	TH1D hNeutronE("hNeutronE","Neutron Energy;Neutron Energy [MeV];Num Events",100,0,25);
	inTree->Draw("neutron_start_energy>>hNeutronE");
	hNeutronE.Write();
	
	// neutron travel distance
	TH1D hNeutronTravelDist("hNeutronTravelDist","Neutron Travel Distance;Distance [cm];Num Events",100,0,125);
	inTree->Draw("neutron_travel_dist>>hNeutronTravelDist");
	hNeutronTravelDist.Write();
	
	// neutron travel distance vs energy
	TH2D hNeutronTravelDistVsEnergy("hNeutronTravelDistVsEnergy","Neutron Travel Distance Vs Energy;Energy [MeV];Distance [cm]",200,0,25,200,0,125);
	inTree->Draw("neutron_travel_dist:neutron_start_energy>>hNeutronTravelDistVsEnergy");
	hNeutronTravelDistVsEnergy.Write();
	
	// neutron travel time
	TH1D hNeutronTravelTime("hNeutronTravelTime","Neutron Travel Time;Time [us];Num Events",100,0,1200);
	inTree->Draw("(neutron_travel_time/1000.)>>hNeutronTravelTime");
	hNeutronTravelTime.Write();
	
	// gamma mulitiplicity
	TH1D hNumGammas("hNumGammas","Gamma Multiplicity (All Nuclides);Num Gammas Emitted;Num Events",100,0,30);
	inTree->Draw("neutron_n_gammas>>hNumGammas");
	hNumGammas.Write();
	
	// gamma energy
	TH1D hGammaE("hGammaE", "Gamma Energy (All Nuclides);Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("gamma_energy>>hGammaE");
	hGammaE.Write();
	
	// total gamma energy from the neutron capture
	TH1D hSumGammaE("hSumGammaE","Total Emitted Gamma Energy (All Nuclides);Sum of Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("neutron_tot_gammaE>>hSumGammaE");
	hSumGammaE.Write();
	
	// plot gamma emission time relative to neutron capture (parent lifetime)
	// (n.b. currently Geant4 does not actually simulate this)
	TH1D hGammaT("hGammaT", "Gamma Emission Time (All Nuclides);Gamma Emission Time [ns];Num Events",100,0,1800E3);
	inTree->Draw("(gamma_time-capt_t)>>hGammaT");
	hGammaT.Write();
	
	// electron mulitiplicity
	TH1D hNumElectrons("hNumElectrons","Electron Multiplicity (All Nuclides);Num Electrons Emitted;Num Events",100,0,30);
	inTree->Draw("neutron_n_electrons>>hNumElectrons");
	hNumElectrons.Write();
	
	// electron energy
	TH1D hElectronE("hElectronE", "Electron Energy (All Nuclides);Electron Energy [MeV];Num Events",100,0,10);
	inTree->Draw("electron_energy>>hElectronE");
	hElectronE.Write();
	
	// total electron energy from the neutron capture
	TH1D hSumElectronE("hSumElectronE","Total Emitted Electron Energy (All Nuclides);Sum of Electron Energy [MeV];Num Events",100,0,10);
	inTree->Draw("neutron_tot_electronE>>hSumElectronE");
	hSumElectronE.Write();
	
	// electron emission time (parent lifetime)
	TH1D hElectronT("hElectronT", "Electron Emission Time (All Nuclides);Electron Emission Time [ns];Num Events",100,0,1800E3);
	inTree->Draw("(electron_time-capt_t)>>hElectronT");
	hElectronT.Write();
	
	// total daughter multiplicity
	TH1D hNumDaughters("hNumDaughters", "Daughter Multiplicity;Num Daughters;Num Events",100,0,30);
	inTree->Draw("neutron_n_daughters>>hNumDaughters");
	hNumDaughters.Write();
	
	// total daughter energy
	TH1D hSumDaughterE("hSumDaughterE", "Total Daughter Energy;Total Energy [MeV];Num Events",100,0,10);
	inTree->Draw("neutron_tot_daughterE>>hSumDaughterE");
	hSumDaughterE.Write();
	
	// pie chart of capture nuclei
	Log(m_unique_name+" making pie chart",v_debug,verbosity);
	// since the daughter nuclides are integers, but are very sparse, we can't uses TTree::Draw
	// to get the set of relevant pdg values (they'll be binned). Just have to do it the old fashioned way.
	std::map<int,int> capture_nuclide_vs_count;
	TBranch* bp = inTree->GetBranch("nuclide_daughter_pdg");
	if(!bp){
		Log(m_unique_name+" Error! No branch 'nuclide_daughter_pdg'!",v_error,verbosity);
	} else {
		TLeaf* lf = (TLeaf*)bp->GetListOfLeaves()->At(0);
		int* nuclidepdgptr = (int*)(lf->GetValuePointer());
		for(int i=0; i<bp->GetEntries(); ++i){
			bp->GetEntry(i);
			capture_nuclide_vs_count[*nuclidepdgptr]++;
		}
	}
	Log(m_unique_name+" found "+toString(capture_nuclide_vs_count.size())+" unique nuclides",v_debug,verbosity);
	TPie pCaptureNuclidePdg = TPie("pCaptureNuclidePdg", "Captures by Nuclide",capture_nuclide_vs_count.size());
	// OK the TPie class sucks - Google sheets is a much better alternative, but one of the biggest
	// problems is that when we have two small slices next to each other, the labels don't have space so overlap.
	// So let's try to fix this by alternating large and small contents. To do this, we need a sorted version
	// of our bin counts. We can't just make a map with swapped key and value because we'd end up merging
	// all nuclides with the same capture count. We can use a vector of pairs.
	std::vector<std::pair<int, int>> bins_by_entries;
	for(auto&& anuclide : capture_nuclide_vs_count){
		bins_by_entries.emplace_back(anuclide.first, anuclide.second);
	}
	std::sort(bins_by_entries.begin(), bins_by_entries.end(),
	    [](std::pair<int,int>& a, std::pair<int,int>& b) -> bool { return a.second < b.second; });
	auto fwit=bins_by_entries.cbegin();
	auto bwit=bins_by_entries.crbegin();
	for(int i=0; i<bins_by_entries.size();++i){
		int anuclide; // pdg
		int count;
		if((i%2)==0){
			anuclide = fwit->first;
			count = fwit->second;
			++fwit;
		} else {
			anuclide = bwit->first;
			count = bwit->second;
			++bwit;
		}
		Log(m_unique_name+": nuclide "+toString(anuclide)+" ("+PdgToString(anuclide)+") had "+toString(count)+" entries",v_debug,verbosity);
		pCaptureNuclidePdg.SetEntryLabel(i,PdgToString(anuclide).c_str());
		pCaptureNuclidePdg.SetEntryVal(i,count);
	}
	// making it look nice
	/*
//	pCaptureNuclidePdg.SetAngularOffset(333); // rotaaate the piiie
	pCaptureNuclidePdg.SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
	pCaptureNuclidePdg.SetValueFormat("%4.0f");
	pCaptureNuclidePdg.SetPercentFormat("%3.0f");
	*/
	for(int i=0; i<pCaptureNuclidePdg.GetEntries(); ++i){
		pCaptureNuclidePdg.SetEntryRadiusOffset(i,0.02); // slightly exploded
		pCaptureNuclidePdg.SetEntryFillColor(i, i+2);    // for some reason default colors are all very similar
		pCaptureNuclidePdg.SetEntryLineStyle(i,0);       // no line
	}
	//pCaptureNuclidePdg.SetLabelFormat("%val"); // too much info and it all overlaps for small slices
	pCaptureNuclidePdg.SetLabelFormat("#splitline{%txt}{(%val)}"); // not sure this gets saved to file, but for reference
	// since we don't draw the labels (isotopes), its then necessary to add a legend, but this is much clearer.
	pCaptureNuclidePdg.SetCircle(0.5, 0.4702026, 0.3302274);
	pCaptureNuclidePdg.SetTextSize(0.03455766);
	// saving to file
	pCaptureNuclidePdg.Draw("goff");
	//TLegend* pieLegend = pCaptureNuclidePdg.MakeLegend(); // does this get saved with the pie? // who owns it?
	pCaptureNuclidePdg.Write();
	
	// ==============================
	// broken down by capture nucleus
	// ==============================
	Log(m_unique_name+" making total gamma energy stack",v_debug,verbosity);
	auto statsboxdefault = gStyle->GetOptStat();
	gStyle->SetOptStat(0); // turn off stats box; overlaps with legends
	// Total Gamma Energy
	// ------------------
	// capture on H
	TH1D hTotGammaE_H("hTotGammaE_H", "Total Gamma Energy (Capture on H);Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("neutron_tot_gammaE>>hTotGammaE_H","(nuclide_daughter_pdg==100045 || nuclide_daughter_pdg==1000010020)");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hTotGammaE_Gd_155("hTotGammaE_Gd_155", "Total Gamma Energy (Capture on Gd-155);Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("neutron_tot_gammaE>>hTotGammaE_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hTotGammaE_Gd_157("hTotGammaE_Gd_157", "Total Gamma Energy (Capture on Gd-157);Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("neutron_tot_gammaE>>hTotGammaE_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hTotGammaE_Stack("hTotGammaE_Stack","Total Gamma E by Capture Nucleus;Gamma Energy [MeV];Num Events");
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
	hTotGammaE_Stack.Draw("goff");
	StackLegend.Draw("goff");
	// add the legend to the list of functions so that it gets saved on Write call
	// a THStack doesn't have a list of functions, so we have to add it to a component histo
	hTotGammaE_H.GetListOfFunctions()->Add(&StackLegend);
	hTotGammaE_Stack.Write();
	// we need to remove it afterwards, though, otherwise the histograms thinks it owns it now,
	// and tries to delete it when the function returns, causing a segfault.
	hTotGammaE_H.GetListOfFunctions()->Clear();
	
	// Gamma Multiplicity
	// -------------------
	Log(m_unique_name+" making gamma multiplicity stack",v_debug,verbosity);
	// capture on H
	TH1D hNumGammas_H("hNumGammas_H", "Gamma Multiplicity (Capture on H);Num Gammas;Num Events",100,0,10);
	inTree->Draw("Length$(gamma_energy[])>>hNumGammas_H","(nuclide_daughter_pdg==100045 || nuclide_daughter_pdg==1000010020)");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hNumGammas_Gd_155("hNumGammas_Gd_155", "Gamma Multiplicity (Capture on Gd-155);Num Gammas;Num Events",100,0,10);
	inTree->Draw("Length$(gamma_energy[])>>hNumGammas_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hNumGammas_Gd_157("hNumGammas_Gd_157", "Gamma Multiplicity (Capture on Gd-157);Num Gammas;Num Events",100,0,10);
	inTree->Draw("Length$(gamma_energy[])>>hNumGammas_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hNumGammas_Stack("hNumGammas_Stack","Gamma Multiplicity by Capture Nucleus;Num Gammas;Num Events");
	hNumGammas_H.SetLineColor(kRed);
	hNumGammas_Gd_155.SetLineColor(kBlue);
	hNumGammas_Gd_157.SetLineColor(kMagenta);
	hNumGammas_Stack.Add(&hNumGammas_H);
	hNumGammas_Stack.Add(&hNumGammas_Gd_155);
	hNumGammas_Stack.Add(&hNumGammas_Gd_157);
	hNumGammas_Stack.Draw("goff");
	StackLegend.Draw("goff");
	hNumGammas_H.GetListOfFunctions()->Add(&StackLegend);
	hNumGammas_Stack.Write();
	hNumGammas_H.GetListOfFunctions()->Clear();
	
	// Gamma Energy Spectrum
	// ---------------------
	Log(m_unique_name+" making gamma spectrum stack",v_debug,verbosity);
	// capture on H
	TH1D hGammaE_H("hGammaE_H", "Gamma Energy (Capture on H);Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("gamma_energy>>hGammaE_H","(nuclide_daughter_pdg==100045 || nuclide_daughter_pdg==1000010020)");
	
	// capture on Gd-155 -> daughter nuclide Gd-156
	TH1D hGammaE_Gd_155("hGammaE_Gd_155", "Gamma Energy (Capture on Gd-155);Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("gamma_energy>>hGammaE_Gd_155","nuclide_daughter_pdg==1000641560");
	
	// capture on Gd-157 -> daughter nuclide Gd-158
	TH1D hGammaE_Gd_157("hGammaE_Gd_157", "Gamma Energy (Capture on Gd-157);Gamma Energy [MeV];Num Events",100,0,10);
	inTree->Draw("gamma_energy>>hGammaE_Gd_157","nuclide_daughter_pdg==1000641580");
	
	// Stack of all of them
	THStack hGammaE_Stack("hGammaE_Stack","Gamma Spectrum by Capture Nucleus;Gamma Energy [MeV];Num Events");
	hGammaE_H.SetLineColor(kRed);
	hGammaE_Gd_155.SetLineColor(kBlue);
	hGammaE_Gd_157.SetLineColor(kMagenta);
	hGammaE_Stack.Add(&hGammaE_Gd_155);
	hGammaE_Stack.Add(&hGammaE_Gd_157);
	hGammaE_Stack.Add(&hGammaE_H); // add H last so that the big H peak is on top of gd continuum
	hGammaE_Stack.Draw("goff");
	StackLegend.Draw("goff");
	hGammaE_H.GetListOfFunctions()->Add(&StackLegend);
	hGammaE_Stack.Write();
	hGammaE_H.GetListOfFunctions()->Clear();
	
	// restore stats box behaviour
	gStyle->SetOptStat(statsboxdefault);
	
	// clean up the friend relationship (apparently it gets written to file...
	if(friendTree){
		inTree->RemoveFriend(friendTree);
	}
	
	return 1;
}

