/* vim:set noexpandtab tabstop=4 wrap */
#include "FitLi9Lifetime.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"
#include "MTreeReader.h"
#include "MTreeSelection.h"

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

// For Fitting
#include "Fit/Fitter.h"
//#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
//#include "Fit/Chi2FCN.h"
//#include "Fit/DataOptions.h"
#include "Fit/FitConfig.h"

// For defining the functions
//#include "TList.h"
#include "Math/WrappedMultiTF1.h"
//#include "HFitInterface.h"

FitLi9Lifetime::FitLi9Lifetime():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

// from 2015 paper Table I
const double li9_lifetime_secs = 0.26;
const double li9_endpoint = 14.5; // MeV

bool FitLi9Lifetime::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("li9_lifetime_dtmin",li9_lifetime_dtmin);
	m_variables.Get("li9_lifetime_dtmax",li9_lifetime_dtmax);
	m_variables.Get("outputFile",outputFile);          // where to save data. If empty, current TFile
	m_variables.Get("treeReaderName",treeReaderName);
	
	myTreeReader = m_data->Trees.at(treeReaderName);
	myTreeSelections = m_data->Selectors.at(treeReaderName);
	
	return true;
}


bool FitLi9Lifetime::Execute(){
	
	GetBranchValues();
	
	// the following cuts are based on muon-lowe pair variables, so loop over muon-lowe pairs
	std::set<size_t> spall_mu_indices = myTreeSelections->GetPassingIndexes("dlt_mu_lowe>200cm");
	Log(toolName+" Looping over "+toString(spall_mu_indices.size())
				+" preceding muons to look for spallation events",v_debug,verbosity);
	for(size_t mu_i : spall_mu_indices){
		// check whether this passed the additional Li9 cuts
		if(not myTreeSelections->GetPassesCut("ntag_FOM>0.995")) continue;
		
		// plot distribution of beta energies from passing triplets, compare to fig 4
		Log(toolName+" filling li9 candidate distributions",v_debug+2,verbosity);
		li9_e_vals.push_back(LOWE->bsenergy); // FIXME weight by num_post_muons
		
		// plot distirbution of mu->beta   dt from passing triplets, compare to fig 6
		li9_muon_dt_vals.push_back(fabs(dt_mu_lowe[mu_i])); // FIXME weight by num_post_muons
	}
	
	return true;
}

bool FitLi9Lifetime::GetBranchValues(){
	// retrieve variables from branches
	bool success = 
	(myTreeReader->Get("LOWE", LOWE)) &&
	(myTreeReader->Get("spadt", dt_mu_lowe));
	
	return success;
}


bool FitLi9Lifetime::Finalise(){
	
	// make a new file if given a filename, or if blank check there is a valid file open
	TFile* fout = nullptr;
	if(outputFile!=""){
		fout = new TFile(outputFile.c_str(),"RECREATE");
		fout->cd();
	} else {
		if(gDirectory->GetFile()==nullptr){
			Log(toolName+" Error! No output file given and no file open!",v_error,verbosity);
			return true;
		}
	}
	
	Log(toolName+" making plot of Li9 candidates Dt distribution",v_debug,verbosity);
	PlotLi9LifetimeDt();
	
	Log(toolName+" making plot of Li9 candidates beta energy distribution",v_debug,verbosity);
	PlotLi9BetaEnergy();
	
	if(fout!=nullptr){
		fout->Close();
		delete fout;
		fout=nullptr;
	}
	
	return true;
}

// =========================================================================
// Li9 lifetime fits
// =========================================================================

bool FitLi9Lifetime::PlotLi9LifetimeDt(){
	
	// make a histogram to bin the data
	std::cout<<"making li9 mu-lowe dt histogram"<<std::endl;
	TH1F li9_muon_dt_hist("li9_muon_dt_hist","Muon to Low-E dt for Li9 triplets",
	                       15,li9_lifetime_dtmin,li9_lifetime_dtmax);
	for(auto&& aval : li9_muon_dt_vals){
		li9_muon_dt_hist.Fill(aval);
	}
	std::cout<<"saving to file"<<std::endl;
	li9_muon_dt_hist.Write();
	
	/*
	// TODO calculate the efficiency of background selection and scale down the number
	// of spallation events to get the expected number of accidental triplets
	// efficiency is based on background acceptance of BDT and of preceding cuts i guess
	// this is subtracted off to estimate the number of li9 events
	
	// get ntag signal and background efficiency
	// ------------------------------------------
	// read the BDT ROC and populate the splines
	fill_roc_ntag();
	// evaluate the splines at the working point used to retrieve the efficiencies
	double ncut_li9=0.995; // this should be read from the cut threshold used in PurewaterLi9Rate
	double ntag_sig_eff = sig->Eval(ncut_li9);
	double ntag_bg_eff = bg->Eval(ncut_li9);
	
	// total background selection efficiency is product of all cuts specific to Li9
	// * dt in li9 lifetime range 0.05-0.5
	// * lowe energy in li9 beta energy range
	// * no other muons within 1ms
	// * ntag bdt fom > 0.995
	// how do we determine the fraction of non-li9 events that pass these? MC? paper doesn't say...
	double li9_lifetime_bg_eff = 1;
	double li9_energy_bg_eff = 1;
	double no_mu_1ms_bg_eff = 1;
	double li9_bg_eff = li9_lifetime_bg_eff * li9_energy_bg_eff * ntag_bg_eff * no_mu_1ms_bg_eff;
	// we multiply this by the number of spallation events to get the total bg contamination
	double num_spall_events = dt_mu_lowe_vals.size();
	double total_li9_bg = num_spall_events * li9_bg_eff;
	// and then we get the number of li9 events as the remainder
	// fig 4 plots the spectrum of li9 energy, and of accidentals - presumably from MC of all bgs?
	// and the expected li9 spectrum - presumably from some suitable reference
	*/
	
	// fit the Li9 lifetime, Fig 6 from the paper. NOT neutron capture times (plot 5!)
	// this fits a combination of 8 exponentials (including li9) with fixed abundances,
	// and uses it to extract the lifetime of li9. The paper only does a binned fit.
	// do we do chi2 fit or binned likelihood fit? why? FIXME how do we do binned likelihood?
	
	// perform a binned chi2 fit
	std::cout<<"doing Li9 lifetime binned chi2 fit"<<std::endl;
	double binned_estimate = BinnedLi9DtChi2Fit(&li9_muon_dt_hist);
	
	// TODO (un)binned extended likelihood fit to also extract number of Li9 events?
	
	return true;
}

double FitLi9Lifetime::BinnedLi9DtChi2Fit(TH1F* li9_muon_dt_hist){
	// this fits the lifetime of li9 as a cross-check, by fitting 8 exponentials
	// based on backgrounds from other isotopes. The amount of other isotopes is based
	// on the previous global fit to muon-lowe dt, and the fraction of those that pass
	// all cuts, from first red, third, lt, li9 dt, li9 energy and ntag
	// To do this we need the efficiency of all of those cuts for all isotopes!
	// finally we fit (unbinned? binned?) the distribution of mu-lowe times leaving only li9 lifetime floating
	
	// currently the above paragraph does not describe what this code does FIXME
	std::cout<<"making TF1 for binned chi2 fit with "<<li9_muon_dt_hist->GetEntries()<<" values"<<std::endl;
	
	// number of li9 left after time t = N = N0*exp(-t/τ)
	// rate of change of li9, i.e. rate of observed li9 events = dN/dt = -N0*τ*exp(-t/τ)
	// we can ignore the sign which just says N is decreasing.
	// simple chi2 fit of 'y = C + A*τ*exp(-dt/τ)' for all values of C, A, τ
	TF1 li9_muon_dt_func("li9_muon_dt_func","[0]+[1]*[2]*exp(-x/[2])",
	                      li9_lifetime_dtmin,li9_lifetime_dtmax);
	
	// set parameter names
	li9_muon_dt_func.SetParName(0,"num bg events");
	li9_muon_dt_func.SetParName(1,"num Li9 events");
	li9_muon_dt_func.SetParName(2,"Li9 lifetime");
	
	// set starting values - we can fix the lifetime
	//li9_muon_dt_func.SetParameters(allparams.data());  // pass an array or set individually
	li9_muon_dt_func.SetParameter(0,0);  // TODO estimate from accidental rate of ntag from MC and num events?
	li9_muon_dt_func.SetParameter(1,li9_muon_dt_hist->GetBinContent(1));
	li9_muon_dt_func.FixParameter(2,li9_lifetime_secs);
	
	// set parameter limits. For minuit at least this is strongly discouraged
	// Maybe ok for chi2 fits but not likelihood fits?
	//li9_muon_dt_func.SetParLimits(4, 0., 200.);
	
	// set num TF1 points for drawing
	li9_muon_dt_func.SetNpx(1000);
	
	// DO THE FIT
	std::cout<<"invoking li9 lifetime chi2 fit"<<std::endl;
	TFitResultPtr fitresult = li9_muon_dt_hist->Fit("li9_muon_dt_func","MRS");
	// options here are  M: better fit, R: use range, S: get resultsptr
	
	// print result
	std::cout<<"li9 lifetime binned chi2 fit parameters = {";
	for(int i=0; i<(li9_muon_dt_func.GetNpar()); i++){
		std::cout<<li9_muon_dt_func.GetParameter(i);
		(i<(li9_muon_dt_func.GetNpar()-1)) ? std::cout<<", " : std::cout<<"};\n";
	}
	//float fitchi2 = fitresult->Chi2();                               // same as below, which
	float fitchi2 = li9_muon_dt_func.GetChisquare();         // doesn't need fitresultptr
	Log(toolName+" li9 mu->lowe dt fit chi2 was "+toString(fitchi2),v_message,verbosity);
	
	// draw result
	/*
	li9_muon_dt_hist->Draw();
	gPad->WaitPrimitive();
	gPad->Clear();
//	//fit->Draw("lsame");   // not necessary, fit is added to histogram's list of functions and drawn automatically
	*/
	li9_muon_dt_hist->GetListOfFunctions()->Clear();
	
	std::cout<<"li9 lifetime binned chi2 fit done"<<std::endl;
	return li9_muon_dt_func.GetParameter(1);
}

// =========================================================================
// Li9 energy spectrum fits
// =========================================================================

bool FitLi9Lifetime::PlotLi9BetaEnergy(){
	
	// plot the distribution to compare to paper Fig 4
	TH1F li9_e_hist("li9_e_hist","Li9 Candidate Beta Energy",7,6,li9_endpoint);
	// add rest mass as reconstructed energy is only kinetic...
	double e_rest_mass = 0.511; // [MeV] TODO replace this with TParticleDatabase lookup
	for(auto&& aval : li9_e_vals) li9_e_hist.Fill(aval+e_rest_mass); // FIXME weight by # post mus & num neutrons
	li9_e_hist.Write();
	
	// to overlay the expected li9 and background plots, we need to take the beta spectra
	// of each isotope and propagate it through the efficiency chain TODO
	
	return true;
}
