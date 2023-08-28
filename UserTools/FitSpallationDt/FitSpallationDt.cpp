/* vim:set noexpandtab tabstop=4 wrap */
#include "FitSpallationDt.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"
#include "MTreeReader.h"
#include "MTreeSelection.h"

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

FitSpallationDt::FitSpallationDt():Tool(){}

const double fiducial_vol = 22.5;               // kton, from paper
const double paper_livetime = 1890;             // days, from paper
const double paper_first_reduction_eff = 99.;   // fraction, from paper
const double paper_dlt_cut_eff = 78.8;          // fraction, from paper

const std::map<std::string,double> lifetimes{  // from 2015 paper, seconds
	{"11Be",19.9},
	{"16N",10.3},
	{"15C",3.53},
	{"8Li",1.21},
	{"8B",1.11},
	{"16C",1.08},
	{"9Li",0.26},
	{"9C",0.18},
	{"8He",0.17},
	{"12Be",0.034},
	{"12B",0.029},
	{"13B",0.025},
	{"14B",0.02},
	{"12N",0.016},
	{"13O",0.013},
	{"11Li",0.012},
	{"8He_9C",0.175},
	{"8Li_8B",1.16},
	{"ncapture",204.8E-6}
};

const std::map<std::string,double> papervals{  // from 2015 paper table II, Rate in events/kton/day
	// NOTE: these have been corrected for the efficiency of 6MeV energy threshold,
	// so are the rate of ALL DECAYS, not just the rate of decays that produce betas with >6MeV!
	{"11Be",5.7},   // read off plot, 16.9 upper limit in table
	{"16N",39.7},
	{"15C",3.5},    // read off plot, 6.7 upper limit in table
	{"8Li",3.9},    // read off plot, combined in table but plot splits them
	{"8B",4.9},     // read off plot, combined in table but plot splits them
//	{"8Li_8B",8.3},
	{"16C",0},
	{"9Li",0.9},
	//{"9C",0.1},   // read off plot, 0.7 upper limit in table
	//{"8He",0.2},  // read off plot, 0.7 upper limit in table (combined in table but plotted separately?)
	{"8He_9C",0.3}, // read off plot, 1.4 upper limit
	{"12Be",0},
	{"12B",19.8},
	{"13B",0},
	{"14B",0},
	{"12N",2.8},
	{"13O",0},
	{"11Li",0},
	{"const",330}   // not in table, use value from plot, uhhhhh units? Events/0.006s/FV?
};

const std::map<std::string,double> papereffs{
	{"11Be",48.84},
	{"16N",57.68},
	{"15C",40.76},
	{"8Li",54.86},
	{"8B",65.76},
	{"16C",0},
	{"9Li",50.25},
	{"9C",64.35},
	{"8He",28.84},
	{"12Be",0},
	{"12B",58.32},
	{"13B",0},
	{"14B",0},
	{"12N",72.04},
	{"13O",0},
	{"11Li",0}
};

bool FitSpallationDt::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);              // how verbose to be
	m_variables.Get("outputFile",outputFile);            // where to save data. If empty, current TFile
	// hacky way to scale the total number of events when comparing to the paper plots
	m_variables.Get("paper_scaling",paper_scaling);      // overall scaling factor of paper
	m_variables.Get("livetime",livetime);                // for when we don't have it from upstream
	
	// this tool may build its own dataset with an upstream ROOT file reader
	m_variables.Get("treeReaderName",treeReaderName);    // reader when getting entries from an upstream source
	// or we can shortcut it by dumping the dataset to a BoostStore and then retrieving in one go.
	m_variables.Get("valuesFileMode",valuesFileMode);    // if not provided by upstream tool
	m_variables.Get("valuesFile",valuesFile);            // if not provided by upstream tool
	
	// for when we're reading files without pre-selection
	m_variables.Get("run_min",run_min);                  // debug, for when we're loading the dt data directly
	m_variables.Get("run_max",run_max);                  // debug, for when we're loading the dt data directly
	
	// various alterations on the fit process while we try to find the source of our discrepancy
	// against previous versions of this study
	m_variables.Get("useHack",useHack);                  // constrain abundances to at least half paper vals
	m_variables.Get("n_dt_bins",n_dt_bins);              // num bins in the dt histogram
	m_variables.Get("binning_type",binning_type);        // '0 = logarithmic' or '1 = linear'
	m_variables.Get("random_subtract",random_subtract);  // subtract pre-mu dt dist from post-mu dt dist before fitting
	m_variables.Get("laurasfile",laurasfile);            // file with lauras dt histogram to fit
	m_variables.Get("use_par_limits",use_par_limits);    // whether to limit isotope abundances to >0 (and <1E7)
	m_variables.Get("fix_const",fix_const);              // whether to include a constant term in the fitting
	m_variables.Get("split_iso_pairs",split_iso_pairs);  // whether to split pairs (e.g. 8Be_8Li) into
	// two expontial terms with a shared amplitude, i.e. (A/2)*{exp(-t/t1)+exp(-t/t2)}
	// or combine them into one term with an average lifetime, i.e. A*(exp(-t/{(t1+t2)*0.5}))
	
	// energy threshold efficiencies, from FLUKA
	m_variables.Get("efficienciesFile",efficienciesFile);
	
	// read efficiencies of the different thresholds of old vs new data
	GetEnergyCutEfficiencies();
	
	// if we're loading data with an upstream ROOT file reader, retrieve the reader
	if(valuesFileMode!="read"){
		myTreeReader = m_data->Trees.at(treeReaderName);
		myTreeSelections = m_data->Selectors.at(treeReaderName);
	}
	
	return true;
}


bool FitSpallationDt::Execute(){
	
	// if retrieving a previously built dataset from a BoostStore, we can skip the Execute loop
	if(valuesFileMode=="read"){
		m_data->vars.Set("StopLoop",1);
		return true;
	}
	
	// this tool supports two types of upstream input file
	// the normal analysis using the 2020 spallation dataset directly
	if(myTreeSelections!=nullptr){
		Analyse();
	} else {
		// or, as part of debugging, files from laura's processing script
		Analyse_Laura();
	}
	
	return true;
}

bool FitSpallationDt::Analyse(){
	// retrieve desired branches
	get_ok = GetBranchValues();
	
	// the following cuts are based on muon-lowe pair variables, so loop over muon-lowe pairs
	std::set<size_t> spall_mu_indices = myTreeSelections->GetPassingIndexes("dlt_mu_lowe>200cm");
	Log(m_unique_name+" Looping over "+toString(spall_mu_indices.size())
				+" preceding muons to look for spallation events",v_debug,m_verbose);
	for(size_t mu_i : spall_mu_indices){
		// record the distribution of dt_mu_lowe
		Log(m_unique_name+" filling mu->lowe dt distribution for spallation entries",v_debug+2,m_verbose);
		// this *should* only contain pre muons with dt < 0:
		if(dt_mu_lowe[mu_i]>0){
			Log(m_unique_name+" error! Spallation muon with time "+toString(dt_mu_lowe[mu_i])+" after lowe event!",
				v_warning,m_verbose);
			continue;
		}
		dt_mu_lowe_vals.push_back(dt_mu_lowe[mu_i]);  // FIXME weight by num_pre_muons
	}
	return true;
}

bool FitSpallationDt::Analyse_Laura(){
	// retrieve desired branches
	get_ok = GetBranchValuesLaura();
	
	// ensure selections are aligned with our analysis for comparison
//	if(dt<-60) return true;
//	if(bsgood_relic<0.5) return true;
//	if(energy_relic<8) return true;
//	if(lt>200) return true;
//	if(muboy_status!=1) return true;
	if((nrunsk<run_min) || (nrunsk>run_max)) return true;
	dt_mu_lowe_vals.push_back(dt);
	return true;
}

bool FitSpallationDt::GetBranchValues(){
	bool success = (
		(myTreeReader->Get("spadt",dt_mu_lowe))
	);
	return success;
}

bool FitSpallationDt::GetBranchValuesLaura(){
	bool success = (
		(myTreeReader->Get("muboy_status",muboy_status)) &&   // for paper version of lt cut
		(myTreeReader->Get("lt",lt)) &&
		(myTreeReader->Get("energy_relic",energy_relic)) &&
		(myTreeReader->Get("bsgood_relic",bsgood_relic)) &&
		(myTreeReader->Get("nrunsk",nrunsk)) &&
		(myTreeReader->Get("dt",dt))
	);
	return success;
}

bool FitSpallationDt::Finalise(){
	
	// if we want to shortcut the file read loop and go straight to finalise,
	// dump any necessary variables to an output file now. Or, if we're skipping
	// the file loop, read that file in now.
	if(valuesFileMode=="write"){
		// set all the values into the BStore
		BStore valueStore(true);
		valueStore.Initnew(valuesFile,uncompressed,true);
		valueStore.Set("livetime",livetime);
		valueStore.Set("dt_mu_lowe_vals",dt_mu_lowe_vals);
		// save BoostStore
		valueStore.Save();
		valueStore.Close(); // necessary to complete the file write!
	} else if(valuesFileMode=="read"){
		BStore valueStore(true);
		valueStore.Initnew(valuesFile,uncompressed,true);
		valueStore.Get("dt_mu_lowe_vals",dt_mu_lowe_vals);
		valueStore.Get("livetime",livetime);
	} else {
		// get any remaining variables from upstream tools
		get_ok = m_data->CStore.Get("livetime",livetime);
	}
	
	Log(m_unique_name+"fitting "+toString(dt_mu_lowe_vals.size())+" spallation dt values",v_debug,m_verbose);
	
	// make a new output file for fit results if given a filename
	// or if no file name is given, check there is a valid ROOT file open, and if so we'll write to it
	TFile* fout = nullptr;
	if(outputFile!=""){
		fout = new TFile(outputFile.c_str(),"RECREATE");
		fout->cd();
	} else {
		if(gDirectory->GetFile()==nullptr){
			Log(m_unique_name+" Error! No output file given and no file open!",v_error,m_verbose);
			return true;
		}
	}
	
	// Fit the data, extract the abundances
	Log(m_unique_name+" Fitting spallation dt distribution",v_debug,m_verbose);
	PlotSpallationDt();
	
	if(fout!=nullptr){
		fout->Close();
		delete fout;
		fout=nullptr;
	}
	
	return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// =====================================================================

bool FitSpallationDt::PlotSpallationDt(){
	/* Main driver function for this analysis. This function fits the distribution of
	   spallation muon to low-e event time differences, using several intermediate fits
	   to obtain initial values, before one final fit with everything floating.
	   The data is fit with a sum of exponentials, one for each contributing spallation isotope.
	   The number of events from each isotope is extracted from the fit and converted
	   (by accounting for livetime, fiducial volume, and selection efficiency)
	   into rate of production of that spallation isotope.
	*/
	
	// dt distribution of muons passing dlt<200cm cut
	// first a version using linear binning, both for spallation candidates (mu-lowe dt<0)
	// and random events (mu-lowe dt>0)
	// 5000 bins (30/0.006) gives chi2/NDOF that matches the 2015 paper.
	int nbins = n_dt_bins;
	binwidth=30./nbins;
	TH1F dt_mu_lowe_hist("dt_mu_lowe_hist","Spallation Muon to Low-E Time Differences",nbins,0,30);
	TH1F dt_mu_lowe_hist_short("dt_mu_lowe_hist_short","Spallation Muon to Low-E Time Differences",500,0,0.25);
	TH1F dt_mu_lowe_rand_hist("dt_mu_lowe_rand_hist","Random Muon to Low-E Time Differences",nbins,0,30);
	for(auto&& aval : dt_mu_lowe_vals){
		if(aval>0){
			dt_mu_lowe_rand_hist.Fill(aval);
		} else {
			dt_mu_lowe_hist.Fill(fabs(aval));
			dt_mu_lowe_hist_short.Fill(fabs(aval));
		}
	}
	dt_mu_lowe_hist.Write();
	dt_mu_lowe_hist_short.Write();
	dt_mu_lowe_rand_hist.Write();
	
	// we can either fit with a constant term to account for (dt-independent) random contamination,
	// or subtract the distribution of random candidates from the spallation candidates
	TH1F* dt_mu_lowe_hist_randsub = (TH1F*)dt_mu_lowe_hist.Clone("dt_mu_lowe_hist_randsub");
	dt_mu_lowe_hist_randsub->Add(&dt_mu_lowe_rand_hist,-1);
	
	// alternative versions of the above with equally spaced bins on a log scale
	// we need to fit the same histogram (with the same binning) for all the intermediate fits,
	// otherwise the bin widths change, the contents change, and the fit parameters change.
	// We can only really achieve suitable binning across the whole dt range with logarithmic binning.
	std::vector<double> binedges = MakeLogBins(0.001, 30, nbins+1);
	TH1F dt_mu_lowe_hist_log("dt_mu_lowe_hist_log","Data;dt(s);Events/0.006 s",
							 nbins, binedges.data());
	TH1F dt_mu_lowe_rand_hist_log("dt_mu_lowe_rand_hist_log","Data;dt(s);Events/0.006 s",
							 nbins, binedges.data());
	for(auto&& aval : dt_mu_lowe_vals){
		if(aval>0) dt_mu_lowe_rand_hist_log.Fill(aval);
		else       dt_mu_lowe_hist_log.Fill(fabs(aval));
	}
	// Since we used different bin widths, to have a consistent y axis
	// (events per fixed time interval) we need to scale each bin's contents by its bin width
	for(int bini=1; bini<dt_mu_lowe_hist_log.GetNbinsX()+1; ++bini){
		dt_mu_lowe_rand_hist_log.SetBinContent(bini,
			dt_mu_lowe_rand_hist_log.GetBinContent(bini)/dt_mu_lowe_rand_hist_log.GetBinWidth(bini));
		dt_mu_lowe_hist_log.SetBinContent(bini,
			dt_mu_lowe_hist_log.GetBinContent(bini)/dt_mu_lowe_hist_log.GetBinWidth(bini));
	}
	// for some reason the 2015 paper plots things in events per 0.006s,
	// rather than events per second, even though our x-axis is in units of seconds!
	// To match it, scale our y-axis (which is in units of events per second) by *0.006s
	dt_mu_lowe_hist_log.Scale(30./nbins);
	dt_mu_lowe_rand_hist_log.Scale(30./nbins);
	
	// random-subtracted log-binned histogram
	TH1F* dt_mu_lowe_hist_log_randsub = (TH1F*)dt_mu_lowe_hist_log.Clone("dt_mu_lowe_hist_log_randsub");
	dt_mu_lowe_hist_log_randsub->Add(&dt_mu_lowe_rand_hist_log,-1);
	
	// select which histogram to fit here:
	TH1* the_hist_to_fit = nullptr;
	if(binning_type==0){                                                    // log binned
		if(random_subtract) the_hist_to_fit = dt_mu_lowe_hist_log_randsub;  // random subtracted
		else the_hist_to_fit = &dt_mu_lowe_hist_log;                        // not random subtracted
	} else {                                                                // linear binned
		if(random_subtract) the_hist_to_fit = dt_mu_lowe_hist_randsub;      // random subtracted
		else the_hist_to_fit = &dt_mu_lowe_hist;                            // not random subtracted
	}
	// one further debug option: fit the histogram of dt values from laura's script directly
	if(laurasfile!=""){
		// pull laura's histogram from her file
		auto current_dir = gDirectory;
		TFile* f = TFile::Open(laurasfile.c_str());
		f->cd();
		TDirectoryFile* spal1 = (TDirectoryFile*)f->Get("spal1");
		spal1->cd();
		the_hist_to_fit = (TH1D*)gDirectory->Get("data_dt-random1_dt");
		Log(m_unique_name+": Fitting dt histogram from file "+laurasfile,v_error,m_verbose);
		current_dir->cd();
	}
	
//	std::cout<<"this plot will be used for all the time distribution fits:"<<std::endl;
//	the_hist_to_fit->Draw();
//	gPad->WaitPrimitive();
	
	// Now we have our histogram, fit it!
	// we do the fitting in 5 stages, initially fitting sub-ranges of the distribution
	for(int i=0; i<5; ++i) FitDtDistribution(*the_hist_to_fit, i);
	
	// the production rate integrated over the whole energy range is given by:
	// Ri = Ni / (FV * T * eff_i)
	// where FV is the fid vol, T is the live time, and eff_i is the total efficiency
	// of the reduction cuts for retaining events from isotope i.
	// (note this is using the paper definiton, not zhangs definition. just changes defn of Ni<->Ni*eff_i)
	std::map<std::string,double> rates;
	for(auto&& anisotope : fit_amps){
		std::string isotope = anisotope.first;
		if(isotope.substr(0,5)=="const") continue;
		// we need to know the efficiency of the selection, which is obtained
		// from the fraction of beta spectrum that's below 8MeV
		if(isotope.find("_")==std::string::npos){
			double energy_cut_eff = reco_effs_8mev.at(anisotope.first)/100.;
			// FIXME for now just assume the same 1st reduction and dlt cut efficiency as 2015 paper
			double efficiency = energy_cut_eff * (paper_first_reduction_eff/100.) * (paper_dlt_cut_eff/100.);
			rates[anisotope.first] = anisotope.second/(efficiency * fiducial_vol * livetime);
			std::cout<<"Num of "<<anisotope.first<<" events is "<<anisotope.second
					 <<", energy cut efficiency is "<<energy_cut_eff<<", total efficiency is "<<efficiency
					 <<" giving a total rate of "<<rates[anisotope.first]<<"/kton/day"<<std::endl;
		} else {
			std::string first_isotope = isotope.substr(0,isotope.find_first_of("_"));
			std::string second_isotope = isotope.substr(isotope.find_first_of("_")+1,std::string::npos);
			double e_cut_eff_1 = reco_effs_8mev.at(first_isotope)/100.;
			double e_cut_eff_2 = reco_effs_8mev.at(second_isotope)/100.;
			// FIXME for now just assume the same 1st reduction and dlt cut efficiency as 2015 paper
			double first_efficiency = e_cut_eff_1 * (paper_first_reduction_eff/100.) * (paper_dlt_cut_eff/100.);
			double second_efficiency = e_cut_eff_2 * (paper_first_reduction_eff/100.) * (paper_dlt_cut_eff/100.);
			rates[first_isotope] = 0.5*anisotope.second/(first_efficiency * fiducial_vol * livetime);
			rates[second_isotope] = 0.5*anisotope.second/(second_efficiency * fiducial_vol * livetime);
			std::cout<<"efficiency of "<<first_isotope<<" is "<<first_efficiency
					 <<", calculated rate of production is "<<rates[first_isotope]
					 <<"/kton/day"<<std::endl
					 <<", efficiency of "<<second_isotope<<" is "<<second_efficiency
					 <<", calculated rate of production is "<<rates[second_isotope]
					 <<"/kton/day"<<std::endl;
		}
	}
	// TODO save this map to an ouput file
	// not worth doing till we re-run and determine the livetime
	
	// the yield across the whole energy range is obtained from the fit amplitude via:
	// Yi = Ni / (Rmu * T * rho * Lmu)
	// where Rmu is the muon rate, T the live time, rho the density of water and Lmu the
	// measured path length of the muon track...? is this event-wise?
	
	return true;
}

std::vector<double> FitSpallationDt::MakeLogBins(double xmin, double xmax, int nbins){
	std::vector<double> binedges(nbins+1);
	double xxmin=log10(xmin);
	double xxmax = log10(xmax);
	for(int i=0; i<nbins; ++i){
		binedges[i] = pow(10,xxmin + (double(i)/(nbins-1.))*(xxmax-xxmin));
	}
	binedges[nbins+1] = xmax; // required
	return binedges;
}

bool FitSpallationDt::FitDtDistribution(TH1& dt_mu_lowe_hist, int rangenum){
	/* Do dt fitting based on section B of the 2015 paper.
	   As described in section B2, we do 4 fits to subsets of the time range,
	   making note of the fit results as we go, then one final fit to the complete time range
	   using our previous results as a starting point for each isotope's abundance.
	*/
	
	switch (rangenum){
	case 0:{
		// fit range 50us -> 0.1s with 12B + 12N only
		// make the function which we'll fit
		Log(m_unique_name+"calling BuildFunction for time range case "+toString(rangenum)
			+", 12B+12N",v_debug,m_verbose);
		TF1 func_sum = BuildFunction({"12B","12N"},50e-6,0.1);
		// do the fit
		Log(m_unique_name+" doing case "+toString(rangenum)+" fit",v_debug,m_verbose);
		std::string fitopt = (m_verbose>2) ? "Rq" : "R";
		dt_mu_lowe_hist.Fit(&func_sum,fitopt.c_str(),"",50e-6,0.1);
		// record the results for the next step
		Log(m_unique_name+" recording results from case "+toString(rangenum)+" fit",v_debug,m_verbose);
		PushFitAmp(func_sum,"12B");
		PushFitAmp(func_sum,"12N");
		fit_amps["const_0"] = func_sum.GetParameter("const");
		if(m_verbose){
			std::cout<<"recorded results were: "<<std::endl
					 <<"12B: "<<fit_amps["12B"]<<std::endl
					 <<"12N: "<<fit_amps["12N"]<<std::endl
					 <<"const_0: "<<fit_amps["const_0"]<<std::endl;
		}
		
		// draw the fit, if we want to check it
		/*
		dt_mu_lowe_hist.Draw();
		gPad->WaitPrimitive();
		*/
		
		// remove this intermediate fit from the list of functions
		// associated with our dt distribution, since we no longer need it
		dt_mu_lowe_hist.GetListOfFunctions()->Clear();
		
		break;
		}
	case 1:{
		// fit range 6-30s with 16N + 11B only
		Log(m_unique_name+"calling BuildFunction for time range case "+toString(rangenum)
			+", 16N+11Be",v_debug,m_verbose);
		TF1 func_sum = BuildFunction({"16N","11Be"},6,30);
		Log(m_unique_name+" doing case "+toString(rangenum)+" fit",v_debug,m_verbose);
		std::string fitopt = (m_verbose>2) ? "Rq" : "R";
		dt_mu_lowe_hist.Fit(&func_sum,fitopt.c_str(),"",6,30);
		// record the results for the next step
		Log(m_unique_name+" recording results from case "+toString(rangenum)+" fit",v_debug,m_verbose);
		PushFitAmp(func_sum,"16N");
		PushFitAmp(func_sum,"11Be");
		fit_amps["const_1"] = func_sum.GetParameter("const");
		if(m_verbose){
			std::cout<<"recorded results were: "<<std::endl
					 <<"16N: "<<fit_amps["16N"]<<std::endl
					 <<"11Be: "<<fit_amps["11Be"]<<std::endl
					 <<"const_1: "<<fit_amps["const_1"]<<std::endl;
		}
		
		// draw the fit, if we want to check it
		/*
		dt_mu_lowe_hist.Draw();
		gPad->WaitPrimitive();
		*/
		
		// remove this intermediate fit from the list of functions
		// associated with our dt distribution, since we no longer need it
		dt_mu_lowe_hist.GetListOfFunctions()->Clear();
		
		break;
		}
	case 2:{
		// fit the range 0.1-0.8s with the components previously fit now fixed,
		// allowing additional components Li9 + a combination of 8He+9C
		Log(m_unique_name+"calling BuildFunction for time range case "+toString(rangenum)
			+", 9Li+8He_9C+8Li_8B+12B+12N+16N+11Be",v_debug,m_verbose);
		TF1 func_sum = BuildFunction({"9Li","8He_9C","8Li_8B","12B","12N","16N","11Be"},0.1,0.8);
		Log(m_unique_name+" retrieving results from past fits in case "+toString(rangenum)+" fit",v_debug,m_verbose);
		// pull fit results from the last two stages
		PullFitAmp(func_sum,"12B");
		PullFitAmp(func_sum,"12N");
		PullFitAmp(func_sum,"16N");
		PullFitAmp(func_sum,"11Be");
		// fit the new components
		Log(m_unique_name+" doing case "+toString(rangenum)+" fit",v_debug,m_verbose);
		std::string fitopt = (m_verbose>2) ? "Rq" : "R";
		dt_mu_lowe_hist.Fit(&func_sum,fitopt.c_str(),"",0.1,0.8);
		// record the results for the next step
		Log(m_unique_name+" recording results from case "+toString(rangenum)+" fit",v_debug,m_verbose);
		PushFitAmp(func_sum,"9Li");
		PushFitAmp(func_sum,"8He_9C");
		PushFitAmp(func_sum,"8Li_8B");
		fit_amps["const_2"] = func_sum.GetParameter("const");
		if(m_verbose){
			std::cout<<"recorded results were: "<<std::endl
					 <<"9Li: "<<fit_amps["9Li"]<<std::endl
					 <<"8He_9C: "<<fit_amps["8He_9C"]<<std::endl
					 <<"8Li_8B: "<<fit_amps["8Li_8B"]<<std::endl
					 <<"const_2: "<<fit_amps["const_2"]<<std::endl;
		}
		
		// draw the fit, if we want to check it
		/*
		dt_mu_lowe_hist.Draw();
		gPad->WaitPrimitive();
		*/
		
		// remove this intermediate fit from the list of functions
		// associated with our dt distribution, since we no longer need it
		dt_mu_lowe_hist.GetListOfFunctions()->Clear();
		
		/*
		// debug check
		std::cout<<"Drawing past fits and this one, to see how they look"<<std::endl;
		TF1 func_early = BuildFunction({"12B","12N"},50E-6,30);
		// pull fit results from the last two stages
		PullFitAmp(func_early,"12B");
		PullFitAmp(func_early,"12N");
		PullFitAmp(func_early,"const_0");
		func_early.SetLineColor(kBlue);
		std::cout<<"early"<<std::endl;
		func_early.Draw();
		//gPad->WaitPrimitive();
		TF1 func_late = BuildFunction({"16N","11Be"},50E-6,30);
		// pull fit results from the last two stages
		PullFitAmp(func_late,"16N");
		PullFitAmp(func_late,"11Be");
		PullFitAmp(func_late,"const_1");
		func_late.SetLineColor(kGreen-1);
		std::cout<<"late"<<std::endl;
		func_late.Draw();
		//gPad->WaitPrimitive();
		
		func_sum.SetRange(50E-6,30);
		std::cout<<"int"<<std::endl;
		func_sum.Draw();
		//gPad->WaitPrimitive();
		
		std::cout<<"everything"<<std::endl;
		dt_mu_lowe_hist.Draw();
		func_sum.Draw("same");
		func_early.Draw("same");
		func_late.Draw("same");
		//gPad->WaitPrimitive();
		*/
		
		break;
		}
	case 3:{
		// fit the range 0.8-6s with the components previously fit now fixed,
		// allowing additional components 15C + 16N
		Log(m_unique_name+"calling BuildFunction for time range case "+toString(rangenum)
			+", 15C+16N+8Li_8B",v_debug,m_verbose);  // FIXME update if we're going to keep everything
		/*
		TF1 func_sum = BuildFunction({"15C","16N","8Li_8B"},0.8,6);
		Log(m_unique_name+" retrieving results from past fits in case "+toString(rangenum)+" fit",v_debug,m_verbose);
		PullFitAmp(func_sum,"8Li_8B");
		PullFitAmp(func_sum,"16N",false); // FIXME add back in
		*/
		TF1 func_sum = BuildFunction({"12B","12N","16N","11Be","9Li","8He_9C","8Li_8B","15C"},0.8,6);
		Log(m_unique_name+" retrieving results from past fits in case "+toString(rangenum)+" fit",v_debug,m_verbose);
		PullFitAmp(func_sum,"12B");
		PullFitAmp(func_sum,"12N");
		PullFitAmp(func_sum,"11Be");
		PullFitAmp(func_sum,"9Li");
		PullFitAmp(func_sum,"8He_9C");
		PullFitAmp(func_sum,"8Li_8B");
		PullFitAmp(func_sum,"16N",false);
		
		// fit the new components
		Log(m_unique_name+" doing case "+toString(rangenum)+" fit",v_debug,m_verbose);
		std::string fitopt = (m_verbose>2) ? "Rq" : "R";
		dt_mu_lowe_hist.Fit(&func_sum,fitopt.c_str(),"",0.8,6);
		// record the results for the next step
		Log(m_unique_name+" recording results from case "+toString(rangenum)+" fit",v_debug,m_verbose);
		PushFitAmp(func_sum,"15C");
		PushFitAmp(func_sum,"16N");
		fit_amps["const_3"] = func_sum.GetParameter("const");
		
		if(m_verbose){
			std::cout<<"recorded results were: "<<std::endl
					 <<"15C: "<<fit_amps["15C"]<<std::endl
					 <<"16N: "<<fit_amps["16N"]<<std::endl
					 <<"const_3: "<<fit_amps["const_3"]<<std::endl;
		}
		
		// draw the fit, if we want to check it
		/*
		dt_mu_lowe_hist.Draw();
		gPad->WaitPrimitive();
		*/
		
		// remove this intermediate fit from the list of functions
		// associated with our dt distribution, since we no longer need it
		dt_mu_lowe_hist.GetListOfFunctions()->Clear();
		
		break;
		}
	case 4:{
		// final case: release all the parameters but keep the previously fit values as starting points.
		Log(m_unique_name+"calling BuildFunction for time range case "+toString(rangenum)
			+", everything",v_debug,m_verbose);
		TF1 func_sum = BuildFunction({"12B","12N","16N","11Be","9Li","8He_9C","8Li_8B","15C"},0,30);
		func_sum.SetLineColor(kRed);
		Log(m_unique_name+" retrieving results from past fits in case "+toString(rangenum)+" fit",v_debug,m_verbose);
		PullFitAmp(func_sum,"12B",false);
		PullFitAmp(func_sum,"12N",false);
		PullFitAmp(func_sum,"16N",false);
		PullFitAmp(func_sum,"11Be",false);
		PullFitAmp(func_sum,"9Li",false);
		PullFitAmp(func_sum,"8He_9C",false);
		PullFitAmp(func_sum,"8Li_8B",false);
		PullFitAmp(func_sum,"15C",false);
		
		Log(m_unique_name+" doing case "+toString(rangenum)+" fit",v_debug,m_verbose);
		std::string fitopt = (m_verbose>2) ? "RSq" : "RS";
		TFitResultPtr fitresult =  dt_mu_lowe_hist.Fit(&func_sum,fitopt.c_str(),"",0,30);
		
		Log(m_unique_name+" recording results from case "+toString(rangenum)+" fit",v_debug,m_verbose);
		if(m_verbose) std::cout<<"recorded results were: "<<std::endl;
		// for some reason TF1::GetParameter() was returning double values truncated to integer
		// Maybe this was something that happens when we fix the sign in PushFitAmp
		// Anyway, using the TFitResultPtr as that seems to work...
		for(auto&& anisotope : fit_amps){
			if(anisotope.first.substr(0,5)=="const") continue; // not a real isotope
			int par_number = func_sum.GetParNumber(("amp_"+anisotope.first).c_str());
			double amp = fitresult->Parameter(par_number);
			PushFitAmp(abs(amp),anisotope.first);
			if(m_verbose) std::cout<<anisotope.first<<": "<<amp<<std::endl;
		}
		fit_amps["const_4"] = func_sum.GetParameter("const");
		
		// Debug check
		Log(m_unique_name+" correcting fit parameter signs",v_debug,m_verbose);
		// we want to constrain the number of each isotope to be >0, which would
		// mean putting a lower limit of 0 on the amplitude in the fit.
		// unfortunately ROOT only lets us put both upper and lower limits (not just one)
		// and Minuit doesn't like suitable limits of e.g. 0 and 10E10 (all sorts of errors).
		// So instead we put no constraint on the parameter, but our function (from BuildFunction)
		// only uses its magnitude. This means we may end up with fit values that are negative,
		// but resulting function should be the same even if we flip the sign to positive.
		// As a debug check, let's correct the signs of all abundances to positive,
		// and then re-draw then function. It should still fit our data nicely.
		for(int pari=0; pari<func_sum.GetNpar(); ++pari){
			std::string parname = func_sum.GetParName(pari);
			if(parname.substr(0,9)=="lifetime_") continue; // skip lifetimes, they're good
			
//			// when adjusting fits interactively they're given a limit,
//			// but it seems like they don't have one in code, so skip this
//			double parmin, parmax;
//			func_sum.GetParLimits(pari,parmin,parmax);
//			double newupper = std::max(abs(parmin),abs(parmax));
//			double parval = func_sum.GetParameter(pari);
//			func_sum.SetParLimits(pari,parmin,newupper);
//			func_sum.SetParameter(pari,abs(parval));
//			func_sum.SetParLimits(pari,0,newupper);
			
			// much simpler in that case
			func_sum.SetParameter(pari,abs(func_sum.GetParameter(pari)));
		}
		dt_mu_lowe_hist.Draw();
		gPad->SetLogx();
		gPad->SetLogy();
//		gPad->Modified();
//		gPad->Update();
//		gPad->WaitPrimitive();
		
		// Debug check
		// we're having issues, so let's see how well the 2015 results fit our data
		TF1 func_paper = BuildFunctionNoHack({"12B","12N","16N","11Be","9Li","8He_9C","8Li","8B","15C"},0,30);
		func_paper.SetName("2015 Paper");
		func_paper.SetTitle("2015 Paper");
		// the 2015 paper had a lower energy threshold, so adjust results
		// based on the different efficiency of selection, to ease comparison
		bool correct_energy_threshold = true;
		Log(m_unique_name+" getting paper parameters for comparsion function",v_debug,m_verbose);
		PullPaperAmp(func_paper,"12B",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"12N",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"16N",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"11Be",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"9Li",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"8He_9C",correct_energy_threshold,paper_scaling);
		//PullPaperAmp(func_paper,"8Li_8B",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"8Li",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"8B",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"15C",correct_energy_threshold,paper_scaling);
		PullPaperAmp(func_paper,"const",correct_energy_threshold,paper_scaling);
		func_paper.SetLineColor(kViolet);
		func_paper.Draw("same");
		
		// Draw the individual isotopic contributions along with the total to reproduce 2015 Figure 3
		std::vector<TF1> indiv_funcs;
		indiv_funcs.reserve(fit_amps.size());
		for(auto&& theisotope : fit_amps){
			if(theisotope.first.substr(0,5)=="const") continue; // not a real isotope
			std::string anisotope = theisotope.first;
			TF1 next_func = BuildFunction({anisotope},0,30);
			PullFitAmp(next_func,anisotope);
			next_func.SetLineColor(colourwheel.GetNextColour());
			
			// we want to fix the y range of the plot to make all the contributions visible,
			// but setrangeuser doesn't seem to be working, so try constraining the x ranges
			// to the time range over which this function will be within our desired y range
			double ltime;
			if(anisotope.find("_")==std::string::npos){
				// not a pair
				ltime = lifetimes.at(anisotope);
			} else {
				// pair of isotopes
				// it's not possible to generically solve the necessary equation
				// to know what time range corresponds to a given y range,
				// but it should be good enough just take the average lifetime
				std::string first_isotope = anisotope.substr(0,anisotope.find_first_of("_"));
				std::string second_isotope = anisotope.substr(anisotope.find_first_of("_")+1,std::string::npos);
				ltime = (lifetimes.at(first_isotope) + lifetimes.at(second_isotope)) / 2.;
			}
			double miny = 1E-2;  // must be less than the plot min or it gets cut off early
			double maxx = -ltime * log(miny*(ltime/theisotope.second));
			if(maxx<0) maxx = 30;
			next_func.SetRange(0.0,maxx);
			
			// the following adds all functions to the plot, but it doesn't add them to the legend!!
			//dt_mu_lowe_hist.GetListOfFunctions()->Add(&indiv_funcs.back());
			
			// instead keep them in a vector
			indiv_funcs.push_back(next_func);
			indiv_funcs.back().SetName(anisotope.c_str());
			indiv_funcs.back().SetTitle(anisotope.c_str());
		}
		// and last but not least the constant background
		TF1 constfunc("constfunc","[0]",0,30);
		constfunc.SetParameter(0,fit_amps["const_4"]);
		indiv_funcs.push_back(constfunc);
		indiv_funcs.back().SetName("const");
		indiv_funcs.back().SetTitle("const");
		
		// draw it all
		dt_mu_lowe_hist.Draw();  // total fit is drawn implicitly as it is owned by the data histo
		dt_mu_lowe_hist.GetYaxis()->SetRangeUser(1.,1E7);
		for(auto&& afunc : indiv_funcs){ std::cout<<"drawing "<<afunc.GetName()<<std::endl; afunc.Draw("same"); }
		func_paper.Draw("same");
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		gPad->GetCanvas()->BuildLegend();
		gPad->Modified();
		gPad->Update();
		gPad->WaitPrimitive();
		
		// reproduce the paper plot including breakdowns, to check our value extraction
		//BuildPaperPlot();
		
		gPad->SetLogx(false);
		gPad->SetLogy(false);
		
		// one last thing: if we've been using BuildFunctionHack, then
		// we've been fitting only the *additional* number of events on top of
		// half the value found by the 2015 paper. Before returning we need to set
		// fit_amps to be the full number of events, including this contribution
		if(useHack){
			for(auto&& anisotope : fit_amps){
				if(anisotope.first.substr(0,5)=="const") continue; // not a real isotope
				// ensure the following matches what's being added by BuildFunctionHack
				double paperval = GetPaperAmp(anisotope.first, true, paper_scaling);
				anisotope.second += (paperval/2.);
			}
		}
		
		break;
		}
	default:
		Log(m_unique_name+" FitSubRangeDt invoked with invalid range "+toString(rangenum),v_error,m_verbose);
	}
	
	return true;
}

// =========================================================================
// =========================================================================

// Helper functions for constructing the TF1s with which we fit our dt distribution

// wrapper around either BuildFunctionNoHack or BuildFunctionHack, depending on whether we want
// to try to coerce the fit result to a number we like better
TF1 FitSpallationDt::BuildFunction(std::vector<std::string> isotopes, double func_min, double func_max){
	if(useHack){  // set via config file parameter
		return BuildFunctionHack(isotopes, func_min, func_max);
	} else {
		return BuildFunctionNoHack(isotopes, func_min, func_max);
	}
	return TF1{}; // to silence compiler warnings
}

TF1 FitSpallationDt::BuildFunctionNoHack(std::vector<std::string> isotopes, double func_min, double func_max){
	/* Construct a TF1 based on the list of isotopes given, over the time range given.
	   We name the function parameters and fix the lifetimes, since they're known. */
	
	std::string total_func="";
	std::string func_name="";
	std::map<std::string,int> parameter_posns;
	int next_par_index=0;
	for(std::string& anisotope : isotopes){
		func_name += anisotope+"_";
		//std::cout<<"adding isotope "<<anisotope<<std::endl;
		// first, this 'isotope' may be a degenerate pair, so try to split it apart
		if((!split_iso_pairs) || anisotope.find("_")==std::string::npos){
			// not a pair
			//std::cout<<"not a pair"<<std::endl;
			int first_index=next_par_index;
			// "([0]/[1])*exp(-x/[1])"
			
			// Or that's the function you'd expect. Fig 3 plots from 0--30s, with x-axis in seconds;
			// which means that F(t) = dN/dt is also in seconds ... but Fig 3's y-axis is in events / 0.006s!!
			// The bin width is variable (it's a log-log plot), so we already need to scale our bin counts
			// by the bin width to get consistent units, so there's no reason not to use events/second.
			// Still, to make a comparable plot, we could scale our histogram bin counts up using TH1::Scale,
			// but then our fit values will be off unless our fit function accounts for it.
			// (this also ensures the paper plot overlay comparison has the correct scaling).
			
			std::string scalestring = std::to_string(binwidth)+"*";  // "0.006*"
			std::string this_func =  scalestring+"(abs(["+toString(next_par_index)+"])"  // original
//			std::string this_func =  scalestring+"(["+toString(next_par_index)+"]"         // laura
									+"/["+toString(next_par_index+1)
									+"])*exp(-x/["+toString(next_par_index+1)+"])";
			next_par_index +=2;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index);
			parameter_posns.emplace("lifetime_"+anisotope,first_index+1);
		} else {
			// it's a pair. For now, only support two isotopes at a time.
			//std::cout<<"a pair, splitting into ";
			int first_index=next_par_index;
			std::string first_isotope = anisotope.substr(0,anisotope.find_first_of("_"));
			std::string second_isotope = anisotope.substr(anisotope.find_first_of("_")+1,std::string::npos);
			//std::cout<<first_isotope<<" and "<<second_isotope<<std::endl;
			// the fit function isn't just the sum of two single isotope functions
			// as they share an amplitude and constant
			//"[0]*0.5*(exp(-x/[1])/[1]+exp(-x/[2])/[2])
			// as above, add in the additional scaling factor to get Y units of events/0.006s
			std::string scalestring = std::to_string(binwidth)+"*";  // "0.006*"
			std::string this_func =  scalestring+"abs(["+toString(next_par_index)+"])*0.5*"
									+"(exp(-x/["+toString(next_par_index+1)+"])/" // don't increment index
									+"["+toString(next_par_index+1)+"]+"
									+"exp(-x/["+toString(next_par_index+2)+"])/"  // don't increment index
									+"["+toString(next_par_index+2)+"])";
			next_par_index += 3;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index++);
			parameter_posns.emplace("lifetime_"+first_isotope,first_index++);
			parameter_posns.emplace("lifetime_"+second_isotope,first_index);
		}
	}
	// add the constant term
	total_func += " + abs(["+toString(next_par_index)+"])";
	parameter_posns.emplace("const",next_par_index);
	
	// build the total function from the sum of all strings
	func_name.pop_back(); // remove trailing '_'
	//std::cout<<"function is "<<total_func<<std::endl;
	TF1 afunc(func_name.c_str(),total_func.c_str(),func_min,func_max);
	
	// OK, propagate parameter names to the function
	for(auto&& next_par : parameter_posns){
		//std::cout<<"setting parname "<<next_par.second<<" to "<<next_par.first<<std::endl;
		afunc.SetParName(next_par.second,next_par.first.c_str());
		if(next_par.first.substr(0,9)=="lifetime_"){
			std::string anisotope = next_par.first.substr(9,std::string::npos);
			//std::cout<<"fixing lifetime of "<<anisotope<<std::endl;
			FixLifetime(afunc,anisotope);
		} else {
			// try setting par limits
			if(use_par_limits) afunc.SetParLimits(next_par.second,0,10000000);
			// for amplitudes / background constant we'll set the lower limit to be 0.
			// This seems to be necessary for some fits otherwise
			// ROOT gives negative amounts of some isotopes.
			//afunc.SetParLimits(next_par.second, 0, 1E9); // FIXME what to use as upper limit????
			// strictly i think we ought to give initial fit values,
			// but thankfully MINUIT manages without them, but it doesn't like having
			// (default) initial values at the limit, so set some small amount
			//afunc.SetParameter(next_par.second,10);
		}
	}
	
	// to try to reproduce Laura's reproduction of the 2015 paper.... do not use the constant.
	if(fix_const){
		int const_par_number = afunc.GetParNumber("const");
		afunc.FixParameter(const_par_number,0);
	}
	
	// return the built function
	return afunc;
}

// this is the hack
TF1 FitSpallationDt::BuildFunctionHack(std::vector<std::string> isotopes, double func_min, double func_max){
	std::string total_func="";
	std::string func_name="";
	std::map<std::string,int> parameter_posns;
	int next_par_index=0;
	for(std::string& anisotope : isotopes){
		func_name += anisotope+"_";
		//std::cout<<"adding isotope "<<anisotope<<std::endl;
		// first, this 'isotope' may be a degenerate pair, so try to split it apart
		if((!split_iso_pairs) || anisotope.find("_")==std::string::npos){
			// not a pair
			//std::cout<<"not a pair"<<std::endl;
			int first_index=next_par_index;
			
			// get paper amplitude, corrected for energy efficiency and scaled by config file scaling
			double paperval = GetPaperAmp(anisotope,true,paper_scaling);
			//std::cout<<"paperval= "<<paperval<<std::endl;
			std::string paperstring = std::to_string(paperval/2.); // fix half the paper val, fit the rest
			
			// "((x.xx + abs([0]))/[1])*exp(-x/[1])"
			std::string scalestring = std::to_string(binwidth)+"*";  // "0.006*"
			std::string this_func =  scalestring+"(("+paperstring+"+abs(["+toString(next_par_index)
									+"]))/["+toString(next_par_index+1)
									+"])*exp(-x/["+toString(next_par_index+1)+"])";
			next_par_index +=2;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index);
			parameter_posns.emplace("lifetime_"+anisotope,first_index+1);
		} else {
			// it's a pair. For now, only support two isotopes at a time.
			//std::cout<<"a pair, splitting into ";
			std::string first_isotope = anisotope.substr(0,anisotope.find_first_of("_"));
			std::string second_isotope = anisotope.substr(anisotope.find_first_of("_")+1,std::string::npos);
			//std::cout<<first_isotope<<" and "<<second_isotope<<std::endl;
			int first_index=next_par_index;
			
			// get paper amplitude, corrected for energy efficiency and scaled by config file scaling
			//std::cout<<"getting paper amp"<<std::endl;
			double paperval = GetPaperAmp(anisotope,true,paper_scaling);
			//std::cout<<"paperval= "<<paperval<<std::endl;
			// fix the abundance to at least half the paper val, fit the rest
			std::string paperstring = std::to_string(paperval/2.);
			
			// the fit function isn't just the sum of two single isotope functions
			// as they share an amplitude and constant
			//"(x.xx + abs([0]))*0.5*(exp(-x/[1])/[1]+exp(-x/[2])/[2])
			std::string scalestring = std::to_string(binwidth)+"*";  // "0.006*"
			std::string this_func =  scalestring+"("+paperstring+"+abs(["+toString(next_par_index)+"]))*0.5*"
									+"(exp(-x/["+toString(next_par_index+1)+"])/" // don't increment index
									+"["+toString(next_par_index+1)+"]+"
									+"exp(-x/["+toString(next_par_index+2)+"])/"  // don't increment index
									+"["+toString(next_par_index+2)+"])";
			next_par_index += 3;
			// add this function to the total function string
			total_func = (total_func=="") ? this_func : total_func+" + "+this_func;
			// add the parameter names to our map
			parameter_posns.emplace("amp_"+anisotope,first_index++);
			parameter_posns.emplace("lifetime_"+first_isotope,first_index++);
			parameter_posns.emplace("lifetime_"+second_isotope,first_index);
		}
	}
	// add the constant term
	total_func += " + abs(["+toString(next_par_index)+"])";
	parameter_posns.emplace("const",next_par_index);
	
	// build the total function from the sum of all strings
	func_name.pop_back(); // remove trailing '_'
	std::cout<<"function is "<<total_func<<std::endl;
	TF1 afunc(func_name.c_str(),total_func.c_str(),func_min,func_max);
	
	// OK, propagate parameter names to the function
	for(auto&& next_par : parameter_posns){
		//std::cout<<"settin parname "<<next_par.second<<" to "<<next_par.first<<std::endl;
		afunc.SetParName(next_par.second,next_par.first.c_str());
		if(next_par.first.substr(0,9)=="lifetime_"){
			std::string anisotope = next_par.first.substr(9,std::string::npos);
			//std::cout<<"fixing lifetime of "<<anisotope<<std::endl;
			FixLifetime(afunc,anisotope);
		} else {
			// for amplitudes / background constant we'll set the lower limit to be 0.
			// This seems to be necessary for some fits otherwise
			if(use_par_limits) afunc.SetParLimits(next_par.second,0,10000000);
			// ROOT gives negative amounts of some isotopes.
			//afunc.SetParLimits(next_par.second, 0, 1E9); // FIXME what to use as upper limit????
			// strictly i think we ought to give initial fit values,
			// but thankfully MINUIT manages without them, but it doesn't like having
			// (default) initial values at the limit, so set some small amount
			//afunc.SetParameter(next_par.second,10);
		}
	}
	
	// to try to reproduce Laura's reproduction of the 2015 paper.... do not use the constant.
	if(fix_const){
		int const_par_number = afunc.GetParNumber("const");
		afunc.FixParameter(const_par_number,0);
	}
	
	// return the built function
	return afunc;
}

// fix lifetime parameter of a TF1
void FitSpallationDt::FixLifetime(TF1& func, std::string isotope){
	int par_number = func.GetParNumber(("lifetime_"+isotope).c_str());
	func.FixParameter(par_number,lifetimes.at(isotope));
}

// extract an abundance paramerer result from a TF1 and record it into our map of amplitudes
void FitSpallationDt::PushFitAmp(TF1& func, std::string isotope){
//	std::cout<<"fixing fit result for func "<<func.GetName()<<", isotope "<<isotope<<std::endl;
//	std::cout<<"available pars are: "<<std::endl;
//	for(int i=0; i<func.GetNpar(); ++i){
//		std::cout<<i<<"="<<func.GetParName(i)<<std::endl;
//	}
//	int par_number = func.GetParNumber(("amp_"+isotope).c_str()); // can do straight by name
	double v = static_cast<double>(func.GetParameter(("amp_"+isotope).c_str()));
	fit_amps[isotope] = func.GetParameter(("amp_"+isotope).c_str());
}

// record an abundance paramerer result, given as a double, into our map of amplitudes
void FitSpallationDt::PushFitAmp(double amp, std::string isotope){
	fit_amps[isotope] = amp;
}


// use our map of amplitudes to set the value of a TF1 abundance parameter
void FitSpallationDt::PullFitAmp(TF1& func, std::string isotope, bool fix){
	std::string parname = (isotope.substr(0,5)=="const") ? "const" : "amp_"+isotope;
	int par_number = func.GetParNumber(parname.c_str());
	if(fix){
		func.FixParameter(par_number,fit_amps.at(isotope));
	} else {
		func.SetParameter(par_number,fit_amps.at(isotope));
	}
}

// =========================================================================
// =========================================================================

// Various functions for retrieving results based on the 2015 paper, for comparison

double FitSpallationDt::GetPaperAmp(std::string isotope, bool threshold_scaling, double fixed_scaling){
	double paperval;
	if(isotope!="const"){
		// the values in papervals are rates in events / kton / day, integrated over all beta energies.
		// To obtain Ni, the initial number of observable events over our live time,
		// we need to multiply by the livetime, the fiducial volume, and the efficiency of observation
		// (i.e. total efficiency all cuts)
		// Everything except the low energy cut is isotope independent.
		// the efficiency of the energy cut we need to look up for each isotope
		// The threshold_scaling bool defines whether we use the 6MeV efficiency from the paper,
		// or an 8MeV efficiency from FLUKA + skdetsim. For the same livetime this would give the number
		// of events we would expect to see with a higher E threshold so that we can compare.
		double energy_cut_eff = 1.;
		if(isotope.find("_")==std::string::npos){
			// not a pairing - can look up efficiency directly
			if(threshold_scaling){
				energy_cut_eff = reco_effs_8mev.at(isotope);
			} else {
				energy_cut_eff = papereffs.at(isotope);
			}
			paperval = papervals.at(isotope);
		} else {
			std::string first_isotope = isotope.substr(0,isotope.find_first_of("_"));
			std::string second_isotope = isotope.substr(isotope.find_first_of("_")+1,std::string::npos);
			// how do we handle efficiency for pairs?
			// best we can do is assume half for each and average the efficiency i think....
			if(threshold_scaling){
				energy_cut_eff = 0.5*(reco_effs_8mev.at(first_isotope) + reco_effs_8mev.at(second_isotope));
			} else {
				energy_cut_eff = 0.5*(papereffs.at(first_isotope) + papereffs.at(second_isotope));
			}
			// for a pair of isotopes we may have one or both
			if(papervals.count(isotope)){
				// if we have a combined amplitude, use that
				paperval = papervals.at(isotope);
			} else {
				// otherwise take the sum
				double paperval1 = papervals.at(first_isotope);
				double paperval2 = papervals.at(second_isotope);
				paperval = (paperval1+paperval2);
			}
		}
		// ok, convert Ri to Ni
//		std::cout<<"calculating paper amplitude for "<<isotope<<", rate = "
//				 <<paperval<<", FV="<<fiducial_vol<<", T="<<paper_livetime
//				 <<" energy cut eff = "<<energy_cut_eff<<"%, total eff = "
//				 <<((paper_first_reduction_eff/100.)*(paper_dlt_cut_eff/100.)*(energy_cut_eff/100.)*100.)
//				 <<std::endl;
		paperval *= fiducial_vol * paper_livetime * (paper_first_reduction_eff/100.)
					 * (paper_dlt_cut_eff/100.) * (energy_cut_eff/100.);
//		std::cout<<"initial abundance = "<<paperval<<std::endl;
	} else {
		// the constant term isn't a rate, it's read straight off the plot so doesn't need conversion.
		// But, if we're comparing across energy thresholds, it should also be scaled
		// down by the relative rate of random backgrounds above 8MeV vs above 6MeV
		// TODO obtain that number...
		paperval = papervals.at(isotope); // for now just use the value read off the plot, no scaling
	}
	
	paperval *= fixed_scaling; // superfluous, just in case
	return paperval;
}

void FitSpallationDt::PullPaperAmp(TF1& func, std::string isotope, bool threshold_scaling, double fixed_scaling){
	double paperval = GetPaperAmp(isotope, threshold_scaling, fixed_scaling);
	std::string parname = (isotope=="const") ? "const" : "amp_"+isotope;
	func.SetParameter(parname.c_str(), paperval);
}

// Debug function
void FitSpallationDt::BuildPaperPlot(){
	/* double check we can reproduce Figure 3 from the paper */
	std::cout<<"BuildPaperPlot reproducing Fig 3"<<std::endl;
	
	// add the individual isotopic contributions to reproduce the old plot
	std::vector<TF1> indiv_funcs;
	std::vector<std::string> isotope_names;
	indiv_funcs.reserve(papervals.size());
	for(auto&& theisotope : papervals){
		if(theisotope.first.substr(0,5)=="const") continue; // not a real isotope
		if(theisotope.second==0) continue; // do not add to plot isotopes with no abundance
		
		std::string anisotope = theisotope.first;
		TF1 next_func = BuildFunctionNoHack({anisotope},0,30);
		
		PullPaperAmp(next_func,anisotope,false);
		next_func.SetLineColor(colourwheel.GetNextColour());
		
		indiv_funcs.push_back(next_func);
		indiv_funcs.back().SetName(anisotope.c_str());
		indiv_funcs.back().SetTitle(anisotope.c_str());
		isotope_names.push_back(anisotope);
	}
	// and last but not least the constant background
	TF1 constfunc("constfunc","[0]",0,30);
	constfunc.SetParameter(0,papervals.at("const"));
	indiv_funcs.push_back(constfunc);
	indiv_funcs.back().SetName("const");
	indiv_funcs.back().SetTitle("const");
	
	// now the sum of everything
	TF1 func_paper = BuildFunctionNoHack(isotope_names,0.001,30);
	func_paper.SetName("total");
	func_paper.SetTitle("total");
	// set the amplitudes
	for(auto&& theisotope : isotope_names){
		PullPaperAmp(func_paper,theisotope,false);
	}
	func_paper.SetParameter("const",papervals.at("const"));
	
	func_paper.SetLineColor(kViolet);
	std::cout<<"Drawing paper total"<<std::endl;
	func_paper.Draw();
	func_paper.GetYaxis()->SetRangeUser(1.,3E7);
//	gPad->Modified();
//	gPad->Update();
//	gPad->WaitPrimitive();
	
	// add all the components
	for(auto&& afunc : indiv_funcs){
		std::cout<<"drawing "<<afunc.GetName()<<std::endl;
		afunc.Draw("same");
//		gPad->Modified();
//		gPad->Update();
//		gPad->WaitPrimitive();
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gPad->GetCanvas()->BuildLegend();
	std::cout<<"Done, drawing paper Fig 3"<<std::endl;
	gPad->Modified();
	gPad->Update();
	gPad->WaitPrimitive();
	
}

// =========================================================================
// =========================================================================

// A function to retrieve the selection efficiency of a 6 or 8 MeV reconstructed energy cut
// for each isotope, based on FLUKA + skdetsim simulations

bool FitSpallationDt::GetEnergyCutEfficiencies(){
	// 2015 paper had a low threshold of 6MeV, newer data has (as of now) 8MeV
	// so to compare yeilds we need to scale the paper values by the fraction
	// of decays that would be above the respective thresholds
	
	std::map<std::string,int> true_events_below_8MeV;
	std::map<std::string,int> true_events_above_8MeV;
	std::map<std::string,int> true_events_below_6MeV;
	std::map<std::string,int> true_events_above_6MeV;
	
	// same but using energy from bonsai
	std::map<std::string,int> reco_events_below_8MeV;
	std::map<std::string,int> reco_events_above_8MeV;
	std::map<std::string,int> reco_events_below_6MeV;
	std::map<std::string,int> reco_events_above_6MeV;
	
	BStore efficiencyStore(true);
	efficiencyStore.Initnew(efficienciesFile, uncompressed, true);
	efficiencyStore.Get("true_events_below_8MeV",true_events_below_8MeV);
	efficiencyStore.Get("true_events_above_8MeV",true_events_above_8MeV);
	efficiencyStore.Get("true_events_below_6MeV",true_events_below_6MeV);
	efficiencyStore.Get("true_events_above_6MeV",true_events_above_6MeV);
	
	efficiencyStore.Get("reco_events_below_8MeV",reco_events_below_8MeV);
	efficiencyStore.Get("reco_events_above_8MeV",reco_events_above_8MeV);
	efficiencyStore.Get("reco_events_below_6MeV",reco_events_below_6MeV);
	efficiencyStore.Get("reco_events_above_6MeV",reco_events_above_6MeV);
	
	for(auto&& anisotope : true_events_below_6MeV){
		std::string isotope = anisotope.first;
		double true_below6 = true_events_below_6MeV.at(isotope);
		double true_above6 = true_events_above_6MeV.at(isotope);
		double true_eff6 = true_above6/(true_below6+true_above6);
		double true_below8 = true_events_below_8MeV.at(isotope);
		double true_above8 = true_events_above_8MeV.at(isotope);
		double true_eff8 = true_above8/(true_below8+true_above8);
		std::cout<<"TRUE ENERGY:"<<std::endl;
		std::cout<<"Isotope "<<anisotope.first<<" had "<<true_below6<<" events below 6MeV vs "
				 <<true_above6<<" above 6 MeV corresponding to an efficiency of "<<(true_eff6*100.)<<"%"
				 <<std::endl<<"compared to "<<true_below8<<" events below 8MeV and "<<true_above8
				 <<" events above 8 MeV corresponding to an efficiency of "<<(true_eff8*100.)<<"%"<<std::endl;
		
		true_effs_6mev.emplace(isotope,true_eff6);
		true_effs_8mev.emplace(isotope,true_eff8);
		true_effs_scaling.emplace(isotope,(true_eff8/true_eff6));
		
		// same with reconstructed energies - the spectra are VERY distorted...!?
		double reco_below6 = reco_events_below_6MeV.at(isotope);
		double reco_above6 = reco_events_above_6MeV.at(isotope);
		double reco_eff6 = reco_above6/(reco_below6+reco_above6);
		double reco_below8 = reco_events_below_8MeV.at(isotope);
		double reco_above8 = reco_events_above_8MeV.at(isotope);
		double reco_eff8 = reco_above8/(reco_below8+reco_above8);
		std::cout<<"RECONSTRUCTED ENERGY:"<<std::endl;
		std::cout<<"Isotope "<<anisotope.first<<" had "<<reco_below6<<" events below 6MeV vs "
				 <<reco_above6<<" above 6 MeV corresponding to an efficiency of "<<(reco_eff6*100.)<<"%"
				 <<std::endl<<"compared to "<<reco_below8<<" events below 8MeV and "<<reco_above8
				 <<" events above 8 MeV corresponding to an efficiency of "<<(reco_eff8*100.)<<"%"<<std::endl;
		
		reco_effs_6mev.emplace(isotope,reco_eff6*100.);
		reco_effs_8mev.emplace(isotope,reco_eff8*100.);
		reco_effs_scaling.emplace(isotope,(reco_eff8/reco_eff6));
	}
	return true;
}
