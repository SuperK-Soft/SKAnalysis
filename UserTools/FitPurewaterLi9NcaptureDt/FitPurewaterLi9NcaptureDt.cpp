/* vim:set noexpandtab tabstop=4 wrap */
#include "FitPurewaterLi9NcaptureDt.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"
#include "MTreeReader.h"
#include "MTreeSelection.h"

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TString.h"

// For Fitting
#include "Fit/Fitter.h"
#include "Fit/UnBinData.h"
#include "Fit/FitConfig.h"
#include "Math/WrappedMultiTF1.h"

FitPurewaterLi9NcaptureDt::FitPurewaterLi9NcaptureDt():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

// from 2015 paper Table I
constexpr double li9_lifetime_secs = 0.26;
constexpr double ncapture_lifetime_secs = 204.8E-6;

bool FitPurewaterLi9NcaptureDt::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("outputFile",outputFile);          // where to save data. If empty, current TFile
	m_variables.Get("li9_ncapture_dtmin",ncap_dtmin);
	m_variables.Get("li9_ncapture_dtmax",ncap_dtmax);
	m_variables.Get("treeReaderName",treeReaderName);
	
	myTreeReader = m_data->Trees.at(treeReaderName);
	myTreeSelections = m_data->Selectors.at(treeReaderName);
	
	return true;
}

bool FitPurewaterLi9NcaptureDt::Execute(){
	
	// retrieve variables from TreeReader
	GetBranchValues();
	
	// the following cuts are based on muon-lowe pair variables, so loop over muon-lowe pairs
	std::set<size_t> spall_mu_indices = myTreeSelections->GetPassingIndexes("dlt_mu_lowe>200cm");
	Log(toolName+" Looping over "+toString(spall_mu_indices.size())
				+" preceding muons to look for spallation events",v_debug,verbosity);
	for(size_t mu_i : spall_mu_indices){
		// now check whether this passed the additional Li9 cuts
		if(not myTreeSelections->GetPassesCut("ntag_FOM>0.995")) continue;
		// plot distribution of beta->ntag dt from passing triplets, compare to fig 5
		// Zhang had no events with >1 ntag candidate: should we only take the first? XXX
		for(size_t neutron_i=0; neutron_i<num_neutron_candidates; ++neutron_i){
			if(myTreeSelections->GetPassesCut("mu_lowe_ntag_triplets",{mu_i, neutron_i})){
				double ncap_time = dt_lowe_n[neutron_i];  // these times are in nanoseconds
				// according to sonias code we need to account for some offset of the AFT trigger timestamps?
				// "Remove 10mus time shift for AFT events and convert time to microseconds"
				// doesn't seem to tie up with what this is actually doing, though
				// adjusted too instead convert presumably ms, to seconds for consistency
				double ncap_time_adjusted = ncap_time < 50000 ? ncap_time : ncap_time - 65000;
				li9_ntag_dt_vals.push_back(ncap_time_adjusted/1E9);
			}
		}
	} // end loop over muons
	
	return true;
}

bool FitPurewaterLi9NcaptureDt::GetBranchValues(){
	// retrieve variables from branches
	bool success = 
	(myTreeReader->Get("np", num_neutron_candidates)) &&
	(myTreeReader->Get("dt", dt_lowe_n));
	
	return success;
}

bool FitPurewaterLi9NcaptureDt::Finalise(){
	
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
	
	Log(toolName+" Fitting Li9 Ntag candidate dt distribution",v_debug,verbosity);
	PlotNcaptureDt();
	
	if(fout!=nullptr){
		fout->Close();
		delete fout;
		fout=nullptr;
	}
	
	return true;
}

// =========================================================================
// Li9 ncapture fits
// =========================================================================

bool FitPurewaterLi9NcaptureDt::PlotNcaptureDt(){
	
	// make a histogram to bin the data
	std::cout<<"making li9 lowe-ncapture dt histogram"<<std::endl;
	TH1F li9_ncap_dt_hist("li9_ncap_dt_hist","Beta to ncapture dt for Li9 triplets",
	                       21,0,500E-6);
	
	// XXX debug investigation - no sign of expl decay of ncapture times... time span is correct,
	// there's no data beyond 500us. What's going on?
	std::cout<<"first 100 ncapture times were: {";
	int ncpi=0;
	for(auto&& aval : li9_ntag_dt_vals){
		// XXX FIXME REMOVE AFTER REPROCESSING IN ANALYSE XXX XXX XXX XXX XXX XXX 
		double ncap_time_adjusted = aval < 50000 ? aval : aval - 65000;
		ncap_time_adjusted /= 1E9;
		if(ncpi<100) std::cout<<ncap_time_adjusted<<", "; ++ncpi;
		li9_ncap_dt_hist.Fill(ncap_time_adjusted);  // FIXME weight by num_post_muons and num neutrons
	}
	std::cout<<"}"<<std::endl;
	std::cout<<"saving to file"<<std::endl;
	// for comparison to the paper, scale the x axis up to us
	li9_ncap_dt_hist.GetXaxis()->SetLimits(0,500);  // changes axis labels but doesn't affect binning
	li9_ncap_dt_hist.Write();
	// set it back for analysis
	li9_ncap_dt_hist.GetXaxis()->SetLimits(0,500E-6);  // XXX this is bad practice
	
	// TODO get expected number of background events
	// this is also needed for the expected background E spectrum, where we also need
	// the energies of the background, so calculate that first and we can just take the number
	
	// fit the ncapture lifetime, Fig 5 from the paper.
	// do we do chi2 fit or binned likelihood fit? why? FIXME how do we do binned likelihood?
	
	// perform an unbinned chi2 fit
	std::cout<<"doing binned chi2 fit"<<std::endl;
	double binned_estimate = BinnedNcapDtChi2Fit(&li9_ncap_dt_hist);
	
	// perform an unbinned likelihood fit
	std::cout<<"doing unbinned likelihood fit"<<std::endl;
	UnbinnedNcapDtLogLikeFit(&li9_ncap_dt_hist, binned_estimate);
	
	return true;
}

double FitPurewaterLi9NcaptureDt::BinnedNcapDtChi2Fit(TH1F* li9_ncap_dt_hist){
	// this fits the lifetime of ncapture to extract the amount of exponential and constant
	std::cout<<"making TF1 for binned chi2 fit with "<<li9_ncap_dt_hist->GetEntries()<<" values"<<std::endl;
	
	// number of neutrons left after time t = N = N0*exp(-t/τ)
	// rate of change of num neutrons, i.e. rate of observed ncapture events = dN/dt = -N0*τ*exp(-t/τ)
	// we can ignore the sign which just says N is decreasing.
	// simple chi2 fit of 'y = C + A*τ*exp(-dt/τ)' for all values of C, A, τ
	TF1 ncap_dt_func("ncap_dt_func","[0]+[1]*[2]*exp(-x/[2])",ncap_dtmin,ncap_dtmax);
	
	// set parameter names
	ncap_dt_func.SetParName(0,"rate of bg events");
	ncap_dt_func.SetParName(1,"num of Li9+n events");
	ncap_dt_func.SetParName(2,"ncapture lifetime");
	
	// set starting values - we can fix the ncapture lifetime
	//ncap_dt_func.SetParameters(allparams.data());  // pass an array or set individually
	ncap_dt_func.SetParameter(0,0);  // TODO estimate from accidental rate of ntag from MC and num events?
	ncap_dt_func.SetParameter(1,li9_ncap_dt_hist->GetBinContent(1));
	ncap_dt_func.FixParameter(2,ncapture_lifetime_secs);
	
	// set num TF1 points for drawing
	ncap_dt_func.SetNpx(1000);
	
	// DO THE FIT
	std::cout<<"invoking li9 ncapture binned chi2 fit"<<std::endl;
	TFitResultPtr fitresult = li9_ncap_dt_hist->Fit("ncap_dt_func","MRS");
	// options here are  M: better fit, R: use range, S: get resultsptr
	
	// print result
	std::cout<<"ncapture dt binned chi2 fit parameters = {";
	for(int i=0; i<(ncap_dt_func.GetNpar()); i++){
		std::cout<<ncap_dt_func.GetParameter(i);
		(i<(ncap_dt_func.GetNpar()-1)) ? std::cout<<", " : std::cout<<"};\n";
	}
	//float fitchi2 = fitresult->Chi2();                 // same as below, which
	float fitchi2 = ncap_dt_func.GetChisquare();         // doesn't need fitresultptr
	Log(toolName+" li9 lowe->ncap dt fit chi2 was "+toString(fitchi2),v_message,verbosity);
	
	// draw result
	li9_ncap_dt_hist->Draw();
	gPad->WaitPrimitive();
	gPad->Clear();
//	//fit->Draw("lsame");   // not necessary, fit is added to histogram's list of functions and drawn automatically
	li9_ncap_dt_hist->GetListOfFunctions()->Clear();
	
	std::cout<<"binned ncapture chi2 fit done"<<std::endl;
	return ncap_dt_func.GetParameter(1);
}

bool FitPurewaterLi9NcaptureDt::UnbinnedNcapDtLogLikeFit(TH1F* li9_ncap_dt_hist, double num_li9_events){
	
	std::cout<<"doing unbinned likelihood fit with "<<li9_ntag_dt_vals.size()<<" values"<<std::endl;
	
	// TF1 of likelihood distribution
	// ROOT provides some basic pdfs: https://root.cern.ch/doc/v610/group__PdfFunc.html
	// but we need our likelihood to be normalized and I don't see how to do it suitably with these
//	TF1 ncapture_dt_unbinned_like("ncapture_dt_unbinned_like",
//	TString::Format("[0]*ROOT::Math::uniform_pdf (x,%.2f,%.2f,0) "
//	                " + [1]*ROOT::Math::exponential_pdf(x,[1],0)",xmin,xmax),xmin,xmax);
	
	// so let's make our own
	int fit_n_pars = 1;
	int func_n_dims = 1;
	std::cout<<"making TF1"<<std::endl;
	TF1 ncapture_dt_unbinned_like("ncapture_dt_unbinned_like",this,&FitPurewaterLi9NcaptureDt::ncap_lifetime_loglike,
	        ncap_dtmin, ncap_dtmax, fit_n_pars, "FitPurewaterLi9NcaptureDt", "ncap_lifetime_loglike");
	
	// set parameter names
	ncapture_dt_unbinned_like.SetParName(0,"fraction of background events");
	ncapture_dt_unbinned_like.SetParLimits(0,0,1.);
	
	//ncapture_dt_unbinned_like.SetParameters(allparams.data());  // pass an array or set individually
	// set starting values to fit val from binned fit
	std::cout<<"setting starting value of background fraction to "
			 <<(1.-(num_li9_events/li9_ntag_dt_vals.size()))<<std::endl;
	ncapture_dt_unbinned_like.SetParameter(0,(1.-(num_li9_events/li9_ntag_dt_vals.size())));
	
	int nbins=10;
	for(int i=nbins; i>0; --i){
		std::vector<double> pars{ncapture_dt_unbinned_like.GetParameter(0)};
		double xval = ncap_dtmin+double(i)*((ncap_dtmax-ncap_dtmin)/double(nbins-1));
		std::cout<<"ncapture unbinned function eval at "<<xval<<" = "
				 <<ncapture_dt_unbinned_like.Eval(xval)<<std::endl
				 <<", base function = "<<ncap_lifetime_loglike(&xval,pars.data())
				 <<std::endl;
	}
	
	// draw for comparison. Since our likelihood function is normalised we need to normalise the data too
	TH1F* li9_ncap_dt_hist_normalised = (TH1F*)li9_ncap_dt_hist->Clone();
	li9_ncap_dt_hist_normalised->Scale(1./li9_ncap_dt_hist->Integral());
	li9_ncap_dt_hist_normalised->Draw();
	ncapture_dt_unbinned_like.Draw("same");
	gPad->WaitPrimitive();
	return true;
	
	// convert the TF1 for use in unbinned fit
	std::cout<<"making TF1 wrapper"<<std::endl;
	ROOT::Math::WrappedMultiTF1 ncapture_dt_unbinned_func(ncapture_dt_unbinned_like,
	                                                      ncapture_dt_unbinned_like.GetNdim());
	
	// build a fitter
	// 'false' says let the fitter calculate derivatives of the function itself (recommended)
	std::cout<<"making fitter"<<std::endl;
	ROOT::Fit::Fitter fitter;
	fitter.SetFunction(ncapture_dt_unbinned_func);  // no bool in ROOT 5 - what's the default?
	// important note: the fit function should be normalized, but simply calculating and dividing by
	// the integral before returning may not be enough as this introduces a correlation
	// between the fit parameters. Apparently the best way is to fit N-1 parameters, then derive
	// the value of the remaining one from the others and a fixed normalization condition
	// https://root-forum.cern.ch/t/root-fit-unbindata-example-is-not-working/28462/3
	// in cases where it is sufficient to simply divide by the integral (which??) 
	// a wrapper is given in https://root-forum.cern.ch/t/unbinned-log-likelihood-fits/18981/2
	
	// can we do it just by fixing the parameter?
	// https://web.pa.msu.edu/people/brock/file_sharing/ATLAS/root/math/mathcore/test/fit/testRooFit.cxx
	//fitter.Config().ParSettings(0).Fix();
	
	// configure the fitter
	// there is a fitter.Config().ParSettings(#) for each parameter #
	// these are available only after calling ROOT::Fit::Fitter::SetFunction
	// the configurable options are:
	// the initial values of the parameters;   (can also be set by input model func)
	// the parameter step sizes;
	// eventual parameter bounds;
	// the minimizer library and the particular algorithm to use;
	// different minimization options (print level, tolerance, max iterations, etc.)
	// the type of parameter errors to compute (parabolic error, Minos errors, re-normalized errors using fitted chi2 values)
	
	// set initial values. Think it also inherits initial values from TF1 on construction?
	//double initialParams[] = {2,1};
	//fitter.Config().SetParamsSettings(2,initialParams);
	
	// set the minimizer algorithm to use
	// https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
	// https://root.cern.ch/root/htmldoc/ROOT__Fit__FitConfig.html
	fitter.Config().SetMinimizer("Minuit2","Migrad");
//	fitter.Config().SetUpdateAfterFit();  // not in root 5
	
	// minimizer options include:
	// Minimizer type (Minuit, Fumili, GSLMultiMin...)
	// Minimizer algorithm (Migrad, Simplex, Scan...)
	// Strategy - Minuit default 1 is to only compute full Hessian matrix after minimization(?)
	// Print level (verbosity 0...)
	// Tolerance to control iterations
	// Max function calls
	// Max iterations (not used by Minuit)
	// Parameter errors value (??) - default 1 for chi2, 0.5 for log-likelihood
	// Precision in evaluation of minimization. Default double. Should we use float?
	
	// MinimizerOptions::SetMaxFunctionCalls(int )
	// MinimizerOptions::SetTolerance(double )
	// MinimizerOptions::SetPrecision(double )
	//ROOT::Math::MinimizerOptions opt = fitter.Config().MinimizerOptions();
	//opt.SetMaxFunctionCalls(5000);  // XXX like this?
	
//	ROOT::Fit::DataOptions opt;
//	ROOT::Fit::DataRange range(0,5);
//	ROOT::Fit::UnBinData ncap_dt_data(opt, range, x.size());
	
	// convert the data into a suitable 'UnBinData' object
	std::cout<<"making unbinned dataset"<<std::endl;
	ROOT::Fit::UnBinData ncap_dt_data(li9_ntag_dt_vals.size());
	for(auto aval : li9_ntag_dt_vals){   // note: use auto NOT auto&&
		ncap_dt_data.Add(aval);   // can we introduce weights here? FIXME
	}
	
	// DO THE FIT
	std::cout<<"doing the fit"<<std::endl;
	fitter.LikelihoodFit(ncap_dt_data);
	
	// print the results
	ROOT::Fit::FitResult r=fitter.Result();
	r.Print(std::cout);
	
	std::cout<<"unbinned likelihood fit done"<<std::endl;
	
	// TODO err, extract and return parameters for drawing
	li9_ncap_dt_hist_normalised->Draw();
	ncapture_dt_unbinned_like.Draw("same");
	gPad->WaitPrimitive();
	gPad->Clear();
	return true;
}

double FitPurewaterLi9NcaptureDt::ncap_lifetime_loglike(double* x, double* par){
	// shouldn't be trying to fit times outside our range
	if(((*x)<ncap_dtmin) || ((*x)>ncap_dtmax)) return 1e10;
	// shouldn't be trying to fit parameter values outside the valid range
	if((par[0]<0.) || (par[0]>1.)) return 1e10;
	
	// likelihood of an event occurring at a given time, using the specified background fraction
	// being that this is a probability the return value should be normalised 0-1
	// since the integral is fixed, all we can vary here is is the fraction
	// of events that are signal vs background
	
	// rate = C + A*exp(-dt/τ)
	// τ = neutron capture lifetime
	// C = rate of accidental backgrounds
	// A = rate of signal events
	
	// we need to make a PDF, so need to normalise. Calculating the integral we get:
	// integral = C*Δt + A*τ*[exp(-tmax/τ) - exp(-tmin/τ)]
	// the normalized function is obtained by scaling by this integral, turning C->c, A->a
	// a = A / integral, c = C / integral, giving:
	// c*Δt + a*τ*[exp(-tmax/τ) - exp(-tmin/τ)] = 1
	// we can relate a and c:
	// a = (c*Δt - 1 ) / τ*[exp(-tmax/τ) - exp(-tmin/τ)]
	// so now c*Δt is the fraction of background events
	// for a to remain poaitive we should restrict c to the range 0 to (1/Δt)
	// if we define par[0] = c' = c*Δt, then we can fit it directly.
	
	// Now, we have been given a value par[0], so
	// calculate a based on our value of c':
	double timespan = ncap_dtmax-ncap_dtmin;
	double ncap_lifetime = li9_lifetime_secs;
	double a =  (par[0] - 1.) / (ncap_lifetime*(exp(-ncap_dtmax/ncap_lifetime) - exp(-ncap_dtmin/ncap_lifetime)));
	
	// calculate likelihood of an event at time (*x) given a background fraction of c'
	// i.e. evaluate the normalised function, c + a*exp(-dt/τ), at this time
	double thelikelihood = (par[0]/timespan) + a*exp(-(*x)/ncap_lifetime);
	
	// return log-likelihood
	return thelikelihood;
	//return log(thelikelihood);    // doesn't this mess with our normalization?
}

