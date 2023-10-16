/* vim:set noexpandtab tabstop=4 wrap */
#include "PlotMuonDtDlt.h"

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
#include "THStack.h"
#include "TString.h"

#include <cmath>  // fabs


PlotMuonDtDlt::PlotMuonDtDlt():Tool(){}

bool PlotMuonDtDlt::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);            // how verbose to be
	m_variables.Get("outputFile",outputFile);          // where to save data. If empty, current TFile
	m_variables.Get("treeReaderName",treeReaderName);
	
	myTreeReader = m_data->Trees.at(treeReaderName);
	myTreeSelections = m_data->Selectors.at(treeReaderName);
	
	return true;
}

bool PlotMuonDtDlt::Execute(){
	
	// retrieve variables from TreeReader
	GetBranchValues();
	
	// pre muons
	// only consider first muboy muon (only for multi-mu events?)
	std::set<size_t> pre_muboy_first_muons = myTreeSelections->GetPassingIndexes("pre_muon_muboy_i==0");
	for(size_t mu_i : pre_muboy_first_muons){
		Log(m_unique_name+" filling spallation dt and dlt distributions",v_debug+2,m_verbose);
		dlt_vals_pre.at(mu_class[mu_i]).push_back(dlt_mu_lowe[mu_i]);   // FIXME weight by num_pre_muons
		dt_vals_pre.at(mu_class[mu_i]).push_back(dt_mu_lowe[mu_i]);     // FIXME weight by num_pre_muons
		
		// to evaluate systematic on lt cut, apply various dt cuts and see how the lt cut efficiency varies
		// since we're interested in the effect on the spallation sample, which is given by
		// the total - post-muon sample, record both pre- and post- muon samples with various dt cuts
		for(int dt_cut_i=0; dt_cut_i<num_dt_cuts; ++dt_cut_i){
			Log(m_unique_name+" checking nominal dlt cut systematic",v_debug+2,m_verbose);
			if(myTreeSelections->GetPassesCut("pre_mu_dt_cut_"+toString(dt_cut_i),mu_i)){
				Log(m_unique_name+" filling spallation dlt distribution for dt cut "
				            +toString(dt_cut_i),v_debug+2,m_verbose);
				dlt_systematic_dt_cuts_pre.at(dt_cut_i).push_back(dt_mu_lowe[mu_i]);
			}
		}
	}
	// post muons
	std::set<size_t> post_muboy_first_muons = myTreeSelections->GetPassingIndexes("post_muon_muboy_i==0");
	for(size_t mu_i : post_muboy_first_muons){
		Log(m_unique_name+" filling spallation dt and dlt distributions",v_debug+2,m_verbose);
		dlt_vals_post.at(mu_class[mu_i]).push_back(dlt_mu_lowe[mu_i]);   // FIXME weight by num_post_muons
		dt_vals_post.at(mu_class[mu_i]).push_back(dt_mu_lowe[mu_i]);     // FIXME weight by num_post_muons
		
		for(int dt_cut_i=0; dt_cut_i<num_dt_cuts; ++dt_cut_i){
			Log(m_unique_name+" checking nominal dlt cut systematic",v_debug+2,m_verbose);
			if(myTreeSelections->GetPassesCut("post_mu_dt_cut_"+toString(dt_cut_i),mu_i)){
				Log(m_unique_name+" filling spallation dlt distribution for dt cut "
				            +toString(dt_cut_i),v_debug+2,m_verbose);
				dlt_systematic_dt_cuts_post.at(dt_cut_i).push_back(dt_mu_lowe[mu_i]);
			}
		}
	}
	// in Finalise we'll substract the two to get dt and dlt distributions for spallation only.
	// we'll also compare across various dt cuts to get the systematic error on the spallation dlt cut.
	
	return true;
}

bool PlotMuonDtDlt::GetBranchValues(){
	// retrieve variables from branches
	bool success = 
	(myTreeReader->Get("mubstatus", mu_class)) &&
	(myTreeReader->Get("spadt", dt_mu_lowe)) &&
	(myTreeReader->Get("spadlt", dlt_mu_lowe));
	
	return success;
}

bool PlotMuonDtDlt::Finalise(){
	
	// make a new file if given a filename, or if blank check there is a valid file open
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
	
	// subtract distribution of post-muons from pre-muons to obtain lt and dt distributions of mu-lowe pairs
	Log(m_unique_name+" making plots of muon-lowe Dt distributions",v_debug,m_verbose);
	PlotMuonDt();
	Log(m_unique_name+" making plots of muon-lowe Dlt distributions",v_debug,m_verbose);
	PlotMuonDlt();
	
	// measure dlt cut systematic TODO
//	MeasureDltSystematic(); // using variation in dlt_systematic_dt_cuts in bin corresponding to dlt=200cm
	// for each entry in dlt_systematic_dt_cuts_pre, dlt_systematic_dt_cuts_post, subtract post from pre.
	// this gives a set of efficiencies of spallation for varying dts.
	// compare efficiency across various dts in each dlt bin to get the dlt systematic error.
	
	if(fout!=nullptr){
		fout->Close();
		delete fout;
		fout=nullptr;
	}
	
	return true;
}

// =====================================================================
// =====================================================================

bool PlotMuonDtDlt::PlotMuonDlt(){
	
	THStack my_spall_dlts;
	// make histograms of transverse distance to muon for all pre- and post-muons
	// and their difference to extract the spallation distributions
	for(auto&& aclass : constants::muboy_class_to_name){  // we have 5 muboy classifications
		int mu_class_i = int(aclass.first);
		const char* mu_class_name = aclass.second.c_str();
		TH1F ahist_pre(TString::Format("dlt_pre_%d",mu_class_i),
					     "All Pre-Muon to Low-E Transverse Distances",8,0,400);
		for(auto&& aval : dlt_vals_pre.at(mu_class_i)) ahist_pre.Fill(aval);
		ahist_pre.Write();
		
		TH1F ahist_post(TString::Format("dlt_post_%d",mu_class_i),
					     "All Post-Muon to Low-E Transverse Distances",8,0,400);
		for(auto&& aval : dlt_vals_post.at(mu_class_i)) ahist_post.Fill(aval);
		ahist_post.Write();
		
		TH1F* dlt_hist_spall = (TH1F*)ahist_pre.Clone(TString::Format("dlt_spall_%s",mu_class_name));
		dlt_hist_spall->Reset(); // clear entries
		// subtract post from pre
		dlt_hist_spall->Add(&ahist_pre,&ahist_post,1,-1);
		std::cout<<"subtracting "<<ahist_post.GetEntries()<<" post entries from "<<ahist_pre.GetEntries()
				 <<" pre entries results in "<<dlt_hist_spall->GetEntries()<<" ("
				 <<(ahist_pre.GetEntries()-ahist_post.GetEntries())<<") entries"<<std::endl;
		dlt_hist_spall->Scale(1./dlt_hist_spall->Integral());
		// XXX we seem to lose half our entries here....??
		dlt_hist_spall->Write();
		dlt_hist_spall->SetDirectory(0);
		my_spall_dlts.Add(dlt_hist_spall);
	}
	
	// create and add the paper versions
	PlotPaperDlt(my_spall_dlts);
	
	my_spall_dlts.Write("spall_dls_stack");
	
	// cleanup
	my_spall_dlts.GetHists()->Delete();
	
	return true;
}

// =====================================================================
// =====================================================================

bool PlotMuonDtDlt::PlotMuonDt(){
	// make histograms of time difference to muon for all pre- and post-muons
	// and their difference to extract spallation distributions.
	
	// colour match to the existing paper
	// muboy_class{ misfit=0, single_thru_going=1, single_stopping=2, multiple_mu=3, also_multiple_mu=4, corner_clipper=5};
	std::map<std::string, EColor> class_colours{
		{"misfit",kBlue},
		{"single_thru_going",kBlack},
		{"single_stopping",kMagenta},
		{"multiple_mu",kRed},
		{"also_multiple_mu",kRed},
		{"corner_clipper",kWhite}
	}; // corner clippers not shown...
	
	// plot over two ranges: full range (30s) and paper plotted range (0.25s)
	std::vector<float> dt_range_full{30,0.25};
	THStack my_spall_dts;
	for(auto&& dt_max : dt_range_full){
		for(auto&& aclass : constants::muboy_class_to_name){  // we have 5 muboy classifications
			int mu_class_i = int(aclass.first);
			const char* mu_class_name = aclass.second.c_str();
			// need to take the fabs of the time so time 0 is in bin 0 for both pre- and post-
			// in order to be able to subtract the bin counts.
			TH1F ahist_pre(TString::Format("dt_pre_%d_%0.2f",mu_class_i,dt_max),
							 "All Pre-Muon to Low-E Time Differences",10,0,dt_max);
			for(auto&& aval : dt_vals_pre.at(mu_class_i)) ahist_pre.Fill(fabs(aval));
			ahist_pre.Write();
			
			TH1F ahist_post(TString::Format("dt_post_%d_%0.2f",mu_class_i,dt_max),
							 "All Post Muon to Low-E Time Differences",10,0,dt_max);
			for(auto&& aval : dt_vals_post.at(mu_class_i)) ahist_post.Fill(aval);
			ahist_post.Write();
			
			TH1F* dt_hist_spall = 
				(TH1F*)ahist_pre.Clone(TString::Format("dt_spall_%s_%0.2f",mu_class_name,dt_max));
			dt_hist_spall->Reset(); // clear entries
			// subtract post from pre
			dt_hist_spall->Add(&ahist_pre,&ahist_post,1,-1);
			dt_hist_spall->SetLineColor(class_colours.at(mu_class_name));
			dt_hist_spall->Scale(1./dt_hist_spall->Integral());
			dt_hist_spall->Write();
			dt_hist_spall->SetDirectory(0);
			if(dt_max<20.f) my_spall_dts.Add(dt_hist_spall);  // add shorter axis ones for comparison
		}
	}
	
	// create and add the paper ones for overlay
	PlotPaperDt(my_spall_dts);
	
	// write the stack to file so we can easily see all histos
	my_spall_dts.Write("spall_dts_stack");
	
	// TH1::Clone creates a copy that we own, so are responsible for cleanup.
	// THStacks do not take ownership of their histos, and we can't make them.
	// The way to cleanup is therefore:
	my_spall_dts.GetHists()->Delete();
	
	return true;
}


bool PlotMuonDtDlt::PlotPaperDt(THStack& ourplots){
	// get the first of our plots to use as a base to define binning
	TH1F* basehist = (TH1F*)ourplots.GetHists()->At(0);
	
	TH1F* h_paper_dt_mu_nonfitted = (TH1F*)basehist->Clone("paper_dt_nonfitted");
	h_paper_dt_mu_nonfitted->Reset();
	// digitized from paper. I didn't do all series' as most overlapped.
	std::vector<double> paper_dt_mu_nonfitted
		{ 0.5458, 0.3111, 0.1014, 0.0403, 0.0014, 0.0806, 0.0694, 0.0111, 0.0014, 0.0056};
	for(int i=0; i<paper_dt_mu_nonfitted.size(); ++i){
		h_paper_dt_mu_nonfitted->SetBinContent(i+1, paper_dt_mu_nonfitted.at(i));
	}
	h_paper_dt_mu_nonfitted->SetLineColor(kBlue);
	h_paper_dt_mu_nonfitted->SetLineStyle(2);
	
	TH1F* h_paper_dt_mu_single = (TH1F*)basehist->Clone("paper_dt_single");
	h_paper_dt_mu_single->Reset();
	std::vector<double> paper_dt_mu_single
		{ 0.5806, 0.2153, 0.0903, 0.0403, 0.0236, 0.0153, 0.0111, 0.0083, 0.0097, 0.0069 };
	for(int i=0; i<paper_dt_mu_single.size(); ++i){
		h_paper_dt_mu_single->SetBinContent(i+1, paper_dt_mu_single.at(i));
	}
	h_paper_dt_mu_single->SetLineColor(kBlack);
	h_paper_dt_mu_single->SetLineStyle(2);
	
	/*
	// put them into a stack
	THStack paper_dts;
	paper_dts.Add(h_paper_dt_mu_nonfitted);
	paper_dts.Add(h_paper_dt_mu_single);
	// write them out to file
	paper_dts.Write("spall_dts_paper");
	// cleanup
	paper_dts.GetHists()->Delete();
	*/
	
	// add to our stack for easy combined plotting
	ourplots.Add(h_paper_dt_mu_nonfitted);
	ourplots.Add(h_paper_dt_mu_single);
	
	return true;
}

bool PlotMuonDtDlt::PlotPaperDlt(THStack& ourplots){
	
	// get first one as a base for defining axes
	TH1F* basehist = (TH1F*)ourplots.GetHists()->At(0);
	
	// digitized from paper
	std::vector<double> paper_dl_multiple
		{ 0.1414, 0.2303, 0.1980, 0.1465, 0.1121, 0.0798, 0.0596, 0.0404 };
	TH1F* h_paper_dl_multiple = (TH1F*)basehist->Clone("paper_dl_multiple");
	h_paper_dl_multiple->Reset();
	for(int i=0; i<paper_dl_multiple.size(); ++i){
		h_paper_dl_multiple->SetBinContent(i+1,paper_dl_multiple.at(i));
	}
	h_paper_dl_multiple->SetLineColor(kRed);
	h_paper_dl_multiple->SetLineStyle(2);
	
	std::vector<double> paper_dl_single
		{ 0.3636, 0.3707, 0.1485, 0.0636, 0.0283, 0.0141, 0.0121, 0.0071 };
	TH1F* h_paper_dl_single = (TH1F*)basehist->Clone("paper_dl_single");
	h_paper_dl_single->Reset();
	for(int i=0; i<paper_dl_single.size(); ++i){
		h_paper_dl_single->SetBinContent(i+1,paper_dl_single.at(i));
	}
	h_paper_dl_single->SetLineColor(kBlack);
	h_paper_dl_single->SetLineStyle(2);
	
	std::vector<double> paper_dl_stopping
		{ 0.4161, 0.3242, 0.0455, 0.0414, 0.0707, 0.0434, 0.0263, 0.0303 };
	TH1F* h_paper_dl_stopping = (TH1F*)basehist->Clone("paper_dl_stopping");
	h_paper_dl_stopping->Reset();
	for(int i=0; i<paper_dl_stopping.size(); ++i){
		h_paper_dl_stopping->SetBinContent(i+1,paper_dl_stopping.at(i));
	}
	h_paper_dl_stopping->SetLineColor(kMagenta);
	h_paper_dl_stopping->SetLineStyle(2);
	
	/*
	// put them into a stack
	THStack paper_dls;
	paper_dls.Add(h_paper_dl_multiple);
	paper_dls.Add(h_paper_dl_single);
	paper_dls.Add(h_paper_dl_stopping);
	// write them out to file
	paper_dls.Write("spall_dls_paper");
	// cleanup
	paper_dls.GetHists()->Delete();
	*/
	
	// add them to our stack for combined plotting
	ourplots.Add(h_paper_dl_multiple);
	ourplots.Add(h_paper_dl_single);
	ourplots.Add(h_paper_dl_stopping);
	
	return true;
}
