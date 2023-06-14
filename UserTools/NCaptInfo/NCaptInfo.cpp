#include "NCaptInfo.h"
#include "Algorithms.h"

#include <thread>
#include <chrono>

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"

// TODO rename this tool

NCaptInfo::NCaptInfo():Tool(){
	// get the name of the tool from its class name
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}

bool NCaptInfo::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("time_match_tolerance",time_match_tolerance);
	m_variables.Get("dist_match_tolerance",dist_match_tolerance);
	// whether to try to match candidates to non-neutron-capture truth vertices
	m_variables.Get("match_mistags",match_mistags);
	m_variables.Get("likelihood_threshold",likelihood_threshold);
	
	// tree reader is used by derived classes to retreive candidates
	std::string treeReaderName="";
	m_variables.Get("treeReaderName",treeReaderName);
	
	// check the TreeReader
	 if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+" failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	// initialise histograms
	MakePlots(0);
	
	return true;
}


bool NCaptInfo::Execute(){
	
	// clear last event so we don't carry over
	std::vector<NCaptCandidate>& candidates = m_data->NCaptureCandidates[m_unique_name];
	candidates.clear();
	
	// call specialised class method to get information for a given algorithm
	Log(m_unique_name+": Getting candidates",v_debug,verbosity);
	get_ok = GetCandidates(candidates);
	if(not get_ok){
		Log(m_unique_name+": Error getting candidates!",v_error,verbosity);
		return false;
	}
	
	// could add any standard (non-algorithm-specific) information
	//for(NCaptCandidate& cand : candidates){
	//	? anything?
	//}
	
	// match to true captures, if we have them
	MatchToTrueCaptures();
	if(match_mistags) MatchMistags();
	
	PrintCandidates();
	
	// fill histograms
	MakePlots(1);
	
	//Log(m_unique_name+": "+toString(candidates.size())+" candidates this entry",v_debug,verbosity);
	
	return true;
}


bool NCaptInfo::Finalise(){
	
	// draw histograms
	MakePlots(2);
	
	return true;
}

bool NCaptInfo::PrintCandidates(){
	
	std::vector<NCaptCandidate>& candidates = m_data->NCaptureCandidates[m_unique_name];
	std::cout<<"This event contained "<<candidates.size()
	         <<" candidate neutron captures from "<<m_unique_name
	         <<"\n==========================================\n";
	for(int i=0; i<candidates.size(); ++i){
		if(i>0) std::cout<<"------------------------------------------\n";
		std::cout<<"Candidate "<<i<<"\n";
		NCaptCandidate& acand = candidates.at(i);
		acand.Print();
	}
	std::cout<<"=========================================="<<std::endl;
	
	return true;
}

bool NCaptInfo::MatchToTrueCaptures(){
	
	std::vector<NCaptCandidate>& candidates = m_data->NCaptureCandidates[m_unique_name];
	
	// make a map of match quality metrics, based on difference in time and position
	std::vector<std::vector<double> > match_merits(candidates.size(), 
	            std::vector<double>(m_data->NCapturesTrue.size(),999999));
	
	// loop over candidates
	Log(m_unique_name+" matching "+toString(candidates.size())+" candidates",v_debug,verbosity);
	for(int candi=0; candi<candidates.size(); ++candi){
		NCaptCandidate& candidate = candidates.at(candi);
		// we could place a cut on likelihood metric first
		if(candidate.likelihood_metric<likelihood_threshold) continue;
		
		// for each candidate, loop over all true captures
		// and match to the best one
		int best_match_index=-1;
		for(int truecapi=0; truecapi<m_data->NCapturesTrue.size(); ++truecapi){
			NCapture& truecap = m_data->NCapturesTrue.at(truecapi);
			double* truetime = truecap.GetTime();
			TVector3* truepos =truecap.GetPos();
			if(truetime==nullptr || truepos==nullptr){
				Log(m_unique_name+" Error! True capture "+toString(truecapi)
				   +" returned nullptr for time or position!",v_error,verbosity);
				return false;
			}
			double timediff = candidate.capture_time - *truetime;
			double posdiff = (*truepos - candidate.capture_pos).Mag();
			if(verbosity>v_debug){
				std::cout<<"candidate "<<candi<<" has position "<<toString(candidate.capture_pos)
					     <<", vs true cap "<<truecapi<<" which has position "<<toString(*truepos)
					     <<" giving distance "<<posdiff<<" compared to tolerance "
					     <<dist_match_tolerance<<std::endl;
				std::cout<<"candidate "<<candi<<" has time "<<toString(candidate.capture_time)
					     <<", vs true cap "<<truecapi<<" which has time "<<toString(*truetime)
					     <<" giving difference "<<timediff<<" compared to tolerance "
					     <<time_match_tolerance<<std::endl;
			}
			
			// always set the best matching true capture for every candidate so that we can
			// make plots of pos/time error, even for candidates without a qualifying match
			// to compare which is the better match, use a metric based on time and distance
			//  (TODO: and likelihood?)
			// since capture times are ~O(20us) and capture distances are ~O(20cm)
			// let's say add them in quadrature, with 1us equivalent to 1cm
			// and say the better match is the one with the lower sum.
			// also, it doesn't make sense for the candidate time to be before the true capture
			// i.e. timediff > 0 (within some uncertainty from PMT timing resolution and
			// a combination of noise + true hits gives a candidate time a little negative)
			// FIXME should we / how do we incorporate this into the match metric?
			double match_metric = std::sqrt(pow(timediff/1000.,2.)+pow(posdiff,2.));
			bool make_match=false;
			if(candidate.GetTrueCapture()!=nullptr){
				// see whether this true capture is a better match for this candidate
				double* old_timediff = candidate.GetTerr();
				double* old_posdiff = candidate.GetPosErr();
				if(old_timediff==nullptr || old_posdiff==nullptr){
					Log(m_unique_name+" Error! Existing true match for candidate "+toString(candi)
					   +" returned nullptr for time or position!",v_error,verbosity);
					return false;
				}
				double old_match_metric = std::sqrt(pow(*old_timediff/1000.,2.)+pow(*old_posdiff,2.));
				if(match_metric<old_match_metric){
					// this is a better match
					make_match=true;
				}
			} else {
				// no current match
				make_match=true;
			}
			
			if(make_match){
				candidate.SetTrueCaptureIdx(truecapi);
			}
			
			// we'll only set the 'matchType' for candidates passing some qualifying limit.
			if(timediff < 1000 || timediff>time_match_tolerance) continue;
			if(posdiff>dist_match_tolerance) continue;
			
			// try to set a specific match type based on isotope, if known
			// otherwise just set to general 'UnknownCapture'
			if(truecap.GetDaughterNuclide()){
				if(truecap.GetDaughterNuclide()->pdg==100045){
					candidate.matchtype = NCaptCandidate::matchType::kHCapture;
				} else if(truecap.GetDaughterNuclide()->pdg==1000641560 || // Gd155
						  truecap.GetDaughterNuclide()->pdg==1000641580){  // Gd157
					candidate.matchtype = NCaptCandidate::matchType::kGdCapture;
				}
			} else {
				candidate.matchtype = NCaptCandidate::matchType::kUnknownCapture;
			}
			
			/*
			// if we wanted to ensure a matching of 1:1 from candidates to true captures,
			// we would need to check whether the best match for this candidate
			// already has a better match with another candidate.
			// to do this, we first build a map of all possible match metrics,
			// then scan over the map, removing candidate+true capture pairs
			// in descending order of match quality.
			Log(m_unique_name+"passed tolerances: match metric is "+toString(match_metric),v_debug,verbosity);
			match_merits.at(candi).at(truecapi)=match_metric;
			*/
			
		}
	}
	
	/*
	// scan over the map to pull out the best matching pairs (optionally with a threshold on FOM)
	Log(m_unique_name+" scanning for best matches",v_debug,verbosity);
	double merit_threshold=999999; // a candidate must have at least this goodness metric to be matched
	while(true){
		double currentmin=999999.;
		int candi=-1, truecapi=-1;
		// double loop to find the minimum element of the 2D array
		for(int rowi=0; rowi<match_merits.size(); rowi++){
			//std::cout<<"scannig candidate "<<rowi<<std::endl;
			std::vector<double> therow = match_merits.at(rowi);
			std::vector<double>::iterator thisrowsminit = std::min_element(therow.begin(), therow.end());
			if((thisrowsminit!=therow.end())&&((*thisrowsminit)<currentmin)){
				currentmin=*thisrowsminit;
				candi = rowi;
				truecapi = std::distance(therow.begin(),thisrowsminit);
				//std::cout<<"next best metric is candidate "<<candi<<", true capture "<<truecapi<<std::endl;
			}
		}
		// found the best (minimum) match goodness
		Log(m_unique_name+" best match is candidate "+toString(candi)
		     +", truecap "+toString(truecapi),v_debug,verbosity);
		// check if we found any match merits
		if(candi<0 || truecapi<0) break;
		// check if the next best match merit passes our minimum requirement
		double minFOM = match_merits.at(candi).at(truecapi);
		Log(m_unique_name+"corresponding metric is "+toString(minFOM),v_debug,verbosity);
		// if the next best match isn't good enough, we're done
		if(minFOM>merit_threshold) break;
		Log(m_unique_name+" passes merit_threshold check",v_debug,verbosity);
		
		// otherwise, make this match
		NCaptCandidate& candidate = candidates.at(candi);
		candidate.SetTrueCaptureIdx(truecapi);
		
		// try to set match type based on isotope, if known
		NCapture& truecap = m_data->NCapturesTrue.at(truecapi);
		if(truecap.GetDaughterNuclide()){
			if(truecap.GetDaughterNuclide()->pdg==100045){
				candidate.matchtype = NCaptCandidate::matchType::kHCapture;
			} else if(truecap.GetDaughterNuclide()->pdg==1000641560 || // Gd155
			          truecap.GetDaughterNuclide()->pdg==1000641580){  // Gd157
				candidate.matchtype = NCaptCandidate::matchType::kGdCapture;
			}
		} else {
			candidate.matchtype = NCaptCandidate::matchType::kUnknownCapture;
		}
		
		// to retain 1:1 matching, remove this candidate and this true capture
		// from the match merit matric
		match_merits.at(candi).assign(match_merits.at(candi).size(),999999);
		for(auto&& avec : match_merits) avec.at(truecapi) = 999999;
	}
	*/
	
	return true;
}

bool NCaptInfo::MatchMistags(){
	// TODO
	return true;
}

bool NCaptInfo::MakePlots(int step){
	
	if(step==0){
		// initialisation; make histograms
		
		h_likelihood = new TH1D("h_likelihood","Likelihood;metric;num events",100,0,1);
		h_tdiff = new TH1D("h_tdiff","Capture Time Error;Terr;num events",100,-10000,1000);
		h_xdiff = new TH1D("h_xdiff","X Position Error;Xerr;num events",100,0,2000);
		h_ydiff = new TH1D("h_ydiff","Y Position Error;Yerr;num events",100,0,2000);
		h_zdiff = new TH1D("h_zdiff","Z Position Error;Zerr;num events",100,0,2000);
		h_rdiff = new TH1D("h_rdiff","Total Position Error;Rerr;num events",100,0,2000);
		h_tdiff_vs_metric = new TH2D("h_tdiff_vs_metric","Terr;Terr;Likelihood;N events",
			                               100,-10000,10000,100,0,1);
		h_rdiff_vs_metric = new TH2D("h_rdiff_vs_metric","Rerr;Rerr;Likelihood;N events",
			                               100,0,2000,100,0,1);
		// would also be good to plot against energy; true energy? but then we can only plot matched captures.
		// reconstructed energy? can't see a suitable general metric... maybe we could calculate our own?
		// e.g. calculate N200 based on the reconstructed time?
		//TH2D h_tdiff_vs_metric("h_tdiff_vs_E","Terr;Terr;Likelihood;N events",100,-10000,10000,100,0,1);
		//TH2D h_rdiff_vs_metric("h_rdiff_vs_E","Rerr;Rerr;Likelihood;N events",100,-10000,10000,100,0,1);
		
	} else if(step==1){
		
		// execution - fill histograms
		std::vector<NCaptCandidate>& candidates = m_data->NCaptureCandidates[m_unique_name];
		
		for(NCaptCandidate& acand : candidates){
			h_likelihood->Fill(acand.likelihood_metric);
			if(acand.GetTerr()){
				h_tdiff->Fill(*acand.GetTerr()/1000.);  // conver to us?
				h_tdiff_vs_metric->Fill(*acand.GetTerr(),acand.likelihood_metric);
			}
			if(acand.GetPosErr()){
				h_xdiff->Fill(acand.capture_pos.X() - acand.GetTrueCapture()->GetPos()->X());
				h_ydiff->Fill(acand.capture_pos.Y() - acand.GetTrueCapture()->GetPos()->Y());
				h_zdiff->Fill(acand.capture_pos.Z() - acand.GetTrueCapture()->GetPos()->Z());
				h_rdiff->Fill(*acand.GetPosErr());
				h_rdiff_vs_metric->Fill(*acand.GetPosErr(), acand.likelihood_metric);
			}
		}
		
	} else if(step==2){
		// finalise - make and save histograms
		
		std::string histfilename = m_unique_name+"_hists.root";
		TFile* histogram_file = new TFile(histfilename.c_str(),"RECREATE");
		
		std::string canvasname = m_unique_name+"_candidates";
		TCanvas* c1 = new TCanvas(canvasname.c_str(),canvasname.c_str(),1200,800);
		c1->Divide(4,2);
		c1->cd(1);
		h_xdiff->Draw();
		c1->cd(2);
		h_ydiff->Draw();
		c1->cd(3);
		h_zdiff->Draw();
		c1->cd(4);
		h_rdiff->Draw();
		c1->cd(5);
		h_tdiff->Draw();
		c1->cd(6);
		h_likelihood->Draw();
		c1->cd(7);
		h_tdiff_vs_metric->Draw("lego");
		c1->cd(8);
		h_rdiff_vs_metric->Draw("lego");
		
		histogram_file->cd();
		h_xdiff->Write();
		h_ydiff->Write();
		h_zdiff->Write();
		h_rdiff->Write();
		h_tdiff->Write();
		h_likelihood->Write();
		h_tdiff_vs_metric->Write();
		h_rdiff_vs_metric->Write();
		
		/*
		std::cout<<"waiting for user to close canvas "<<canvasname<<std::endl;
		while(gROOT->FindObject(canvasname.c_str())!=nullptr){
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			std::this_thread::sleep_for(std::chrono::milliseconds(200));
		}
		*/
		// c1 is deleted when canvas is closed
		
		histogram_file->Close();
		// hsitograms are deleted when file is closed
		
		delete histogram_file;
		
	}
	
	return true;
}

