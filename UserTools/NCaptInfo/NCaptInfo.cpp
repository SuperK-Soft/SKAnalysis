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
	m_variables.Get("time_match_tolerance",time_match_tolerance); // [ns]
	m_variables.Get("dist_match_tolerance",dist_match_tolerance); // [cm]
	// whether to try to match candidates to non-neutron-capture truth vertices
	m_variables.Get("match_mistags",match_mistags);
	m_variables.Get("likelihood_threshold",likelihood_threshold);
	m_variables.Get("outfilename",outfilename);
	
	// initialisation required to read in
	InitCandidateReader();
	
	// get name of MC truth file used for matching, if applicable
	m_data->CStore.Get("mcparticlesfile",mctruth_file);
	
	// initialise histograms
	if(outfilename!="") MakePlots(0);
	
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
	
	//PrintCandidates();
	
	// fill histograms
	if(outfilename!="") MakePlots(1);
	
	//Log(m_unique_name+": "+toString(candidates.size())+" candidates this entry",v_debug,verbosity);
	
	return true;
}


bool NCaptInfo::Finalise(){
	
	// draw histograms
	if(outfilename!="") MakePlots(2);
	
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
		if(candidate.capture_likelihood_metric<likelihood_threshold) continue;
		
		// for each candidate, loop over all true captures
		// and match to the best one
		int best_match_index=-1;
		//std::cout<<"scanning "<<m_data->NCapturesTrue.size()<<" true captures"<<std::endl;
		for(int truecapi=0; truecapi<m_data->NCapturesTrue.size(); ++truecapi){
			NCapture& truecap = m_data->NCapturesTrue.at(truecapi);
			double* truetime = truecap.GetTime();
			TVector3* truepos =truecap.GetPos();
			if(truetime==nullptr || truepos==nullptr){
				Log(m_unique_name+" Error! True capture "+toString(truecapi)
				   +" returned nullptr for time or position!",v_error,verbosity);
				return false;
			}
			double timediff = (candidate.capture_time - *truetime)/1000.;
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
			double match_metric = std::sqrt(pow(timediff,2.)+pow(posdiff,2.));
			bool make_match=false;
			if(candidate.GetTrueCapture()!=nullptr){
				// see whether this true capture is a better match for this candidate
				double* old_timediff = candidate.GetCaptTerr();
				double* old_posdiff = candidate.GetCaptPosErr();
				if(old_timediff==nullptr || old_posdiff==nullptr){
					Log(m_unique_name+" Error! Existing true match for candidate "+toString(candi)
					   +" returned nullptr for time or position!",v_error,verbosity);
					return false;
				}
				double old_match_metric = std::sqrt(pow(*old_timediff/1000.,2.)+pow(*old_posdiff,2.));
				if(match_metric<old_match_metric){
					// this is a better match
					make_match=true;
					//std::cout<<"found better match metric "<<match_metric
					//         <<" than old one "<<old_match_metric<<std::endl;
				}
			} else {
				// no current match
				make_match=true;
			}
			
			if(make_match){
				//std::cout<<"matching candidate "<<candi<<" to true capture "<<truecapi<<std::endl;
				//std::cout<<"which is:"<<std::endl;
				//m_data->NCapturesTrue.at(truecapi).Print();
				candidate.SetTrueCaptureIdx(truecapi);
				//std::cout<<"this is now matched to true capture:"<<std::endl;
				//m_data->NCapturesTrue.at(truecapi).Print();
				//std::cout<<"which is at "<<&m_data->NCapturesTrue.at(truecapi)<<std::endl;
				//std::cout<<"compared to "<<candidate.GetTrueCapture()<<std::endl;
				//std::cout<<"and printing one more time "<<std::endl;
				//candidate.GetTrueCapture()->Print();
			}
			
			// we'll only set the 'matchType' for candidates passing some qualifying limit.
			Log(m_unique_name+"checking if candidate passes cuts:\n"+
			    +"timediff: "+toString(timediff)+" vs tolerance "+toString(time_match_tolerance)
			    +"\nposdiff: "+toString(posdiff)+" vs tolerance "+toString(dist_match_tolerance),
			    v_debug,verbosity);
			if(timediff>time_match_tolerance) continue;
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
				// skdetsim doesn't record daughter Gds, so probably Gd if it's skdetsim...
				candidate.matchtype = NCaptCandidate::matchType::kUnknownCapture;
			}
			
			/*
			// XXX Note 1
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
	// XXX see Note 1
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
		out_file = new TFile(outfilename.c_str(),"RECREATE");
		if(out_file==nullptr || out_file->IsZombie()){
			Log(m_unique_name+" Error making histogram file "+outfilename,v_error,verbosity);
			return false;
		}
		out_file->cd();
		// a tree for information about candidates
		candtree = new TTree("candtree",m_unique_name.c_str());
		// a tree for information on how well they match their best true n capture, if available
		matchtree = new TTree("matchtree",m_unique_name.c_str());
		
		// note the input files in the TTree header
		TNamed* candfile = new TNamed("candidates_file", candidates_file);
		candtree->GetUserInfo()->Add(candfile);
		
		TNamed* mcfilename = new TNamed("mc_file", mctruth_file.c_str());
		matchtree->GetUserInfo()->Add(mcfilename);
		
		// info about candidates
		std::vector<std::string> civariables{"nrunsk", "nsubsk", "nevsk", "cand_num"};
		std::vector<std::string> cdvariables{"prompt_t","prompt_x", "prompt_y", "prompt_z",
		                                     "cap_t","cap_x", "cap_y", "cap_z",
		                                     "n_travel_t", "n_travel_d","prompt_e",
		                                     "prompt_goodness", "cap_likelihood"};
		// info about matches
		std::vector<std::string> mivariables{"matchtype"};
		std::vector<std::string> mdvariables{"prompt_terr", "prompt_derr", 
		                                     "prompt_xerr", "prompt_yerr", "prompt_zerr",
		                                     "prompt_eerr", "cap_terr", "cap_derr",
		                                     "cap_xerr", "cap_yerr", "cap_zerr",
		                                     "n_travel_terr", "n_travel_derr"};
		
		// make branches in candidate tree
		cibranchvars.reserve(civariables.size());
		for(auto&& avar : civariables){
			cibranchvars[avar] = 0;
			candtree->Branch(avar.c_str(), &cibranchvars[avar]);
		}
		cdbranchvars.reserve(cdvariables.size());
		for(auto&& avar : cdvariables){
			cdbranchvars[avar] = 0;
			candtree->Branch(avar.c_str(), &cdbranchvars[avar]);
		}
		
		// make branches in match tree
		mibranchvars.reserve(mivariables.size());
		for(auto&& avar : mivariables){
			mibranchvars[avar] = 0;
			matchtree->Branch(avar.c_str(), &mibranchvars[avar]);
		}
		mdbranchvars.reserve(mdvariables.size());
		for(auto&& avar : mdvariables){
			mdbranchvars[avar] = 0;
			matchtree->Branch(avar.c_str(), &mdbranchvars[avar]);
		}
		
	} else if(step==1){
		
		// execution - fill branches/histograms
		cibranchvars["nrunsk"] = skhead_.nrunsk;
		cibranchvars["nsubsk"] = skhead_.nsubsk;
		cibranchvars["nevsk"] = skhead_.nevsk;
		
		std::vector<NCaptCandidate>& candidates = m_data->NCaptureCandidates[m_unique_name];
		Log(m_unique_name+" Filling "+toString(candidates.size())+" candidates",v_debug,verbosity);
		for(int cand_i=0; cand_i<candidates.size(); ++cand_i){
			cibranchvars["cand_num"] = cand_i;
			NCaptCandidate& acand = candidates.at(cand_i);
			
			// candidate info
			cdbranchvars["cap_likelihood"] = acand.capture_likelihood_metric;
			cdbranchvars["cap_t"] = acand.capture_time;       // [ns]
			cdbranchvars["cap_x"] = acand.capture_pos.X();    // these have peaks at edges
			cdbranchvars["cap_y"] = acand.capture_pos.Y();    // 
			cdbranchvars["cap_z"] = acand.capture_pos.Z();    // especially z
			// prompt event info (used for start of search, usually)
			LoweCandidate* prompt_event = acand.GetPromptEvent();
			if(prompt_event){
				cdbranchvars["prompt_goodness"] = prompt_event->goodness_metric;
				cdbranchvars["prompt_t"] = prompt_event->event_time;
				cdbranchvars["prompt_x"] = prompt_event->event_pos.X();
				cdbranchvars["prompt_y"] = prompt_event->event_pos.Y();
				cdbranchvars["prompt_z"] = prompt_event->event_pos.Z();
				cdbranchvars["prompt_e"] = prompt_event->event_energy;
				cdbranchvars["n_travel_t"] = (acand.capture_time - prompt_event->event_time);
				cdbranchvars["n_travel_d"] = (acand.capture_pos - prompt_event->event_pos).Mag();
			} else {
				cdbranchvars["prompt_goodness"] = -1;
				cdbranchvars["prompt_t"] = 9999;  // ignore 9999
				cdbranchvars["prompt_x"] = 9999;  // ignore 9999
				cdbranchvars["prompt_y"] = 9999;  // ignore 9999
				cdbranchvars["prompt_z"] = 9999;  // ignore 9999
				cdbranchvars["prompt_e"] = 9999;
				cdbranchvars["n_travel_t"] = 9999;
				cdbranchvars["n_travel_d"] = 9999;
			}
			
			candtree->Fill();
			
			// match info if available
			NCapture* truecap = acand.GetTrueCapture();
			//std::cout<<"candidate "<<cand_i<<" has true cap at: "<<truecap<<std::endl;
			if(truecap){
				mibranchvars["matchtype"] = int(acand.matchtype);
				mdbranchvars["cap_terr"] = (acand.GetCaptTerr() ? *acand.GetCaptTerr() : 0);
				mdbranchvars["cap_derr"] = (acand.GetCaptPosErr() ? *acand.GetCaptPosErr() : 0);
				TVector3* truecaptpos = truecap->GetPos();
				if(truecaptpos){
					mdbranchvars["cap_xerr"] = acand.capture_pos.X() - truecaptpos->X();
					mdbranchvars["cap_yerr"] = acand.capture_pos.Y() - truecaptpos->Y();
					mdbranchvars["cap_zerr"] = acand.capture_pos.Z() - truecaptpos->Z();
				} else {
					mdbranchvars["cap_xerr"] = 0;
					mdbranchvars["cap_yerr"] = 0;
					mdbranchvars["cap_zerr"] = 0;
				}
				double truetraveld = 0;
				if(truecap->NeutronTravelDist(truetraveld)){
					mdbranchvars["n_travel_derr"] = cdbranchvars["n_travel_d"] - truetraveld;
				} else {
					mdbranchvars["n_travel_derr"] = 0;
				}
				double truetravelt = 0;
				if(truecap->NeutronTravelTime(truetravelt)){
					mdbranchvars["n_travel_terr"] = cdbranchvars["n_travel_t"] - truetravelt;
				} else {
					mdbranchvars["n_travel_terr"] = 0;
				}
				
				// also compare prompt events
				// first init to defaults
				mdbranchvars["prompt_terr"] = 0;
				mdbranchvars["prompt_derr"] = 0;
				mdbranchvars["prompt_xerr"] = 0;
				mdbranchvars["prompt_yerr"] = 0;
				mdbranchvars["prompt_zerr"] = 0;
				mdbranchvars["prompt_eerr"] = 0;
				MParticle* truepositron = truecap->GetIBDPositron();
				if(prompt_event!=nullptr && truepositron!=nullptr){
					TVector3* cand_prompt_pos = &prompt_event->event_pos;
					TVector3* true_prompt_pos = truepositron->GetStartPos();
					if(cand_prompt_pos && true_prompt_pos){
						TVector3 diffvec = (*cand_prompt_pos - *true_prompt_pos);
						mdbranchvars["prompt_derr"] = diffvec.Mag();
						mdbranchvars["prompt_xerr"] = diffvec.X();
						mdbranchvars["prompt_yerr"] = diffvec.Y();
						mdbranchvars["prompt_zerr"] = diffvec.Z();
					}
					double* cand_prompt_t = &prompt_event->event_time;
					double* true_prompt_t = truepositron->GetStartTime();
					if(cand_prompt_t && true_prompt_t){
						mdbranchvars["prompt_terr"] = (*cand_prompt_t - (*true_prompt_t));
					}
					double* cand_prompt_e = &prompt_event->event_energy;
					double* true_prompt_e = truepositron->GetStartE();
					if(cand_prompt_e && true_prompt_e){
						mdbranchvars["prompt_eerr"] = (*cand_prompt_e - *true_prompt_e);
					}
				}
			}
			
			matchtree->Fill();
			
		}
		
	} else if(step==2){
		// finalise - write to file
		out_file->cd();
		candtree->Write();
		matchtree->Write();
		
		//make and save histograms
		std::string canvasname = m_unique_name+"_candidates";
		TCanvas* c1 = new TCanvas(canvasname.c_str(),canvasname.c_str(),1200,800);
		std::string savename;
		
		///*
		TH1D h_cap_likelihood("h_cap_likelihood","cap_likelihood;metric;num events",100,0,1);
		TH1D h_cap_t("h_cap_t","Capture Time;T;num events",100,-1000,60000);
		TH1D h_cap_x("h_cap_x","Capture X Position;X;num events",100,-2000,2000);
		TH1D h_cap_y("h_cap_y","Capture Y Position;Y;num events",100,-2000,2000);
		TH1D h_cap_z("h_cap_z","Capture Z Position;Z;num events",100,-2000,2000);
		
		TH1D h_prompt_goodness("h_prompt_goodness","prompt_goodness;metric;num events",100,0,1);
		TH1D h_prompt_t("h_prompt_t","Prompt Time;T;num events",100,-100,2000);
		TH1D h_prompt_x("h_prompt_x","Prompt X Position;X;num events",100,-2000,2000);
		TH1D h_prompt_y("h_prompt_y","Prompt Y Position;Y;num events",100,-2000,2000);
		TH1D h_prompt_z("h_prompt_z","Prompt Z Position;Z;num events",100,-2000,2000);
		TH1D h_prompt_e("h_prompt_e","Prompt Energy;Energy [MeV];num events",100,0,20);
		
		TH1D h_n_travel_d("h_n_travel_d","Neutron Travel D;D;num events",100,0,2000);
		TH1D h_n_travel_t("h_n_travel_t","Neutron Travel T;T;num events",100,0,2000);
		
		TH1D h_cap_terr("h_cap_terr","Capture Time Error;Terr;num events",100,-1000,1000);
		TH1D h_cap_xerr("h_cap_xerr","Capture X Position Error;Xerr;num events",100,0,2000);
		TH1D h_cap_yerr("h_cap_yerr","Capture Y Position Error;Yerr;num events",100,0,2000);
		TH1D h_cap_zerr("h_cap_zerr","Capture Z Position Error;Zerr;num events",100,0,2000);
		TH1D h_cap_derr("h_cap_derr","Capture Total Position Error;Rerr;num events",100,0,2000);
		
		TH1D h_prompt_terr("h_prompt_terr","Prompt Time Error;Terr;num events",100,-1000,2000);
		TH1D h_prompt_xerr("h_prompt_xerr","Prompt X Position Error;Xerr;num events",100,0,2000);
		TH1D h_prompt_yerr("h_prompt_yerr","Prompt Y Position Error;Yerr;num events",100,0,2000);
		TH1D h_prompt_zerr("h_prompt_zerr","Prompt Z Position Error;Zerr;num events",100,0,2000);
		TH1D h_prompt_derr("h_prompt_derr","Prompt Total Position Error;Rerr;num events",100,0,2000);
		TH1D h_prompt_eerr("h_prompt_eerr","Prompt Energy Error;Eerr;num events",100,0,10);
		
		TH2D h_capt_terr_vs_metric("h_capt_terr_vs_metric","Capture Terr vs Likelihood;cap_likelihood;Terr",200,0,1,200,-300,300);
		TH2D h_capt_derr_vs_metric("h_capt_derr_vs_metric","Capture Rerr vs Likelihood;cap_likelihood;Rerr",200,0,1,200,0,2000);
		
		TH2D h_prompt_terr_vs_metric("h_prompt_terr_vs_metric","Prompt Terr vs Goodness;prompt_goodness;Terr",200,0,1,200,-300,300);
		TH2D h_prompt_derr_vs_metric("h_prompt_derr_vs_metric","Prompt Rerr vs Goodness;prompt_goodness;Rerr",200,0,1,200,0,2000);
		
		// also plot against visible energy? need a suitable metric... N50? N200?
		//TH2D h_capt_terr_vs_E("h_capt_terr_vs_E","Terr;Terr;cap_likelihood;N events",100,-10000,10000,100,0,1);
		//TH2D h_capt_derr_vs_E("h_capt_derr_vs_E","Rerr;Rerr;cap_likelihood;N events",100,-10000,10000,100,0,1);
		
		/*
		std::cout<<"candtree: "<<std::endl;
		candtree->Print();
		std::cout<<"matchtree: "<<std::endl;
		matchtree->Print();
		*/
		
		candtree->Draw("cap_likelihood>>h_cap_likelihood","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_t>>h_cap_t","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_x>>h_cap_x","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_y>>h_cap_y","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_z>>h_cap_z","prompt_t!=9999 && cap_t!=9999");
		
		candtree->Draw("prompt_goodness>>h_prompt_goodness","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_t>>h_prompt_t","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_x>>h_prompt_x","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_y>>h_prompt_y","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_z>>h_prompt_z","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_e>>h_prompt_e","prompt_e!=9999 && cap_t!=9999");
		
		candtree->Draw("n_travel_t>>h_n_travel_t","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("n_travel_d>>h_n_travel_d","prompt_t!=9999 && cap_t!=9999");
		
		candtree->AddFriend(matchtree);
		
		candtree->Draw("cap_terr>>h_cap_terr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_xerr>>h_cap_xerr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_yerr>>h_cap_yerr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_zerr>>h_cap_zerr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_derr>>h_cap_derr","prompt_t!=9999 && cap_t!=9999");
		
		candtree->Draw("prompt_terr>>h_prompt_terr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_xerr>>h_prompt_xerr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_yerr>>h_prompt_yerr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_zerr>>h_prompt_zerr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_derr>>h_prompt_derr","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_eerr>>h_prompt_eerr","prompt_t!=9999 && cap_t!=9999");
		
		candtree->Draw("cap_terr:cap_likelihood>>h_capt_terr_vs_metric","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("cap_derr:cap_likelihood>>h_capt_derr_vs_metric","prompt_t!=9999 && cap_t!=9999");
		candtree->Draw("prompt_terr:prompt_goodness>>h_prompt_terr_vs_metric","prompt_t!=9999 && cap_t!=9999 && prompt_goodness >0 && prompt_goodness<1");
		candtree->Draw("prompt_derr:prompt_goodness>>h_prompt_derr_vs_metric","prompt_t!=9999 && cap_t!=9999 && prompt_goodness >0 && prompt_goodness<1");
		
		c1->Divide(2,2);
		c1->cd(1);
		h_cap_t.Draw("goff");
		c1->cd(2);
		h_cap_x.Draw("goff");
		c1->cd(3);
		h_cap_y.Draw("goff");
		c1->cd(4);
		h_cap_z.Draw("goff");
		savename = m_unique_name + "_captures.png";
		c1->SaveAs(savename.c_str());
		
		c1->cd(1);
		h_prompt_t.Draw("goff");
		c1->cd(2);
		h_prompt_x.Draw("goff");
		c1->cd(3);
		h_prompt_y.Draw("goff");
		c1->cd(4);
		h_prompt_z.Draw("goff");
		savename = m_unique_name + "_prompts.png";
		c1->SaveAs(savename.c_str());
		
		c1->cd(1);
		h_cap_likelihood.Draw("goff");
		c1->cd(2);
		h_prompt_goodness.Draw("goff");
		c1->cd(3);
		h_n_travel_d.Draw("goff");
		c1->cd(4);
		h_n_travel_t.Draw("goff");
		savename = m_unique_name + "_extras.png";
		c1->SaveAs(savename.c_str());
		
		c1->Divide(4,2);
		c1->cd(1);
		h_cap_xerr.Draw("goff");
		c1->cd(2);
		h_cap_yerr.Draw("goff");
		c1->cd(3);
		h_cap_zerr.Draw("goff");
		c1->cd(4);
		h_cap_derr.Draw("goff");
		c1->cd(5);
		h_cap_terr.Draw("goff");
		c1->cd(6);
		h_cap_likelihood.Draw("goff");
		c1->cd(7);
		h_capt_terr_vs_metric.Draw("goff colz");
		c1->cd(8);
		h_capt_derr_vs_metric.Draw("goff colz");
		savename = m_unique_name + "_capture_errs.png";
		c1->SaveAs(savename.c_str());
		
		c1->cd();
		c1->Clear();
		c1->Divide(4,2);
		c1->cd(1);
		h_prompt_xerr.Draw("goff");
		c1->cd(2);
		h_prompt_yerr.Draw("goff");
		c1->cd(3);
		h_prompt_zerr.Draw("goff");
		c1->cd(4);
		h_prompt_derr.Draw("goff");
		c1->cd(5);
		h_prompt_terr.Draw("goff");
		c1->cd(6);
		h_prompt_terr_vs_metric.Draw("goff");
		c1->cd(7);
		h_prompt_derr_vs_metric.Draw("goff");
		savename = m_unique_name + "_prompt_errs.png";
		c1->SaveAs(savename.c_str());
		
		h_cap_likelihood.Write();
		h_cap_t.Write();
		h_cap_x.Write();
		h_cap_y.Write();
		h_cap_z.Write();
		h_prompt_goodness.Write();
		h_prompt_t.Write();
		h_prompt_x.Write();
		h_prompt_y.Write();
		h_prompt_z.Write();
		h_n_travel_d.Write();
		h_n_travel_t.Write();
		h_cap_xerr.Write();
		h_cap_yerr.Write();
		h_cap_zerr.Write();
		h_cap_derr.Write();
		h_cap_terr.Write();
		h_capt_terr_vs_metric.Write();
		h_capt_derr_vs_metric.Write();
		h_prompt_xerr.Write();
		h_prompt_yerr.Write();
		h_prompt_zerr.Write();
		h_prompt_derr.Write();
		h_prompt_terr.Write();
		h_prompt_terr_vs_metric.Write();
		h_prompt_derr_vs_metric.Write();
		//*/
		
		/*
		std::cout<<"waiting for user to close canvas "<<canvasname<<std::endl;
		while(gROOT->FindObject(canvasname.c_str())!=nullptr){
			c1->Modified();
			c1->Update();
			gSystem->ProcessEvents();
			std::this_thread::sleep_for(std::chrono::milliseconds(200));
		}
		*/
		// c1 is deleted when canvas is closed, so only delete if not shown & closed
		Log(m_unique_name+": Cleanup",v_debug,verbosity);
		if(gROOT->FindObject(canvasname.c_str())!=nullptr) delete c1;
		
		matchtree->RemoveFriend(candtree);
		candtree->ResetBranchAddresses();
		matchtree->ResetBranchAddresses();
		out_file->Close();
		// hsitograms are deleted when file is closed
		
		delete out_file;
		
	}
	
	return true;
}

