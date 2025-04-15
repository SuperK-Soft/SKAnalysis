#include "NCaptInfo_BDT.h"
#include "LoweCandidate.h"

bool NCaptInfo_BDT::InitCandidateReader(){
	
	// tree reader to retreive candidates
	std::string treeReaderName="";
	m_variables.Get("treeReaderName",treeReaderName);
	 if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+" failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	candidates_file = myTreeReader->GetFile()->GetName();
	
	return true;
}

bool NCaptInfo_BDT::GetCandidates(std::vector<NCaptCandidate>& candidates){
	
	// get variables from Tree branches
	// 1D arrays of size np (number of candidates)
	basic_array<float*> metrics;
	basic_array<float*> times;
	basic_array<float*> prompt_times;
	basic_array<float*> xs;
	basic_array<float*> ys;
	basic_array<float*> zs;
	basic_array<float*> bses;
	
	get_ok  = myTreeReader->Get("neutron5",metrics);
	get_ok &= myTreeReader->Get("nvx",xs);
	get_ok &= myTreeReader->Get("nvy",ys);
	get_ok &= myTreeReader->Get("nvz",zs);
	// dt= time of prompt event, dtn=time of capture event
	// hmm, actually both seem to have an offset of about 5-10ns from the true capture time.
	// 'dt' is marginally better, but more in spread than in central value.
	get_ok &= myTreeReader->Get("dt",times);
	get_ok &= myTreeReader->Get("bse",bses);  // is this bonsai energy for the ncapture candidate?
	
	
	// prompt information is from upstream LOWE branch which gets propagated
	LoweInfo* myLowE=nullptr;
	get_ok &= myTreeReader->Get("LOWE",myLowE);
	
	if(!get_ok){
		Log(m_unique_name+": error getting candidates!",v_error,m_verbose);
		return false;
	}
	
	// populate candidates structure
	candidates.resize(times.size());
	for(int icand=0; icand<times.size(); ++icand){
		
		Log(m_unique_name+": Get BDT candidate "+toString(icand),v_debug,m_verbose);
		NCaptCandidate& cand = candidates.at(icand);
		cand.algo = "BDT";
		
		// likelihood metric should be defined such that 1 is high confidence and 0 is no confidence.
		// didn't NTag and the BDT define metrics in the opposite way...? one might need `1-metrics.at(icand)`
		cand.capture_likelihood_metric = metrics.at(icand);
		cand.capture_time = times.at(icand);
		cand.capture_pos = TVector3(xs.at(icand),
		                            ys.at(icand),
		                            zs.at(icand));
		cand.capture_E = bses.at(icand);
		
		// we ought to save some info about the prompt event too, as this is used for the search starting point
		if(myLowE){
			
			//LoweCandidate prompt;
			m_data->LoweCandidates.resize(m_data->LoweCandidates.size()+1);
			LoweCandidate& prompt = m_data->LoweCandidates.back();
			
			prompt.algo = "bonsai";
			prompt.event_pos = TVector3(myLowE->bsvertex[0],
			                            myLowE->bsvertex[1],
			                            myLowE->bsvertex[2]);
			prompt.event_time = myLowE->bsvertex[3];
			prompt.event_energy = myLowE->bsenergy;
			prompt.goodness_metric = myLowE->bsgood[1]; // seems to be the main 'goodness' metric
			std::vector<float> bsgood_vec{myLowE->bsgood, myLowE->bsgood+3};
			prompt.recoVars->Set("bsgood",bsgood_vec);
			prompt.recoVars->Set("bsdirks",myLowE->bsdirks);
			prompt.recoVars->Set("ovaQ",myLowE->linfo[26]);
			
			//m_data->LoweCandidates.push_back(prompt);
			cand.SetPromptEvent(m_data->LoweCandidates.size()-1);
			
		}
		
		// TODO
		// populate feature variables used by the tagging algorithm?
		//cand.featureMap = BStore(true, BSTORE_BINARY_FORMAT);
		/*
		 // feature branches used by bdt
		 N10                     : N10[np]/I
		 N10d                    : N10d[np]/I
		 Nc                      : Nc[np]/I
		 Nback                   : Nback[np]/I
		 N300                    : N300[np]/I
		 trms                    : trms[np]/F
		 trmsdiff                : trmsdiff[np]/F
		 fpdist                  : fpdist[np]/F
		 bpdist                  : bpdist[np]/F
		 fwall                   : fwall[np]/F
		 bwall                   : bwall[np]/F
		 bse                     : bse[np]/F
		 mintrms_3               : mintrms_3[np]/F
		 mintrms_6               : mintrms_6[np]/F
		 Qrms                    : Qrms[np]/F
		 Qmean                   : Qmean[np]/F
		 thetarms                : thetarms[np]/F
		 NLowtheta               : NLowtheta[np]/I
		 phirms                  : phirms[np]/F
		 thetam                  : thetam[np]/F
		 NhighQ                  : NhighQ[np]/I
		 Nlow1                   : Nlow1[np]/I
		*/
		
	}
	
	return get_ok;
	
}
