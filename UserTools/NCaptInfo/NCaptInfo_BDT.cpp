#include "NCaptInfo_BDT.h"

bool NCaptInfo_BDT::InitCandidateReader(){
	
	// tree reader to retreive candidates
	std::string treeReaderName="";
	m_variables.Get("treeReaderName",treeReaderName);
	 if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+" failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
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
	basic_array<float*> xs;
	basic_array<float*> ys;
	basic_array<float*> zs;
	
	get_ok  = myTreeReader->Get("neutron5",metrics);
	get_ok &= myTreeReader->Get("nvx",xs);
	get_ok &= myTreeReader->Get("nvy",ys);
	get_ok &= myTreeReader->Get("nvz",zs);
	// dt= time of prompt event, dtn=time of capture event
	get_ok &= myTreeReader->Get("dtn",times);
	
	// prompt information is from upstream LOWE branch which gets propagated
	LoweInfo* myLowE=nullptr;
	get_ok &= myTreeReader->Get("LOWE",myLowE);
	
	if(!get_ok){
		Log(m_unique_name+": error getting candidates!",v_error,verbosity);
		return false;
	}
	
	// populate candidates structure
	candidates.resize(times.size());
	for(int icand=0; icand<times.size(); ++icand){
		
		Log(m_unique_name+": Get BDT candidate "+toString(icand),v_debug,verbosity);
		NCaptCandidate& cand = candidates.at(icand);
		cand.algo = "BDT";
		
		// likelihood metric should be defined such that 1 is high confidence and 0 is no confidence.
		// didn't NTag and the BDT define metrics in the opposite way...? one might need `1-metrics.at(icand)`
		cand.likelihood_metric = metrics.at(icand);
		cand.capture_time = times.at(icand);
		cand.capture_pos = TVector3(xs.at(icand),
		                            ys.at(icand),
		                            zs.at(icand));
		if(myLowE){
			cand.prompt_pos = TVector3(myLowE->bsvertex[0],
			                           myLowE->bsvertex[1],
			                           myLowE->bsvertex[2]);
			cand.prompt_time = myLowE->bsvertex[3];
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
