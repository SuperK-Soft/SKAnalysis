#include "NCaptInfo_NTag.h"
#include "LoweCandidate.h"

bool NCaptInfo_NTag::InitCandidateReader(){
	
	// tree reader to retreive candidates
	std::string candidateTreeReaderName="";
	m_variables.Get("candidateTreeReaderName",candidateTreeReaderName);
	 if(m_data->Trees.count(candidateTreeReaderName)==0){
		Log(m_unique_name+" failed to find TreeReader "+candidateTreeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		candidatesTreeReader = m_data->Trees.at(candidateTreeReaderName);
	}
	
	candidates_file = candidatesTreeReader->GetFile()->GetName();
	
	// we need another tree reader to access the prompt vertex info
	std::string variablesTreeReaderName="";
	m_variables.Get("variablesTreeReaderName",variablesTreeReaderName);
	 if(m_data->Trees.count(variablesTreeReaderName)==0){
		Log(m_unique_name+" failed to find TreeReader "+variablesTreeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		variablesTreeReader = m_data->Trees.at(variablesTreeReaderName);
	}
	
	return true;
}

bool NCaptInfo_NTag::GetCandidates(std::vector<NCaptCandidate>& candidates){
	
	// get variables from Tree branches
	std::vector<float>* times = nullptr;
	std::vector<float>* metrics = nullptr;
	std::vector<float>* xs = nullptr;
	std::vector<float>* ys = nullptr;
	std::vector<float>* zs = nullptr;
	std::vector<float>* qsums = nullptr;
	
	get_ok  = candidatesTreeReader->Get("TMVAOutput",metrics);
	get_ok &= candidatesTreeReader->Get("TrmsFitVertex_X",xs);
	get_ok &= candidatesTreeReader->Get("TrmsFitVertex_Y",ys);
	get_ok &= candidatesTreeReader->Get("TrmsFitVertex_Z",zs);
	get_ok &= candidatesTreeReader->Get("TRMS",times);       // RMS  of candidate hit times
	//get_ok &= candidatesTreeReader->Get("ReconCT",times);  // mean of candidate hit times
	get_ok &= candidatesTreeReader->Get("QSum",qsums);
	
	TVector3* prompt_vtx=nullptr;
	float prompt_t;
	get_ok &= variablesTreeReader->Get("prompt_vertex",prompt_vtx);
	get_ok &= variablesTreeReader->Get("geant_t0",prompt_t);
	
	if(!get_ok){
		Log(m_unique_name+": error getting candidates!",v_error,verbosity);
		return false;
	}
	
	// populate candidates structure
	candidates.resize(times->size());
	for(int icand=0; icand<times->size(); ++icand){
		
		Log(m_unique_name+": Get NTag candidate "+toString(icand),v_debug,verbosity);
		NCaptCandidate& cand = candidates.at(icand);
		cand.algo = "NTag";
		
		cand.capture_likelihood_metric = metrics->at(icand);
		cand.capture_time = times->at(icand)*1000.;     // XXX why *1000?
		cand.capture_pos = TVector3(xs->at(icand),
		                            ys->at(icand),
		                            zs->at(icand));
		cand.capture_E = qsums->at(icand); // XXX suitable variable?
		
		// save information about the reconstructed prompt event used as a starting point?
		// TODO
		LoweCandidate prompt;
		prompt.event_pos = *prompt_vtx;
		prompt.event_time = prompt_t;
		/*
		// where are these? surely they're stored?
		prompt.event_energy = 0;
		prompt.goodness_metric = 0;
		
		// any other variables we want to save
		cand.recoVars.Set("bsgood",std::vector<float>{myLowE->bsgood,myLowE->bsgood+3});
		cand.recoVars.Set("bsdirks",myLowE->bsdirks);
		cand.recoVars.Set("ovaQ",myLowE->linfo[26]);
		*/
		m_data->LoweCandidates.push_back(prompt);
		cand.SetPromptEvent(m_data->LoweCandidates.size()-1);
		
		// TODO
		// populate feature variables used by the tagging algorithm?
		//cand.featureMap = BStore(true, BSTORE_BINARY_FORMAT);
		/*
		 // other branches in this tree
		 AngleMean       = (vector<float>*)0x17fb680
		 AngleSkew       = (vector<float>*)0x17fb860
		 AngleStdev      = (vector<float>*)0x234d440
		 Beta1           = (vector<float>*)0x240b140
		 Beta2           = (vector<float>*)0x24137e0
		 Beta3           = (vector<float>*)0x241be90
		 Beta4           = (vector<float>*)0x2424540
		 Beta5           = (vector<float>*)0x242cbf0
		 DWall           = (vector<float>*)0x243d950
		 DWallMeanDir    = (vector<float>*)0x2446000
		 DWall_n         = (vector<float>*)0x244e6b0
		 N1300           = (vector<float>*)0x2456d60
		 N200            = (vector<float>*)0x245f410
		 N50             = (vector<float>*)0x2467ac0
		 NHits           = (vector<float>*)0x2470170
		 QSum            = (vector<float>*)0x2478820
		 ThetaMeanDir    = (vector<float>*)0x249a2e0
		 decay_e_like    = (vector<float>*)0x24a2990
		 prompt_nfit     = (vector<float>*)0x24ab040
		*/
		
	}
	
	return get_ok;
	
}
