#include "PostLoweReconstructionCuts.h"

PostLoweReconstructionCuts::PostLoweReconstructionCuts():Tool(){}

// see also $SKOFL_ROOT/examples/lowe/lomufit_gd.F or /home/skofl/analysis/lowe/lomufit_gd.F (same file)

bool PostLoweReconstructionCuts::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	
	// see if recording cuts, and if so make the selector
	get_ok = m_variables.Get("selectorName",selectorName);
	if(get_ok){
		// make note of all the cuts we're going to make in the order we're going to apply them
		// AddCut(selectorName, cutname, description)
		m_data->AddCut(selectorName, "BSenergy", "relic candidate bonsai energy",true,8,100);
		m_data->AddCut(selectorName, "ClusfitGoodness", "relic candidate clusfit goodness",true,0.3,10);
		m_data->AddCut(selectorName, "BSgoodness", "relic candidate bonsai goodness",true,0.5,10);
		m_data->AddCut(selectorName, "BSovaQ", "relic candidate bonsai vertex+direction combined goodness",true,0.25,100);
		m_data->AddCut(selectorName, "WallDistance", "direct distance to closest wall",true,200,1E9);
		m_data->AddCut(selectorName, "SNR", "signal to noise ratio",true,0,0.55);
	}
	
	// whether we need to read the lowe reconstruction variables from file via skroot_get_lowe_
	// or whether it's already populated by an upstream Tool (e.g. lfallfit) within this ToolChain
	get_ok = m_variables.Get("getLoweVarsFromFile",getLoweVarsFromFile);
	
	// make_precut.F sets:
	//skoptn_("31,30');  // read header and TQreal
	
	return true;
}


bool PostLoweReconstructionCuts::Execute(){
	
	// only apply these to lowe events before matching
	EventType eventType;
	get_ok = m_data->vars.Get("eventType", eventType);
	if(eventType!=EventType::LowE){
		Log(m_unique_name+" skipping non-lowE event",v_debug,m_verbose);
		return true;
	}
	
	// if required get lowe reconstruction variables into common blocks
	if(getLoweVarsFromFile){
		skroot_get_lowe_(&skheadf_.root_id,
		                 &get_ok,
		                 skroot_lowe_.bsvertex,
		                 skroot_lowe_.bsresult,
		                 skroot_lowe_.bsdir,
		                 skroot_lowe_.bsgood,
		                 &skroot_lowe_.bsdirks,
		                 skroot_lowe_.bseffhit,
		                 &skroot_lowe_.bsenergy,
		                 &skroot_lowe_.bsn50,
		                 &skroot_lowe_.bscossun,
		                 skroot_lowe_.clvertex,
		                 skroot_lowe_.clresult,
		                 skroot_lowe_.cldir,
		                 &skroot_lowe_.clgoodness,
		                 &skroot_lowe_.cldirks,
		                 skroot_lowe_.cleffhit,
		                 &skroot_lowe_.clenergy,
		                 &skroot_lowe_.cln50,
		                 &skroot_lowe_.clcossun,
		                 &skroot_lowe_.latmnum,
		                 &skroot_lowe_.latmh,
		                 &skroot_lowe_.lmx24,
		                 &skroot_lowe_.ltimediff,
		                 &skroot_lowe_.lnsratio,
		                 skroot_lowe_.lsdir,
		                 &skroot_lowe_.spaevnum,
		                 &skroot_lowe_.spaloglike,
		                 &skroot_lowe_.sparesq,
		                 &skroot_lowe_.spadt,
		                 &skroot_lowe_.spadll,
		                 &skroot_lowe_.spadlt,
		                 &skroot_lowe_.spamuyn,
		                 &skroot_lowe_.spamugdn,
		                 skroot_lowe_.posmc,
		                 skroot_lowe_.dirmc,
		                 skroot_lowe_.pabsmc,
		                 skroot_lowe_.energymc,
		                 &skroot_lowe_.darkmc,
		                 &skroot_lowe_.islekeep,
		                 &skroot_lowe_.bspatlik,
		                 &skroot_lowe_.clpatlik,
		                 &skroot_lowe_.lwatert,
		                 &skroot_lowe_.lninfo,
		                 skroot_lowe_.linfo);
		
		if(get_ok!=0){
			Log(m_unique_name+" Error "+toString(get_ok)+" calling skroot_get_lowe_!",v_error,m_verbose);
			return false;
		}
	}
	
	// 4.16 "pre-energy" cut v2; bonsai energy low threshold or reconstruction failed
	// n.b. this is a looser cut than the run-wise energy cut in make_precut, so skip it
	/*
	if((skroot_lowe_.bsenergy < 4.) || (skroot_lowe_.bsenergy>9000)){
		if(!selectorName.empty()) m_data->ApplyCut(selectorName, "BSenergyLoose",skroot_lowe_.bsenergy);
		m_data->vars.Set("Skip", true);
		return true;
	}
	*/
	
	// 1. energy cut
	// n.b. relic_sk4_ana lowfit_sk4 cut is different before run 68671, which had a higher SHE threshold,
	// but since we'll be processing SK-VI+ we will only need the one cut.
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "BSenergy",skroot_lowe_.bsenergy);
	if(skroot_lowe_.bsenergy < 8 || skroot_lowe_.bsenergy>100 ){
		// note the different placement of ApplyCut vs AddPassingEvent
		// the former records the value even if it fails so we can see the cutoff in distribution
		// so needs to be called before we `return`.
		// the latter is for boolean cuts; nothing to record if it fails, its absence is enough.
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed BSenergy cut with value "+toString(skroot_lowe_.bsenergy),v_debug,m_verbose);
		return true;
	}
	
	// 4.13 remove bad clusfit goodness events
	// XXX what about events with qismsk > NHITCUT that bypass clusfit???
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "ClusfitGoodness",skroot_lowe_.clgoodness);
	if(skroot_lowe_.clgoodness<0.3){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed clusfit goodness cut with value "+toString(skroot_lowe_.clgoodness),v_debug,m_verbose);
		return true;
	}
	
	// remove bad bonsai goodness events
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "BSgoodness",skroot_lowe_.bsgood[1]);
	if(skroot_lowe_.bsgood[1]<0.5){
		Log(m_unique_name+": event failed BSgoodness cut with value "+toString(skroot_lowe_.bsgood[1]),v_debug,m_verbose);
		m_data->vars.Set("Skip", true);
		return true;
	}
	
	// 5. bonsai ovaQ fit quality cut (lowfit_sk4)
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "BSovaQ",skroot_lowe_bsovaq);
	if(skroot_lowe_bsovaq < 0.25){    // alias for skroot_lowe_.linfo[25]. calculated in lfallfit.
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed BSovaQ cut with value "+toString(skroot_lowe_bsovaq),v_debug,m_verbose);
		return true;
	}
	
	// 4.14 FV cut
	// XXX FV distance can be configured in WallCut Tool
	skroot_lowe_bswallsk = wallsk_(skroot_lowe_.bsvertex); // XXX seems to already be calculated in lfallfit
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "WallDistance",skroot_lowe_bswallsk);
	if(skroot_lowe_bswallsk < 200){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed WallDistance cut with value "+toString(skroot_lowe_bswallsk),v_debug,m_verbose);
		return true;
	}
	
	// 4.9 remove events with high signal-to-noise ratioq
	// SNR is defined as N(Q<0.5pe)/N(anyQ) - i.e. fraction of low-charge hits
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "SNR",skroot_lowe_.lnsratio);
	if(skroot_lowe_.lnsratio > 0.55){  // aka LoweInfo::lnsratio; calculated in lfallfit
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed SNR cut with value "+toString(skroot_lowe_.lnsratio),v_debug,m_verbose);
		return true;
	}
	
	// 4.10 remove "clustered ATM hits" events...?
	// this cut is commented out in both SKOFL and relic_sk4_ana versions of lf_1st_reduction
	/*
	float atmratio=-1;   // skroot_lowe_.latmh and skroot_lowe_.latmnum are calculated in lfallfit
	if(skroot_lowe_.latmnum > 0) atmratio = float(skroot_lowe_.latmh) / float(skroot_lowe_.latmnum);
	if(atmratio > 0.95){
		if(!selectorName.empty()) m_data->ApplyCut(selectorName, "ATMcluster",atmratio);
		m_data->vars.Set("Skip", true);
		return true;
	}
	*/
	
	// 4.11 remove flasher events
	// again commented out in lf_1st_reduction, "flasher cut removed" since 2006
	/*
	if(skroot_lowe_nflf>1){  // alias for skroot_lowe_.linfo[0]. calculated in lfallfit
		if(!selectorName.empty()) m_data->ApplyCut(selectorName, "Flasher",atmratio);
		m_data->vars.Set("Skip", true);
		return true;
	}
	// see also $SKOFL_ROOT/sklowe/lfflash.F which defines `lfflash_` which appears to be a newer version
	// of `lfflasher_` defined in lf_1st_reduction.... but is not used by anything in $SKOFL_ROOT.
	*/
	
	// 6. effective (back-projected) distance from wall
	// commented out in relic_sk4_ana make_precut.F
	/*
	int fromIDwall=1; // 1=distance from ID wall, anything else= distance from OD wall
	skroot_lowe_poswal = effwallf_(&one, &skroot_lowe_.bsvertex[0], &skroot_lowe_.bsdir[0], skroot_lowe_poswal);
	// skroot_lowe_poswal == skroot_lowe_.linfo[33]. XXX check if we need to recalculate after lfallfitg/33
	
	// n.b. this represents the back-projected position at the wall.
	if(skroot_lowe_poswal < ???){  // seems they weren't sure what cut to place
		if(!selectorName.empty()) m_data->ApplyCut(selectorName, "BSeffwall",skroot_lowe_poswal);
		m_data->vars.Set("Skip", true);
		return true;
	}
	*/
	Log(m_unique_name+": event passed Post-Lowe Reconstruction cuts",v_debug,m_verbose);
	
	return true;
}


bool PostLoweReconstructionCuts::Finalise(){
	
	return true;
}
