#include "DataQualityCuts.h"
#include <bitset>

DataQualityCuts::DataQualityCuts():Tool(){}


bool DataQualityCuts::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	
	// these cuts are based on some in lowfit_sk4.F, and others in make_precut.F
	// order is different, but there's no point doing lowe reconstruction if we're
	// just going to cut the event at a later stage
	
	// lowfit_sk4.F sets:
	//skoptn_("31,30,26,25");
	//skbadopt_(23);
	
	// see if recording cuts, and if so make the selector
	get_ok = m_variables.Get("selectorName",selectorName);
	if(get_ok){
		// make note of all the cuts we're going to make in the order we're going to apply them
		// AddCut(selectorName, cutname, description)
		m_data->AddCut(selectorName, "Incomplete", "reject incomplete events",false);
		m_data->AddCut(selectorName, "ID_Off", "reject events with the ID off",false);
		m_data->AddCut(selectorName, "OD_Off", "reject events with the OD off",false);
		m_data->AddCut(selectorName, "SlowData", "reject slowcontrol entries",false);
		m_data->AddCut(selectorName, "RunInfo", "reject runinfo entries",false);
		m_data->AddCut(selectorName, "Spacer", "reject spacer entries",false);
		m_data->AddCut(selectorName, "LEDburst", "reject LED burst entries",false);
	}
	
	return true;
}


bool DataQualityCuts::Execute(){
	
	// 3. bad run cut - now part of TreeReader (FIXME TODO)
	
	// cuts from lf_1st_reduction, called by make_precut, which don't depend on lowe reco variables.
	// if we're going to skip these events, may as well do it early.
	
	// 4. 1st reduction:
	// note: relic_sk4_ana has its own lf_1st_reduction_hz, which only has 2 differences from
	// the standard lf_1st_reduction in $SKOFL_ROOT/lowe/sklowe/lf_1st_reduction.F; it comments out
	// the FV cut and the 50us time cut. However, it then makes the exact same cuts in make_precut.F,
	// which is where it calls lf_1st_reduction_hz. So: no difference.
	
	// the following cuts seem to be related to outright data quality;
	// these events probably just shouldn't be used
	std::bitset<32> triggerID(skhead_.idtgsk);
	std::bitset<32> eventFlags(skhead_.ifevsk);
	
	// 4.4 remove incomplete events
	if(eventFlags.test(EventFlagSKIV::INCOMPLETE_TQ)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed Incomplete cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "Incomplete");
	
	// 4.5 remove ID off events
	if(eventFlags.test(EventFlagSKIV::INNER_DETECTOR_OFF)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed ID_Off cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "ID_Off");
	
	// 4.5 remove OD off events
	if(eventFlags.test(EventFlagSKIV::ANTI_DETECTOR_OFF)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed OD_Off cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "OD_Off");
	
	// skip 'slow data' events (? slow control data?)
	if(eventFlags.test(EventFlagSKIV::INNER_SLOW_DATA)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed SlowData cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "SlowData");
	
	// remove 'run info' events
	if(eventFlags.test(EventFlagSKIV::RUN_INFORMATION)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed RunInfo cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "RunInfo");
	
	// remove spacer events
	if(eventFlags.test(EventFlagSKIV::SPACER_BLOCK)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed Spacer cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "Spacer");
	
	// remove LED burst
	if(eventFlags.test(EventFlagSKIV::LED_BURST_ON)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed LEDburst cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "LEDburst");
	
	// TODO: from relic sk4:
	// reject "test runs"
	// reject "calibration runs" (LINAC?)
	// reject runs < 5mins
	// reject runs started <15 mins after HV recovery
	// reject runs with hardware problems?
	// reject runs with "odd event distriubtions"
	// reject "badly processed events"?
	Log(m_unique_name+": event passed Data Quality cuts",v_debug,m_verbose);
	
	return true;
}


bool DataQualityCuts::Finalise(){
	
	return true;
}
