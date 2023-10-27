#include "PreLoweReconstructionCuts.h"
#include <bitset>
#include <limits>

PreLoweReconstructionCuts::PreLoweReconstructionCuts():Tool(){}


bool PreLoweReconstructionCuts::Initialise(std::string configfile, DataModel &data){
	
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
		m_data->AddCut(selectorName, "SLE_trigger", "reject SLE triggers",false);
		m_data->AddCut(selectorName, "Periodic_trigger", "reject Periodic triggers",false);
		m_data->AddCut(selectorName, "PeriodicSimple_trigger", "reject Periodic simple triggers",false);
		m_data->AddCut(selectorName, "Calib_trigger", "reject Calibration triggers",false);
		m_data->AddCut(selectorName, "IDlaser_trigger", "reject ID laser triggers",false);
		m_data->AddCut(selectorName, "ODlaser_trigger", "reject OD laser triggers",false);
		m_data->AddCut(selectorName, "T2K_trigger", "reject T2K triggers",false);
		m_data->AddCut(selectorName, "OD_trigger", "reject OD triggers",false);
		m_data->AddCut(selectorName, "LE_trigger", "require LE trigger",false);
		m_data->AddCut(selectorName, "Pedestal_flag", "reject pedestal events",false);
		m_data->AddCut(selectorName, "IQBCalMode_flag", "reject IQBCALMODE==0 events",false);
		m_data->AddCut(selectorName, "NQISK", "reject events with ID charge>2000",true, 0, 2000);
		m_data->AddCut(selectorName, "50usTimeCut", "reject events less than 50us from a previous trigger",true, 50E3, std::numeric_limits<double>::max());
		m_data->AddCut(selectorName, "ODActivity", "reject events with >20 OD hits in 500-->1300ns window",true, 0, 20);
	}
	
	return true;
}


bool PreLoweReconstructionCuts::Execute(){
	
	// basic cuts for relics. For events with untagged muons we'll keep the events
	// to maximise our spallation muon time matching coverage
	EventType eventType;
	get_ok = m_data->vars.Get("eventType", eventType);
	if(eventType==EventType::Muon || eventType==EventType::AFT){
		// what about relic events in a event with a following untagged muon...?
		// reject those as we can't do ntag search?
		Log(m_unique_name+" skipping "+((eventType==EventType::AFT) ? "AFT" : "Muon")+" event",
		    v_debug,m_verbose);
		return true;
	}
	
	// get trigger bits and event flags
	// XXX we can't move these cuts to the SkipTriggers and SkipEventFlags tools unless those tools
	// also acknowledge the 'muonFlag'
	std::bitset<32> triggerID(skhead_.idtgsk);
	std::bitset<32> eventFlags(skhead_.ifevsk);
	
	// 4.1 reject SLE triggers, if slekeep is 2 or 3? XXX what sets skroot_lowe_.islekeep? (aka. LoweInfo::islekeep)
	// pretty sure SLE triggers actually get removed upstream of rfm files, so this is probably redundant
	if(triggerID.test(TriggerType::SLE) && (skroot_lowe_.islekeep==2 || skroot_lowe_.islekeep==3)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed SLE trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "SLE_trigger");
	
	// 4.2 remove periodic trigger
	if(triggerID.test(TriggerType::Periodic)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed periodic trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "Periodic_trigger");
	
	// remove simple periodic trigger
	if(triggerID.test(TriggerType::Periodic_simple)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed simple periodic trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "PeriodicSimple_trigger");
	
	// 4.3 remove calibration trigger
	if(triggerID.test(TriggerType::AFT_or_Cal)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed AFT trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "Calib_trigger");
	
	// remove ID laser trigger
	if(triggerID.test(TriggerType::ID_Laser)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed ID laser trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "IDlaser_trigger");
	
	// remove OD laser trigger
	if(triggerID.test(TriggerType::OD_Laser)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed OD laser trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "ODlaser_trigger");
	
	// remove T2K trigger
	if(triggerID.test(TriggerType::T2K)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed T2K trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "T2K_trigger");
	
	// 4.7 remove OD trigger
	if(triggerID.test(TriggerType::OD_or_Fission)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed OD trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "OD_trigger");
	
	// *require* (at least) LE
	if(!triggerID.test(TriggerType::LE)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed LE trigger requirement cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "LE_trigger");
	
	// 4.4 remove pedestal events
	if(eventFlags.test(EventFlagSKIV::PEDESTAL_ON)){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed Pedestal flag cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "Pedestal_flag");
	
	// reject events with IQBCALMODE==0, whatever that is, from skpdst.h. No C equivalent?
	if(skpdst2_.iqbcalmode==0){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed IQBCalMode trigger cut",v_debug,m_verbose);
		return true;
	}
	if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "IQBCalMode_flag");
	
	// remove very high energy events
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "NQISK", skq_.nqisk);
	if(skq_.nqisk > 2000){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed NQISK>2000 trigger cut",v_debug,m_verbose);
		return true;
	}
	
	// 4.6 remove events within 50us from another previous event
	// skroot_lowe_.ltimediff is not populated in RFM files, but we can calculate it.
	// tdiff gets the PrevT0 branch (which holds the counter_32 and IT0SK values for the last event
	// of each trigger type) and returns the time elapsed (in ns) since the most recent trigger
	// of any kind except SLE and pedestal (trigger IDs 2, 18, 30 - fortran indexing remember)
	// also tdiff of AFT triggers doesn't work as it0xsk is not properly saved (is 0) for them?? XXX
	// tdiff_sub appears to be the same as tdiff, except that it does not call get_prev_t0 first,
	// and does not account for rollover of the SK-VI+ clock.
	// tdiff_muon only supports SK-VI+ root files and, rather than scanning PrevT0 branch,
	// accepts a nevhwsk (aka counter_32) and it0xsk value, returning the time since that event.
	// It also accepts a flag to indicate whether the target event is expected to be before
	// (mode 0) or after (mode 1) the current event (important to account for rollover properly).
	// returned differences should be <0 for target events *following* the current event.
	int unused[3];
	sk::tdiff_(&unused[0], &skroot_lowe_.ltimediff);
	// TODO ltimediff of relics should be updated after MuonSearch to t0_sub(i)/count_per_nsec
	// where t0_sub(i) is subtrigger time of untagged muon
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "50usTimeCut",skroot_lowe_.ltimediff);
	if(skroot_lowe_.ltimediff < 50E3){  // aka LoweInfo::ltimediff, but not set in RFM files
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed 50us time cut",v_debug,m_verbose);
		return true;
	}
	
	// 4.8 remove events with >20 OD hits in [500,1300] ns time window from primary trigger
	// (n.b. time window is hard-coded. For a configurable equivalent use ODCut Tool)
	// whats the justification for this cut?
	lfnhita_(&skroot_lowe_lnahit);  // skroot_lowe_lnahit == linfo[1]; but doesn't seem set in RFM files
	if(!selectorName.empty()) m_data->ApplyCut(selectorName, "ODActivity",skroot_lowe_lnahit);
	if(skroot_lowe_lnahit>20){
		m_data->vars.Set("Skip", true);
		Log(m_unique_name+": event failed OD activity cut",v_debug,m_verbose);
		return true;
	}
	
	// mark this flag as LowE candidate
	m_data->vars.Set("eventType", EventType::LowE);
	
	Log(m_unique_name+": event passed Pre-Lowe Reconstruction cuts",v_debug,m_verbose);
	
	return true;
}


bool PreLoweReconstructionCuts::Finalise(){
	
	return true;
}
