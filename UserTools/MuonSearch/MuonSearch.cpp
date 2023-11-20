#include "MuonSearch.h"
#include "Constants.h"

#include "fortran_routines.h"
#include "softtrg_tbl.h" // FIXME we can't put this in `fortran_routines.h`
// because someone forgot to declare the structs in it as 'extern'!!
// TODO run fh2h.pl on $SKOFL_ROOT/inc/softtrg_tblF.h and put it in $SKOFL_ROOT/inc

#include <bitset>

MuonSearch::MuonSearch():Tool(){}

// Applies the software trigger to search for coincident HE+OD triggers
// if found, marks the event as containing one or muons.

bool MuonSearch::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	m_variables.Get("verbosity",m_verbose);
	m_variables.Get("coincidence_threshold",coincidence_threshold);
	
	// add a cut to the selector if being used
	get_ok = m_variables.Get("selectorName", selectorName);
	if(get_ok){
		std::string description="accept events with a coincident HE+OD trigger pair anywhere"
		                        " in the readout window. Coincidence defined as within "
		                        +toString(coincidence_threshold)+" ns";
		m_data->AddCut(selectorName, m_unique_name, description, false);
	}
	
	// convert coincidence_threshold from ns to clock ticks for it0sk and swtrgt0ctr
	coincidence_threshold *= COUNT_PER_NSEC;
	
	return true;
}


bool MuonSearch::Execute(){
	
	
	// skip AFTs after a relic
	EventType lastEventType = eventType;
	m_data->vars.Get("eventType", eventType);
	/* don't see a compelling reason we shouldn't scan AFTs after LowE's for untagged muons
	if(eventType==EventType::AFT && lastEventType==EventType::LowE){
		Log(m_unique_name+" skipping AFT after LowE event",v_debug,m_verbose);
		return true;
	}
	*/
	
	// get trigger settings from file (why bother?)
	int idetector [32], ithr [32], it0_offset [32],ipret0 [32],ipostt0 [32];
	softtrg_get_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);
	
	// disable all triggers except 1 (HE) and 3 (OD) by setting threshold to 100k and window size to 0
	for(int i = 0; i < 32; i++){
		if(i != 1 && i != 3){
			ithr[i] = 100000;
			it0_offset[i]=0;
			ipret0[i]=0;
			ipostt0[i]=0;
		}
	}
	// pass to the software trigger algorithm
	softtrg_set_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);
	
	// call softtrg_inittrgtbl_ to populate the swtrgtbl_ common block.
	int max_qb = 1280;
	int one = 1;
	int zero = 0;
	int ntrg = softtrg_inittrgtbl_(&skhead_.nrunsk, &zero, &one, &max_qb);
	
	Log(m_unique_name+" found "+toString(ntrg)+" software triggers...",v_debug,m_verbose);
	
	std::vector<int> untaggedMuonTime;
	
	// search for pairs of HE+OD within a 100ns window - consider these muons
	for(int i = 0; i < ntrg; i++){                                                   // loop over triggers found
		Log(m_unique_name+" trigger "+toString(i)+" is of type "
		    +toString(swtrgtbl_.swtrgtype[i]),v_debug,m_verbose);
		if(swtrgtbl_.swtrgtype[i] == 1){                                             // for each HE trigger...
			Log(m_unique_name+" found HE trigger, looking for coincident OD trigger",v_debug,m_verbose);
			for(int j = 0; j < ntrg; j++){                                           // loop over triggers again
				Log("\t"+m_unique_name+" trigger "+toString(i)+" is of type "
				    +toString(swtrgtbl_.swtrgtype[i]),v_debug,m_verbose);
				if(swtrgtbl_.swtrgtype[j] == 3){                                     // looking for OD triggers...
					Log(m_unique_name+" SHE+OD pair with Δt="                        // in time coincidence
					    +toString(swtrgtbl_.swtrgt0ctr[j] - swtrgtbl_.swtrgt0ctr[i]),v_debug,m_verbose);
					if(abs(swtrgtbl_.swtrgt0ctr[j] - swtrgtbl_.swtrgt0ctr[i])< coincidence_threshold){
						// swtrgt0ctr is t0_sub, so time from it0sk. We would need to add it0sk to get it0xsk.
						untaggedMuonTime.push_back(swtrgtbl_.swtrgt0ctr[i]);
					}
				}
			}
		}
	}
	
	// seems redundant, but we can also check the primary trigger
	std::bitset<sizeof(int)*8> triggerID = skhead_.idtgsk;
	if(triggerID.test(1) && triggerID.test(3)){
		if(untaggedMuonTime.empty()){
			// worrying if we did not find it...
			Log(m_unique_name+" Warning! software trigger scan did not pick up primary muon event!",
			    v_error,m_verbose);
			Log(m_unique_name+" Trigger bits in this event were: "+GetTriggerNames(skhead_.idtgsk),
			    v_error,m_verbose);
			untaggedMuonTime.push_back(0);
		} else {
			// this is probably the first entry from our scan
			if(std::abs(untaggedMuonTime.front()) > coincidence_threshold){
				// time difference of >100ns from the primary trigger??
				Log(m_unique_name+" Warning! software trigger scan first HE+OD event at "
				   +toString(untaggedMuonTime.front())+" does not line up with primary trigger HE+OD event?"
				   " (why was this not picked up by online trigger?)",v_warning,m_verbose);
				// i guess we can scan for a more robust check...?
				// see if another match lines up with the primary trigger
				bool foundit=false;
				for(auto&& atime : untaggedMuonTime){
					if(std::abs(atime) < coincidence_threshold){
						foundit=true;
						Log(m_unique_name+" found match (Δt="+toString(std::abs(atime - skheadqb_.it0sk))
						    +") in later subtrigger",v_warning,m_verbose);
					}
				}
				// add it to our list i guess
				if(!foundit) untaggedMuonTime.push_back(0);
			}
		}
	} else if(!untaggedMuonTime.empty()){
		// primary trigger says not a muon event, but we found one in subtriggers
		Log(m_unique_name+": Untagged muon found!",v_debug,m_verbose);
	}
	
	// if we found any muons
	if(!untaggedMuonTime.empty()){
		// flag it for downstream tools (this will not be a relic candidate)
		eventType = EventType::Muon;
		m_data->vars.Set("eventType", eventType);
		// and pass their times
		m_data->CStore.Set("muonTimes", untaggedMuonTime);
		// mark this event as 'passing' the cut
		if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, m_unique_name);
	}
	
	Log(m_unique_name+" Found "+toString(untaggedMuonTime.size())+" muons in this event",v_debug,m_verbose);
	
	
	return true;
	
}


bool MuonSearch::Finalise(){
	
	return true;
}
