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
	
	// times of muons (HE+OD occurrances) in this event
	std::vector<int> untaggedMuonTime;
	bool prim_mu=false;
	
	// check the primary trigger bits for HE+OD first, since we already have them
	std::bitset<sizeof(int)*8> triggerID = skhead_.idtgsk;
	if(triggerID.test(1) && triggerID.test(3)){
		untaggedMuonTime.push_back(0);
		prim_mu=true;
	}
	
	// FIXME - for now downstream toolchain can only pair a relic to one muon per event
	// with 2Hz of muons and a -+60s matching window, it's probably sufficient just to match
	// to any muon on the event...? So don't bother searching for more if we already have one
	// TODO switch to config variable to optionally do this search anyway when we support >1 muon per relic
	else {
		
		// get trigger settings from file (get those we're interested in - HE, OD thresholds for this run)
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
		
		// search for pairs of HE+OD within a 100ns window - consider these muons
		int untagged_count=0;
		bool found_prim=false;
		for(int i = 0; i < ntrg; i++){                                                   // loop over triggers found
			Log(m_unique_name+" trigger "+toString(i)+" is of type "
			    +toString(swtrgtbl_.swtrgtype[i])+" at time "+toString(swtrgtbl_.swtrgt0ctr[i]),v_debug,m_verbose);
			if(swtrgtbl_.swtrgtype[i] == 1){                                             // for each HE trigger...
				Log(m_unique_name+" found HE trigger, looking for coincident OD trigger",v_debug,m_verbose);
				for(int j = 0; j < ntrg; j++){                                           // loop over triggers again
					Log("\t"+m_unique_name+" trigger "+toString(i)+" is of type "
					    +toString(swtrgtbl_.swtrgtype[i]),v_debug,m_verbose);
					if(swtrgtbl_.swtrgtype[j] == 3){                                     // looking for OD triggers...
						Log(m_unique_name+" SHE+OD pair with Δt="                        // in time coincidence
						    +toString(swtrgtbl_.swtrgt0ctr[j] - swtrgtbl_.swtrgt0ctr[i]),v_debug,m_verbose);
						if(std::abs(swtrgtbl_.swtrgt0ctr[j] - swtrgtbl_.swtrgt0ctr[i])< coincidence_threshold){
							// skip it if this is the primary trigger
							// assumes primary trigger time is matched to SHE time not OD time...?
							if( prim_mu && (swtrgtbl_.swtrgt0ctr[i]< coincidence_threshold) ){
								found_prim=true;
							} else {
								// n.b. swtrgt0ctr is t0_sub, so time from it0sk. We would need to add it0sk to get it0xsk.
								untaggedMuonTime.push_back(swtrgtbl_.swtrgt0ctr[i]);
								++untagged_count;
							}
							break;
						}
					}
				}
			}
		}
		
		// if the primary trigger had HE+OD bits set, we would expect to have found the primary trigger in our scan
		if(prim_mu && !found_prim){
			// is 100ns a tighter window than used by the normal trigger algorithm perhaps?
			// we haven't accounted for trigger offsets either...
			Log(m_unique_name+" Warning! software trigger scan did not pick up primary muon event!",
			    v_error,m_verbose);
			//Log(m_unique_name+" Trigger bits in this event were: "+GetTriggerNames(skhead_.idtgsk),
			//    v_error,m_verbose);
		}
		
		if(untagged_count!=0){
			// primary trigger says not a muon event, but we found one in subtriggers
			Log(m_unique_name+": "+std::to_string(untagged_count)+" Untagged muons found!",v_debug,m_verbose);
		}
		
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
