#include "MuonSearch.h"
#include "Constants.h"

#include "fortran_routines.h"
#include "SK_helper_functions.h"

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
		std::string description="accept events with a coincident HE+OD trigger pair anywhere in the readout window. Coincidence defined as within "+toString(coincidence_threshold)+" ns";
		m_data->AddCut(selectorName, m_unique_name, description);
	}
	
	// convert coincidence_threshold from ns to clock ticks for it0sk and swtrgt0ctr
	coincidence_threshold *= COUNT_PER_NS;
	
	return true;
}


bool MuonSearch::Execute(){
	
	std::bitset<sizeof(int)*8> triggerID = skhead_.idtgsk;
	
	int idetector [32], ithr [32], it0_offset [32],ipret0 [32],ipostt0 [32];
	
	// get trigger settings from file (why bother?)
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
	int ntrg = softtrg_inittrgtbl_(&myHeader->nrunsk, &zero, &one, &max_qb);
	
	std::vector<float> untaggedMuonTime;
	
	// search for pairs of HE+OD within a 100ns window - consider these muons
	for(int i = 0; i < ntrg; i++){                                                   // loop over triggers found
		if(swtrgtbl_.swtrgtype[i] == 1){                                             // for each HE trigger...
			for(int j = 0; j < ntrg; j++){                                           // loop over triggers again
				if( (swtrgtbl_.swtrgtype[j] == 3) &&                                 // looking for OD triggers...
					(abs(swtrgtbl_.swtrgt0ctr[j] - swtrgtbl_.swtrgt0ctr[i])          // within e.g. ~100ns
					     < coincidence_threshold)){
						untaggedMuonTime.push_back(swtrgtbl_.swtrgt0ctr[i]);
					}
				}
			}
		}
	}
	
	// seems redundant, but we can also check the primary trigger
	if(triggerID.test(1) && triggerID.test(3)){
		if(untaggedMuonTime.empty()){
			// worrying if we did not find it...
			Log(m_unique_name+" Warning! software trigger scan did not pick up primary muon event!",
			    v_error,verbosity);
			untaggedMuonTime.push_back(skheadqb_.it0sk);
		} else {
			// this is probably the first entry from our scan
			if(std::abs(skheadqb_.it0sk - untaggedMuonTime.front()) > coincidence_threshold){
				// time difference of >100ns??
				Log(m_unique_name+" Warning! software trigger scan first HE+OD event at "
				   +toString(untaggedMuonTime.front())+" does not line up with primary trigger HE+OD event at "
				   +toString(skheadqb_.it0sk),v_warning,verbosity);
				// i guess we can scan for a more robust check...?
				bool foundit=false;
				for(auto&& atime : untaggedMuonTime){
					if(std::abs(atime - skheadqb_.it0sk) < coincidence_threshold){
						foundit=true;
						Log(m_unique_name+" found match (Δt="+toString(std::abs(atime - skheadqb_.it0sk))
						    +") in later subtrigger",v_warning,verbosity);
					}
				}
				// add it to our list i guess
				if(!foundit) untaggedMuonTime.push_back(skheadqb_.it0sk);
			}
		}
	} else if(!untaggedMuonTime.empty()){
		// primary trigger says not a muon event, but we found one in subtriggers
		Log(m_unique_name+": Untagged muon found!",v_debug,verbosity);
	}
	
	// if we found any muons
	if(!untaggedMuonTime.empty()){
		// flag it for downstream tools (this will not be a relic candidate)
		m_data->vars.Set("newMuon", true);
		// and pass their times
		m_data->vars.Set("muonTimes", untaggedMuonTime);
		// mark this event as 'passing' the cut
		if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, m_unique_name);
	} else {
		m_data->vars.Set("newMuon", false);
		m_data->vars.Set("muonTimes", std::vector<float>{});
	}
	
	return true;
	
}


bool MuonSearch::Finalise(){
	
	return true;
}
