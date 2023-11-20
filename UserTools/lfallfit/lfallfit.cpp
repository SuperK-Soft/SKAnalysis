/* vim:set noexpandtab tabstop=4 wrap */
#include "lfallfit.h"

#include "Algorithms.h"
#include "Constants.h"
#include "MTreeReader.h"

#include <string>
#include <vector>
#include <iostream>

// declarations and #includes for SK fortran routines
#include "fortran_routines.h"

lfallfit::lfallfit():Tool(){}

bool lfallfit::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);
	m_variables.Get("readerName",readerName);
	m_variables.Get("StepsToPerform",flag_skip);
	m_variables.Get("checkEventType",checkEventType);
	
	// check the upstream TreeReader name
	if(m_data->Trees.count(readerName)==0){
		Log("Failed to find TreeReader "+readerName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	// we need the associated LUN number for skroot_* functions
	lun = m_data->GetLUN(readerName);
	// and to check if this file is MC or not, so we can use the appropriate lf_allfit version
	MC = m_data->Trees.at(readerName)->GetMCFlag();
	
	// skip clusfit if we have too many hits. get the limit for 'too many hits'
	get_ok = m_variables.Get("hitLimitForClusfit",max_nqisk_for_clusfit);
	if(!get_ok){
		// fall back to relic sk4 defaults...
		max_nqisk_for_clusfit = (MC) ? 800 : 1000;
	}
	
	// initialize bonsai
	m_data->BonsaiInit();
	
	return true;
}


bool lfallfit::Execute(){
	
	// if applicable, check for a lowe event flag
	if(checkEventType){
		EventType eventType;
		get_ok = m_data->vars.Get("eventType",eventType);
		if(!get_ok || eventType!=EventType::LowE){
			return true;
		}
	}
	
	Log(m_unique_name+" Passed event type check, executing...",v_debug,m_verbose);
	if((nread%10000)==0){
		Log(m_unique_name+" read loop "+toString(nread)+", current run "
		    +toString(skhead_.nrunsk),v_message,m_verbose);
	}
	++nread;
	
	// clear variables for low & mu fitters
	Log(m_unique_name+" clearall",v_debug,m_verbose);
	// XXX n.b. lfnewfit would seem to indicate we may only need to do this on change of run/subrun
	// and that (perhaps after doing so??) we may need to re-set the darkrate, watert and bad channels?
	lfclear_all_();
	
	// fetch latest water transparency, which may have changed if we changed run.
	float watert;
	m_data->vars.Get("watert",watert);
	
	/* the variants of lfallfit have mostly the same signature:
		inputs:
			watert     water transparency for fit (will be stored in lwatert)
			nhitcut    max. nqisk for clusfit
			flag_skip  ==0  clear variables, then do everything
				       ==1,others  just recalculate energy and related variables
				       ==2  recalculate energy, clik, msg and related variables
			flag_log   LOGLV in findconst
				       1) Lots of output
				       2) only prints filenames
				       3) only prints when not found
				       4) do not print
				       5) only prints when found
		outputs:
			lfflag is output: 0=fit success, -1=too many hits for lowe, -2: no hits (nqisk==0)
	*/
	
	int lfflag;
	int flag_log;
	switch(m_verbose){
		case 0:
		case 1: flag_log = 3; break; // error or warn: print only on error
		case 2:
		case 3: flag_log = 2; break; // message/debug: print filenames (of what?)
		default: flag_log = 1;       // high verbosity: enable lfallfit high verb.
	}
	
	Log(m_unique_name+" calling lfallfit with MC "+toString(MC)
	    +", geometry "+toString(skheadg_.sk_geometry)+", flag_skip "+toString(flag_skip)
	    +", nhitcut "+toString(max_nqisk_for_clusfit)+" and water transp "
	    +toString(watert), v_debug,m_verbose);
	if(MC){
		switch (skheadg_.sk_geometry) {
			case 1:  //lfallfit_sk1_mc_(&watert, &max_nqisk_for_clusfit, &lfflag);    // does not exist
			case 2:  //lfallfit_sk2_mc_(&watert, &max_nqisk_for_clusfit, &lfflag);    // does not exist
			case 3:  //lfallfit_sk3_mc_(&watert, &max_nqisk_for_clusfit, &lfflag);    // does not exist
				Log(m_unique_name+": Error! lfallfit does not exist for SK-"+toString(skheadg_.sk_geometry)
				    +" MC!",v_error,m_verbose);
				return false;
			case 4: {
				// which lfallfit_sk4 would you prefer...?
				//lfallfit_sk4_mc_(&watert, &max_nqisk_for_clusfit, &lfflag);
				//lfallfit_sk4_final_qe41_mc_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				lfallfit_sk4_final_qe43_mc_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				//lfallfit_sk4_gain_corr_mc_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			case 5: {
				lfallfit_sk5_mc_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			case 6: {
				lfallfit_sk6_mc_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			case 7: {
				// FIXME doesn't yet exist; for now use SK-6
				// not sure if this will work: it doesn't for data....
				Log(m_unique_name+" Warning! Calling lfallfit_sk6_mc for SK geometry 7"
				    ", because no lfallfit_sk7_mc exists at time of writing",v_warning,m_verbose);
				lfallfit_sk6_mc_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			default: {
				Log(m_unique_name+": Error! lfallfit does not exist for SK-"+toString(skheadg_.sk_geometry)
				    +" MC!",v_error,m_verbose);
				return false;
			}
		}
	} else {
		switch (skheadg_.sk_geometry) {
			case 1: {
				lfallfit_sk1_data_(&watert, &max_nqisk_for_clusfit, &lfflag);
				break;
			}
			case 2: {
				lfallfit_sk2_data_(&watert, &max_nqisk_for_clusfit, &lfflag);
				break;
			}
			case 3: {
				// lfallfit_sk3_data_(&watert, &max_nqisk_for_clusfit, &lfflag); // does not exist
				Log(m_unique_name+": Error! lfallfit does not exist for SK-III!",v_error,m_verbose);
				return false;
			}
			case 4: {
				//lfallfit_sk4_data_(&watert, &max_nqisk_for_clusfit, &lfflag);
				//lfallfit_sk4_final_qe41_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				lfallfit_sk4_final_qe43_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				//lfallfit_sk4_gain_corr_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			case 5: {
				lfallfit_sk5_data_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			case 6: {
				lfallfit_sk6_data_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			/*
			case 7: {
				// FIXME can't do this, it crashes with "trabs_sk6: absfitfac is negative!"
				lfallfit_sk6_data_(&watert, &max_nqisk_for_clusfit, &flag_skip, &flag_log, &lfflag);
				break;
			}
			*/
			default: {
				Log(m_unique_name+": Error! lfallfit does not exist for SK-"+toString(skheadg_.sk_geometry)
				    +" data!",v_error,m_verbose);
				return false;
			}
		}
	}
	Log(m_unique_name+" fitting done with return value lfflag "+toString(lfflag),v_debug,m_verbose);
	
	if(lfflag==-1){
		Log(m_unique_name+" Warning! lfallfit failed with status -1: 'too many hits for lowe'!\n"
		    "Maybe an upstream nqisk cut is required? NQISK for this event: "+toString(skq_.nqisk),
		    v_warning,m_verbose);
		//return true;   // XXX i guess maybe just record what we have?
	} else if(lfflag==-2){
		Log(m_unique_name+" Warning! lfallfit failed with status -2: 'no hits'!\n",v_error,m_verbose);
		return false;    // XXX we really should not have events with no hits...!
	}
	
	Log(m_unique_name+" setting results into lowe branch",v_debug,m_verbose);
	// pass reconstructed variables from skroot_lowe_ common block (populated by lfallfit) to skroot file
	skroot_set_lowe_(&lun,                      &skroot_lowe_.bsvertex[0], &skroot_lowe_.bsresult[0],
	                 &skroot_lowe_.bsdir[0],    &skroot_lowe_.bsgood[0],   &skroot_lowe_.bsdirks,
	                 &skroot_lowe_.bseffhit[0], &skroot_lowe_.bsenergy,    &skroot_lowe_.bsn50,
	                 &skroot_lowe_.bscossun,    &skroot_lowe_.clvertex[0], &skroot_lowe_.clresult[0],
	                 &skroot_lowe_.cldir[0],    &skroot_lowe_.clgoodness,  &skroot_lowe_.cldirks,
	                 &skroot_lowe_.cleffhit[0], &skroot_lowe_.clenergy,    &skroot_lowe_.cln50,
	                 &skroot_lowe_.clcossun,    &skroot_lowe_.latmnum,     &skroot_lowe_.latmh,
	                 &skroot_lowe_.lmx24,       &skroot_lowe_.ltimediff,   &skroot_lowe_.lnsratio,
	                 &skroot_lowe_.lsdir[0],    &skroot_lowe_.spaevnum,    &skroot_lowe_.spaloglike,
	                 &skroot_lowe_.sparesq,     &skroot_lowe_.spadt,       &skroot_lowe_.spadll,
	                 &skroot_lowe_.spadlt,      &skroot_lowe_.spamuyn,     &skroot_lowe_.spamugdn,
	                 &skroot_lowe_.posmc[0],    &skroot_lowe_.dirmc[0],    &skroot_lowe_.pabsmc[0],
	                 &skroot_lowe_.energymc[0], &skroot_lowe_.darkmc,      &skroot_lowe_.islekeep,
	                 &skroot_lowe_.bspatlik,    &skroot_lowe_.clpatlik,    &skroot_lowe_.lwatert,
	                 &skroot_lowe_.lninfo,      &skroot_lowe_.linfo[0]);
	
	return true;
	
}


bool lfallfit::Finalise(){
	
	return true;
}

