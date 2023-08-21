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
	m_variables.Get("reference_watert_run", reference_watert_run);
	
	// check the upstream TreeReader name
	if(m_data->Trees.count(readerName)==0){
		Log("Failed to find TreeReader "+readerName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	// we need the associated LUN number for skroot_* functions
	lun = m_data->GetLUN(readerName);
	// and to check if this file is MC or not, so we can use the appropriate lf_allfit version
	MC = m_data->Trees.at(readerName)->GetMCFlag();
	
	// initialize water transparency table
	skrunday_();
	skwt_gain_corr_();
	
	// initialize bonsai
	m_data->BonsaiInit();
	
	if(MC && skhead_.nrunsk==999999){
		Log(m_unique_name+" using reference run "+toString(reference_watert_run)
		    +" for water transparency",v_warning,m_verbose);
	}
	
	return true;
}


bool lfallfit::Execute(){
	
	if((nread%10000)==0){
		Log(m_unique_name+" read loop "+toString(nread)+", current run "
		    +toString(skhead_.nrunsk),v_message,m_verbose);
	}
	++nread;
	
	// FIXME should this whole water transparency / bad channel thing be part of e.g. the TreeReader?
	
	// lfallfit requires a water transparency
	// we can look up a suitable value with lfwater_, but we need a run number
	// for MC, if the run number hasn't already been overridden, use a suitably representative run
	if(MC && skhead_.nrunsk==999999){
		// nrunsk of 999999 is expected for MC...
		// this gets printed every time as SKREAD overwrites nrunsk with each new event
		//Log(m_unique_name+" using reference run "+toString(reference_watert_run)
		//    +" for water transparency",v_debug,m_verbose);
		skhead_.nrunsk = reference_watert_run;
	}
	
	// retreive new water transparency when changing runs
	if(skhead_.nrunsk!=nrunsk_last){
		int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];  // defined in skdayC.h
		lfwater_(&days_to_run_start, &watert);
		Log(m_unique_name+" loaded new water transparency value "+toString(watert)
			+" for run "+toString(skhead_.nrunsk),v_debug,m_verbose);
		nrunsk_last = skhead_.nrunsk;
	}
	
	// update the bad channel list when changing runs or subruns (MC only? FIXME does this make sense?)
	if(MC){
		if(skhead_.nrunsk!=nrunsk_last || skhead_.nsubsk!=nsubsk_last){
			int ierr;
			skbadch_(&skhead_.nrunsk,&skhead_.nsubsk,&ierr);
			nrunsk_last = skhead_.nrunsk;
			nsubsk_last = skhead_.nsubsk;
			if(skhead_.nrunsk!=nrunsk_last) darklf_(&skhead_.nrunsk);
		}
	}
	
	// clear variables for low & mu fitters
	lfclear_all_();
	
	// apply lowfit
	int NHITCUT = (MC) ? 800 : 1000;  //  number of hit limit for clusfit
	
	int lfflag;
	int log_level;
	int flag_skip=0;
	
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
			lfflag is output: 0=fit success, -1=too many hits for lowe, -2: no hits
	*/
	if(MC){
		switch (skheadg_.sk_geometry) {
			case 1:  //lfallfit_sk1_mc_(&watert, &NHITCUT, &lfflag);    // does not exist
			case 2:  //lfallfit_sk2_mc_(&watert, &NHITCUT, &lfflag);    // does not exist
			case 3:  //lfallfit_sk3_mc_(&watert, &NHITCUT, &lfflag);    // does not exist
				Log(m_unique_name+": Error! lfallfit does not exist for SK-"+toString(skheadg_.sk_geometry)
				    +" MC!",v_error,m_verbose);
				return false;
			case 4: {
				// which lfallfit_sk4 would you prefer...?
				//lfallfit_sk4_mc_(&watert, &NHITCUT, &lfflag);
				//lfallfit_sk4_final_qe41_mc_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				lfallfit_sk4_final_qe43_mc_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				//lfallfit_sk4_gain_corr_mc_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				break;
			}
			case 5: {
				lfallfit_sk5_mc_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				break;
			}
			case 6: {
				lfallfit_sk6_mc_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
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
				lfallfit_sk1_data_(&watert, &NHITCUT, &lfflag);
				break;
			}
			case 2: {
				lfallfit_sk2_data_(&watert, &NHITCUT, &lfflag);
				break;
			}
			case 3: {
				// lfallfit_sk3_data_(&watert, &NHITCUT, &lfflag); // does not exist
				Log(m_unique_name+": Error! lfallfit does not exist for SK-III!",v_error,m_verbose);
				return false;
			}
			case 4: {
				//lfallfit_sk4_data_(&watert, &NHITCUT, &lfflag);
				//lfallfit_sk4_final_qe41_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				lfallfit_sk4_final_qe43_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				//lfallfit_sk4_gain_corr_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				break;
			}
			case 5: {
				lfallfit_sk5_data_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				break;
			}
			case 6: {
				lfallfit_sk6_data_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
				break;
			}
			default: {
				Log(m_unique_name+": Error! lfallfit does not exist for SK-"+toString(skheadg_.sk_geometry)
				    +" data!",v_error,m_verbose);
				return false;
			}
		}
	}
	
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

