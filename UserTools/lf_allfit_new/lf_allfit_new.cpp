/* vim:set noexpandtab tabstop=4 wrap */
#include "lf_allfit_new.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <string>
#include <vector>
#include <iostream>

// declarations and #includes for SK fortran routines
#include "fortran_routines.h"

lf_allfit_new::lf_allfit_new():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool lf_allfit_new::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("readerName",readerName);          // name given to the TreeReader used for file handling
	m_variables.Get("writeout",writeout);              // whether to write out to a new file or not (just reconstruction vs inline)
	
	
	// use the readerName to find the LUN associated with this file
	std::map<std::string,int> lunlist;
	m_data->CStore.Get("LUNList",lunlist);
	if(lunlist.count(readerName)==0){
		Log(toolName+" error! No LUN associated with readerName "+readerName,v_error,verbosity);
		return false;
	}
	lun = lunlist.at(readerName);
	
	// determine if this file is MC or not, so we can use the appropriate lf_allfit version
	TTree* t = skroot_get_tree(&lun);
	if(t) MC = (t->FindBranch("MC")!=nullptr);
	
	// initialize water transparency table
	skrunday_();
	skwt_gain_corr_();
	
	// initialize bonsai
	int MAXPM_var = MAXPM;
	float* xyzpm = &geopmt_.xyzpm[0][0];
	cfbsinit_(&MAXPM_var, xyzpm);
	
	return true;
}


bool lf_allfit_new::Execute(){
	
	if((nread%10000)==0){
		Log(toolName+" read loop "+toString(nread)+", current run "+toString(skhead_.nrunsk),v_message,verbosity);
	}
	++nread;
	
	if(MC && skhead_.nrunsk==999999){
		Log(toolName+" warning: no run number!!",v_warning,verbosity);
		skhead_.nrunsk = 75000;
	}
	
	// once per run update the water transparency
	if(skhead_.nrunsk!=nrunsk_last){
		int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];  // defined in skdayC.h
		lfwater_(&days_to_run_start, &watert);
		Log(toolName+" loaded new water transparency value "+toString(watert)
			+" for run "+toString(skhead_.nrunsk),v_debug,verbosity);
		nrunsk_last = skhead_.nrunsk;
	}
	
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
	int NHITCUT = (MC) ? 800 : 1000;  //  number of hit limit for clusfit (changed for relic analysis)
	
	int lfflag;
	int log_level;
	int flag_skip=0;
	if(MC){
		lfallfit_sk4_mc_(&watert, &NHITCUT, &lfflag);
		//lfallfit_sk4_gain_corr_mc_(&watert, &NHITCUT, 0, &log_level, &lfflag);
	} else {
		//lfallfit_sk4_data_(&watert, &NHITCUT, &lfflag);
		//lfallfit_sk4_gain_corr(&watert, &NHITCUT, 0, &log_level, &lfflag);
		lfallfit_sk4_final_qe43_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
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
	
	// if the writeout variable is set then pass reconstructed information
	// from the skroot_lowe_ common block to the output skroot file
	if(writeout){
		// remove hits outside 1.3 microsec
		delete_outside_hits_();
		
		// store header & TQ info.
		// skroot_set_tree.F is just a wrapper around another couple of calls like skroot_set_lowe,
		// that simply pull variables from fortran common blocks and pass them to the TreeManager
		skroot_set_tree_(&lun);
		
		// invokes TTree::Fill. Only use it in SKROOT mode WRITE or COPY!
		skroot_fill_tree_(&lun);
	}
	
	return true;
	
}


bool lf_allfit_new::Finalise(){
	
	// terminate bonsai
	cfbsexit_();
	
	return true;
}

