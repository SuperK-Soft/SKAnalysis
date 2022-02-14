/* vim:set noexpandtab tabstop=4 wrap */
#include "lf_allfit.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <string>
#include <vector>
#include <iostream>

// declarations and #includes for SK fortran routines
#include "fortran_routines.h"

lf_allfit::lf_allfit():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool lf_allfit::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	m_variables.Get("fname_in",fname_in);
	m_variables.Get("fname_out",fname_out);
	
	lun = 10;  // TODO we should track these in ToolAnalysis to ensure uniqueness
	
	// set up branches we don't need to read on input or write out
	std::vector<std::string> in_branches_to_skip
		{"SPACERS", "QBEESTATUS", "DBSTATUS", "MISMATCHEDHITS", "ATMPD", "UPMU"};
	std::vector<std::string> out_branches_to_skip
		{"TQLIST","ODTQLIST","SOFTWARETRG","EVENTTRAILER","HWTRGLIST",
		"PEDESTALS","EVENTHEADER","GPSLIST","PREVT0","SLE","T2KGPSLIST"};
	
	// open input file
	skroot_open_(&lun, fname_out.c_str(), fname_out.size());
	skroot_set_input_file_(&lun, fname_in.c_str(), fname_in.size());
	
	// disable unused input branches
	int io_dir = 0; // 0 for input branch, 1 for output branch
	for(auto&& abranch : in_branches_to_skip){
		skroot_zero_branch_(&lun, &io_dir, abranch.c_str(), abranch.size());
	}
	// disable output branches...
	io_dir = 1;
	for(auto&& abranch : out_branches_to_skip){
		skroot_zero_branch_(&lun, &io_dir, abranch.c_str(), abranch.size());
	}
	
	// initialize the tree manager
	skroot_init_(&lun);
	// skroot_init.F is mostly just a wrapper around skroot_initialize_,
	// but also sets the variable `SK_FILE_FORMAT = 1` in common block /SKHEADF/
	
	// initialize data structure (zbs)
	kzinit_();
	std::string options = "31,30,26,25";
	skoptn_(const_cast<char*>(options.c_str()), options.size());
	int badopt = 23;
	skbadopt_(&badopt);
	// need to set skgeometry in skheadg common block
	skheadg_.sk_geometry = 4;
	geoset_();
	
	// initialize water transparency table
	skrunday_();
	skwt_gain_corr_();
	
	// MAXPM is a #defined constant in tqrealroot.h, but we need to pass it by reference, so need a copy.
	int MAXPM_var = MAXPM;
	// XYZPM is a float[MAXPM][3] defined in the geompmt common block in geopmtC.h
	float* xyzpm = &geopmt_.xyzpm[0][0];
	// initialize bonsai
	cfbsinit_(&MAXPM_var, xyzpm);
	
	return true;
}


bool lf_allfit::Execute(){
	
	if((nread%10000)==0) std::cout<<"nrunsk/nread = "<<skhead_.nrunsk<<"/"<<nread<<std::endl;
	++nread;
	
	// read input
	int ierr;
	skcrawread_(&lun, &ierr);
	if(ierr==1){
		std::cerr<<"read error"<<std::endl;
		m_data->vars.Set("StopLoop",1);
		return true;
	} else if(ierr==2){
		m_data->vars.Set("StopLoop",1);
		return true;
	} else if(ierr!=0) {
		//std::cout<<"possibly some recoverable error? continuing"<<std::endl;
		// this happens a lot...
		return true;
	}
	// same now for SKREAD...
	int neglun = -lun;  // is used to indicate to SKREAD that it's a root file not zbs
	skcread_(&neglun, &ierr);
	if(ierr==1){
		std::cerr<<"read error"<<std::endl;
		m_data->vars.Set("StopLoop",1);
		return true;
	} else if(ierr==2){
		std::cerr<<"end of file"<<std::endl;
		m_data->vars.Set("StopLoop",1);
		return true;
	} else if(ierr!=0) {
		//std::cout<<"possibly some recoverable error? continuing"<<std::endl;
		// this happens a lot (or, one of them does)
		return true;
	}
	
	// once per run water transparency
	if(skhead_.nrunsk!=nrunsk_last){
		int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];  // defined in skdayC.h
		// function can be found at /usr/local/skofl/src/wtlib.obsolete/  to calculate watert
		lfwater_(&days_to_run_start, &watert);
		std::cout<<"nrunsk/watert = "<<skhead_.nrunsk<<"/"<<watert<<std::endl;
		nrunsk_last = skhead_.nrunsk;
	}
	
	// clear variables for low & mu fitters
	lfclear_all_();
	
	int lfflag;
	int log_level;
	int NHITCUT = 1000;  //  number of hit limit for clusfit (changed for relic analysis)
	int flag_skip=0;
	// apply lowfit
	/*
	C         call lfallfit_sk4_data(watert, NHITCUT, lfflag)
	C         call lfallfit_sk4_gain_corr(watert, NHITCUT, 0, log_level, lfflag)
	*/
	lfallfit_sk4_final_qe43_(&watert, &NHITCUT, &flag_skip, &log_level, &lfflag);
	
	/*
	c***         read all branches into memory for output
	c***         an error occur when filling output tree if this is not called.
	c***         in skrawread(), only head, tq and pedestal branches are read.
	c***         to copy entire event, call this function.
	c***         however, if all irrelevant branches are skipped, then it's ok
	c***         w/o calling this function.
	c            call skroot_get_entry(lun)
	*/
	
	// store LOWE branch to skroot file
	// at this point we're back to calling fortran interface C functions
	// and can drop the trailing underscore
	// (though we still need to make everything pointers)
	// of course we're passing with fortan variables from common blocks...
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
	
	// remove hits outside 1.3 microsec
	delete_outside_hits_();
	
	// store header & TQ info.
	// skroot_set_tree.F is just a wrapper around another couple of calls like skroot_set_lowe,
	// that simply pull variables from fortran common blocks and pass them to the TreeManager
	skroot_set_tree_(&lun);
	
	// output root file - another fortran interface function we can call directly.
	skroot_fill_tree_(&lun);
	
	return true;
}


bool lf_allfit::Finalise(){
	
	// close input skroot files, delete the TTreeManager
	skroot_close_(&lun);
	// delete the SuperManager
	skroot_end_();
	
	// terminate bonsai
	cfbsexit_();
	
	return true;
}

