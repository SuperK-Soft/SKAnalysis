/*==============================================================================
C  
C   lowfit_sk4.F 
C    - apply lowfit to rfm data for relic analysis
C    - /usr/local/sklib_g77/skofl/lowe/sklowe/lfallfit_sk4_data.F
C
C    [output]
C     (arg1) .root lo file  
C
C    [input]
C     (arg2-) .root rfm file
C
==============================================================================*/
#include "lowfit_sk4_stripped.h"
#include <stdint.h>
#include "skheadC.h"
#include "skparmC.h"
#include "sktqC.h"
#include "skbadcC.h"
#include "geopmtC.h"
//#include "skruninfC.h"
// following are in $SKOFL_ROOT/inc/lowe, no C conversion
// automatic conversion done with fh2h, headers in local inc:
#include "skdayC.h"
#include "skwtC.h"
// following had a lot of manual work as well to resolve 'EQUIVALENCE' use
#include "skroot_loweC.h"  // in inc, manual conversion...
// other errors seen during automatic conversion, so use with caution

//#include "skroot.h"
#include "fortran_interface.h"

#include <string>
#include <vector>
#include <iostream>

// these seem ok
extern "C" void skroot_init_(int*);
extern "C" void kzinit_();
extern "C" void skoptn_(char*, int);
extern "C" void skbadopt_(int*);
extern "C" void geoset_();
extern "C" void delete_outside_hits_();
extern "C" void skcrawread_(int*, int*);
extern "C" void skcread_(int*, int*);
extern "C" void skroot_set_tree_(int*);

// the following are provided by libwtlib_5.1.a
extern "C" void skrunday_();
extern "C" void skwt_gain_corr_();
extern "C" void lfwater_(int*, float*);
// skday_data_, common block
// but we have a new undefined ref in `skwt_gain_corr_` to `sortzv_`

// the following are provided by libbonsai_3.3.a
extern "C" void cfbsinit_(int*, float*);
extern "C" void cfbsexit_();

// these are provided by libsklowe_7.0.a
extern "C" void lfclear_all_();
extern "C" void lfallfit_sk4_final_qe43_(float*, int*, int*, int*, int*);
// skroot_lowe_ common block
// but we have maaaaany new undefined references.
// had to manually add a bunch of libs to the makefile, cernlibs in particular, possibly others

int do_lowfit(int argc, const char* argv[]){

	// need to set skgeometry in skheadg common block
	skheadg_.sk_geometry = 4;
	
	if(argc<3){
		std::cout<<"usage: "<<argv[0] <<" filename_out filename_in..."<<std::endl;
		return 0;
	}
	std::string fname_out = argv[1];
	std::string fname_in = argv[2];
	
	// ok so actually all this initial SKROOT stuff is really just obfuscation
	// around standard ROOT file IO performed by the TreeManager.
	// note: "skrooth.h" provides #defines that make some fortran interface methods
	// accessible without trailing underscores, but not all of them.
	// We'll reproduce it here to minimize changes, but there are easier ways to handle ROOT files.
	
	std::vector<std::string> in_branches_to_skip
		{"SPACERS", "QBEESTATUS", "DBSTATUS", "MISMATCHEDHITS", "ATMPD", "UPMU"};
	std::vector<std::string> out_branches_to_skip
		{"TQLIST","ODTQLIST","SOFTWARETRG","EVENTTRAILER","HWTRGLIST",
		"PEDESTALS","EVENTHEADER","GPSLIST","PREVT0","SLE","T2KGPSLIST"};
	
	// open input file
	int lun = 10;                // logic unit number
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
	// skroot_init.F sets the variable `SK_FILE_FORMAT = 1`
	// in common block /SKHEADF/, so perhaps this is important.
	
	// c***  read runinf (y/m/d, start time, end time, etc)
	// c     call runinfsk
	
	kzinit_();  // initialize data structure(zbs)
	std::string options = "31,30,26,25";
	skoptn_(const_cast<char*>(options.c_str()), options.size());
	int badopt = 23;
	skbadopt_(&badopt);
	geoset_();
	
	/*
	c***  initialize water transparency table
	call skrunday
C	call skwt
	call skwt_gain_corr
	*/
	
	skrunday_();
	skwt_gain_corr_();
	
	/*
	c***  initialize bonsai 
	call cfbsinit(MAXPM,xyzpm)
	*/
	// MAXPM is a #defined constants in tqrealroot.h
	int MAXPM_var = MAXPM;
	// XYZPM is a 2D array of floats defined in the GEOPMT common block in geopmt.h OR geopmtC.h
	float* xyzpm = &geopmt_.xyzpm[0][0];  /// xyzpm is a float[MAXPM][3]
	cfbsinit_(&MAXPM_var, xyzpm);
	
	// main loop 
	int nread=0;
	int  nrunsk_last=0;
	float watert;
	
	// event loop
	while(true) {
		++nread;
		if((nread%10000)==0) std::cout<<"nrunsk/nread = "<<skhead_.nrunsk<<"/"<<nread<<std::endl;
		
		// read input
		/*
		call SKRAWREAD(lun, *1002, *1001, *10, *10)   !  read data from bank to commom
		call SKREAD(-lun, *1002, *1001, *10, *10)     !  read event when lun is positive
		*/
		int rawret;
		skcrawread_(&lun, &rawret);
		if(rawret==1){
			std::cerr<<"read error"<<std::endl;
			break;
		} else if(rawret==2){
			std::cerr<<"end of file"<<std::endl;
			break;
		} else if(rawret!=0) {
			//std::cout<<"possibly some recoverable error? continuing"<<std::endl;
			// this happens a lot...
			continue;
		}
		// same now for SKREAD...
		int neglun = -lun;  // is used to indicate to SKREAD that it's a root file not zbs
		skcread_(&neglun, &rawret);
		if(rawret==1){
			std::cerr<<"read error"<<std::endl;
			break;
		} else if(rawret==2){
			std::cerr<<"end of file"<<std::endl;
			break;
		} else if(rawret!=0) {
			//std::cout<<"possibly some recoverable error? continuing"<<std::endl;
			// this happens a lot (or, one of them does)
			continue;
		}
		
		// once per run water transparency
		if(skhead_.nrunsk!=nrunsk_last){
			int days_to_run_start = skday_data_.relapse[skhead_.nrunsk]; // defined in skdayC.h
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
		          call lfallfit_sk4_final_qe43(watert, NHITCUT, 0, log_level, lfflag)
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
		skroot_set_lowe_(&lun, &skroot_lowe_.bsvertex[0], &skroot_lowe_.bsresult[0],
		                &skroot_lowe_.bsdir[0], &skroot_lowe_.bsgood[0], &skroot_lowe_.bsdirks,
		                &skroot_lowe_.bseffhit[0], &skroot_lowe_.bsenergy, &skroot_lowe_.bsn50,
		                &skroot_lowe_.bscossun, &skroot_lowe_.clvertex[0], &skroot_lowe_.clresult[0],
		                &skroot_lowe_.cldir[0], &skroot_lowe_.clgoodness, &skroot_lowe_.cldirks,
		                &skroot_lowe_.cleffhit[0], &skroot_lowe_.clenergy, &skroot_lowe_.cln50,
		                &skroot_lowe_.clcossun, &skroot_lowe_.latmnum, &skroot_lowe_.latmh,
		                &skroot_lowe_.lmx24, &skroot_lowe_.ltimediff, &skroot_lowe_.lnsratio,
		                &skroot_lowe_.lsdir[0], &skroot_lowe_.spaevnum, &skroot_lowe_.spaloglike,
		                &skroot_lowe_.sparesq, &skroot_lowe_.spadt, &skroot_lowe_.spadll,
		                &skroot_lowe_.spadlt, &skroot_lowe_.spamuyn, &skroot_lowe_.spamugdn,
		                &skroot_lowe_.posmc[0], &skroot_lowe_.dirmc[0], &skroot_lowe_.pabsmc[0],
		                &skroot_lowe_.energymc[0], &skroot_lowe_.darkmc, &skroot_lowe_.islekeep,
		                &skroot_lowe_.bspatlik, &skroot_lowe_.clpatlik, &skroot_lowe_.lwatert,
		                &skroot_lowe_.lninfo, &skroot_lowe_.linfo[0]);
		
		// remove hits outside 1.3 microsec
		delete_outside_hits_();
		
		// store header & TQ info.
		// skroot_set_tree.F is just a wrapper around another couple of calls like skroot_set_lowe,
		// that simply pull variables from fortran common blocks and pass them to the TreeManager
		skroot_set_tree_(&lun);
		
		// output root file - another fortran interface function we can call directly.
		skroot_fill_tree_(&lun);
		
	}
	
	// close input skroot files, delete the TTreeManager
	skroot_close_(&lun);
	// delete the SuperManager
	skroot_end_();
	
	// terminate bonsai
	cfbsexit_();
	
	return 0;
}
