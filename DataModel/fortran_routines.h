// (c-versions of) headers that define fortran common blocks
// these are offficial in $SKOFL_ROOT/inc/
#include "skheadC.h"
#include "skparmC.h"
#include "sktqC.h"
#include "skbadcC.h"
#include "geopmtC.h"
#include "skruninfC.h"  // commented out in original lowfit_sk4
//#include "geotnkC.h"  // this conflicts with SK2p2MeV.h because someone reused the #define'd macro names as variables (e.g. `const Float_t RINTK`, when geotnkC.h does `#define RINTK`)

// these are from $SKOFL_ROOT/inc/lowe and have no official C version
// the following two were done via automatic conversion done with fh2h
#include "skdayC.h"
#include "skwtC.h"
// this last one was manually edited to reflect 'EQUIVALENCE' use (use with caution)
#include "skroot_loweC.h"

#include "geopmaC.h"
#include "skwaterlenC.h"
#include "odmaskflagC.h"
#include "skdbstatC.h"
#include "skqbstatC.h"
#include "skspacerC.h"
#include "skgpsC.h"
#include "skprevt0C.h"
#include "sktrighitC.h"
#include "sktrgC.h"
#include "skvetoC.h"
#include "softtrg_listC.h"
#include "vcvrtxC.h"
#include "vcworkC.h"

// header for skroot_* functions. These are actually C functions.
#include "fortran_interface.h"

// we need to declare all external fortran routines used
// the following were fine without any additional libraries being added to the makefile
extern "C" void skroot_init_(int*);
extern "C" void kzinit_();
extern "C" void skoptn_(char*, int);
extern "C" void skbadopt_(int*);
extern "C" void skbadch_(int*, int*, int*);
extern "C" void geoset_();
extern "C" void delete_outside_hits_();
extern "C" void skcrawread_(int*, int*);
extern "C" void skcread_(int*, int*);
extern "C" void skroot_set_tree_(int*);
extern "C" void skroot_get_entry_(int*);

// additional functions from lfallfit
extern "C" void skatmmap2_(int*,int*,int*,int*,int*,int*);
extern "C" void lfneihit_(int*,int*,int*);
extern "C" void lfqhit_(int*,float*); 
extern "C" void lfflasher_(int*,float*,int*);
extern "C" float wallsk_(float*);
extern "C" void lfdir2_(float*,float*,float*);
extern "C" void lfneweff_sk3_final_(float*, float*, Float_t*);
extern "C" void lfneweff_sk4_final_qe43_(float*,float*,float*,float*,int*,int*);
extern "C" void lfeffwt_sk3_(float*, float*,float*);
extern "C" float enelf_sk4_(float*,int*);
extern "C" void lfdir4_(float*,float*,float*,float*);
extern "C" void lfdirks_(float*,float*,float*);
extern "C" float patliklf_sk4_(float*,float*,float*,int*);
extern "C" float effwallf_(int*,float*,float*,float*);
extern "C" void slsundir3_(int*,int*,float*,float*,float*,float*);
extern "C" void lfnhit_tisk_(int*,float*,int*,int*,int*);
extern "C" void lfnhit2_(int*,float*,int*,int*,int*);
extern "C" void lfnhit_(int*,float*,int*,int*);
extern "C" float lfcal_r02_(int*,int*,float*);
extern "C" void lfariadne_(float*,int*,float*,float*,float*,float*,int*,float*);

// from $ATMPD_ROOT/src/programs/TreeBuilder/examples/fort_fopen.F
extern "C" void fort_fopen_(int*, const char*, char*, int* ,int);

// read runinf (y/m/d, start time, end time, etc)
extern "C" void runinfsk_();

//
extern "C" void softtrg_set_cond_(int*, int*, int*, int*, int*);
extern "C" void get_sub_triggers_(int*, int*, int*, int*);
extern "C" void get_sub_triggers_add_(int*, int*, int*, int*, int*);
extern "C" void get_sub_triggers2_(int*, int*, int*, int*, int*);
extern "C" void set_timing_gate_(int*);
extern "C" void set_timing_gate_nsec_(float*);
extern "C" int make_swtrg_list(int*,int*);
extern "C" void softtrg_get_cond_(int*, int*, int*, int*, int*);
// for these two, see $SKOFL_ROOT/src/softtrg/softtrg_ofldata.c
extern "C" int softtrg_inittrgtbl_(int*, int*, int*,int*); // probably call this
extern "C" int softtrg_inittrgtbl(int , int, int, int);    // which internally calls this

// the following are provided by libwtlib_5.1.a
extern "C" void skrunday_();
extern "C" void skwt_();
extern "C" void skwt_gain_corr_();
extern "C" void lfwater_(int*, float*);
// skday_data_, common block

// the following are provided by libbonsai_3.3.a
extern "C" void cfbsinit_(int*, float*);
extern "C" void cfbsexit_();
extern "C" int bonsaifit_(float*,float*,float*,int*,int*,int*,float*,float*);

// the following are provided by libsklowe_7.0.a
extern "C" void lfclear_all_();
extern "C" void lfallfit_sk1_data_(float* watert, int*nhitcut, int*lfflag);
//extern "C" void lfallfit_sk1_mc_(float* watert, int*nhitcut, int*lfflag);    // does not exist
extern "C" void lfallfit_sk2_data_(float* watert, int*nhitcut, int*lfflag);
//extern "C" void lfallfit_sk2_mc_(float* watert, int*nhitcut, int*lfflag);    // does not exist
// extern "C" void lfallfit_sk3_data_(float* watert, int*nhitcut, int*lfflag); // does not exist
// extern "C" void lfallfit_sk3_mc_(float* watert, int*nhitcut, int*lfflag); // does not exist
extern "C" void lfallfit_sk4_data_(float* watert, int*nhitcut, int*lfflag);
extern "C" void lfallfit_sk4_mc_(float* watert, int*nhitcut, int*lfflag);
extern "C" void lfallfit_sk5_mc_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk5_data_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk6_mc_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk6_data_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
// for some reason SK-IV has a bunch of extra versions...
extern "C" void lfallfit_sk4_final_qe41_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk4_final_qe41_mc_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk4_final_qe43_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk4_final_qe43_mc_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk4_gain_corr_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);
extern "C" void lfallfit_sk4_gain_corr_mc_(float* watert, int* nhitcut, int* flag_skip, int* flag_log, int* lfflag);

// skroot_lowe_ common block

extern "C" void rluxgo_(int*, int*, int*, int*);
extern "C" float rlu_();  // src/monlib/rlu.F

extern "C" void slmcmklow_(int*, int*, int*, int*, float*);  // $SKOFL_ROOT/lowe/sollib/slmcmklow.F
extern "C" void slredtimev_(int*, int*, int*, int*, int*, float*, float*, int*); // $SKOFL_ROOT/lowe/sollib/slredtimev.f

extern "C" int lfbadrun_(int*, int*); // $SKOFL_ROOT/lowe/sklowe/lfbadrun.F

extern "C" void darklf_(int*);  // $SKOFL_ROOT/lowe/sklowe/darklf.F

// after that there were many undefined references, e.g. `sortzv_`, `hf1_`, `hf2_`...
// after some trial and error these are resolved, but i lost track of which provided what.
// cernlibs in particular resolved a lot of repeated undefined issues, they may be the main culprit.

// muon reconstruction
extern "C" void fix_maxqisk_();
extern "C" void lfmufit_sk4_();
extern "C" void newmufit_(float (*)[3], float (*)[3], float*);
extern "C" void makededx_(float (*)[4], float (*)[3], int (*)[11146], float (*)[11146], float (*)[11146], float (*)[11146][3], int*, float (*)[200]);
extern "C" void makededx_intg_(float (*)[4], float (*)[3], float*, int (*)[11146], float (*)[11146], float (*)[11146], float (*)[11146][3], int*, int*, float (*)[200], int (*)[334380], int*);
extern "C" void mfmuselect_(float (*)[3], float (*)[3], float*, int*);
extern "C" void mffastfast_(float (*)[3], float (*)[3], int*);
extern "C" void muboy_zbs_(int*, int*, float (*)[4], float (*)[3], float*, float*, int*, float (*)[36], int*);

// vectgen
extern "C" void spectrum_(int*,int*,int*,float*,int*,float*,float*,int*);

// SK I/O
extern "C" {
//    void kzinit_();
//    void geoset_();
    void zbsinit_();
    void kzwrit_(int&);
    void kzeclr_();
    // void set_rflist_(int&, const char*, const char*, const char*, const char*,
                     // const char*, const char*, const char*, const char*, const char*,
                     // long, long,  long,  long,  long,  long,  long,  long,  long  );
    // void skopenf_(int&, int&, const char*, int&, long);
    void set_rflist_(int*, const char*, const char*, const char*, const char*,
                     const char*, const char*, const char*, const char*, const char*,
                     int, int, int, int, int, int, int, int, int);
    void skopenf_(int*, int*, const char*, int*, int*);
//    void skoptn_(const char*, int);
//    void skbadopt_(int*);
//    void skbadch_(int*, int*, int*);
    int  skread_(int*);
    int  skrawread_(int*);
    void skclosef_(int*);
//    void skroot_init_(int*);
    void cclose_(int*); // $SKOFL_ROOT/src/iolib/cclose.f
}

// data control
extern "C" {
    void  nerdnebk_(float*);
    void  skgetv_();
    void  apflscndprt_();
    void  trginfo_(float*);
    void  aprstbnk_(int*);
    void  odpc_2nd_s_(int*);
    void  inpmt_(float*, int&);
    float wallsk_(float*);
}

// BONSAI
extern "C" {
// from SKLibs - maybe incorrect signature?
//    void bonsai_ini_(int*);
//    void bonsai_fit_(int*, float*, float*, float*, int*, int*, float*, float*, float*,
//                     float*, float*, float*, float*, float*, float*);
    // from SK2p2MeV
    void bonsai_ini_();
    void bonsai_fit_(float *t0, float *tisksend, float *qisksend, int *cabsend, int *ntisksend,
                     float *tenergy, float *tvx, float *tvy, float *tvz, float *tbsgood);
    void bonsai_end_();
     void bonsai_energy_(float *,float*,float*,float*,int*,int*,float*);
    // also from SK2p2MeV - but commented out
    //void bonsai_combined_ini_();
    //void bonsai_fit_combined_(float *t0, float *tiskp, float *qiskp, int *cabp, int *ntiskp,
    //                          float *tisksend, float *qisksend, int *cabsend, int *ntisksend,
    //                          float *tvx, float *tvy, float *tvz, float *tbsgood);
    //void bonsai_combined_end_();
}

// stopmu fit
extern "C" {
    void stmfit_(float*, float*, float&, float&);
    void sparisep_(int&, int&, int&, int&);
    void pfdodirfit_(int&);
    void sppang_(int&, float&, float&);
    void spfinalsep_();
//	float pttruewaterlen_(float&);
}
