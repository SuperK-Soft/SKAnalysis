/*******************************************
*
* @file SKLibs.hh
*
* @brief Defines SKOFL and ATMPD functions
* to use in NTag.
*
********************************************/

#ifndef SKLIBS_HH
#define SKLIBS_HH 1

// SK I/O
extern "C" {
    void kzinit_();
    void geoset_();
//XXX    void set_rflist_(int*, const char*, const char*, const char*, const char*,
//                     const char*, const char*, const char*, const char*, const char*,
//                     int, int, int, int, int, int, int, int, int);
//XXX    void skopenf_(int*, int*, const char*, int*);
//XXX    void skoptn_(const char*, int);
    void skbadopt_(int*);
    void skbadch_(int*, int*, int*);
    int  skread_(int*);
    void skclosef_(int*);
    void skroot_init_(int*);
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
//XXX    void bonsai_ini_(int*);
//XXX    void bonsai_fit_(int*, float*, float*, float*, int*, int*, float*, float*, float*,
//                     float*, float*, float*, float*, float*, float*);
    void bonsai_end_();
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

// right mess here; see ReadMCInfo/atmpd_skvectC.h,
// but it's used by multiple NTag Tools so move it here.
extern struct atmpd_vcwork_common {
  float  pos[3];
  int    ip[50];
  float  pin[50][3];
  float  pabs[50];
  int    nvect;
} skvect_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct skvect_common *skvect;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct skvect_common *skvect = &skvect_;
#endif

#endif
