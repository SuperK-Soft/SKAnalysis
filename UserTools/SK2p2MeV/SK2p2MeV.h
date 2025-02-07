/*****************************************************
 * SK2p2MeV.h                                        *
 *                                                   *
 * Author : Haibing ZHANG                            *
 *          <zhanghb02@gmail.com>                    *
 * Date   : 08/02/10                                 *
 * Comment: Modified from old class SK2p2MeV         *
 *****************************************************/

#ifndef __SK2P2MEV_H__
#define __SK2P2MEV_H__

///////////////////////////////////////////////////////
//                                                   //
// This is the base class, providing basic functions //
// for 2.2 MeV analysis.                             //
//                                                   //
///////////////////////////////////////////////////////

#include <vector>

#include "MTreeReader.h"
#include "SkrootHeaders.h"
#include "fortran_routines.h"

// use anonymous namespace to keep these local to SK2p2MeV
namespace {
    // Relevant constants
    const Float_t   RINTK = 1690.;     // Radius of inner tank
    const Float_t   ZPINTK = 1810.;    // Half height of inner tank
    const Float_t   C_WATER = 21.5833; // speed of light in water (cm/ns), see skofl/inc/lowe/lflight.par
    const Int_t     MAXNP   = 500;     // Maximum N10 peaks in 535 us data
    const Int_t     MAXPMT   = 11146;     // Maximum number of PMTs
    const Int_t     MAXN10  = 200;
    const Int_t     N200CUT = 100;
    const Float_t   T0TH = -10;
    const Int_t     MAXHITS = 100000;
    const Float_t   cut_window = 12000; //continuous dark noise cut
    const Float_t   cut_window2 = 6000; //continuous dark noise cut
    //// CHANGE!!!!!!!!!
    //const Float_t cut_window2 = 0; //continuous dark noise cut
    //const Float_t cut_window = 0; //continuous dark noise cut
    const Float_t   primary_window = 1000; //safe window for primary
    const Float_t   she_tmax = 35.e3;
    
    // Experimental
    const Int_t     MAXCOMB = 10000;  // Maximum combinations of any two hits
    const Int_t     MAXTOT  = 10000;

}

// 2.2 MeV search results
typedef struct
{
    Int_t    nhits;
    Int_t    np;          // Number of N10 peaks
    Int_t    N200M;       // Maximum N200 in 535 us
    Float_t  T200M;       // timing of N200M peak, sometimes not calculated. used in daughter function
    Int_t    N10[MAXNP];
    Int_t    N200[MAXNP];
    Int_t    N10d[MAXNP];
    Int_t    Nc[MAXNP];
    Int_t    Nback[MAXNP];
    Int_t    N300[MAXNP];
    Float_t  trms[MAXNP];
    Float_t  trmsdiff[MAXNP];
    Float_t  phirms[MAXNP];
    Float_t  bsdirks[MAXNP];
    Float_t  thetam[MAXNP];
    Float_t  dt[MAXNP];    // Time diff from delayed signal to prompt signal
    Float_t  dtn[MAXNP];   // Time diff from delayed signal to prompt signal
    Float_t  nvx[MAXNP], nvy[MAXNP], nvz[MAXNP];
    Float_t  cvertexx[MAXNP], cvertexy[MAXNP], cvertexz[MAXNP];
    Float_t  pvertexx[MAXNP], pvertexy[MAXNP], pvertexz[MAXNP];
    Float_t  nvertexx[MAXNP], nvertexy[MAXNP], nvertexz[MAXNP];
    Float_t goodness_combined[MAXNP], goodness_prompt[MAXNP], goodness_neutron[MAXNP];
    Float_t goodness_prompt1[MAXNP], goodness_neutron1[MAXNP], goodness_window[MAXNP];
    Float_t  fpdist[MAXNP];   // force fit to initial
    Float_t  bpdist[MAXNP];   // bonsai fit to initial
    Float_t  fwall[MAXNP];    // force fit to wall
    Float_t  bwall[MAXNP];    // bonsai fit to wall
    Float_t  bse[MAXNP];      // bonsai energy
    Float_t  pvx[MAXNP], pvy[MAXNP], pvz[MAXNP];   // initial vertex
    Int_t    tindex[MAXNP];
    Int_t    mctruth_neutron[MAXNP];
    Int_t    n40index[MAXNP];
    Float_t  mintrms_6[MAXNP], mintrms_3[MAXNP];
    //*********new added by Yang Zhang*****
    Float_t  Q10[MAXNP];   //total p.e. in 10 ns
    Float_t  Qmean[MAXNP];   //mean p.e. in 10 ns
    Float_t  Qrms[MAXNP];   //p.e. rms in 10 ns
    Float_t  thetarms[MAXNP];
    Int_t    NLowtheta[MAXNP];
    Int_t    cable[MAXCOMB]; //use the same number
    //************************************
    Int_t    Neff[MAXNP];  // obselete, equal to N10-Nback
    Float_t  ratio[MAXNP]; // obselete, weight ratio
    Int_t    Nc1[MAXNP];   // obselete, previous definition
    Int_t    NhighQ[MAXNP];// obselete
    Int_t    NlowQ[MAXNP];// obselete
    
    // # of PMTs w/ low hit probability
    Int_t    Nlow1[MAXNP]; // 0.50
    Int_t    Nlow2[MAXNP]; // 0.55
    Int_t    Nlow3[MAXNP]; // 0.60
    Int_t    Nlow4[MAXNP]; // 0.65
    Int_t    Nlow5[MAXNP]; // 0.70
    Int_t    Nlow6[MAXNP]; // 0.75
    Int_t    Nlow7[MAXNP]; // 0.80
    Int_t    Nlow8[MAXNP]; // 0.85
    Int_t    Nlow9[MAXNP]; // 0.90
    
    // Truth level signal fraction
    Float_t sig_frac_peak[MAXNP];
    Float_t sig_frac_tot[MAXNP];
    Int_t nsignal;
    
    // Experimental
    Int_t    ncomb;            // Number of combinations of any two hits
    Int_t    ncomb3;           // Number of hit vectors
    Float_t  dalpha[MAXCOMB];  // Angles between any two hits.
    Int_t    ncomb1;           // after cluster cut
    Float_t  dalpha1[MAXCOMB]; // after cluster cut
    Int_t    ncomb2;           // after backward hits search
    Float_t  dalpha2[MAXCOMB]; // after backward hits search
    Float_t  theta[MAXCOMB];   // theta distritution
    
    Float_t EneThre;
    
    // Store SHE hits for dark noise cut
    std::vector<float> shet;
    std::vector<int> shecab;
    std::vector<float> aftt;
    std::vector<float> aftq;
    std::vector<int> aftcab;
    
    void Clear()
        {
            np   = 0;
            N200M = 0;
            T200M = -9999.;
            for (Int_t i=0; i<MAXNP; i++) {
                N10[i]   = 0;
                Nc[i]    = 0;
                Nback[i] = 0;
                N300[i]  = 0;
                //*********new added by Yang Zhang*****
                Q10[i] = 0;
                Qmean[i] = 0;
                Qrms[i] = 0;
                thetarms[i] = 0;
                NLowtheta[i] = 0;
                //************************************
                trms[i]  = 0.;
                trmsdiff[i]  = 0.;
                fpdist[i]  = 0.;
                phirms[i] = 0.;
                bsdirks[i] = 0.;
                thetam[i] = 0.;
                dt[i]    = 0.;
                
                Neff[i]  = 0;
                ratio[i] = 0.;
                NhighQ[i] = 0;
                NlowQ[i] = 0;
                Nc1[i] = 0;
                
                Nlow1[i] = Nlow2[i] = Nlow3[i] = Nlow4[i] = 0;
                Nlow5[i] = Nlow6[i] = Nlow7[i] = Nlow8[i] = Nlow9[i]= 0;
            }
            
            // Experimental
            ncomb = ncomb1 = ncomb2 = ncomb3 = 0;
            for (Int_t i=0; i<MAXCOMB; i++) dalpha[i]= dalpha1[i] = dalpha2[i] = theta[i]  = 0;
            for(Int_t i=0; i<MAXCOMB; i++) cable[i] = 0;
            
            // Clear vectors
            shet.clear();
            shecab.clear();
            aftt.clear();
            aftq.clear();
            aftcab.clear();
        }
} results_t;

typedef struct prompt_hits{
    Int_t nhits;
    Int_t *cab = new Int_t[MAXHITS];
    Float_t *t = new Float_t[MAXHITS];
    Float_t *q = new Float_t[MAXHITS];
    
    ~prompt_hits(){
        delete[] cab;
        delete[] t;
        delete[] q;
    }
} prompt_hits;

class TChain;
class TFile;
class TTree;

class Header;
class TQReal;
class LoweInfo;

class SK2p2MeV: public TObject {
  public:
    SK2p2MeV (const Float_t (*xyz)[3]=0);
    virtual ~SK2p2MeV ();
    virtual void Analyze (long entry, bool last_entry) = 0;
    virtual void Print ();
    virtual bool GetBranchValues();
    TTree* MakeOutputTree(bool is_wit);
    Bool_t CheckBadMis (const Int_t cab);
    void SetPrompt (Float_t, Float_t);
    void SetAFTGate (const Float_t gate);
    void SetWindowSize (const Float_t time);
    void SetExtra (const Int_t extr);
    void SetBadch (const Int_t nbad, const Int_t *ibad);
    void SetMisch (const Int_t nmis, const Int_t *imis);
    void SetN10Threshold (const Int_t n10);
    void SetN10CutThreshold (const Int_t n10);
    void SetVerbosity (const Int_t v=1);
    void SetNcCutFlag (const Bool_t flag = kTRUE);
    void SetHighQCutFlag (const Bool_t flag = kTRUE);
    void SetBackCutFlag (const Bool_t flag = kTRUE);
    void SetVertex (const Float_t x, const Float_t y, const Float_t z);
    void SetDarkRate(Int_t run);
    void Clear();
  
  private:
    SK2p2MeV (const SK2p2MeV& rhs);
    SK2p2MeV& operator= (const SK2p2MeV& rhs);
    
  protected:
    
    // for reading inputs
    MTreeReader* myTreeReader = nullptr;
    bool isWIT = false;
    
    // Geometry of ID PMTs
    const Float_t (*xyz)[3];
    
    // Input tree branches
    Header    *HEADER;
    TQReal    *TQI;
    TQReal    *TQA;
    LoweInfo  *LOWE;
    int *is_signal;
    
    //Input for WIT
    Float_t BONSAI[13];
    UShort_t swtrigger[3];
    Float_t mc[11];
    UShort_t cable[100000];
    Float_t t[100000];
    Float_t q[100000];
    int nhit;
    int nrunsk;
    int nsubsk;
    int nevsk;
    int ltimediff;
    int idtgsk;
    int trigid;
    float mct;
    float bsdirks;
    
    // File and tree for output
    TFile *fout;
    TTree *theOTree;
    
    // Primary vertex
    Float_t VX;
    Float_t VY;
    Float_t VZ;
    
    // Results
    results_t res;
    
    // Prompt hits
    //prompt_hits prompt;
    
    // After trigger gate, in ns
    Float_t AFT_GATE;
    
    // N10 threshold
    Int_t N10TH;
    Int_t N10cutTH;
    
    // Window size
    Float_t twin;
    
    // Extra wiggle room when reoptimizing window
    Int_t extra;
    
    // flag indicating whether cutting bad hits or not
    Bool_t removeCluster;
    Bool_t removeHighQ;
    Bool_t removeBack;
    
    // Verbosity
    Int_t verbosity;
    
    // Bad channel info
    Int_t NBAD;        // # of back channels
    const Int_t *IBAD; // Cable number of bad channels
    
    // Missing channel info
    Int_t NMIS;        // # of missing channels
    const Int_t *IMIS; // Calbe number of missing channels
    
    // Private functions
    Int_t   N200Max  (Float_t tstart, Float_t tend);
    Int_t   N200Max  (Float_t tstart, Float_t tend, Float_t &t200m);
    void    NeutronSearch (Float_t tstart, Float_t tend, Float_t TOFFSET=0.);
 
    Int_t   GetNXX   (Int_t nhits, Float_t *t, Float_t twin, Float_t tcenter);
    Int_t   GetNhits (Float_t *v, Int_t start_index, Float_t width, Int_t nhits);
    Int_t GetNhits_flag(Float_t *v, Int_t *flag, Int_t flagcut, Int_t start_index, Float_t width, Int_t nhits);
    Float_t EffCos   ( Float_t costh );
    Float_t GetWeight (const Float_t xyz[3], const Float_t v[3]);
    Float_t GetWeightThreshold (const Float_t *w, const Float_t frac=0.75);
    
    Bool_t  CheckHighQ (Float_t *qi, Int_t *flag, Int_t N10, Int_t &ncut, Float_t qth=3.);
    //by Yang Zhang
    Bool_t  CheckLowQ (Float_t *qi, Int_t *flag, Int_t N10, Int_t &ncut, Float_t qth=0.5);
    //
    Bool_t  CheckBackHits (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
                           Int_t N10, Int_t &ncut, Float_t angle=90.);
    Int_t   GetCluster (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
                        Int_t N10, Int_t &ncut, Int_t ncth=3, Float_t thr=0.97);
    Int_t   GetLowHits (Int_t *ci, Int_t N10, Float_t *wt, Float_t acceptance=0.7);
    Float_t GetThetaMean (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10);
    Float_t GetPhiRms (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10);
    Float_t GetDirKS (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10);
    Float_t GetRatio (Int_t *ci, Float_t *wt, Int_t *flag, Int_t N10);
    Float_t GetTrms (Float_t *ti, Int_t *flag, Int_t N10);
    //*********new added by Yang Zhang*****
    Float_t GetQ10 (Float_t *qi, Int_t *flag, Int_t N10);
    Float_t GetQrms (Float_t *qi, Int_t *flag, Int_t N10);
    Float_t GetQmean (Float_t *qi, Int_t *flag, Int_t N10);
    Float_t GetThetaRms (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10);
    Int_t GetNLowTheta (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10, Float_t thetaTh);
    //***************************************
    // Obsolete
    Bool_t  CheckCherenkovLike (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
                                Int_t N10, Int_t &ncut);
    Bool_t  CheckCluster (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
                          Int_t N10, Int_t &ncut, Float_t angle=10.);
    // new
     Float_t MinimizeTrms(Float_t* tiskz, Int_t* cabiz, Int_t startindex, Int_t* index, Int_t nhits, Float_t& CVX, Float_t& CVY, Float_t& CVZ, Float_t pVX, Float_t pVY, Float_t pVZ, Float_t discut);
     Float_t BasicTof(Float_t* tiskz, Int_t* cabiz, Float_t VX, Float_t VY, Float_t VZ, Int_t nhits, Int_t* index);
     void GetMinTrms (Float_t *ti, Int_t N10, Float_t &mintrms_6, Float_t &mintrms_5, Float_t &mintrms_4, Float_t &mintrms_3, Int_t *flag);
    
};

#endif // __SK2P2MEV_H__
