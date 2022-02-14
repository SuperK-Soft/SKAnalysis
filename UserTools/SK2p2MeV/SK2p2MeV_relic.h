/*****************************************************
 * SK2p2MeV_relic.h                                   *
 *                                                   *
 * Author : Haibing ZHANG                            * 
 *          <zhanghb02@gmail.com>                    *
 * Date   : 02/03/10                                 *
 * Comment: Derived from SK2p2MeV                    * 
 *****************************************************/

#ifndef __SK2P2MEV_RELIC_H__
#define __SK2P2MEV_RELIC_H__

///////////////////////////////////////////////////////
//                                                   //  
//                                                   //
///////////////////////////////////////////////////////

#include "SK2p2MeV.h"

class MuInfo;
class ThirdRed;

class SK2p2MeV_relic: public SK2p2MeV {
  public:
    SK2p2MeV_relic (const Float_t (*xyz)[3]=0);
    ~SK2p2MeV_relic();
    bool Initialise(MTreeReader* reader, bool fake, int seed, const char* t2kinfo, const char* t2kdir, int timebin);
    bool GetBranchValues();
    void Analyze (long entry, bool last_entry);
    void SetEneThreshold (const Double_t ethre);
    void Print ();

  private:
    SK2p2MeV_relic (const SK2p2MeV_relic& rhs);
    SK2p2MeV_relic& operator= (const SK2p2MeV_relic& rhs);
    
  protected:
    
    // additional input variables
    MuInfo *MU;         // Parent muon info
    ThirdRed *thirdred; // Third reduction observables
    
    // additional output variables
    Header   *head0 = nullptr;
    LoweInfo *lowe0 = nullptr;
    MuInfo   *mu0   = nullptr;
    ThirdRed *third0 = nullptr;
    //TQReal *tq0    = nullptr;
    //TQReal *tqa0   = nullptr;
    
    // configuration / processing variables
    Float_t totpe;
    Int_t nnhits;
    int multispa;
    Int_t type;   // Event type (0=primary only, 1=350 us, 2=500 us, 3=800 us, 4=other)
    
    // SK2p2MeV::SetMisch takes a pointer to a noisy PMT array, but it does not
    // make a new array from it! It just takes the pointer and expects that to remain valid!
    int* MISCH = nullptr;
    int NMIS = 0;
    
    // variables that need to be carried over between Analyze calls
    Int_t  pre_nevsk = -1;
    Bool_t pre_saved = kFALSE;
    Int_t  n_pair = 0;  // Number of good SHE+AFT pairs
    int nrun=0, nsub=0;
    unsigned int tstart, tend;
    
    // not sure if these are carried over between calls but should be fine either way
    float tbegin = -1; // Timing of the first event in SHE trigger
    float tmax = -1;   // Timing of the first event in SHE trigger
    float deltat;
    
    // additional variables if using fake data
    bool fake = false;
    TChain *t2kch = nullptr;
    Header *t2khead = nullptr;
    TQReal *t2ktq = nullptr;
    
    // for loading SHE / AFT entries
    TTree* ch = nullptr;
};

#endif // __SK2P2MEV_RELIC_H__
