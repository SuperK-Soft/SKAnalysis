/*****************************************************
 * SK2p2MeV_merged.h                                   *
 *                                                   *
 * Author : Haibing ZHANG                            * 
 *          <zhanghb02@gmail.com>                    *
 * Date   : 02/03/10                                 *
 * Comment: Derived from SK2p2MeV                    * 
 *****************************************************/

#ifndef __SK2P2MEV_MERGED_H__
#define __SK2P2MEV_MERGED_H__

///////////////////////////////////////////////////////
//                                                   //  
//                                                   //
///////////////////////////////////////////////////////

#include "SK2p2MeV.h"

class MuInfo;
class ThirdRed;

class SK2p2MeV_merged: public SK2p2MeV {
  public:
    SK2p2MeV_merged (const Float_t (*xyz)[3]=0);
    ~SK2p2MeV_merged();
    bool Initialise(MTreeReader* reader);
    bool GetBranchValues();
    void Analyze (long entry, bool last_entry);
    void SetEneThreshold(const Double_t ethre);
    void Print ();

  private:
    SK2p2MeV_merged (const SK2p2MeV_merged& rhs);
    SK2p2MeV_merged& operator= (const SK2p2MeV_merged& rhs);
    
  protected:
    
    // additional input variables
    MuInfo *MU;         // Parent muon info
    ThirdRed *thirdred; // Third reduction observables
    
    // additional output variables
    Header   *head0 = nullptr;
    LoweInfo *lowe0 = nullptr;
    MuInfo   *mu0   = nullptr;
    ThirdRed *third0 = nullptr;
    
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
    int nrun=0, nsub=0;
    unsigned int tstart, tend;
    
    // for loading input data
    TTree* ch = nullptr;
};

#endif // __SK2P2MEV_MERGED_H__
