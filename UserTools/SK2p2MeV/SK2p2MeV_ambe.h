/*****************************************************
 * SK2p2MeV_ambe.h                                   *
 *                                                   *
 * Author : Haibing ZHANG                            * 
 *          <zhanghb02@gmail.com>                    *
 * Date   : 03/25/10                                 *
 * Comment: Derived from SK2p2MeV                    * 
 *****************************************************/

#ifndef __SK2P2MEV_AMBE_H__
#define __SK2P2MEV_AMBE_H__

///////////////////////////////////////////////////////
//                                                   //  
//  For AmBe analysis.                               //
//                                                   //
///////////////////////////////////////////////////////

#include "SK2p2MeV.h"

class SK2p2MeV_ambe: public SK2p2MeV {
  public:
    SK2p2MeV_ambe (const Float_t (*xyz)[3]=0);
    ~SK2p2MeV_ambe();
    bool Initialise(MTreeReader* reader);
    bool GetBranchValues();
    void Analyze(long entry, bool last_entry);
    void Print ();

  private:
    SK2p2MeV_ambe (const SK2p2MeV_ambe& rhs);
    SK2p2MeV_ambe& operator= (const SK2p2MeV_ambe& rhs);
    
  protected:
    
    // additional output variables
    Header   *head0 = nullptr;
    LoweInfo *lowe0 = nullptr;
    Float_t totpe;
    
    // variables that need to be carried over between Analyze calls
    Int_t  pre_nevsk = -1;
    Bool_t pre_she_good = kFALSE;
    Int_t  n_pair = 0; // Number of good SHE+AFT pairs
    unsigned int tstart, tend;
    float deltat;
    
    // for loading SHE / AFT commons
    TTree* ch = nullptr;
    
};

#endif // __SK2P2MEV_AMBE_H__
