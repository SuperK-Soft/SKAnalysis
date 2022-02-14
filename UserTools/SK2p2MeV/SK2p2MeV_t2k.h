/*****************************************************
 * SK2p2MeV_t2k.h                                   *
 *                                                   *
 * Author : Haibing ZHANG                            * 
 *          <zhanghb02@gmail.com>                    *
 * Date   : 03/24/10                                 *
 * Comment: Derived from SK2p2MeV                    * 
 *****************************************************/

#ifndef __SK2P2MEV_T2K_H__
#define __SK2P2MEV_T2K_H__

///////////////////////////////////////////////////////
//                                                   //  
//  To evaluate BG using T2K dummay data.            //
//                                                   // 
///////////////////////////////////////////////////////

#include "SK2p2MeV.h"

class SK2p2MeV_t2k: public SK2p2MeV {
  public:
    SK2p2MeV_t2k (const Float_t (*xyz)[3]=0);
    ~SK2p2MeV_t2k();
    bool Initialise(MTreeReader* reader, bool random);
    bool GetBranchValues();
    void Analyze (long entry, bool last_entry);
    void Print ();

  private:
    SK2p2MeV_t2k (const SK2p2MeV_t2k& rhs);
    SK2p2MeV_t2k& operator= (const SK2p2MeV_t2k& rhs);
    
  protected:
    
    // configuration variables
    bool random = false;
    
    // additional output variables
    Header *head0 = nullptr;
    LoweInfo *lowe0 = nullptr;
    Float_t totpe;
    Int_t nnhits;
    
    // variables that need to be carried over between Analyze calls
    int nrun=0,nsub=0;
    
};

#endif // __SK2P2MEV_T2K_H__
