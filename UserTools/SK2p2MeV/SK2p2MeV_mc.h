/*****************************************************
 * SK2p2MeV_mc.h                                     *
 *                                                   *
 * Author : Haibing ZHANG                            * 
 *          <zhanghb02@gmail.com>                    *
 * Date   : 04/22/10                                 *
 * Comment: Derived from SK2p2MeV                    * 
 *****************************************************/

#ifndef __SK2P2MEV_MC_H__
#define __SK2P2MEV_MC_H__

///////////////////////////////////////////////////////
//                                                   //  
//                                                   //
///////////////////////////////////////////////////////

#include "SK2p2MeV.h"

class MCInfo;
class ThirdRed;

class SK2p2MeV_mc: public SK2p2MeV {
  public:
    SK2p2MeV_mc (const Float_t (*xyz)[3]=0);
    ~SK2p2MeV_mc ();
    bool Initialise(MTreeReader* reader);
    bool GetBranchValues();
    void Analyze (long entry, bool last_entry);
    void Print ();
    void SetSmearFlag (Bool_t flag=kTRUE);
    void SetMCFlag (Bool_t flag=kTRUE);
    void SetVertexResolution (const Float_t x, const Float_t y, const Float_t z);
    void SetShiftFlag (Bool_t flag=kTRUE);
    void SetShiftDistance (const Float_t dis);
    void SetSeed(int seed);

  private:
    SK2p2MeV_mc (const SK2p2MeV_mc& rhs);
    SK2p2MeV_mc& operator= (const SK2p2MeV_mc& rhs);

  protected:
    
    // additional input variables
    MCInfo *MC;
    ThirdRed *thirdred;
    
    // configuration variables
    Bool_t  fSmear;
    Bool_t  useMC;
    Float_t sigma_xyz[3]; // vertex resolution
    Bool_t  fShift;
    Float_t rshift;
    int rseed;
    
    // additional output variables
    Header *head0 = nullptr;
    LoweInfo *lowe0 = nullptr;
    MCInfo *mc0 = nullptr;
    ThirdRed *third0 = nullptr;
    float smearedvertex[3];
    
    // variables that need to be carried over between Analyze calls
    int prevrun = 0;
    int prevsub = 0;
};

#endif // __SK2P2MEV_MC_H__
