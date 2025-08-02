#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"

#include "thirdredvars.h"

#include "SK2p2MeV_merged.h"

#include <iostream>
#include <cstdlib>

SK2p2MeV_merged::SK2p2MeV_merged (const Float_t (*geomxyz)[3]) 
    : SK2p2MeV (geomxyz)
{
    MU = new MuInfo;
    thirdred = new ThirdRed;
    
    // additional branches
    head0   = new Header;
    lowe0   = new LoweInfo;
    mu0     = new MuInfo;
    third0  = new ThirdRed;
    
}

SK2p2MeV_merged::~SK2p2MeV_merged(){
    theOTree->ResetBranchAddresses();
    
    delete MU;
    delete thirdred;
    
    delete head0;
    delete lowe0;
    delete mu0;
    delete third0;
    
}

void SK2p2MeV_merged::Print()
{
    // Show current status
    std::cout << " Class: SK2p2MeV_merged" << std::endl;
}

void SK2p2MeV_merged::SetEneThreshold(const Double_t ethre)
{
    // Set bsenergy threshold
    res.EneThre = ethre;
}

bool SK2p2MeV_merged::GetBranchValues(){
    // get base class branches
    bool get_ok = SK2p2MeV::GetBranchValues();
    
    // get additional branch variables - MU and ThirdRed propagated to output TTree
    get_ok &= (myTreeReader->Get("MU", MU));
  //get_ok &= (myTreeReader->Get("multispa_dist", &multispa));
    //get_ok &= (myTreeReader->Get("ThirdRed", thirdred));
    
    return get_ok;
}

bool SK2p2MeV_merged::Initialise(MTreeReader* reader){
    
    myTreeReader = reader;
    ch = myTreeReader->GetTree();
    
    // Additional output branches
    // ==========================
    theOTree->Branch("HEADER", "Header", &head0, 1024*1024, 0);
    theOTree->Branch("LOWE", "LoweInfo", &lowe0, 1024*1024, 0);
    theOTree->Branch("MU", "MuInfo", &mu0, 1024*1024, 0);
    theOTree->Branch("ThirdRed", "ThirdRed", &third0, 1024*1024, 0);
    //theOTree->Branch("multispa_dist", &multispa, "multispa_dist/I");
    theOTree->Branch("totpe", &totpe,  "totpe/F");
    theOTree->Branch("nnhits", &nnhits,  "nnhits/I");
    theOTree->Branch("type", &type, "type/I");
    
    //multispa = 1000000;
    
    return true;
    
}

////////////////////////////////////////////////////////
// =====================================================
////////////////////////////////////////////////////////

void SK2p2MeV_merged::Analyze(long entry, bool last_entry)
{
    
    // Clear results of previous event
    res.Clear();
    type = 0;
    
    // get event data
    GetBranchValues();
    
    // Set energy threshold and AFT gate
    // these are part of args passed to NeutronSearch,
    // and then subsequently used as a cut on lowe energy
    // (NeutronSearch calls bonsai internally... do
    // we cut on the right energy though? previously it cut on LOWE->bsenergy)
    double enethre, AFT_gate;
    if( 61525 <= nrun && nrun <= 62427){
      enethre=10.0;
      AFT_gate=350e3;
    }
    if(62428 <= nrun && nrun <= 68670){
      enethre=10.0;
      AFT_gate=500e3;
    }
    if(68671 <= nrun ){//&& nrun <= 77958){
      enethre=8.0;
      AFT_gate=500e3;
    }
    // FIXME this was before the OTree->Fill call, so acted as a cut
    // although bonsai is called within NeutronSearch, i don't see it updating the LOWE class member 
    // FIXME FIXME move to upstream RunWiseCut tool
    if(LOWE->bsenergy < enethre) return;
    SetEneThreshold(enethre);
    // likewise this doesn't appear to do anything, AFTGate is passed below to NeutronSearch
    // but SetAFTGate just sets an unused internal variable.
    SetAFTGate(AFT_gate);
    
    /*
    // FIXME FIXME what is this about? was only for AFT? not all 'type's written, so cut?
    // these thresholds probably need changing since we moved to merged SHE+AFT
    if ( TQI->T[TQI->nhits - 1] < 400000 ) type = 1;      // 300 us
    else if ( TQI->T[TQI->nhits - 1] < 600000 ) type = 2; // 500 us
    else type = 3;                                        // 800 us
    if(type==3) return true;  // FIXME FIXME FIXME previously in relic OTree->Fill was bypassed for type 3
    */
    
    // debug prints
    if ( verbosity == 1 ) {
        if ( entry%1000 == 0 ) {
            std::cout << "entry/nrunsk/nevsk/idtgsk/nqiskz: " << entry
                 << " " << HEADER->nrunsk
                 << " " << HEADER->nevsk
                 << " " << std::hex << HEADER->idtgsk
                 << " " << std::dec << TQI->nhits << std::endl;
        }
    }
    else if ( verbosity == 2) {
        std::cout << "entry/nrunsk/nevsk/idtgsk/nqiskz/tdiff: " << entry
             << " " << HEADER->nrunsk 
             << " " << HEADER->nevsk 
             << " " << std::hex << HEADER->idtgsk 
             << " " << std::dec << TQI->nhits 
             << " " << LOWE->ltimediff << std::endl;
    }
    
    // FIXME FIXME FIXME we shouldn't need to do this
    // replace with just TTree::CloneTree so that branches point at the same object
    *head0 = *HEADER;
    *lowe0 = *LOWE;
    *mu0   = *MU;
    *third0 = *thirdred;
    
    // bad ch file import (Takeuchi-code import by H.Ito, August, 2019)
    if(nrun != HEADER->nrunsk || nsub != HEADER->nsubsk ){
        nrun = HEADER->nrunsk;
        nsub = HEADER->nsubsk;
        int imask = 23;
        int log_level = 0;
        int istat;
        combad_.log_level_skbadch = log_level;
        skbadopt_(&imask);
        skbadch_(&nrun, &nsub, &istat);
        SetBadch(combad_.nbad, combad_.ibad);
        
        if(MISCH==nullptr){
            NMIS = combad_.nbad;
            MISCH = new int[NMIS];
        } else if(combad_.nbad > NMIS){
            delete[] MISCH;
            NMIS = combad_.nbad;    /// FIXME typo? is this supposed to be nbad again? maybe.
            MISCH = new int[NMIS];
        }
        for(int j=0;j<NMIS;j++){
            MISCH[j]=combad_.isqbad[j];
        }
        SetMisch (NMIS, MISCH);
    }
    
    nnhits = TQI->nhits;
    
    // Calculate totpe
    totpe = 0.;
    for ( Int_t j=0; j<TQI->nhits; j++) {
        Int_t flag = TQI->cables[j]&0xFFFF0000; // Get upper 16 bit
        if ( (flag>>16)&0x1 ) { // in 1.3us
            totpe += TQI->Q[j];
        }
    }
    
    // Set primary vertex
    VX = LOWE->bsvertex[0];
    VY = LOWE->bsvertex[1];
    VZ = LOWE->bsvertex[2];
    
    // "Some AFT events are observed to have no SHE" << comment in other wrapper classe. 
    // Well, that doesn't seem right, but we probably don't have those as AFT-only probably got skipped by upstream ToolChain.
    // SHE without AFT - will surely happen and will be kept... do we need to handle these differently? FIXME
    
    // FIXME set appropriate time window...
    // perviously in SK2p2MeV_relic:
    // 1) N200Max: SHE searched 6k-40k, then AFT searched 0-AFT_gate, with 35k added to AFT time
    // 2) NeutronSearch: SHE searched 6k-40k, then AFT searched 0-(AFT_gate-cut_window), with AFT hit times shifted by
    //    35k-100k+time between last SHE hit and first AFT hit... XXX huh???
    // since we merge both sets of hits, i don't think any time shift makes sense for NeutronSearch
    // so we only need to decide a window length... and it probably just needs to be "all of ths hits" right??
    // well, we start from 6k to exclude muon afterpulsing, IIRC
    res.N200M = N200Max(6000., 40000 + AFT_GATE, res.T200M);
    NeutronSearch (1050, 40000 + AFT_GATE - cut_window);
    
    // Fill the output tree 
    theOTree->Fill();
    return;
    
}
