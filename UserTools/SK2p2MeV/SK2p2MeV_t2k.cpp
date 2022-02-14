#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"

#include "tqrealroot.h"
#include "loweroot.h"

#include "SK2p2MeV_t2k.h"

#include <iostream>
#include "skheadC.h"
#include "skparmC.h"
#include "geopmtC.h"
#include "skbadcC.h"

SK2p2MeV_t2k::SK2p2MeV_t2k (const Float_t (*geomxyz)[3]) 
    : SK2p2MeV (geomxyz)
{
    head0 = new Header;
    lowe0 = new LoweInfo;
}

SK2p2MeV_t2k::~SK2p2MeV_t2k(){
    
    theOTree->ResetBranchAddresses();
    
    delete head0;
    delete lowe0;
}

void SK2p2MeV_t2k::Print()
{
    // Show current status
    
    std::cout << " Class: SK2p2MeV_t2k" << std::endl;
    std::cout << "    AFT_GATE: " << AFT_GATE << std::endl;
}

bool SK2p2MeV_t2k::Initialise(MTreeReader* reader, bool random){
    
    myTreeReader = reader;
    
    // Add additional branches
    theOTree->Branch("HEADER", "Header", &head0, 1024*1024, 0);
    theOTree->Branch("LOWE", "LoweInfo", &lowe0, 1024*1024, 0);
    theOTree->Branch("totpe", &totpe,  "totpe/F");
    theOTree->Branch("nnhits", &nnhits,  "nnhits/I");
    
    if (random){
        UInt_t seed = HEADER->nevsk + HEADER->nrunsk;
        gRandom->SetSeed(seed);
    }
    
    return true;
}

bool SK2p2MeV_t2k::GetBranchValues(){
    return SK2p2MeV::GetBranchValues();
}

void SK2p2MeV_t2k::Analyze (long entry, bool last_entry)
{
    
    // bad ch file import (Takeuchi-code import by H.Ito, August, 2019)
    if(nrun != HEADER->nrunsk || nsub != HEADER->nsubsk ){
        nrun = HEADER->nrunsk;
        nsub = HEADER->nsubsk;
        int imask = 23;
        int log_level = 0;//1;
        skheadg_.sk_geometry = 4; // SK-IV
        combad_.log_level_skbadch = log_level; 
        skbadopt_(&imask);
        int istat;
        skbadch_(&nrun, &nsub, &istat);
        SetBadch(combad_.nbad, combad_.ibad);
        
        const int NMIS=combad_.nbad;
        int MISCH[NMIS];
        for(int j=0;j<NMIS;j++){
            MISCH[j]=combad_.isqbad[j];
        }
        SetMisch (NMIS, MISCH);
    }
    // output info
    if ( entry < 0 ) verbosity = 2; // 2011-04-27 for Nlow check
    else verbosity = 1;
    
    if( entry%100==0 ) std::cout << "processing " << entry << " th events" << std::endl;
    
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
    
    if ( ! (HEADER->idtgsk & 0x80000000) ) return; // T2K trig.
    // NOTE: It is observed that some events have very large number of
    //       hits, e.g. 110000. These are very likely non-physical events.
    //       They may cause program crash.
    if ( TQI->nhits > 70000 ) return;
    
    if(random){
        // Using random vertex for backgrund study
        VZ = gRandom->Uniform(-1610., 1610.);
        do {
            VX = gRandom->Uniform(-1490., 1490.);
            VY = gRandom->Uniform(-1490., 1490.);
        } while ( sqrt(VX*VX + VY*VY) > 1490.);
    }
    
    // Clear results of previous event
    res.Clear();
    // New results
    *head0 = *HEADER;
    *lowe0 = *LOWE;
    lowe0->bsvertex[0] = VX;
    lowe0->bsvertex[1] = VY;
    lowe0->bsvertex[2] = VZ;
    
    // Cal totpe
    totpe = 0.;
    for ( Int_t j=0; j<TQI->nhits; j++) {
        Int_t flag = TQI->cables[j]&0xFFFF0000; // Get upper 16 bit
        if ( (flag>>16)&0x1 ) { // in 1.3us
            totpe += TQI->Q[j];
        }
    }
    
    nnhits = TQI->nhits ;
    // Split window in 2 to reproduce timing of real data
    float tmin = 0;
    float aft_gate = 488000;
    for(int sign = -1; sign <= 1; sign+=2){
        float tbegin = sign > 0 ? tmin : -aft_gate;
        float tend = sign > 0 ? aft_gate : tmin;
        // Check N200 in AFT data
        res.Clear();
        res.N200M = N200Max(tbegin - cut_window, tend + cut_window);
        
        // Eliminate events with "big" peak, e.g. muons
        if ( res.N200M > 10000 ) return;
        
        // Do 2.2 MeV search
        std::cout << "t2k total hits " << TQI->nhits << std::endl;
        NeutronSearch (tbegin, tend);
        
        // Fill the output tree 
        theOTree->Fill();
    }
    
}
