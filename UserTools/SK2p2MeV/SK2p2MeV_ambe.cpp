#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"

#include "tqrealroot.h"
#include "loweroot.h"

#include "SK2p2MeV_ambe.h"

#include <iostream>

SK2p2MeV_ambe::SK2p2MeV_ambe (const Float_t (*geomxyz)[3]) 
    : SK2p2MeV (geomxyz)
{
    // additional branches
    head0   = new Header;
    lowe0   = new LoweInfo;
}

SK2p2MeV_ambe::~SK2p2MeV_ambe(){
    //    std::cout << "ntot/npair: " << nentries << " " << n_pair << std::endl;
    theOTree->ResetBranchAddresses();
    
    delete head0;
    delete lowe0;
}

void SK2p2MeV_ambe::Print()
{
    // Show current status

    std::cout << " Class: SK2p2MeV_ambe" << std::endl;
    std::cout << "    AFT_GATE: " << AFT_GATE << std::endl;
}

bool SK2p2MeV_ambe::GetBranchValues(){
    // get base class branches
    return SK2p2MeV::GetBranchValues();
}

bool SK2p2MeV_ambe::Initialise(MTreeReader* reader){
    
    myTreeReader = reader;
    ch = myTreeReader->GetTree();
    
    // Add additional branches
    theOTree->Branch("HEADER", "Header", &head0, 1024*1024, 0);
    theOTree->Branch("LOWE", "LoweInfo", &lowe0, 1024*1024, 0);
    theOTree->Branch("totpe", &totpe,  "totpe/F");
    
    return true;
}

void SK2p2MeV_ambe::Analyze(long entry, bool last_entry)
{
    
    // printouts
    // =========
    
    //if (entry % 10 == 0) std::cout << "Entry number " << entry << "/" << nentries << std::endl;
    
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
    
    // analysis
    // ========
    
    if ( HEADER->idtgsk & 0x10000000 ) {// Primary SHE
        // Clear results of previous event
        res.Clear();
        
        // New results
        *head0 = *HEADER;
        *lowe0 = *LOWE;
        
        // Cal totpe
        totpe = 0.;
        for ( Int_t j=0; j<TQI->nhits; j++) {
            Int_t flag = TQI->cables[j]&0xFFFF0000; // Get upper 16 bit
            if ( (flag>>16)&0x1 ) { // in 1.3us
                totpe += TQI->Q[j];
            }
        }
        
        // Check N200 in primary event
        res.N200M = N200Max(6000., 40000., res.T200M);
        
        // Do 2.2MeV search in primary event (~35us)
        // Save AFT hit info for dark noise search
        if (!last_entry){
            ch->GetEntry(entry + 1);
            GetBranchValues();
            if ( HEADER->idtgsk & 0x20000000 ){
                for(int k = 0; k < TQI->nhits; k++){
                    if (TQI->T[k] > cut_window) break;
                    int cb = TQI->cables[k]&0xFFFF;
                    if (cb < 0 || cb > MAXPM) return;
                    if ( ! ((TQI->cables[k]&0xFFFF0000)&0x20000 ) ) return;
                    if ( CheckBadMis(cb) ) return;
                    res.aftt.push_back(TQI->T[k] + she_tmax);
                    res.aftcab.push_back(cb);
                }
            }
            ch->GetEntry(entry);
            GetBranchValues();
        }
        // Perform neutron search for SHE
        NeutronSearch (1500, 40000.);
        // Clear AFT hit info used for dark noise search
        res.aftt.clear();
        res.aftcab.clear();
        
        // Save SHE hit info for dark noise search with AFT trigger
        for(int k = 0; k < TQI->nhits; k++){
            int cb = TQI->cables[k]&0xFFFF;
            if (cb < 0 || cb > MAXPM) return;
            if ( ! ((TQI->cables[k]&0xFFFF0000)&0x20000 ) ) return;
            if ( CheckBadMis(TQI->cables[k]&0xFFFF) ) return;
            res.shet.push_back(TQI->T[k]);
            res.shecab.push_back(cb);
        }
        pre_she_good = kTRUE;
        tstart = HEADER->t0;
    }
    else if ( (HEADER->idtgsk & 0x20000000) && pre_she_good ) { // AFT
        // SHE and AFT should be two successive events.
        // Some AFT events are observed to have no SHE accomponied.
        if ( HEADER->nevsk - pre_nevsk != 1 ) { 
            pre_nevsk = HEADER->nevsk;
            return;
        }
        
        n_pair ++;
        
        // Check N200 in AFT data
        Float_t t200m;
        Int_t n200m = N200Max(0., AFT_GATE, t200m);
        if ( n200m > res.N200M ) {
            res.N200M = n200m;
            res.T200M = t200m + 35000.; // shift 35 us
        }
        
        // Do 2.2 MeV search
        tend = HEADER->t0;
        deltat = (tend - tstart)/1.92;
        NeutronSearch (0., AFT_GATE, 100000. + deltat - 35000.); // shift t0 in AFT by 100us
        
        // Fill the output tree 
        theOTree->Fill();
        res.Clear();
        
        pre_she_good = kFALSE;
    }
    else if ((HEADER->idtgsk & 0x4000) && ! (HEADER->idtgsk & 0x08)) {//Ni trig.
        // Clear results of previous event
        res.Clear();
        
        // New results
        *head0 = *HEADER;
        *lowe0 = *LOWE;
        
        // Check N200 in AFT data
        res.N200M = N200Max(0., AFT_GATE);
        
        if(verbosity == 2) std::cout << "hits Ni" << std::endl;
        // Do 2.2 MeV search
        SetPrompt(0, 1050);
        NeutronSearch (0., AFT_GATE);

        // Fill the output tree 
        theOTree->Fill();
        pre_she_good = kFALSE;
    }
    else if ( HEADER->idtgsk & 0x00000800 ) {// Random Wide Trigger (BG run)
        // Clear results of previous event
        res.Clear();
        
        // New results
        *head0 = *HEADER;
        *lowe0 = *LOWE;
        
        // Cal totpe
        totpe = 0.;
        for ( Int_t j=0; j<TQI->nhits; j++) {
            Int_t flag = TQI->cables[j]&0xFFFF0000; // Get upper 16 bit
            if ( (flag>>16)&0x1 ) { // in 1.3us
                totpe += TQI->Q[j];
            }
        }
        
        // Check N200 in random wide trigger
        res.N200M = N200Max(0., 1000000., res.T200M);
        // Do 2.2MeV search in randowm wide trigger event (1000us)
        NeutronSearch (0, 1000000.);
        
        // Fill the output tree 
        theOTree->Fill();
        pre_she_good = kFALSE;
    }
    else {
        // other events.
        pre_she_good = kFALSE;
        res.Clear();
    }
    
    pre_nevsk = HEADER->nevsk;
    
    return;
}
