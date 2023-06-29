#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"

#include "thirdredvars.h"

#include "SK2p2MeV_relic.h"

#include <iostream>
#include <cstdlib>

SK2p2MeV_relic::SK2p2MeV_relic (const Float_t (*geomxyz)[3]) 
    : SK2p2MeV (geomxyz)
{
    MU = new MuInfo;
    //    thirdred = new ThirdRed;
    
    // additional branches
    head0   = new Header;
    lowe0   = new LoweInfo;
    mu0     = new MuInfo;
    //    third0  = new ThirdRed;
    //tq0   = new TQReal;
    //tqa0  = new TQReal;
    
}

SK2p2MeV_relic::~SK2p2MeV_relic(){
    theOTree->ResetBranchAddresses();
    
    delete MU;
    //delete thirdred;
    
    delete head0;
    delete lowe0;
    delete mu0;
    //delete third0;
    //delete tq0;
    //delete tqa0;
    
    if(t2khead) delete t2khead;
    if(t2ktq) delete t2ktq;
    if(t2kch) delete t2kch;
}

void SK2p2MeV_relic::Print()
{
    // Show current status

    std::cout << " Class: SK2p2MeV_relic" << std::endl;
    std::cout << "    AFT_GATE: " << AFT_GATE << std::endl;
}

void SK2p2MeV_relic::SetEneThreshold (const Double_t ethre)
{
    // Set bsenergy threshold
    res.EneThre = ethre;
}

bool SK2p2MeV_relic::GetBranchValues(){
    // get base class branches
    bool get_ok = SK2p2MeV::GetBranchValues();
    
    // get additional branch variables - MU and ThirdRed propagated to output TTree
    get_ok &= (myTreeReader->Get("MU", MU));
  //get_ok &= (myTreeReader->Get("multispa_dist", &multispa));
    //    get_ok &= (myTreeReader->Get("ThirdRed", thirdred));
    
    return get_ok;
}

bool SK2p2MeV_relic::Initialise(MTreeReader* reader, bool fake, int seed, const char* t2kinfo, const char* t2kdir, int timebin){
    
    myTreeReader = reader;
    ch = myTreeReader->GetTree();
    
    // Additional output branches
    // ==========================
    theOTree->Branch("HEADER", "Header", &head0, 1024*1024, 0);
    theOTree->Branch("LOWE", "LoweInfo", &lowe0, 1024*1024, 0);
    theOTree->Branch("MU", "MuInfo", &mu0, 1024*1024, 0);
    //    theOTree->Branch("ThirdRed", "ThirdRed", &third0, 1024*1024, 0);
    //theOTree->Branch("multispa_dist", &multispa, "multispa_dist/I");
    //theOTree->Branch("TQREAL", "TQReal", &tq0, 1024*1024, 0);
    //theOTree->Branch("TQAREAL", "TQReal", &tqa0, 1024*1024, 0);
    theOTree->Branch("totpe", &totpe,  "totpe/F");
    theOTree->Branch("nnhits", &nnhits,  "nnhits/I");
    theOTree->Branch("type", &type, "type/I");
    
    //multispa = 1000000;
    
    // Here add processing for fake data
    if (fake){
        if(verbosity > 0) std::cout << "Running with fake data" << std::endl;
        
        if(t2kch==nullptr) t2kch = new TChain("data");
        if(t2khead==nullptr) t2khead = new Header;
        if(t2ktq==nullptr) t2ktq = new TQReal;
        
        // Find seed and possibly the time period
        int fixed_time = (timebin > 0) ? 1 : 0;
        // Find run time bin and the corresponding T2K runs to use
        // Info in a text file
        std::map<int,int> runs;
        FILE *fruns = fopen(t2kinfo, "r");
        int run, bin, dum;
        int startbin = -1, oldbin = -1;
        int count = 0;
        // Fill run table and find time bin for current run
//        ch->GetEntry(0);
        if (timebin < 0){
            while(!feof(fruns)){
                fscanf(fruns, "%d %d %d\n", &run, &bin, &dum);
                if (bin != oldbin){
                    if (timebin != -1) break;
                    startbin = count;
                    oldbin = bin;
                }
                if (run >= HEADER->nrunsk){
                    timebin = bin;
                }
                std::pair<int,int> p(run,dum);
                runs.insert(p);
                count++;
            }
            // Remove everything not in the right time bin
            // (for the first bin, nothing to remove)
            // So at the end we have a hash table that, for a given run, tells us
            // whether it is dummy or real bin
            if (startbin > 0)
                for(int q = 0; q < startbin; q++)
                    runs.erase(runs.begin());
        }
        else{
            while(!feof(fruns)){
                fscanf(fruns, "%d %d %d\n", &run, &bin, &dum);
                if (bin == timebin){
                    std::pair<int,int> p(run,dum);
                    runs.insert(p);
                }
                if (bin > timebin) break;
            }
        }
        fclose(fruns);
        if(verbosity > 1) std::cout << "Time bin is " << timebin  << " with " << runs.size() << " runs" << std::endl;
        if(verbosity > 1) std::cout << "Seed is " << seed << std::endl;
        
        // Load T2K data for the run period we want
        //dummy trigger chain
        if(verbosity > 1) std::cout << "Loading T2K events..." << std::endl;
        for(std::map<int,int>::iterator iter = runs.begin(); iter != runs.end(); iter++){
            int irun = iter->first;
            char inname[500];
            for (int k=0; k < 9; k++){
                sprintf(inname, "%s/inbtwRun%d_%d/t2k.0%d*.root", t2kdir, k, k+1, irun);
                t2kch->Add(inname);
            }
        }
        if(verbosity > 1) std::cout << "T2K events: " << t2kch->GetEntries() << "entries" << std::endl;
        t2kch->SetBranchAddress("TQREAL", &t2ktq);
        t2kch->SetBranchAddress("HEADER", &t2khead);
        srand(seed);
    }
    
    pre_nevsk = -1;
    pre_saved = kFALSE;
    
    return true;
    
}

////////////////////////////////////////////////////////
// =====================================================
////////////////////////////////////////////////////////

void SK2p2MeV_relic::Analyze(long entry, bool last_entry)
{
    
    // bad ch file import (Takeuchi-code import by H.Ito, August, 2019)
    if(nrun != HEADER->nrunsk || nsub != HEADER->nsubsk ){
        nrun = HEADER->nrunsk;
        nsub = HEADER->nsubsk;
        int imask = 23;
        int log_level = 0;//1;
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
    double enethre, AFT_gate;
    // Set energy threshold and AFT gate
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
    SetAFTGate(AFT_gate);
    SetEneThreshold(enethre);
    
    // If fake data fill in T2K dummy instead of real data
    if(fake){
        type = 0;
        //[>*****************************************
        // Uncomment this part and comment out the following tmax 
        // assignement if you want to merge SHE and AFT for fakes.
        // First get SHE events that are followed by an AFT event
#ifdef MERGED
        if ( !(HEADER->idtgsk & 0x10000000) ) return;
        if (!last_entry){
            int nev_current = HEADER->nevsk;
            ch->GetEntry(entry+1);
            GetBranchValues();
            int nev_new = HEADER->nevsk;
            if ( !(HEADER->idtgsk & 0x20000000) || (nev_new - nev_current != 1)) return;
            tend = HEADER->t0;
        }
        else return;
        deltat = (tend - tstart)/1.92;
        ch->GetEntry(entry);
        GetBranchValues();
        res.aftt.clear();
        *head0 = *HEADER;
        *lowe0 = *LOWE;
        *mu0   = *MU;
	//        *third0 = *thirdred;
        float tmax = AFT_GATE + 35000;
#endif
        //*******************************/
        // Find maximum hit time in real data
        // This sets the limit for the SHE trigger
        // Reproduce SHE window of the real event
        // If AFT, hits have already been stored
        if (res.aftt.size() > 0){
            ch->GetEntry(entry);
            GetBranchValues();
            tend = HEADER->t0;
            deltat = (tend - tstart)/1.92;
            (TQI->cables).clear();
            (TQI->T).clear();
            (TQI->Q).clear();
            TQI->nhits = 0;
            for(int q = 0; q < res.aftt.size(); q++){
                (TQI->cables).push_back((res.aftcab[q]|0x20000));
                (TQI->T).push_back(res.aftt[q] - deltat);
                (TQI->Q).push_back(res.aftq[q]);
                TQI->nhits++;
            }
            res.aftt.clear();
            res.aftq.clear();
            res.aftcab.clear();
        }
        else{
            // If SHE, fill for the first time
            // Use begin and end times of the real SHE data
            tbegin = TQI->T[0];
#ifndef MERGED
            tmax = TQI->T[TQI->nhits - 1]; // Comment out for SHE+AFT together mode
#endif
            (TQI->cables).clear();
            (TQI->T).clear();
            (TQI->Q).clear();
            TQI->nhits = 0;
            // Inject hits
            int nt2k = t2kch->GetEntries();
            int t2kentry = (int) (random()/((double) RAND_MAX) * nt2k);
            do{
                t2kch->GetEntry(entry);
                t2kentry = (int) (random()/((double) RAND_MAX) * nt2k);
            }while(t2kch->GetEntries() < 70000 && t2kch->GetEntries() > 30000);
            std::cout << "T2K event: " << t2khead->nevsk << std::endl;
            for(int q = 0; q < t2ktq->nhits; q++){
                (TQI->cables).push_back(t2ktq->cables[q]);
                float newtime = t2ktq->T[q] - t2ktq->T[0] + tbegin;
                (TQI->T).push_back(newtime);
                (TQI->Q).push_back(t2ktq->Q[q]);
                if (newtime > tmax) break;
                TQI->nhits++;
            }
        }
    }
    
    // output info
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
    
    // Set primary vertex
    VX = LOWE->bsvertex[0];
    VY = LOWE->bsvertex[1];
    VZ = LOWE->bsvertex[2];
    
#ifdef MERGED
    // Uncomment to merge SHE and AFT for fake data
    if (fake){
        // Set primary vertex
        VX = LOWE->bsvertex[0];
        VY = LOWE->bsvertex[1];
        VZ = LOWE->bsvertex[2];
        res.N200M = N200Max(6000., AFT_GATE+deltat - cut_window, res.T200M);
        NeutronSearch (1050, AFT_GATE + deltat - cut_window);
        theOTree->Fill();
        res.Clear();
        return;
    }
#endif
    
    if ( HEADER->idtgsk & 0x10000000 ) {// Primary SHE
    //if ( ! (HEADER->idtgsk & 0x20000000) ) {// Some events do not have SHE trigger
        // Check if previous results saved, which may happen when two
        // successive primary events exist.
        //if ( entry>0 && (! pre_saved) ) theOTree->Fill(); // Be careful with the first event
        if(entry > 0 && !pre_saved) std::cout << "badev " << HEADER->nrunsk << " " << HEADER->nevsk << std::endl;
        
        // Clear results of previous event
        res.Clear();
        
        // New results
        *head0 = *HEADER;
        *lowe0 = *LOWE;
        *mu0   = *MU;
	//        *third0 = *thirdred;
        type = 0;
        
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
        //NeutronSearch (2000, 40000.);
        //SetPrompt(0,1050);
        
        // Save AFT hit info for dark noise search
        if (!last_entry){
            if (!fake){
                // If real data look for AFT events
                ch->GetEntry(entry+1);
                GetBranchValues();
                if ( HEADER->idtgsk & 0x20000000 ){
                    for(int k = 0; k < TQI->nhits; k++){
                        if (TQI->T[k] > cut_window) break;
                        int cb = TQI->cables[k]&0xFFFF;
                        if (cb < 0 || cb > MAXPM) continue;
                        if ( ! ((TQI->cables[k]&0xFFFF0000)&0x20000 ) ) continue;
                        if ( CheckBadMis(cb) ) continue;
                        res.aftt.push_back(TQI->T[k] + she_tmax);
                        res.aftcab.push_back(cb);
                    }
                }
                ch->GetEntry(entry);
                GetBranchValues();
            }
            else{
                // Otherwise look for T2K dummy data
                int nhits = TQI->nhits; // End of the fake SHE trigger, beginning of fake AFT
                for(int k = TQI->nhits; k < t2ktq->nhits; k++){
                    int cb = t2ktq->cables[k]&0xFFFF;
                    if (cb < 0 || cb > MAXPM) continue;
                    if ( ! ((t2ktq->cables[k]&0xFFFF0000)&0x20000 ) ) continue;
                    if ( CheckBadMis(cb) ) continue;
                    float newtime = t2ktq->T[k] - t2ktq->T[0] + tbegin;
                    if (newtime > AFT_GATE + she_tmax) break;
                    res.aftt.push_back(newtime);
                    res.aftq.push_back(t2ktq->Q[k]);
                    res.aftcab.push_back(cb);
                }
            }
        }
        NeutronSearch (1050, 40000.); //dt peaks around 950 ns
        if(!fake){
            res.aftt.clear();
            res.aftcab.clear();
        }
        
        // Save SHE hit info for dark noise search
        for(int k = 0; k < TQI->nhits; k++){
            int cb = TQI->cables[k]&0xFFFF;
            if (cb < 0 || cb > MAXPM) continue;
            if ( ! ((TQI->cables[k]&0xFFFF0000)&0x20000 ) ) continue;
            if ( CheckBadMis(cb) ) continue;
            res.shet.push_back(TQI->T[k]);
            res.shecab.push_back(cb);
        }
        // Do not save, wait for AFT results.
        pre_saved = kFALSE;
        tstart = HEADER->t0;
    }
    else if ( HEADER->idtgsk & 0x20000000 ) { // AFT
        // SHE and AFT should be two successive events.
        // Some AFT events are observed to have no SHE accomponied.
        if ( HEADER->nevsk - pre_nevsk != 1 ) { 
            pre_nevsk = HEADER->nevsk;
            return;
        }
        
        // NOW THIS IS A GOOD AFT EVENT.
        if ( TQI->T[TQI->nhits - 1] < 400000 ) type = 1;      // 300 us
        else if ( TQI->T[TQI->nhits - 1] < 600000 ) type = 2; // 500 us
        else type = 3;                                        // 800 us
        nnhits = TQI->nhits ;
        n_pair ++;
        
        // Check N200 in AFT data
        Float_t t200m;
        Int_t n200m = N200Max(0., AFT_gate, t200m); 
        if ( n200m > res.N200M ) {
            res.N200M = n200m;
            res.T200M = t200m + 35000.; // shift 35 us
        }
        
        // Do 2.2 MeV search
        tend = HEADER->t0;
        deltat = (tend - tstart)/1.92;
        NeutronSearch (0., AFT_gate - cut_window, 100000. + deltat - 35000.); // shift t0 in AFT by 100us
        
        // Fill the output tree 
        if(LOWE->bsenergy > enethre && type < 3)
            theOTree->Fill();
        pre_saved = kTRUE;
        res.Clear();
    }
    else {
        std::cerr << "  Unexpected event: idtgsk = " << HEADER->idtgsk << std::endl;
        type = 4;
        pre_nevsk = HEADER->nevsk;
        std::cout << "badev2 " << HEADER->nrunsk << " " << HEADER->nevsk << std::endl;
        return;
    }
    
    pre_nevsk = HEADER->nevsk;
    
    return;
}
