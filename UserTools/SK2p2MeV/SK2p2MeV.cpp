/* vim:set expandtab tabstop=4 wrap */
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

#include "SK2p2MeV.h"

#include <iostream>

SK2p2MeV::SK2p2MeV (const Float_t (*geomxyz)[3])
{
    
    // Geometry
    xyz = geomxyz;
    
    // After trig. gate
    AFT_GATE = 500000; // ns
    
    // N10 threshold
    N10TH = 8;
    N10cutTH = 7;
    
    // Window size
    twin = 10.0;
    
    // Wiggle room for reoptimizing window
    extra = 12;
    
    // cut bad hits or not
    removeCluster = kFALSE;
    removeHighQ = kFALSE;
    removeBack = kFALSE;
    
    // Verbosity
    verbosity = 1;
    
    // Bad & missing ch.
    NBAD=NMIS=0;
    IBAD=IMIS=0;
    
    // bonsai
    bonsai_ini_();
    //bonsai_combined_ini_();
}

SK2p2MeV::~SK2p2MeV ()
{
    // cleanup
    delete TQI;
    delete LOWE;
    delete HEADER;
    if(isWIT){
        delete TQA;
        delete[] is_signal;
    }
    
    TQI = 0;
    LOWE = 0;
    HEADER = 0;
    TQA = 0;
    is_signal = 0;
    
    if(!isWIT){
       nhit = 0;
       ltimediff = 0;
    }
    
    bonsai_end_();
    //bonsai_combined_end_();
}

void SK2p2MeV::Clear(){
	res.Clear();
}

void SK2p2MeV::Print()
{
    // Show current status of SK2p2MeV
    std::cout << "\nThe After trigger gate is " << AFT_GATE << std::endl;
}

Bool_t SK2p2MeV::CheckBadMis (const Int_t cab)
{
    //std::cout << NBAD << " " << NMIS << std::endl;
    // Check if cab is bad or missing. Retrun true if yes.
    for (Int_t i=0; i<NBAD; i++) {
        if ( IBAD[i] == cab ) return kTRUE;
    }
    for (Int_t i=0; i<NMIS; i++) {
        if ( IMIS[i] == cab ) return kTRUE;
    }
    return kFALSE;
}

void SK2p2MeV::SetBadch (const Int_t nb, const Int_t *ib)
{
    // Set bad channel info
    NBAD = nb;
    IBAD = ib;
}

void SK2p2MeV::SetMisch (const Int_t nm, const Int_t *im)
{
    // Set missing channel info
    NMIS = nm;
    IMIS = im;
}

// Fill in prompt events
void SK2p2MeV::SetPrompt (Float_t tmin, Float_t tmax){
    prompt_hits candidates;
    const Int_t nhits = TQI->nhits;
    Int_t sel_hits = 0;
    for(int i = 0; i < nhits; i++){
        if (i >= MAXHITS) break;
        if ((TQI->T[i] < tmax) && (TQI->T[i] > 0)){
            //prompt.cab[sel_hits] = TQI->cables[i]&0xFFFF;
            //prompt.t[sel_hits] = TQI->T[i];
            //prompt.q[sel_hits] = TQI->Q[i];
            //prompt.t[prompt.cab[i] - 1] = TQI->T[i];
            //prompt.q[prompt.cab[i] - 1] = TQI->Q[i];
            sel_hits++;
            if (sel_hits >= MAXHITS){
                std::cout << "Error!! More than " << MAXHITS << " hits for prompt vertex!" << std::endl;
            }
        }
    }
    //prompt.nhits = sel_hits;
}

void SK2p2MeV::SetN10Threshold (const Int_t n10)
{
    // Set N10 threshold
    
    N10TH = n10;
}

void SK2p2MeV::SetN10CutThreshold (const Int_t n10)
{
    // Set N10 threshold
    
    N10cutTH = n10;
}

void SK2p2MeV::SetWindowSize (const Float_t t_win)
{
    // Set window size
    
    twin = t_win;
}

void SK2p2MeV::SetExtra (const Int_t ext)
{
    // Set wiggle room for new window
    
    extra = ext;
}

void SK2p2MeV::SetAFTGate(const Float_t gate)
{
    // Set After trigger gate
    
    AFT_GATE = gate;
}

void SK2p2MeV::SetVerbosity (const Int_t v)
{
    // verbosity = 0: output nothing
    // verbosity = 1: output necessary info
    // verbosity = 2: output details, mainly for debugging
    
    verbosity = v;
}

void SK2p2MeV::SetNcCutFlag (const Bool_t flag)
{
    removeCluster = flag;
}

void SK2p2MeV::SetHighQCutFlag (const Bool_t flag)
{
    removeHighQ = flag;
}

void SK2p2MeV::SetBackCutFlag (const Bool_t flag)
{
    removeBack = flag;
}

void SK2p2MeV::SetVertex (const Float_t x, const Float_t y, const Float_t z)
{
    VX = x;
    VY = y;
    VZ = z;
}

bool SK2p2MeV::GetBranchValues(){
    // retrieve variables from branches
    if(myTreeReader==nullptr){
        std::cout << "  No tree specified for analysis. " << std::endl;
        return false;
    }
    
    bool get_ok = true;
    
    if(isWIT){
      //get_ok &= (myTreeReader->Get("HEADER",  HEADER));
        get_ok &= (myTreeReader->Get("BONSAI",    BONSAI));
        get_ok &= (myTreeReader->Get("cable",     cable));
        get_ok &= (myTreeReader->Get("t",         t));
        get_ok &= (myTreeReader->Get("q",         q));
        get_ok &= (myTreeReader->Get("nhit",      nhit));
        get_ok &= (myTreeReader->Get("nevsk",     nevsk));
        get_ok &= (myTreeReader->Get("nrunsk",    nrunsk));
        get_ok &= (myTreeReader->Get("nsubsk",    nsubsk));
        get_ok &= (myTreeReader->Get("idtgsk",    idtgsk));
        get_ok &= (myTreeReader->Get("ltimediff", ltimediff));
        get_ok &= (myTreeReader->Get("mct",       mct));
        get_ok &= (myTreeReader->Get("trigid",    trigid));
        get_ok &= (myTreeReader->Get("bsdirks",   bsdirks));
        get_ok &= (myTreeReader->Get("swtrigger", swtrigger));
        get_ok &= (myTreeReader->Get("mc",        mc));
    } else {
        get_ok &= (myTreeReader->Get("HEADER",    HEADER));
        get_ok &= (myTreeReader->Get("TQREAL",    TQI));
        get_ok &= (myTreeReader->Get("TQAREAL",   TQA));
        get_ok &= (myTreeReader->Get("LOWE",      LOWE));;
        //get_ok &= (myTreeReader->Get("issignal",  is_signal)); //
    }
    
    return get_ok;
}


TTree* SK2p2MeV::MakeOutputTree (bool is_wit)
{
    isWIT = is_wit;
    
    theOTree = new TTree("sk2p2", "SK 2.2 MeV");
    
    // 2.2 MeV info
    theOTree->Branch("np",                 &res.np,                  "np/I");
    theOTree->Branch("N200M",              &res.N200M,               "N200M/I");
    theOTree->Branch("T200M",              &res.T200M,               "T200M/F");
    theOTree->Branch("N10",                 res.N10,                 "N10[np]/I");
    theOTree->Branch("N200",                res.N200,                "N200[np]/I");
    theOTree->Branch("N10d",                res.N10d,                "N10d[np]/I");
    theOTree->Branch("Nc",                  res.Nc,                  "Nc[np]/I");
    theOTree->Branch("Nback",               res.Nback,               "Nback[np]/I");
    theOTree->Branch("N300",                res.N300,                "N300[np]/I");
    theOTree->Branch("trms",                res.trms,                "trms[np]/F");
    theOTree->Branch("trmsdiff",            res.trmsdiff,            "trmsdiff[np]/F");
    theOTree->Branch("fpdist",              res.fpdist,              "fpdist[np]/F");
    theOTree->Branch("bpdist",              res.bpdist,              "bpdist[np]/F");
    theOTree->Branch("fwall",               res.fwall,               "fwall[np]/F");
    theOTree->Branch("bwall",               res.bwall,               "bwall[np]/F");
    theOTree->Branch("pvx",                 res.pvx,                 "pvx[np]/F");
    theOTree->Branch("pvy",                 res.pvy,                 "pvy[np]/F");
    theOTree->Branch("pvz",                 res.pvz,                 "pvz[np]/F");
    theOTree->Branch("bse",                 res.bse,                 "bse[np]/F");
    theOTree->Branch("mintrms_3",           res.mintrms_3,           "mintrms_3[np]/F");
    theOTree->Branch("mintrms_6",           res.mintrms_6,           "mintrms_6[np]/F");
    //*********new added by Yang Zhang****
    theOTree->Branch("Q10",                 res.Q10,                 "Q10[np]/F");
    theOTree->Branch("Qrms",                res.Qrms,                "Qrms[np]/F");
    theOTree->Branch("Qmean",               res.Qmean,               "Qmean[np]/F");
    theOTree->Branch("thetarms",            res.thetarms,            "thetarms[np]/F");
    theOTree->Branch("NLowtheta",           res.NLowtheta,           "NLowtheta[np]/I");
    //************************************
    theOTree->Branch("phirms",              res.phirms,              "phirms[np]/F");
    theOTree->Branch("bsdirks",             res.bsdirks,             "bsdirks[np]/F");
    theOTree->Branch("thetam",              res.thetam,              "thetam[np]/F");
    theOTree->Branch("dt",                  res.dt,                  "dt[np]/F");
    theOTree->Branch("dtn",                 res.dtn,                 "dtn[np]/F");
    theOTree->Branch("nvx",                 res.nvx,                 "nvx[np]/F");
    theOTree->Branch("nvy",                 res.nvy,                 "nvy[np]/F");
    theOTree->Branch("nvz",                 res.nvz,                 "nvz[np]/F");
    theOTree->Branch("tindex",              res.tindex,              "tindex[np]/I");
    theOTree->Branch("n40index",            res.n40index,            "n40index[np]/I");
    
    theOTree->Branch("Neff",                res.Neff,                 "Neff[np]/I");
    theOTree->Branch("ratio",               res.ratio,                "ratio[np]/F");
    theOTree->Branch("Nc1",                 res.Nc1,                  "Nc1[np]/I");
    theOTree->Branch("NhighQ",              res.NhighQ,               "NhighQ[np]/I");
    theOTree->Branch("NlowQ",               res.NlowQ,                "NlowQ[np]/I");
    
    theOTree->Branch("Nlow1",               res.Nlow1,                "Nlow1[np]/I");
    theOTree->Branch("Nlow2",               res.Nlow2,                "Nlow2[np]/I");
    theOTree->Branch("Nlow3",               res.Nlow3,                "Nlow3[np]/I");
    theOTree->Branch("Nlow4",               res.Nlow4,                "Nlow4[np]/I");
    theOTree->Branch("Nlow5",               res.Nlow5,                "Nlow5[np]/I");
    theOTree->Branch("Nlow6",               res.Nlow6,                "Nlow6[np]/I");
    theOTree->Branch("Nlow7",               res.Nlow7,                "Nlow7[np]/I");
    theOTree->Branch("Nlow8",               res.Nlow8,                "Nlow8[np]/I");
    theOTree->Branch("Nlow9",               res.Nlow9,                "Nlow9[np]/I");
    //theOTree->Branch("ncomb3", &res.ncomb3,  "ncomb3/I");
    //theOTree->Branch("theta", res.theta, "theta[ncomb3]/F");
    //theOTree->Branch("cable", res.cable, "cable[ncomb3]/I");
    /////////////////////////////////////////////////////////////////////
    //     theOTree->Branch("ncomb", &res.ncomb,  "ncomb/I");
    //     theOTree->Branch("dalpha", res.dalpha, "dalpha[ncomb]/F");
    //     theOTree->Branch("ncomb1", &res.ncomb1,  "ncomb1/I");
    //     theOTree->Branch("dalpha1", res.dalpha1, "dalpha1[ncomb1]/F");
    //     theOTree->Branch("ncomb2", &res.ncomb2,  "ncomb2/I");
    //     theOTree->Branch("dalpha2", res.dalpha2, "dalpha2[ncomb2]/F");
    /////////////////////////////////////////////////////////////////////
    
    if(not isWIT){
    theOTree->Branch("nhits",              &res.nhits,               "nhits/I");
    theOTree->Branch("goodness_combined",   res.goodness_combined,   "goodness_combined[np]/F");
    theOTree->Branch("goodness_prompt1",    res.goodness_prompt1,    "goodness_prompt1[np]/F");
    theOTree->Branch("goodness_neutron1",   res.goodness_neutron1,   "goodness_neutron1[np]/F");
    theOTree->Branch("goodness_prompt",     res.goodness_prompt,     "goodness_prompt[np]/F");
    theOTree->Branch("goodness_neutron",    res.goodness_neutron,    "goodness_neutron[np]/F");
    theOTree->Branch("goodness_window",     res.goodness_window,     "goodness_window[np]/F");
    theOTree->Branch("cvertexx",            res.cvertexx,            "cvertexx[np]/F");
    theOTree->Branch("cvertexy",            res.cvertexy,            "cvertexy[np]/F");
    theOTree->Branch("cvertexz",            res.cvertexz,            "cvertexz[np]/F");
    theOTree->Branch("pvertexx",            res.pvertexx,            "pvertexx[np]/F");
    theOTree->Branch("pvertexy",            res.pvertexy,            "pvertexy[np]/F");
    theOTree->Branch("pvertexz",            res.pvertexz,            "pvertexz[np]/F");
    theOTree->Branch("nvertexx",            res.nvertexx,            "nvertexx[np]/F");
    theOTree->Branch("nvertexy",            res.nvertexy,            "nvertexy[np]/F");
    theOTree->Branch("nvertexz",            res.nvertexz,            "nvertexz[np]/F");
    theOTree->Branch("sig_frac_peak",       res.sig_frac_peak,       "sig_frac_peak[np]/F");
    theOTree->Branch("sig_frac_tot",        res.sig_frac_tot,        "sig_frac_tot[np]/F");
    }
    
    TQI = new TQReal;
    LOWE = new LoweInfo;
    HEADER = new Header;
    if(!isWIT){
        TQA    = new TQReal;
        is_signal = new int[MAXHITS];
    }
    
    return theOTree;
}

Int_t SK2p2MeV::N200Max (Float_t tstart, Float_t tend)
{
    // Search N200 peak. Return the maximum N200.
    
    //const Int_t MAXHITS = 100000;
    Int_t   nhits;
    Int_t *cabiz = new Int_t[MAXHITS];
    Int_t *cabiz2 = new Int_t[MAXHITS];
    Float_t *tiskz = new Float_t[MAXHITS];
    Float_t *tiskz2 = new Float_t[MAXHITS];
    Int_t   *index = new Int_t[MAXHITS];
    
    nhits = 0;
    Int_t i;
    for (i=0; i<TQI->nhits; i++) {
        if (i >= MAXHITS) break;
        // Exclude non-ingate hit
        // Exclude bad ch.
        // Although bad & missing ch are masked in extract.cc, it's
        // needed here for MC
        if ( ! ((TQI->cables[i]&0xFFFF0000)&0x20000 ) ) continue;
        if ( CheckBadMis(TQI->cables[i]&0xFFFF) ) {
            //std::cout << " bad ch found: " << " " << CheckBadMis(TQI->cables[i]&0xFFFF) << std::endl;
            continue;
        }
        cabiz2[nhits] = TQI->cables[i]&0xFFFF;
        tiskz2[nhits] = TQI->T[i];
        nhits ++;
    }
    
    // Sort hits by raw time
    TMath::Sort(nhits, tiskz2, index, kFALSE); // In increasing order
    for (i=0; i<nhits; i++){
        if (i >= MAXHITS) break;
        cabiz[i] = cabiz2[ index[i] ];
        tiskz[i] = tiskz2[ index[i] ];
    }
    
    // Search N200 peak in sorted time arry tiskz2
    Int_t n200, n200max;
    n200max = 0;
    for (i=0; i<nhits; i++) {
        if(i >= MAXHITS) break;
        if ( tiskz[i] < tstart ) continue;
        if ( tiskz[i] > tend-200.) break;
        n200  = GetNhits(tiskz, i, 200., nhits);
        if ( n200 > n200max ) {
            n200max = n200;
        }
    }
    delete[] cabiz;
    delete[] cabiz2;
    delete[] tiskz;
    delete[] tiskz2;
    delete[] index;
    return n200max;
}

Int_t SK2p2MeV::N200Max (Float_t tstart, Float_t tend, Float_t & t200m)
{
    // Search N200 peak. Return the maximum N200 and timing of the peak.
    
    const Int_t MAXHITS = 100000;
    Int_t   nhits;
    Int_t *cabiz = new Int_t[MAXHITS];
    Int_t *cabiz2 = new Int_t[MAXHITS];
    Float_t *tiskz = new Float_t[MAXHITS];
    Float_t *tiskz2 = new Float_t[MAXHITS];
    Int_t   *index = new Int_t[MAXHITS];
    
    nhits = 0;
    Int_t i;
    for (i=0; i<TQI->nhits; i++) {
        if (i >= MAXHITS) break;
        // Exclude non-ingate hit
        // Exclude bad ch.
        // Although bad & missing ch are masked in extract.cc, it's
        // needed here for MC
        if ( ! ((TQI->cables[i]&0xFFFF0000)&0x20000 ) ) continue;
        if ( CheckBadMis(TQI->cables[i]&0xFFFF) ) {
            //std::cout << " bad ch found: " << std::endl;
            continue;
        }
        cabiz2[nhits] = TQI->cables[i]&0xFFFF;
        tiskz2[nhits] = TQI->T[i];
        nhits ++;
    }
    
    // Sort hits by raw time
    TMath::Sort(nhits, tiskz2, index, kFALSE); // In increasing order
    for (i=0; i<nhits; i++){
        if (i >= MAXHITS) break;
        cabiz[i] = cabiz2[ index[i] ];
        tiskz[i] = tiskz2[ index[i] ];
    }
    
    // Search N200 peak in sorted time arry tiskz2
    Int_t n200, n200max;
    n200max = 0;
    for (i=0; i<nhits; i++) {
        if (i >= MAXHITS) break;
        if ( tiskz[i] < tstart ) continue;
        if ( tiskz[i] > tend-200.) break;
        n200  = GetNhits(tiskz, i, 200., nhits);
        if ( n200 > n200max ) {
            n200max = n200;
            t200m = tiskz[i] + 100.; // output!!
        }
    }
    delete[] cabiz;
    delete[] cabiz2;
    delete[] tiskz;
    delete[] tiskz2;
    delete[] index;
    return n200max;
}

void SK2p2MeV::NeutronSearch (Float_t tstart, Float_t tend, Float_t TOFFSET)
{
    // Search N10 peak and cal N10, Neff, etc., for each peak.
    
    // Calculate weight for each PMT
    const Int_t MAXPMT = 11146;
    Float_t wt[MAXPMT];
    for (Int_t i=0; i<MAXPMT; i++) {
        Float_t v[3];
        v[0] = VX;
        v[1] = VY;
        v[2] = VZ;
        try {
          wt[i] = GetWeight(xyz[i], v);
        }
        catch (const std::exception& e){
          std::cerr << "SK2p2MeV::GetWeight exception!"<<std::endl;
          std::cout << e.what() << "\n";
          std::cout << "pmt: " << i << "\n";
          for (int j = -5; j < 5; ++j){ 
            std::cout << "xyz[pmt][0]: " << xyz[i+j][0] << " ";
            std::cout << "xyz[pmt][1]: " << xyz[i+j][1] << " ";
            std::cout << "xyz[pmt][1]: " << xyz[i+j][2] << "\n\n";
          }
          exit(0);
        }
    }
    Float_t maxwt = TMath::MaxElement(MAXPMT, wt);
    Float_t tot_wt = 0.;
    for (Int_t i=0; i<MAXPMT; i++){
        wt[i] = wt[i] / maxwt; //calculate relative wt
        tot_wt += wt[i];
    }
    
    //
    //const Int_t MAXHITS = 100000;
    // FIXME yo this is crazy inefficient, we should not be invoking new for every call
    Int_t   nhits;
    Int_t *cabiz = new Int_t[MAXHITS];
    Int_t *cabiz2 = new Int_t[MAXHITS];
    Int_t *cabiz3 = new Int_t[MAXHITS];
    Float_t *tiskz = new Float_t[MAXHITS];
    Float_t *tiskz2 = new Float_t[MAXHITS];
    Float_t *tiskz3 = new Float_t[MAXHITS];
    Float_t *qiskz = new Float_t[MAXHITS];
    Float_t *qiskz2 = new Float_t[MAXHITS];
    Int_t *index = new Int_t[MAXHITS];
    Int_t *is_signal2 = new Int_t[MAXHITS];
    Int_t nindex[MAXN10];
    Float_t *hitv_x = new Float_t[MAXHITS]; //hit vector
    Float_t *hitv_y = new Float_t[MAXHITS];
    Float_t *hitv_z = new Float_t[MAXHITS];
    Int_t dark_flag[MAXHITS]={0}, dark_flag0[MAXHITS]={0};
    Int_t ndark=0;
    
    nhits = 0;
    Int_t i;
    res.nhits = TQI->nhits;
    for (i=0; i<TQI->nhits; i++) {  //nhits in defined to be the # of total hits(not only in 1.3us in extract)
        if (i > MAXHITS) break;
        // Exclude non-ingate hit
        // Exclude bad ch.
        // Although bad & missing ch are masked in extract.cc, it's needed here for MC
        int cb = TQI->cables[i]&0xFFFF;
        if (cb < 0 || cb > MAXPM) continue;
        if ( ! ((TQI->cables[i]&0xFFFF0000)&0x20000 ) ) continue;
        if ( CheckBadMis(TQI->cables[i]&0xFFFF) ) {
            //std::cout << " bad ch found: " << " " << CheckBadMis(TQI->cables[i]&0xFFFF) << " " << TQI->nhits << std::endl;
            //std::cout << " bad ch found: " << std::endl;
            continue;
        }
        cabiz2[nhits] = TQI->cables[i]&0xFFFF;
        cabiz3[nhits] = TQI->cables[i];
        tiskz2[nhits] = TQI->T[i];
        tiskz3[nhits] = TQI->T[i];
        qiskz2[nhits] = TQI->Q[i];
        is_signal2[nhits] = is_signal[i];
        
        // Look for continuous dark noise
        // this scans through preceding SHE hits, looking for another hit on the same PMT within a given time window.
        // The scan is performed with two window sizes, setting the dark_flag0 to 1 if a hit is found in the first window,
        // then setting (potentially overriding) it to 2 if a hit is found in the second window.
        int tmpl=nhits;
        // First, for AFT in data look in previous SHE
        if (res.shet.size() > 0 && res.shecab.size() > 0 && tiskz2[nhits] < cut_window){
            tmpl = res.shet.size() - 1;
            while(TMath::Abs((tiskz2[nhits]+she_tmax-res.shet[tmpl]))<cut_window && tmpl>=0 && res.shet[tmpl] > primary_window){
                if(cabiz2[nhits]==res.shecab[tmpl]){
                    dark_flag0[nhits]=1;
                    break;
                }
                tmpl--;
            }
            tmpl = res.shet.size() - 1;
            while(TMath::Abs((tiskz2[nhits]+she_tmax-res.shet[tmpl]))<cut_window2 && tmpl>=0 && res.shet[tmpl] > primary_window){
                if(cabiz2[nhits]==res.shecab[tmpl]){
                    dark_flag0[nhits]=2;
                    break;
                }
                tmpl--;
            }
        }
        
        // Now perform regular search (hits in the same trigger window)
        tmpl = nhits;
        while(TMath::Abs((tiskz2[nhits]-tiskz2[tmpl]))<cut_window && tmpl>=0){
            if(cabiz2[nhits]==cabiz2[tmpl]&&tmpl!=nhits&&(TMath::Abs(tiskz2[tmpl])>primary_window&&TMath::Abs(tiskz2[nhits])>primary_window)){
                dark_flag0[tmpl]=1;
                dark_flag0[nhits]=1;
                break;
            }
            tmpl--;
        }
        
        tmpl=nhits;
        while(TMath::Abs((tiskz2[nhits]-tiskz2[tmpl]))<cut_window2 && tmpl>=0){
            if(cabiz2[nhits]==cabiz2[tmpl]&&tmpl!=nhits&&(TMath::Abs(tiskz2[tmpl])>primary_window&&TMath::Abs(tiskz2[nhits])>primary_window)){
                dark_flag0[tmpl]=2;
                dark_flag0[nhits]=2;
                break;
            }
            tmpl--;
        }
        
        // Finally, for SHE in data look in following AFT
        if (res.aftt.size() > 0 && res.aftcab.size() > 0 && tiskz2[nhits] > she_tmax - cut_window){
            tmpl = 0;
            while(TMath::Abs((tiskz2[nhits]-res.aftt[tmpl]))<cut_window && tmpl < res.aftt.size()){
                if(cabiz2[nhits]==res.aftcab[tmpl]){
                    dark_flag0[nhits]=1;
                    break;
                }
                tmpl++;
            }
            tmpl = 0;
            while(TMath::Abs((tiskz2[nhits]-res.aftt[tmpl]))<cut_window2 && tmpl < res.aftt.size()){
                if(cabiz2[nhits]==res.aftcab[tmpl]){
                    dark_flag0[nhits]=2;
                    break;
                }
                tmpl++;
            }
        }
        nhits++;
    }
    if(verbosity > 1) std::cout<<"TOTAL: "<<nhits<<std::endl;
    res.nsignal = 0;
    for (i=0; i<nhits; i++) {
        if (tiskz3[i] > 19000)
            res.nsignal += is_signal2[i];
    }
    if(verbosity > 1) std::cout << "nsignal is " << res.nsignal << std::endl;
    
    // TOF
    for (i=0; i<nhits; i++) {
        if (i > MAXHITS) break;
        Float_t tof;
        tof = TMath::Sqrt((VX - xyz[cabiz2[i]-1][0]) * (VX - xyz[cabiz2[i]-1][0])
                +(VY - xyz[cabiz2[i]-1][1]) * (VY - xyz[cabiz2[i]-1][1])
                +(VZ - xyz[cabiz2[i]-1][2]) * (VZ - xyz[cabiz2[i]-1][2])) / C_WATER;
        tiskz2[i] -= tof;
    }
    
    // Sort hits by TOF-corrected time
    TMath::Sort(nhits, tiskz2, index, kFALSE); // In increasing order
    for (i=0; i<nhits; i++){
        if (i > MAXHITS) break;
        cabiz[i] = cabiz2[ index[i] ];
        tiskz[i] = tiskz2[ index[i] ];
        qiskz[i] = qiskz2[ index[i] ];
        is_signal[i] = is_signal2[index[i]];
        dark_flag[i]=dark_flag0[index[i]];
    }
    
    // Calculate hit vectors
    for (i=0; i<nhits; i++) {
        if (i > MAXHITS) break;
        Float_t pmt_r;
        pmt_r = TMath::Sqrt((VX - xyz[cabiz[i]-1][0]) * (VX - xyz[cabiz[i]-1][0])
                +(VY - xyz[cabiz[i]-1][1]) * (VY - xyz[cabiz[i]-1][1])
                +(VZ - xyz[cabiz[i]-1][2]) * (VZ - xyz[cabiz[i]-1][2]));
        hitv_x[i] = (xyz[cabiz[i]-1][0] - VX)/pmt_r;
        hitv_y[i] = (xyz[cabiz[i]-1][1] - VY)/pmt_r;
        hitv_z[i] = (xyz[cabiz[i]-1][2] - VZ)/pmt_r;
    }
    
    // Use a 10 ns window to search 2.2MeV candidate
    Float_t uvx[MAXN10], uvy[MAXN10], uvz[MAXN10];
    Float_t qi[MAXN10], ti[MAXN10];
    Float_t qn[MAXPM], tn[MAXPM];
    Float_t qp[MAXPM], tp[MAXPM];
    Int_t   ci[MAXN10]; //cable
    Int_t   cn[MAXPM]; //cable
    Int_t   cp[MAXPM]; //cable
    Int_t   uvf[MAXN10]; //flag: 0=not cut, 1=first cut, 2=second cut, etc.
    Int_t   ncut;
    float sig_frac_peak, sig_frac_tot;
    
    Bool_t   pre_t0_set = kFALSE;
    Float_t  t0 = -1, t0n=0;
    Int_t    N10i, N10 = 0, N10in, N10d, N10n=0, N40index=0;
    Int_t pN10=0, nN10=0;
    Int_t n40indexlow = N40index, n40indexhigh, N10send;
    Float_t bsenergy, bsvertexx, bsvertexy, bsvertexz, bsgood;
    Float_t bsgood_combined[3];
    Float_t qisend[MAXPM], tisend[MAXPM];
    Int_t   cisend[MAXPM*2];
    Int_t    Nc, Nback, Nlow[9], Neff, Nc1, NhighQ, NlowQ, n200;
    Float_t  mintrms_3,mintrms_4,mintrms_5,mintrms_6;
    Float_t  trms, newtrms, fpdist, fwall, trmsdiff, phirms, bsdirks, thetam, ratio;
    Float_t  bpdist, bwall, bse, bvr, bvz;
    Float_t  fvr, fvz;
    //*********new added by Yang Zhang*****
    Float_t Q10, Qmean, Qrms, thetarms;
    Float_t nvx, nvy, nvz;
    Float_t tnvx, tnvy, tnvz;
    Float_t cvertexx, cvertexy, cvertexz;
    Float_t pvertexx, pvertexy, pvertexz;
    Float_t nvertexx, nvertexy, nvertexz;
    Float_t goodness_combined, goodness_prompt, goodness_neutron;
    Float_t goodness_prompt1, goodness_neutron1, goodness_window;
    Int_t    NLowtheta;
    Int_t  tindex=0, n40hits=0;
    //************************************
    for ( i=0; i<nhits; i++) {
        if (i > MAXHITS) break;
        if ( tiskz[i] < tstart ) continue;
        if ( tiskz[i] > tend-twin ) continue;
        
        // EXCLUDE LE, HE or SHE PEAK REGION
        //if ( tiskz[i]>500. && tiskz[i]<1500.) continue;
        
        // Calculate hits in 10 ns window
        N10i = GetNhits(tiskz, i, twin, nhits);
        pN10=N10i;
        int darkcut_flag;
        if(N10i>N10cutTH){
            darkcut_flag=2;
            if(dark_flag[i]==2)continue;
        }
        else{
            darkcut_flag=1;
            if(dark_flag[i]!=0)continue;
        }
        int N10iold = N10i;
        N10i=GetNhits_flag(tiskz,dark_flag,darkcut_flag,i,10.,nhits);
        //if (N10iold != N10i) std::cout << "ndiff " << N10iold << " " << N10i << std::endl;
        
        // Only consider candidates with N10 >= N10TH && N10 <= 50
        if ((N10i < N10TH) || (N10i > 50)) continue;
        
        // Store info of previous peak when encountering a new peak.
        // We regard time diff of t0 of two N10 peaks greater than
        // 20 ns as two seperated peaks.
        // dt is t0 of previous N10 peak
        if ( pre_t0_set && (tiskz[i] - t0 > twin * 2) ) {
            // Save result to array
            res.N10[res.np]    = N10;
            res.N10d[res.np]    = N10d;
            res.goodness_combined[res.np] = goodness_combined;
            res.goodness_prompt1[res.np] = goodness_prompt1;
            res.goodness_prompt[res.np] = goodness_prompt;
            res.goodness_neutron1[res.np] = goodness_neutron1;
            res.goodness_neutron[res.np] = goodness_neutron;
            res.goodness_window[res.np] = goodness_window;
            res.Nc[res.np]     = Nc;
            res.Nback[res.np]  = Nback;
            res.N300[res.np]   = GetNXX(nhits, tiskz, 300., t0+twin/2.);
            res.N200[res.np]   = n200;
            res.trms[res.np]   = trms;
            res.trmsdiff[res.np]   = trmsdiff;
            res.fpdist[res.np]   = fpdist;
            res.bpdist[res.np]   = bpdist;
            res.fwall[res.np]   = fwall;
            res.bwall[res.np]   = bwall;
            res.pvx[res.np]   = VX;
            res.pvy[res.np]   = VY;
            res.pvz[res.np]   = VZ;
            res.cvertexx[res.np]   = cvertexx;
            res.cvertexy[res.np]   = cvertexy;
            res.cvertexz[res.np]   = cvertexz;
            res.pvertexx[res.np]   = pvertexx;
            res.pvertexy[res.np]   = pvertexy;
            res.pvertexz[res.np]   = pvertexz;
            res.nvertexx[res.np]   = nvertexx;
            res.nvertexy[res.np]   = nvertexy;
            res.nvertexz[res.np]   = nvertexz;
            res.bse[res.np]   = bse;
            res.mintrms_6[res.np]   = mintrms_6;
            res.mintrms_3[res.np]   = mintrms_3;
            res.nvx[res.np] = nvx;
            res.nvy[res.np] = nvy;
            res.nvz[res.np] = nvz;
            //*********new added by Yang Zhang*****
            res.Q10[res.np]   = Q10;
            res.Qrms[res.np]   = Qrms;
            res.Qmean[res.np]   = Qmean;
            res.thetarms[res.np] = thetarms;
            res.NLowtheta[res.np] = NLowtheta;
            //************************************
            res.phirms[res.np] = phirms;
            res.bsdirks[res.np] = bsdirks;
            res.thetam[res.np] = thetam;
            res.dt[res.np]     = t0 + TOFFSET; // TOFFSET only useful in SHE + AFT case
            res.dtn[res.np]    = t0n+TOFFSET;
            
            res.Neff[res.np]   = Neff;
            res.ratio[res.np]  = ratio;
            res.Nc1[res.np]    = Nc1;
            res.NhighQ[res.np] = NhighQ;
            res.NlowQ[res.np] = NlowQ;
            
            res.Nlow1[res.np]  = Nlow[0];
            res.Nlow2[res.np]  = Nlow[1];
            res.Nlow3[res.np]  = Nlow[2];
            res.Nlow4[res.np]  = Nlow[3];
            res.Nlow5[res.np]  = Nlow[4];
            res.Nlow6[res.np]  = Nlow[5];
            res.Nlow7[res.np]  = Nlow[6];
            res.Nlow8[res.np]  = Nlow[7];
            res.Nlow9[res.np]  = Nlow[8];
            res.sig_frac_peak[res.np] = sig_frac_peak;
            res.sig_frac_tot[res.np] = sig_frac_tot;
            
            res.np ++;
            if ( res.np >= MAXNP - 1 ) break; // This should not happen. If it does happen, it's
            // a bad event.
            //***********1**added by Yang Zhang, for experimental test of theta
            //Float_t dir_tmp[3];
            //dir_tmp[0] = dir_tmp[1] = dir_tmp[2] = 0.;
            //for (Int_t i=0; i<N10; i++) {
                //dir_tmp[0] += uvx[i];
                //dir_tmp[1] += uvy[i];
                //dir_tmp[2] += uvz[i];
            //}
            //Float_t vr_tmp = sqrt(dir_tmp[0]*dir_tmp[0] + dir_tmp[1]*dir_tmp[1]
                    //+ dir_tmp[2]*dir_tmp[2]);
            //dir_tmp[0] = dir_tmp[0] / vr_tmp;
            //dir_tmp[1] = dir_tmp[1] / vr_tmp;
            //dir_tmp[2] = dir_tmp[2] / vr_tmp;
            //Float_t the;
            //for (Int_t i=0; i<N10; i++) {
                //the = dir_tmp[0]*uvx[i] + dir_tmp[1]*uvy[i] + dir_tmp[2]*uvz[i];
                //the = TMath::ACos(the)*180./TMath::Pi();
                //// all hits
                //res.theta[res.ncomb3] = the;
                ////cable
                //res.cable[res.ncomb3] = ci[i];
                //res.ncomb3 ++;
                //if ( res.ncomb3 > 10000 ) {
                    //std::cout << "ncomb3: " << res.ncomb3 << std::endl;
                    //break;
                //}
            //}
            N10 = 0; //after a peak is stored, re-initialize N10
        }
        //std::cout << std::endl << "Candidate total = " << res.np << std::endl << std::endl;
        
        if ( N10i <= N10 ) continue;
        //if(is_signal2[i] > 0 && tiskz[i] > 19000 && N10i > 6) std::cout << "sssig " << N10i << " " << tiskz[i] << " " << tiskz3[i] << " " << i << std::endl;
        
        // Now this is a "better" peak.
        // Calculate discriminating variables for the new peak.
        N10 = N10i;
        tindex = i;
        t0 = tiskz[i];
        ndark=0;
        n200 = GetNXX(TQI->nhits, tiskz, 200., t0+twin/2.);
        int totsig = 0;
        pre_t0_set = kTRUE;
        for (Int_t j=0; j<N10; j++) {
            if(dark_flag[i+j+ndark]==darkcut_flag){
                ndark++;
                j--;
                continue;
            }
            uvx[j] = hitv_x[i+j+ndark]; // hit vector in 10 ns
            uvy[j] = hitv_y[i+j+ndark];
            uvz[j] = hitv_z[i+j+ndark];
            uvf[j] = 0; // NOT CUT YET
            ti[j] =  tiskz[i+j+ndark];
            qi[j]  = qiskz[i+j+ndark];
            ci[j]  = cabiz[i+j+ndark];
            totsig += is_signal[i+j+ndark];
        }
        sig_frac_peak = totsig * 1.0/N10;
        //std::cout << "ssig " << totsig << " " << N10 << " " << sig_frac_peak << std::endl;
        sig_frac_tot = res.nsignal > 0 ? totsig * 1./res.nsignal : 0;
        
        // trms
        trms = GetTrms (ti, uvf, N10);
        //*********new added by Yang Zhang*****
        
        ndark=0;
        for (Int_t j=0; j<N10; j++) {
            if(dark_flag0[index[i+j+ndark]]==darkcut_flag){
                ndark++;
                j--;
                continue;
            }
            //std::cout<<"Reset t "<<j<<" "<<ti[j]<<" "<<ci[j]<<" "<<cabiz2[index[i+j]]<<std::endl;
            ti[j]  = tiskz3[index[i+j+ndark]];
            //TQI->T[index[tindex+j]];
            qi[j]  = qiskz2[index[i+j+ndark]];
            ci[j]  = cabiz2[index[i+j+ndark]];
            tisend[ci[j]-1] = ti[j];
            qisend[ci[j]-1] = qi[j];
            //qi[j]  = TQI->Q[index[tindex+j]];
            //ci[j]  = TQI->cables[index[tindex+j]]&0xFFFF;
            //std::cout<<"Final t "<<j<<" "<<ti[j]<<" "<<ci[j]<<std::endl;
        }
        // Run Bonsai for hits in the small window, extract goodness of fit
        //bonsai_fit_(&t0, tisend, qisend, ci, &N10, &bsenergy, &bsvertexx, &bsvertexy, &bsvertexz, &bsgood);
        //goodness_window = bsgood;
        
        //now look for minimum trms point, requiring it within 200cm of pvx.
        //this still allows for good improvement of neutron signal, but background is not affected so much.
        //cannot use original ti ci because ti has to been w/o TOF here
        Float_t newtrms = MinimizeTrms(ti, ci, 0, nindex , N10, nvx, nvy, nvz, VX, VY, VZ, 10000);
        ndark=0;
        for (Int_t j=0; j<N10; j++) {
            if (i+j+ndark >= nhits || i+j+ndark < 0) continue;
            if(dark_flag0[index[i+j+ndark]]==darkcut_flag){
                ndark++;
                j--;
                continue;
            }
            ti[j]  = tiskz3[index[i+j+ndark]];
            qi[j]  = qiskz2[index[i+j+ndark]];
            ci[j]  = cabiz2[index[i+j+ndark]];
        }
        newtrms = MinimizeTrms(ti, ci, 0, nindex , N10, tnvx, tnvy, tnvz, VX, VY, VZ, 200);
        
        Int_t nhitaim = N10+extra;
        n40hits = 0;
        Int_t low = extra/2;
        
        for (Int_t j=0; j<nhitaim; j++)
        {
            if (i+j-low >= nhits || i+j-low < 0) continue;
            if(dark_flag0[index[i+j-low]]==darkcut_flag){
                continue;
            }
            ti[n40hits]  = tiskz3[index[i+j-low]];
            qi[n40hits]  = qiskz2[index[i+j-low]];
            ci[n40hits]  = cabiz2[index[i+j-low]];
            n40hits++;
        }
        //Time of flight correct these 20 hits to the new vertex, found by minimizing trms
        BasicTof(ti,ci,tnvx, tnvy, tnvz, n40hits,nindex);
        N10n = 0;
        
        //search for a new best N10 (N10n) from these new ToF corrected hits
        for (Int_t j = 0; j < n40hits; j++)
        {
            //if(dark_flag[i+j]==darkcut_flag)continue;
            N10in = GetNhits_flag (ti, dark_flag, 100, j, 10, n40hits);
            if (N10in > N10n)
            {
                N10n = N10in;
                N40index = j;
                t0n = ti[j];
            }
        }
        //std::cout<<t0n<<" "<<N10i<<std::endl;
        
        for (Int_t j=0; j<n40hits; j++) {
            Float_t pmt_r;
            pmt_r = TMath::Sqrt((VX - xyz[ci[j]-1][0]) * (VX - xyz[ci[j]-1][0])
                    +(VY - xyz[ci[j]-1][1]) * (VY - xyz[ci[j]-1][1])
                    +(VZ - xyz[ci[j]-1][2]) * (VZ - xyz[ci[j]-1][2]));
            hitv_x[j] = (xyz[ci[j]-1][0] - VX)/pmt_r;
            hitv_y[j] = (xyz[ci[j]-1][1] - VY)/pmt_r;
            hitv_z[j] = (xyz[ci[j]-1][2] - VZ)/pmt_r;
            
            if(hitv_x[j]!=hitv_x[j])
            {
                //nan check for events with vertex at pmt location
                hitv_x[j]=0;
                hitv_y[j]=0;
                hitv_z[j]=0;
            }
        }
        for (Int_t j=0; j<N10n; j++) {
            uvx[j] = hitv_x[N40index+j];
            uvy[j] = hitv_y[N40index+j];
            uvz[j] = hitv_z[N40index+j];
            uvf[j] = 0;
            ti[j] =  ti[N40index+j];
            qi[j]  = qi[N40index+j];
            ci[j]  = ci[N40index+j];
        }
        // Check low hits
        for (Int_t j=0; j<9; j++) {
            Nlow[j] = GetLowHits (ci, N10n, wt, 0.5+j*0.05); //last number is accep
        }
        // cluster search
        ncut = 0;
        Nc = ncut;
        GetCluster(uvx, uvy, uvz, uvf, N10n, ncut, 3, 0.97); // 14.1 degree
        Nc = ncut - Nc;
        /////////////////////////////////////////////////////
        // do not cut hits
        if (!removeCluster){
            for (Int_t j=0; j<N10n; j++) uvf[j] = 0;
            ncut = 0;
        }
        /////////////////////////////////////////////////////
        // check high Q
        NhighQ = ncut;
        CheckHighQ (qi, uvf, N10n, ncut, 3.); // 5 is qth
        NhighQ = ncut - NhighQ;
        /////////////////////////////////////////////////////
        // do not cut hits
        //Linyan
        if (!removeHighQ){
            for (Int_t j=0; j<N10; j++) uvf[j] = 0;
            ncut = 0;
        }
        /////////////////////////////////////////////////////
        // check backward hits
        Nback = ncut;
        CheckBackHits (uvx, uvy, uvz, uvf, N10n, ncut, 90.);
        Nback = ncut - Nback;
        if (!removeBack){
            for (Int_t j=0; j<N10; j++) uvf[j] = 0;
            ncut = 0;
        }
        
        // Neff, equal to N10-Nback
        Neff = N10n - ncut;
        // ratio
        ratio = GetRatio (ci, wt, uvf, N10n);
        //Q10
        Q10 = GetQ10 (qi, uvf, N10n);
        //Qrms
        Qrms = GetQrms (qi, uvf, N10n);
        //Qmean
        Qmean = GetQmean (qi, uvf, N10n);
        //thetarms
        thetarms = GetThetaRms (uvx, uvy, uvz, uvf, N10n);
        //NlowTheta
        NLowtheta = GetNLowTheta (uvx, uvy, uvz, uvf, N10n, 20.);
        //************************************
        // phirms
        phirms = GetPhiRms (uvx, uvy, uvz, uvf, N10n);
        // dirks
        bsdirks = GetDirKS (uvx, uvy, uvz, uvf, N10n);
        // theta mean
        thetam = GetThetaMean (uvx, uvy, uvz, uvf, N10n);
        
        // min-trms (different from minimizing trms we did earlier)
        GetMinTrms (ti, N10n, mintrms_6, mintrms_5, mintrms_4, mintrms_3, uvf);
        newtrms=GetTrms(ti, uvf, N10n);
        trmsdiff=trms-newtrms;
        N10d=N10n-N10;
        fpdist=TMath::Sqrt(TMath::Power(nvx-VX,2)+TMath::Power(nvy-VY,2)+TMath::Power(nvz-VZ,2));
        fvr=RINTK - TMath::Sqrt(nvx*nvx+nvy*nvy);
        fvz = ZPINTK - TMath::Abs(nvz);
        fwall=(fvr < fvz) ? fvr : fvz;
        
        N10send = 0;
        bool out_pmt = false;
        for (Int_t j = 0; j < nhits; j++)
        {
            if ((tiskz3[j] > t0 - 520.8) && (tiskz3[j] < t0 + 779.2))
            {
                if (cabiz2[j]>MAXPMT) out_pmt = true;
                if (N10send >= MAXPM) continue;
                tisend[cabiz2[j]-1] = tiskz3[j];
                qisend[cabiz2[j]-1] = qiskz2[j];
                cisend[N10send] = cabiz2[j];
                tn[N10send] = tiskz3[j];
                qn[N10send] = qiskz2[j];
                cn[N10send] = cabiz2[j];
                N10send++;
            }
        }
        if(out_pmt) continue;
        if (N10send< 10 || N10send>1000)
        {
            std::cout << "N1300 out of range = " << N10send << std::endl;
            std::cout << "BS1300 assigned to truncated BS values..." << std::endl;
            bsenergy = 0;
            bsvertexx = nvx;
            bsvertexy = nvy;
            bsvertexz = nvz;
            bsgood = 0;
        }
        else
        {
            bsvertexx = 1e5;
            bsvertexy = 1e5;
            bsvertexz = 1e5;
            bonsai_fit_(&t0, tisend, qisend, cisend, &N10send, &bsenergy, &bsvertexx, &bsvertexy, &bsvertexz, &bsgood);
            goodness_neutron = bsgood;
        }
        bpdist=TMath::Sqrt(TMath::Power(bsvertexx-nvx,2)+TMath::Power(bsvertexy-nvy,2)+TMath::Power(bsvertexz-nvz,2));
        bvr=RINTK - TMath::Sqrt(bsvertexx*bsvertexx+bsvertexy*bsvertexy);
        bvz = ZPINTK - TMath::Abs(bsvertexz);
        bwall=(bvr < bvz) ? bvr : bvz;
        bse=bsenergy;
    }
    
    // Fill the last peak. IMPORTANT!!!
    if ( N10 > 0) { //then N10 is stored by N10i, and N10i has passed the N10TH
        //    if ( N10 > 0 && n200 < N200CUT) { //then N10 is stored by N10i, and N10i has passed the N10TH
        // Save result to array
        res.nvx[res.np] = nvx;
        res.nvy[res.np] = nvy;
        res.nvz[res.np] = nvz;
        res.goodness_combined[res.np] = goodness_combined;
        res.goodness_prompt1[res.np] = goodness_prompt1;
        res.goodness_prompt[res.np] = goodness_prompt;
        res.goodness_neutron1[res.np] = goodness_neutron1;
        res.goodness_neutron[res.np] = goodness_neutron;
        res.goodness_window[res.np] = goodness_window;
        res.N10[res.np]    = N10;
        res.N10d[res.np]    = N10d;
        res.Nc[res.np]     = Nc;
        res.Nback[res.np]  = Nback;
        res.N300[res.np]   = GetNXX(nhits, tiskz, 300., t0+twin/2.);
        res.trms[res.np]   = trms;
        res.trmsdiff[res.np]   = trmsdiff;
        res.fpdist[res.np]   = fpdist;
        res.bpdist[res.np]   = bpdist;
        res.fwall[res.np]   = fwall;
        res.bwall[res.np]   = bwall;
        res.pvx[res.np]   = VX;
        res.pvy[res.np]   = VY;
        res.pvz[res.np]   = VZ;
        res.cvertexx[res.np]   = cvertexx;
        res.cvertexy[res.np]   = cvertexy;
        res.cvertexz[res.np]   = cvertexz;
        res.pvertexx[res.np]   = pvertexx;
        res.pvertexy[res.np]   = pvertexy;
        res.pvertexz[res.np]   = pvertexz;
        res.nvertexx[res.np]   = nvertexx;
        res.nvertexy[res.np]   = nvertexy;
        res.nvertexz[res.np]   = nvertexz;
        res.bse[res.np]   = bse;
        res.mintrms_6[res.np]   = mintrms_6;
        res.mintrms_3[res.np]   = mintrms_3;
        res.Q10[res.np]   = Q10;
        res.Qrms[res.np]   = Qrms;
        res.Qmean[res.np]   = Qmean;
        res.thetarms[res.np]   = thetarms;
        res.NLowtheta[res.np]   = NLowtheta;
        res.phirms[res.np] = phirms;
        res.bsdirks[res.np] = bsdirks;
        res.thetam[res.np] = thetam;
        res.dt[res.np]     = t0 + TOFFSET; // TOFFSET only useful in SHE + AFT case
        res.dtn[res.np]     = t0n + TOFFSET; // TOFFSET only useful in SHE + AFT case
        res.N200[res.np]   = n200;
        res.n40index[res.np] = N40index;
        res.tindex[res.np] = tindex;
        
        res.Neff[res.np]   = Neff;
        res.ratio[res.np]  = ratio;
        res.Nc1[res.np]    = Nc1;
        res.NhighQ[res.np] = NhighQ;
        res.NlowQ[res.np] = NlowQ;
        
        res.Nlow1[res.np]  = Nlow[0];
        res.Nlow2[res.np]  = Nlow[1];
        res.Nlow3[res.np]  = Nlow[2];
        res.Nlow4[res.np]  = Nlow[3];
        res.Nlow5[res.np]  = Nlow[4];
        res.Nlow6[res.np]  = Nlow[5];
        res.Nlow7[res.np]  = Nlow[6];
        res.Nlow8[res.np]  = Nlow[7];
        res.Nlow9[res.np]  = Nlow[8];
        res.sig_frac_peak[res.np] = sig_frac_peak;
        res.sig_frac_tot[res.np] = sig_frac_tot;
        
        res.np ++;
        //***********added by Yang Zhang, for experimental test of theta
        //Float_t dir_tmp[3];
        //dir_tmp[0] = dir_tmp[1] = dir_tmp[2] = 0.;
        //for (Int_t i=0; i<N10; i++) {
            //dir_tmp[0] += uvx[i];
            //dir_tmp[1] += uvy[i];
            //dir_tmp[2] += uvz[i];
        //}
        //Float_t vr_tmp = sqrt(dir_tmp[0]*dir_tmp[0] + dir_tmp[1]*dir_tmp[1]
                //+ dir_tmp[2]*dir_tmp[2]);
        //dir_tmp[0] = dir_tmp[0] / vr_tmp;
        //dir_tmp[1] = dir_tmp[1] / vr_tmp;
        //dir_tmp[2] = dir_tmp[2] / vr_tmp;
        //Float_t the;
        //for (Int_t i=0; i<N10; i++) {
            //the = dir_tmp[0]*uvx[i] + dir_tmp[1]*uvy[i] + dir_tmp[2]*uvz[i];
            //the = TMath::ACos(the)*180./TMath::Pi();
            //// all hits
            //res.theta[res.ncomb3] = the;
            ////cable
            //res.cable[res.ncomb3] = ci[i];
            //res.ncomb3 ++;
            //if ( res.ncomb3 > 10000 ) {
                //std::cout << "ncomb3: " << res.ncomb3 << std::endl;
                //break;
            //}
        //}
        N10 = 0;
    }
    delete[] cabiz ;
    delete[] cabiz2;
    delete[] cabiz3;
    delete[] tiskz ;
    delete[] tiskz2;
    delete[] tiskz3;
    delete[] is_signal2;
    delete[] qiskz ;
    delete[] qiskz2;
    delete[] index ;
    delete[] hitv_x; //hit vector
    delete[] hitv_y;
    delete[] hitv_z;
}

Int_t SK2p2MeV::GetNhits_flag(Float_t *v, Int_t *flag, Int_t flagcut, Int_t start_index, Float_t width, Int_t nhits)
    // Get nhits from v[0] in time window = width
    // v[] - hit timing, sorted in time
    // width - time window width
    // start_index - index of start pmt hit
    // nhits - total number  of v[] hits
{
    Int_t i = start_index;
    Int_t dark=0;
    while (1) {
        if(flag[i]>=flagcut)dark++;
        i++;
        if((i > nhits-1 ) || (TMath::Abs((v[i]-v[start_index])) > width)) break;
    }
    return TMath::Abs(i - start_index)-dark;
}

Int_t SK2p2MeV::GetNXX(Int_t nhits, Float_t *t, Float_t twin, Float_t tcenter)
    // Cal hits in twin which centered at tcenter
    // nhits: number of hits in t[]
    // t[]  : ordered timings
    // twin : time window width
    // tcenter: center of time window
{
    Int_t nxx = 0;
    for (Int_t i=0; i<nhits; i++){
        if (i > MAXHITS) break;
        if ( t[i] < tcenter - twin/2. ) continue;
        if ( t[i] > tcenter + twin/2. ) break;
        nxx ++;
    }
    return nxx;
}

Int_t SK2p2MeV::GetNhits(Float_t *v, Int_t start_index, Float_t width, Int_t nhits)
    // Get nhits from v[0] in time window = width
    // v[] - hit timing, sorted in time
    // width - time window width
    // start_index - index of start pmt hit
    // nhits - total number  of v[] hits
{
    Int_t i = start_index;
    while (1) {
        i++;
        if((i > nhits-1 ) || (TMath::Abs((v[i]-v[start_index])) > width)) break;
    }
    return TMath::Abs(i - start_index);
}

Float_t SK2p2MeV::EffCos ( Float_t costh )
{
    // Cal effective cos theta
    // See also: /usr/local/sklib_g77/skofl_12a/src/sklib/coseffsk.F
    
    //if ( costh < -0.001 || costh > 1.001 ) {
      if ( costh < -1.001 || costh > 1.001 ) {
        std::cerr << " Invalid cos theta: " << costh << std::endl;
        throw std::runtime_error("bad angle");
        //exit (0);
    }
    
    return 0.205349 + 0.523981*costh + 0.389951*costh*costh
         - 0.131959 *costh*costh*costh;
}

Float_t SK2p2MeV::GetWeight(const Float_t pmt_coords[3], const Float_t vertex_coords[3]){
  const Float_t ATT_LEN = 9000.;

  const Float_t pmt_x = pmt_coords[0], pmt_y = pmt_coords[1], pmt_z = pmt_coords[2];
  const Float_t vertex_x = vertex_coords[0], vertex_y = vertex_coords[1], vertex_z = vertex_coords[2];

  const Float_t radius = sqrt(pow(pmt_x - vertex_x,2) +
                              pow(pmt_y - vertex_y,2) +
                              pow(pmt_z - vertex_z,2));
  
    const auto is_on_barrel = [](Float_t z){
    const double barrel_pmt_max_z = 1800;
    return (fabs(z) < barrel_pmt_max_z);
  };
  
  Float_t cos_theta = 0;
  
  if (is_on_barrel(pmt_z)){
    cos_theta = ( pmt_x * (pmt_x - vertex_x) +
                  pmt_y * (pmt_y - vertex_y) )
      / (sqrt(pow(pmt_x,2) + pow(pmt_y, 2)) * radius);
  }
  else {
    cos_theta  = (fabs(pmt_z - vertex_z) / radius);
  }
  
  return EffCos(cos_theta) * exp(-radius / ATT_LEN) / pow(radius,2);
}

Float_t SK2p2MeV::GetWeightThreshold (const Float_t *w, const Float_t frac)
{
    Float_t w1[MAXPM], w2[MAXPM];
    Int_t   index[MAXPM];
    Float_t tot = 0.;
    for (Int_t i=0; i<MAXPM; i++) {
        w2[i] = w[i];
        tot += w[i];
    }
    TMath::Sort(MAXPM, w2, index, kTRUE); // In decreasing order
    for (Int_t i=0; i<MAXPM; i++) {
        w1[i] = w2[ index[i] ];
    }
    Float_t t = 0;
    for (Int_t i=0; i<MAXPM; i++) {
        t += w1[i];
        if ( t/tot > frac ) {
            if ( verbosity > 2 ) {
                std::cout << " X/Y/Z= " << VX << " " << VY << " " << VZ << std::endl;
                std::cout << " # of high weight PMTs (frac=" << frac << "): " << i+1 << std::endl;
            }
            return w1[i];
        }
    }
    
    return -1.;
}

Bool_t SK2p2MeV::CheckCluster (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
        Int_t N10, Int_t &ncut, Float_t angle)
{
    // Search a cluster and cut hits in the cluster.
    // ux[], uy[], uz[] are unit vectors
    // flag[] indicates if this hit has already been cut.
    // Return true if cluster found.
    
    // Find the direction of the cluster
    const Float_t PI    = TMath::Pi();
    const Float_t COS10 = TMath::Cos(angle*PI/180.);
    Int_t theta_c, phi_c; // direction of the cluster
    Int_t nc = 0;
    for (Int_t theta=0; theta<180; theta+=2) {
        for (Int_t phi=0; phi<360; phi+=4) {
            Float_t x, y, z, cross;
            x = TMath::Sin(theta*PI/180.) * TMath::Cos(phi*PI/180.);
            y = TMath::Sin(theta*PI/180.) * TMath::Sin(phi*PI/180.);
            z = TMath::Cos(theta*PI/180.);
            Int_t n = 0;
            for (Int_t i=0; i<N10; i++) {
                if ( flag[i] != 0 ) continue;
                cross = x * ux[i] + y * uy[i] + z * uz[i];
                if ( cross > COS10 ) n ++;
            }
            if ( n > nc ) {
                nc  = n;
                theta_c = theta;
                phi_c   = phi;
            }
        }
    }
    
    if ( nc < 3 ) return kFALSE; // No cluster at all
    
    // Cut hits inside the cluster
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        Float_t cx = TMath::Sin(theta_c*PI/180.) * TMath::Cos(phi_c*PI/180.);
        Float_t cy = TMath::Sin(theta_c*PI/180.) * TMath::Sin(phi_c*PI/180.);
        Float_t cz = TMath::Cos(theta_c*PI/180.);
        if ( cx*ux[i] + cy*uy[i] + cz*uz[i] > COS10 ) {
            ncut ++;
            flag[i] = ncut;
        }
    }
    
    return kTRUE;
    
    /*
    //
    // solution 2
    //
    
    // Sort hits by theta
    Float_t uxt[N10], uxt2[N10], uyt[N10], uyt2[N10], uzt[N10], uzt2[N10];
    Int_t   oid[N10], oid2[N10]; // original index
    
    Int_t neff = 0;
    for (Int_t i=0; i<N10; i++) {
    if ( flag[i] != 0 ) continue;
    uxt2[neff] = ux[i];
    uyt2[neff] = uy[i];
    uzt2[neff] = uz[i];
    oid2[neff] = i;
    
    neff ++;
    }
    
    if ( neff < 3 ) return kFALSE; // No cluster at all
    
    Int_t index[neff];
    TMath::Sort(neff, uzt2, index, kFALSE); // Sort according to cos(theta)
    for (Int_t i=0; i<neff; i++){
    uxt[i] = uxt2[ index[i] ];
    uyt[i] = uyt2[ index[i] ];
    uzt[i] = uzt2[ index[i] ];
    oid[i] = oid2[ index[i] ];
    }
    
    // Search cluster, start from the lowest one (maximum theta)
    const Float_t PI    = TMath::Pi();
    const Float_t COS10 = TMath::Cos(10.*PI/180.);
    const Float_t COS20 = TMath::Cos(20.*PI/180.);
    Int_t ncm = 0;
    Int_t idm[N10];
    for (Int_t i=0; i<neff - 1; i++) {
    if ( uxt[i]*uxt[i+1] + uyt[i]*uyt[i+1] + uzt[i]*uzt[i+1] < COS20 ) continue;
    Float_t ucx, ucy, ucz;
    ucx = (uxt[i] + uxt[i+1])/2.;
    ucy = (uyt[i] + uyt[i+1])/2.;
    ucz = (uzt[i] + uzt[i+1])/2.;
    ucx = ucx / sqrt(ucx*ucx+ucy*ucy+ucz*ucz); // unit, unnecessary, only dir is important
    ucy = ucy / sqrt(ucx*ucx+ucy*ucy+ucz*ucz);
    ucz = ucz / sqrt(ucx*ucx+ucy*ucy+ucz*ucz);
    
    Int_t id[N10];
    Int_t nc = 0; // # of other hits fits in the small circle(r=COS10)
    for (Int_t j=0; j<neff; j++) {
    if ( j==i || j==i+1 ) continue;
    if ( ucx*uxt[j]+ucy*uyt[j]+ucz*uzt[j] > COS10 ) {
    id[nc] = j;
    nc ++;
    }
    }
    if ( nc+2 > ncm ) {
    for (Int_t k=0; k<nc; k++) idm[k] = id[k];
    idm[nc] = i;
    idm[nc+1] = i+1;
    ncm = nc + 2;
    }
    }
    
    if ( ncm < 3 ) return kFALSE; // No cluster at all
    
    // Cut hits in the cluster
    for (Int_t i=0; i<ncm; i++) {
    ncut ++;
    flag[ oid[ idm[i] ] ] = ncut;
    }
    
    return kTRUE;
    */
        
        /*
        //
        // solution 3
        //
        
        const Float_t DELTA_THETA = 20.; // degree
        const Float_t DELTA_PHI   = 20.;
        const Int_t N10MAX = N10 * 2;
        Float_t theta[N10MAX], theta2[N10MAX], theta3[N10MAX];
        Float_t phi[N10MAX],   phi2[N10MAX],   phi3[N10MAX];
        Int_t   oldi[N10MAX],  oldi2[N10MAX],  oldi3[N10MAX], index[N10MAX];
        
        Int_t neff = 0;
        for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        theta2[neff] = TMath::ACos(uz[i]);
        //if ( TMath::Abs(TMath::Abs(uz[i])-1.) < 0.000001 ) phi2[neff] = 0.;
        //else phi2[neff] = TMath::ATan2 (uy[i], ux[i]);
        phi2[neff] = TMath::ATan2 (uy[i], ux[i]);
        if ( phi2[neff] < 0. ) phi2[neff] += 2.*TMath::Pi();
        
        theta2[neff] = theta2[neff] * 180./TMath::Pi();
        phi2[neff]   = phi2[neff] * 180./TMath::Pi(); // in degree
        oldi2[neff]  = i;
        
        neff ++;
        }
        
        // NOTE: Different from GetNhits in tisk, the hits with phi near 2 pi
        //       should be calculated together with those near 0.
        Int_t oldn = neff;
        for ( Int_t i=0; i<oldn; i++ ) {
        if ( phi2[i] < DELTA_PHI ) {
        theta2[neff] = theta2[i];
        phi2[neff]   = phi2[i] + 360.; // Add 360 deg
        oldi2[neff]  = oldi2[i];
        
        neff ++;
        }
        }
        
        // Sort according to phi
        TMath::Sort(neff, phi2, index, kFALSE); // In increasing order
        for (Int_t i=0; i<neff; i++){
        theta[i] = theta2[ index[i] ];
        phi[i]   = phi2[ index[i] ];
        oldi[i]  = oldi2[ index[i] ];
        }
        
        // Find cluster in phi direction
        Int_t ncluster, n, id;
        ncluster = 0;
        for (Int_t i=0; i<neff; i++) {
        n = GetNhits(phi, i, DELTA_PHI, neff); // 20 deg
        if ( n > ncluster ) {
        ncluster = n;
        id = i;
        }
        }
        
        if ( ncluster < 2 ) return kFALSE; // No cluster at all
        
        for (Int_t i=0; i<ncluster; i++) {
        theta3[i] = theta[i+id];
        phi3[i] = phi[i+id];
        oldi3[i] = oldi[i+id];
        }
        
        neff = ncluster;
        // Sort according to theta
        TMath::Sort(neff, theta3, index, kFALSE); // In increasing order
        for (Int_t i=0; i<neff; i++){
            theta[i] = theta3[ index[i] ];
            phi[i]   = phi3[ index[i] ];
            oldi[i]  = oldi3[ index[i] ];
        }
    
    // Find cluster in theta direction
    ncluster = 0;
    for (Int_t i=0; i<neff; i++) {
        n = GetNhits(theta, i, DELTA_THETA, neff); // 20 deg
        if ( n > ncluster ) {
            ncluster = n;
            id = i;
        }
    }
    
    //
    if ( ncluster < 3 ) return kFALSE; // No cluster at all
    
    for ( Int_t i=0; i<ncluster; i++ ) {
        if ( flag[ oldi[i+id] ] == 0 ) {
            ncut ++;
            flag[ oldi[i+id] ] = ncut;
        }
    }
    
    return kTRUE;
    */
}

Int_t SK2p2MeV::GetCluster (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
        Int_t N10, Int_t &ncut, Int_t ncth, Float_t thr)
{
    // Search cluster.
    
    // Ordering is not necessary, keep the orginal order
    // Sort hits according to uz
    //     Float_t ux1[N10], uy1[N10], uz1[N10];
    //     Int_t   index[N10];
    //     Int_t   neff = 0;
    //     for (Int_t i=0; i<N10; i++) {
    //         if ( flag[i] != 0 ) continue; // skip hits that are already in clusters
    //         ux1[i] = ux[i];
    //         uy1[i] = uy[i];
    //         uz1[i] = uz[i];
    //         neff ++;
    //      }
    //     TMath::Sort(N10, uz1, index, kFALSE); // In increasing order
    //     for (Int_t i=0; i<N10; i++) {
    //         ux[i] = ux1[ index[i] ];
    //         uy[i] = uy1[ index[i] ];
    //         uz[i] = uz1[ index[i] ];
    //     }
    Int_t index[N10];
    Int_t nc, nc_m;
    for (Int_t i=0; i<N10; i++) {
        if (  flag[i] != 0 ) continue; // skip hits that are already in clusters
        for (Int_t j=0; j<N10; j++) index[j] = 0;
        index [i] = 1; // first hit in a cluster
        nc_m = 0;
        while (1) {
            for (Int_t j=0; j<N10; j++) {
                if (  flag[j] != 0 ) continue; // skip hits that are already in clusters
                // scan index
                for (Int_t k=0; k<N10; k++) {
                    if ( k == j ) continue;       // do not compute angle to itself
                    if ( index[k] == 0 ) continue;// skip non-candidate hit
                    Float_t dalpha = ux[j]*ux[k]+uy[j]*uy[k]+uz[j]*uz[k];
                    if ( dalpha > thr ) index[j] = 1;
                }
            }
            // count hits in cluster
            nc = 0;
            for (Int_t j=0; j<N10; j++) {
                if ( index[j] == 1 ) nc ++;
            }
            if ( nc > nc_m ) nc_m = nc;
            else break; // no more hits belong to this cluster, stop
        }
        if ( nc_m >= ncth ) { // this is a cluster
            for (Int_t j=0; j<N10; j++) {
                if ( index[j] == 1 ) {
                    ncut ++; //the initial value of ncut = 0
                    flag[j] = ncut;
                }
            }
        }
    }
    return 1;
}

Bool_t SK2p2MeV::CheckHighQ (Float_t *qi, Int_t *flag, Int_t N10, Int_t &ncut, Float_t qth)
{
    // Cut hits with Q>QTH
    // Return true if bad Q found.
    
    Int_t nbad = 0;
    for ( Int_t i=0; i<N10; i++ ) {
        if ( flag[i] != 0 ) continue;
        if ( qi[i] > qth) {
            ncut ++;
            flag[i] = ncut;
            nbad ++;
        }
    }
    
    return ( nbad > 0 ) ? kTRUE : kFALSE;
}

Bool_t SK2p2MeV::CheckLowQ (Float_t *qi, Int_t *flag, Int_t N10, Int_t &ncut, Float_t qth)
{
    // Cut hits with Q<QTH
    // Return true if bad Q found.
    
    Int_t nbad = 0;
    for ( Int_t i=0; i<N10; i++ ) {
        if ( flag[i] != 0 ) continue;
        if ( qi[i] < qth) {
            ncut ++;
            flag[i] = ncut;
            nbad ++;
        }
    }
    
    return ( nbad > 0 ) ? kTRUE : kFALSE;
}

Float_t SK2p2MeV::GetRatio (Int_t *ci, Float_t *wt, Int_t *flag, Int_t N10)
{
    Float_t tot_wt = 0.;
    for (Int_t i=0; i<MAXPM; i++) {
        tot_wt += wt[i];
    }
    Int_t   neff = 0;
    Float_t nwt = 0.;
    Float_t ratio = 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        nwt += wt[ ci[j] - 1];
        neff ++;
    }
    if ( neff == 0 ) ratio = 0.;
    else ratio = nwt * MAXPM / tot_wt / neff;
    
    return ratio;
}

//*********new added by Yang Zhang*****
Float_t SK2p2MeV::GetQ10 (Float_t *qi, Int_t *flag, Int_t N10)
{
    Float_t q10 = 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        q10 += qi[j];
    }
    return q10 ;
}

Float_t SK2p2MeV::GetQmean (Float_t *qi, Int_t *flag, Int_t N10)
{
    Int_t   neff = 0;
    Float_t qmean= 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        qmean += qi[j];
        neff ++;
    }
    if ( neff < 2 ) return -1.;
    qmean = qmean/neff;
    return qmean;
}

Float_t SK2p2MeV::GetQrms (Float_t *qi, Int_t *flag, Int_t N10)
{
    Int_t   neff = 0;
    Float_t qmean= 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        qmean += qi[j];
        neff ++;
    }
    if ( neff < 2 ) return -1.;
    qmean = qmean/neff;
    Float_t qrms = 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        qrms += (qi[j] - qmean) * (qi[j] - qmean);
    }
    qrms = TMath::Sqrt(qrms/neff);
    return qrms;
}

Float_t SK2p2MeV::GetThetaRms (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10)
{
    Int_t neff = 0;
    Float_t dir[3];
    dir[0] = dir[1] = dir[2] = 0.;
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        dir[0] += ux[i];
        dir[1] += uy[i];
        dir[2] += uz[i];
        neff ++;
    }
    if ( neff < 2 ) return 0.;
    Float_t vr = sqrt(dir[0]*dir[0] + dir[1]*dir[1]
            + dir[2]*dir[2]);
    dir[0] = dir[0] / vr;
    dir[1] = dir[1] / vr;
    dir[2] = dir[2] / vr;
    
    Float_t theta;
    Float_t mean = 0.;
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        theta = dir[0]*ux[i] + dir[1]*uy[i] + dir[2]*uz[i];
        theta = TMath::ACos(theta)*180./TMath::Pi();
        mean += theta;
    }
    mean = mean / neff;
    Float_t rms = 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        theta = dir[0]*ux[j] + dir[1]*uy[j] + dir[2]*uz[j];
        theta = TMath::ACos(theta)*180./TMath::Pi();
        rms += (theta - mean) * (theta - mean);
    }
    rms = TMath::Sqrt(rms/neff);
    return rms;
}

Int_t SK2p2MeV::GetNLowTheta (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10, Float_t thetaTh)
{
    Int_t nlow = 0;
    Float_t dir[3];
    dir[0] = dir[1] = dir[2] = 0.;
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        dir[0] += ux[i];
        dir[1] += uy[i];
        dir[2] += uz[i];
    }
    Float_t vr = sqrt(dir[0]*dir[0] + dir[1]*dir[1]
            + dir[2]*dir[2]);
    dir[0] = dir[0] / vr;
    dir[1] = dir[1] / vr;
    dir[2] = dir[2] / vr;
    
    Float_t theta;
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        theta = dir[0]*ux[i] + dir[1]*uy[i] + dir[2]*uz[i];
        theta = TMath::ACos(theta)*180./TMath::Pi();
        if(theta<thetaTh) nlow ++;
    }
    
    return nlow;
}
//**********************************

Float_t SK2p2MeV::GetTrms (Float_t *ti, Int_t *flag, Int_t N10)
{
    Int_t   neff = 0;
    Float_t tmean= 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        tmean += ti[j];
        neff ++;
    }
    if ( neff < 2 ) return -1.;
    tmean = tmean/neff;
    Float_t trms = 0.;
    for (Int_t j=0; j<N10; j++) {
        if ( flag[j] != 0 ) continue;
        trms += (ti[j] - tmean) * (ti[j] - tmean);
    }
    trms = TMath::Sqrt(trms/neff);
    return trms;
}

Bool_t SK2p2MeV::CheckBackHits (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
        Int_t N10, Int_t &ncut, Float_t angle)
{
    const Float_t THR = TMath::Cos(angle*3.14159265/180.);
    
    // estimate dir
    Float_t dir[3];
    dir[0] = dir[1] = dir[2] = 0.;
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        dir[0] += ux[i];
        dir[1] += uy[i];
        dir[2] += uz[i];
    }
    Float_t r = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0] = dir[0] / r;
    dir[1] = dir[1] / r;
    dir[2] = dir[2] / r;
    // cut hits going backward
    for ( Int_t i=0; i< N10; i++ ) {
        if ( flag[i] != 0 ) continue;
        Float_t cross = dir[0] * ux[i] + dir[1] * uy[i] + dir[2] * uz[i];
        if ( cross < THR ) {
            ncut ++;
            flag[i] = ncut;
        }
    }
    return kTRUE;
}

Int_t SK2p2MeV::GetLowHits (Int_t *ci, Int_t N10, Float_t *wt, Float_t acceptance)
{
    // Cut hits with low hit probability
    
    // get weight threshold
    Float_t wlow  = GetWeightThreshold(wt, acceptance);
    
    Int_t nlow = 0;
    for (Int_t j=0; j<N10; j++) {
        if ( wt[ ci[j] - 1] < wlow  ) nlow ++;
    }
    
    return nlow;
}

Bool_t SK2p2MeV::CheckCherenkovLike (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
        Int_t N10, Int_t &ncut)
{
    // Search maximum # of hits in a Cherenkov cone (angle=70 deg) and
    // cut hits outside this cone.
    // ux[], uy[], uz[] are unit vectors
    // flag[] indicates if this hit has already been cut.
    // Return true if no hits found outside Cherenkov cone
    
    // Find the direction of Cherenkov cone
    const Float_t PI    = TMath::Pi();
    const Float_t COS70 = TMath::Cos(70.*PI/180.);
    Int_t theta_c, phi_c; // Direction of Cherenkov cone
    Int_t nclike = 0;
    for (Int_t theta=0; theta<180; theta+=2) {
        for (Int_t phi=0; phi<360; phi+=4) {
            Float_t x, y, z, cross;
            x = TMath::Sin(theta*PI/180.) * TMath::Cos(phi*PI/180.);
            y = TMath::Sin(theta*PI/180.) * TMath::Sin(phi*PI/180.);
            z = TMath::Cos(theta*PI/180.);
            Int_t n = 0;
            for (Int_t i=0; i<N10; i++) {
                if ( flag[i] != 0 ) continue;
                cross = x * ux[i] + y * uy[i] + z * uz[i];
                if ( cross > COS70 ) n ++;
            }
            if ( n > nclike ) {
                nclike  = n;
                theta_c = theta;
                phi_c   = phi;
            }
        }
    }
    
    if ( nclike == N10 - ncut ) return kTRUE;
    
    // Cut hits outside the Chrenkov cone
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        Float_t cx = TMath::Sin(theta_c*PI/180.) * TMath::Cos(phi_c*PI/180.);
        Float_t cy = TMath::Sin(theta_c*PI/180.) * TMath::Sin(phi_c*PI/180.);
        Float_t cz = TMath::Cos(theta_c*PI/180.);
        if ( cx*ux[i] + cy*uy[i] + cz*uz[i] < COS70 ) {
            ncut ++;
            flag[i] = ncut;
        }
    }
    
    return kFALSE;
}

Float_t SK2p2MeV::GetDirKS (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10)
{
    // Calculate standard deviation of Phi of hit vectors
    Float_t vect_sum[3], angle_x[N10], angle_y[N10], angle_z[N10];
    Int_t neff = 0;
    //std::cout << "C event " << N10 << std::endl;
    for (Int_t j=0; j<N10; j++){
        if ( flag[j] != 0 ) continue;
        angle_x[neff] = ux[j];
        angle_y[neff] = uy[j];
        angle_z[neff] = uz[j];
        //std::cout << ux[j] << " " << uy[j] << " " << uz[j] << std::endl;
        neff ++;
    }
    if ( neff < 2 ) return -1.;
    // calculate sum of hit vectors (estimate dir)
    vect_sum[0] = vect_sum[1] = vect_sum[2] = 0.;
    for (Int_t i=0; i<neff; i++) {
        vect_sum[0] += ux[i];
        vect_sum[1] += uy[i];
        vect_sum[2] += uz[i];
    }
    Float_t vr = sqrt(vect_sum[0]*vect_sum[0] + vect_sum[1]*vect_sum[1]
            + vect_sum[2]*vect_sum[2]);
    vect_sum[0] = vect_sum[0] / vr;
    vect_sum[1] = vect_sum[1] / vr;
    vect_sum[2] = vect_sum[2] / vr;
    // calculate phi_div
    TVector3 hitv[neff], pv[neff], dirv, hitvnew[neff];
    Float_t  phiv[neff];
    dirv.SetXYZ(vect_sum[0], vect_sum[1], vect_sum[2]); //sum vector
    for (Int_t j=0; j<neff; j++) {
        Float_t cross;
        Float_t scale;
        cross = vect_sum[0] * angle_x[j] + vect_sum[1] * angle_y[j]
            + vect_sum[2] * angle_z[j];
        scale = 1 / cross;
        hitv[j].SetXYZ (angle_x[j] * scale, angle_y[j] * scale, angle_z[j] * scale);
        hitvnew[j].SetXYZ (angle_x[j], angle_y[j], angle_z[j]);
        //perpendicular vector
        pv[j] = hitv[j] - dirv;
        // Rotate vdir to z axis
        pv[j].Rotate(-dirv.Phi(), TVector3(0., 0., 1.));
        pv[j].Rotate(-dirv.Theta(), TVector3(0., 1., 0.));
        phiv[j] = pv[j].Phi();
        if ( phiv[j] < 0. ) phiv[j] += 2.*TMath::Pi();
        hitvnew[j].Rotate(-dirv.Phi(), TVector3(0., 0., 1.));
        hitvnew[j].Rotate(-dirv.Theta(), TVector3(0., 1., 0.));
        
        phiv[j]=hitvnew[j].Phi();
        if(phiv[j]<0) phiv[j]+=2.*TMath::Pi();
    }
    Float_t phiv2[neff];
    Int_t index[neff];
    for (Int_t j=0; j<neff; j++) phiv2[j] = phiv[j];
    TMath::Sort(neff, phiv2, index, kFALSE); // Sort phi in increasing order
    for (Int_t j=0; j<neff; j++){
        phiv[j] = phiv2[ index[j] ];
    }
    // SONIA: hacks function to compute dirks
    Float_t avephi = 2.*TMath::Pi() / neff;
    float pksmax = -2 * TMath::Pi();
    float pksmin = 2 * TMath::Pi();
    for (Int_t j=1; j<neff; j++) {
        float pks = phiv[j] - avephi * j;
        if (pks > pksmax) pksmax = pks;
        if (pks < pksmin) pksmin = pks;
    }
    Float_t dks = (pksmax - pksmin)/(2 * TMath::Pi());
    return dks;
}

Float_t SK2p2MeV::GetPhiRms (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10)
{
    // Calculate standard deviation of Phi of hit vectors
    Float_t vect_sum[3], angle_x[N10], angle_y[N10], angle_z[N10];
    Int_t neff = 0;
    for (Int_t j=0; j<N10; j++){
        if ( flag[j] != 0 ) continue;
        angle_x[neff] = ux[j];
        angle_y[neff] = uy[j];
        angle_z[neff] = uz[j];
        neff ++;
    }
    if ( neff < 2 ) return -1.;
    // calculate sum of hit vectors (estimate dir)
    vect_sum[0] = vect_sum[1] = vect_sum[2] = 0.;
    for (Int_t i=0; i<neff; i++) {
        vect_sum[0] += ux[i];
        vect_sum[1] += uy[i];
        vect_sum[2] += uz[i];
    }
    Float_t vr = sqrt(vect_sum[0]*vect_sum[0] + vect_sum[1]*vect_sum[1]
            + vect_sum[2]*vect_sum[2]);
    vect_sum[0] = vect_sum[0] / vr;
    vect_sum[1] = vect_sum[1] / vr;
    vect_sum[2] = vect_sum[2] / vr;
    // calculate phi_div
    TVector3 hitv[neff], pv[neff], dirv, hitvnew[neff];
    Float_t  phiv[neff];
    dirv.SetXYZ(vect_sum[0], vect_sum[1], vect_sum[2]); //sum vector
    for (Int_t j=0; j<neff; j++) {
        Float_t cross;
        Float_t scale;
        cross = vect_sum[0] * angle_x[j] + vect_sum[1] * angle_y[j]
            + vect_sum[2] * angle_z[j];
        scale = 1 / cross;
        hitv[j].SetXYZ (angle_x[j] * scale, angle_y[j] * scale, angle_z[j] * scale);
        hitvnew[j].SetXYZ (angle_x[j], angle_y[j], angle_z[j]);
        //perpendicular vector
        pv[j] = hitv[j] - dirv;
        // Rotate vdir to z axis
        pv[j].Rotate(-dirv.Phi(), TVector3(0., 0., 1.));
        pv[j].Rotate(-dirv.Theta(), TVector3(0., 1., 0.));
        phiv[j] = pv[j].Phi();
        if ( phiv[j] < 0. ) phiv[j] += 2.*TMath::Pi();
        hitvnew[j].Rotate(-dirv.Phi(), TVector3(0., 0., 1.));
        hitvnew[j].Rotate(-dirv.Theta(), TVector3(0., 1., 0.));
        
        phiv[j]=hitvnew[j].Phi();
        if(phiv[j]<0) phiv[j]+=2.*TMath::Pi();
    }
    Float_t phiv2[neff];
    Int_t index[neff];
    for (Int_t j=0; j<neff; j++) phiv2[j] = phiv[j];
    TMath::Sort(neff, phiv2, index, kFALSE); // Sort phi in increasing order
    for (Int_t j=0; j<neff; j++){
        phiv[j] = phiv2[ index[j] ];
    }
    Float_t avephi = 2.*TMath::Pi() / neff;
    Float_t delta_phi;
    Float_t phi_div = 0.;
    for (Int_t j=1; j<neff; j++) {
        delta_phi = phiv[j] - phiv[j-1];
        phi_div += (delta_phi - avephi) * (delta_phi - avephi);
    }
    delta_phi = phiv[0] + 2.*TMath::Pi() - phiv[neff-1];
    phi_div += (delta_phi - avephi) * (delta_phi - avephi);
    phi_div = TMath::Sqrt(phi_div/neff);
    return phi_div*180/TMath::Pi();
}

Float_t SK2p2MeV::GetThetaMean (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N10)
{
    Int_t neff = 0;
    Float_t dir[3];
    dir[0] = dir[1] = dir[2] = 0.;
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        dir[0] += ux[i];
        dir[1] += uy[i];
        dir[2] += uz[i];
        neff ++;
    }
    if ( neff < 2 ) return 0.;
    Float_t vr = sqrt(dir[0]*dir[0] + dir[1]*dir[1]
            + dir[2]*dir[2]);
    dir[0] = dir[0] / vr;
    dir[1] = dir[1] / vr;
    dir[2] = dir[2] / vr;
    
    Float_t theta;
    Float_t mean = 0.;
    for (Int_t i=0; i<N10; i++) {
        if ( flag[i] != 0 ) continue;
        theta = dir[0]*ux[i] + dir[1]*uy[i] + dir[2]*uz[i];
        theta = TMath::ACos(theta)*180./TMath::Pi();
        mean += theta;
    }
    mean = mean / neff;
    
    return mean;
}

Float_t SK2p2MeV::MinimizeTrms(Float_t* tiskz, Int_t* cabiz, Int_t startindex, Int_t* index, Int_t nhits, Float_t& CVX, Float_t& CVY, Float_t& CVZ, Float_t pVX, Float_t pVY, Float_t pVZ, Float_t discut)
{
    Float_t inc;
    (discut > 200) ? inc = 100 : inc = discut / 2;
    Float_t tiskzmin[nhits], tiskz2[nhits];
    Int_t cabizmin[nhits];
    for (int i = 0; i < nhits; i++)
    {
        //std::cout << "hit at " << tiskz[i] << std::endl;
        tiskz2[i] = tiskz[i];
    }
    Int_t maxscanxy, maxscanz;
    maxscanz = (Int_t)(2*ZPINTK/(float)inc);
    maxscanxy = (Int_t)(2*RINTK/(float)inc);
    Float_t VX = 0;
    Float_t VY = 0;
    Float_t VZ = 0;
    Float_t tVX, tVY, tVZ, vx, vy, vz;
    //added for when cutdis is small
    // inc = inc / 4.;
    
    Float_t mintrms = 9999;
    Float_t mintmean = 0;
    Float_t tmean, trms, trms1;
    // VX = pVX;
    // VY = pVY;
    // VZ = pVZ;
    while (inc >0.5)
    {
        //  std::cout << "inc = " << inc << std::endl;
        for (Float_t x = 0; x < maxscanxy; x++)
        {
            //std::cout << x << " ";
            vx = inc*(x-maxscanxy/2.) + VX;
            for (Float_t y = 0; y < maxscanxy; y++)
            {
                //std::cout << y << " ";
                vy = inc*(y-maxscanxy/2.) + VY;
                if (TMath::Sqrt(vx*vx + vy*vy) > RINTK) continue;
                for (Float_t z = 0; z < maxscanz; z++)
                {
                    //  std::cout << z << std::endl;
                    vz = inc*(z-maxscanz/2.) + VZ;
                    if (vz > ZPINTK || vz < -ZPINTK) continue;
                    Float_t dis = TMath::Sqrt((vx - pVX)*(vx - pVX) + (vy - pVY)*(vy - pVY) + (vz - pVZ) * (vz - pVZ));
                    if (dis > discut) continue;
                    //          std::cout << "passed discut" << std::endl;
                    tmean = 0;
                    trms = 0;
                    for (int i=0; i<nhits; i++)
                    {
                        tiskz[i] = tiskz2[i];
                        Float_t tof;
                        //std::cout<<vx<<" "<<cabiz[i]-1<<" "<<xyz[cabiz[i]-1][0]<<std::endl;
                        tof = TMath::Sqrt((vx - xyz[cabiz[i]-1][0]) * (vx - xyz[cabiz[i]-1][0])
                                +(vy - xyz[cabiz[i]-1][1]) * (vy - xyz[cabiz[i]-1][1])
                                +(vz - xyz[cabiz[i]-1][2]) * (vz - xyz[cabiz[i]-1][2])) / C_WATER;
                        tiskz[i] -= tof;
                        //std::cout << "tiskz in min: " << tof << std::endl;
                        tmean += tiskz[i];
                    }
                    tmean = tmean/nhits;
                    for (int i=0; i<nhits; i++)
                    {
                        trms += (tiskz[i] - tmean)*(tiskz[i] - tmean);
                    }
                    trms = TMath::Sqrt(trms/nhits);
                    //std::cout << "trms = " << trms << std::endl;
                    if (trms < mintrms)
                    {
                        mintrms = trms;
                        mintmean = tmean;
                        tVX = vx;
                        tVY = vy;
                        tVZ = vz;
                        //        std::cout << "new best: " << mintrms << std::endl;
                        for (int i=0; i<nhits; i++)
                        {
                            //std::cout << tiskz[i] << std::endl;
                            
                            tiskzmin[i] = tiskz[i];
                            cabizmin[i] = cabiz[i];
                        
                        }
                    
                    }
                }
            }
        }
        VX = tVX;
        VY = tVY;
        VZ = tVZ;
        inc = inc / 2.;
    }
    //std::cout<<VX<<" "<<pVX<<" "<<mintrms<<std::endl;
    TMath::Sort(nhits, tiskzmin, index, kFALSE); // In increasing order
    for (int i = 0; i < nhits; i++)
    {
        tiskz[i] = tiskzmin[index[i]];
        cabiz[i] = cabizmin[index[i]];
        //std::cout<<tiskz[i]<<" "<<cabiz[i]<<std::endl;
    }
    /*delete [] tiskz2;
      delete [] tiskzmin;
      delete [] cabizmin;
      tiskz2 = NULL;
      tiskzmin = NULL;
      cabizmin = NULL;
      */
    CVX = VX;
    CVY = VY;
    CVZ = VZ;
    return mintrms;
}

Float_t SK2p2MeV::BasicTof(Float_t* tiskz, Int_t* cabiz, Float_t VX, Float_t VY, Float_t VZ, Int_t nhits, Int_t* index)
{
    Float_t tmean = 0;
    Float_t trms = 0;
    for (int i=0; i<nhits; i++)
    {
        if (i > MAXHITS)  break;
        Float_t tof;
        //std::cout << cabiz[i] << std::endl;
        tof = TMath::Sqrt((VX - xyz[cabiz[i]-1][0]) * (VX - xyz[cabiz[i]-1][0])
                +(VY - xyz[cabiz[i]-1][1]) * (VY - xyz[cabiz[i]-1][1])
                +(VZ - xyz[cabiz[i]-1][2]) * (VZ - xyz[cabiz[i]-1][2])) / C_WATER;
        //std::cout << "hits: " << tiskz2[i] << " - tof: " << tiskz2[i]-tof << std::endl;
        //       std::cout << "tof " << tof << std::endl;
        tiskz[i] -= tof;
        tmean += tiskz[i];
    }
    tmean = tmean/nhits;
    Float_t tiskz1[nhits];
    Int_t cabiz1[nhits];
    TMath::Sort(nhits, tiskz, index, kFALSE);
    
    for (Int_t i = 0; i<nhits; i++)
    {
        if (i > MAXHITS) break;
        tiskz1[i] = tiskz[index[i]];
        cabiz1[i] = cabiz[index[i]];
        trms += (tiskz1[i] - tmean)*(tiskz1[i]-tmean);
    }
    for (Int_t i = 0; i<nhits; i++)
    {
        tiskz[i] = tiskz1[i];
        cabiz[i] = cabiz1[i];
    }
    trms = TMath::Sqrt(trms/nhits);
    //  delete [] xyz;
    return trms;

}

void SK2p2MeV::GetMinTrms (Float_t *ti, Int_t N10, Float_t &mintrms_6, Float_t &mintrms_5, Float_t &mintrms_4, Float_t &mintrms_3, Int_t *flag)
{
    Int_t usedHits[6];
    Float_t trms_3, trms_4, trms_5, trms_6, tmean_3, tmean_4, tmean_5, tmean_6;
    mintrms_3 = -1.;
    mintrms_4 = -1.;
    mintrms_5 = -1.;
    mintrms_6 = -1.;
    
    Float_t tis[N10];
    Int_t t = 0;
    for (Int_t j = 0; j < N10; j++)
    {
        //if (flag[j] != 0) continue;
        tis[t] = ti[j];
        t++;
    }
    
    //for blocks of 6 hit
    if (t >= 6)
    {
        for (Int_t j = 0; j < (N10-5); j++)
        {
            usedHits[0] = j;
            for (Int_t k = j+1; k < (N10-4); k++)
            {
                usedHits[1] = k;
                for (Int_t l = k + 1; l < (N10-3); l++)
                {
                    usedHits[2] = l;
                    for (Int_t m = l + 1; m < (N10-2); m++)
                    {
                        usedHits[3] = m;
                        for (Int_t n = m + 1; n < (N10 - 1); n++)
                        {
                            usedHits[4] = n;
                            for (Int_t o = n + 1; o < N10; o++)
                            {
                                usedHits[5] = o;
                                //             std::cout << usedHits[0] << " " << usedHits[1] << " " << usedHits[2] << " " << usedHits[3] << std::endl;
                                trms_6 = 0.;
                                tmean_6 = 0.;
                                for (Int_t p = 0; p < 6; p++)
                                {
                                    //std::cout << tmean_6 << std::endl;
                                    tmean_6 += tis[usedHits[p]];
                                }
                                tmean_6 = tmean_6/6.;
                                
                                for (Int_t p = 0; p < 6; p++)
                                {
                                    trms_6 += (tis[usedHits[p]] - tmean_6)*(tis[usedHits[p]] - tmean_6);
                                }
                                
                                trms_6 = TMath::Sqrt(trms_6/6.);
                                
                                if ((trms_6 < mintrms_6) || (mintrms_6 < 0)) mintrms_6 = trms_6;
                            
                            }
                        }
                    
                    }
                }
            }
        }
    }
    
    //blocks of 5 hits
    
    if (t >= 5)
    {
        for (Int_t k = 0; k < (N10-4); k++)
        {
            usedHits[0] = k;
            for (Int_t l = k + 1; l < (N10-3); l++)
            {
                usedHits[1] = l;
                for (Int_t m = l + 1; m < (N10-2); m++)
                {
                    usedHits[2] = m;
                    for (Int_t n = m + 1; n < (N10 - 1); n++)
                    {
                        usedHits[3] = n;
                        for (Int_t o = n + 1; o < N10; o++)
                        {
                            usedHits[4] = o;
                            //             std::cout << usedHits[0] << " " << usedHits[1] << " " << usedHits[2] << " " << usedHits[3] << std::endl;
                            trms_5 = 0;
                            tmean_5 = 0;
                            for (Int_t p = 0; p < 5; p++)
                            {
                                tmean_5 += tis[usedHits[p]];
                            }
                            tmean_5 = tmean_5/5.;
                            for (Int_t p = 0; p < 5; p++)
                            {
                                trms_5 += (tis[usedHits[p]] - tmean_5)*(tis[usedHits[p]] - tmean_5);
                            }
                            trms_5 = TMath::Sqrt(trms_5/5.);
                            if ((trms_5 < mintrms_5) || (mintrms_5 < 0)) mintrms_5 = trms_5;
                        }
                    }
                }
            }
        }
    }
    // blocks of 4 hits
    if (t >= 4)
    {
        for (Int_t l = 0; l < (N10-3); l++)
        {
            usedHits[0] = l;
            for (Int_t m = l + 1; m < (N10-2); m++)
            {
                usedHits[1] = m;
                for (Int_t n = m + 1; n < (N10 - 1); n++)
                {
                    usedHits[2] = n;
                    for (Int_t o = n + 1; o < N10; o++)
                    {
                        usedHits[3] = o;
                        //             std::cout << usedHits[0] << " " << usedHits[1] << " " << usedHits[2] << " " << usedHits[3] << std::endl;
                        trms_4 = 0;
                        tmean_4 = 0;
                        for (Int_t p = 0; p < 4; p++)
                        {
                            tmean_4 += tis[usedHits[p]];
                        }
                        tmean_4 = tmean_4/4.;
                        for (Int_t p = 0; p < 4; p++)
                        {
                            trms_4 += (tis[usedHits[p]] - tmean_4)*(tis[usedHits[p]] - tmean_4);
                        }
                        trms_4 = TMath::Sqrt(trms_4/4.);
                        if ((trms_4 < mintrms_4) || (mintrms_4 < 0)) mintrms_4 = trms_4;
                    }
                }
            
            }
        }
    }
    
    //blocks of 3 hits
    if (t >= 3)
    {
        for (Int_t m = 0; m < (N10-2); m++)
        {
            usedHits[0] = m;
            for (Int_t n = m + 1; n < (N10 - 1); n++)
            {
                usedHits[1] = n;
                for (Int_t o = n + 1; o < N10; o++)
                {
                    usedHits[2] = o;
                    //             std::cout << usedHits[0] << " " << usedHits[1] << " " << usedHits[2] << " " << usedHits[3] << std::endl;
                    trms_3 = 0;
                    tmean_3 = 0;
                    for (Int_t p = 0; p < 3; p++)
                    {
                        tmean_3 += tis[usedHits[p]];
                    }
                    tmean_3 = tmean_3/3.;
                    for (Int_t p = 0; p < 3; p++)
                    {
                        trms_3 += (tis[usedHits[p]] - tmean_3)*(tis[usedHits[p]] - tmean_3);
                    }
                    trms_3 = TMath::Sqrt(trms_3/3.);
                    if ((trms_3 < mintrms_3) || (mintrms_3 < 0)) mintrms_3 = trms_3;
                }
            }
        
        }
    }
}

void SK2p2MeV::SetDarkRate(Int_t run)
{
    std::cout << " *** This is obsolete now! Please don't use it. *** " << std::endl;
    
    //     TString fname;
    //     fname.Form("/skofl/const/darkr/darkr.%06d.root", run);
    //     std::cout << " Reading dark rate from " << fname.Data() << std::endl;
    //     TChain *ch = new TChain("skdark");
    //     ch->Add(fname.Data());
    //     ch->SetBranchAddress ("dark_rate",  dark_rate);
    //     ch->GetEntry(0);
    //     for (Int_t i=0; i<11146; i++) {
    //         std::cout << i << " " << dark_rate[i] << std::endl;
    //     }
}
