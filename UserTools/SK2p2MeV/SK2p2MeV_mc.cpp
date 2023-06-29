#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"

#include "tqrealroot.h"
#include "loweroot.h"
#include "mcinfo.h"
#include "thirdredvars.h"

#include "SK2p2MeV_mc.h"

#include <iostream>
#include <fstream>

SK2p2MeV_mc::SK2p2MeV_mc (const Float_t (*geomxyz)[3]) 
    : SK2p2MeV (geomxyz)
{
    // additional input branches
    MC = new MCInfo;
    //    thirdred = new ThirdRed;
    
    // additional output branches
    head0 = new Header;
    lowe0 = new LoweInfo;
    mc0 = new MCInfo;
    //    third0 = new ThirdRed;
    
}

SK2p2MeV_mc::~SK2p2MeV_mc (){
    theOTree->ResetBranchAddresses();
    
    delete MC;
    delete thirdred;
    
    delete head0;
    delete lowe0;
    //    delete third0;
    delete mc0;
}

void SK2p2MeV_mc::Print()
{
    // Show current status
    
    std::cout << " Class: SK2p2MeV_mc" << std::endl;
    std::cout << "    AFT_GATE: " << AFT_GATE << std::endl;
    std::cout << "      fSmear: " << fSmear << std::endl;
    std::cout << "       sigma: x=" << sigma_xyz[0] << " y=" << sigma_xyz[1]
         << " z=" << sigma_xyz[2] << std::endl;
    std::cout << "      fShift: " << fShift << std::endl;
    std::cout << "      rshift: " << rshift << std::endl;
}

void SK2p2MeV_mc::SetSeed (int seed) 
{
    rseed = seed;
}

void SK2p2MeV_mc::SetSmearFlag (Bool_t flag) 
{
    fSmear = flag;
}

void SK2p2MeV_mc::SetMCFlag (Bool_t flag) 
{
    useMC = flag;
}

void SK2p2MeV_mc::SetVertexResolution (const Float_t x, const Float_t y, const Float_t z)
{
    sigma_xyz[0] = x;
    sigma_xyz[1] = y;
    sigma_xyz[2] = z;
}

void SK2p2MeV_mc::SetShiftFlag (Bool_t flag) 
{
    fShift = flag;
}

void SK2p2MeV_mc::SetShiftDistance (const Float_t dis)
{
    rshift = dis;
}

bool SK2p2MeV_mc::GetBranchValues(){
    // get base class branches
    bool get_ok = SK2p2MeV::GetBranchValues();
    if(not get_ok){
        return false;
    }
    
    // get additional branch variables - MU and ThirdRed propagated to output TTree
    get_ok  = (myTreeReader->Get("MC", MC));                   //&&
    //get_ok &= (myTreeReader->Get("ThirdRed", thirdred));
    
    return get_ok;
}

bool SK2p2MeV_mc::Initialise(MTreeReader* reader){
    
    myTreeReader = reader;
    
    theOTree->Branch("HEADER", "Header", &head0, 1024*1024, 0);
    theOTree->Branch("LOWE", "LoweInfo", &lowe0, 1024*1024, 0);
    theOTree->Branch("MC", "MCInfo", &MC, 1024*1024, 0);
    // theOTree->Branch("ThirdRed", "ThirdRed", &third0, 1024*1024, 0);
    theOTree->Branch("smearedvertex", smearedvertex, "smearedvertex[3]/F");
    
    prevrun = 0;
    prevsub = 0;
    fSmear = kFALSE;
    sigma_xyz[0] = sigma_xyz[1] = sigma_xyz[2] = 0.;
    fShift = kFALSE;
    rshift = 0.;
    
    // Set seed
    gRandom->SetSeed(rseed);
    
    return true;
    
}

void SK2p2MeV_mc::Analyze (long entry, bool last_entry)
{
    // print output info
    // =================
    if ( verbosity == 2 ) std::cout << "processing " << entry << " th event" << std::endl;
    
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
    
    // Bad channels
    int imask = 0;
    int istat;
    if ((prevrun != HEADER->nrunsk || prevsub != HEADER->nsubsk) && HEADER->nrunsk > 60000){
        prevrun = HEADER->nrunsk;
        prevsub = HEADER->nsubsk;
        int nsub = HEADER->nsubsk < 0 ? 1 : HEADER->nsubsk;
        skbadopt_(&imask);
        skbadch_(&(HEADER->nrunsk), &nsub, &istat);
        if(verbosity == 2) std::cout<<"skbadch: istat = "<<istat<<std::endl;
        /*
*        istat  ; -1 : error
*        istat  ;  0 : normal end  no modification
*        istat  ;  1 : normal end  read only badch.XXXXXX
*        istat  ;  2 : normal end  read only badch.XXXXXX.XXXXXX
*        istat  ;  3 : normal end  read only /skam/const/badch.dat
*        istat  ;  4 : normal end  read only badch2/badch.XXXXXX
*        istat  ;  5 : normal end  read only badch2/badch.XXXX00.XXXX99
*        istat  ;+10 : normal end  additional read /skam/const/badch.dat
         */
        SetBadch (combad_.nbad, combad_.isqbad);
    }
    
    // Set primary vertex
    if (useMC){
        VX = MC->pvtxvc[0][0];
        VY = MC->pvtxvc[0][1];
        VZ = MC->pvtxvc[0][2];
    }
    else{
        VX = LOWE->bsvertex[0];
        VY = LOWE->bsvertex[1];
        VZ = LOWE->bsvertex[2];
    }
    std::cout << VX << " " << VY << " " << VZ << std::endl;
    // Fiducial volume check
    if (sqrt(VX*VX + VY*VY) > 1690 || fabs(VZ) > 1810) return;
    if ( useMC && fSmear ) {
        // Smear the true vertex. 
        Double_t dx, dy, dz;
        do {
            dx = gRandom->Gaus(VX, sigma_xyz[0]);
            dy = gRandom->Gaus(VY, sigma_xyz[1]);
            dz = gRandom->Gaus(VZ, sigma_xyz[2]);
        } while ( sqrt(dx*dx+dy*dy) > 1690 || fabs(dz) > 1810 ); // fiducial check
        VX = dx;
        VY = dy;
        VZ = dz;
    }
       // END CHANGE!!
    
    if ( fShift ) {
        // Shift rshift from the true vertex. 
        Double_t dx, dy, dz;
        do {
            gRandom->Sphere(dx, dy, dz, rshift);
        } while ( sqrt((VX+dx)*(VX+dx)+(VY+dy)*(VY+dy)) > 1640 || fabs(VZ+dz) > 1760 ); // inside detector
        std::cout << VX << " " << VY << " " << VZ << std::endl;
        VX += dx;
        VY += dy;
        VZ += dz;
        std::cout << VX << " " << VY << " " << VZ << std::endl;
    }
    
    // Clear previous results
    res.Clear();
    
    // New results
    *head0 = *HEADER;
    *lowe0 = *LOWE;
    //    *third0 = *thirdred;
    smearedvertex[0]= VX;
    smearedvertex[1]= VY;
    smearedvertex[2]= VZ;
    
    // Check N200
    res.N200M = N200Max (18e3+1e3, 535e3+1e3);  //for 2.2 
    
    // Do 2.2MeV search
    NeutronSearch (18e3, 535e3 - cut_window); // for 2.2
    
    // Fill the output tree 
    theOTree->Fill();
    
}
