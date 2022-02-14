#include "ReadMCInfo.h"

#include <cmath>

#include <skroot.h>
#undef MAXHWSK
#undef MAXPM
#undef MAXPMA
#include "skvectC.h"
#include "geotnkC.h"
#include "skheadC.h"

#include "apscndryC.h"

#include "SKLibs.h"

#include "TrueCapture.h"
#include "Particle.h"
#include "MTreeReader.h"
#include "SuperWrapper.h"

bool ReadMCInfo::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    m_data= &data;
    m_data->tool_configs[name] = &m_variables;
    
    std::string treeReaderName;
    m_variables.Get("treeReaderName",treeReaderName);
    int LUN = m_data->GetLUN(treeReaderName);
    if(LUN<0){
        Log("ReadMCInfo failed to find TreeReader "+treeReaderName+" in DataModel",0,0);
        return false;
    }
    // use our hacked workaround to override the corresponding SuperManager method
    TreeManager* manager = GetTreeManager(LUN);
    // if there's a TreeManager associated with it, we're reading a ROOT file.
    // if not, must be reading a zbs file.
    inputIsSKROOT = (manager!=nullptr);
    
    //inputIsSKROOT = (inFilePath.find(".root")!=std::string::npos);
    //std::cout<<"ReadMCInfo will process file "<<inFilePath<<" with inputIsSKROOT ="<<inputIsSKROOT<<std::endl;
    return true;
}

bool ReadMCInfo::Execute()
{
    m_data->eventPrimaries.Clear();
    m_data->eventSecondaries.Clear();
    m_data->eventTrueCaptures.Clear();
    
    // double check we're processing MC
    bool isMC = false;
    m_data->vars.Get("inputIsMC",isMC);
    if(!isMC){
        Log("ReadMCInfo in toolchain but MC flag is not set in datamodel! Is the input file MC???",0,0);
        exit(-1);
        return true;
    }

    float geantT0;
    trginfo_(&geantT0);
    m_data->eventVariables.Set("geant_t0", geantT0);
    Log(Form("Simulation T0: %3.2f ns", geantT0));

    // Primaries
    skgetv_();

    for (int iVec = 0; iVec < skvect_.nvect; iVec++) {
        Particle primary(skvect_.ip[iVec], geantT0, TVector3(skvect_.pos), TVector3(skvect_.pin[iVec]));
TVector3 vp(skvect_.pos);
std::cout<<"primary vertex at "<<vp.X()<<", "<<vp.Y()<<", "<<vp.Z()<<std::endl;
        m_data->eventPrimaries.Append(primary);
    }

    // NEUT

    // Read secondary bank
    double vtxscale=1.;
    if (inputIsSKROOT) {
        int lun = 10;

        TreeManager* mgr  = skroot_get_mgr(&lun);
        SecondaryInfo* SECONDARY = mgr->GetSECONDARY();
        mgr->GetEntry();

        secndprt_.nscndprt = SECONDARY->nscndprt;

        std::copy(std::begin(SECONDARY->iprtscnd), std::end(SECONDARY->iprtscnd), std::begin(secndprt_.iprtscnd));
        std::copy(std::begin(SECONDARY->iprntprt), std::end(SECONDARY->iprntprt), std::begin(secndprt_.iprntprt));
        std::copy(std::begin(SECONDARY->lmecscnd), std::end(SECONDARY->lmecscnd), std::begin(secndprt_.lmecscnd));
        std::copy(std::begin(SECONDARY->tscnd), std::end(SECONDARY->tscnd), std::begin(secndprt_.tscnd));

        // is there a unit mismatch or sth? All the vertices from SKG4 are outside the FV....
        // FIXME: this probably needs to be fixed in SKG4, not here! Will probably BREAK skdetsim use! FIXME
        vtxscale=0.1;
        std::copy(&SECONDARY->vtxscnd[0][0], &SECONDARY->vtxscnd[0][0] + 3*SECMAXRNG, &secndprt_.vtxscnd[0][0]);
        std::copy(&SECONDARY->pscnd[0][0], &SECONDARY->pscnd[0][0] + 3*SECMAXRNG, &secndprt_.pscnd[0][0]);
    }
    else {
        apflscndprt_();
    }

    int nAllSec = secndprt_.nscndprt;

    // Loop over all secondaries in secondary common block
    for (int iSec = 0; iSec < nAllSec; iSec++) {
        TVector3 startvtx(secndprt_.vtxscnd[iSec]);
        startvtx *= vtxscale;
        std::cout<<"secondary at "<<startvtx.X()<<", "<<startvtx.Y()<<", "<<startvtx.Z()<<std::endl;

        Particle secondary(secndprt_.iprtscnd[iSec],
                           secndprt_.tscnd[iSec] + geantT0,
                           startvtx,
                           TVector3(secndprt_.pscnd[iSec]),
                           secndprt_.iprntprt[iSec],
                           secndprt_.lmecscnd[iSec]);

        // neutrons
        if (secondary.PID() == 2112)  // neutron
            m_data->eventSecondaries.Append(secondary);

        // deuteron, gamma, electrons
        else if (secondary.PID() == 100045 ||     // skdetsim deuteron
                 secondary.PID() == 1000010020 || // SKG4 deuteron
                 secondary.PID() == 1000641580 || // SKG4 Gd157
                 secondary.PID() == 100045 ||     // SKG4 Gd155
                 secondary.PID() == 22 ||         // gamma
                 // electrons over Cherenkov threshold momentum, from interaction other than multiple scattering
                 (fabs(secondary.PID()) == 11 && secondary.Momentum().Mag() > 0.579 && secondary.IntID() != 2)) {

            int inPMT; inpmt_(secndprt_.vtxscnd[iSec], inPMT);

            // Check if the capture is within ID volume
            if (secondary.Vertex().Perp() < RINTK && fabs(secondary.Vertex().z()) < ZPINTK && !inPMT) {

                // Save secondary (deuteron, gamma, electrons)
                m_data->eventSecondaries.Append(secondary);

                bool isNewCapture = true;

                // particle produced by n-capture
                if (secondary.IntID() == 18) {

                    // Check saved capture stack
                    auto nCaptures = m_data->eventTrueCaptures.GetSize();
                    for (int iCapture = 0; iCapture < nCaptures; iCapture++) {
                        TrueCapture* capture = &(m_data->eventTrueCaptures.At(iCapture));
                        if (fabs((double)(secondary.Time()-capture->Time()))<1.e-7) {
                            isNewCapture = false;
                            if (secondary.PID() == 22) {
                                capture->Append(secondary);
                            }
                        }
                    }

                    if (isNewCapture && secondary.PID() == 22) {
                        TrueCapture capture;
                        capture.Append(secondary);
                        m_data->eventTrueCaptures.Append(capture);
                    }
                }
            }
//        } else {
//            std::cout<<"Skipping non-neutron, non daughter secondary with pdg "<<secondary.PID()
//                     <<" and creation process "<<secondary.IntID()<<std::endl;
        }
    }
    Log("True capture information:");
    m_data->eventTrueCaptures.Sort();
    m_data->eventTrueCaptures.DumpAllElements();
    
    m_data->eventVariables.Set("true_neutron_count", m_data->eventTrueCaptures.GetSize());

    return true;
}

bool ReadMCInfo::Finalise()
{
    return true;
}
