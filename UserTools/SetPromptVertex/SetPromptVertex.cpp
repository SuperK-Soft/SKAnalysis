#include "SetPromptVertex.h"
#include <TRandom.h>

#include <skroot.h>
#undef MAXHWSK
#include <apmringC.h>
#include <apmueC.h>
#include <apmsfitC.h>
#include <appatspC.h>
#include <loweroot.h>
#include <geotnkC.h>
#include <skvectC.h>
#undef MAXPM
#undef MAXPMA

#include "SKLibs.h"

#include "Calculator.h"

bool SetPromptVertex::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	m_data= &data;
	m_data->tool_configs[name] = &m_variables;
	m_variables.Get("m_verbose",m_verbose);
	
    std::string vertexModeCmd(""); 
    m_variables.Get("vertex_mode", vertexModeCmd);
    m_variables.Get("smearing_resolution", smearingResolution);

    cmdToModeMap["APFIT"] = mAPFIT;
    cmdToModeMap["BONSAI"] = mBONSAI;
    cmdToModeMap["CUSTOM"] = mCUSTOM;
    cmdToModeMap["TRUE"] = mTRUE;
    cmdToModeMap["STMU"] = mSTMU;

    if (cmdToModeMap.count(vertexModeCmd) > 0) {
        vertexMode = cmdToModeMap[vertexModeCmd];
        if (vertexMode == mCUSTOM){
            std::string customPromptVertexString;
            if(!m_variables.Get("prompt_vertex", customPromptVertexString)) {
                Log("Custom prompt vertex not specified!", pWARNING, m_verbose);
                Log("Aborting program...", pERROR, m_verbose);
                return false;
            } else {
                // parse the string and convert to TVector3
                char* pend=nullptr;
                double xtemp = strtod(customPromptVertexString.c_str(),&pend);
                double ytemp = strtod(++pend,&pend);
                double ztemp = strtod(++pend,&pend);
                customPromptVertex = TVector3(xtemp,ytemp,ztemp);
            }
        }

        if (vertexMode == mSTMU) {
            // initialize PDG muon range if STMU
            muMom = {0., 339.6, 1301., 2103., 3604., 4604., 5605., 7105., 8105., 9105.,
                       10110., 12110., 14110., 17110., 20110.};
            muRange  = {0., 103.9, 567.8, 935.3, 1595., 2023., 2443., 3064., 3472., 3877.,
                       4279., 5075., 5862., 7030., 8183.};
        }

        return true;
    }
    else {
        std::string allCmdKeys;
        for (auto const& pair: cmdToModeMap)
            allCmdKeys += (", " + pair.first);

        Log("Vertex mode not correctly specified. Specified command: " + vertexModeCmd, pWARNING, m_verbose);
        Log("Please specify correct vertex mode from " + allCmdKeys, pWARNING, m_verbose);
        Log("Aborting program...", pERROR, m_verbose);
        return false;
    }

    return true;
}

bool SetPromptVertex::Execute()
{
    TVector3 promptVertex;

    switch (vertexMode) {
        case mAPFIT: {
            int bank = 0;
            aprstbnk_(&bank);
            promptVertex = TVector3(apcommul_.appos); break;
        }
        case mBONSAI: {
            int lun = 10;
            TreeManager* mgr  = skroot_get_mgr(&lun);
            LoweInfo*    LOWE = mgr->GetLOWE();
            mgr->GetEntry();
            promptVertex = TVector3(LOWE->bsvertex); break;
        }
        case mCUSTOM: {
            promptVertex = customPromptVertex; break;
        }
        case mTRUE: {
            skgetv_();
            float dx=0, dy=0, dz=0;
            if(smearingResolution>0){
                    dx = 2*RINTK, dy = 2*RINTK, dz = 2*ZPINTK;
                    while (Norm(dx, dy, dz) > 150.) {
                        dx = gRandom->BreitWigner(0, smearingResolution);
                        dy = gRandom->BreitWigner(0, smearingResolution);
                        dz = gRandom->BreitWigner(0, smearingResolution);
                    }
            }
            promptVertex = TVector3(skvect_.pos) + TVector3(dx, dy, dz); break;
        }
        case mSTMU: {
            float initPoint[3], momDir[3];
            float muonMom;
            float goodness, entryQ; // probably dummy
            int iRing = 0; // first ring
            enum {iGamma, iElectron, iMuon};
            stmfit_(initPoint, momDir, goodness, entryQ);
            if (goodness < 0)
                Log("STMUFIT error occurred.", pWARNING, m_verbose);

            // 1-ring muon
            apcommul_.apnring = 1; apcommul_.apip[iRing] = 13;
            for (int dim = 0; dim < 3; dim++) {
                apcommul_.appos[dim] = initPoint[dim];
                apcommul_.apdir[iRing][dim] = momDir[dim];
            }

            // Fix Cherenkov angle to 42 deg
            apcommul_.apangcer[iRing] = 42.;
            int iCall = 0, iTrCor = 0, iPAng = 0; // probably dummy
            int nRing = 1;
            sparisep_(iCall, iTrCor, iPAng, nRing); // momentum

            // MS-fit
            pffitres_.pffitflag = 1; // 0: normal fit, 1: fast fit
            int iPID = iMuon + 1; // muon (+1 for fortran)
            pfdodirfit_(iPID); // direction fit

            // Replace APFit direction with MS-fit direction
            for (int dim = 0; dim < 3; dim++)
                apcommul_.apdir[iRing][dim] = pffitres_.pfdir[2][iRing][dim];
            sparisep_(iCall, iTrCor, iPAng, nRing); // momentum
            sppang_(apcommul_.apip[iRing], apcommul_.apamom[iRing], apcommul_.apangcer[iRing]); // Cherenkov angle
            // Calculate momentum
            appatsp_.approb[iRing][iElectron] = -100.; // e-like probability
            appatsp_.approb[iRing][iMuon]     = 0.;    // mu-like probability (set this higher than e-like!)
            spfinalsep_();

            for (int dim = 0; dim < 3; dim++)
                momDir[dim] = apcommul_.apdir[iRing][dim];

            // final mu-like momentum
            muonMom = appatsp2_.apmsamom[iRing][iMuon];

            // Muon range as a function of momentum (almost linear, PDG2020)
            int iBin = 0;
            float range = 0;
            for (iBin = 0; muMom[iBin] < muonMom && iBin < MU_RANGE_NBINS; iBin++);
            if (iBin < MU_RANGE_NBINS)
                // linear interpolation for muon momentum within data
                range = muRange[iBin-1]
                        + (muRange[iBin]-muRange[iBin-1]) * (muonMom-muMom[iBin-1])/(muMom[iBin]-muMom[iBin-1]);
            else
                // for muon momentum larger than 20 GeV/c,
                // assume constant stopping power 2.3 MeV/cm
                range = muonMom / 2.3;
            TVector3 stopPoint;
            stopPoint = TVector3(initPoint) + range * TVector3(momDir);
            // for estimated stop point outside tank, force it to be within tank
            float r = stopPoint.Perp();
            if (r > RINTK)
                stopPoint *= RINTK/r;
            float z = abs(stopPoint.z());
            if (z > ZPINTK)
                stopPoint *= ZPINTK/z;
            Log(Form("Initial position : (%3.0f, %3.0f, %3.0f) cm", initPoint[0], initPoint[1], initPoint[2]));
            Log(Form("Momentum direction : (%3.3f, %3.3f, %3.3f)", momDir[0], momDir[1], momDir[2]));
            Log(Form("Momentum : %3.0f MeV/c", muonMom));
            Log(Form("Stopping point : (%3.0f, %3.0f, %3.0f) cm", stopPoint.x(), stopPoint.y(), stopPoint.z()));

            promptVertex = stopPoint; break;
        }
    }

    float pv[3] = {(float)promptVertex.x(), (float)promptVertex.y(), (float)promptVertex.z()};
    float dWall = wallsk_(pv);

    Log(Form("Prompt vertex: (%3.0f, %3.0f, %3.0f) cm", promptVertex.x(), promptVertex.y(), promptVertex.z()));
    Log(Form("dwall: %3.0f cm", dWall));
    m_data->eventVariables.Set("prompt_vertex", promptVertex);
    m_data->eventVariables.Set("d_wall", dWall);

    return true;
}

bool SetPromptVertex::Finalise()
{
    return true;
}
