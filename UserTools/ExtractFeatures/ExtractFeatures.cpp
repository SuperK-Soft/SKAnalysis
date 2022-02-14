#include "ExtractFeatures.h"

#include "skheadC.h"
#include "geotnkC.h"

#include "Calculator.h"
#include "SK_helper_functions.h"

bool ExtractFeatures::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	m_data= &data;
	m_data->tool_configs[name] = &m_variables;
	m_variables.Get("tool_verbosity",m_verbose);
	
    // read ntag options
    if (!m_variables.Get("TWIDTH", tWidth))
        Log("TWIDTH not specified!", pERROR,m_verbose);
    // TODO if the TWIDTH value needs to match that loaded by the SearchCandidates tool,
    // it's probably better for the SearchCandidates tool to put this into m_vata->vars
    // and have this tool read that, so that it's only specified in one place. Currently
    // we don't have 'global' configuration variables, but they might be something we'll
    // add as they seem potentially useful

    // read trms-fit options
    if (!m_variables.Get("INITGRIDWIDTH", initGridWidth)) initGridWidth = 800;
    if (!m_variables.Get("MINGRIDWIDTH", minGridWidth)) minGridWidth = 50;
    if (!m_variables.Get("GRIDSHRINKRATE", gridShrinkRate)) gridShrinkRate = 0.5;
    if (!m_variables.Get("VTXSRCRANGE", vertexSearchRange)) vertexSearchRange = 5000;

    // read MC true capture match options
    m_data->vars.Get("inputIsMC",inputIsMC);
    if (inputIsMC) {
        if (!m_variables.Get("TMATCHWINDOW", tMatchWindow)) tMatchWindow = 50;
    }
    
    // define key list to branch
    // TODO again this needs to match the list of features in the ApplyTMVA Tool,
    // so we should ensure that somehow.
    std::vector<std::string>keyList = {"NHits", "N50", "N200", "N1300", "ReconCT", "TRMS", "QSum", 
                                       "Beta1", "Beta2", "Beta3", "Beta4", "Beta5", 
                                       "AngleMean", "AngleSkew", "AngleStdev", "CaptureType",
                                       "DWall", "DWallMeanDir", "ThetaMeanDir", "DWall_n", "prompt_nfit",
                                       "decay_e_like", "TrmsFitVertex_X", "TrmsFitVertex_Y", "TrmsFitVertex_Z"};
    
    // register keys to the feature map
    m_data->eventCandidates.RegisterFeatureNames(keyList);

    return true;
}

bool ExtractFeatures::Execute()
{
    // EventPMTHits must be filled
    if (m_data->eventPMTHits.IsEmpty()) {
        Log("EventPMTHits class is empty. Skipping...", pWARNING,m_verbose);
        //m_data->vars.Set("Skip",true);
        return true;
    }
    else if (m_data->eventCandidates.IsEmpty()) {
        Log("No candidates found in the event! Skipping...", pWARNING,m_verbose);
        //m_data->vars.Set("Skip",true);
        return true;
    }
    
    TVector3 promptVertex(0,0,0);
    float dWall = 0.0;
    m_data->eventVariables.Get("prompt_vertex", promptVertex);
    m_data->eventVariables.Get("d_wall", dWall);

    EventPMTHits* eventHits = &(m_data->eventPMTHits);
    EventCandidates* eventCans = &(m_data->eventCandidates);
    EventTrueCaptures* eventCaps = &(m_data->eventTrueCaptures);
    EventParticles* eventSecs = &(m_data->eventSecondaries);

    unsigned int nCandidates = eventCans->GetSize();

    if (!eventHits->HasVertex())
        eventHits->SetVertex(promptVertex);

    // loop over candidates
    for (unsigned int i = 0; i < nCandidates; i++) {
        Candidate* candidate = &(eventCans->At(i));
        int firstHitID = candidate->HitID();

        PMTHitCluster hitsInTWIDTH = eventHits->Slice(firstHitID, tWidth);
        PMTHitCluster hitsIn50ns   = eventHits->Slice(firstHitID, tWidth/2.- 50, tWidth/2.+ 50);
        PMTHitCluster hitsIn200ns  = eventHits->Slice(firstHitID, tWidth/2.-100, tWidth/2.+100);
        PMTHitCluster hitsIn1300ns = eventHits->Slice(firstHitID, tWidth/2.-520, tWidth/2.+780);

        // Number of hits
        candidate->Set("NHits", hitsInTWIDTH.GetSize());
        candidate->Set("N50",   hitsIn50ns.GetSize());
        candidate->Set("N200",  hitsIn200ns.GetSize());
        candidate->Set("N1300", hitsIn1300ns.GetSize());

        // Time
        float reconCT = hitsInTWIDTH.Find(HitFunc::T, Calc::Mean) * 1e-3;
        candidate->Set("ReconCT", reconCT);
        candidate->Set("TRMS", hitsInTWIDTH.Find(HitFunc::T, Calc::RMS));

        // Charge
        candidate->Set("QSum", hitsInTWIDTH.Find(HitFunc::Q, Calc::Sum));

        // Beta's
        std::array<float, 6> beta = hitsInTWIDTH.GetBetaArray();
        candidate->Set("Beta1", beta[1]);
        candidate->Set("Beta2", beta[2]);
        candidate->Set("Beta3", beta[3]);
        candidate->Set("Beta4", beta[4]);
        candidate->Set("Beta5", beta[5]);

        // DWall
        auto dirVec = hitsInTWIDTH[HitFunc::Dir];
        auto meanDir = GetMean(dirVec).Unit();
        candidate->Set("DWall", dWall);
        candidate->Set("DWallMeanDir", GetDWallInDirection(promptVertex, meanDir));

        // Mean angle formed by all hits and the mean hit direction
        std::vector<float> angles;
        for (auto const& dir: dirVec) {
            angles.push_back((180/M_PI)*meanDir.Angle(dir));
        }
        float meanAngleWithMeanDirection = GetMean(angles);
        candidate->Set("ThetaMeanDir", meanAngleWithMeanDirection);

        // Opening angle stats
        OpeningAngleStats openingAngleStats = hitsInTWIDTH.GetOpeningAngleStats();
        candidate->Set("AngleMean",  openingAngleStats.mean);
        candidate->Set("AngleStdev", openingAngleStats.stdev);
        candidate->Set("AngleSkew",  openingAngleStats.skewness);

        // TRMS-fit
        TVector3 trmsFitVertex = hitsInTWIDTH.FindTRMSMinimizingVertex(/* TRMS-fit options */
                                                                   initGridWidth, minGridWidth, gridShrinkRate, vertexSearchRange);
        candidate->Set("TrmsFitVertex_X", trmsFitVertex.X());
        candidate->Set("TrmsFitVertex_Y", trmsFitVertex.Y());
        candidate->Set("TrmsFitVertex_Z", trmsFitVertex.Z());
        candidate->Set("DWall_n", GetDWall(trmsFitVertex));
        candidate->Set("prompt_nfit", (promptVertex-trmsFitVertex).Mag());

        int passDecayECut = 0;
        if ((candidate->Get("N50") > 50) && reconCT < 20) {
            passDecayECut = 1;
        }
        candidate->Set("decay_e_like", passDecayECut);

        // BONSAI

        // MC info
        if (inputIsMC) {
            Log("MC!",5,m_verbose);
            // default: not a capture
            int captureType = 0;
            int nTrueCaptures = eventCaps->GetSize();
            int nSecondaries = eventSecs->GetSize();

            // label candidates as signal if a true capture with matching capture time exists
            for (int iCapture = 0; iCapture < nTrueCaptures; iCapture++) {
                TrueCapture& capture = eventCaps->At(iCapture);
                if (fabs(capture.Time() - reconCT*1e3) < tMatchWindow) {
                    if (capture.Energy() > 6.) captureType = 2; // Gd
                    else                       captureType = 1; // H
                }
            }
            
            // search for a matching decay electron and overwrite label
            for (int iSec = 0; iSec < nSecondaries; iSec++) {
                Particle& secondary = eventSecs->At(iSec);

                // decay electron from a parent muon
                if (fabs(secondary.PID()) == 11 &&
                    fabs(secondary.ParentPID()) == 13 &&
                    secondary.IntID() == 5) {
                    // check if ReconCT matches the decay within TMATCHWINDOW
                    if (fabs(secondary.Time() - reconCT*1e3) < tMatchWindow) {
                        captureType = 3; // decay electron
                    }
                }
            }
            candidate->Set("CaptureType", captureType);
        }
    }

    return true;
}

bool ExtractFeatures::Finalise()
{
    return true;
}
