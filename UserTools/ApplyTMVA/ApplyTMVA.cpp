#include "ApplyTMVA.h"
#include "TMVA/Reader.h"

#include "PathGetter.h"

#include "Candidate.h"

bool ApplyTMVA::Initialise(std::string configfile, DataModel &data)
{
    
    if(configfile!="")  m_variables.Initialise(configfile);
    m_data= &data;
    m_data->tool_configs[name] = &m_variables;
    
    m_variables.Get("verbosity",m_verbose);
    // defaults
    likelihoodThreshold = 0.7f;
    mvaMethodName="MLP";
    weightFilePath = GetENV("NTAGPATH") + std::string("weights/MLP_Gd0.011p_calibration.xml");
    // update from config file
    m_variables.Get("likelihood_threshold", likelihoodThreshold);
    m_variables.Get("mva_method_name", mvaMethodName);
    m_variables.Get("weight_file_path", weightFilePath);

    featureContainer["NHits"] = 0;
    featureContainer["N200"] = 0;
    featureContainer["AngleMean"] = 0.;
    featureContainer["AngleSkew"] = 0.;
    featureContainer["AngleStdev"] = 0.;
    featureContainer["Beta1"] = 0.;
    featureContainer["Beta2"] = 0.;
    featureContainer["Beta3"] = 0.;
    featureContainer["Beta4"] = 0.;
    featureContainer["Beta5"] = 0.;
    featureContainer["DWall"] = 0.;
    featureContainer["DWall_n"] = 0.;
    featureContainer["DWallMeanDir"] = 0.;
    featureContainer["prompt_nfit"] = 0.;
    featureContainer["ThetaMeanDir"] = 0.;
    featureContainer["TRMS"] = 0.;


    tmvaReader = new TMVA::Reader();

    for (auto& pair: featureContainer)
        tmvaReader->AddVariable(pair.first, &(pair.second));

    tmvaReader->AddSpectator("CaptureType", &(captureType));

    tmvaReader->BookMVA(mvaMethodName, weightFilePath);
    
    m_data->eventCandidates.RegisterFeatureName("TMVAOutput");

    return true;
}

bool ApplyTMVA::Execute()
{
    int taggedNeutronCount = 0;

    EventCandidates* eventCans = &(m_data->eventCandidates);
    unsigned int nCandidates = eventCans->GetSize();

    // candidate loop
    for (unsigned int i = 0; i < nCandidates; i++) {
        Candidate* candidate = &(eventCans->At(i));
        float tmvaOutput = GetClassifierOutput(candidate);
        if (tmvaOutput > likelihoodThreshold) taggedNeutronCount++;
        candidate->Set("TMVAOutput", tmvaOutput);
    }

    if(m_verbose>2){
        Log("Candidates information:");
        eventCans->Print();
    }
    
    m_data->eventVariables.Set("tagged_neutron_count", taggedNeutronCount);
    std::cout << std::endl;
    Log(Form("%d neutron-like signals tagged in this event.", taggedNeutronCount));

    return true;
}

bool ApplyTMVA::Finalise()
{
    return true;
}

float ApplyTMVA::GetClassifierOutput(Candidate* candidate)
{
    // get features from candidate and fill feature container
    for (auto const& pair: featureContainer) {
        float value = candidate->Get(pair.first);
        featureContainer[pair.first] = value;
    }

    // get spectator
    captureType = candidate->Get("CaptureType");

    return tmvaReader->EvaluateMVA(mvaMethodName);
}
