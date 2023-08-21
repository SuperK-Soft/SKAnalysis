#include "NTupleMatcher.h"
#include "NTupleReader.h"

#include "TFile.h"
#include "TTree.h"

#include "skheadC.h"

bool NTupleMatcher::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	m_data= &data;
	m_data->tool_configs[name] = &m_variables;
    iEntry = 0; nEntries = 0; eventNo = 1;

    // get ntuple file name and file id
    std::string ntupleFilePath; int fileID;
    m_variables.Get("ntuple_file_path", ntupleFilePath);
    m_variables.Get("mc_file_id", fileID);
    
    // get ntuple and set branch address of event number
    ntupleFile = TFile::Open(ntupleFilePath.c_str());
    ntuple = (TTree*)ntupleFile->Get("h1");
    //ntuple->SetBranchAddress("nev", &eventNo);
    nEntries = ntuple->GetEntries();
    ntupleReader = new NTupleReader(ntuple);
    ntupleReader->nev = 1;
    
    // fill the first event number
    int prevEventNo = -1;
    
    // skip file
    for (int i=0; i<fileID; i++) {
        while (prevEventNo <= eventNo) {
            prevEventNo = eventNo;
            GetNewEntry();
        }
        prevEventNo = eventNo;
    }
    Log(Form("NTuple iEntry: %d, nEntries: %d", iEntry, nEntries));
    std::cout << "\n";
    
    return true;
}

bool NTupleMatcher::Execute()
{
    // get event number from skread
    Log(Form("skhead_.nevsk: %d, eventNo: %d", skhead_.nevsk, eventNo));
    if (skhead_.nevsk < eventNo) {
        Log(Form("No matching event in the ntuple (skhead_.nevsk %d, eventNo: %d). Skipping...", skhead_.nevsk, eventNo));
        m_data->vars.Set("Skip",true);
    }
    else if (skhead_.nevsk == eventNo) {
        SetNTupleVariables();
        Log(Form("A matching event exists for event %d in the ntuple! Continuing toolchain execution...", eventNo));
        //try { GetNewEntry(); }
        //catch(...){}
        //catch (ExceptionBehavior& e) { throw e; }
        GetNewEntry();
    }
    else {
        Log(Form("skhead_.nevsk %d is larger than eventNo %d! Ending toolchain execution...", skhead_.nevsk, eventNo));
        m_data->vars.Set("Skip",true);
        m_data->vars.Set("StopLoop",true);
    }

    return true;
}

bool NTupleMatcher::Finalise()
{
    ntupleFile->Close();
    delete ntupleReader;
    return true;
}

void NTupleMatcher::GetNewEntry()
{
    iEntry++;
    if (iEntry < nEntries) {
        ntupleReader->GetEntry(iEntry);
        eventNo = ntupleReader->nev;
    }
    else {
        Log("Reached the end of ntuple!");
        m_data->vars.Set("Skip",true);
        m_data->vars.Set("StopLoop",true);
    }
}

void NTupleMatcher::SetNTupleVariables()
{
    m_data->eventVariables.Set("nev", eventNo);
    m_data->eventVariables.Set("nring", ntupleReader->nring);
    m_data->eventVariables.Set("nhitac", ntupleReader->nhitac);
    m_data->eventVariables.Set("evis", ntupleReader->evis);
    m_data->eventVariables.Set("wall", ntupleReader->wall);
    m_data->eventVariables.Set("ip", ntupleReader->ip[0]);
    m_data->eventVariables.Set("ipnu", ntupleReader->ipnu[0]);
    m_data->eventVariables.Set("amome", ntupleReader->amome[0]);
    m_data->eventVariables.Set("amomm", ntupleReader->amomm[0]);
    m_data->eventVariables.Set("potot", ntupleReader->potot);
    m_data->eventVariables.Set("nn", ntupleReader->ntag_nn);
    m_data->eventVariables.Set("mctruth_nn", ntupleReader->ntag_mctruth_nn);
    m_data->eventVariables.Set("ndcy", ntupleReader->ndcy);
    m_data->eventVariables.Set("nmue", ntupleReader->nmue);
}
