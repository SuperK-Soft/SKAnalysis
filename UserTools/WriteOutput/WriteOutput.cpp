#include "WriteOutput.h"
#include <unistd.h>      // system
#include <cstdlib>       // getenv

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSysEvtHandler.h"

#include "Algorithms.h"  // getOutputFromFunctionCall
#include "Printer.h"     // PrintBlock


bool WriteOutput::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    m_data= &data;
    m_data->tool_configs[name] = &m_variables;
    m_variables.Get("verbosity",m_verbose);
    
    // get info about the current application git revision
    std::string gitcommit = getOutputFromFunctionCall(system, "git rev-parse HEAD");
    gitcommit.pop_back();      // remove trailing newline
    m_variables.Set("ntag_commit_id", gitcommit);
    
    // if being verbose, print some additional info
    if(m_verbose>1){
        std::string gitcommitdate = getOutputFromFunctionCall(system, "git log -1 --format=%cd");
        gitcommitdate.pop_back();  // remove trailing newline
        PrintNTag(gitcommit.c_str(), gitcommitdate.c_str());
    }
    if(m_verbose>2){
        char* installPath = getenv("NTAGPATH");
        std::string empty = "";
        if(installPath==nullptr) installPath = const_cast<char*>(empty.c_str());
        char* pwd = getenv("PWD");
        if (strcmp(pwd,installPath)!=0) {
            PrintBlock("NTAGPATH");
            std::cout << installPath << std::endl;
        }
    }
    // pass on high verbosity to the StoreConverter for debugging
    if(m_verbose>2) m_data->StoreConverter.SetVerbosity(10);
    
    // XXX ???
    TClass::GetClass("TVector3")->IgnoreTObjectStreamer();
    
    // make (or open) output file
    TString outFilePath;
    m_variables.Get("output_file_path", outFilePath);
    m_variables.Get("output_mode", outputMode);
    Log("Output mode: " + outputMode);
    outFile = new TFile(outFilePath, outputMode);
    
    // check if we need to create the MC branches
    m_data->vars.Get("inputIsMC",inputIsMC);
    
    // create or retrieve the output TTrees
    if (outputMode == "update") {
        if (CheckTreesExist()) {
            Log("Trees exist!");
            GetTrees();
        }
        else {
            Log("Trees don't exist!");
            CreateTrees();
            outputMode = "recreate";
        }
    }
    else
        CreateTrees();
    
    // if processing MC, get pointers to the MC truth info
    if (inputIsMC) {
        Log("MC!",v_debug,m_verbose);
        primaries = &(m_data->eventPrimaries);
        secondaries = &(m_data->eventSecondaries);
        trueCaptures = &(m_data->eventTrueCaptures);
    } else {
        Log("Data!",v_debug,m_verbose);
    }
    
    // the ntagInfo tree stores the ToolChain configuration scraped from tool config files
    // this may mean we have variables we do not need to save (e.g. verbose/verbosity)
    // and there may also be duplicates, which will be ambiguous in the output file
    // (and will result in multiple TTree Fill calls when there should only be one)
    // so consolidate our data here.
    Store all_variables;
    std::map<std::string,std::string> duplicates;
    // loop over tools
    for(std::pair<const std::string, Store*> atoolconfig : m_data->tool_configs){
        // loop over tool config variables
        for (const std::pair<const std::string, std::string>& avariable : *(atoolconfig.second->GetMap())){
            std::string key = avariable.first;
            if(key=="verbose" || key=="verbosity" || key=="tool_verbosity") continue;
            std::string value = avariable.second;
            if(all_variables.Has(key)){
                // found a duplicate. remove it for now.
                all_variables.Erase(key);
                // instead add it to the list of duplicates
                duplicates.emplace(atoolconfig.first, key);
            } else {
                // unique so far
                all_variables.Set(key,value);
            }
        }
    }
    // now let's re-add any duplicates, this time making their name unique
    // by prepending the tool name
    for(std::pair<const std::string, std::string>& avariable : duplicates){
        std::string unique_key = avariable.first + "_" + avariable.second;
        // get the value from the appropriate store
        std::string value;
        m_data->tool_configs.at(avariable.first)->Get(avariable.second,value);
        // add to our full list
        all_variables.Set(unique_key, value);
    }
    
    m_data->StoreConverter.MakeBranches(ntagInfoTree, &all_variables);
    m_data->StoreConverter.FillBranches(ntagInfoTree, &all_variables);
    
    // some of the other output trees have a dynamic structure that depends on
    // the contents of DataModel stores. Those won't be filled yet,
    // so we'll need to delay initialization of those trees until the first write call.
    firstWrite=true; // flag for setting up output TTree branches

    return true;
}

bool WriteOutput::Execute()
{
    
    m_data->eventCandidates.FillVectorMap();
    
    if(firstWrite){
        // only setup TTree branches on first write
        // we need to delay at least some of these to the first Execute call
        // as the branches are based on Stores that aren't filled until then.
        if (outputMode == "update") {
            // updating an existing file, branches already exist.
            for (auto& pair: m_data->eventCandidates.featureVectorMap)
                candidateTree->SetBranchAddress(pair.first.c_str(), &(pair.second));
            if (inputIsMC) {
                mcTree->SetBranchAddress("primaries", &primaries);
                mcTree->SetBranchAddress("secondaries", &secondaries);
                mcTree->SetBranchAddress("captures", &trueCaptures);
            }
            // variableTree branch address setting is handled by the StoreConverter
        } else {
            // make the TTree branch structure
            m_data->StoreConverter.MakeBranches(variableTree, &(m_data->eventVariables));
            for (auto& pair: m_data->eventCandidates.featureVectorMap)
                candidateTree->Branch(pair.first.c_str(), &(pair.second));
            if (inputIsMC) {
                Int_t bufsize = 32000;
                mcTree->Branch("primaries", &primaries,bufsize,0);
                mcTree->Branch("secondaries", &secondaries, bufsize, 3);
                // FIXME - is this related to IgnoreTObjectStreamer set above...??
                // There's a bizarre issue here: the written TVector3 of trueCaptures gets corrupted!
                // e.g. a TrueCapture with vertex (33,55,77) becomes a file entry with vertex (55,77,00).
                // This doesn't happen with splitlevel=0.
                // meanwhile setting splitlevel=2 or less for the 'secondaries' branch results in
                // a segfault when invoking mcTree->Fill()???
                mcTree->Branch("captures", "EventTrueCaptures", &trueCaptures,bufsize,0);
            }
        }
        firstWrite=false;
    }
    
    candidateTree->Fill();
    
    m_data->StoreConverter.FillBranches(variableTree, &(m_data->eventVariables));
    if(m_verbose>2) m_data->eventVariables.Print();
    
    if (inputIsMC) {
        if(m_verbose>2){
            primaries->DumpAllElements();
            secondaries->DumpAllElements();
            trueCaptures->DumpAllElements();
        }
        
        mcTree->Fill();
    }
    
    return true;
}

bool WriteOutput::Finalise()
{
    if (m_verbose > pDEFAULT) PrintTrees();
    
    outFile->cd();
    
    int option = outputMode == "update" ? TObject::kWriteDelete : TObject::kOverwrite;
    WriteTrees(option);
    
    outFile->Close();
    
    return true;
}

void WriteOutput::CreateTrees()
{
    outFile->cd();
    
    variableTree = new TTree("variables", "Event variables");
    candidateTree = new TTree("candidates", "Neutron-signal candidate features");
    
    // NTag
    ntagInfoTree = new TTree("ntaginfo", "NTag options and information");
    
    // MC
    if (inputIsMC)
        mcTree = new TTree("mc", "MC truth information");
}

bool WriteOutput::CheckTreesExist()
{
    bool variableTreeExists = (outFile->Get("variables") != nullptr);
    bool candidateTreeExists = (outFile->Get("candidates") != nullptr);
    bool ntagInfoTreeExists = (outFile->Get("ntaginfo") != nullptr);
    bool mcTreeExists = (outFile->Get("mc") != nullptr);

    
    bool dataTreesExist = variableTreeExists && candidateTreeExists && ntagInfoTreeExists;
    
    if (!inputIsMC)
        return dataTreesExist;
    else
        return dataTreesExist && mcTreeExists;
}

void WriteOutput::GetTrees()
{
    variableTree = (TTree*)outFile->Get("variables");
    candidateTree = (TTree*)outFile->Get("candidates");
     
    // NTag
    ntagInfoTree = (TTree*)outFile->Get("ntaginfo");
    
    // MC
    if (inputIsMC)
        mcTree = (TTree*)outFile->Get("mc");
    
}

void WriteOutput::PrintTrees()
{
    ntagInfoTree->Print();
    variableTree->Print();
    candidateTree->Print();
    
    if (inputIsMC)
        mcTree->Print();
}

void WriteOutput::WriteTrees(int option)
{
    ntagInfoTree->Write(0, option, 0);
    variableTree->Write(0, option, 0);
    candidateTree->Write(0, option, 0);
    
    if (inputIsMC)
        mcTree->Write(0, option, 0);
}

