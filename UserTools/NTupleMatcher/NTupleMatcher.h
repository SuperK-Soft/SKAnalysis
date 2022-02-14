#ifndef NTUPLEMATCHER_HH
#define NTUPLEMATCHER_HH

#include "Tool.h"

class TFile;
class TTree;

class NTupleReader;

class NTupleMatcher : public Tool
{
    public:
        NTupleMatcher() { name = "NTupleMatcher"; }
        
        virtual bool Initialise(std::string configfile, DataModel &data);
        virtual bool Execute();
        virtual bool Finalise();
        
    private:
        std::string name;
        void GetNewEntry();
        void SetNTupleVariables();
        
        TFile* ntupleFile;
        TTree* ntuple;
        
        NTupleReader* ntupleReader;
        
        int iEntry, nEntries;
        int eventNo;
};

#endif
