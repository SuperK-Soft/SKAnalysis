#ifndef WRITEOUTPUT_HH
#define WRITEOUTPUT_HH

class TFile;
class TTree;
class TBranch;

#include "Tool.h"

#include "EventParticles.h"
#include "EventTrueCaptures.h"

#include "TString.h"
#include "TVector3.h"
#include "TSysEvtHandler.h"

class TInterruptHandler : public TSignalHandler
{
   public:
        TInterruptHandler(Tool* tool):TSignalHandler(kSigInterrupt, kFALSE)
        { WriteOutput = tool; }

        virtual Bool_t Notify()
        {
            std::cerr << "Received SIGINT. Writing output..." << std::endl;
            WriteOutput->Finalise();
            _exit(2);
            return kTRUE;
        }

    private:
        Tool* WriteOutput;
};

class WriteOutput : public Tool
{
    public:
        WriteOutput() { name = "WriteOutput"; }

        virtual bool Initialise(std::string configfile, DataModel &data);
        virtual bool Execute();
        virtual bool Finalise();

    private:
        std::string name;
        void CreateTrees();
        bool CheckTreesExist();
        void GetTrees();
        void PrintTrees();
        void WriteTrees(int option);

        void MakeParticleTrees(TTree* particleTree, EventParticles* eventParticles);
        void MakeTrueCaptureTrees(TTree* particleTree, EventTrueCaptures* eventTrueCaptures);
    
        TFile* outFile;
        TTree* variableTree;
        TTree* candidateTree;
        bool firstWrite;

        // MC
        bool inputIsMC;
        TTree* mcTree;
        
        // NTag
        TTree* ntagInfoTree;
        
        TString outputMode;
        
        EventParticles* primaries;
        EventParticles* secondaries;
        EventTrueCaptures* trueCaptures;
        
        TInterruptHandler* handler;
        
        // verbosity levels: if 'verbosity' < this level, the message type will be logged.
        int v_error=0;
        int v_warning=1;
        int v_message=2;
        int v_debug=3;
};

#endif
