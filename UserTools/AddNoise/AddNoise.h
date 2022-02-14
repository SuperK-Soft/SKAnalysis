#ifndef ADDNOISE_HH
#define ADDNOISE_HH

#include "Tool.h"

#include "TChain.h"


class TChain;
class TQReal;

class AddNoise : public Tool
{
    public:
        AddNoise(): noiseChain("data")
        { name = "AddNoise"; }
        
        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();
        
    private:
        std::string name;
        void GetNewEntry();
        void SetNoiseHitCluster();
        
        TChain noiseChain;
        
        TQReal* tqi;
        
        float noiseEventLength;
        float noiseStartTime, noiseEndTime, noiseTimeWindowWidth;
        
        int iEntry, nEntries;
        int iPart, nParts;
        
        unsigned int iHit, nHits;
        float noiseT0;
        std::vector<float> t;
        std::vector<float> q;
        std::vector<int>   i;
        
        PMTHitCluster noiseEventHits;
};

#endif
