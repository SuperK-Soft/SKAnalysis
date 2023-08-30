#ifndef READHITS_HH
#define READHITS_HH

#include "Tool.h"

class ReadHits : public Tool
{
    enum TrigType
    {
        tELSE,
        tSHE,
        tAFT
    };
    
    public:
        ReadHits(): prevEvTrigType(tELSE), currentEvTrigType(tELSE), isPrevEvProcessed(true)
        { name = "ReadHits"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();

    private:
        std::string name;
        void SetTriggerType();
        TrigType prevEvTrigType;
        TrigType currentEvTrigType;
        bool isPrevEvProcessed;

};

#endif
