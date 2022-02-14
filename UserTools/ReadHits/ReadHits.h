#ifndef READHITS_HH
#define READHITS_HH

#include "Tool.h"

enum TriggerType
{
    tELSE,
    tSHE,
    tAFT
};

class ReadHits : public Tool
{
    public:
        ReadHits(): prevEvTrigType(tELSE), currentEvTrigType(tELSE), isPrevEvProcessed(true)
        { name = "ReadHits"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();

    private:
        std::string name;
        void SetTriggerType();
        TriggerType prevEvTrigType;
        TriggerType currentEvTrigType;
        bool isPrevEvProcessed;

};

#endif
