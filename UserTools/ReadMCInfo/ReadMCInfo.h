#ifndef READMCINFO_HH
#define READMCINFO_HH

#include "Tool.h"

class ReadMCInfo : public Tool
{
    public:
        ReadMCInfo() { name = "ReadMCInfo"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();

    private:
        std::string name;
        bool inputIsSKROOT;
};

#endif
