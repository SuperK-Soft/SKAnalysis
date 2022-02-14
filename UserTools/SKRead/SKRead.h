#ifndef SKREAD_HH
#define SKREAD_HH

#include "Tool.h"

enum ReadStatus
{
    readOK,
    readError,
    readEOF
};

class SKRead : public Tool
{
    public:
        SKRead() { name = "SKRead"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();

        inline int GetReadStatus() { return readStatus; }

    private:
        std::string name;
        bool inputIsSKROOT;
        int readStatus;
        int exeCounter=0;
};

#endif
