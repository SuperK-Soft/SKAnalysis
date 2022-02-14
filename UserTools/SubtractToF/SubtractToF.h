#ifndef SUBTRACTTOF_HH
#define SUBTRACTTOF_HH

#include "Tool.h"

class SubtractToF : public Tool
{
    public:
        SubtractToF() { name = "SubtractToF"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();
        
   private:
        std::string name;
};

#endif
