#ifndef APPLYTMVA_HH
#define APPLYTMVA_HH

#include "Tool.h"

#include <memory>

#include "TMVA/Reader.h"

class Candidate;

class ApplyTMVA : public Tool
{
    public:
        ApplyTMVA() { name = "ApplyTMVA"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();

        float GetClassifierOutput(Candidate* candidate);

    private:
        std::string name;
        std::map<std::string, float> featureContainer;
        int captureType;
        float likelihoodThreshold;

        std::string mvaMethodName;
        std::string weightFilePath;
        TMVA::Reader* tmvaReader=nullptr;
};

#endif
