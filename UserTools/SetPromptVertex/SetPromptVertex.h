#ifndef SETPROMPTVERTEX_HH
#define SETPROMPTVERTEX_HH

#define MU_RANGE_NBINS 15

#include <array>
#include <map>

#include "Tool.h"
#include "TVector3.h"

enum VertexMode
{
    mAPFIT,
    mBONSAI,
    mCUSTOM,
    mTRUE,
    mSTMU
};

class SetPromptVertex : public Tool
{
    public:
        SetPromptVertex()
        :vertexMode(mCUSTOM), customPromptVertex(0, 0, 0), smearingResolution(7.)
        { name = "SetPromptVertex"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();

    private:
        std::string name;
        VertexMode vertexMode;
        TVector3 customPromptVertex;
        float smearingResolution;

        std::array<float, MU_RANGE_NBINS> muMom;
        std::array<float, MU_RANGE_NBINS> muRange;

        std::map<std::string, VertexMode> cmdToModeMap;
};

#endif
