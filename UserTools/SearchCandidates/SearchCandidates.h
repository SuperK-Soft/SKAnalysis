#ifndef SEARCHCANDIDATES_HH
#define SEARCHCANDIDATES_HH

#include "Tool.h"
#include "TVector3.h"

class SearchCandidates : public Tool
{
    public:
        SearchCandidates():
        T0TH(2e3), T0MX(535e3), TWIDTH(14), TMINPEAKSEP(0),
        NHITSTH(7), NHITSMX(70), N200TH(0), N200MX(200) 
        { name = "SearchCandidates"; }
    
        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();
        
    private:
        std::string name;
        TVector3 promptVertex;
        float T0TH, T0MX, TWIDTH, TMINPEAKSEP;
        int   NHITSTH, NHITSMX, N200TH, N200MX;
};

#endif
