#ifndef EXTRACTFEATURES_HH
#define EXTRACTFEATURES_HH

#include "Tool.h"

class ExtractFeatures : public Tool
{
    public:
        ExtractFeatures():
        tWidth(14), tMatchWindow(50),
        initGridWidth(800), minGridWidth(50), gridShrinkRate(0.5), vertexSearchRange(5000)
        { name = "ExtractFeatures"; }

        bool Initialise(std::string configfile, DataModel &data);
        bool Execute();
        bool Finalise();
	
    private:
        std::string name;
        float tWidth;
        float tMatchWindow;
        float initGridWidth, minGridWidth, gridShrinkRate, vertexSearchRange;
        
        bool inputIsMC;
};

#endif
