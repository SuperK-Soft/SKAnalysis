#ifndef TRUECAPTURE_HH
#define TRUECAPTURE_HH

#include "Particle.h"

class TrueCapture
{
    public:
        TrueCapture();
        virtual ~TrueCapture();
        void Append(const Particle& particle);

        inline float Time() const { return t; }
        inline float Energy() const { return E; }
        inline TVector3 Vertex() const { return v; }

        void Dump();

    private:
        TVector3 v;
        float t;
        float E;
        int nGamma;
        
    ClassDef(TrueCapture, 1)
};

#endif
