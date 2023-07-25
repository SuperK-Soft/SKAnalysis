#ifndef EVENTPARTICLES_HH
#define EVENTPARTICLES_HH

#include "Particle.h"
#include "Cluster.h"

class EventParticles : public Cluster<Particle> 
{
    public:
        EventParticles() {}
        virtual ~EventParticles();
        
        void DumpAllElements();
    
    ClassDef(EventParticles, 1)
};

#endif
