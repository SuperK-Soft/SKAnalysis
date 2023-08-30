#ifndef EVENTTRUECAPTURES_HH
#define EVENTTRUECAPTURES_HH

#include <memory>

#include "Cluster.h"
#include "TrueCapture.h"

class EventTrueCaptures : public Cluster<TrueCapture>
{
    public:
        virtual ~EventTrueCaptures();
        void Sort();
        void DumpAllElements();
        
    ClassDef(EventTrueCaptures, 1)
};

#endif
