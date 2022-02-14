#ifndef THIRDREDVARS_H
#define THIRDREDVARS_H

#include "TObject.h"
#include "TNamed.h"
#include <vector>

// in this file we define the classes that
// contain the values of the observables used for third reduction


class ThirdRed : public TNamed {
public:
    
    int odn51; // OD hits within 5m of lowe hit (use OD window)
    int odn52; // OD hits within 5m of lowe hit (use tof sub window)
    int odclus1; // OD clusters (use OD window)
    int odclus2; // OD clusters (use tof sub window)
    int maxpre, maxpost; // pre/post activity (all hits)
    int maxpregate, maxpostgate; // pre/post activity (in-gate hits)
    float ring_angle; // angle between two rings (if two rings)
    float ari_msg; // multiple scattering goodness
    float angle; // cherenkov angle
    float pilike; // fuzzyness of ring (pion identification)
    float q50; // charge in 50ns window with maximal number of hits
    int n50; // maximal number of hits in 50ns window
    float ovaq; // goodness**2 - dirks**2
    int nmue; // reconstructed number of decay electrons
    float dwall; // distance to wall
    float effwall; // effective distance to wall

    void Clear() {
        odn51 = 0;
        odn52 = 0;
        odclus1 = 0;
        odclus2 = 0;
        maxpre = 0;
        maxpost = 0;
        ring_angle = 0;
        ari_msg = -5;
        angle = -5;
        pilike = -5;
        q50 = 0;
        n50 = -1;
        ovaq = 0;
        nmue = 0;
        dwall = 1000000;
        effwall = 10000000;
    }
  
    ClassDef(ThirdRed,3) // increase version number when structure changed. 
};

#endif
