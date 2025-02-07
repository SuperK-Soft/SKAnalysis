#ifndef GetSubTrg_SetTimingGate_H
#define GetSubTrg_SetTimingGate_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class GetSubTrg_SetTimingGate
*
* Please fill out the descripton and author information.
*
* $Author: ?.????? $
* $Date: ????/??/?? $
* $Contact: ???@km.icrr.u-tokyo.ac.jp
*/

class GetSubTrg_SetTimingGate: public Tool {

 public:

  GetSubTrg_SetTimingGate();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

};


#endif
