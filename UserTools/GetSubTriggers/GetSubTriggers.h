#ifndef GetSubTriggers_H
#define GetSubTriggers_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class GetSubTriggers
*
* Please fill out the descripton and author information.
*
* $Author: ?.????? $
* $Date: ????/??/?? $
* $Contact: ???@km.icrr.u-tokyo.ac.jp
*/

class GetSubTriggers: public Tool {

 public:

  GetSubTriggers();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

};


#endif
