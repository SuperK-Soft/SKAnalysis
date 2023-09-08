#ifndef SLESearch_H
#define SLESearch_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class SLESearch
*
* Please fill out the descripton and author information.
*
* $Author: ?.????? $
* $Date: ????/??/?? $
* $Contact: ???@km.icrr.u-tokyo.ac.jp
*/

class SLESearch: public Tool {

 public:

  SLESearch();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  bool previous_entry_was_muon = false;
  
};


#endif
