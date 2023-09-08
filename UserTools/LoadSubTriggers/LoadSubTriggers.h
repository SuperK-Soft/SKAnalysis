#ifndef LoadSubTriggers_H
#define LoadSubTriggers_H

#include <string>
#include <iostream>

#include "Tool.h"

class LoadSubTriggers: public Tool {

 public:

  LoadSubTriggers();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

};


#endif
