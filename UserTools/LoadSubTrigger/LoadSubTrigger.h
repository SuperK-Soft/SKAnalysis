#ifndef LoadSubTrigger_H
#define LoadSubTrigger_H

#include <string>
#include <iostream>

#include "Tool.h"

class LoadSubTrigger: public Tool {

 public:

  LoadSubTrigger();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  int trigger_idx = 0;
  int GetReaderLUN();
  
};


#endif
