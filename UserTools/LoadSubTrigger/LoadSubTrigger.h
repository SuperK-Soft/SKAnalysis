#ifndef LoadSubTrigger_H
#define LoadSubTrigger_H

#include "TH1D.h"

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

  std::string trigger_time_names = "";
  int trigger_idx = 0;
  
  int GetReaderLUN();
  int TreeReaderLUN;
  int neglun;
  
};


#endif
