#ifndef CalculatePreactivityObservables_H
#define CalculatePreactivityObservables_H

#include <string>
#include <iostream>

#include "Tool.h"

class CalculatePreactivityObservables: public Tool {

 public:

  CalculatePreactivityObservables();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

};


#endif
