#ifndef CalculatePreactivityObservables_H
#define CalculatePreactivityObservables_H

#include <string>
#include <iostream>

#include "Tool.h"

struct Hit {
  double time = 0;
  double goodness = 0;
  Hit(const double& t, const double& g) : time{t}, goodness{g} {}
};

class CalculatePreactivityObservables: public Tool {

 public:

  CalculatePreactivityObservables();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  double dark_threshold = 4;
  double fraction = 0.4;
  ConnectionTable* connection_table = nullptr;

  double TimeOfFlight(const float*, const float*) const;
  bool CalculateGoodness(const double&, const double&) const;
  
};


#endif
