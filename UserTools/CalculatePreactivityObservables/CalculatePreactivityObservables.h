#ifndef CalculatePreactivityObservables_H
#define CalculatePreactivityObservables_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"

struct Hit {
  double time = 0;
  double goodness = 0;
  double charge = 0;
  Hit(const double& t, const double& g, const double& q) : time{t}, goodness{g}, charge{q} {}
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
  MTreeReader* LOWE_tree_reader;

  void GetTreeReader();
  double TimeOfFlight(const float*, const float*) const;
  bool CalculateGoodness(const double&, const double&) const;
  
};


#endif
