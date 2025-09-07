#ifndef CalculatePreactivityObservables_H
#define CalculatePreactivityObservables_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "TH1.h"

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
  double q50n50_window_size = 50; //think this is in ns
  double preact_window_size = 15;
  double preact_window_cutoff = 12;
  ConnectionTable* connection_table = nullptr;
  MTreeReader* LOWE_tree_reader;
  TH1F h_maxpre;
  TH1F h_maxpregate;
  TH1F h_bsvertex_t;
  TH1F h_q50n50;

  void GetTreeReader();
  double TimeOfFlight(const float*, const float*) const;
  double CalculateGoodness(const double&, const double&) const;
  
};


#endif
