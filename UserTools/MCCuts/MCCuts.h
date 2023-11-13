#ifndef MCCuts_H
#define MCCuts_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

class MCCuts: public Tool {

 public:

  MCCuts();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  TH1D nqisk_plot;
  TH1D bsenergy_plot;
  TH1D bsgoodness_plot;
  TH1D ovaq_plot;
  TH1D clusfit_goodness_plot;
  TH1D d_wall_plot;
  TH1D n_od_plot;
  
};

#endif
