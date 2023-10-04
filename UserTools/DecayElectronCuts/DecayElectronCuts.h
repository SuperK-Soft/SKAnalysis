#ifndef DecayElectronCuts_H
#define DecayElectronCuts_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

class DecayElectronCuts: public Tool {

 public:

  DecayElectronCuts();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
  
  TH1D pre_q50n50_ratio_cut;
  TH1D pre_nmue_cut;
  TH1D pre_maxpre_cut;
  TH1D pre_maxpregate_cut;

  void SkipEntry();
  
};


#endif
