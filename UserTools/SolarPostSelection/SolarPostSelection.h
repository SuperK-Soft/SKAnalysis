#ifndef SolarPostSelection_H
#define SolarPostSelection_H

#include <string>
#include <array>
#include <iostream>

#include "Tool.h"
#include "TH1D.h"
//#include "SolarRelic.h"

class SolarPostSelection: public Tool {

public:

  SolarPostSelection();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  TFile* plotfile = nullptr;
  TH1D bsok_cut;
  TH1D pre_bsenergy_cut;
  TH1D pre_ovaQ_cut;
  TH1D pre_FV_cut;
  TH1D pre_effwall_cut;
  
  TH1D post_bsenergy_cut;
  TH1D post_ovaQ_cut;
  TH1D post_FV_cut;
  TH1D post_effwall_cut;

  
  
};

#endif
