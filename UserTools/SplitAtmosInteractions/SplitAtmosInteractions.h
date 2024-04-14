#ifndef SplitAtmosInteractions_H
#define SplitAtmosInteractions_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "MTreeReader.h"

#include "TH1D.h"

// const static std::map<int, std::string> interaction_modes = {
  
// };

class SplitAtmosInteractions: public Tool {

 public:

  SplitAtmosInteractions();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  MTreeReader* tree_reader_ptr = nullptr;
  void GetReader();

  TH1D interaction_mode_plot;
  TH1D interaction_mode_zero;
  TH1D NCQE_plot;
  TH1D non_NCQE_plot;

  std::map<std::string, int> interactions_for_pi_chart;

};

#endif


