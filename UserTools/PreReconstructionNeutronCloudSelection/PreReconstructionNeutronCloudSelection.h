#ifndef PreReconstructionNeutronCloudSelection_H
#define PreReconstructionNeutronCloudSelection_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

class PreReconstructionNeutronCloudSelection: public Tool {

 public:

  PreReconstructionNeutronCloudSelection();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
  
  TFile* outfile=nullptr;
  TH1D pre_dt_cut;
  TH1D post_dt_cut;
  
  
  
};


#endif
