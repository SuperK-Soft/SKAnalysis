#ifndef PostReconstructionNeutronCloudSelection_H
#define PostReconstructionNeutronCloudSelection_H

#include <string>
#include <array>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

#include "NeutronInfo.h"

class PostReconstructionNeutronCloudSelection: public Tool {

public:

  PostReconstructionNeutronCloudSelection();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  TH1D pre_bsgood_cut;
  TH1D pre_bsdirks_cut;
  TH1D pre_bsn50_cut;
  TH1D pre_ldt_cut;

  TH1D post_bsgood_cut;
  TH1D post_bsdirks_cut;
  TH1D post_bsn50_cut;
  TH1D post_ldt_cut;

  MTreeReader* tree_reader_ptr = nullptr;
  void GetReader();
  
  std::vector<NeutronInfo> neutrons = {}; //make struct
  
  float CalculateLongitudinalDistance(float*, float*, float*) const;
};

#endif
