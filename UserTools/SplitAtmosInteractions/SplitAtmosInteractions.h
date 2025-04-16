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

  MTreeReader* tree_reader_root_ptr = nullptr;
  MTreeReader* tree_reader_ZBS_ptr = nullptr;

  TFile* output_file_ptr =  nullptr;
  TTree* output_tree_ptr = nullptr;

  LoweInfo* lowe_ptr = nullptr;
  MCInfo* mcinfo_ptr = nullptr;
  Header* header_ptr = nullptr;
  TQReal* tqreal_ptr = nullptr;
  TQReal* tqareal_ptr = nullptr;
  SecondaryInfo* secondary_info_ptr = nullptr;

  std::vector<float>*  neutron_likelihoods_ptr;
  int interaction_mode = 0;
  
  void GetReaders();
  void GetValues();
  void UpdateBranches();
  void SetOutputBranches(TTree*);
  
  bool first = true;
  
  // TH1D interaction_mode_plot;
  // TH1D interaction_mode_zero;
  // TH1D NCQE_plot;
  // TH1D non_NCQE_plot;

  // std::map<std::string, int> interactions_for_pi_chart;

};

#endif


