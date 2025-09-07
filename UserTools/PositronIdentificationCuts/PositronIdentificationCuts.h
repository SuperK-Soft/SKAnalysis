#ifndef PositronIdentificationCuts_H
#define PositronIdentificationCuts_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

class PositronIdentificationCuts: public Tool {

 public:

  PositronIdentificationCuts();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
  int total=0;
  int nmue_rejects=0;
  int max_pre_rejects=0;
  int max_pre_gate_rejects=0;
  int q50n50_rejects=0;
  int accepted=0;
  
  bool rejected=false; // current event
  
  double q50n50_ratio = 0;
  int nmue = 0;
  std::vector<double> mue_times;
  int max_pre = 0;
  int max_pregate = 0;
  
  int q50n50_threshold;
  int nmue_thresh;
  int max_pre_threshold;
  int max_pregate_threshold;
 
  TFile* outfile=nullptr; 
  TTree* out_tree=nullptr;
  
  TFile* plotfile=nullptr; 
  TH1D pre_q50n50_ratio_cut;
  TH1D pre_nmue_cut;
  TH1D pre_maxpre_cut;
  TH1D pre_maxpregate_cut;

  void SkipEntry();
  
};


#endif
