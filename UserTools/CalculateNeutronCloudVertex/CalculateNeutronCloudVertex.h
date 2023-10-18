#ifndef CalculateNeutronCloudVertex_H
#define CalculateNeutronCloudVertex_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

class CalculateNeutronCloudVertex: public Tool {

 public:

  CalculateNeutronCloudVertex();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();
  
 private:

  TH1D mult_plot;
  TH1D dist_to_mu_plot;
  
  int mult = 0;
  std::vector<double> neutron_cloud_vertex = {};
  MTreeReader* MU_tree_reader = nullptr;

  TFile* nvc_file_ptr = nullptr;
  TTree* nvc_tree_ptr = nullptr;

  void GetTreeReader();
  void CreateOutputFile();
  double ClosestApproach(const std::vector<double>&) const;
  
};

#endif
