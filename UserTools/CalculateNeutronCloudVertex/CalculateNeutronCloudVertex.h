#ifndef CalculateNeutronCloudVertex_H
#define CalculateNeutronCloudVertex_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TH1D.h"

class CalculateNeutronCloudVertex: public Tool {

 public:

  CalculateNeutronCloudVertex();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  TH1D mult_plot;
  
  int mult = 0;
  MTreeReader* LOWE_tree_reader = nullptr;

  void GetTreeReader();
  
};

#endif
