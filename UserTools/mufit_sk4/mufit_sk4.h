#ifndef mufit_sk4_H
#define mufit_sk4_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

class mufit_sk4: public Tool {

public:

  mufit_sk4(); 
  bool Initialise(std::string configfile,DataModel &data); 
  bool Execute(); 
  bool Finalise(); 

private:

  MTreeReader* tree_reader_ptr = nullptr;
  void GetReader();

  std::string tree_reader_str = "";
  
  int nread = 0;
  int nmuon = 0;

  TH1D charge_plot;
  
};

#endif
