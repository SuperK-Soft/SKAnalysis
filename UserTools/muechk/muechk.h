#ifndef muechk_H
#define muechk_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

class muechk: public Tool {

 public:

  muechk();
  bool Initialise(std::string configfile,DataModel &data); 
  bool Execute();
  bool Finalise();

 private:

  MTreeReader* tree_reader_ptr = nullptr;

  TH1D nmue_plot;
  
};

#endif
