#ifndef SplitAtmosInteractions_H
#define SplitAtmosInteractions_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "MTreeReader.h"

class SplitAtmosInteractions: public Tool {

 public:

  SplitAtmosInteractions();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  MTreeReader* tree_reader_ptr = nullptr;
  void GetReader();
};


#endif
