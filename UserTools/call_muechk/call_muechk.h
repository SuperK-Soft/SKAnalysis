#ifndef call_muechk_H
#define call_muechk_H

#include <string>
#include <iostream>

#include "Tool.h"

class call_muechk: public Tool {

 public:

  call_muechk();
  bool Initialise(std::string configfile,DataModel &data); 
  bool Execute();
  bool Finalise();

 private:

  MTreeReader* tree_reader_ptr = nullptr;
  
};

#endif
