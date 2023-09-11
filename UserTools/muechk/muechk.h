#ifndef muechk_H
#define muechk_H

#include <string>
#include <iostream>

#include "Tool.h"

class muechk: public Tool {

 public:

  muechk();
  bool Initialise(std::string configfile,DataModel &data); 
  bool Execute();
  bool Finalise();

 private:

  MTreeReader* tree_reader_ptr = nullptr;
  
};

#endif
