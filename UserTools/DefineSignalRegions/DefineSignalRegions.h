#ifndef DefineSignalRegions_H
#define DefineSignalRegions_H

#include <string>
#include <iostream>

#include "Tool.h"

class DefineSignalRegions: public Tool {

public:

  DefineSignalRegions(); 
  bool Initialise(std::string configfile,DataModel &data); 
  bool Execute(); 
  bool Finalise(); 

private:

};

#endif
