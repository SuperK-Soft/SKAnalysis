#ifndef CalculateNeutronCloudVertex_H
#define CalculateNeutronCloudVertex_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "NeutronInfo.h"

class CalculateNeutronCloudVertex: public Tool {

 public:

  CalculateNeutronCloudVertex();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  double GetWeighting(const NeutronInfo) const;
  
};

#endif
