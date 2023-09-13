#ifndef NeutCloudCorrelationCuts_H
#define NeutCloudCorrelationCuts_H

#include <string>
#include <iostream>

#include "Tool.h"

class NeutCloudCorrelationCuts: public Tool {

 public:

  NeutCloudCorrelationCuts();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  std::string relic_reader_name = "";
  std::vector<TVector3> GetTensor(const std::vector<double>&, const std::vector<double>&) const;
  
};

#endif
