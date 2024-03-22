#ifndef SpectralFit_H
#define SpectralFit_H

#include <string>
#include <vector>
#include <iostream>
#include <map>

#include "TGraph.h"

#include "Tool.h"

static const int N_regions = 6;
static const int N_distributions = 6;

struct Distribution {
  std::array<TGraph, N_regions> regions;
};

class SpectralFit: public Tool {

public:

  SpectralFit();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

private:

  void GetPDFs();

  const std::array<std::string, N_distributions> distribution_names = {
    "SRN signal"
    "Invisible muons and pions",
    "nu_e CC interactions",
    "mu/pi-producing interactions",
    "NCQE interactions",
    "Spallation backgrounds"};

  const std::array<std::string, N_regions> region_names = {
    "20-38deg, N_tagged = 1",
    "38-50deg, N_tagged = 1",
    "78-90deg, N_tagged = 1",
    "20-38deg, N_tagged != 1",
    "38-50deg, N_tagged != 1",
    "78-90deg, N_tagged != 1"};
  
  std::map<std::string, std::array<TGraph, N_regions>> distributions;
  std::array<TGraph, N_distributions> global_distributions;
  
};

#endif
