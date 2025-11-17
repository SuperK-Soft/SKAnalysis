#ifndef WriteSolarMatches_H
#define WriteSolarMatches_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SolarRelic.h"
#include "MTreeReader.h"


class WriteSolarMatches: public Tool {

 public:

  WriteSolarMatches();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
 
  std::deque<SolarRelic>* relics_this_run=nullptr;
  
  void GetReader();
 
  MTreeReader* rfmReader = nullptr;
  TFile* outfile = nullptr;
  TFile* solarFile = nullptr;
  TTree* solarTree = nullptr;
 
  std::vector<int> matched_relics;
  std::vector<double> matched_tdiffs;
  std::vector<double> matched_dists;
  std::set<int> relics_to_write;
 
  std::string solarSelectorName;
  int outentry;
  double match_dist_limit=400.; // [cm] (SK-VI tech note, SK-IV used 490)

};


#endif
