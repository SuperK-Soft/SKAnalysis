#ifndef DefineSignalRegions_H
#define DefineSignalRegions_H

#include <string>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"

#include "MTreeReader.h"

#include "Tool.h"

class DefineSignalRegions: public Tool {

public:

  DefineSignalRegions(); 
  bool Initialise(std::string configfile,DataModel &data); 
  bool Execute(); 
  bool Finalise(); 

private:

  MTreeReader* LOWE_TreeReader = nullptr;
  MTreeReader* BDT_TreeReader = nullptr;

  bool need_to_get_treereaders;
  int events_matched = 0;
  int current_lowe_event_number = 0;
  int current_bdt_event_number = 0;

  TH1D r0;
  TH1D r1;
  TH1D r2;
  TH1D r3;
  TH1D r4;
  TH1D r5;
    
  void GetTreeReaders();
  //bool HasExactlyOneNeutron(const std::vector<double>&);
  bool HasExactlyOneNeutron(const basic_array<float*>);
  void MakePlot();
};

#endif
