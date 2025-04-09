#ifndef SLESearch_H
#define SLESearch_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TH1D.h"

class SLESearch: public Tool {

 public:

  SLESearch();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  bool previous_entry_was_muon = false;

  int max_triggers = -999;
  
  ConnectionTable* connection_table = nullptr;
  bool include_offset = false;
  const double SLE_t0_offset = -885.417;
  
  double TimeOfFlight(const float*, const float*) const;

  TH1D hit_times_plot;
  TH1D event_hit_times_plot;
  
};


#endif
