#ifndef MakeSpectralFitHistos_H
#define MakeSpectralFitHistos_H

#include <string>
#include <iostream>

#include "MTreeReader.h"

#include "Tool.h"

#include "TH1D.h"

class MakeSpectralFitHistos: public Tool {

 public:
  
  MakeSpectralFitHistos();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  MTreeReader* rw_tree_reader_ptr = nullptr;
  MTreeReader* bdt_tree_reader_ptr = nullptr;
  std::string tree_reader_str = "";
  std::string weight_name = "";
  
  TH1D mc_energy_full;
  TH1D bs_energy_full;
  TH1D bs_energy_bdt_full;
  TH1D bs_energy_bdt_mccuts_full;

  TH1D r0;
  TH1D r1;
  TH1D r2;
  TH1D r3;
  TH1D r4;
  TH1D r5;

  void GetReaders();
  void GetWeightName();
  bool HasExactlyOneNeutron(const basic_array<float*> likelihoods);

  
};


#endif
