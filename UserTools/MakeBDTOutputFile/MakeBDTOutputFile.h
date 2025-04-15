#ifndef MakeBDTOutputFile_H
#define MakeBDTOutputFile_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TFile.h"
#include "TTree.h"

#include "MTreeReader.h"

class MakeBDTOutputFile: public Tool {

 public:

  MakeBDTOutputFile();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  bool ok = false;
  
  int entry_number_tmp = 0;
  
  TFile* output_file_ptr = nullptr;
  TTree* output_tree_ptr = nullptr;
  
  MTreeReader* reader_input_ptr;
    MTreeReader* reader_ntagbdt_ptr;

  Header* header_ptr = nullptr;
  LoweInfo* lowe_ptr = nullptr;
  MCInfo* mcinfo_ptr = nullptr;
  TQReal* tqreal_ptr = nullptr;

  float current_first_neutron_likelihood = 0;
  
  std::vector<float> neutron_likelihoods = {};
  
  void SetOutputBranches(TTree* tree);
  
};


#endif
