#ifndef MakeBDTOutputFile_H
#define MakeBDTOutputFile_H

#include <string>
#include <iostream>
#include <map>

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
  bool carry_ibd_weights = false;
  int entry_number_tmp = 0;
  std::string outfile_str = "";
  
  TFile* output_file_ptr = nullptr;
  TTree* output_tree_ptr = nullptr;
  
  MTreeReader* reader_input_ptr;
  MTreeReader* reader_ntagbdt_ptr;

  Header* header_ptr = nullptr;
  LoweInfo* lowe_ptr = nullptr;
  MCInfo* mcinfo_ptr = nullptr;
  SecondaryInfo* secondary_info_ptr = nullptr;
  TQReal* tqreal_ptr = nullptr;
  TQReal* tqareal_ptr = nullptr;
  TClonesArray* tqlist_ptr = nullptr;
  TClonesArray* odtqlist_ptr = nullptr;
  float current_first_neutron_likelihood = 0;

  float* neutron5_temp = nullptr;
  int np_temp = 0;
  
  std::vector<float> neutron_likelihoods = {};

  float* x_temp;
  float* y_temp;
  float* z_temp;
  float* bse_temp;
  float* dtn_temp;
  
  void SetOutputBranches(TTree* tree);

  std::map<std::string, float*> srn_weights = {
    {"weight_tabriziNH", nullptr},
    {"weight_galaisNH", nullptr},
    {"weight_horiuchi3", nullptr},
    {"weight_horiuchi2", nullptr},
    {"weight_galaisIH", nullptr},
    {"weight_kresseNHhigh", nullptr},
    {"weight_nakazato2", nullptr},
    {"weight_malaney", nullptr},
    {"weight_nakazato1", nullptr},
    {"weight_totani", nullptr},
    {"weight_lunardini", nullptr},
    {"weight_hartmann", nullptr},
    {"weight_ando", nullptr},
    {"weight_li9", nullptr},
    {"weight_reactor", nullptr}};
  
};


#endif
