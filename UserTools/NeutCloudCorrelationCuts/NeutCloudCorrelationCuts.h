#ifndef NeutCloudCorrelationCuts_H
#define NeutCloudCorrelationCuts_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"

#include "TH1D.h"

class NeutCloudCorrelationCuts: public Tool {

 public:

  NeutCloudCorrelationCuts();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  std::string relic_reader_name = "";
  std::vector<TVector3> GetTensor(const std::vector<double>&, const std::vector<double>&) const;
  void SkipEntry();
  void GetTreeReaders();

  TH1D pre_sample_total_dt;
  TH1D pre_sample_m2_dt;
  TH1D pre_sample_m3_dt;
  TH1D pre_sample_m45_dt;
  TH1D pre_sample_m69_dt;
  TH1D pre_sample_m10_dt;

  TH1D post_sample_total_dt;
  TH1D post_sample_m2_dt;
  TH1D post_sample_m3_dt;
  TH1D post_sample_m45_dt;
  TH1D post_sample_m69_dt;
  TH1D post_sample_m10_dt;

  TH1D pre_sample_total_dl;
  TH1D pre_sample_m2_dl;
  TH1D pre_sample_m3_dl;
  TH1D pre_sample_m45_dl;
  TH1D pre_sample_m69_dl;
  TH1D pre_sample_m10_dl;

  TH1D post_sample_total_dl;
  TH1D post_sample_m2_dl;
  TH1D post_sample_m3_dl;
  TH1D post_sample_m45_dl;
  TH1D post_sample_m69_dl;
  TH1D post_sample_m10_dl;

  MTreeReader* relic_tree_reader = nullptr;
  MTreeReader* cloud_tree_reader = nullptr;
  std::string cloud_tree_reader_str = "";
  
  std::vector<int>* relicMatchedEntryNums = nullptr;
  std::vector<float>* relicTimeDiffs = nullptr;
 
};

#endif
