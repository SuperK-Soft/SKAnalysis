#include "PreReconstructionNeutronCloudSelection.h"

#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"

PreReconstructionNeutronCloudSelection::PreReconstructionNeutronCloudSelection():Tool(){}


bool PreReconstructionNeutronCloudSelection::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  pre_dt_cut = TH1D("pre_dt_cut", "pre_dt_cut;dt", 100, 0, 1000);
  post_dt_cut = TH1D("post_dt_cut", "post_dt_cut;dt", 100, 0, 1000);
  
  return true;
}


bool PreReconstructionNeutronCloudSelection::Execute(){

  /* 
     We can cut on time difference between muon and SLE before we go through reconstruction:
     Cut all SLE events where dt is < 35us or > 535us;
  */

  std::vector<int> SLE_times;
  m_data->CStore.Get("SLE_times", SLE_times);

  for (const auto& dt : SLE_times){
    pre_dt_cut.Fill(dt);
  }

  std::vector<int> post_cut_SLE_times = {};
  std::copy_if(SLE_times.begin(), SLE_times.end(), std::back_inserter(post_cut_SLE_times),
	       [](int i){
		 double ticks = i / (COUNT_PER_NSEC * 1000);
		 return (ticks >= 35 && ticks <= 535);});

  for (const auto& dt: post_cut_SLE_times){
    post_dt_cut.Fill(dt);
  }

  m_data->CStore.Set("SLE_times", SLE_times);
  
  return true;
}


bool PreReconstructionNeutronCloudSelection::Finalise(){

  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "pre_recon_neut_cloud_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "UPDATE");
  if (outfile == nullptr){
    throw std::runtime_error("PreReconstructionNeutronCloudSelection::Finalise - Couldn't open output file");
  }

  pre_dt_cut.Write();
  post_dt_cut.Write();

  return true;
}
