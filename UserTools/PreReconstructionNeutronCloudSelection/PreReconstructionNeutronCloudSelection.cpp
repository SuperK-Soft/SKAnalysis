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

  pre_dt_cut = TH1D("pre_dt_cut", "pre_dt_cut;dt", 50, 0, 0);
  post_dt_cut = TH1D("post_dt_cut", "post_dt_cut;dt", 50, 0, 0);
  
  return true;
}


bool PreReconstructionNeutronCloudSelection::Execute(){

  /* 
     We can cut on time difference between muon and SLE before we go through reconstruction:
     Cut all SLE events where dt is < 35us or > 535us;
  */
  Log("PreReconstructionNeutronCloudSelection: cut out all SLE triggers that have dt < 35us || dt > 535us", v_debug, m_verbose);
  
  std::vector<double> SLE_times;
  m_data->CStore.Get("SLE_times", SLE_times); // in nsec remember

  Log(m_unique_name+": SLE times pre cut:",v_debug,m_verbose);
  for (const auto& dt : SLE_times){
    Log(toString(dt)+"ns or , "+toString(dt / 1000)+"us",v_debug,m_verbose);
    pre_dt_cut.Fill(dt);
  }

  std::vector<double> post_cut_SLE_times = {};
  std::copy_if(SLE_times.begin(), SLE_times.end(), std::back_inserter(post_cut_SLE_times),
	       [](double time_ns){
		 return (time_ns >= 35E3 && time_ns <= 535E3);});

  for (const auto& dt: post_cut_SLE_times){
    post_dt_cut.Fill(dt);
    Log(m_unique_name+" trigger making the cut at "+toString(dt),v_debug,m_verbose);
  }
  
  Log(m_unique_name+": number of SLE triggers making the dt cut: "+toString(post_cut_SLE_times.size()),v_debug,m_verbose);
  
  m_data->CStore.Set("SLE_times", post_cut_SLE_times);
  const int N_SLE = post_cut_SLE_times.size();
  m_data->CStore.Set("N_SLE", N_SLE); //need this for the subtoolchain                                                                                      

  
  return true;
}


bool PreReconstructionNeutronCloudSelection::Finalise(){

  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "pre_recon_neut_cloud_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()){
    throw std::runtime_error("PreReconstructionNeutronCloudSelection::Finalise - Couldn't open output file");
  }

  pre_dt_cut.Write();
  post_dt_cut.Write();

  return true;
}
