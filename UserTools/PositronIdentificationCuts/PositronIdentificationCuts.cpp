#include "PositronIdentificationCuts.h"

#include "TFile.h"
#include "TH1D.h"

PositronIdentificationCuts::PositronIdentificationCuts():Tool(){}


bool PositronIdentificationCuts::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  pre_q50n50_ratio_cut = TH1D("pre_q50n50_ratio_cut", "pre_q50n50_ratio_cut;q50/n50", 100, 0, 0);
  pre_nmue_cut = TH1D("pre_nmue_cut", "pre_nmue_cut;nmue", 10, 0, 0);
  pre_maxpre_cut = TH1D("pre_maxpre_cut", "pre_maxpre_cut;maxpre", 20, 0, 0);
  pre_maxpregate_cut = TH1D("pre_maxpregate_cut", "pre_maxpregate_cut;maxpregate", 20, 0, 0);
  
  return true;
}


bool PositronIdentificationCuts::Execute(){

  double q50n50_ratio = 0;
  m_data->CStore.Get("q50n50_ratio", q50n50_ratio);
  pre_q50n50_ratio_cut.Fill(q50n50_ratio);
  
  int nmue = 0;
  m_data->CStore.Get("nmue", nmue);
  pre_nmue_cut.Fill(nmue);

  int max_pre = 0;
  m_data->CStore.Get("max_pre", max_pre);
  pre_maxpre_cut.Fill(max_pre);
  
  int max_pregate = 0;
  m_data->CStore.Get("max_pregate", max_pregate);
  pre_maxpregate_cut.Fill(max_pregate);
  
  if (nmue != 0){
    SkipEntry();
    return true;
  }

  if (max_pre >= 12){
    SkipEntry();
    return true;
  }
  
  if (max_pregate >= 5){
    SkipEntry();
    return true;
  }
 

  if (q50n50_ratio > 2){
    SkipEntry();
    return true;
  }
  
  return true;
}


bool PositronIdentificationCuts::Finalise(){

  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "decay_electrons_cuts_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "UPDATE");
  if (outfile == nullptr){
    throw std::runtime_error("PositronIdentificationCuts::Finalise - Couldn't open output file");
  }

  pre_nmue_cut.Write();
  pre_maxpre_cut.Write();
  pre_maxpregate_cut.Write();
  pre_q50n50_ratio_cut.Write();
  
  return true;
}

void PositronIdentificationCuts::SkipEntry(){
  bool skip = true;
  m_data->CStore.Set("Skip", skip);
  return;
}
