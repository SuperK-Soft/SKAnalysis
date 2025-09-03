#include "PositronIdentificationCuts.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"

PositronIdentificationCuts::PositronIdentificationCuts():Tool(){}


bool PositronIdentificationCuts::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  
  nmue_thresh=0;
  max_pre_threshold=12;
  max_pregate_threshold=5;
  q50n50_threshold=2;
  m_variables.Get("nmue_thresh",nmue_thresh);
  m_variables.Get("max_pre_threshold",max_pre_threshold);
  m_variables.Get("max_pregate_threshold",max_pregate_threshold);
  m_variables.Get("q50n50_threshold",q50n50_threshold);
  
  //TDirectory* old_directory = gDirectory;
  TDirectory::TContext ctxt(); // restores gDirectory at end of scope
  // constructor may also take old and new directories, effectively cd's when appropriate
  // note this is better as it handles exceptions and early returns as well.
  
  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "decay_electrons_cuts_out.root";}

  outfile = TFile::Open(outfile_name.c_str(), "UPDATE");
  if (outfile == nullptr){
    throw std::runtime_error("PositronIdentificationCuts::Initialise - Couldn't open output file");
  }

  pre_q50n50_ratio_cut = TH1D("pre_q50n50_ratio_cut", "pre_q50n50_ratio_cut;q50/n50", 100, 0, 0);
  pre_nmue_cut = TH1D("pre_nmue_cut", "pre_nmue_cut;nmue", 10, 0, 0);
  pre_maxpre_cut = TH1D("pre_maxpre_cut", "pre_maxpre_cut;maxpre", 20, 0, 0);
  pre_maxpregate_cut = TH1D("pre_maxpregate_cut", "pre_maxpregate_cut;maxpregate", 20, 0, 0);
  
  //old_directory->cd();
  
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
  
  ++total;
  
  if (nmue != nmue_thresh){
    SkipEntry();
    ++nmue_rejects;
    return true;
  }

  if (max_pre >= max_pre_threshold){
    SkipEntry();
    ++max_pre_rejects;
    return true;
  }
  
  if (max_pregate >= max_pregate_threshold){
    SkipEntry();
    ++max_pre_gate_rejects;
    return true;
  }

  if (q50n50_ratio > q50n50_threshold){
    SkipEntry();
    ++q50n50_rejects;
    return true;
  }
  
  ++accepted;
  
  return true;
}


bool PositronIdentificationCuts::Finalise(){

  Log(m_unique_name+" Accepted "+toString(accepted)+"/"+toString(total)+" events\n"+
      toString(nmue_rejects)+" events rejected by nmue==0 cut\n"+
      toString(max_pre_rejects)+" events rejected by max_pre cut\n"+
      toString(max_pre_gate_rejects)+" events rejected by max_pre_gate cut\n"+
      toString(q50n50_rejects)+" events rejected by q50/n50 cut",v_debug,m_verbose);
  
  if(outfile){
    outfile->cd();
    pre_nmue_cut.Write();
    pre_maxpre_cut.Write();
    pre_maxpregate_cut.Write();
    pre_q50n50_ratio_cut.Write();
    outfile->Close();
  }
  
  return true;
}

void PositronIdentificationCuts::SkipEntry(){
  bool skip = true;
  m_data->CStore.Set("Skip", skip);
  return;
}
