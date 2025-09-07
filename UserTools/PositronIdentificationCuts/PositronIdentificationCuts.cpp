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
  if (!ok || outfile_name.empty()){ outfile_name = "decay_electrons_cut_out.root";}
  
  outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (outfile == nullptr){
    throw std::runtime_error("PositronIdentificationCuts::Initialise - Couldn't open output file");
  }
  
  out_tree = new TTree("muechk","muechk");
  out_tree->Branch("rejected",&rejected);
  out_tree->Branch("nevsk",&skhead_.nevsk); // to check event alignment
  out_tree->Branch("nmue",&nmue);
  out_tree->Branch("mue_times",&mue_times);
  out_tree->Branch("q50n50",&q50n50_ratio);
  out_tree->Branch("maxpre",&max_pre);
  out_tree->Branch("maxpregate",&max_pregate);
  
  std::string plotfile_name = "";
  ok = m_variables.Get("plotfile_name", plotfile_name);
  if (!ok || plotfile_name.empty()){ plotfile_name = "decay_electrons_cut_plots.root";}
  
  plotfile = TFile::Open(plotfile_name.c_str(), "RECREATE");
  if (plotfile == nullptr){
    throw std::runtime_error("PositronIdentificationCuts::Initialise - Couldn't open plot file");
  }
  
  pre_q50n50_ratio_cut = TH1D("pre_q50n50_ratio_cut", "pre_q50n50_ratio_cut;q50/n50", 100, 0, 0);
  pre_nmue_cut = TH1D("pre_nmue_cut", "pre_nmue_cut;nmue", 10, 0, 0);
  pre_maxpre_cut = TH1D("pre_maxpre_cut", "pre_maxpre_cut;maxpre", 20, 0, 0);
  pre_maxpregate_cut = TH1D("pre_maxpregate_cut", "pre_maxpregate_cut;maxpregate", 20, 0, 0);
  
  //old_directory->cd();
  
  return true;
}


bool PositronIdentificationCuts::Execute(){
  
  rejected=true;
  
  m_data->CStore.Get("q50n50_ratio", q50n50_ratio);
  pre_q50n50_ratio_cut.Fill(q50n50_ratio);
  
  m_data->CStore.Get("nmue", nmue);
  m_data->CStore.Get("mue_times", mue_times);
  pre_nmue_cut.Fill(nmue);

  m_data->CStore.Get("max_pre", max_pre);
  pre_maxpre_cut.Fill(max_pre);
  
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
  
  rejected=false;
  ++accepted;
  out_tree->Fill();
  
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
    out_tree->Write();
    outfile->Close();
  }
  
  if(plotfile){
    plotfile->cd();
    pre_nmue_cut.Write();
    pre_maxpre_cut.Write();
    pre_maxpregate_cut.Write();
    pre_q50n50_ratio_cut.Write();
    plotfile->Close();
  }
  
  return true;
}

void PositronIdentificationCuts::SkipEntry(){
  bool skip = true;
  m_data->CStore.Set("Skip", skip);
  out_tree->Fill();
  return;
}
