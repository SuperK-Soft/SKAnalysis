#include "DecayElectronCuts.h"

#include "TFile.h"
#include "TH1D.h"

DecayElectronCuts::DecayElectronCuts():Tool(){}


bool DecayElectronCuts::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  pre_nmue_cut = TH1D("pre_nmue_cut", "pre_nmue_cut;nmue", 10, 0, 10);
  pre_maxpre_cut = TH1D("pre_maxpre_cut", "pre_maxpre_cut;maxpre", 20, 0, 20);
  pre_maxpregate_cut = TH1D("pre_maxpregate_cut", "pre_maxpregate_cut;maxpregate", 20, 0, 20);
  
  return true;
}


bool DecayElectronCuts::Execute(){

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
  
  if(max_pregate >= 5){
    SkipEntry();
    return true;
  }
  
  

  
  return true;
}


bool DecayElectronCuts::Finalise(){

  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "decay_electrons_cuts_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "UPDATE");
  if (outfile == nullptr){
    throw std::runtime_error("DecayElectronCuts::Finalise - Couldn't open output file");
  }

  pre_nmue_cut.Write();
  pre_maxpre_cut.Write();
  pre_maxpregate_cut.Write();
  
  return true;
}

void DecayElectronCuts::SkipEntry(){
  bool skip = true;
  m_data->CStore.Set("Skip", skip);
  return;
}
