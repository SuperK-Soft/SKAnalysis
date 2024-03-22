#include "MakeSpectralFitHistos.h"

#include "MTreeReader.h"
#include "TROOT.h"

#include "TFile.h"
#include <cstdlib>
MakeSpectralFitHistos::MakeSpectralFitHistos():Tool(){}

bool MakeSpectralFitHistos::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReaders();
  GetWeightName();

  gROOT->cd();

  mc_energy_full = TH1D("mc", "mc", 100, 0, 0);
  bs_energy_full = TH1D("bs", "bs", 100, 0, 0);
  bs_energy_bdt_full = TH1D("bs_bdt", "bs_bdt", 100, 0, 0);
  bs_energy_bdt_mccuts_full = TH1D("bs_bdt_mccuts", "bs_bdt_mccuts", 100, 0, 0);
  
  r0 = TH1D("r0", "20#circ < #theta_{c} < 38#circ : N_{tagged} != 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 0);
  r1 = TH1D("r1", "38#circ < #theta_{c} < 53#circ : N_{tagged} != 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 0);
  r2 = TH1D("r2", "70#circ < #theta_{c} < 90#circ : N_{tagged} != 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 0);
  r3 = TH1D("r3", "20#circ < #theta_{c} < 38#circ : N_{tagged} == 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 0);
  r4 = TH1D("r4", "38#circ < #theta_{c} < 53#circ : N_{tagged} == 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 0);
  r5 = TH1D("r5", "70#circ < #theta_{c} < 90#circ : N_{tagged} == 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 0);
  
  return true;
}


bool MakeSpectralFitHistos::Execute(){

  bool ok = false;
  
  float weight = 1;
  if (weight_name != "unweighted"){
    std::cout << "current entry: " << rw_tree_reader_ptr->GetEntryNumber() << std::endl;
    rw_tree_reader_ptr->GetTree()->Show();
    std::cout << "weight name: " << weight_name << std::endl;
    ok = rw_tree_reader_ptr->GetBranchValue(weight_name, weight);
    if (!ok){throw std::runtime_error("MakeSpectralFitHistos::Execute bad branch value");}
  }

  if (std::isnan(weight) || std::isinf(weight)){
    return true;
  }
  
  std::cout << "weight:  " << weight << std::endl;
  
  // if (weight > 400e6{
  //   return true;
  // }
     
  double mcenergy = 0;
  m_data->CStore.Get("mcenergy", mcenergy);
  std::cout << "mcenergy:  " << mcenergy << std::endl;
  mc_energy_full.Fill(mcenergy, weight);

  double bsenergy = 0;
  m_data->CStore.Get("bsenergy", bsenergy);
  std::cout << "bsenergy:  " << bsenergy << std::endl;

  if( bsenergy > 1000){
    return true;
  }
  
  bs_energy_full.Fill(bsenergy, weight);

  basic_array<float*> neutron5;
  ok = bdt_tree_reader_ptr->GetBranchValue("neutron5", neutron5);
  if (!ok){
    throw std::runtime_error("MakeSpectralFitHistos::Execute - Couldn't Get() variable neutron5");
  }

  const bool has_neutron = HasExactlyOneNeutron(neutron5);
  
  if (has_neutron){
    bs_energy_bdt_full.Fill(bsenergy, weight);
  }

  bool pass = false;
  m_data->CStore.Get("pass", pass);

  if (!pass){
    return true;
  }

  if (has_neutron){
    bs_energy_bdt_mccuts_full.Fill(bsenergy, weight);  
  }
  
  double cherenkov_angle = 0;
  m_data->CStore.Get("c_angle", cherenkov_angle);
  if (cherenkov_angle > 20 && cherenkov_angle < 38){
    has_neutron ? r3.Fill(bsenergy, weight) : r0.Fill(bsenergy, weight);
  }

  if (cherenkov_angle > 38 && cherenkov_angle < 53){
    has_neutron ? r4.Fill(bsenergy, weight) : r1.Fill(bsenergy, weight);
  }
    
  if (cherenkov_angle > 70 && cherenkov_angle < 90){
    has_neutron ? r5.Fill(bsenergy, weight) : r2.Fill(bsenergy, weight);
  }
  
  
  return true;
}


bool MakeSpectralFitHistos::Finalise(){
  
  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "region_plot_"+weight_name+".root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (outfile == nullptr){
    throw std::runtime_error("PositronIdentificationCuts::Finalise - Couldn't open output file");
  }

  mc_energy_full.Write();
  bs_energy_bdt_full.Write();
  bs_energy_full.Write();
  bs_energy_bdt_mccuts_full.Write();

  r0.Write();
  r1.Write();
  r2.Write();
  r3.Write();
  r4.Write();
  r5.Write();
  
  return true;
}

void MakeSpectralFitHistos::GetReaders(){
  m_variables.Get("bdt_tree_reader_str", tree_reader_str);
  if (tree_reader_str.empty() || m_data->Trees.count(tree_reader_str) == 0){
    throw std::invalid_argument("no valid bdt_treereader specified!");
  }
  bdt_tree_reader_ptr = m_data->Trees.at(tree_reader_str);
  if (bdt_tree_reader_ptr == nullptr){
    throw std::runtime_error("couldn't get treereader");
  }
  tree_reader_str = "";
  m_variables.Get("rw_tree_reader_str", tree_reader_str);
  if (tree_reader_str.empty() || m_data->Trees.count(tree_reader_str) == 0){
    throw std::invalid_argument("no valid rw_treereader specified!");
  }
  rw_tree_reader_ptr = m_data->Trees.at(tree_reader_str);
  if (rw_tree_reader_ptr == nullptr){
    throw std::runtime_error("couldn't get treereader");
  }
  return;
}

void MakeSpectralFitHistos::GetWeightName(){
  m_variables.Get("weight_name", weight_name);
  if (weight_name == ""){weight_name = "unweighted";}
  return;
}

bool MakeSpectralFitHistos::HasExactlyOneNeutron(const basic_array<float*> likelihoods){
  double cut = -1;
  bool ok = m_variables.Get("likelihood_cut", cut);
  if (!ok || cut == -1){
    throw std::runtime_error("DefineSignalRegions::HasOneNeutron - Couldn't get neutron likelihood cut value!");
  }
  const int n_neutrons = std::count_if(likelihoods.begin(), likelihoods.end(), [cut](double l){return l > cut;});
  return n_neutrons == 1;
}
