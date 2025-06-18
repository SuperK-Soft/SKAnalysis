#include "PostReconstructionNeutronCloudSelection.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TFile.h"

PostReconstructionNeutronCloudSelection::PostReconstructionNeutronCloudSelection():Tool(){}


bool PostReconstructionNeutronCloudSelection::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  bsenergy_plot = TH1D("bonsai_energy", "bonsai_energy;bsenergy", 100, 0, 0);
  
  pre_bsgood_cut = TH1D("pre_bsgood_cut", "pre_bsgood_cut;bsgood", 100, 0, 1);
  pre_bsdirks_cut = TH1D("pre_bsdirks_cut", "pre_bsdirks_cut;bsdirks", 100, 0, 1);
  pre_bsn50_cut = TH1D("pre_bsn50_cut", "pre_bsn50_cut;bsn50", 100, 0, 100);
  pre_ldt_cut = TH1D("pre_ldt_cut", "pre_ldt_cut;ldt", 100, 0, 1000);

  post_bsgood_cut = TH1D("post_bsgood_cut", "post_bsgood_cut;bsgood", 100, 0, 1);
  post_bsdirks_cut = TH1D("post_bsdirks_cut", "post_bsdirks_cut;bsdirks", 100, 0, 1);
  post_bsn50_cut = TH1D("post_bsn50_cut", "post_bsn50_cut;bsn50", 100, 0, 100);
  post_ldt_cut = TH1D("post_ldt_cut", "post_ldt_cut;ldt", 100, 0, 1000);

  return true;
}


bool PostReconstructionNeutronCloudSelection::Execute(){

  /* 
     Criteria on Reconstructed SLE Triggers are as follows:

     1. Time difference between SLE and muon is within [35us, 535us] (done pre-recon)
     2. bsgood > 0.4 
     3. dirks < 0.4
     4. 24 <= bsn50 <= 50
     5. Distance between SLE and muon is < 500cm
     
  */

  std::vector<int> SLE_times;
  m_data->CStore.Get("SLE_times", SLE_times);
  N = SLE_times.size();
  
  
  std::cout << "this neutron candidate has:" << std::endl;
  std::cout << "skroot_lowe_.bsenergy: " << skroot_lowe_.bsenergy << std::endl;
  std::cout << "skroot_lowe_.bsgood[1]: " << skroot_lowe_.bsgood[1] << std::endl;
  std::cout << "skroot_lowe_.bsdirks: " << skroot_lowe_.bsdirks << std::endl;
  std::cout << "skroot_lowe_.bsn50: " << skroot_lowe_.bsn50 << std::endl;

  // MuInfo* my_muinfo = new MuInfo();
  // tree_reader_ptr->Get("MU", my_muinfo);
  // if (my_muinfo == nullptr){
  //   throw std::runtime_error("PostReconstructionNeutronCloudSelection::Execute: failed to get MU branch from tree reader");
  // }
  // const int muboy_idx = my_muinfo->muinfo[7];

  
  // std::cout << "muboy_idx: " << muboy_idx << std::endl;
  // const float ldt = CalculateLongitudinalDistance(skroot_mu_.muboy_entpos[muboy_idx], skroot_mu_.muboy_dir, skroot_lowe_.bsvertex);
  // pre_ldt_cut.Fill(ldt);
  // std::cout << "ldt: " << ldt << std::endl;

  if (skroot_lowe_.bsenergy > 9000){
    std::cout << "neutron candidate was not reconstructed properly - rejected" << std::endl;
    --N;
    return true;
  }
  bsenergy_plot.Fill(skroot_lowe_.bsenergy);
  
  pre_bsgood_cut.Fill(skroot_lowe_.bsgood[1]);
  if (skroot_lowe_.bsgood[1] < 0.4){
    --N;
    std::cout << "neutron candidate does not meet bsgood cut of > 0.4" << std::endl;
    return true;
  }
  post_bsgood_cut.Fill(skroot_lowe_.bsgood[1]);

  pre_bsdirks_cut.Fill(skroot_lowe_.bsdirks);
  if (skroot_lowe_.bsdirks > 0.4){
    --N;
    std::cout << "neutron candidate does not meet bsdirks cut of < 0.4" << std::endl;
    return true;
  }
  post_bsdirks_cut.Fill(skroot_lowe_.bsdirks);

  pre_bsn50_cut.Fill(skroot_lowe_.bsn50);
  if (skroot_lowe_.bsn50 < 24 || skroot_lowe_.bsn50 > 50){
    --N;
    std::cout << "neutron candidate does not meet bsn50 > 24 && bsn50 < 50" << std::endl;
    return true;
  }
  post_bsn50_cut.Fill(skroot_lowe_.bsn50);

  
  
  // if (ldt > 500){
  //   std::cout << "neutron candidate does not meet ldt cut of < 500" << std::endl;
  //   return true;
  // }
  // post_ldt_cut.Fill(ldt);

  
  
  
  if (neutrons.size() == N){
    m_data->CStore.Set("event_neutrons", neutrons);
    neutrons.clear();
  } else {
    NeutronInfo n;
    n.bs_goodness = skroot_lowe_.bsgood[1];
    n.bs_dirks = skroot_lowe_.bsdirks;
    n.bsn50 = skroot_lowe_.bsn50;
    n.bs_vertex = {skroot_lowe_.bsvertex[0], skroot_lowe_.bsvertex[1], skroot_lowe_.bsvertex[2]};
    neutrons.push_back(n);
  }

  //  delete my_muinfo; my_muinfo = nullptr;
  
  return true;
}


bool PostReconstructionNeutronCloudSelection::Finalise(){
  
  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "post_recon_neut_cloud_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (outfile == nullptr){
    throw std::runtime_error("PostReconstructionNeutronCloudSelection::Finalise - Couldn't open output file");
  }

  bsenergy_plot.Write();
  
  pre_bsgood_cut.Write();
  pre_bsdirks_cut.Write();
  pre_bsn50_cut.Write();
  pre_ldt_cut.Write();

  post_bsgood_cut.Write();
  post_bsdirks_cut.Write();
  post_bsn50_cut.Write();
  post_ldt_cut.Write();
  
  return true;
}

float PostReconstructionNeutronCloudSelection::CalculateLongitudinalDistance(float* muon_entry_point,
									     float* muon_direction,
									     float* bs_vertex) const {
  float result = 0;
  for (int i = 0; i < 3; ++i){
    result += pow( ( muon_entry_point[i] - bs_vertex[i] ) * muon_direction[(i+1)%3] -
		   ( muon_entry_point[(i+1)%3] - bs_vertex[(i+1)%3] ) * muon_direction[i], 2);
  }
  return sqrt(result);
}

void PostReconstructionNeutronCloudSelection::GetReader(){
  std::string tree_reader_str = "";
  m_variables.Get("tree_reader_str", tree_reader_str);
  if (tree_reader_str.empty() || m_data->Trees.count(tree_reader_str) == 0){
    throw std::invalid_argument("no valid treereader specified!");
  }
  tree_reader_ptr = m_data->Trees.at(tree_reader_str);
  if (tree_reader_ptr == nullptr){
    throw std::runtime_error("couldn't get treereader");
  }
  return;
}
