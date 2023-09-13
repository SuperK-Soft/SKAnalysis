#include "PostReconstructionNeutronCloudSelection.h"

#include "TFile.h"

PostReconstructionNeutronCloudSelection::PostReconstructionNeutronCloudSelection():Tool(){}


bool PostReconstructionNeutronCloudSelection::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

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

  pre_bsgood_cut.Fill(skroot_lowe_.bsgood[1]);
  if (skroot_lowe_.bsgood[1] < 0.4){
    return true;
  }
  post_bsgood_cut.Fill(skroot_lowe_.bsgood[1]);

  pre_bsdirks_cut.Fill(skroot_lowe_.bsdirks);
  if (skroot_lowe_.bsdirks > 0.4){
    return true;
  }
  post_bsdirks_cut.Fill(skroot_lowe_.bsdirks);

  pre_bsn50_cut.Fill(skroot_lowe_.bsn50);
  if (skroot_lowe_.bsn50 < 24 || skroot_lowe_.bsn50 > 50){
    return true;
  }
  post_bsn50_cut.Fill(skroot_lowe_.bsn50);

  const float ldt  = CalculateLongitudinalDistance(skroot_mu_.muentpoint, skroot_mu_.mudir, skroot_lowe_.bsvertex);
  pre_ldt_cut.Fill(ldt);
  if (ldt  > 500){
    return true;
  }
  post_ldt_cut.Fill(ldt);

  std::vector<int> SLE_times;
  m_data->CStore.Get("SLE_times", SLE_times);
  
  if (neutrons.size() == SLE_times.size()){
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
  
  return true;
}


bool PostReconstructionNeutronCloudSelection::Finalise(){
  
  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "pre_recon_neut_cloud_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "UPDATE");
  if (outfile == nullptr){
    throw std::runtime_error("PostReconstructionNeutronCloudSelection::Finalise - Couldn't open output file");
  }

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
