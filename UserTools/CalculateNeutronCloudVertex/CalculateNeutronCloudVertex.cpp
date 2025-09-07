#include "CalculateNeutronCloudVertex.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

#include <algorithm>

CalculateNeutronCloudVertex::CalculateNeutronCloudVertex():Tool(){}

bool CalculateNeutronCloudVertex::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetTreeReader();

  CreateOutputFile();
  
  return true;
}

bool CalculateNeutronCloudVertex::Execute(){
  
  int N_SLE = 0;
  m_data->CStore.Get("N_SLE", N_SLE);
  // just for record keeping, number of SLE triggers (before neutron reconstruction and cuts)
  N_SLE_plot.Fill(N_SLE);
  
  std::vector<NeutronInfo> neutrons = {};
  m_data->CStore.Get("event_neutrons", neutrons);
  if (N_SLE==0 || neutrons.empty()){
    // note that if SLE_Search Tool finds no triggers, the subtoolchain won't have run,
    // so neutrons will not have been cleared by PostReconstructionNeutronCloudSelection.
    // so it's ok if neutrons is not empty but N_SLE=0
    
    // multiplicity 0, set neutron cloud vertex to (0,0,0)
    mult = 0;
    mult_plot.Fill(mult);
    std::fill(neutron_cloud_vertex.begin(), neutron_cloud_vertex.end(), 0);
    nvc_tree_ptr->Fill();
    return true;
  }
  
  // record number of reconstructed neutrons after cuts
  mult = neutrons.size();
  mult_plot.Fill(mult);
  
  // calculate mean position of all neutron vertices
  for (const auto& neutron : neutrons){
    for (int dim = 0; dim < 3; ++dim){
      neutron_cloud_vertex.at(dim) += neutron.bs_vertex.at(dim) / neutrons.size();
    }
  }
  
  // distance from neutron cloud to ... muon?
  // i suppose this could be interesting...
  double closestappr=ClosestApproach(neutron_cloud_vertex);
  dist_to_mu_plot.Fill(closestappr);
  //std::cout<<"distance muon to neutron cloud: "<<closestappr<<std::endl;
  
  // fill output tree
  nvc_tree_ptr->Fill();

  // intermittent write so we don't lose everything if toolchain crashes
  if (MU_tree_reader->GetEntryNumber() % 1000 == 0){
    nvc_file_ptr->cd();
    nvc_tree_ptr->Write("*",TObject::kOverwrite);
  }
  
  return true;
}

bool CalculateNeutronCloudVertex::Finalise(){

  plotfile->cd();
  mult_plot.Write();
  N_SLE_plot.Write();
  dist_to_mu_plot.Write();
  plotfile->Close();

  nvc_file_ptr->cd();
  nvc_tree_ptr->Write("*",TObject::kOverwrite);
  nvc_file_ptr->Close();
  
  return true;
}

void CalculateNeutronCloudVertex::GetTreeReader(){
  std::string tree_reader_str = "";
  m_variables.Get("MU_TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("CalculateNeutronCloudVertex::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  MU_tree_reader = m_data->Trees.at(tree_reader_str);
}


double CalculateNeutronCloudVertex::ClosestApproach(const std::vector<double>& vertex) {
  // skroot_mu_ not populated unless we do so, just use MU branch
  //const std::vector<double> muon_ent(skroot_mu_.muentpoint, skroot_mu_.muentpoint + 3);
  //const std::vector<double> muon_dir(skroot_mu_.mudir, skroot_mu_.mudir + 3);
  
  MuInfo* MU_ptr = nullptr;
  MU_tree_reader->Get("MU", MU_ptr);
  if(MU_ptr == nullptr){
    throw std::runtime_error("CalculateNeutronCloudVertex::Execute: failed to get MU branch from tree reader");
  }
  const int muboy_idx = MU_ptr->muinfo[7];
  const std::vector<double> muon_dir(MU_ptr->muboy_dir, MU_ptr->muboy_dir+3);
  const std::vector<double> muon_ent(MU_ptr->muboy_entpos[muboy_idx],MU_ptr->muboy_entpos[muboy_idx]+3);
  //std::cout<<"muon enters at "<<muon_ent[0]<<", "<<muon_ent[1]<<", "<<muon_ent[2]
  //         <<" with direction "<<muon_dir[0]<<", "<<muon_dir[1]<<", "<<muon_dir[2]<<std::endl;

  std::vector<double> diff = std::vector<double>(3, 0);
  for (int i = 0; i < 3; ++i){diff.at(i) = vertex.at(i) - muon_ent.at(i);}

  std::vector<double> proj = std::vector<double>(3, 0);
  std::vector<double> dist_vec = std::vector<double>(3, 0);

  const double diff_mudir_ip = std::inner_product(diff.begin(), diff.end(), muon_dir.begin(), 0);
  const double mudir_muent_ip = std::inner_product(muon_dir.begin(), muon_dir.end(), muon_ent.begin(), 0);
  
  for (int i = 0; i < 3; ++i){
    proj.at(i) = (diff_mudir_ip / mudir_muent_ip) * muon_dir.at(i);
    dist_vec.at(i) = vertex.at(i) - proj.at(i) - muon_ent.at(i);
  }
  
  return sqrt(std::inner_product(dist_vec.begin(), dist_vec.end(), dist_vec.begin(), 0));
}

void CalculateNeutronCloudVertex::CreateOutputFile(){
  
  // output tree of results - one neutron cloud vertex for each muon tree entry
  // we will leave applying cuts to the relic tree for a subsequent toolchain
  // which will process the relic tree and for each relic, check the distance
  // to each matched muon's neutron cloud vertex.
  std::string nvc_file_str = "";
  m_variables.Get("nvc_file_str", nvc_file_str);
  if (nvc_file_str.empty()){
    throw std::runtime_error("CalculateNeutronCloudVertex::CreateOutputFile - no output file specified!");
  }
  nvc_file_ptr = TFile::Open(nvc_file_str.c_str(), "RECREATE");
  
  nvc_tree_ptr = new TTree("neutron_cloud_info", "neutron_cloud_info");
  nvc_tree_ptr->Branch("neutron_cloud_multiplicity", &mult);
  nvc_tree_ptr->Branch("neutron_cloud_vertex", &neutron_cloud_vertex);
  nvc_tree_ptr->Branch("nevsk",&skhead_.nevsk); // to check Tree alignment
  
  
  // for plots
  std::string plotfile_name = "";
  bool ok = m_variables.Get("plotfile_name", plotfile_name);
  if (!ok || plotfile_name.empty()){ plotfile_name = "calculateneutroncloudvertex_out.root";}

  plotfile = TFile::Open(plotfile_name.c_str(), "RECREATE");
  if (plotfile == nullptr){
    throw std::runtime_error("CalculateNeutronCloudVertex::Finalise - Couldn't open output file");
  }
  
  N_SLE_plot = TH1D("N_SLE_plot", "number of SLE triggers after muon", 20, 0, 20);
  mult_plot = TH1D("mult_plot", "multiplcity of neutron cloud;multiplcity", 20, 0, 20);
  dist_to_mu_plot = TH1D("dist_to_mu_plot", "distance to muon plot; distance [cm]", 100, 100, 100);
  
  return;
}
