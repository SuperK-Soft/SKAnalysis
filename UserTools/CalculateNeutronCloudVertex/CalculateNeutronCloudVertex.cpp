#include "CalculateNeutronCloudVertex.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TH1D.h"

#include <algorithm>

CalculateNeutronCloudVertex::CalculateNeutronCloudVertex():Tool(){}

bool CalculateNeutronCloudVertex::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetTreeReader();

  mult_plot = TH1D("mult_plot", "multiplcity of neutron cloud;multiplcity", 20, 0, 20);
  dist_to_mu_plot = TH1D("dist_to_mu_plot", "distance to muon plot; distance [cm]", 100, 100, 100);
  
  return true;
}

bool CalculateNeutronCloudVertex::Execute(){

  std::vector<NeutronInfo> neutrons = {};
  m_data->CStore.Get("event_neutrons", neutrons);

  if (neutrons.empty()){
    bool skip = true;
    m_data->CStore.Set("Skip", skip);
    return true;
  }

  std::vector<double> neutron_cloud_vertex = {}; 
  
  for (const auto& neutron : neutrons){
    for (int dim = 0; dim < 3; ++dim){
      neutron_cloud_vertex.at(dim) = neutron.bs_vertex.at(dim) / neutrons.size();
    }
  }
  
  mult = neutrons.size();
  mult_plot.Fill(mult);


  dist_to_mu_plot.Fill(ClosestApproach(neutron_cloud_vertex));
  
  LOWE_tree_reader->GetTree()->Branch("neutron_cloud_vertex", &neutron_cloud_vertex);
  LOWE_tree_reader->GetTree()->Branch("neutron_multiplicity", &mult);
  
  return true;
}

bool CalculateNeutronCloudVertex::Finalise(){

  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "calculateneutroncloudvertex_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "UPDATE");
  if (outfile == nullptr){
    throw std::runtime_error("CalculateNeutronCloudVertex::Finalise - Couldn't open output file");
  }

  mult_plot.Write();
  dist_to_mu_plot.Write();
  
  return true;
}

void CalculateNeutronCloudVertex::GetTreeReader(){
  std::string tree_reader_str = "";
  m_variables.Get("LOWE_TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("CalculateNeutronCloudVertex::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  LOWE_tree_reader = m_data->Trees.at(tree_reader_str);
}


double CalculateNeutronCloudVertex::ClosestApproach(const std::vector<double>& vertex) const {
  const std::vector<double> muon_ent = std::vector<double>(skroot_mu_.muentpoint, skroot_mu_.muentpoint + 3);
  const std::vector<double> muon_dir = std::vector<double>(skroot_mu_.mudir, skroot_mu_.mudir + 3);

  std::vector<double> diff = std::vector<double>(3, 0);
  for (int i = 0; i < 3; ++i){diff.at(i) = vertex.at(i) - muon_ent.at(i);}

  std::vector<double> proj = std::vector<double>(3, 0);
  std::vector<double> dist_vec = std::vector<double>(3, 0);

  const double diff_mudir_ip = std::inner_product(diff.begin(), diff.end(), muon_dir.begin(), 0);
  const double mudir_muent_ip = std::inner_product(muon_dir.begin(), muon_dir.end(), muon_dir.end(), 0);
  
  for (int i = 0; i < 3; ++i){
    proj.at(i) = (diff_mudir_ip / mudir_muent_ip) * muon_dir.at(i);
    dist_vec.at(i) = vertex.at(i) - proj.at(i) - muon_ent.at(i);
  }
  
  return sqrt(std::inner_product(dist_vec.begin(), dist_vec.begin(), dist_vec.end(), 0));
}
