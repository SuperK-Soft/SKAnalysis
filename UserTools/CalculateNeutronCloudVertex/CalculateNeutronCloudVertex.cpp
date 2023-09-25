#include "CalculateNeutronCloudVertex.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TH1D.h"

CalculateNeutronCloudVertex::CalculateNeutronCloudVertex():Tool(){}

bool CalculateNeutronCloudVertex::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetTreeReader();

  mult_plot = TH1D("mult_plot", "multiplcity of neutron cloud;multiplcity", 20, 0, 20);
  
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
