#include "CalculateNeutronCloudVertex.h"

#include "NeutronInfo.h"

CalculateNeutronCloudVertex::CalculateNeutronCloudVertex():Tool(){}


bool CalculateNeutronCloudVertex::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

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
    double weighting = GetWeighting(neutron);
    for (int dim = 0; dim < 3; ++dim){
      neutron_cloud_vertex.at(dim) = neutron.bs_vertex.at(dim) * weighting / neutrons.size(); //this is mean position, but we could do median etc...
    }
  }

  m_data->CStore.Set("neutron_cloud_vertex", neutron_cloud_vertex);
  int mult = neutrons.size();
  m_data->CStore.Set("neutron_cloud_multiplicity", mult);
  
  return true;
}


bool CalculateNeutronCloudVertex::Finalise(){

  return true;
}

double CalculateNeutronCloudVertex::GetWeighting(const NeutronInfo n) const {
  return 1; // can we think of a good weighting?
}
