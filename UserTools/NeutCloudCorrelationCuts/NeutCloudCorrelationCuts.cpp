#include "NeutCloudCorrelationCuts.h"

#include "TVector3.h"

NeutCloudCorrelationCuts::NeutCloudCorrelationCuts():Tool(){}


bool NeutCloudCorrelationCuts::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  if (!m_variables.Get("relic_reader_name", relic_reader_name)) return false;

  return true;
}


bool NeutCloudCorrelationCuts::Execute(){

  /*
    for each relic in a 60s window around muon, do the following:
    
    - calculate dx, dy, dz - dt already calculated
    - make appropriate cuts based on multiplicity of neutrons in the cloud
   */

  std::vector<double> muon_dir = std::vector<double>(skroot_mu_.mudir, skroot_mu_.mudir +3);
  std::vector<double> neutron_cloud_vertex;
  m_data->CStore.Get("neutron_cloud_vertex", neutron_cloud_vertex);
  
  std::vector<TVector3> coord_change_tensor = GetTensor(muon_dir, neutron_cloud_vertex);
  
  std::vector<int> matched_relic_entries = {};
  std::vector<double> matched_relic_tdiffs = {};
  
  int multiplicity = 0;
  m_data->CStore.Get("neutron_cloud_multiplicity", multiplicity);

  for (int i = 0; i < matched_relic_entries.size(); ++i){
    const double dt = matched_relic_tdiffs.at(i);

    m_data->getTreeEntry(relic_reader_name, matched_relic_entries.at(i));

    // old - regular sk coordinate system
    // new - azimuthal coordinate system with z aligning with muon track
    const TVector3 dr_old = TVector3(neutron_cloud_vertex.data()) - TVector3(skroot_lowe_.bsvertex);
    const std::vector<double> dr_squared_new = {
      pow(dr_old.Dot(coord_change_tensor.at(0)), 2),
      pow(dr_old.Dot(coord_change_tensor.at(1)), 2),
      pow(dr_old.Dot(coord_change_tensor.at(2)),2)
    };

    
    
  }
							 
  return true;
}


bool NeutCloudCorrelationCuts::Finalise(){

  return true;
}

std::vector<TVector3> NeutCloudCorrelationCuts::GetTensor(const std::vector<double>& mu_dir,
							  const std::vector<double>& nc_vert) const {

  TVector3 mu_dir_v3 = TVector3(mu_dir.data()),  nc_vert_v3 = TVector3(nc_vert.data());

  if (mu_dir_v3.Angle(nc_vert_v3) == 0){
    mu_dir_v3.SetX(0);
  }
  
  TVector3 z = mu_dir_v3;
  z.SetMag(1);

  TVector3 y = z.Cross(nc_vert_v3);
  y.SetMag(1);

  TVector3 x = y.Cross(z);
  
  return {x,y,z};
}
