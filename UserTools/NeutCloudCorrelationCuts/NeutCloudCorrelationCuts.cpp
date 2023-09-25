#include "NeutCloudCorrelationCuts.h"

#include "TVector3.h"
#include "TH1D.h"
#include "TFile.h"

NeutCloudCorrelationCuts::NeutCloudCorrelationCuts():Tool(){}

bool NeutCloudCorrelationCuts::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  if (!m_variables.Get("relic_reader_name", relic_reader_name)) return false;

  pre_sample_total_dt = TH1D("pre_sample_total_dt", "total;abs. time from relic candidate [s]", 100, 0, 60);
  pre_sample_m2_dt = TH1D("pre_sample_m2_dt", "m = 2;abs. time from relic candidate [s]", 100, 0, 60);
  pre_sample_m3_dt = TH1D("pre_sample_m3_dt", "m = 3;abs. time from relic candidate [s]", 100, 0, 60);
  pre_sample_m45_dt = TH1D("pre_sample_m45_dt", "m = 4,5;abs. time from relic candidate [s]", 100, 0, 60);
  pre_sample_m69_dt = TH1D("pre_sample_m69_dt", "m = 6-9;abs. time from relic candidate [s]", 100, 0, 60);
  pre_sample_m10_dt = TH1D("pre_sample_m10_dt", "m = 10+;abs. time from relic candidate [s]", 100, 0, 60);  

  post_sample_total_dt = TH1D("post_sample_total_dt", "total;abs. time from relic candidate [s]", 100, 0, 60);
  post_sample_m2_dt = TH1D("post_sample_m2_dt", "m = 2;abs. time from relic candidate [s]", 100, 0, 60);
  post_sample_m3_dt = TH1D("post_sample_m3_dt", "m = 3;abs. time from relic candidate [s]", 100, 0, 60);
  post_sample_m45_dt = TH1D("post_sample_m45_dt", "m = 4,5;abs. time from relic candidate [s]", 100, 0, 60);
  post_sample_m69_dt = TH1D("post_sample_m69_dt", "m = 6-9;abs. time from relic candidate [s]", 100, 0, 60);
  post_sample_m10_dt = TH1D("post_sample_m10_dt", "m = 10+;abs. time from relic candidate [s]", 100, 0, 60);  

  pre_sample_total_dl = TH1D("pre_sample_total_dl", "total;distance from relic candidate [cm]", 100, 0, 5000);
  pre_sample_m2_dl = TH1D("pre_sample_m2_dl", "m = 2;distance from relic candidate [cm]", 100, 0, 5000);
  pre_sample_m3_dl = TH1D("pre_sample_m3_dl", "m = 3;distance from relic candidate [cm]", 100, 0, 5000);
  pre_sample_m45_dl = TH1D("pre_sample_m45_dl", "m = 4,5;distance from relic candidate [cm]", 100, 0, 5000);
  pre_sample_m69_dl = TH1D("pre_sample_m69_dl", "m = 6-9;distance from relic candidate [cm]", 100, 0, 5000);
  pre_sample_m10_dl = TH1D("pre_sample_m10_dl", "m = 10+;distance from relic candidate [cm]", 100, 0, 5000);  

  post_sample_total_dl = TH1D("post_sample_total_dl", "total;distance from relic candidate [cm]", 100, 0, 5000);
  post_sample_m2_dl = TH1D("post_sample_m2_dl", "m = 2;distance from relic candidate [cm]", 100, 0, 5000);
  post_sample_m3_dl = TH1D("post_sample_m3_dl", "m = 3;distance from relic candidate [cm]", 100, 0, 5000);
  post_sample_m45_dl = TH1D("post_sample_m45_dl", "m = 4,5;distance from relic candidate [cm]", 100, 0, 5000);
  post_sample_m69_dl = TH1D("post_sample_m69_dl", "m = 6-9;distance from relic candidate [cm]", 100, 0, 5000);
  post_sample_m10_dl = TH1D("post_sample_m10_dl", "m = 10+;distance from relic candidate [cm]", 100, 0, 5000);  

  return true;
}

bool NeutCloudCorrelationCuts::Execute(){

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
    // new - cloud coordinate system with z aligning with muon track
    const TVector3 dr_old = TVector3(neutron_cloud_vertex.data()) - TVector3(skroot_lowe_.bsvertex);
    const std::vector<double> dr_squared_new = {
      pow(dr_old.Dot(coord_change_tensor.at(0)), 2),
      pow(dr_old.Dot(coord_change_tensor.at(1)), 2),
      pow(dr_old.Dot(coord_change_tensor.at(2)),2)
    };

    const double dL = std::inner_product(dr_squared_new.begin(), dr_squared_new.end(), dr_squared_new.begin(), 0);

    if (dt < 0){
      pre_sample_total_dl.Fill(dL);
      pre_sample_total_dt.Fill(dt);
    } else {
      post_sample_total_dl.Fill(dL);
      post_sample_total_dt.Fill(dt);
    }
       
    const auto ellipse = [dr_squared_new](double a, double b, double c){
      return
	((dr_squared_new.at(0) / a) +
	 (dr_squared_new.at(1) / b) +
	 (dr_squared_new.at(2) / c));
    };
      
    //ellipse cuts
    if (multiplicity == 2){
      if (dt < 0){
	pre_sample_m2_dl.Fill(dL);
	pre_sample_m2_dt.Fill(dt);
      } else {
	post_sample_m2_dl.Fill(dL);
	post_sample_m2_dt.Fill(dt);
      }
      if ((std::abs(dt) < 30) && ellipse(40000, 40000, 160000) < 1.2){
	SkipEntry();
      }
    }

    if (multiplicity == 3){
      if (dt < 0){
	pre_sample_m3_dl.Fill(dL);
	pre_sample_m3_dt.Fill(dt);
      } else {
	post_sample_m3_dl.Fill(dL);
	post_sample_m3_dt.Fill(dt);
      }      
      if ((std::abs(dt) < 60) && ellipse(60000, 60000, 250000) < 1.2){
	SkipEntry();
      }
    }
  
    if((multiplicity == 4) || (multiplicity == 5)){
      if (dt < 0){
	pre_sample_m45_dl.Fill(dL);
	pre_sample_m45_dt.Fill(dt);
      } else {
	post_sample_m45_dl.Fill(dL);
	post_sample_m45_dt.Fill(dt);
      }
      if ((std::abs(dt) < 60) && ellipse(120000, 120000, 302500) < 1.2){
	SkipEntry();
      }
    }
    
    if ((multiplicity > 6) && (multiplicity < 9)){
      if (dt < 0){
	pre_sample_m69_dl.Fill(dL);
	pre_sample_m69_dt.Fill(dt);
      } else {
	post_sample_m69_dl.Fill(dL);
	post_sample_m69_dt.Fill(dt);
      }
      if ((std::abs(dt) < 60) && ellipse(200000, 200000, 422500) < 1.2){
	SkipEntry();
      }
    }
    
    if (multiplicity >= 10){
      if (dt < 0){
	pre_sample_m10_dl.Fill(dL);
	pre_sample_m10_dt.Fill(dt);
      } else {
	post_sample_m10_dl.Fill(dL);
	post_sample_m10_dt.Fill(dt);
      }
      if ((std::abs(dt) < 60) && ellipse(250000, 250000, 490000) < 1.2){
	SkipEntry();
      }
    }

    //box cuts
    if ((multiplicity > 2) && ((std::abs(dt) < 0.1 && dL < 1200) || (std::abs(dt) < 1 && dL < 800))){
      SkipEntry();
    }
    
  }
  
  return true;
}


bool NeutCloudCorrelationCuts::Finalise(){
  
  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "pre_recon_neut_cloud_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "UPDATE");
  if (outfile == nullptr){
    throw std::runtime_error("PostReconstructionNeutronCloudSelection::Finalise - Couldn't open output file");
  }
  
  pre_sample_total_dt.Write();
  pre_sample_m2_dt.Write();
  pre_sample_m3_dt.Write();
  pre_sample_m45_dt.Write();
  pre_sample_m69_dt.Write();
  pre_sample_m10_dt.Write();

  post_sample_total_dt.Write();
  post_sample_m2_dt.Write();
  post_sample_m3_dt.Write();
  post_sample_m45_dt.Write();
  post_sample_m69_dt.Write();
  post_sample_m10_dt.Write();

  pre_sample_total_dl.Write();
  pre_sample_m2_dl.Write();
  pre_sample_m3_dl.Write();
  pre_sample_m45_dl.Write();
  pre_sample_m69_dl.Write();
  pre_sample_m10_dl.Write();

  post_sample_total_dl.Write();
  post_sample_m2_dl.Write();
  post_sample_m3_dl.Write();
  post_sample_m45_dl.Write();
  post_sample_m69_dl.Write();
  post_sample_m10_dl.Write();
  
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
 
void NeutCloudCorrelationCuts::SkipEntry(){
  bool skip = true;
  m_data->CStore.Set("Skip", skip);
}
