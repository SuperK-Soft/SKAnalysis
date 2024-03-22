#include "NeutCloudCorrelationCuts.h"

#include "TVector3.h"
#include "TH1D.h"
#include "TFile.h"

#include "MTreeReader.h"

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

  GetTreeReaders();
  
  return true;
}

bool NeutCloudCorrelationCuts::Execute(){

  bool ok = relic_tree_reader->Get("MatchedOutEntryNums", relicMatchedEntryNums);
  if (!ok){throw std::runtime_error("NeutCloudCorrelationCuts::Execute - couldn't retrieve matched entries");}
  ok = relic_tree_reader->Get("MatchedTimeDiffs", relicTimeDiffs);
  if (!ok){throw std::runtime_error("NeutCloudCorrelationCuts::Execute - couldn't retrieve matched time differences");}  

  if (relicMatchedEntryNums->empty()){return true;}

  for (size_t i = 0; i < relicMatchedEntryNums->size(); ++i){

    const int match_idx = relicMatchedEntryNums->at(i);
    std::cout << "match_idx: " << match_idx << std::endl;
    const double dt = relicTimeDiffs->at(i);
    std::cout << "matched time diff" << std::endl;

    ok = m_data->getTreeEntry(cloud_tree_reader_str, match_idx);
    if (!ok){throw std::runtime_error("NeutCloudCorrelationCuts::Execute - failed to retrieve cloud file entry");} 

    std::vector<double> muon_dir = {}; //std::vector<double>(skroot_mu_.mudir, skroot_mu_.mudir +3);
    std::vector<double> neutron_cloud_vertex = {};
    int* multiplicity = nullptr;

    cloud_tree_reader->Get("neutron_cloud_vertex", neutron_cloud_vertex);
    cloud_tree_reader->Get("neutron_cloud_multiplicity", *multiplicity);
    cloud_tree_reader->Get("muon_dir", muon_dir);
  
    std::vector<TVector3> coord_change_tensor = GetTensor(muon_dir, neutron_cloud_vertex);
     
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
    if (*multiplicity == 2){
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

    if (*multiplicity == 3){
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
  
    if((*multiplicity == 4) || (*multiplicity == 5)){
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
    
    if ((*multiplicity > 6) && (*multiplicity < 9)){
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
    
    if (*multiplicity >= 10){
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
    if ((*multiplicity > 2) && ((std::abs(dt) < 0.1 && dL < 1200) || (std::abs(dt) < 1 && dL < 800))){
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

void NeutCloudCorrelationCuts::GetTreeReaders(){
  std::string tree_reader_str = "";
  m_variables.Get("relic_TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("CalculateNeutronCloudVertex::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  relic_tree_reader = m_data->Trees.at(tree_reader_str);

  m_variables.Get("cloud_TreeReader", tree_reader_str);
  cloud_tree_reader_str = tree_reader_str;
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("CalculateNeutronCloudVertex::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  cloud_tree_reader = m_data->Trees.at(tree_reader_str);

}
