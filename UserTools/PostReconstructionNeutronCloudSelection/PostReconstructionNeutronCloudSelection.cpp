#include "PostReconstructionNeutronCloudSelection.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TFile.h"

PostReconstructionNeutronCloudSelection::PostReconstructionNeutronCloudSelection():Tool(){}

// helper
namespace{
        //enum reject { fail } reject;
        class reject {
                public:
                reject(const std::string& m) : msg(m){};
                std::string msg;
                const std::string what(){ return msg; };
        };
};

bool PostReconstructionNeutronCloudSelection::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  
  GetReader();
  
  // TODO make cut thresholds input variables
  
  // add a cut to the selector if being used
  get_ok = m_variables.Get("selectorName", selectorName);
  if(get_ok){
    // n.b. if you use PrintMTreeSelections on the cuts file the printed number of 'passing events'
    // only accounts for *unique input Tree entries*. Since each entry has many SLE triggers
    // this is not the number of passing neutrons. The entry count of the distros is, though.
    // type 1/2 cuts aren't quite right as we have in input branch to use as an index.
    // honestly this is redundant with the histograms below.
    m_data->AddCut(selectorName, "bsvertex", "bsvertex[0] != 9999",false);
    m_data->AddCut(selectorName, "bsgood", "bsgood > 0.4",true,0.4,1.0);
    m_data->AddCut(selectorName, "bsdirks", "bsdirks < 0.4",true,0,0.4);
    m_data->AddCut(selectorName, "bsn50", "24 > bsn50 > 50",true,24,50);
    m_data->AddCut(selectorName, "dist", "distance to muon point of max energy deposition < 500",true,0,500);
  }
  
    
  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "post_recon_neut_cloud_out.root";}
  
  outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (outfile == nullptr){
    throw std::runtime_error("PostReconstructionNeutronCloudSelection::Initialise - Couldn't open output file");
  }
  
  bsenergy_plot = TH1D("bonsai_energy", "bonsai_energy;bsenergy", 100, 0, 0);
  h_n_neutrons = TH1D("h_n_neutrons","events:n_neutrons",50,0,50);
  pre_ldt_cut2 = TH1D("pre_ldt_cut2", "pre_ldt_cut2;ldt", 100, 0, 5000); // transverse distance only
  
  pre_bsgood_cut = TH1D("pre_bsgood_cut", "pre_bsgood_cut;bsgood", 100, 0, 1);
  pre_bsdirks_cut = TH1D("pre_bsdirks_cut", "pre_bsdirks_cut;bsdirks", 100, 0, 1);
  pre_bsn50_cut = TH1D("pre_bsn50_cut", "pre_bsn50_cut;bsn50", 100, 0, 100);
  pre_ldt_cut = TH1D("pre_ldt_cut", "pre_ldt_cut;ldt", 100, 0, 1000);  // total distance

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
  
  if(N==0){
    // new event from parent toolchain
    std::vector<int> SLE_times;
    m_data->CStore.Get("SLE_times", SLE_times);
    N = SLE_times.size();
    Log(m_unique_name+": next event had "+toString(N)+" subtriggers",v_debug,m_verbose);
    if(N==0){
      // prevent carry-over of previous event neutrons
      // this might happen if pre-reconstruction selection nukes all SLE triggers
      m_data->CStore.Set("event_neutrons", neutrons);
      return true;
    }
  }
  
  try {
    
    Log(m_unique_name+": this neutron candidate has\nskroot_lowe_.bsvertex: ("
       +toString(skroot_lowe_.bsvertex[0])+", "+toString(skroot_lowe_.bsvertex[1])+", "+toString(skroot_lowe_.bsvertex[2])+")"
       +"\nskroot_lowe_.bsenergy: "+toString(skroot_lowe_.bsenergy)
       +"\nskroot_lowe_.bsgood[1]: "+toString(skroot_lowe_.bsgood[1])
       +"\nskroot_lowe_.bsdirks: "+toString(skroot_lowe_.bsdirks)
       +"\nskroot_lowe_.bsn50: "+toString(skroot_lowe_.bsn50),v_debug,m_verbose);

    /* energy is reconstructed independently of vertex, both can fail independently.
    // but SK-4/6 analyses did not cut on bsenergy, so  maybe we don't care
    if (skroot_lowe_.bsenergy > 9000){
      throw reject("neutron energy was not reconstructed properly");
    }
    */
    bsenergy_plot.Fill(skroot_lowe_.bsenergy);
    
    // cut on bad vertex reconstruction
    if(skroot_lowe_.bsvertex[0]==9999){
      throw reject("neutron vertex was not reconstructed properly");
    }
    if(!selectorName.empty()) m_data->AddPassingEvent(selectorName, "bsvertex");
    
    // bsgood
    pre_bsgood_cut.Fill(skroot_lowe_.bsgood[1]);
    if(!selectorName.empty()) m_data->ApplyCut(selectorName, "bsgood", skroot_lowe_.bsgood[1]);
    if (skroot_lowe_.bsgood[1] < 0.4){
      throw reject("neutron candidate does not meet bsgood cut of > 0.4");
    }
    post_bsgood_cut.Fill(skroot_lowe_.bsgood[1]);
    
    // bsdirks
    pre_bsdirks_cut.Fill(skroot_lowe_.bsdirks);
    if(!selectorName.empty()) m_data->ApplyCut(selectorName, "bsdirks", skroot_lowe_.bsdirks);
    if (skroot_lowe_.bsdirks > 0.4){
      throw reject("neutron candidate does not meet bsdirks cut of < 0.4");
    }
    post_bsdirks_cut.Fill(skroot_lowe_.bsdirks);
    
    // bsn50
    pre_bsn50_cut.Fill(skroot_lowe_.bsn50);
    if(!selectorName.empty()) m_data->ApplyCut(selectorName, "bsn50", skroot_lowe_.bsn50);
    if (skroot_lowe_.bsn50 < 24 || skroot_lowe_.bsn50 > 50){
      throw reject("neutron candidate does not meet bsn50 > 24 && bsn50 < 50");
    }
    post_bsn50_cut.Fill(skroot_lowe_.bsn50);
    
    // distance to muon
    // n.b. skroot_mu_ is not populated unless we call skroot_get_mu_ so just use MU branch pointer.
    MuInfo* MU_ptr = nullptr;
    tree_reader_ptr->Get("MU", MU_ptr);
    if(MU_ptr == nullptr){
      throw std::runtime_error("PostReconstructionNeutronCloudSelection::Execute: failed to get MU branch from tree reader");
    }
    const int muboy_idx = MU_ptr->muinfo[7];
    Log(m_unique_name+": muboy_idx: "+toString(muboy_idx),v_debug,m_verbose);
    float* muon_dedx = MU_ptr->muboy_dedx;
    basic_array<float> muboy_dir(MU_ptr->muboy_dir);                   // direction from muboy
    basic_array<float[10][4]> muboy_entrypoint(MU_ptr->muboy_entpos);  // entry points from muboy
    float* muon_entrypoint = const_cast<float*>(muboy_entrypoint[muboy_idx].data());
    float* muon_direction = const_cast<float*>(muboy_dir.data());
    //std::cout<<"muon enters at "<<muon_entrypoint[0]<<", "<<muon_entrypoint[1]<<", "<<muon_entrypoint[2]
    //         <<" with direction "<<muon_direction[0]<<", "<<muon_direction[1]<<", "<<muon_direction[2]<<std::endl;
    
    const float dist = CalculateDistanceToMuon(muon_entrypoint, muon_direction, skroot_lowe_.bsvertex, muon_dedx);
    Log(m_unique_name+": dist to muon: "+toString(dist),v_debug,m_verbose);
    
    pre_ldt_cut.Fill(dist);
    if(!selectorName.empty()) m_data->ApplyCut(selectorName, "dist", dist);
    if (dist > 500){
      throw reject("neutron candidate does not meet dist to muon cut of < 500");
    }
    post_ldt_cut.Fill(abs(dist));
    
    Log(m_unique_name+": Passing Neutron!",v_debug,m_verbose);
    NeutronInfo n;
    n.bs_goodness = skroot_lowe_.bsgood[1];
    n.bs_dirks = skroot_lowe_.bsdirks;
    n.bsn50 = skroot_lowe_.bsn50;
    n.bs_vertex = {skroot_lowe_.bsvertex[0], skroot_lowe_.bsvertex[1], skroot_lowe_.bsvertex[2]};
    neutrons.push_back(n);
      
  } catch(reject& r){
    Log(m_unique_name+r.what(),v_message,m_verbose);
    --N;
  }
  
  if(neutrons.size() == N){
    m_data->CStore.Set("event_neutrons", neutrons);
    //std::cout<<"fill h_n_neutrons with "<<N<<std::endl;
    h_n_neutrons.Fill(N);
    neutrons.clear();
    N=0;
  }
  
  return true;
}


bool PostReconstructionNeutronCloudSelection::Finalise(){
  
  outfile->cd();
  
  bsenergy_plot.Write();
  h_n_neutrons.Write();
  
  pre_bsgood_cut.Write();
  pre_bsdirks_cut.Write();
  pre_bsn50_cut.Write();
  pre_ldt_cut.Write();
  pre_ldt_cut2.Write();

  post_bsgood_cut.Write();
  post_bsdirks_cut.Write();
  post_bsn50_cut.Write();
  post_ldt_cut.Write();
  
  outfile->Close();
  
  return true;
}

float PostReconstructionNeutronCloudSelection::CalculateDistanceToMuon(float* muon_entrypoint,
                                                                       float* muon_direction,
                                                                       float* bs_vertex,
                                                                       float* muon_dedx) {
/* cross product of muon direction and vector from muon to neutron = transverse distance...
  float result = 0;
  for (int i = 0; i < 3; ++i){
    result += pow( ( muon_entry_point[i] - bs_vertex[i] ) * muon_direction[(i+1)%3] -
                 ( muon_entry_point[(i+1)%3] - bs_vertex[(i+1)%3] ) * muon_direction[i], 2);
  }
  return sqrt(result);
*/
  
  /*
  std::cout<<"muon enters at ("<<muon_entrypoint[0]<<", "<<muon_entrypoint[1]<<", "<<muon_entrypoint[2]<<") "
           <<"with direction ("<<muon_direction[0]<<", "<<muon_direction[1]<<", "<<muon_direction[2]<<"), "
           <<"neutron candidate is at ("<<bs_vertex[0]<<", "<<bs_vertex[1]<<", "<<bs_vertex[2]<<")"<<std::endl;
  */
  
  // FIXME tech notes are vague on how it defines 'distance to muon' - do we just take dlt?
  // what if the muon stops before the point of closest approach?
  // spallation calcultion determines both dlt and dll, where latter is calculated from
  // point of maximum energy deposition along the muon track. Do we do the same here?
  
  // FIXME would be better not to duplicate this functionality already in RelicMuonPlots/CalculateSpallationVariables
  // get transverse distance and distance along muon track of point of closest approach
  float dlt = 0, appr = 0;
  getdl_(muon_direction,
         &bs_vertex[0],
         &bs_vertex[1],
         &bs_vertex[2],
         muon_entrypoint,
         &dlt,
         &appr);
  Log(m_unique_name+": getdl suggests transvese distance of "+toString(dlt)+" cm with point of closest approach at "
           +toString(appr)+" cm along muon track",v_debug,m_verbose);
  pre_ldt_cut2.Fill(dlt);
  
  // perhaps transverse distance is a better metric than total distance...?
  // but we don't necessarily care about getting 'all' the neutrons - we just want the cloud vertex,
  // so perhaps a smaller sample of 'better' neutrons gives a more accurate vertex...?
  //return dlt;
  
  // find point of maximum energy deposition along muon track
  double max_edep = 0;
  int max_edep_bin=0;
  for(int i=0;i<111;i++){
    double e_dep_in_window = 0 ;
    for(int j=0;j<9;j++){
      e_dep_in_window = e_dep_in_window + muon_dedx[i+j];
    }
    if(e_dep_in_window > max_edep){
      max_edep_bin = i+4;
      max_edep = e_dep_in_window;
    }
  }
  double max_energy_dep_pos = 50.*max_edep_bin;
  float dll = max_energy_dep_pos - appr;
  Log(m_unique_name+": max e dep at "+toString(max_energy_dep_pos)+" cm along track, or "+toString(dll)
     +" cm from point of closest approach",v_debug,m_verbose);
  
  return sqrt(pow(dlt,2.)+pow(dll,2.));

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
