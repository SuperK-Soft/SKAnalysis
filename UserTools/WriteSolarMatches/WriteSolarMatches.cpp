#include "WriteSolarMatches.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TFile.h"

WriteSolarMatches::WriteSolarMatches():Tool(){}

bool WriteSolarMatches::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  // SK-IV: 490cm, SK-VI: 400cm
  m_variables.Get("match_dist_limit",match_dist_limit);
  
  // get solar tree
  GetReader();
  // add branch to hold match distances, which are calculated in this Tool
  if(solarTree) solarTree->Branch("MatchedDistance",&matched_dists);
  
  // get the relics from this run. These are loaded into memory
  // by the SolarPreSelection Tool, which populates time match info
  intptr_t solar_relics_ptr=0;
  bool ok = m_data->CStore.Get("solar_relics",solar_relics_ptr);
  if(!ok || solar_relics_ptr==0){
    throw std::runtime_error("WriteSolarMatches::Initialise - couldn't get solar relics. Is SolarPreSelection in ToolChain?");
  }
  relics_this_run=reinterpret_cast<std::deque<SolarRelic>*>(solar_relics_ptr);
  
  return true;
}


bool WriteSolarMatches::Execute(){

  m_data->CStore.Get("solar_relic_matches",matched_relics); // indices in relics_this_run deque
  m_data->CStore.Get("solar_relic_tdiffs",matched_tdiffs);
  
  matched_dists.clear();
  
  // record this as a matched solar for all the relics it was time matched to
  for(int& relic_num : matched_relics){
    
    SolarRelic& a_relic = relics_this_run->at(relic_num);
    
    // calculate distance to this relic
    double dist = sqrt(pow(a_relic.vtx[0]-skroot_lowe_.bsvertex[0],2)+
                       pow(a_relic.vtx[1]-skroot_lowe_.bsvertex[1],2)+
                       pow(a_relic.vtx[2]-skroot_lowe_.bsvertex[2],2)
                      );
    
    // cut tracking
    if(!solarSelectorName.empty()){
            m_data->ApplyCut(solarSelectorName, "relic_solar_dist", dist);
    }
    
    // add this match to the relic
    a_relic.matched_ev_nums.push_back(skhead_.nevsk);
    a_relic.matched_in_entry_nums.push_back(rfmReader->GetEntryNumber());
    outentry = (solarTree) ? solarTree->GetEntries() : 0;
    a_relic.matched_out_entry_nums.push_back(outentry); // no -1 as we've not yet added this solar to the tree
    a_relic.matched_tdiffs.push_back(-matched_tdiffs.at(relic_num));   // flip sign to make it relative to the relic
    a_relic.matched_dists.push_back(dist);
    if(dist<match_dist_limit){
      a_relic.rejected=true;
      a_relic.rejected_by.push_back(a_relic.matched_ev_nums.size()-1);
    }
    
  }
  
  if(solarTree){
    Log(m_unique_name+" Filling solar tree",v_message,m_verbose);
    solarTree->Fill();
  }
  
  return true;
}


bool WriteSolarMatches::Finalise(){
  
  if(solarTree){
    Log(m_unique_name+" writing out solar Tree",v_message,m_verbose);
    solarFile->cd();
    solarFile->Write("*",TObject::kOverwrite);
  }
  
  return true;
}

void WriteSolarMatches::GetReader(){
  std::string tree_reader_str = "";
  m_data->vars.Get("SolarRfmReader", tree_reader_str);
  if (tree_reader_str.empty() || m_data->Trees.count(tree_reader_str) == 0){
    throw std::invalid_argument("no valid treereader specified!");
  }
  rfmReader = m_data->Trees.at(tree_reader_str);
  if (rfmReader == nullptr){
    throw std::runtime_error("couldn't get treereader");
  }
  int lun = m_data->GetLUN(tree_reader_str);
  TreeManager* mgr = skroot_get_mgr(&lun);
  if(mgr->GetMode()==0){
    solarTree = mgr->GetOTree();
    solarFile = mgr->GetOTree()->GetCurrentFile();
  }
  
  return;
}
