#include "mufit_sk4.h"

#include "MTreeReader.h"

#include <bitset>

#include "skroot_loweC.h"
#include "sktqC.h"

#include "TH1D.h"

mufit_sk4::mufit_sk4():Tool(){}

bool mufit_sk4::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  GetReader();

  charge_plot = TH1D("r1", "r1", 100, 0, 2000);
  
  return true;
}

bool mufit_sk4::Execute(){
  
  //  if (tree_reader_ptr->GetEntryNumber() % 1000 == 0){
    std::cout << "nread: " << nread << ", nrunsk: " << skhead_.nrunsk << ", nevsk: " << skhead_.nevsk << ", nmuon: " << nmuon << "\n";
    //}

  // //ignore LED burst events
  // std::bitset<sizeof(int)*8> bits = skhead_.nevsk;
  // if (bits.test(26)){
  //   m_data->vars.Set("Skip", true);
  //   return true;
  // }

  //ignore incomplete events
  // if ( (bits & std::bitset(pow(2,20)) ) != 0){
  //   m_data->vars.Set("Skip", true);
  //   return true;    
  // }

  std::cout << "calling lfclear_all_()\n";
  charge_plot.Fill(skq_.qismsk);

  lfclear_all_();

  
  
  //skip muon events - judged by total number of PMT hits
  if(skq_.nqisk < 1000){
    m_data->vars.Set("Skip", true);
    std::cout << "skipping event\n";
    return true;
  } 

  ++nmuon;
  skroot_mu_.muinfo[0] = skq_.qismsk;
  skroot_mu_.muninfo = 1;
  
  fix_maxqisk_();

  skroot_mu_.muqismsk = skq_.qismsk;

  lfmufit_sk4_();

  int lun = m_data->GetLUN(tree_reader_str);
  skroot_set_mu_(&lun, skroot_mu_.muentpoint, skroot_mu_.mudir, &skroot_mu_.mutimediff,
		 &skroot_mu_.mugoodness, &skroot_mu_.muqismsk, &skroot_mu_.muyn, &skroot_mu_.mufast_flag,
		 &skroot_mu_.muboy_status, &skroot_mu_.muboy_ntrack, skroot_mu_.muboy_entpos, skroot_mu_.muboy_dir,
		 &skroot_mu_.muboy_goodness, &skroot_mu_.muboy_length, skroot_mu_.muboy_dedx,
		 skroot_mu_.mubff_entpos, skroot_mu_.mubff_dir, &skroot_mu_.mubff_goodness, &skroot_mu_.muninfo, skroot_mu_.muinfo);

  delete_outside_hits_();

  skroot_set_tree_(&lun);
  skroot_fill_tree_(&lun);
  
  return true;
}

bool mufit_sk4::Finalise(){

  charge_plot.SaveAs("charge_plot.root");
  
  return true;
}

void mufit_sk4::GetReader(){
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
