#include "SolarPostSelection.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TFile.h"

SolarPostSelection::SolarPostSelection():Tool(){}

bool SolarPostSelection::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  // TODO make cut thresholds input variables
  
  std::string plotfile_name = "";
  bool ok = m_variables.Get("plotfile_name", plotfile_name);
  if (!ok || plotfile_name.empty()){ plotfile_name = "solar_post_plots.root";}
  plotfile = TFile::Open(plotfile_name.c_str(), "RECREATE");
  if (plotfile == nullptr){
    throw std::runtime_error("SolarPostSelection::Initialise - Couldn't open plot file '"+plotfile_name+"'");
  }
  
  bsok_cut = TH1D("bsok_cut", "bsok_cut;bs_ok", 2, 0, 2); // did bonsai vertex and energy reconstruct correctly
  // remember these pre/post histograms should fix their binning to be able to overlay them
  // as well as hadd results from multiple jobs
  pre_bsenergy_cut = TH1D("pre_bsenergy_cut", "pre_bsenergy_cut;bsenergy", 100, 0, 100);
  pre_ovaQ_cut = TH1D("pre_ovaQ_cut","pre_ovaQ_cut:ovaQ",100,0,1);
  pre_FV_cut = TH1D("pre_FV_cut", "pre_FV_cut;wallsk", 100, 0, 1800);     // FV is wallsk>200
  pre_effwall_cut = TH1D("pre_effwall_cut", "pre_effwall_cut;d_effwall", 100, 0, 5000);
  
  post_bsenergy_cut = TH1D("post_bsenergy_cut", "post_bsenergy_cut;bsenergy", 100, 0, 100);
  post_ovaQ_cut = TH1D("post_ovaQ_cut","post_ovaQ_cut:ovaQ",100,0,1);
  post_FV_cut = TH1D("post_FV_cut", "post_FV_cut;wallsk", 100, 0, 1800);  // FV is wallsk>200
  post_effwall_cut = TH1D("post_effwall_cut", "post_effwall_cut;d_effwall", 100, 0, 5000);
  
  return true;
}


bool SolarPostSelection::Execute(){

  /* 
     Criteria on Reconstructed Solar Triggers are as follows:

     1. 6 < bonsai energy  < 25 MeV (energy successfully reconstructed)
     2. ovaQ > 0.2 (for bonsai energy > 7.5), ovaQ > 0.25 (for bonsai energy < 7.5)
     3. within FV == wallsk > 200 (vertex successfully reconstructed)
     4. deffwall > 400 (for bonsai energy > 8), deffwall > 650 (for bonsai energy < 8)
     
  */
  std::string rejectedby="";
  
  Log(m_unique_name+": this solar event has\nbsvertex: ("
     +toString(skroot_lowe_.bsvertex[0])+", "+toString(skroot_lowe_.bsvertex[1])+", "+toString(skroot_lowe_.bsvertex[2])+")"
     +"\nbsenergy: "+toString(skroot_lowe_.bsenergy)
     +"\novaQ: "+toString(skroot_lowe_bsovaq)  // alias for skroot_lowe_.linfo[25]
     +"\nwallsk: "+toString(skroot_lowe_bswallsk) // alias for skroot_lowe_.linfo[9]
     +"\ndeffwall: "+toString(skroot_lowe_bseffwal),v_debug,m_verbose); // alias for skroot_lowe_.info[5]
  
  // cut on bad vertex reconstruction
  if(skroot_lowe_.bsvertex[0]==9999 || skroot_lowe_.bsenergy == 9999){
    bsok_cut.Fill(0);
    rejectedby="candidate vertex or energy was not reconstructed properly";
    goto reject;
  }
  bsok_cut.Fill(1);
  
  // bsenergy cut
  pre_bsenergy_cut.Fill(skroot_lowe_.bsenergy);
  if (skroot_lowe_.bsenergy < 6 || skroot_lowe_.bsenergy > 25 ){
    rejectedby="candidate does not meet bsenergy cut of 6 < bse < 25";
    goto reject;
  }
  post_bsenergy_cut.Fill(skroot_lowe_.bsenergy);
  
  // ovaQ cut
  pre_ovaQ_cut.Fill(skroot_lowe_bsovaq);
  if (skroot_lowe_bsovaq < (skroot_lowe_.bsenergy > 7.5 ? 0.2 : 0.25)){
    rejectedby="candidate does not meet ovaQ cut of > 0.2 / 0.25";
    goto reject;
  }
  post_ovaQ_cut.Fill(skroot_lowe_bsovaq);
  
  // FV cut
  // superscan and lf_1st_reduction.F this is just dwall < 200cm
  pre_FV_cut.Fill(skroot_lowe_bswallsk);
  if (skroot_lowe_bswallsk < 200){
    rejectedby="candidate does not meet FV cut of wallsk > 200";
    goto reject;
  }
  post_FV_cut.Fill(skroot_lowe_bswallsk);
  
  // effwall cut
  pre_effwall_cut.Fill(skroot_lowe_bseffwal);
  if (skroot_lowe_bseffwal < (skroot_lowe_.bsenergy > 8 ? 400 : 650)){
    rejectedby="candidate does not meet effwall > 400 / 650";
    goto reject;
  }
  post_effwall_cut.Fill(skroot_lowe_bseffwal);
  
  Log(m_unique_name+": Passing solar event!",v_debug,m_verbose);
  return true;
  
  reject:
  Log(m_unique_name+"solar rejected: "+rejectedby,v_debug,m_verbose);
  m_data->vars.Set("Skip",1);
  return true;  
  
}


bool SolarPostSelection::Finalise(){
  
  plotfile->cd();
  
  bsok_cut.Write();
  pre_bsenergy_cut.Write();
  pre_ovaQ_cut.Write();
  pre_FV_cut.Write();
  pre_effwall_cut.Write();
  
  post_bsenergy_cut.Write();
  post_ovaQ_cut.Write();
  post_FV_cut.Write();
  post_effwall_cut.Write();
  
  plotfile->Close();
  
  return true;
}
