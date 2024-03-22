#include "MCCuts.h"

#include "TH1D.h"
#include "TROOT.h"

MCCuts::MCCuts():Tool(){}

bool MCCuts::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();
  
  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReader();

  LUN = m_data->GetLUN(tree_reader_str);
  TreeManager* tree_manager = skroot_get_mgr(&LUN);
  if (tree_manager == nullptr){
    throw std::runtime_error("MCCuts::Initialise: Couldn't get tree manager!");
  }
  TTree* out_tree = tree_manager->GetOTree();
  if (out_tree == nullptr){
    throw std::runtime_error("MCCuts::Initialise: Null output tree");
  }

  out_tree->Branch("passed_mc_cuts", &pass);

  gROOT->cd();
    
  nqisk_plot = TH1D("nqisk_plot", "nqisk", 100, 20, 20);
  bsenergy_plot = TH1D("bsenergy_plot", "bsenergy", 100, 20, 20);
  bsgoodness_plot = TH1D("bsgoodness_plot", "bsgoodness", 100, 20, 20);
  ovaq_plot = TH1D("ovaq_plot", "ovaq", 100, 20, 20);
  clusfit_goodness_plot = TH1D("clusfit_goodness_plot", "clusfit_goodness", 100, 20, 20);
  d_wall_plot = TH1D("d_wall_plot", "d_wall", 100, 20, 20);
  n_od_plot = TH1D("n_od_plot", "n_od", 100, 20, 20);
  
  return true;
}


bool MCCuts::Execute(){
  
  skroot_get_lowe_(&skheadf_.root_id,
		   &get_ok,
		   skroot_lowe_.bsvertex,
		   skroot_lowe_.bsresult,
		   skroot_lowe_.bsdir,
		   skroot_lowe_.bsgood,
		   &skroot_lowe_.bsdirks,
		   skroot_lowe_.bseffhit,
		   &skroot_lowe_.bsenergy,
		   &skroot_lowe_.bsn50,
		   &skroot_lowe_.bscossun,
		   skroot_lowe_.clvertex,
		   skroot_lowe_.clresult,
		   skroot_lowe_.cldir,
		   &skroot_lowe_.clgoodness,
		   &skroot_lowe_.cldirks,
		   skroot_lowe_.cleffhit,
		   &skroot_lowe_.clenergy,
		   &skroot_lowe_.cln50,
		   &skroot_lowe_.clcossun,
		   &skroot_lowe_.latmnum,
		   &skroot_lowe_.latmh,
		   &skroot_lowe_.lmx24,
		   &skroot_lowe_.ltimediff,
		   &skroot_lowe_.lnsratio,
		   skroot_lowe_.lsdir,
		   &skroot_lowe_.spaevnum,
		   &skroot_lowe_.spaloglike,
		   &skroot_lowe_.sparesq,
		   &skroot_lowe_.spadt,
		   &skroot_lowe_.spadll,
		   &skroot_lowe_.spadlt,
		   &skroot_lowe_.spamuyn,
		   &skroot_lowe_.spamugdn,
		   skroot_lowe_.posmc,
		   skroot_lowe_.dirmc,
		   skroot_lowe_.pabsmc,
		   skroot_lowe_.energymc,
		   &skroot_lowe_.darkmc,
		   &skroot_lowe_.islekeep,
		   &skroot_lowe_.bspatlik,
		   &skroot_lowe_.clpatlik,
		   &skroot_lowe_.lwatert,
		   &skroot_lowe_.lninfo,
		   skroot_lowe_.linfo);

  double cherenkov_angle = acos(skroot_lowe_.bsresult[3]) * 180 / 3.1415926535;
  m_data->CStore.Set("c_angle", cherenkov_angle);
  
  const double nqisk = skq_.nqisk;
  nqisk_plot.Fill(nqisk);

  double bsenergy = skroot_lowe_.bsenergy;
  if (bsenergy < 1000){
    bsenergy_plot.Fill(bsenergy);
  }
    m_data->CStore.Set("bsenergy", bsenergy);
  

  double positron_energy=-1;
  positron_energy = skroot_lowe_.energymc[1];
  // for(int i=0; i<m_data->eventParticles.size(); ++i){
  //   MParticle& aparticle = m_data->eventParticles.at(i);
  //   if(aparticle.pdg==-11){
  //     if(aparticle.GetStartE()!=nullptr) positron_energy = *aparticle.GetStartE();
  //     break;
  //   }
  // }
  m_data->CStore.Set("mcenergy", positron_energy);
  
  const double bsgoodness = skroot_lowe_.bsgood[1];
  if (bsgoodness < 100){bsgoodness_plot.Fill(bsgoodness);}

  const double ovaq = skroot_lowe_bsovaq;
  ovaq_plot.Fill(ovaq);

  const double clusfit_goodness = skroot_lowe_.clgoodness;
  clusfit_goodness_plot.Fill(clusfit_goodness);

  const double wall = skroot_lowe_bswallsk;
  d_wall_plot.Fill(wall);
  
  lfnhita_(&skroot_lowe_lnahit);
  n_od_plot.Fill(skroot_lowe_lnahit);
  
  pass = ((nqisk < 2000) &&
		     (bsenergy > 8 && bsenergy < 100) &&
		     (bsgoodness > 0.5) &&
		     (ovaq > 0.25) &&
		     (clusfit_goodness > 0.3) &&
		     (wall > 200) &&
		     (skroot_lowe_lnahit < 20));

  m_data->CStore.Set("pass", pass);
  
  skroot_fill_tree_(&LUN);
  
  return true;
}


bool MCCuts::Finalise(){

  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "MC_cuts_plots_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()){
    throw std::runtime_error("MCCuts::Finalise - Couldn't open output file");
  }

  outfile->cd();
  
  nqisk_plot.Write();
  bsenergy_plot.Write();
  bsgoodness_plot.Write();
  ovaq_plot.Write();
  clusfit_goodness_plot.Write();
  d_wall_plot.Write();
  n_od_plot.Write();
  
  return true;
}

void MCCuts::GetReader(){
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
