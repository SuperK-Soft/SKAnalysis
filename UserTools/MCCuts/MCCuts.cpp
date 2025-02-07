#include "MCCuts.h"

#include "TH1D.h"
#include "TROOT.h"

#include "fortran_routines.h"

MCCuts::MCCuts():Tool(){}

bool MCCuts::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();
  
  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReader();

  std::string type_str = "";
  m_variables.Get("type", type_str);
  if (type_str == "IBD"){
    event_type = IBD;
  } else if (type_str == "ATMOS"){
    event_type = ATMOS;
  } else {
    throw std::runtime_error("MCCuts::Initialise() valid values for type are IBD and ATMOS");
  }

  SetUpPlots();
    
  return true;
}

bool MCCuts::Execute(){
  
  GetValues();
  
  if (lowe_ptr->bsenergy > 9000){
    return true;
  }

  xy.Fill(lowe_ptr->bsvertex[0], lowe_ptr->bsvertex[1]);
  rz.Fill(sqrt(pow(lowe_ptr->bsvertex[0],2) + pow(lowe_ptr->bsvertex[1],2)), lowe_ptr->bsvertex[2]);
 
  double dwall = wallsk_(lowe_ptr->bsvertex);
  std::cout << "dwall: " << dwall << std::endl;
  double ovaq = pow(lowe_ptr->bsgood[1],2) - pow(lowe_ptr->bsdirks, 2);
  std::cout << "ovaq: " << ovaq << std::endl;


  // count number of tagged neutrons
  int n_neutrons = 0;

  double likelihood_cut = 0;
  bool ok = m_variables.Get("likelihood_cut", likelihood_cut);
  if (!ok){throw std::runtime_error("MCCuts::Execute: Couldn't retrieve likelihood_cut!");}
  
  for (const auto& likelihood : *neutron_likelihoods_ptr){
      if (likelihood >= likelihood_cut){++n_neutrons;}
  }

  // cut efficiencies:
  if (event_type == IBD){
    if (n_neutrons == 1){

      for (const auto& likelihood : *neutron_likelihoods_ptr){
	likelihood_plot_1n.Fill(likelihood);
      }

      nqisk_plot_1n.Fill(skq_.nqisk);
      ovaq_plot_1n.Fill(ovaq);
      dwall_plot_1n.Fill(dwall);
      bsgood_plot_1n.Fill(lowe_ptr->bsgood[1]);
      clgood_plot_1n.Fill(lowe_ptr->clgoodness);

    } else if (n_neutrons == 0) {
    
      for (const auto& likelihood : *neutron_likelihoods_ptr){
	likelihood_plot_0n.Fill(likelihood);
      }
  
      nqisk_plot_0n.Fill(skq_.nqisk);
      ovaq_plot_0n.Fill(ovaq);
      dwall_plot_0n.Fill(dwall);
      bsgood_plot_0n.Fill(lowe_ptr->bsgood[1]);
      clgood_plot_0n.Fill(lowe_ptr->clgoodness);
    
    } else {

      for (const auto& likelihood : *neutron_likelihoods_ptr){
	likelihood_plot_n1o0n.Fill(likelihood);
      }
  
      nqisk_plot_n1o0n.Fill(skq_.nqisk);
      ovaq_plot_n1o0n.Fill(ovaq);
      dwall_plot_n1o0n.Fill(dwall);
      bsgood_plot_n1o0n.Fill(lowe_ptr->bsgood[1]);
      clgood_plot_n1o0n.Fill(lowe_ptr->clgoodness);
    
    }
  } else if (event_type == ATMOS){

    int n_true_neutrons = 0;
    m_data->CStore.Get("n_true_neutrons", n_true_neutrons);
    std::cout << "number of true neutrons in the event: " << n_true_neutrons << std::endl;
    std::cout << "number of tagged neutrons in the event: " << n_neutrons << std::endl;

    true_to_tagged.Fill(n_true_neutrons, n_neutrons); //x, y

    atmos_nqisk_neutron_plots.at(indexer(n_true_neutrons, n_neutrons))->Fill(skq_.nqisk);
    atmos_dwall_neutron_plots.at(indexer(n_true_neutrons, n_neutrons))->Fill(dwall);
    atmos_bsgood_neutron_plots.at(indexer(n_true_neutrons, n_neutrons))->Fill(lowe_ptr->bsgood[1]);
    atmos_clgood_neutron_plots.at(indexer(n_true_neutrons, n_neutrons))->Fill(lowe_ptr->clgoodness);
    atmos_ovaq_neutron_plots.at(indexer(n_true_neutrons, n_neutrons))->Fill(ovaq);
    
    for (const auto& l : *neutron_likelihoods_ptr){
      atmos_ntag_neutron_plots.at(indexer(n_true_neutrons, n_neutrons))->Fill(l);
    }

    
      // for (const auto& likelihood : *neutron_likelihoods_ptr){likelihood_plot_0tr_0ta.Fill(likelihood);}
      // nqisk_plot_0tr_0ta.Fill(skq_.nqisk);
      // ovaq_plot_0tr_0ta.Fill(ovaq);
      // dwall_plot_0tr_0ta.Fill(dwall);
      // bsgood_plot_0tr_0ta.Fill(lowe_ptr->bsgood[1]);
      // clgood_plot_0tr_0ta.Fill(lowe_ptr->clgoodness);
    }
      
  // actually do the cutting:
      if (event_type == ATMOS){
    if (interaction_mode ==  51 || interaction_mode ==  52 || 
	interaction_mode == -51 || interaction_mode == -52){
      // NCQE
      NCQE_nocut_plot.Fill(lowe_ptr->bsenergy);
      if (skq_.nqisk > 2000){
	return true;
      }
      NCQE_nqisk_plot.Fill(lowe_ptr->bsenergy);
      if (dwall < 200){
	return true;
      }
      NCQE_nqisk_dwall_plot.Fill(lowe_ptr->bsenergy);
      if (n_neutrons != 1){
	return true;
      }
      NCQE_nqisk_dwall_ntag_plot.Fill(lowe_ptr->bsenergy);
      if (lowe_ptr->bsgood[1] <= 0.5){
	return true;
      }
      NCQE_nqisk_dwall_ntag_bsgood_plot.Fill(lowe_ptr->bsenergy);
      if (lowe_ptr->clgoodness < 0.3){
	return true;
      }
      NCQE_nqisk_dwall_ntag_bsgood_clgood_plot.Fill(lowe_ptr->bsenergy);
      if (ovaq < 0.25){
	return true;
      }
      NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot.Fill(lowe_ptr->bsenergy);
      if (lowe_ptr->bsenergy < 8 || lowe_ptr->bsenergy > 100){
	return true;
      }
      NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot.Fill(lowe_ptr->bsenergy);
    
    } else {
      // non-NCQE
      nonNCQE_nocut_plot.Fill(lowe_ptr->bsenergy);
      if (skq_.nqisk > 2000){
	return true;
      }
      nonNCQE_nqisk_plot.Fill(lowe_ptr->bsenergy);
      if (dwall < 200){
	return true;
      }
      nonNCQE_nqisk_dwall_plot.Fill(lowe_ptr->bsenergy);
      if (n_neutrons != 1){
	return true;
      }
      nonNCQE_nqisk_dwall_ntag_plot.Fill(lowe_ptr->bsenergy);
      if (lowe_ptr->bsgood[1] <= 0.5){
	return true;
      }
      nonNCQE_nqisk_dwall_ntag_bsgood_plot.Fill(lowe_ptr->bsenergy);
      if (lowe_ptr->clgoodness < 0.3){
	return true;
      }
      nonNCQE_nqisk_dwall_ntag_bsgood_clgood_plot.Fill(lowe_ptr->bsenergy);
      if (ovaq < 0.25){
	return true;
      }
      nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot.Fill(lowe_ptr->bsenergy);
      if (lowe_ptr->bsenergy < 8 || lowe_ptr->bsenergy > 100){
	return true;
      }
      nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot.Fill(lowe_ptr->bsenergy);
    }
  } else {
    //IBD
    for (int i = 0; i < weight_names.size(); ++i){
      IBD_nocut_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));

      if (skq_.nqisk > 2000){continue;}
      IBD_nqisk_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));

      if (dwall < 200){continue;}
      IBD_nqisk_dwall_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));

      if (n_neutrons != 1){continue;}
      IBD_nqisk_dwall_ntag_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));

      if (lowe_ptr->bsgood[1] <= 0.5){continue;}
      IBD_nqisk_dwall_ntag_bsgood_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));

      if (lowe_ptr->clgoodness < 0.3){continue;}
      IBD_nqisk_dwall_ntag_bsgood_clgood_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));

      if (ovaq < 0.25){continue;}
      IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));

      if (lowe_ptr->bsenergy < 8 || lowe_ptr->bsenergy > 100){continue;}
      IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plots.at(i)->Fill(lowe_ptr->bsenergy, weight_values.at(i));
    }
  }

  // double cherenkov_angle = acos(skroot_lowe_.bsresult[3]) * 180 / 3.1415926535;
  // m_data->CStore.Set("c_angle", cherenkov_angle);
  
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

  NCQE_cut_efficiencies.Fill("nocut", NCQE_nocut_plot.GetEntries());
  NCQE_cut_efficiencies.Fill("nqisk", NCQE_nqisk_plot.GetEntries());
  NCQE_cut_efficiencies.Fill("dwall", NCQE_nqisk_dwall_plot.GetEntries());
  NCQE_cut_efficiencies.Fill("ntag", NCQE_nqisk_dwall_ntag_plot.GetEntries());
  NCQE_cut_efficiencies.Fill("bsgood", NCQE_nqisk_dwall_ntag_bsgood_plot.GetEntries());
  NCQE_cut_efficiencies.Fill("clgood", NCQE_nqisk_dwall_ntag_bsgood_clgood_plot.GetEntries());
  NCQE_cut_efficiencies.Fill("ovaq", NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot.GetEntries());
  NCQE_cut_efficiencies.Fill("bsenergy", NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot.GetEntries());

  non_NCQE_cut_efficiencies.Fill("nocut", nonNCQE_nocut_plot.GetEntries());
  non_NCQE_cut_efficiencies.Fill("nqisk", nonNCQE_nqisk_plot.GetEntries());
  non_NCQE_cut_efficiencies.Fill("dwall", nonNCQE_nqisk_dwall_plot.GetEntries());
  non_NCQE_cut_efficiencies.Fill("ntag", nonNCQE_nqisk_dwall_ntag_plot.GetEntries());
  non_NCQE_cut_efficiencies.Fill("bsgood", nonNCQE_nqisk_dwall_ntag_bsgood_plot.GetEntries());
  non_NCQE_cut_efficiencies.Fill("clgood", nonNCQE_nqisk_dwall_ntag_bsgood_clgood_plot.GetEntries());
  non_NCQE_cut_efficiencies.Fill("ovaq", nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot.GetEntries());
  non_NCQE_cut_efficiencies.Fill("bsenergy", nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot.GetEntries());

  IBD_cut_efficiencies.Fill("nocut", IBD_nocut_plots.at(0)->GetEntries());
  IBD_cut_efficiencies.Fill("nqisk", IBD_nqisk_plots.at(0)->GetEntries());
  IBD_cut_efficiencies.Fill("dwall", IBD_nqisk_dwall_plots.at(0)->GetEntries());
  IBD_cut_efficiencies.Fill("ntag", IBD_nqisk_dwall_ntag_plots.at(0)->GetEntries());
  IBD_cut_efficiencies.Fill("bsgood", IBD_nqisk_dwall_ntag_bsgood_plots.at(0)->GetEntries());
  IBD_cut_efficiencies.Fill("clgood", IBD_nqisk_dwall_ntag_bsgood_clgood_plots.at(0)->GetEntries());
  IBD_cut_efficiencies.Fill("ovaq", IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_plots.at(0)->GetEntries());
  IBD_cut_efficiencies.Fill("bsenergy", IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plots.at(0)->GetEntries());
  
  xy.Write();
  rz.Write();
  true_to_tagged.Write();
  
  likelihood_plot_1n.Write();
  bsgood_plot_1n.Write();
  clgood_plot_1n.Write();
  ovaq_plot_1n.Write();
  dwall_plot_1n.Write();
  nqisk_plot_1n.Write();

  likelihood_plot_0n.Write();
  bsgood_plot_0n.Write();
  clgood_plot_0n.Write();
  ovaq_plot_0n.Write();
  dwall_plot_0n.Write();
  nqisk_plot_0n.Write();

  likelihood_plot_n1o0n.Write();
  bsgood_plot_n1o0n.Write();
  clgood_plot_n1o0n.Write();
  ovaq_plot_n1o0n.Write();
  dwall_plot_n1o0n.Write();
  nqisk_plot_n1o0n.Write();
    
  if (event_type == ATMOS){
    NCQE_nocut_plot.Write();
    NCQE_nqisk_plot.Write();
    NCQE_nqisk_dwall_plot.Write();
    NCQE_nqisk_dwall_ntag_plot.Write();
    NCQE_nqisk_dwall_ntag_bsgood_plot.Write();
    NCQE_nqisk_dwall_ntag_bsgood_clgood_plot.Write();
    NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot.Write();
    NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot.Write();
    NCQE_cut_efficiencies.Write();
    
    nonNCQE_nocut_plot.Write();
    nonNCQE_nqisk_plot.Write();
    nonNCQE_nqisk_dwall_plot.Write();
    nonNCQE_nqisk_dwall_ntag_plot.Write();
    nonNCQE_nqisk_dwall_ntag_bsgood_plot.Write();
    nonNCQE_nqisk_dwall_ntag_bsgood_clgood_plot.Write();
    nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot.Write();
    nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot.Write();
    non_NCQE_cut_efficiencies.Write();

    for (int i = 0; i < 9; ++i){
      atmos_nqisk_neutron_plots.at(i)->Write();
      atmos_dwall_neutron_plots.at(i)->Write();
      atmos_ntag_neutron_plots.at(i)->Write();
      atmos_clgood_neutron_plots.at(i)->Write();
      atmos_bsgood_neutron_plots.at(i)->Write();
      atmos_ovaq_neutron_plots.at(i)->Write();

    }
    
  } else if (event_type == IBD){
    for (int i = 0; i < weight_names.size(); ++i){
    IBD_nocut_plots.at(i)->Write();
    IBD_nqisk_plots.at(i)->Write();
    IBD_nqisk_dwall_plots.at(i)->Write();
    IBD_nqisk_dwall_ntag_plots.at(i)->Write();
    IBD_nqisk_dwall_ntag_bsgood_plots.at(i)->Write();
    IBD_nqisk_dwall_ntag_bsgood_clgood_plots.at(i)->Write();
    IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_plots.at(i)->Write();
    IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plots.at(i)->Write();
    }
    IBD_cut_efficiencies.Write();
    
  }

  outfile->Write();
   
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

void MCCuts::GetValues(){

  bool ok = tree_reader_ptr->Get("neutron_likelihoods", neutron_likelihoods_ptr);
  if (!ok){throw std::runtime_error("couldn't get neutron_likelihoods branch");}
  
  ok = tree_reader_ptr->Get("LOWE", lowe_ptr);
  if (!ok){throw std::runtime_error("couldn't get LOWE branch");}

  if (event_type == ATMOS){
    ok = tree_reader_ptr->Get("interaction_mode", interaction_mode);
    if (!ok){throw std::runtime_error("couldn't get interaction_mode branch");};
    return;
  } else {
    for (int i = 0; i < weight_names.size(); ++i){
      ok = tree_reader_ptr->Get(weight_names.at(i), weight_values.at(i));
      if (!ok){throw std::runtime_error("couldn't retrieve "+weight_names.at(i)+" value");};
    }
  }
    return;
}

void MCCuts::SetUpPlots(){

  xy = TH2D("xy","xy;x;y", 10,0,0,10,0,0);
  rz = TH2D("rz","rz;r;z", 10,0,0,10,0,0);
  
  likelihood_plot_1n = TH1D("likelihood_plot_1n", "bdt likelihoods;likelihood_metric;numb. candidates", 100, 0, 0);
  bsgood_plot_1n = TH1D("bsgood_plot_1n", "bsgood;bsgood;n events", 100, 0, 0);
  clgood_plot_1n = TH1D("clgood_plot_1n", "clgood;clgood;n events", 100, 0, 0);
  ovaq_plot_1n = TH1D("ovaq_plot_1n", "ovaq;ovaq;n events", 100, 0, 0);
  dwall_plot_1n = TH1D("wall_plot_1n", "wall;wall; nevents", 100, 0, 0);
  nqisk_plot_1n = TH1D("nqisk_plot_1n", "nqisk;nqisk; nevents", 100, 0, 0);

  likelihood_plot_0n = TH1D("likelihood_plot_0n", "bdt likelihoods;likelihood_metric;numb. candidates", 100, 0, 0);
  bsgood_plot_0n = TH1D("bsgood_plot_0n", "bsgood;bsgood;n events", 100, 0, 0);
  clgood_plot_0n = TH1D("clgood_plot_0n", "clgood;clgood;n events", 100, 0, 0);
  ovaq_plot_0n = TH1D("ovaq_plot_0n", "ovaq;ovaq;n events", 100, 0, 0);
  dwall_plot_0n = TH1D("wall_plot_0n", "wall;wall; nevents", 100, 0, 0);
  nqisk_plot_0n = TH1D("nqisk_plot_0n", "nqisk;nqisk; nevents", 100, 0, 0);
  
  likelihood_plot_n1o0n = TH1D("likelihood_plot_n1o0n", "bdt likelihoods;likelihood_metric;numb. candidates", 100, 0, 0);
  bsgood_plot_n1o0n = TH1D("bsgood_plot_n1o0n", "bsgood;bsgood;n events", 100, 0, 0);
  clgood_plot_n1o0n = TH1D("clgood_plot_n1o0n", "clgood;clgood;n events", 100, 0, 0);
  ovaq_plot_n1o0n = TH1D("ovaq_plot_n1o0n", "ovaq;ovaq;n events", 100, 0, 0);
  dwall_plot_n1o0n = TH1D("wall_plot_n1o0n", "wall;wall; nevents", 100, 0, 0);
  nqisk_plot_n1o0n = TH1D("nqisk_plot_n1o0n", "nqisk;nqisk; nevents", 100, 0, 0);
  
  true_to_tagged = TH2D("true_to_tagged_plot", "true_to_tagged_plot", 10, 0, 10, 10, 0, 10); 
  
  NCQE_nocut_plot = TH1D("nocut", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  NCQE_nqisk_plot = TH1D("nqisk", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  NCQE_nqisk_dwall_plot = TH1D("dwall", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  NCQE_nqisk_dwall_ntag_plot = TH1D("ntag", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  NCQE_nqisk_dwall_ntag_bsgood_plot = TH1D("bsgood", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  NCQE_nqisk_dwall_ntag_bsgood_clgood_plot = TH1D("clgood", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot = TH1D("ovaq", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot = TH1D("bsenergy", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  
  nonNCQE_nocut_plot = TH1D("nonNCQE_nocut", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  nonNCQE_nqisk_plot = TH1D("nonNCQE_nqisk", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  nonNCQE_nqisk_dwall_plot = TH1D("nonNCQE_dwall", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  nonNCQE_nqisk_dwall_ntag_plot = TH1D("nonNCQE_ntag", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  nonNCQE_nqisk_dwall_ntag_bsgood_plot = TH1D("nonNCQE_bsgood", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  nonNCQE_nqisk_dwall_ntag_bsgood_clgood_plot = TH1D("nonNCQE_clgood", ";bsenergy [MeV];Events/MeV", 100, 0, 100);
  nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot = TH1D("nonNCQE_ovaq", ";bsenergy [MeV];Events/MeV", 100, 0, 100);  
  nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot = TH1D("nonNCQE_bsenergy", ";bsenergy [MeV];Events/MeV", 100, 0, 100);

  for (int i = 0; i < weight_names.size(); ++i){
    IBD_nocut_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_nocut").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
    IBD_nqisk_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_nqisk").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
    IBD_nqisk_dwall_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_dwall").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
    IBD_nqisk_dwall_ntag_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_ntag").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
    IBD_nqisk_dwall_ntag_bsgood_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_bsgood").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
    IBD_nqisk_dwall_ntag_bsgood_clgood_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_clgood").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
    IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_ovaq").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
    IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plots.emplace_back( new TH1D((weight_names.at(i)+"_IBD_bsenergy").c_str(), ";bsenergy [MeV];Events/MeV", 100, 0, 100));
  }  
  NCQE_cut_efficiencies = TH1D("NCQE_cut_efficiencies", "NCQE_cut_efficiencies", 8, 0, 7);
  non_NCQE_cut_efficiencies = TH1D("non_NCQE_cut_efficiencies", "non_NCQE_cut_efficiencies", 8, 0, 7);
  IBD_cut_efficiencies = TH1D("IBD_cut_efficiencies", "IBD_cut_efficiencies", 8, 0, 7);

  for (int i = 0; i < 9; ++i){
    atmos_nqisk_neutron_plots.emplace_back( new TH1D(("atmos_nqisk_"+std::to_string(i)).c_str(), ";nqisk; n_events", 100, 0, 0));
    atmos_dwall_neutron_plots.emplace_back( new TH1D(("atmos_dwall_"+std::to_string(i)).c_str(), ";dwall; n_events", 100, 0, 0));
    atmos_ntag_neutron_plots.emplace_back( new TH1D(("atmos_ntag_"+std::to_string(i)).c_str(), ";ntag; n_events", 100, 0, 0));
    atmos_bsgood_neutron_plots.emplace_back( new TH1D(("atmos_bsgood_"+std::to_string(i)).c_str(), ";bsgood; n_events", 100, 0, 0));
    atmos_clgood_neutron_plots.emplace_back( new TH1D(("atmos_clgood_"+std::to_string(i)).c_str(), ";clgood; n_events", 100, 0, 0));
    atmos_ovaq_neutron_plots.emplace_back( new TH1D(("atmos_ovaq_"+std::to_string(i)).c_str(), ";ovaq; n_events", 100, 0, 0));
  }
  
  return;  
}
