#ifndef MCCuts_H
#define MCCuts_H

#include <string>
#include <iostream>

#include "MTreeReader.h"

#include "Tool.h"

#include "TH1D.h"
#include "TH2D.h"


const auto indexer = [](int x, int y){
  if (x >= 2){x = 2;}
  if (y >= 2){y = 2;}
  return 3*x + y;
 };


class MCCuts: public Tool {

 public:

  MCCuts();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  enum EventType {IBD, ATMOS};

 private:

  void GetReader();
  void GetValues();
  void SetUpPlots();
  
  MTreeReader* tree_reader_ptr = nullptr;
  std::string tree_reader_str = "";
  int LUN = 0;
  bool pass = false;

  int interaction_mode = 0;
  LoweInfo* lowe_ptr = nullptr;
  std::vector<float>* neutron_likelihoods_ptr = nullptr;

  std::vector<std::string> weight_names = {
    "weight_hartmann",
    "weight_horiuchi1",
    "weight_horiuchi3_min",
    "weight_kaplinghat",
    "weight_horiuchi3",
    "weight_nakazato1",
    "weight_malaney",
    "weight_nakazato2",
    "weight_horiuchi2",
    "weight_totani",
    "weight_ando",
    "weight_lunardini",
    "weight_horiuchi3_center",
    "weight_li9",
    "weight_reactor"};
  
  TH1D likelihood_plot_1n;
  TH1D bsgood_plot_1n;
  TH1D clgood_plot_1n;
  TH1D ovaq_plot_1n;
  TH1D dwall_plot_1n;

  TH1D likelihood_plot_0n;
  TH1D bsgood_plot_0n;
  TH1D clgood_plot_0n;
  TH1D ovaq_plot_0n;
  TH1D dwall_plot_0n;

  TH1D likelihood_plot_n1o0n;
  TH1D bsgood_plot_n1o0n;
  TH1D clgood_plot_n1o0n;
  TH1D ovaq_plot_n1o0n;
  TH1D dwall_plot_n1o0n;
  
  TH1D NCQE_nocut_plot;
  TH1D NCQE_nqisk_plot;
  TH1D NCQE_nqisk_dwall_plot;
  TH1D NCQE_nqisk_dwall_ntag_plot;
  TH1D NCQE_nqisk_dwall_ntag_bsgood_plot;
  TH1D NCQE_nqisk_dwall_ntag_bsgood_clgood_plot;
  TH1D NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot;
  TH1D NCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot;

  TH1D nonNCQE_nocut_plot;
  TH1D nonNCQE_nqisk_plot;
  TH1D nonNCQE_nqisk_dwall_plot;
  TH1D nonNCQE_nqisk_dwall_ntag_plot;
  TH1D nonNCQE_nqisk_dwall_ntag_bsgood_plot;
  TH1D nonNCQE_nqisk_dwall_ntag_bsgood_clgood_plot;
  TH1D nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_plot;
  TH1D nonNCQE_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plot;


  // std::vector<TH1D> IBD_nocut_plots;
  // std::vector<TH1D> IBD_nqisk_plots;
  // std::vector<TH1D> IBD_nqisk_dwall_plots;
  // std::vector<TH1D> IBD_nqisk_dwall_ntag_plots;
  // std::vector<TH1D> IBD_nqisk_dwall_ntag_bsgood_plots;
  // std::vector<TH1D> IBD_nqisk_dwall_ntag_bsgood_clgood_plots;
  // std::vector<TH1D> IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_plots;
  // std::vector<TH1D> IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plots;

  std::vector<TH1D*> IBD_nocut_plots;
  std::vector<TH1D*> IBD_nqisk_plots;
  std::vector<TH1D*> IBD_nqisk_dwall_plots;
  std::vector<TH1D*> IBD_nqisk_dwall_ntag_plots;
  std::vector<TH1D*> IBD_nqisk_dwall_ntag_bsgood_plots;
  std::vector<TH1D*> IBD_nqisk_dwall_ntag_bsgood_clgood_plots;
  std::vector<TH1D*> IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_plots;
  std::vector<TH1D*> IBD_nqisk_dwall_ntag_bsgood_clgood_ovaq_bsenergy_plots;

  
  TH1D NCQE_cut_efficiencies;
  TH1D non_NCQE_cut_efficiencies;
  TH1D IBD_cut_efficiencies;

  TH2D xy;
  TH2D rz;
  TH2D true_to_tagged;
  
  TH1D nqisk_plot_1n;
  TH1D nqisk_plot_0n;
  TH1D nqisk_plot_n1o0n;
  TH1D bsenergy_plot;
  TH1D bsgoodness_plot;
  TH1D clusfit_goodness_plot;
  TH1D d_wall_plot;
  TH1D n_od_plot;

  double weight_value = 0;

  EventType event_type;

  std::vector<float> weight_values = std::vector<float>(weight_names.size(), 0);

  std::vector<TH1D*> atmos_nqisk_neutron_plots;
  std::vector<TH1D*> atmos_dwall_neutron_plots;
  std::vector<TH1D*> atmos_ntag_neutron_plots;
  std::vector<TH1D*> atmos_clgood_neutron_plots;
  std::vector<TH1D*> atmos_bsgood_neutron_plots;
  std::vector<TH1D*> atmos_ovaq_neutron_plots;
      
};

#endif
