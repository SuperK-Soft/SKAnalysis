#include "SplitAtmosInteractions.h"
#include "neworkC.h"

#include "loweroot.h"
#include "MTreeReader.h"

#include "Constants.h"

SplitAtmosInteractions::SplitAtmosInteractions():Tool(){}


bool SplitAtmosInteractions::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReaders();

  std::string outfile_str = "";
  m_variables.Get("outfile", outfile_str);
  if (outfile_str.empty()){throw std::runtime_error("SplitAtmosInteractions: no output file specified!");}

  output_file_ptr = new TFile(outfile_str.c_str(), "RECREATE");
  output_tree_ptr = new TTree("data", "data");

  // interaction_mode_plot = TH1D("interaction_mode", "interaction_mode", 100, 0, 0);
  // interaction_mode_zero = TH1D("interaction_mode_zero", "interaction_mode_zero;bsenergy", 100, 0, 0);
  // NCQE_plot = TH1D("NCQE_plot", "NCQE_plot;bsenergy", 100, 0, 0);
  // non_NCQE_plot = TH1D("non_NCQE_plot", "non_NCQE_plot;bsenergy", 100, 0, 0);
  
  return true;
}


bool SplitAtmosInteractions::Execute(){

  // LoweInfo* my_lowe_ptr = new LoweInfo();
  // tree_reader_root_ptr->Get("LOWE", my_lowe_ptr);
  // const double bsenergy = my_lowe_ptr->bsenergy;
  // std::cout << "reconstructed energy: " << bsenergy << std::endl;

  GetValues();
  if (first){
    SetOutputBranches(output_tree_ptr);
    first = false;
  }
  UpdateBranches();
  
  interaction_mode = nework_.modene;
  std::cout << "interaction mode: " << interaction_mode << std::endl;
  std::cout << NEUTInteractionModeToString(interaction_mode) << std::endl;
  std::cout << "neutron candidates: " << std::endl;
  std::cout << neutron_likelihoods_ptr->size() << std::endl;
  // if (bsenergy > 9000){
  //   return true;
  // }

  // interaction_mode_plot.Fill(interaction_mode);
  
  // if (interaction_mode == 0){
  //   interaction_mode_zero.Fill(bsenergy);
  // }
  // else if (interaction_mode ==  51 ||
  // 	   interaction_mode ==  52 ||
  // 	   interaction_mode == -51 ||
  // 	   interaction_mode == -52){
  //   NCQE_plot.Fill(bsenergy);
  // } else {
  //   non_NCQE_plot.Fill(bsenergy);
  // }

  // //interactions_for_pi_chart[NEUTInteractionModeToString(interaction_mode)]++;
  // interactions_for_pi_chart[GetNEUTModeProcess(interaction_mode)]++;

  output_tree_ptr->Fill();
  
  return true;
}


bool SplitAtmosInteractions::Finalise(){

  output_file_ptr->cd();
  output_file_ptr->Write();
  
  // TFile* output_file = TFile::Open("interaction_mode.root", "UPDATE");
  // if (output_file == nullptr){
  //   throw std::runtime_error("SplitAtmosInteractions::Finalise: Couldn't open output file!");
  // }
  
  // interaction_mode_zero.Write();
  // interaction_mode_plot.Write();
  // NCQE_plot.Write();
  // non_NCQE_plot.Write();

  // std::cout << "interaction values are\n";
  // for (const auto& [x,y] : interactions_for_pi_chart){
  //   std::cout << x << " = " << y << std::endl;
  // }
  
  return true;
}

void SplitAtmosInteractions::GetReaders(){
 std::string tree_reader_str = "";
  m_variables.Get("TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("SplitAtmosInteractions::Execute - Failed to get root treereader "+tree_reader_str+"!");
  }
  tree_reader_root_ptr = m_data->Trees.at(tree_reader_str);
  
  return;
} 

void SplitAtmosInteractions::GetValues(){

  bool ok = tree_reader_root_ptr->Get("LOWE", lowe_ptr);
  if (!ok){throw std::runtime_error("couldn't get LOWE branch");}

  ok = tree_reader_root_ptr->Get("HEADER", header_ptr);
  if (!ok){throw std::runtime_error("couldn't get HEADER branch");}
  
  ok = tree_reader_root_ptr->Get("MC", mcinfo_ptr);
  if (!ok){throw std::runtime_error("couldn't get MC branch");}
  
  ok = tree_reader_root_ptr->Get("TQREAL", tqreal_ptr);
  if (!ok){throw std::runtime_error("couldn't get TQREAL branch");}

  ok = tree_reader_root_ptr->Get("TQAREAL", tqareal_ptr);
  if (!ok){throw std::runtime_error("couldn't get TQAREAL branch");}
  
  ok = tree_reader_root_ptr->Get("neutron_likelihoods", neutron_likelihoods_ptr);
  if (!ok){throw std::runtime_error("couldn't get neutron_likelihoods branch");}

  ok = tree_reader_root_ptr->Get("SECONDARY", secondary_info_ptr);
  if (!ok){throw std::runtime_error("couldn't get SECONDARY branch!");}

  return ;
}

void SplitAtmosInteractions::UpdateBranches(){

  //output_tree_ptr->Print();
  //std::cout << "before UpdatesBranches" << std::endl;

  output_tree_ptr->SetBranchAddress("LOWE", &lowe_ptr);
  output_tree_ptr->SetBranchAddress("HEADER", &header_ptr);
  output_tree_ptr->SetBranchAddress("MC", &mcinfo_ptr);
  output_tree_ptr->SetBranchAddress("TQREAL", &tqreal_ptr);
  output_tree_ptr->SetBranchAddress("TQAREAL", &tqareal_ptr);
  output_tree_ptr->SetBranchAddress("neutron_likelihoods", &neutron_likelihoods_ptr);
  output_tree_ptr->SetBranchAddress("interaction_mode", &interaction_mode);
  output_tree_ptr->SetBranchAddress("SECONDARY", &secondary_info_ptr);
  
  return;
}

void SplitAtmosInteractions::SetOutputBranches(TTree* tree){

  tree->Branch("LOWE", lowe_ptr, 1024*1024, 0);
  tree->Branch("HEADER", header_ptr, 1024*1024, 0);
  tree->Branch("MC", mcinfo_ptr, 1024*1024, 0);
  tree->Branch("TQREAL", tqreal_ptr, 1024*1024, 0);
  tree->Branch("TQAREAL", tqareal_ptr, 1024*1024, 0);
  tree->Branch("neutron_likelihoods", neutron_likelihoods_ptr, 1024*1024, 0);
  tree->Branch("interaction_mode", &interaction_mode, 1024*1024, 0);
  tree->Branch("SECONDARY", "SecondaryInfo", secondary_info_ptr, 1024*1024, 0);

  
  std::cout << "after set output branches" << std::endl;
  //output_tree_ptr->Print();
  
  return;
}
