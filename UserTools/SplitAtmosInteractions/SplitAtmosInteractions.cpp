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

  GetReader();
  
  
  interaction_mode_plot = TH1D("interaction_mode", "interaction_mode", 100, 0, 0);
  interaction_mode_zero = TH1D("interaction_mode_zero", "interaction_mode_zero;bsenergy", 100, 0, 0);
  NCQE_plot = TH1D("NCQE_plot", "NCQE_plot;bsenergy", 100, 0, 0);
  non_NCQE_plot = TH1D("non_NCQE_plot", "non_NCQE_plot;bsenergy", 100, 0, 0);
  
  return true;
}


bool SplitAtmosInteractions::Execute(){

  LoweInfo* my_lowe_ptr = new LoweInfo();
  tree_reader_ptr->Get("LOWE", my_lowe_ptr);
  const double bsenergy = my_lowe_ptr->bsenergy;
  std::cout << "reconstructed energy: " << bsenergy << std::endl;
  
  int interaction_mode = nework_.modene;
  std::cout << "interaction mode: " << interaction_mode << std::endl;

  std::cout << NEUTInteractionModeToString(interaction_mode) << std::endl;
  
  if (bsenergy > 9000){
    return true;
  }

  interaction_mode_plot.Fill(interaction_mode);
  
  if (interaction_mode == 0){
    interaction_mode_zero.Fill(bsenergy);
  }
  else if (interaction_mode ==  51 ||
	   interaction_mode ==  52 ||
	   interaction_mode == -51 ||
	   interaction_mode == -52){
    NCQE_plot.Fill(bsenergy);
  } else {
    non_NCQE_plot.Fill(bsenergy);
  }

  //interactions_for_pi_chart[NEUTInteractionModeToString(interaction_mode)]++;
  interactions_for_pi_chart[GetNEUTModeProcess(interaction_mode)]++;
  
  return true;
}


bool SplitAtmosInteractions::Finalise(){

  TFile* output_file = TFile::Open("interaction_mode.root", "UPDATE");
  if (output_file == nullptr){
    throw std::runtime_error("SplitAtmosInteractions::Finalise: Couldn't open output file!");
  }
  
  interaction_mode_zero.Write();
  interaction_mode_plot.Write();
  NCQE_plot.Write();
  non_NCQE_plot.Write();

  std::cout << "interaction values are\n";
  for (const auto& [x,y] : interactions_for_pi_chart){
    std::cout << x << " = " << y << std::endl;
  }
  
  return true;
}

void SplitAtmosInteractions::GetReader(){
 std::string tree_reader_str = "";
  m_variables.Get("TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("SplitAtmosInteractions::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  tree_reader_ptr = m_data->Trees.at(tree_reader_str);
  return;
} 
 
