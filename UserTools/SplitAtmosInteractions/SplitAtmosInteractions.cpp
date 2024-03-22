#include "SplitAtmosInteractions.h"
//#include "neworkC.h"

#include "loweroot.h"
#include "MTreeReader.h"

SplitAtmosInteractions::SplitAtmosInteractions():Tool(){}


bool SplitAtmosInteractions::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReader();
  
  return true;
}


bool SplitAtmosInteractions::Execute(){

  LoweInfo* my_lowe_ptr = nullptr;
  tree_reader_ptr->Get("LOWE", my_lowe_ptr);
  const double bsenergy = my_lowe_ptr->bsenergy;
  std::cout << "reconstructed energy: " << bsenergy << std::endl;
  
  // int interaction_mode = nework_.modene;
  // std::cout << "interaction mode: " << interaction_mode << std::endl;
  
  return true;
}


bool SplitAtmosInteractions::Finalise(){
 std::string tree_reader_str = "";
  m_variables.Get("TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("SplitAtmosInteractions::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  tree_reader_ptr = m_data->Trees.at(tree_reader_str);
  return true;
}

void SplitAtmosInteractions::GetReader(){
  return;
}
