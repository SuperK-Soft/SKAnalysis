#include "LookForSHEAFT.h"

#include <bitset>

LookForSHEAFT::LookForSHEAFT():Tool(){}

bool LookForSHEAFT::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  return true;
}


bool LookForSHEAFT::Execute(){

  std::string input_fname = "/home/mattnich/disk3stor/prepared_data/sk2p2/rfm_run086613.001132_sk2p2.root";
  //m_variables.Get("fname", input_fname);

  std::cout << "got file" << std::endl;
  
  TFile* file_ptr = TFile::Open(input_fname.c_str(), "open");
  TTree* data_tree_ptr = static_cast<TTree*>(file_ptr->Get("sk2p2"));

  std::cout << "got tree" << std::endl;
  
  Header* header_ptr = new Header();
  data_tree_ptr->SetBranchAddress("HEADER", &header_ptr);
  if (header_ptr == nullptr){
    throw std::runtime_error("couldn't get header");
  }

  std::cout << "got header" << std::endl;
  
  for(int i=0; i<data_tree_ptr->GetEntries(); ++i){
  data_tree_ptr->GetEntry(i);
  std::bitset<32> bits{header_ptr->idtgsk};
  if(bits.test(28) || bits.test(29)){
    std::cout<<"Found one in entry "<<i<<" with nevsk "<< header_ptr->nevsk <<  std::endl;  
   }
  else {
    std::cout << "Dint find nowt in this one (aka, entry " << i << " chief" << std::endl;
  }
}					     

  std::cout << "done" << std::endl;

  delete header_ptr; header_ptr = nullptr;
  
  return true;
}


bool LookForSHEAFT::Finalise(){

  return true;
}
