#include "MakeBDTOutputFile.h"

#include "TFile.h"
#include "TTree.h"

#include "loweroot.h"
#include "mcinfo.h"
#include "tqrealroot.h"

MakeBDTOutputFile::MakeBDTOutputFile():Tool(){}


bool MakeBDTOutputFile::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  std::string reader_str = "";
  m_variables.Get("reader_input", reader_str);
  if (m_data->Trees.count(reader_str) == 0 || reader_str == ""){
    throw std::runtime_error("MakeBDTOutputFile: failed to get input treereader!");
  }
  reader_input_ptr = m_data->Trees.at(reader_str);

  if (m_data->Trees.count("ntag_BDT_OutTree") == 0){throw std::runtime_error("MakeBDTOutputFile: couldn't get bdt tree");}
  reader_ntagbdt_ptr = m_data->Trees.at("ntag_BDT_OutTree");
  
  std::string outfile_str = "";
  m_variables.Get("outfile", outfile_str);
  if (outfile_str.empty()){throw std::runtime_error("MakeBDTOutputFile: no output file specified!");}

  output_file_ptr = new TFile(outfile_str.c_str(), "RECREATE");
  output_tree_ptr = new TTree("data", "data");

  SetOutputBranches(output_tree_ptr);
  
  return true;  
}


bool MakeBDTOutputFile::Execute(){

  std::cout << "reader_ntagbdt_ptr->GetTree(): " << reader_ntagbdt_ptr->GetTree() << std::endl;
  std::cout << "reader_ntagbdt_ptr->GetTree()->GetEntries(): " << reader_ntagbdt_ptr->GetTree()->GetEntries() << std::endl;
  
  neutron_likelihoods.clear();
  int np_temp = 0;
  if (reader_ntagbdt_ptr->GetTree()->GetEntries() == entry_number_tmp){
    std::cout << "didn't find any neutrons on this loop" << std::endl;                                                                                         
  } else {
    
    ok = reader_ntagbdt_ptr->Get("np", np_temp);
    if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve np branch");}
    
  }
  std::cout << "np_temp: " << np_temp << std::endl;
   
  if(np_temp!=0){
    float neutron5_temp[500]{0};
    //float neutron5_temp = nullptr; 
    ok = reader_ntagbdt_ptr->Get("neutron5", neutron5_temp);
    if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve neutron5 branch");}
    std::cout << "neutron5_temp[0]: " << neutron5_temp[0] << std::endl;
    neutron_likelihoods = std::vector<float>(neutron5_temp, neutron5_temp+np_temp);
  }  
  entry_number_tmp = reader_ntagbdt_ptr->GetTree()->GetEntries();
  
  output_tree_ptr->Fill();
  
  return true;
}


bool MakeBDTOutputFile::Finalise(){

  output_file_ptr->cd();
  output_tree_ptr->Write();
  
  return true;
}

void MakeBDTOutputFile::SetOutputBranches(TTree* tree) {
  reader_input_ptr->Get("LOWE", lowe_ptr);
  tree->Branch("LOWE", "LoweInfo", lowe_ptr, 1024*1024, 0);
  
  reader_input_ptr->Get("HEADER", header_ptr);
  tree->Branch("HEADER", "Header", header_ptr, 1024*1024, 0);
  
  reader_input_ptr->Get("MC", mcinfo_ptr);
  tree->Branch("MC", "MCInfo", mcinfo_ptr, 1024*1024, 0);

  reader_input_ptr->Get("TQREAL", tqreal_ptr);
  tree->Branch("TQREAL", "TQReal", tqreal_ptr, 1024*1024, 0);

  // reader_ptr->Get("neutron5", neutron5);
  // tree.Branch("neutron5", neutron5, "neutron5[np]/F");

  tree->Branch("neutron_likelihoods", &neutron_likelihoods, 1024*1024, 0);
  
  return;
}
