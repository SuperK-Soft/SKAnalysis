#include "MakeBDTOutputFile.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

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
  
  std::string running_ibd_str = "";
  m_variables.Get("ibd", running_ibd_str);
  if (running_ibd_str == "yes"){carry_ibd_weights = true;}
  
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
  np_temp = 0;
  if (reader_ntagbdt_ptr->GetTree()->GetEntries() == entry_number_tmp){
    std::cout << "didn't find any neutrons on this loop" << std::endl;                                                                                         
  } else {
    
    ok = reader_ntagbdt_ptr->Get("np", np_temp);
    if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve np branch");}
    
  }
  std::cout << "np_temp: " << np_temp << std::endl;
   
  if(np_temp!=0){
    //float neutron5_temp[5]{0};
    neutron5_temp = nullptr; 
    ok = reader_ntagbdt_ptr->Get("neutron5", neutron5_temp);
    if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve neutron5 branch");}

    // copying over bdt variables arrays (np, nvx, etc) doesn't work - some sort of type mismatch
    // need to use clone tree so that the output branches  track the arrays as we get them
    // good luck 
    
    // ok = reader_ntagbdt_ptr->Get("nvx", x_temp);
    // if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve nvx branch");}
    // ok = reader_ntagbdt_ptr->Get("nvy", y_temp);
    // if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve nvy branch");}
    // ok = reader_ntagbdt_ptr->Get("nvz", z_temp);
    // if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve nvz branch");}
    // ok = reader_ntagbdt_ptr->Get("dtn", dtn_temp);
    // if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve dtn branch");}
    // ok= reader_ntagbdt_ptr->Get("bse", bse_temp);
    // if(!ok){throw std::runtime_error("MakeBDTOutputFile: Couldn't retrieve bse branch");}

    std::cout << "neutron5_temp[0]: " << neutron5_temp[0] << std::endl;
    neutron_likelihoods = std::vector<float>(neutron5_temp, neutron5_temp+np_temp);
  }  
  entry_number_tmp = reader_ntagbdt_ptr->GetTree()->GetEntries();

  std::cout << "filling output tree" << std::endl;
  output_tree_ptr->Fill();
  
  return true;
}


bool MakeBDTOutputFile::Finalise(){

  std::cout << "saving to output file: " << outfile_str << std::endl;
  
  output_file_ptr->cd();
  output_tree_ptr->Write();
  
  return true;
}

void MakeBDTOutputFile::SetOutputBranches(TTree* tree){

  bool ok = true;
  
  ok = reader_input_ptr->Get("LOWE", lowe_ptr);
  if (!ok){throw std::runtime_error("couldn't get LOWE branch");}
  tree->Branch("LOWE", "LoweInfo", lowe_ptr, 1024*1024, 0);
  
  ok = reader_input_ptr->Get("HEADER", header_ptr);
  if (!ok){throw std::runtime_error("couldn't get HEADER branch");}
  tree->Branch("HEADER", "Header", header_ptr, 1024*1024, 0);
  
  ok = reader_input_ptr->Get("MC", mcinfo_ptr);
  if (!ok){throw std::runtime_error("couldn't get MC branch");}
  tree->Branch("MC", "MCInfo", mcinfo_ptr, 1024*1024, 0);

  ok = reader_input_ptr->Get("SECONDARY", secondary_info_ptr);
  if (!ok){throw std::runtime_error("couldn't get Secondary branch");}
  tree->Branch("SECONDARY", "SecondaryInfo", secondary_info_ptr, 1024*1024, 0);

  ok = reader_input_ptr->Get("TQREAL", tqreal_ptr);
  if (!ok){throw std::runtime_error("couldn't get TQREAL branch");}
  tree->Branch("TQREAL", "TQReal", tqreal_ptr, 1024*1024, 0);

  ok = reader_input_ptr->Get("TQAREAL", tqareal_ptr);
  if (!ok){throw std::runtime_error("couldn't get TQAREAL  branch");}
  tree->Branch("TQAREAL", "TQAReal", tqareal_ptr, 1024*1024, 0);

  // ok = reader_input_ptr->Get("TQLIST", tqlist_ptr);
  // if (!ok){throw std::runtime_error("couldn't get TQLIST branch");}
  // tree->Branch("TQLIST", "TQList", tqlist_ptr, 1024*1024, 0);
  
  // ok = reader_input_ptr->Get("ODTQLIST", odtqlist_ptr);
  // if (!ok){throw std::runtime_error("couldn't get ODTQLIST branch");}
  // tree->Branch("ODTQLIST", "ODTQList", odtqlist_ptr, 1024*1024, 0);
  
  tree->Branch("neutron_likelihoods", &neutron_likelihoods, 1024*1024, 0);
  // tree->Branch("neutron5", neutron5_temp, "neutron5[np]/F", 0);
  // tree->Branch("np", &np_temp);
  // tree->Branch("nvx", x_temp, "nvx[np]/F", 0);
  // tree->Branch("nvy", y_temp, "nvy[np]/F", 0);
  // tree->Branch("nvz", z_temp, "nvz[np]/F", 0);
  // tree->Branch("dtn", dtn_temp, "dtn[np]/F", 0);
  // tree->Branch("bse", bse_temp, "bse[np]/F", 0);

  if (!carry_ibd_weights){return;}

  for (auto&& [name, val] : srn_weights){
    reader_input_ptr->Get(name.c_str(), val);
    
    tree->Branch(name.c_str(), val);
  }
  
  return;
  }
