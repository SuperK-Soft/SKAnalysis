#include "DefineSignalRegions.h"

#include "MTreeReader.h"

DefineSignalRegions::DefineSignalRegions():Tool(){}


bool DefineSignalRegions::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
  
  r0 = TH1D("r0", "20#circ < #theta_{c} < 38#circ : N_{tagged} != 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 100);
  r1 = TH1D("r1", "38#circ < #theta_{c} < 53#circ : N_{tagged} != 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 100);
  r2 = TH1D("r2", "70#circ < #theta_{c} < 90#circ : N_{tagged} != 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 100);
  r3 = TH1D("r3", "20#circ < #theta_{c} < 38#circ : N_{tagged} == 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 100);
  r4 = TH1D("r4", "38#circ < #theta_{c} < 53#circ : N_{tagged} == 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 100);
  r5 = TH1D("r5", "70#circ < #theta_{c} < 90#circ : N_{tagged} == 1;Reconstructed Energy [Mev];#Events / 0.01MeV", 100, 0, 100);

  need_to_get_treereaders = true;
  
  return true;
}


bool DefineSignalRegions::Execute(){

  try{

    if (need_to_get_treereaders){
      GetTreeReaders();
    }
    
    m_data->vars.Set("Skip", false);
    
    LoweInfo* lowe_ptr = nullptr;
    bool ok = LOWE_TreeReader->Get("LOWE", lowe_ptr);
    if (!ok || lowe_ptr == nullptr){
      throw std::runtime_error("DefineSignalRegions::Execute - Failed to get LOWE branch!");
    }

    Header* header_ptr = nullptr;
    ok = LOWE_TreeReader->Get("HEADER", header_ptr);
    if (!ok || header_ptr == nullptr){
      throw std::runtime_error("DefineSignalRegions::Execute - Failed to get HEADER branch!");
    }
    
    current_lowe_event_number = header_ptr->nevsk;
    
    ok = BDT_TreeReader->Get("HEADER", header_ptr);
    if (!ok){
      throw std::runtime_error("DefineSignalRegions::Execute - Failed to get HEADER branch!");
    }

    current_bdt_event_number = header_ptr->nevsk;

    std::cout << "current_lowe_event_number: " << current_lowe_event_number << "\n";
    std::cout << "current_bdt_event_number:  " << current_bdt_event_number << "\n";
    std::cout << "LOWE_TreeReader->GetEntryNumber(): " << LOWE_TreeReader->GetEntryNumber() << "\n";
    std::cout << "BDT_TreeReader->GetEntryNumber():  " << BDT_TreeReader->GetEntryNumber() << "\n";
    //std::cout << "LOWE_TreeReader->GetEntries(): " << LOWE_TreeReader->GetEntries() << "\n";
    //std::cout << "BDT_TreeReader->GetEntries():  " << BDT_TreeReader->GetEntries() << "\n";
    std::cout << "events_matched: " << events_matched << "\n";

    
    if (current_bdt_event_number != current_lowe_event_number){
      std::cout << "DefineSignalRegions:Execute - not yet matched , skipping rest of toolchain\n";

      
      m_data->vars.Set("Skip", true);
      /*
      if (LOWE_TreeReader->GetEntries() == LOWE_TreeReader->GetEntryNumber() + 1 &&
	  BDT_TreeReader->GetEntries() == BDT_TreeReader->GetEntryNumber() + 1 &&
	  BDT_TreeReader->GetEntries() != events_matched){
	throw std::runtime_error("DefineSignalRegions:Execute - didn't match all events");
      }
      */


      return true;
    }

    std::cout << "DefineSignalRegions:Execute - matched event! filling histograms!\n";
    ++events_matched;
    std::cout << "events_matched: " << events_matched << "\n";
    
    // if (LOWE_TreeReader->GetEntries() - 1 == LOWE_TreeReader->GetEntryNumber() && BDT_TreeReader->GetEntryNumber() + 1 != events_matched){
    //   std::cout << "current_bdt_event_number: " << current_bdt_event_number << "\n";
    //   std::cout << "current_lowe_event_number: " << current_lowe_event_number << "\n";
    //   std::cout << "events found: " << events_matched << "\n";
    //   throw std::runtime_error("DefineSignalRegions:Execute - didn't match all events");
    // }
    
    const double cherenkov_angle = acos(lowe_ptr->bsresult[3]) * 180 / 3.1415926535;
    std::cout << "angle: " << cherenkov_angle << "\n";

    //    double neutron5[np];
    basic_array<float*> neutron5;
    ok = BDT_TreeReader->GetBranchValue("neutron5", neutron5);
    if (!ok){
      throw std::runtime_error("DefineSignalRegions::Execute - Couldn't Get() variable neutron5");
    }

    //    const std::vector<double> likelihoods(neutron5);
    //const bool has_neutron = HasExactlyOneNeutron(likelihoods);
    const bool has_neutron = HasExactlyOneNeutron(neutron5);
    
    if (lowe_ptr->bsenergy > 1000){
      return true;
    }

    if (cherenkov_angle > 20 && cherenkov_angle < 38){
      has_neutron ? r3.Fill(lowe_ptr->bsenergy) : r0.Fill(lowe_ptr->bsenergy);
    }

    if (cherenkov_angle > 38 && cherenkov_angle < 53){
      has_neutron ? r4.Fill(lowe_ptr->bsenergy) : r1.Fill(lowe_ptr->bsenergy);
    }
    
    if (cherenkov_angle > 70 && cherenkov_angle < 90){
      has_neutron ? r5.Fill(lowe_ptr->bsenergy) : r2.Fill(lowe_ptr->bsenergy);
    }
    
  }
  catch(const std::exception& e){
    std::cout << e.what() << "\n";
    return false;
  }
  
  return true;
}


bool DefineSignalRegions::Finalise(){

  TFile* out = new TFile("regions.root", "recreate");
  out->cd();

  r0.Write();
  r1.Write();
  r2.Write();
  r3.Write();
  r4.Write();
  r5.Write();

  MakePlot();
  
  out->Close();
  delete out; out = nullptr;
  
  return true;
}

void DefineSignalRegions::GetTreeReaders(){
  std::string tree_reader_str = "";
  m_variables.Get("LOWE_TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("DefineSignalRegions::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  LOWE_TreeReader = m_data->Trees.at(tree_reader_str);

  m_variables.Get("BDT_TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("DefineSignalRegions::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
    
  BDT_TreeReader = m_data->Trees.at(tree_reader_str);
  need_to_get_treereaders = false;    
}

bool DefineSignalRegions::HasExactlyOneNeutron(const basic_array<float*> likelihoods){
  double cut = -1;
  bool ok = m_variables.Get("likelihood_cut", cut);
  if (!ok || cut == -1){
    throw std::runtime_error("DefineSignalRegions::HasOneNeutron - Couldn't get neutron likelihood cut value!");
  }
  const int n_neutrons = std::count_if(likelihoods.begin(), likelihoods.end(), [cut](double l){return l > cut;});
  return n_neutrons == 1;
}

void DefineSignalRegions::MakePlot(){
  TFile* out = new TFile("regions_plot.root", "recreate");
  out->cd();
 
  TCanvas c1("c1", "c1", 3000, 2000);
  std::vector<TH1D> rs = {r0, r1, r2, r3, r4, r5};
  c1.Divide(3,2);
  c1.SetGrid();
  for (int i = 0; i < rs.size(); ++i){
    c1.cd(i+1);

    rs.at(i).Draw();
  }

  c1.Write();
  
  out->Close();
  delete out; out = nullptr;

  return;
}
