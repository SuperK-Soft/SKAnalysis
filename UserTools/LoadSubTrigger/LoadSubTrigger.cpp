#include "LoadSubTrigger.h"

#include "fortran_routines.h"
#include <bitset>

extern "C" void set_timing_gate_proc_(int*);

LoadSubTrigger::LoadSubTrigger():Tool(){}

bool LoadSubTrigger::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  // get a negated version of the logic unit number for the relevant file / TreeReader
  TreeReaderLUN = GetReaderLUN();
  neglun = -std::abs(TreeReaderLUN);

  bool ok = m_variables.Get("trigger_time_names", trigger_time_names);
  if (!ok || trigger_time_names.empty()){throw std::runtime_error("LoadSubTrigger::Execute - on trigger time variables name selected!");}
  
  return true;

}


bool LoadSubTrigger::Execute(){

  if (!made_plots){
    hits_before_inwindow = TH1D("hits_before_inwindow", "hits_before_inwindow", 100,
		      *std::min_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz),
		       *std::max_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz));
    hits_before_outwindow = TH1D("hits_before_outwindow", "hits_before_outwindow", 100,
		      *std::min_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz),
		      *std::max_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz));
    
    for (int i = 0; i < sktqz_.nqiskz; ++i){
      if((sktqz_.ihtiflz[i] & 0x01)==0){
	hits_before_outwindow.Fill(sktqz_.tiskz[i]);
      } else {
	hits_before_inwindow.Fill(sktqz_.tiskz[i]);
      }
    }
    std::cout << "skheadqb_.it0xsk: " << skheadqb_.it0xsk << std::endl;
  }
  
  std::vector<double> SLE_times = {};
  bool ok = m_data->CStore.Get(trigger_time_names, SLE_times); // ns
  if(!ok){throw std::runtime_error("LoadSubTrigger::Execute couldn't retrieve trigger times");}
  
  double this_subtrigger_ticks = (SLE_times.at(trigger_idx) * COUNT_PER_NSEC) + skheadqb_.it0sk; // ticks
  //  skheadqb_.it0xsk = this_subtrigger_ticks;
  std::cout << "this_subtrigger ticks:  " << this_subtrigger_ticks << std::endl;
  
  // int icabbf_counter = 0;
  // int itiskz_counter = 0;
  // int id_start_counter = 0;
  // int id_end_counter = 0;
  // int total_counter = 0;

  // std::cout << "nrunsk: "  << skhead_.nrunsk << std::endl;
  // std::cout << "SKGATE_START_COUNT: " << SKGATE_START_COUNT << std::endl;
  // std::cout << "SKGATE_END_COUNT: " << SKGATE_END_COUNT << std::endl;
  
  // for (int i = 0; i < 5; ++i){
  //   std::cout << "sktqz_.itiskz["<<i<<"]: " << sktqz_.itiskz[i] << std::endl;
  //   std::cout << "rawtqinfo_.itiskz_raw["<<i<<"]: " << rawtqinfo_.itiskz_raw[i] << std::endl;
  // }
  
  // for (int i = 0; i < sktqz_.nqiskz; ++i){
  //   if (sktqz_.icabiz[i] > MAXPM){
  //     continue;
  //   }
  //   int count = 0;
  //   std::bitset<16> bits(rawtqinfo_.icabbf_raw[i] >> 16);
    
  //   if (bits.test(2)){
  //     ++icabbf_counter;
  //     ++count;
  //   }
  //   if (sktqz_.tiskz[i] >= SKGATE_START_COUNT + (SLE_times.at(trigger_idx) * COUNT_PER_NSEC)){
  //     ++id_start_counter;
  //     ++count;
  //   }
  //   if (sktqz_.tiskz[i] <= SKGATE_END_COUNT + (SLE_times.at(trigger_idx) * COUNT_PER_NSEC)){
  //     ++id_end_counter;
  //     ++count;
  //   }
  //   if (count == 3){++total_counter;}
  // }

  // std::cout << "if checks:" << std::endl;
  // std::cout << "icabbf: " << icabbf_counter << std::endl;
  // std::cout << "id start: " << id_start_counter << std::endl;
  // std::cout << "id end: " << id_end_counter << std::endl;
  // std::cout << "total: " << total_counter << "\\" << sktqz_.nqiskz << std::endl;
  
  
  int this_subtrigger_ticks_tmp = int(this_subtrigger_ticks);
  set_timing_gate_proc_(&this_subtrigger_ticks_tmp); // ticks

  //set_timing_gate_(&skheadqb_.it0sk);

  // call `skcread` to re-load common blocks for this subtrigger
  int get_ok = 0;
  skcread_(&neglun, &get_ok);

  // for (int i = 0; i < sktqz_.nqiskz; ++i){
  //   hits_before.Fill(sktqz_.tiskz[skchnl_.ihcab[i]]);
  // }

  // for (int j = 0; j < sktqz_after.nqiskz; ++j){
  //   hits_after.Fill(sktqz_after.tiskz[skchnl_after.ihcab[j]]);
  // }
  
  // check for errors
  // get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
  if(get_ok!=0){
    Log("LoadSubTrigger::Execute:: Error! skcread returned "+std::to_string(get_ok)
	+ " when reloading SLE subtrigger!\n", 0, 0);
    return false;
  }
  
  // increment the trigger idx unless we've loaded all subtriggers, in which case set to zero, ready for the next event.
  trigger_idx < SLE_times.size() - 1 ? ++trigger_idx : trigger_idx = 0;


  // std::cout << "COMPARING SKTQZ" << std::endl;
  // CompareSKTQZ(sktqz_, skchnl_, sktqz_after, skchnl_after);
  // std::cout << std::endl << "COMPARING SKQ" << std::endl;
  // Compare_skq(skq_, skchnl_, m_data->skq_common_dupl, m_data->skchnl_common_dupl);
  // std::cout << std::endl << "COMPARING SKT" << std::endl;
  // Compare_skt(skt_, skchnl_, m_data->skt_common_dupl, m_data->skchnl_common_dupl);

  // if (!saved){
  //   saved = true;
  //   hits_before.SaveAs("hits_before.root");
  //   hits_after.SaveAs("hits_after.root");
  // }

  if (!made_plots){
    hits_after_inwindow = TH1D("hits_after_inwindow", "hits_after_inwindow", 100,
		      *std::min_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz),
		      *std::max_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz));
    hits_after_outwindow = TH1D("hits_after_outwindow", "hits_after_outwindow", 100,
		      *std::min_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz),
		      *std::max_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz));

    for (int i = 0; i < sktqz_.nqiskz; ++i){
      if((sktqz_.ihtiflz[i] & 0x01)==0){
      hits_after_outwindow.Fill(sktqz_.tiskz[i]);
      } else {
	hits_after_inwindow.Fill(sktqz_.tiskz[i]);
      }
    }
    std::cout << "skheadqb_.it0xsk: " << skheadqb_.it0xsk << std::endl;
    made_plots = true;
    
    hits_before_inwindow.SaveAs("hits_before_inwindow.root");
    hits_after_inwindow.SaveAs("hits_after_inwindow.root");
    hits_before_outwindow.SaveAs("hits_before_outwindow.root");
    hits_after_outwindow.SaveAs("hits_after_outwindow.root");

  }

  
  return true;
}


bool LoadSubTrigger::Finalise(){
  
  return true;
}

int LoadSubTrigger::GetReaderLUN(){
  std::string reader_name = "";
  m_variables.Get("reader_name", reader_name);
  if (m_data->Trees.count(reader_name) != 1){
    throw std::runtime_error("LoadSubTrigger::GetReaderLUN couldn't get the LUN of a valid TreeReader");
  }
  return m_data->GetLUN(reader_name);
}

void LoadSubTrigger::CompareSKTQZ(const sktqz_common& before, const skchnl_common& before_chnl,
				  const sktqz_common& after, const skchnl_common& after_chnl){
  
  if (before.nqiskz == after.nqiskz){
    std::cout << "sktqz_.nqiskz left unchanged" << std::endl;
  } else {
    std::cout << "sktqz_.nqiskz has changed:" << std::endl;
    std::cout << "before: " << before.nqiskz << ", after: " << after.nqiskz << std::endl;
  }

  if (CompareArray(before.ihtiflz, after.ihtiflz, 100)){
    std::cout << "sktqz_.ihtiflz left unchanged" << std::endl;
  } else {
    std::cout << "sktqz_.ihtiflz has changed: " << std::endl;
    PrintPMTArray(before.ihtiflz, before_chnl, after.ihtiflz, after_chnl);
  }

  if (CompareArray(before.icabiz, after.icabiz, 100)){
    std::cout << "sktqz_.icabiz left unchanged" << std::endl;
  } else {
    std::cout << "sktqz_.icabiz has changed: " << std::endl;
    PrintArray(before.icabiz, after.icabiz);
  }
  
  if (CompareArray(before.itiskz, after.itiskz, 100)){
    std::cout << "sktqz_.itiskz left unchanged" << std::endl;
  } else {
    std::cout << "sktqz_.itiskz has changed: " << std::endl;
    PrintPMTArray(before.itiskz, before_chnl, after.itiskz, after_chnl);
  }

  if (CompareArray(before.iqiskz, after.iqiskz, 100)){
    std::cout << "sktqz_.iqiskz left unchanged" << std::endl;
  } else {
    std::cout << "sktqz_.iqiskz has changed: " << std::endl;
    PrintPMTArray(before.iqiskz, before_chnl, after.iqiskz, after_chnl);
  }

  if (CompareArray(before.tiskz, after.tiskz, 100)){
    std::cout << "sktqz_.tiskz left unchanged" << std::endl;
  } else {
    std::cout << "sktqz_.tiskz has changed: " << std::endl;
    PrintPMTArray(before.tiskz, before_chnl, after.tiskz, after_chnl);
  }

  if (CompareArray(before.qiskz, after.qiskz, 100)){
    std::cout << "sktqz_.qiskz left unchanged" << std::endl;
  } else {
    std::cout << "sktqz_.qiskz has changed: " << std::endl;
    PrintPMTArray(before.qiskz, before_chnl, after.qiskz, after_chnl);
  }

  return;
}

template <typename T>
bool LoadSubTrigger::CompareArray(T* a, T* b, const int& N){
  bool match = true;
  for (int i = 0; i < N; ++i){match *= (a[i] == b[i]);}
  return match;
}

template <typename T>
void LoadSubTrigger::PrintPMTArray(T* a, const skchnl_common& a_chnl,
				   T* b, const skchnl_common& b_chnl){
  for (int i = 0; i < 5; ++i){
    int a_cabN = a_chnl.ihcab[i];
    int b_cabN = b_chnl.ihcab[i];
    std::cout << "before["<<a_cabN<<"]: "<<a[a_cabN]<<", after["<<b_cabN<<"]: "<<b[b_cabN]<<std::endl;
  }
  return;
}

template <typename T>
void LoadSubTrigger::PrintArray(T* a, T* b){
  for (int i = 0; i < 5; ++i){
    std::cout << "before["<<i<<"]: "<<a[i]<<", after["<<i<<"]: "<<b[i]<<std::endl;
  }
  return;
}

void LoadSubTrigger::Compare_skq(const skq_common& orig, const skchnl_common& origchnl, const skq_common& dupl, const skchnl_common& duplchnl) const {
  // nqisk:                                                                                                                                                                                                                                                                       
  if (orig.nqisk != dupl.nqisk){
    std::cout << "skq_.nqisk does not match after copy!" << std::endl;
    std::cout << "orig.nqisk: " << orig.nqisk << ", dupl.nqisk: " << dupl.nqisk << std::endl;
  }
  // qismsk:                                                                                                                                                                                                                                                                      
  if (orig.qismsk != dupl.qismsk){
    std::cout << "skq_.qismsk does not match after copy!" << std::endl;
    std::cout << "photoelectrons" << std::endl;
    std::cout << "orig.qismsk: " << orig.qismsk << ", dupl.qismsk: " << dupl.qismsk << std::endl;
  }

  // qimxsk:                                                                                                                                                                                                                                                                      
  if (orig.qimxsk != dupl.qimxsk){
    std::cout << "skq_.qimxsk does not match after copy!" << std::endl;
    std::cout << "photoelectrons" << std::endl;
    std::cout << "orig.qimxsk: " << orig.qimxsk << ", dupl.qimxsk: " << dupl.qimxsk << std::endl;
  }

  // mxqisk:                                                                                                                                                                                                                                                                      
  if (orig.mxqisk != dupl.mxqisk){
    std::cout << "skq_.mxqisk does not match after copy!" << std::endl;
    std::cout << "pmt number that acquired max charge" << std::endl;
    std::cout << "orig.mxqisk: " << orig.mxqisk << ", dupl.mxqisk: " << dupl.mxqisk << std::endl;

  }
  // qisk[MAXPM]                                                                                                                                                                                                                                                                  
  // if (!CompareArray(orig.qisk, dupl.qisk, MAXPM)){                                                                                                                                                                                                                             
  //   std::cout << "skq_.qisk[MAXPM] does not match after copy!" << std::endl;                                                                                                                                                                                                   
  std::cout << "some values from the arrays: " << std::endl;
  std::cout << "photoelectrons" << std::endl;
  for (int i = 0; i < 20/*std::min(std::min(20, orig.nqisk), dupl.nqisk)*/; ++i){
    int orig_cableNumber = origchnl.ihcab[i];
    int dupl_cableNumber = duplchnl.ihcab[i];
    std::cout << "orig.qisk["<<orig_cableNumber<<"-1] " <<  orig.qisk[orig_cableNumber-1] << ", dupl.qisk["<<dupl_cableNumber<<"-1] " << dupl.qisk[dupl_cableNumber-1] << std::endl;
  }
  // }                                                                                                                                                                                                                                                                            
  return;
}


void LoadSubTrigger::Compare_skt(const skt_common& orig, const skchnl_common& origchnl, const skt_common& dupl, const skchnl_common& duplchnl) const {
  // timnsk:                                                                                                                                                                                                                                   
  if (orig.timnsk != dupl.timnsk){
    std::cout << "skt_.timnsk does not match after copy!" << std::endl;
    std::cout << "ns:  ";
    std::cout << "orig.timnsk: " << orig.timnsk << ", dupl.timnsk: " << dupl.timnsk << std::endl;
  }
  // timxsk:                                                                                                                                                                                                                                   
  if (orig.timxsk != dupl.timxsk){
    std::cout << "ns:  ";
    std::cout << "skt_.timxsk does not match after copy!" << std::endl;
    std::cout << "orig.timxsk: " << orig.timxsk << ", dupl.timxsk: " << dupl.timxsk << std::endl;
  }
  // mntisk:                                                                                                                                                                                                                                   
  if (orig.mntisk != dupl.mntisk){
    std::cout << "ns:  ";
    std::cout << "skt_.mntisk does not match after copy!" << std::endl;
    std::cout << "orig.mntisk: " << orig.mntisk << ", dupl.mntisk: " << dupl.mntisk << std::endl;
  }
  // mxtisk:                                                                                                                                                                                                                                   
  if (orig.mxtisk != dupl.mxtisk){
    std::cout << "ns:  ";
    std::cout << "skt_.mxtisk does not match after copy!" << std::endl;
    std::cout << "orig.mxtisk: " << orig.mxtisk << ", dupl.mxtisk: " << dupl.mxtisk << std::endl;
  }
  // tisk[MAXPM]                                                                                                                                                                                                                               
  // if (!CompareArray(orig.tisk, dupl.tisk, MAXPM)){                                                                                                                                                                                          
  //   std::cout << "skt_.tisk[MAXPM] does not match after copy!" << std::endl;                                                                                                                                                                
  std::cout << "some values from the arrays: " << std::endl;
  std::cout << "ns:  ";
  for (int i = 0; i < 20; i++){
    int orig_cableNumber = origchnl.ihcab[i];
    int dupl_cableNumber = duplchnl.ihcab[i];
    std::cout << "orig.tisk["<<orig_cableNumber<<"-1] " <<  orig.tisk[orig_cableNumber-1] << ", dupl.tisk["<<dupl_cableNumber<<"-1] " << dupl.tisk[dupl_cableNumber-1] << std::endl;}
  // }                                                                                                                                                                                                                                         
  return;
}
