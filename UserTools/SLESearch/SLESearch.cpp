#include "SLESearch.h"

#include "fortran_routines.h"
#include "softtrg_tbl.h"

#include <bitset>

SLESearch::SLESearch():Tool(){}

bool SLESearch::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  return true;
}


bool SLESearch::Execute(){

  const int SHE_idx = 28, OD_idx = 3, AFT_idx = 29;
  std::bitset<32> trigger_id{skhead_.idtgsk};
  
  if (!previous_entry_was_muon){
    std::cout << "last entry was not muon\n";
    previous_entry_was_muon = (trigger_id.test(SHE_idx) &&
			       trigger_id.test(OD_idx) &&
			       !trigger_id.test(AFT_idx));
  }
    
  previous_entry_was_muon = (trigger_id.test(SHE_idx) &&
			     trigger_id.test(OD_idx) &&
			     !trigger_id.test(AFT_idx));

  if (!trigger_id.test(AFT_idx)){
    return true; 
  }
  
  std::cout << "running software trigger for SLE triggers\n";
  int idetector[32], ithr[32], it0_offset[32] ,ipret0[32], ipostt0[32];
  softtrg_get_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);

  const int SLE_idx = 2;
  for (int bit_idx = 0; bit_idx < 32; ++bit_idx){
    if (bit_idx != SLE_idx){
      ithr[bit_idx] = 100000;
      it0_offset[bit_idx] = 0;
      ipret0[bit_idx] = 0;
      ipostt0[bit_idx] = 0;
    }
  }

  softtrg_set_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);

  int max_qb = 1280;
  int one = 1;
  int zero = 0;
  int ntrg = softtrg_inittrgtbl_(&skhead_.nrunsk, &zero, &one, &max_qb);

  std::cout << "SLESearch found " << std::to_string(ntrg) << " SLE triggers\n";
  
  std::vector<int> SLE_times;
  for (int trig_idx = 0; trig_idx < ntrg; ++trig_idx){
    SLE_times.push_back(swtrgtbl_.swtrgt0ctr[trig_idx]);
  }

  m_data->SLE_times = SLE_times;
  
  return true;
}


bool SLESearch::Finalise(){

  return true;
}
