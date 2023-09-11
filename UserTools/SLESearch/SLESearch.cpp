#include "SLESearch.h"

#include "fortran_routines.h"
#include "softtrg_tbl.h"
#include "Constants.h"

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

  std::bitset<32> trigger_id{skhead_.idtgsk};
  
  if (!previous_entry_was_muon){
    Log("SLESearch::Execute: Last entry was not muon\n", 0, 0);
    previous_entry_was_muon = (trigger_id.test(static_cast<int>(TriggerType::SHE)) &&
			       trigger_id.test(static_cast<int>(TriggerType::OD_or_Fission)) &&
			       !trigger_id.test(static_cast<int>(TriggerType::AFT)));
  }
  
  previous_entry_was_muon = (trigger_id.test(static_cast<int>(TriggerType::SHE)) &&
			     trigger_id.test(static_cast<int>(TriggerType::OD_or_Fission)) &&
			     !trigger_id.test(static_cast<int>(TriggerType::AFT)));
  
  if (!trigger_id.test(static_cast<int>(TriggerType::AFT))){
    return true; 
  }

  Log("SLESearch::Execute: Running software trigger to look for SLE", 0, 0);
  int idetector[32], ithr[32], it0_offset[32] ,ipret0[32], ipostt0[32];
  softtrg_get_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);

  for (int bit_idx = 0; bit_idx < 32; ++bit_idx){
    if (bit_idx != static_cast<int>(TriggerType::SLE)){
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

  Log("SLESearch::Execute: Found " + std::to_string(ntrg) + " SLE triggers\n", 0, 0);
  
  std::vector<int> SLE_times;
  for (int trig_idx = 0; trig_idx < ntrg; ++trig_idx){
    SLE_times.push_back(swtrgtbl_.swtrgt0ctr[trig_idx]);
  }

  m_data->CStore.Set("SLE_times", SLE_times);

  int N_SLE = SLE_times.size();
  m_data->CStore.Set("N_SLE", N_SLE); //need this for the subtoolchain
  
  return true;
}


bool SLESearch::Finalise(){

  return true;
}
