#include "GetSubTrg_SetTimingGate.h"

GetSubTrg_SetTimingGate::GetSubTrg_SetTimingGate():Tool(){}


bool GetSubTrg_SetTimingGate::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  return true;
}


bool GetSubTrg_SetTimingGate::Execute(){

  int ntrigsfound=0;
  int MAX_SUBTRIGS=100;
  std::vector<int> t0_sub(MAX_SUBTRIGS,-1);  // relative time of subtrigger to IT0SK                                                             
  int SLE_idx = TriggerType::SLE;
  
  get_sub_triggers_(&SLE_idx, &ntrigsfound, t0_sub.data(), &MAX_SUBTRIGS);
  int newt0 = t0_sub.at(0) + skheadqb_.it0sk;

  t0_sub.resize(ntrigsfound);

  std::cout << "from GetSubTrg_SetTimingGate: " << std::endl;
  std::cout << "found " << ntrigsfound << " subtriggers " << std::endl;
  for (const auto& t : t0_sub){std::cout << "subtrigger: " << t + skheadqb_.it0sk << std::endl;}
  
  m_data->CStore.Set("trigger_times", t0_sub);
  m_data->CStore.Set("n_subtriggers", ntrigsfound);
  
  set_timing_gate_(&newt0);
  
  return true;
}


bool GetSubTrg_SetTimingGate::Finalise(){

  return true;
}
