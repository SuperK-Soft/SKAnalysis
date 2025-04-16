#include "GetSubTriggers.h"

GetSubTriggers::GetSubTriggers():Tool(){}


bool GetSubTriggers::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  return true;
}


bool GetSubTriggers::Execute(){

  int ntrigsfound=0;
  int MAX_SUBTRIGS=32;
  std::vector<int> t0_sub(MAX_SUBTRIGS,-1);  // relative time of subtrigger to IT0SK                                                                           
  int SLE_idx = TriggerType::SLE;

  get_sub_triggers_(&SLE_idx, &ntrigsfound, t0_sub.data(), &MAX_SUBTRIGS);

  t0_sub.resize(ntrigsfound);

  std::cout << "GetSubTriggers found " << ntrigsfound << " triggers, they were:" << std::endl;
  for (const auto& t: t0_sub){std::cout << t << std::endl;}
  
  m_data->CStore.Set("trigger_times_soft", t0_sub);
  m_data->CStore.Set("n_subtriggers_soft", ntrigsfound);
  
  return true;
}


bool GetSubTriggers::Finalise(){

  return true;
}
