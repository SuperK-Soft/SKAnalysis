#include "LoadSubTrigger.h"

#include "fortran_routines.h"

LoadSubTrigger::LoadSubTrigger():Tool(){}

bool LoadSubTrigger::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  return true;
}


bool LoadSubTrigger::Execute(){
   
  std::vector<int> SLE_times = {};
  m_data->CStore.Get("SLE_times", SLE_times);
  int this_subtrigger_ticks = SLE_times.at(trigger_idx);
  
  set_timing_gate_(&this_subtrigger_ticks);

  // get a negated version of the logic unit number for the relevant file / TreeReader
  int TreeReaderLUN = m_data->GetLUN("reader");
  int neglun = -std::abs(TreeReaderLUN);

  // call `skcread` to re-load common blocks for this subtrigger
  int get_ok = 0;
  skcread_(&neglun, &get_ok);

  // check for errors
  // get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
  if(get_ok!=0){
    std::cout << " Error! skcread returned "+std::to_string(get_ok)
	      << " when reloading SLE subtrigger!\n";
    return false;
  }

  ++trigger_idx;
  
  return true;
}


bool LoadSubTrigger::Finalise(){

  return true;
}
