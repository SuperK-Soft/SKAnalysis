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
  int TreeReaderLUN = GetReaderLUN();
  int neglun = -std::abs(TreeReaderLUN);

  // call `skcread` to re-load common blocks for this subtrigger
  int get_ok = 0;
  skcread_(&neglun, &get_ok);

  // check for errors
  // get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
  if(get_ok!=0){
    Log("LoadSubTrigger::Execute:: Error! skcread returned "+std::to_string(get_ok)
	+ " when reloading SLE subtrigger!\n", 0, 0);
    return false;
  }

  ++trigger_idx;
  
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
