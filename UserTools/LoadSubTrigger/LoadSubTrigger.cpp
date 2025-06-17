#include "LoadSubTrigger.h"

#include "fortran_routines.h"
#include <bitset>

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

  std::vector<double> SLE_times = {};
  bool ok = m_data->CStore.Get(trigger_time_names, SLE_times); // ns
  if(!ok){throw std::runtime_error("LoadSubTrigger::Execute couldn't retrieve trigger times");}
  
  double this_subtrigger_ticks_dbl = (SLE_times.at(trigger_idx) * COUNT_PER_NSEC) + skheadqb_.it0sk; // ticks
  int this_subtrigger_ticks = int(this_subtrigger_ticks_dbl);
  skheadqb_.it0xsk = this_subtrigger_ticks;
  std::cout << "this_subtrigger ticks:  " << this_subtrigger_ticks << std::endl;
  
  set_timing_gate_m_(&this_subtrigger_ticks); // ticks

  // call `skcread` to re-calculate charge in 1.3us, min and max hit time and charges and their PMTs
  int get_ok = 0;
  skcread_(&neglun, &get_ok);  // FIXME maybe skip if we can?

  // check for errors
  // get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
  if(get_ok!=0){
    Log("LoadSubTrigger::Execute:: Error! skcread returned "+std::to_string(get_ok)
	+ " when reloading SLE subtrigger!\n", 0, 0);
    return false;
  }
  
  // increment the trigger idx unless we've loaded all subtriggers, in which case set to zero, ready for the next event.
  trigger_idx < SLE_times.size() - 1 ? ++trigger_idx : trigger_idx = 0;

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
