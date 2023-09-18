#include "CalculatePreactivityObservables.h"

// #include <deque>

CalculatePreactivityObservables::CalculatePreactivityObservables():Tool(){}

bool CalculatePreactivityObservables::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  return true;
}

bool CalculatePreactivityObservables::Execute(){

  /*
    - get arrays of hit times 
    - from first hit count how many subsequent hits are less than 15ns away
    - construct a deque from these hits
    - * maxpre = deque.size()
    - if first hit in deque is within 1.3us from main trigger, maxpregate = deque.size()
    - find the time, dt, between the last hit in the window and the next one in the total readout
    - add the new hit to end of the deque
    - ** remove the hits from start that are within dt from the first hit.
    - loop from * to ** until the last hit in the deque is 12ns before the main trigger
  */

  std::vector<double> window;

  int max_pre = 0;
  int max_pregate = 0;

  int last_hit_idx = 0;
  
  for (int i = 0; i < sktqz_.nqiskz; ++i){
    if (sktqz_.tiskz[i] < 15){
      window.push_back(sktqz_.tiskz[i]);
    } else {
      last_hit_idx = i;
      break;
    }
  }

  while ((skheadqb_.it0sk - window.back()) < 12){
  
    if (window.size() > max_pre){
      max_pre = window.size();
    }

    if ((window.size() > max_pregate) && ((window.front() - skheadqb_.it0sk < 1300))){
      max_pregate = window.size();
    }

    double dt_to_next_hit = sktqz_.tiskz[last_hit_idx + 1] - window.back();
    window.push_back(sktqz_.tiskz[last_hit_idx + 1]);

    double current_first_hit = window.front();
  
    for (auto hit_it = window.begin(); hit_it != window.end(); ++hit_it){
      if (*hit_it - current_first_hit < dt_to_next_hit){
	hit_it = window.erase(hit_it);
      } else {
	break;
      }
    }

    ++last_hit_idx;
    
  }

  m_data->CStore.Set("max_pre", max_pre);
  m_data->CStore.Set("max_pregate", max_pregate);

  return true;
}

bool CalculatePreactivityObservables::Finalise(){

  return true;
}
