#include "CalculatePreactivityObservables.h"

#include <algorithm>

#include "TableReader.h"
#include "TableEntry.h"

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
    - construct a vector from these hits
    - * maxpre = vector.size()
    - if first hit in vector is within 1.3us from main trigger, maxpregate = vector.size()
    - find the time, dt, between the last hit in the window and the next one in the total readout
    - add the new hit to end of the vector
    - ** remove the hits from start that are within dt from the first hit.
    - loop from * to ** until the last hit in the vector is 12ns before the main trigger
  */
  
  const auto tof = [](const float* x, const float* y){
    double dist = 0;
    for (int i = 0; i < 3; ++i){
      dist += pow(x[i] - y[i], 2);
    }
    return dist / (3 * pow(10, 8));
  };  
  
  ConnectionTable* connection_table = m_data->GetConnectionTable();

  std::vector<double> tof_sub_times = {};
  for (int pmt_idx = 0; pmt_idx < sktqz_.nqiskz; ++pmt_idx){
    const int cable_number = sktqz_.icabiz[pmt_idx];
    float* pmt_loc = nullptr;
    connection_table->GetTubePosition(cable_number, pmt_loc);
    tof_sub_times.emplace_back(sktqz_.tiskz[pmt_idx] - tof(skroot_lowe_.bsvertex, pmt_loc));
  }
  std::sort(tof_sub_times.begin(), tof_sub_times.end());
 
  std::vector<double> window = {};
  int last_hit_idx = 0;
  
  int max_pre = 0;
  int max_pregate = 0;
    
  for (int i = 0; i < sktqz_.nqiskz; ++i){
    if (sktqz_.tiskz[i] < 15){ //tof
      float* pmt_loc = nullptr;
      window.emplace_back(sktqz_.tiskz[i] - tof(skroot_lowe_.bsvertex, pmt_loc));
    } else {
      last_hit_idx = i - 1;
      break;
    }
  }

  while (((skheadqb_.it0sk - window.back()) < 12 && (last_hit_idx != tof_sub_times.size()))){
  
    if (window.size() > max_pre){
      max_pre = window.size();
    }

    if ((window.size() > max_pregate) && ((window.front() - skheadqb_.it0sk < 1300))){
      max_pregate = window.size();
    }

    const double dt_to_next_hit = tof_sub_times.at(last_hit_idx + 1) - window.back();
    window.push_back(tof_sub_times.at(last_hit_idx + 1));

    const double current_first_hit = window.front();
  
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
