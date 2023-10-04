#include "CalculatePreactivityObservables.h"


#include <algorithm>

#include "TableReader.h"
#include "TableEntry.h"
#include "Constants.h"

CalculatePreactivityObservables::CalculatePreactivityObservables():Tool(){}

bool CalculatePreactivityObservables::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  m_variables.Get("dark_threshold", dark_threshold);
  m_variables.Get("fraction", fraction);

  connection_table = m_data->GetConnectionTable();
  
  return true;
}

bool CalculatePreactivityObservables::Execute(){

  /*
    - get arrays of hit times 
    - from first hit count how many subsequent hits are less than 15ns away
    - construct a vector from these hits
    - calculate goodness values and remove the hits that don't make goodness cut
    - * maxpre = vector.size()
    - if first hit in vector is within 1.3us from main trigger, maxpregate = vector.size()
    - find the time, dt, between the last hit in the window and the next one in the total readout
    - add the new hit to end of the vector
    - ** remove the hits from start that are within dt from the first hit.
    - loop from * to ** until the last hit in the vector is 12ns before the main trigger
  */
  
  double lowest_in_gate_time = 9999;
  
  std::vector<Hit> tof_sub_hits = std::vector<Hit>(sktqz_.nqiskz, {0,0,0});
  
  for (int pmt_idx = 0; pmt_idx < sktqz_.nqiskz; ++pmt_idx){
    const int cable_number = sktqz_.icabiz[pmt_idx];
    float pmt_loc[3] = {};
    connection_table->GetTubePosition(cable_number, pmt_loc);
    const double new_time = sktqz_.tiskz[pmt_idx] - skroot_lowe_.bsvertex[3] - TimeOfFlight(skroot_lowe_.bsvertex, pmt_loc);
    if (((sktqz_.ihtiflz[pmt_idx] & 0x01)==1) && (new_time < lowest_in_gate_time)){
      lowest_in_gate_time = new_time;
    }
    tof_sub_hits.emplace_back(new_time, 0, sktqz_.qiskz[pmt_idx]); //calculate goodness in the next loop
  }

  std::sort(tof_sub_hits.begin(), tof_sub_hits.end(), [](const Hit& h1, const Hit& h2){return h1.time < h2.time;});

  /*
    Since a lot of this code would be duplicated, we'll also calculate the q50/n50 variables whilst we're at it.
   */

  std::deque<Hit> q50n50_window = {};
  size_t last_hit_idx = 0;
  
  // prepopulate q50n50_window with first 50ns worth of hits
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    if (tof_sub_hits.at(i).time - tof_sub_hits.front().time < 50){
      q50n50_window.push_back(tof_sub_hits.at(i));
    } else {
      last_hit_idx = i - 1;
      break;
    }
  }

  double q50n50_ratio = 0;
  
  while (last_hit_idx != tof_sub_hits.size()){

    const double current_q50 = std::accumulate(q50n50_window.begin(), q50n50_window.end(), 0, [](double a, const Hit& h){return h.charge;});
    if (q50n50_ratio > current_q50 / q50n50_window.size()){
      q50n50_ratio = current_q50 / q50n50_window.size();
    }
    
    // add the next subsequent hit not already the q50n50_window                                                                                                                                                    
    const double dt_to_next_hit = tof_sub_hits.at(last_hit_idx + 1).time - q50n50_window.back().time;
    q50n50_window.push_back(tof_sub_hits.at(last_hit_idx + 1));
    
    const double current_first_hit = q50n50_window.front().time;
    
    //remove all hits that are within dt from the front                                                                                                                                                             
    for (auto hit_it = q50n50_window.begin(); hit_it != q50n50_window.end(); ++hit_it){
      if (hit_it->time - current_first_hit > dt_to_next_hit){
        q50n50_window.erase(q50n50_window.begin(), hit_it - 1);
        break;
      }
    }
    
    ++last_hit_idx;
    
  }

  // back to calculating the preactivity...
  
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    for (size_t j = i; j < tof_sub_hits.size(); ++j){
      tof_sub_hits.at(i).goodness += 2 * CalculateGoodness(tof_sub_hits.at(i).time, tof_sub_hits.at(j).time);
    }
  }
  
  const auto max_it = std::max_element(tof_sub_hits.begin(), tof_sub_hits.end(), [](const Hit& h1, const Hit& h2){return (h1.goodness < h2.goodness);});
  const double max_goodness = max_it->goodness;

  //erase hits that fall below goodness threshold
  tof_sub_hits.erase(std::remove_if(tof_sub_hits.begin(), tof_sub_hits.end(),
				    [this, max_goodness](const Hit& h){
				      return (h.goodness < dark_threshold + fraction * max_goodness * exp(-h.time / 60));
				    }), tof_sub_hits.end());
  
  std::deque<double> preact_window = {};
  last_hit_idx = 0;
  
  int max_pre = 0;
  int max_pregate = 0;

  // prepopulate preact_window with first 15ns worth of hits
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    if (tof_sub_hits.at(i).time - tof_sub_hits.front().time < 15){
      preact_window.push_back(tof_sub_hits.at(i).time);
    } else {
      last_hit_idx = i - 1;
      break;
    }
  }
  
  while ((preact_window.back() < -12) && (last_hit_idx != tof_sub_hits.size() - 1)){

    if (preact_window.size() > static_cast<size_t>(max_pre)){
      max_pre = preact_window.size();
    }

    if ((preact_window.size() > static_cast<size_t>(max_pregate)) && ((preact_window.front() >= lowest_in_gate_time))){
      max_pregate = preact_window.size();
    }

    // add the next subsequent hit not already the preact_window
    const double dt_to_next_hit = tof_sub_hits.at(last_hit_idx + 1).time - preact_window.back();
    preact_window.push_back(tof_sub_hits.at(last_hit_idx + 1).time);

    const double current_first_hit = preact_window.front();

    //remove all hits that are within dt from the front
    for (auto hit_it = preact_window.begin(); hit_it != preact_window.end(); ++hit_it){
      if (*hit_it - current_first_hit > dt_to_next_hit){
	preact_window.erase(preact_window.begin(), hit_it - 1);
	break;
      }
    }
    
    ++last_hit_idx;
  }

  m_data->CStore.Set("q50n50_ratio", q50n50_ratio);
  m_data->CStore.Set("max_pre", max_pre);
  m_data->CStore.Set("max_pregate", max_pregate);

  return true;
}

bool CalculatePreactivityObservables::Finalise(){

  return true;
}

double CalculatePreactivityObservables::TimeOfFlight(const float* x, const float* y) const {
  double dist = 0;
  for (int i = 0; i < 3; ++i){
    dist += pow(x[i] - y[i], 2);
  }
  return (dist / SOL_IN_CM_PER_NS_IN_WATER); // speed of light in cm/ns
};  

bool CalculatePreactivityObservables::CalculateGoodness(const double& t1, const double& t2) const {
  const double dt2 = pow(0.2*(t1-t2), 2);
  return dt2 < 25 ? exp(-0.5 * dt2) : 0;
};
