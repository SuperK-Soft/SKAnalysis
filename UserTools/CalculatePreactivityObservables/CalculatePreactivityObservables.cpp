#include "CalculatePreactivityObservables.h"

#include <algorithm>

#include "MTreeReader.h"
#include "TableReader.h"
#include "TableEntry.h"
#include "Constants.h"

CalculatePreactivityObservables::CalculatePreactivityObservables():Tool(){}

bool CalculatePreactivityObservables::Initialise(std::string configfile, DataModel &data){

  std::cout << "HELLLLLOOOO!" << std::endl;

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  m_variables.Get("dark_threshold", dark_threshold);
  m_variables.Get("fraction", fraction);

  GetTreeReader();
  
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

  // set some constants
  double lowest_in_gate_time = 9999999;
  const double q50n50_window_size = 50; //think this is in ns
  const double preact_window_size = 15; 

  //std::vector<Hit> tof_sub_hits = std::vector<Hit>(sktqz_.nqiskz, {0,0,0});
  std::cout << "number of hits: " << sktqz_.nqiskz << std::endl;
  std::vector<Hit> tof_sub_hits;

  // fill hits into vector, doing time of flight subtraction as we go.
  for (int hit_idx = 0; hit_idx < sktqz_.nqiskz; ++hit_idx){
    // we need to seperate out the cable numbers from sktqz_.icabiz as usual:
    const int cable_number = sktqz_.icabiz[hit_idx];
    float pmt_loc[3] = {};
    connection_table->GetTubePosition(cable_number, pmt_loc);
    
    const double new_time = sktqz_.tiskz[hit_idx] - skroot_lowe_.bsvertex[3] - TimeOfFlight(skroot_lowe_.bsvertex, pmt_loc);
    if (((sktqz_.ihtiflz[hit_idx] & 0x01)==1) && (new_time < lowest_in_gate_time)){
      lowest_in_gate_time = new_time;
    }
    tof_sub_hits.emplace_back(new_time, 0, sktqz_.qiskz[hit_idx]); //calculate goodness in the next loop

    //    if (hit_idx < 5){
      std::cout << "hit_idx: " << hit_idx << " has location: (" << pmt_loc[0] << ", " << pmt_loc[1] << ", " << pmt_loc[2] << ")" << std::endl; //cm
      std::cout << "hit time: " << sktqz_.tiskz[hit_idx] << std::endl; // nsec
      std::cout << "charge: " << sktqz_.qiskz[hit_idx] << std::endl;
      std::cout << "skroot_lowe_.bsvertex[3]: " << skroot_lowe_.bsvertex[3] << std::endl;  // index [0-2] is mm, index [3] is ns
      std::cout << "TimeOfFlight(skroot_lowe_.bsvertex, pmt_loc): " << TimeOfFlight(skroot_lowe_.bsvertex, pmt_loc) << std::endl;
      std::cout << "new time: " << new_time << std::endl;
      ///    }

  }
  
  std::cout << "lowest_in_gate_time was: " << lowest_in_gate_time << std::endl;

  // print a few hits before the sorting
  std::cout << "some of the " << tof_sub_hits.size() << " hits before sort:" << std::endl;
  for (size_t hit_idx = 0; hit_idx < tof_sub_hits.size(); ++hit_idx){
    if (hit_idx < 10){std::cout << tof_sub_hits.at(hit_idx).time << std::endl;}
  }

  //sort hits in time
  std::sort(tof_sub_hits.begin(), tof_sub_hits.end(), [](const Hit& h1, const Hit& h2){return h1.time < h2.time;});

  //print a few hits after sorting 
  std::cout << "some of the " << tof_sub_hits.size() << " hits after sort:" << std::endl;
  for (size_t hit_idx = 0; hit_idx < tof_sub_hits.size(); ++hit_idx){
    if (hit_idx < 10){std::cout << tof_sub_hits.at(hit_idx).time << std::endl;}
  }

  /*
    Since a lot of this code would be duplicated, we'll also calculate the q50/n50 variables whilst we're at it.
   */

  //rolling window for the q50n50 calculation 
  std::deque<Hit> q50n50_window = {};
  size_t next_hit_idx = 0;
  
  // prepopulate q50n50_window with first 50ns worth of hits
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    if (i < 5){
      std::cout << "tof_sub_hits.at(i).time: " << tof_sub_hits.at(i).time << std::endl;
      std::cout << "tof_sub_hits.front().time: " << tof_sub_hits.front().time << std::endl;
      std::cout << "tof_sub_hits.at(i).time - tof_sub_hits.front().time: " << tof_sub_hits.at(i).time - tof_sub_hits.front().time << std::endl;
      std::cout << "q50n50_window_size: "  << q50n50_window_size << std::endl;
    }
    if (tof_sub_hits.at(i).time - tof_sub_hits.front().time < q50n50_window_size){
      q50n50_window.push_back(tof_sub_hits.at(i));
    } else {
      next_hit_idx = i;
      break;
    }
  }

  std::cout << "q50n50_window.size(): " << q50n50_window.size() << std::endl;
  std::cout << "next_hit_idx: " << next_hit_idx << std::endl;
  
  double q50n50_ratio = 0;

  std::cout << "tof_sub_hits.size(): " << tof_sub_hits.size() << std::endl;
  
  //go through hits
  std::cout << "go through hits" << std::endl;
  while (next_hit_idx != tof_sub_hits.size() - 1){

    
    std::cout << "accumulate charge" << std::endl;
    std::cout << "q50n50_window.size():" << q50n50_window.size() << std::endl;
    double current_q50 = 0;
    for (const auto& h : q50n50_window){
      current_q50 += h.charge;
    }
    std::cout << "current_q50: " << current_q50 << std::endl;

    std::cout << "calculate q50n50 ratio" << std::endl;
    //const double current_q50 = std::accumulate(q50n50_window.begin(), q50n50_window.end(), 0, [](double a, const Hit& h){return a + h.charge;});
    if (q50n50_ratio > current_q50 / q50n50_window.size()){
      q50n50_ratio = current_q50 / q50n50_window.size();
    }

    std::cout << "drop the first hit" << std::endl;
    // drop the first hit in the q50n50_window
    if (!q50n50_window.empty()){
      q50n50_window.pop_front();
    } else {
      ++next_hit_idx;
    }

    std::cout << "add new hits to the window" << std::endl;
    std::cout << "next_hit_idx: " << next_hit_idx << std::endl;
    std::cout << "tof_sub_hits.size() " << tof_sub_hits.size() << std::endl;
    // add new hits until the newly truncated window is 50ns long again
    while(next_hit_idx < tof_sub_hits.size()){
      if(tof_sub_hits.at(next_hit_idx).time - q50n50_window.front().time < q50n50_window_size){
	q50n50_window.push_back(tof_sub_hits.at(next_hit_idx));
	++next_hit_idx;
      } else {
	break;
      }
    }
  }


  std::cout << "give each hit a goodness" << std::endl;
  // back to calculating the preactivity...
  // by looping over all pairs of hits, we give each hit a goodness 
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    for (size_t j = i; j < tof_sub_hits.size(); ++j){
      tof_sub_hits.at(i).goodness += 2 * CalculateGoodness(tof_sub_hits.at(i).time, tof_sub_hits.at(j).time);
    }
  }

  //get max goodness
  std::cout << "calculate max goodness " << std::endl;
  const auto max_it = std::max_element(tof_sub_hits.begin(), tof_sub_hits.end(), [](const Hit& h1, const Hit& h2){return (h1.goodness < h2.goodness);});
  const double max_goodness = max_it->goodness;

  std::cout << "erase hits below goodness threshold " << std::endl;
  //erase hits that fall below goodness threshold
  tof_sub_hits.erase(std::remove_if(tof_sub_hits.begin(), tof_sub_hits.end(),
				    [this, max_goodness](const Hit& h){
				      return (h.goodness < dark_threshold + fraction * max_goodness * exp(-h.time / 60));
				    }), tof_sub_hits.end());

  //create preactivity window
  std::deque<double> preact_window = {};
  next_hit_idx = 0;
  
  int max_pre = 0;
  int max_pregate = 0;

  // prepopulate preact_window with first 15ns worth of hits
  std::cout << "prepopulate preact_window with first 15ns worth of hits" << std::endl;
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    if (tof_sub_hits.at(i).time - tof_sub_hits.front().time < preact_window_size){
      preact_window.push_back(tof_sub_hits.at(i).time);
    } else {
      next_hit_idx = i;
      break;
    }
  }

  std::cout << "go through the hits in preactivity window" << std::endl;
  while ((preact_window.back() < -12) && (next_hit_idx != tof_sub_hits.size() - 1)){

    std::cout <<"get max_pre for this window iteration" << std::endl;
    if (preact_window.size() > static_cast<size_t>(max_pre)){
      max_pre = preact_window.size();
    }

    std::cout << " and the same for max_pregate" << std::endl;
    if ((preact_window.size() > static_cast<size_t>(max_pregate)) && ((preact_window.front() >= lowest_in_gate_time))){
      max_pregate = preact_window.size();
    }

    std::cout << "drop the first hit in the preact_window" << std::endl;
    if (!preact_window.empty()){
      preact_window.pop_front();
    } else {
      ++next_hit_idx;
    }

    std::cout << " add new hits until the newly truncated window is 12ns long again" << std::endl;
    while (next_hit_idx < tof_sub_hits.size()){
      if (tof_sub_hits.at(next_hit_idx).time - preact_window.front() < preact_window_size){
	preact_window.push_back(tof_sub_hits.at(next_hit_idx).time);
	++next_hit_idx;
      } else {
	break;
      } 
    }
  }

  std::cout << "q50n50_ratio: " << q50n50_ratio << std::endl;
  std::cout << "max_pre: " << max_pre << std::endl;
  std::cout << "max_pregate: " << max_pregate << std::endl;
  
  m_data->CStore.Set("q50n50_ratio", q50n50_ratio);
  m_data->CStore.Set("max_pre", max_pre);
  m_data->CStore.Set("max_pregate", max_pregate);

  return true;
}

bool CalculatePreactivityObservables::Finalise(){

  return true;
}

double CalculatePreactivityObservables::TimeOfFlight(const float* x, const float* y) const {
  // we need to convert bsvertex from mm to cm 
  double dist = 0;
  for (int i = 0; i < 3; ++i){
    dist += pow((x[i]/1000) - y[i], 2);
  }
  return (dist / SOL_IN_CM_PER_NS_IN_WATER); // speed of light in cm/ns
};  

bool CalculatePreactivityObservables::CalculateGoodness(const double& t1, const double& t2) const {
  const double dt2 = pow(0.2*(t1-t2), 2);
  return dt2 < 25 ? exp(-0.5 * dt2) : 0;
};

void CalculatePreactivityObservables::GetTreeReader(){
  std::string tree_reader_str = "";
  m_variables.Get("reader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("CalculatePreactivityObservables::GetTreeReader(): - Failed to get treereader "+tree_reader_str+"!");
  }
  LOWE_tree_reader = m_data->Trees.at(tree_reader_str);
}
