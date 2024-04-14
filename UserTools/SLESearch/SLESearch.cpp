#include "SLESearch.h"

#include <algorithm>

#include "MTreeReader.h"
#include "fortran_routines.h"
#include "softtrg_tbl.h"
#include "Constants.h"
#include "TableReader.h"
#include "TableEntry.h"

#include <algorithm>

#include "TH1D.h"

#include <bitset>

SLESearch::SLESearch():Tool(){}

bool SLESearch::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  if(!m_variables.Get("include_offset", include_offset)) include_offset=false;

  connection_table = m_data->GetConnectionTable();

  hit_times_plot = TH1D("hit_times", "hit times;time us; number of hits", 1000, 0, 0);
  
  return true;
}


bool SLESearch::Execute(){

  std::bitset<32> trigger_id{skhead_.idtgsk};

  /* 
     SLE search algorithm:
     1. Construct vector of hits from event
     2. sort hits in time order
     3. Starting at the first hit, construct a deque of the first [window size] hits
     4. is number of hits in the window greater than SLE threshold?
     5. Yes to 4? Use "algorithm" to determine the timing of the SLE trigger (eg median or quantile calculation)
     6. move start of trigger search to time of SLE trigger + trigger deadtime
     7. Store new SLE time
     8. Drop the first hit from the beginning of the window
     9. Add on hits not already in the deque until the window size (last hit time - first hit time) is back up to [window_length]
     10. Loop steps 4 to 9 until all last hit in window == last hit in event

  */
  std::vector<double> SLE_times;
  double current_SLE_time = 0;
      
  // 1. Construct vector of hits from event
  std::cout << "constructing vector of hits from event" << std::endl;
  std::vector<double> hits = std::vector<double>(sktqz_.tiskz, sktqz_.tiskz + sktqz_.nqiskz);
  std::vector<double> first_hit_times;

  event_hit_times_plot = TH1D("event_hit_times_plot", "event_hit_times_plot", 200, hits.front(), hits.back());
  
  // 2. sort hits in time order
  std::cout << "sorting hits into time order" << std::endl;
  std::sort(hits.begin(), hits.end());

  for (const auto& time : hits){event_hit_times_plot.Fill(time);}
  
  // 3. Starting from the first hit, construct a deque of the first [window size] hits
  const double window_size = 200;

  std::deque<double> window = {};

  std::cout << "hits.size(): " << hits.size() << std::endl;
  std::cout << "hits.front(): " << hits.front() << std::endl;
  std::cout << "hits.back(): " << hits.back() << std::endl;
  std::cout<<"event window length: "<<(hits.back()-hits.front())<<" ns"<<std::endl;
  std::cout<<"window length: "<<window_size<<" ns"<<std::endl;
  size_t next_hit_idx = 0;
  for (size_t i = 0; i < hits.size(); ++i){
    if (hits.at(i) - hits.front() < window_size){
      window.push_back(hits.at(i));
    } else {
      next_hit_idx = i;
      break;
    }
  }

  const int SLE_threshold = 34;
  const double SLE_deadtime = 0; // assume no deadtime?
  const double SLE_readout_length = 1501.56; //nope

  std::cout<<"initial hits window had "<<window.size()<<" hits vs threshold "<<SLE_threshold<<std::endl;

  //  while (window.back() != hits.back()){
  std::cout<<"window front index: "<<next_hit_idx<<std::endl;
  while ( (next_hit_idx != hits.size()) && (window.back() != hits.back()) ){
    //  std::cout << "next_hit_idx: " << next_hit_idx << std::endl;

    // 4. is number of hits in the window greater than SLE threshold?
    //const int SLE_threshold = GetTriggerThreshold(TriggerType::SLE);

    if (window.size() > SLE_threshold){
      //std::cout << "number of hits in window exceeds SLE threshold!" << std::endl;
      std::cout << " threshold exceeded with window.size(): " << window.size() << std::endl;
      //std::cout << "hits size: " << hits.size() << " next hit idx: " << next_hit_idx << std::endl;
      //std::cout << "before adding deadtime" << std::endl;
      //      std::cout<< "at time: " << window.front() << ", window.back() at: " << window.back() << std::endl;

      std::cout << "window.front() at in ticks: " << (window.front() * COUNT_PER_NSEC) + skheadqb_.it0sk;
      std::cout << "window.back() at in ticks: " << (window.back() * COUNT_PER_NSEC) + skheadqb_.it0sk;
      
      first_hit_times.push_back(window.front() + (include_offset ? SLE_t0_offset : 0));

      // 5. Yes to 4? Use "algorithm" to determine the timing of the SLE trigger (eg median or quantile calculation)
      // since this may not be accurate for where the hit cluster is, maybe also try to better estimate it
      current_SLE_time = window.at(static_cast<size_t>(window.size() / 2)); //shit - fix this
      //std::cout << "sum: " << std::accumulate(window.begin(), window.end(), 0.0) << std::endl;
      // double this_SLE_time = 0;
      //      for (int i = ; 
      // 7. Store new SLE time
      std::cout << "estimated SLE time: " << current_SLE_time << std::endl;
      if (current_SLE_time > 0){
	hit_times_plot.Fill(current_SLE_time);
	SLE_times.push_back(current_SLE_time);
      } else {
	std::cout << "SLE time was before primary trigger - so we ignore" << std::endl; 
      }
	
      std::cout << std::endl << std::endl;

      /*std::cout << "some hits: " << std::endl;
	for (int i = 0; i < window.size(); ++i){
	if (i % 50 == 0){std::cout << window.at(i) << std::endl;}
	}*/
      
      // 6. move start of trigger search to start of trigger window + trigger window length + trigger deadtime

      // 6a. add on SLE_deadtime worth of hits
      double new_start_time = window.front() + SLE_readout_length + SLE_deadtime;
      //std::cout << "last_hit_in_window_time: " << last_hit_in_window_time << std::endl;
      //std::cout << "last hit in event time " << hits.back() << std::endl;
      std::cout<<"scanning forward by readout lenght "<<SLE_readout_length<<" + deadtime "<<SLE_deadtime
	       <<" from end of window "<<new_start_time<<" ns -> new start time "<<new_start_time << " ns"<<std::endl;
      for (size_t i = next_hit_idx; i < hits.size(); ++i){
	++next_hit_idx; // next hit to consider
	//if (hits.at(i) - new_start_time < SLE_deadtime){
        if(hits.at(i) < new_start_time){
	  //if (i % 1000 == 0){ std::cout << "hits.at(i): " << hits.at(i) << std::endl;}
	  //window.push_back(hits.at(i));
	} else {
	  break;
	}
      }
      //      std::cout << "added SLE deadtime: " << SLE_deadtime << std::endl;
      //     std::cout << "after adding on dead time" << std::endl;
      //      std::cout<< "window.front(): " << window.front() << ", window.back(): " << window.back() << std::endl;
      
      // 6b. remove hits from the front of the window so that the window returns back to its proper size
      // empty the window
      window.clear();
      // re-fill the new trigger window
      for (size_t i = next_hit_idx; i < hits.size(); ++i){
	if (hits.at(i) - hits.at(next_hit_idx) < window_size){
	  window.push_back(hits.at(i));
	} else {
	  next_hit_idx = i;
	  break;
	}
      }
      /*
	last_hit_in_window_time = window.back();
	window.erase(std::remove_if(window.begin(), window.end(),
	[ last_hit_in_window_time, window_size](const double& hit){
	return last_hit_in_window_time - hit > window_size;
	}), window.end());
	std::cout << "after removing hits so the window has its original size" << std::endl;
	std::cout<< "window.front(): " << window.front() << ", window.back(): " << window.back() << std::endl;
      */

      // // 6b*. construct new window, with same length, with start time equal to old first time + deadtime
      // const double last_hit_in_window_time = window.back();
      // for (size_t i = next_hit_idx; i < hits.size()){
      // 	if (hits.at(i) - last_hit_in_window_time < SLE_deadtime){
      // 	  ++next_hit_idx;
      // 	} else {
      // 	  next_hit_idx = i;
      // 	  break;
      // 	}
      // }
      
      // window.clear();
      // for (size_t i = next_hit_idx; i < hits.size(); ++i){
      // 	if (hits.at(i) - hits.front() < window_size){
      // 	  window.push_back(hits.at(i));
      // 	} else {
      // 	  next_hit_idx = i;
      // 	  break;
      // 	}
      // }
            
    }

    // 8. Drop the first hit from the beginning of the window
    else {
      if (window.size() > 0){window.pop_front();}
      // 9. Add on hits not already in the deque until the window size (last hit time - first hit time) is back up to [window_size]
      //std::cout << "e0" << std::endl;
      for (size_t i = next_hit_idx; i < hits.size(); ++i){
        //std::cout << "i: " << i << "   next hit idx: " << next_hit_idx << std::endl;
        if ((window.size() == 0) || (window.back() - window.front() < window_size)){
	  ++next_hit_idx;
	  window.push_back(hits.at(i));
        } else {
	  break;
        }
      }
    }
  }

  //m_data->CStore.Set("SLE_times", first_hit_times);
  m_data->CStore.Set("SLE_times", SLE_times);

  int N_SLE = SLE_times.size();
  std::cout << "found " << N_SLE << " SLE times, they are:" << std::endl;

  for (const auto& time : SLE_times){std::cout << "time in ns: " << time << std::endl;}
  std::cout << "then" << std::endl;
  for (const auto& time : SLE_times){std::cout << "time in ticks: " << (time * COUNT_PER_NSEC) + skheadqb_.it0sk << std::endl;}


  // we do this again in Pre recon cuts
  m_data->CStore.Set("N_SLE", N_SLE); //need this for the subtoolchain

  //  hit_times_plot.SaveAs("hit_times.root");
  event_hit_times_plot.SaveAs("event_hit_times_plot.root");
  
  return true;
}


bool SLESearch::Finalise(){
  
  return true;
}
