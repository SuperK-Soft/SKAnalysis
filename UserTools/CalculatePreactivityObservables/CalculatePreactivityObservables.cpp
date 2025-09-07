#include "CalculatePreactivityObservables.h"

#include <algorithm>
#include <cmath>

#include "MTreeReader.h"
#include "TableReader.h"
#include "TableEntry.h"
#include "Constants.h"
#include "TGraph.h"
#include "TCanvas.h"

//#define PREACTIVITY_DEBUG   << to make histograms for every event

CalculatePreactivityObservables::CalculatePreactivityObservables():Tool(){}

bool CalculatePreactivityObservables::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  m_variables.Get("dark_threshold", dark_threshold);
  m_variables.Get("fraction", fraction);
  
  // FIXME make configurable
  q50n50_window_size = 50; //think this is in ns
  preact_window_size = 15;
  preact_window_cutoff = 12;

  GetTreeReader();
  
  connection_table = m_data->GetConnectionTable();
  
#ifdef PREACTIVITY_DEBUG
  h_maxpre = TH1F("h_maxpre","h_maxpre",100,0,0);
  h_q50n50 = TH1F("h_q50n50","h_q50n50",100,0,0);
  h_maxpregate = TH1F("h_maxpregate","h_maxpregate",100,0,0);
  h_bsvertex_t = TH1F("h_bsvertex_t","h_bsvertex_t",100,0,0);
#endif
  
  return true;
}

bool CalculatePreactivityObservables::Execute(){

  /*
    - get array of hit times
    - from first hit count how many subsequent hits are less than 15ns away
    - construct a vector from these hits
    - calculate goodness values and remove the hits that don't make goodness cut
    - * maxpre = vector.size()
    - if first hit in vector is within 1.3us from main trigger, maxpregate = vector.size()
    - find the next hit in the total readout
    - calculate goodness of this hit
    - if fails goodness test, move to next hit and try again - repeat until one passes
    - find the time, dt, between the last hit in the window and the new hit
    - add the new hit to end of the vector
    - ** remove the hits from start that are within dt from the first hit.
    - loop from * to ** until the last hit in the vector is 12ns before the main trigger
  */

  // used in calculating max_pregate
  double lowest_in_gate_time = 9999999;
  
  //std::vector<Hit> tof_sub_hits = std::vector<Hit>(sktqz_.nqiskz, {0,0,0});
  Log(m_unique_name+": total number of hits: "+std::to_string(sktqz_.nqiskz),v_debug,m_verbose);
  std::vector<Hit> tof_sub_hits;
  
  LoweInfo* lowe_ptr = nullptr;
  bool ok = LOWE_tree_reader->Get("LOWE", lowe_ptr);
  if (!ok || lowe_ptr == nullptr){
    throw std::runtime_error("couldn't get lowe branch");
  }
  TQReal* TQREAL = nullptr;
  ok = LOWE_tree_reader->Get("TQREAL", TQREAL);
  if (!ok || TQREAL == nullptr){
    throw std::runtime_error("couldn't get TQREAL branch");
  }
  
#ifdef PREACTIVITY_DEBUG
  TCanvas cpre;
  std::string titlet="h_times_"+std::to_string(skhead_.nevsk);
  double mint=*std::min_element(TQREAL->T,TQREAL->T+sktqz_.nqiskz);
  double maxt=*std::max_element(TQREAL->T,TQREAL->T+sktqz_.nqiskz);
  TH1F h_times(titlet.c_str(),titlet.c_str(),200,-50E3,50e3);
  std::string titlett="h_tofs_"+std::to_string(skhead_.nevsk);
  TH1F h_tofs(titlett.c_str(),titlett.c_str(),200,-50E3,50E3);
  std::vector<std::pair<float,float>> times_and_corrs;
  std::pair<float,int> min_time{999,-1};
  std::pair<float,int> min_time_corr{999,-1};
#endif
  
  // fill hits into vector, doing time of flight subtraction as we go.
  for (int hit_idx = 0; hit_idx < sktqz_.nqiskz; ++hit_idx){
    // get the cable numbers from sktqz_.icabiz as usual:
    const int cable_number = sktqz_.icabiz[hit_idx]; // n.b. no need for bitmask here
    if(cable_number>MAXPM || cable_number<=0) continue;
    // skip out-of-gate hits (i.e. not within trigger window at all)
    if((sktqz_.ihtiflz[hit_idx] & 0x02)==0) continue;
    
    //float raw_time = sktqz_.tiskz[hit_idx]; // something in muechk seems to change tiskz - perhaps it's internally calling set_timing_gate_
    float raw_time = TQREAL->T[hit_idx];      // get the times from TQREAL T branch directly (FIXME maybe call set_timing_gate(0) for others)
    
#ifdef PREACTIVITY_DEBUG
    if(raw_time<min_time.first){ min_time.first=raw_time; min_time.second=hit_idx; }
#endif
    
    // this Tool is to look for *pre-activity* - we probably don't need to process hits way after the trigger (e.g. AFT hits)
    if(raw_time>50E3) continue;
    
    float pmt_loc[3] = {};
    connection_table->GetTubePosition(cable_number, pmt_loc);
    
    double tof = TimeOfFlight(lowe_ptr->bsvertex, pmt_loc);
    const double new_time = raw_time - lowe_ptr->bsvertex[3] - tof;
    if (((sktqz_.ihtiflz[hit_idx] & 0x01)==1) && (new_time < lowest_in_gate_time)){
      lowest_in_gate_time = new_time; // 'in-gate' here refers to in 1.3us window
    }
    tof_sub_hits.emplace_back(new_time, 0, TQREAL->Q[hit_idx]); //calculate goodness in the next loop
    
#ifdef PREACTIVITY_DEBUG
    h_times.Fill(raw_time);
    h_tofs.Fill(new_time);
    times_and_corrs.push_back(std::pair<float,float>(raw_time,lowe_ptr->bsvertex[3]+tof));
    if(new_time<min_time_corr.first){ min_time_corr.first=new_time; min_time_corr.second=hit_idx; }
#endif
    
    /*
    if (hit_idx < 5){
      std::cout << "hit_idx: " << hit_idx << " on PMT "<< cable_number << "(cf MAXPM: "<<MAXPM<<")"
                <<" has location: (" << pmt_loc[0] << ", " << pmt_loc[1] << ", " << pmt_loc[2] << ")" << std::endl //cm
                <<" bonsai vtx: (" << lowe_ptr->bsvertex[0] << ", " << lowe_ptr->bsvertex[1] << ", " << lowe_ptr->bsvertex[2] << ")" << std::endl; //cm
      //std::cout << "hit time: " << sktqz_.tiskz[hit_idx] << std::endl; // nsec
      //std::cout << "charge: " << sktqz_.qiskz[hit_idx] << std::endl;
      //std::cout << "bsvertex[3]: " << lowe_ptr->bsvertex[3] << std::endl;  // index [0-2] is cm, index [3] is ns
      std::cout << "TimeOfFlight(lowe_ptr->bsvertex, pmt_loc): " << TimeOfFlight(lowe_ptr->bsvertex, pmt_loc) << std::endl;
      std::cout << "new time: " << new_time << std::endl;
    }
    */
  }
  
#ifdef PREACTIVITY_DEBUG
  h_times.Draw("");
  cpre.Modified();
  std::string h_times_fname=titlet+".png";
  cpre.SaveAs(h_times_fname.c_str());
  h_tofs.Draw("");
  cpre.Modified();
  std::string h_tofs_fname=titlett+".png";
  cpre.SaveAs(h_tofs_fname.c_str());
#endif
  
  /*
  std::cout << "tof_sub_hits.size(): (before goodness check) " << tof_sub_hits.size() << std::endl;
  std::cout << "lowest_in_gate_time was: " << lowest_in_gate_time << std::endl;
  std::cout<<"min ID hit time: hit "<<min_time.second<<" at "<<min_time.first<<std::endl;
  std::cout<<"min corr time: "<<min_time_corr.second<<" at "<<min_time_corr.first<<std::endl;
  */
  
  // print a few hits before the sorting
  /*
  std::cout << "some of the " << tof_sub_hits.size() << " hits before sort:" << std::endl;
  for (size_t hit_idx = 0; hit_idx < tof_sub_hits.size(); ++hit_idx){
    if (hit_idx < 10){ printf("%.3f\n",tof_sub_hits.at(hit_idx).time); }
  }
  */

  //sort hits in time
  std::sort(tof_sub_hits.begin(), tof_sub_hits.end(), [](const Hit& h1, const Hit& h2){return h1.time < h2.time;});

  /*
  //print a few hits after sorting
  std::cout << "some of the " << tof_sub_hits.size() << " hits after sort:" << std::endl;
  for (size_t hit_idx = 0; hit_idx < tof_sub_hits.size(); ++hit_idx){
    if (hit_idx < 10){ printf("%.3f\n",tof_sub_hits.at(hit_idx).time); }
  }
  */
  
  /*
  // debug checks
  if(skhead_.nevsk==1129864270){
          for(int hit_idx=0,j=0; j<10; ++hit_idx){
                  if((sktqz_.icabiz[hit_idx])<=0 || (sktqz_.icabiz[hit_idx]>MAXPM)) continue;
                  if((sktqz_.ihtiflz[hit_idx] & 0x02)==0) continue;
                  if(sktqz_.tiskz[hit_idx]>-4000) continue;
                  ++j;
                  std::cout<<"event "<<skhead_.nevsk<<" hit "<<hit_idx<<" at time: "<<sktqz_.tiskz[hit_idx]<<" ("<<TQREAL->T.at(hit_idx)<<")"
                           <<" PMT "<<sktqz_.icabiz[hit_idx]<<" ("<<(TQREAL->cables.at(hit_idx) & 0xFFFF)<<")"
                           <<", flags "<<sktqz_.ihtiflz[hit_idx]<<" ("<<(TQREAL->cables.at(hit_idx) >> 16)<<")"
                           <<", <MAXPM: "<<(sktqz_.icabiz[hit_idx]<MAXPM && sktqz_.icabiz[hit_idx]>0)
                           <<", in-gate: "<<(sktqz_.ihtiflz[hit_idx] & 0x02)<<std::endl;
                  
          }
          for(int i=0, j=0; j<10; ++i){
                  int hit_idx=TQREAL->T.size()-1-i;
                  if((sktqz_.icabiz[hit_idx])<=0 || (sktqz_.icabiz[hit_idx]>MAXPM)) continue;
                  if((sktqz_.ihtiflz[hit_idx] & 0x02)==0) continue;
                  ++j;
                  std::cout<<"event "<<skhead_.nevsk<<" hit "<<hit_idx
                           <<" PMT "<<sktqz_.icabiz[hit_idx]<<" ("<<(TQREAL->cables.at(hit_idx) & 0xFFFF)<<")"
                           <<", flags "<<sktqz_.ihtiflz[hit_idx]<<" ("<<(TQREAL->cables.at(hit_idx) >> 16)<<")"
                           <<", <MAXPM: "<<(sktqz_.icabiz[hit_idx]<MAXPM && sktqz_.icabiz[hit_idx]>0)
                           <<", in-gate: "<<(sktqz_.ihtiflz[hit_idx] & 0x02)<<std::endl;
          }
          
          std::sort(times_and_corrs.begin(), times_and_corrs.end(),
                  [](const std::pair<float,float>& h1, const std::pair<float,float>& h2) {return h1.first < h2.first;}
          );
          for(int i=0, j=0; j<10; ++i){
                  if(times_and_corrs.at(i).first>-4000) continue;
                  ++j;
                  std::cout<<"hit "<<i<<" TOF,vtx T corrected time: "<<tof_sub_hits.at(i).time<<std::endl;
                  std::cout<<"raw time, correction: "<<times_and_corrs.at(i).first<<", "<<times_and_corrs.at(i).second<<std::endl;
          }
          m_data->vars.Set("StopLoop",1);
          return true;
  }
  */

  /*
    Since a lot of this code would be duplicated, we'll also calculate the q50/n50 variables whilst we're at it.
   */
  /////////////////
  //std::cout << "calculate max q50n50" << std::endl;
  //rolling window for the q50n50 calculation 
  std::deque<Hit> q50n50_window = {};
  size_t next_hit_idx = 0;
  
  // prepopulate q50n50_window with first 50ns worth of hits
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    if (tof_sub_hits.at(i).time - tof_sub_hits.front().time < q50n50_window_size){
      q50n50_window.push_back(tof_sub_hits.at(i));
    } else {
      next_hit_idx = i;
      break;
    }
  }
  //std::cout << "initial q50n50_window.size(): " << q50n50_window.size() << std::endl;
  //std::cout << "next_hit_idx: " << next_hit_idx << std::endl;
  int n50=0;
  double q50n50_ratio=0;
  //go through hits
  while (next_hit_idx != tof_sub_hits.size() - 1){
    
    //std::cout << "calculate q50n50 ratio" << std::endl;
    //const double current_q50 = std::accumulate(q50n50_window.begin(), q50n50_window.end(), 0, [](double a, const Hit& h){return a + h.charge;});
    if (q50n50_window.size() > n50){
      n50 = q50n50_window.size();
      //std::cout << "new n50:" << n50 << std::endl;
      //std::cout << "accumulate charge" << std::endl;
      double current_q50 = 0;
      for (const auto& h : q50n50_window){
        current_q50 += h.charge;
        //if(q50n50_ratio==0) std::cout<<" hit charge "<<h.charge<<std::endl;
      }
      //std::cout << "new q50: " << current_q50 << std::endl;
      q50n50_ratio = current_q50 / n50;
      
      /* debug
      std::cout<<"new q50n50: "<<q50n50_ratio<<" based on hits: "<<std::endl;
      current_q50=0;
      n50=0;
      for(const auto& h : q50n50_window){
              current_q50+=h.charge;
              ++n50;
              q50n50_ratio=current_q50/n50;
              printf("hit at time: %.2f, dt: %.2f, charge: %.2f, total: %.2f, n50: %d, ratio: %.2f\n", h.time, (h.time - q50n50_window.front().time), h.charge, current_q50, n50, q50n50_ratio);
      }
      */
    }
    //std::cout<<"q50n50 ratio: "<<q50n50_ratio<<std::endl;
    
    if(next_hit_idx == tof_sub_hits.size()) break; // no more hits to grab

    //std::cout << "shifting q50n50 window" << std::endl;
    // drop hits from the the front of the q50n50_window until we exceed the dt to the next hit
    // (otherwise we're just removing hits and q50 is just going to be falling)
    double dt_to_next_hit = tof_sub_hits.at(next_hit_idx).time - q50n50_window.back().time;
    //std::cout<<"dt to next hit "<<dt_to_next_hit<<" vs q50n50 window size: "<<q50n50_window_size<<std::endl;
    // shortcut
    if(std::abs(dt_to_next_hit)>q50n50_window_size){
        q50n50_window.clear();
    } else {
        double current_first_hit = q50n50_window.front().time;
        while (q50n50_window.size() && (q50n50_window.front().time - current_first_hit) < std::abs(dt_to_next_hit)){
            //std::cout<<"pop with q50n50 size "<<q50n50_window.size()<<std::endl;
            q50n50_window.pop_front();
        }
    }
    if(q50n50_window.empty()){
      q50n50_window.push_back(tof_sub_hits.at(next_hit_idx));
      ++next_hit_idx;
    }

    //std::cout << "add new hits to the window" << std::endl;
    //std::cout << "next_hit_idx: " << next_hit_idx << std::endl;
    //std::cout << "tof_sub_hits.size() " << tof_sub_hits.size() << std::endl;
    // add new hits until the newly truncated window is 50ns long again
    while(next_hit_idx < tof_sub_hits.size()){
      if(std::abs(tof_sub_hits.at(next_hit_idx).time - q50n50_window.front().time) < q50n50_window_size){
        q50n50_window.push_back(tof_sub_hits.at(next_hit_idx));
        ++next_hit_idx;
      } else {
        break;
      }
    }
  }
  //std::cout<<"final q50n50_ratio: "<<q50n50_ratio<<std::endl;
  //////////////////////


  // back to calculating the preactivity...
  //std::cout << "give each hit a goodness" << std::endl;
  // by looping over all pairs of hits, we give each hit a goodness
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    for (size_t j = i+1; j < tof_sub_hits.size(); ++j){
      double g = CalculateGoodness(tof_sub_hits.at(i).time, tof_sub_hits.at(j).time);
      tof_sub_hits.at(i).goodness += g;
      tof_sub_hits.at(j).goodness += g;
    }
  }

  //get max goodness
  const auto max_it = std::max_element(tof_sub_hits.begin(), tof_sub_hits.end(), [](const Hit& h1, const Hit& h2){return (h1.goodness < h2.goodness);});
  const double max_goodness = max_it->goodness;
  //std::cout<<"max goodness: "<<max_goodness<<std::endl;
  
  // print a few hits before the goodness cut
  /*
  std::cout<<"first 10 hits before goodness cut:"<<std::endl;
  for (size_t hit_idx = 0; hit_idx < tof_sub_hits.size(); ++hit_idx){
    if (hit_idx < 10){ printf("%.3f\n",tof_sub_hits.at(hit_idx).time); }
  }
  */
  
#ifdef PREACTIVITY_DEBUG
  //TH1D h_goodness("h_goodness","goodness dist",100,0,6);
  std::vector<double> gs, ts;
  for(auto&& h : tof_sub_hits){
          //h_goodness.Fill(h.goodness);
          gs.push_back(h.goodness);
          ts.push_back(h.time);
  }
  //h_goodness.SaveAs("h_goodness.root");
  TGraph g_goodness(tof_sub_hits.size(),gs.data(),ts.data());
  g_goodness.Draw("a*");
  std::string fname="g_goodness_"+std::to_string(skhead_.nevsk)+".png";
  cpre.SaveAs(fname.c_str());
#endif
  
  //std::cout << "erase hits below goodness threshold " << std::endl;
  //erase hits that fall below goodness threshold
  //std::cout<<"hits before goodness cut: "<<tof_sub_hits.size();
  tof_sub_hits.erase(std::remove_if(tof_sub_hits.begin(), tof_sub_hits.end(),
                                    [this, max_goodness](const Hit& h){
                                     return (h.goodness < dark_threshold) ||
                                            (h.time > 0 && h.goodness < (dark_threshold + fraction*max_goodness*exp(-h.time/60.)) );
                                    }), tof_sub_hits.end());
  //std::cout<<", after goodness cut: "<<tof_sub_hits.size()<<std::endl;
  
#ifdef PREACTIVITY_DEBUG
  gs.clear(); ts.clear();
  for(auto&& h : tof_sub_hits){
          gs.push_back(h.goodness);
          ts.push_back(h.time);
  }
  TGraph g_goodness_pass(tof_sub_hits.size(),gs.data(),ts.data());
  g_goodness_pass.Draw("a*");
  fname="g_goodness_pass_"+std::to_string(skhead_.nevsk)+".png";
  cpre.SaveAs(fname.c_str());
#endif
  
  /*
  std::cout<<"first 10 hits after goodness cut:"<<std::endl;
  for (size_t hit_idx = 0; hit_idx < tof_sub_hits.size(); ++hit_idx){
    if (hit_idx < 10){ printf("%.3f\n",tof_sub_hits.at(hit_idx).time); }
  }
  */
  
  //create preactivity window
  std::deque<double> preact_window = {};
  // TODO the deque is redundant: the original lecompte.F does it more efficiently, just tracking indices in tof_sub_hits.
  // replace instances of preact_window.push_back with increment an end index,
  // replace instances of preact_window.pop_front with increment a start index
  // and then replace references to preact_window entries with ones to tof_sub_hits.
  next_hit_idx = 0;
  
  size_t max_pre = 0;
  size_t max_pregate = 0;

  // prepopulate preact_window with first 15ns worth of hits
  //std::cout << "prepopulate preact_window with first "<<preact_window_size<<"ns worth of hits" << std::endl;
  for (size_t i = 0; i < tof_sub_hits.size(); ++i){
    if ( std::abs(tof_sub_hits.at(i).time - tof_sub_hits.front().time) < preact_window_size && tof_sub_hits.at(i).time<-preact_window_cutoff){
      preact_window.push_back(tof_sub_hits.at(i).time);
    } else {
      next_hit_idx = i;
      break;
    }
  }
  //std::cout<<"initial preact_window size: "<<preact_window.size()<<std::endl;
  
  if(preact_window.size()==0){
    // i think this is prooobably fine...
    Log(m_unique_name+": no hits passing goodness cut with t<preact_window_cutoff",v_message,m_verbose);
  }

  //std::cout << "go through the hits in preactivity window" << std::endl;
  //std::cout<<"last preact hit time "<<preact_window.back()<<" vs cutoff "<<-preact_window_cutoff
  //         <<", next hit_idx: "<<next_hit_idx<<" vs tof_sub_hits.size(): "<<tof_sub_hits.size()<<std::endl;
  
  if(preact_window.size()){
    while ((preact_window.back() < -preact_window_cutoff) && (next_hit_idx != tof_sub_hits.size() - 1)){

      //std::cout<<"preact_window.size() after adding new hits "<<preact_window.size()<<std::endl;
      
      //std::cout <<"get max_pre for this window iteration" << std::endl;
      if (preact_window.size() > max_pre){
        max_pre = preact_window.size();
      }

      //std::cout << " and the same for max_pregate" << std::endl;
      if ((preact_window.size() > max_pregate) && ((preact_window.front() >= lowest_in_gate_time))){
        max_pregate = preact_window.size();
      }

      if(next_hit_idx == tof_sub_hits.size()) break; // no more hits to grab
      
      //std::cout << "shifting preact_window" << std::endl;
      // drop hits from the the front of the preact_window_size until we exceed the dt to the next hit
      // (otherwise we're just removing hits and max_pre/max_pregate are just going to be falling)
      double dt_to_next_hit = tof_sub_hits.at(next_hit_idx).time - preact_window.back();
      //std::cout<<"dt to next hit: "<<dt_to_next_hit<<std::endl;
      // shortcut
      if(std::abs(dt_to_next_hit)>preact_window_size){
          preact_window.clear();
      } else {
          double current_first_hit = preact_window.front();
          while (preact_window.size() && (preact_window.front() - current_first_hit) < std::abs(dt_to_next_hit)){
              preact_window.pop_front();
          }
      }
      if(preact_window.empty()){
        preact_window.push_back(tof_sub_hits.at(next_hit_idx).time);
        ++next_hit_idx;
      }
      //std::cout<<"preact_window.size() before adding new hits "<<preact_window.size()<<std::endl;

      //std::cout << " add new hits until the newly truncated window is 12ns long again" << std::endl;
      while (next_hit_idx < tof_sub_hits.size()){
        //std::cout<<"next hit at "<<tof_sub_hits.at(next_hit_idx).time<<" vs current preact first hit "<<preact_window.front()<<std::endl;
        if (std::abs(tof_sub_hits.at(next_hit_idx).time - preact_window.front()) < preact_window_size){
          preact_window.push_back(tof_sub_hits.at(next_hit_idx).time);
          //std::cout<<"added"<<std::endl;
          ++next_hit_idx;
        } else {
          break;
        } 
      }
      
    }
  }

  //std::cout << "max_pre: " << max_pre << std::endl;
  //std::cout << "max_pregate: " << max_pregate << std::endl;
#ifdef PREACTIVITY_DEBUG
  h_maxpre.Fill(max_pre);
  h_maxpregate.Fill(max_pregate);
  h_bsvertex_t.Fill(lowe_ptr->bsvertex[3]);
  h_q50n50.Fill(q50n50_ratio);
#endif
  
  m_data->CStore.Set("q50n50_ratio", q50n50_ratio);
  m_data->CStore.Set("max_pre", max_pre);
  m_data->CStore.Set("max_pregate", max_pregate);

  return true;
}

bool CalculatePreactivityObservables::Finalise(){
  
#ifdef PREACTIVITY_DEBUG
  h_maxpre.SaveAs("h_maxpre.root");
  h_maxpregate.SaveAs("h_maxpregate.root");
  h_q50n50.SaveAs("h_q50n50.root");
  TCanvas cpre("cpre","cpre");
  h_bsvertex_t.Draw();
  cpre.SaveAs("h_bsvertext.png");
#endif
  
  return true;
}

double CalculatePreactivityObservables::TimeOfFlight(const float* x, const float* y) const {
  double dist = 0;
  for (int i = 0; i < 3; ++i){
    dist += pow(x[i] - y[i], 2);
  }
  return (sqrt(dist) / SOL_IN_CM_PER_NS_IN_WATER); // speed of light in cm/ns
};  

double CalculatePreactivityObservables::CalculateGoodness(const double& t1, const double& t2) const {
  const double dt2 = pow(t2-t1, 2.)/25.;
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
