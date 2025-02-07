#include "CopyHits.h"

#include <algorithm>

#include "tqrealroot.h"

#include "MTreeReader.h"

CopyHits::CopyHits():Tool(){}


bool CopyHits::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReader();
  
  return true;
}

bool CopyHits::Execute(){

  
  /* 
     To run bonsai we need the skt_ and skq_ commons to actually have hits.
     There is probably some deep clandestine reason why set_timing_gate() followed by skcread() is removing these hits
     Rather than finding this, we'll just do things quick and dirty - copying the hits from the TQREAL branch to the required commons.
  */

  // TQReal* tqreal_ptr = nullptr;
  // if (!tree_reader_ptr->Get("TQREAL", tqreal_ptr) || tqreal_ptr == nullptr){
  //   throw std::runtime_error("CopyHits::Execute: failed to get TQREAL branch");
  // }

  // TQReal* tqareal_ptr = nullptr;
  // if (!tree_reader_ptr->Get("TQAREAL", tqareal_ptr) || tqareal_ptr == nullptr){
  //   throw std::runtime_error("CopyHits::Execute: failed to get TQAREAL branch");
  // }
  
  // std::cout << "got TQREAL and TQAREAL branches, here are some hits:" << std::endl;
  //std::cout << "TQREAL::nhits = " << tqreal_ptr->nhits << std::endl;
  //for (int i = 0; i < std::min(tqreal_ptr->nhits, 5); ++i){
  //   std::cout << "TQREAL hit no: " << i
  // 	      << " has time " << tqreal_ptr->T.at(i)
  // 	      << " and charge " << tqreal_ptr->Q.at(i)
  // 	      << std::endl;
  // }
  // for (int i = 0; i < std::min(tqareal_ptr->nhits, 5); ++i){
  //   std::cout << "TQAREAL hit no: " << i
  // 	      << " has time " << tqareal_ptr->T.at(i)
  // 	      << " and charge " << tqareal_ptr->Q.at(i)
  // 	      << std::endl;
  // }

  /*
    I _think_ we need to remove the previous event's hits else the, for example, skt_.tisk array will have more hits than skt_.ntisk
  */

  for (int i = 0; i < MAXPM; ++i){
    m_data->skt_common_dupl.tisk[i] = 0;
    m_data->skq_common_dupl.qisk[i] = 0;
  }
  
  /* 
     Now let's copy the info over to the m_data->skt_common_dupl and m_data->skq_common_dupl commons
     
     TQReal::T() has indexing running from 0 to TQReal::nhits-1
     m_data->skt_common_dupl.tisk[] is indexed with the (PMT number - 1) given by skchnl_.ihcab

     For example, we need m_data->skt_common_dupl.tisk[skchnl_.ihcab[i] - 1] == TQReal::T(i) for i : 0 -> TQReal::T().size()-1

  */

  // let's zero/max things where appropriate
  m_data->skq_common_dupl.nqisk = 0;
  m_data->skq_common_dupl.qismsk = 0;
  m_data->skq_common_dupl.qimxsk = 0;
  m_data->skq_common_dupl.mxqisk = 0;
  m_data->skt_common_dupl.timnsk = 9999999;
  m_data->skt_common_dupl.timxsk = 0;

  // skqa_.nqask = 0;
  // skqa_.qasmsk = 0;
  // skta_.tamnsk = 9999999;
  // skta_.mxtask = 0;

  int just_in_window_hits_idx = 0;
  const unsigned int in_gate_mask = ~0 ^ 1; 

  std::cout << "m_data->sktqz_common_dupl.nqiskz: " << m_data->sktqz_common_dupl.nqiskz << std::endl;


  //need to offset the 1.3us window boundaryies with the time of the found SLE entry
  std::vector<double> SLE_times = {};                                                                                   
  if (!m_data->CStore.Get("SLE_times", SLE_times)){throw std::runtime_error("CompareCommons::Execute couldn't retrieve SLESearch subtrigger times");}

  std::cout << "SLE_times.at(0): " << SLE_times.at(0) << std::endl;
  std::cout << "skheadqb_.it0sk: " << skheadqb_.it0sk << std::endl;
							 
  double SLE_time_in_ticks = SLE_times.at(0) - static_cast<double>(skheadqb_.it0sk);
  std::cout << "SLE_time_in_ticks: " << SLE_time_in_ticks << std::endl;
  
  // ID:
  //for (int all_hits_idx = 0; all_hits_idx < tqreal_ptr->nhits; ++all_hits_idx){
  for (int all_hits_idx = 0; all_hits_idx < m_data->sktqz_common_dupl.nqiskz; ++all_hits_idx){
	
      // m_data->sktqz_common_dupl.ihtiflz[just_in_window_hits_idx] = tqreal_ptr->cables.at(all_hits_idx) >> 16;
    
    // std::cout << "m_data->sktqz_common_dupl.tiskz["<<all_hits_idx<<"] * COUNT_PER_NSEC: " << m_data->sktqz_common_dupl.tiskz[all_hits_idx] * COUNT_PER_NSEC << std::endl;
    // std::cout << "skheadqb_.it0xsk: " << skheadqb_.it0xsk << std::endl;
    // std::cout << "one_p_three_gate_pre_t0: " << one_p_three_gate_pre_t0 << std::endl;
    // std::cout << "one_p_three_gate_post_t0: " << one_p_three_gate_post_t0 << std::endl;
    // std::cout << "skheadqb_.it0xsk + one_p_three_gate_pre_t0: " << skheadqb_.it0xsk + one_p_three_gate_pre_t0 << std::endl;
    // std::cout << "skheadqb_.it0xsk + one_p_three_gate_post_t0: " << skheadqb_.it0xsk + one_p_three_gate_post_t0 << std::endl;

    m_data->sktqz_common_dupl.tiskz[all_hits_idx] -= (SLE_time_in_ticks / COUNT_PER_NSEC); 
    
    if (((m_data->sktqz_common_dupl.tiskz[all_hits_idx] * COUNT_PER_NSEC) < /*SLE_time_in_ticks +*/ one_p_three_gate_pre_t0) ||
        ((m_data->sktqz_common_dupl.tiskz[all_hits_idx] * COUNT_PER_NSEC) > /*SLE_time_in_ticks +*/ one_p_three_gate_post_t0)){
      m_data->sktqz_common_dupl.ihtiflz[all_hits_idx] = m_data->sktqz_common_dupl.ihtiflz[all_hits_idx] & in_gate_mask;
    } else {
      m_data->sktqz_common_dupl.ihtiflz[all_hits_idx] = m_data->sktqz_common_dupl.ihtiflz[all_hits_idx] | ~in_gate_mask;

      // we need to populate the m_data->skchnl_common_dupl.ichab (PMT cable numbers) and IHTIFLZ (the hit flags)
      // the upper 16 bits are the cable numbers

      //m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx] = (tqreal_ptr->cables.at(all_hits_idx) & 0x0000FFFF) - 1;
      m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx] = m_data->sktqz_common_dupl.icabiz[all_hits_idx];
      // the lower 16 are the hit flags
      
      // increment the number of ID hits:
      ++m_data->skq_common_dupl.nqisk;
      // and the total ID charge:
      //m_data->skq_common_dupl.qismsk += tqreal_ptr->Q.at(all_hits_idx);
      // raw charge is in 0-10 bits of m_data->skq_common_dupl.qismsk
      m_data->skq_common_dupl.qismsk += (m_data->sktqz_common_dupl.iqiskz[all_hits_idx] & 0x7FF);
      // max charge deposited on ID PMT:
      //if (tqreal_ptr->Q.at(all_hits_idx) > m_data->skq_common_dupl.qimxsk){
      if ((m_data->sktqz_common_dupl.iqiskz[all_hits_idx] & 0x7FF) > m_data->skq_common_dupl.qimxsk){
	//m_data->skq_common_dupl.qimxsk = tqreal_ptr->Q.at(all_hits_idx);
	m_data->skq_common_dupl.qimxsk = (m_data->sktqz_common_dupl.iqiskz[all_hits_idx] & 0x7FF);
	// the PMT number with that charge:
	m_data->skq_common_dupl.mxqisk = m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx];
      }
      // min hit time in ID:
      //if (m_data->skt_common_dupl.timnsk > tqreal_ptr->T.at(all_hits_idx)){
      if (m_data->skt_common_dupl.timnsk > m_data->sktqz_common_dupl.tiskz[all_hits_idx]){
	//m_data->skt_common_dupl.timnsk = tqreal_ptr->T.at(all_hits_idx);
	m_data->skt_common_dupl.timnsk = m_data->sktqz_common_dupl.tiskz[all_hits_idx];
	// the PMT number with that time:
	m_data->skt_common_dupl.mntisk = m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx];
      }
      // max hit time in ID:
      //if (m_data->skt_common_dupl.timxsk < tqreal_ptr->T.at(all_hits_idx)){
      if (m_data->skt_common_dupl.timxsk < m_data->sktqz_common_dupl.tiskz[all_hits_idx]){
	//m_data->skt_common_dupl.timxsk = tqreal_ptr->T.at(all_hits_idx);
	m_data->skt_common_dupl.timxsk = m_data->sktqz_common_dupl.tiskz[all_hits_idx];
	// the PMT number with that time:
	m_data->skt_common_dupl.mxtisk = m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx];
      }
      // then the actual hits in the ID:
      //m_data->skt_common_dupl.tisk[m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx] - 1] = tqreal_ptr->T.at(all_hits_idx);
      m_data->skt_common_dupl.tisk[m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx] - 1] = m_data->sktqz_common_dupl.tiskz[all_hits_idx];
      //m_data->skq_common_dupl.qisk[m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx] - 1] = tqreal_ptr->Q.at(all_hits_idx);
      m_data->skq_common_dupl.qisk[m_data->skchnl_common_dupl.ihcab[just_in_window_hits_idx] - 1] = m_data->sktqz_common_dupl.iqiskz[all_hits_idx] & 0x7FF;

      ++just_in_window_hits_idx;
    }
    
    
  }

  std::cout << "for our boy smy, may we deck him:" << std::endl;
  std::cout << "nqisk: " << m_data->skq_common_dupl.nqisk << std::endl;
  std::cout << "qismsk: " << m_data->skq_common_dupl.qismsk << std::endl;
  
  // BONSAI supposedly doesn't need the OD hits to reconstruct but there might be other sections of the reduction that do require them, so this is commented out for now.
  // FIXME: charges and times were not being read in properly - the 1.3us selection was cutting out everything
  
  // // OD:
  // for (int j = 0; j < tqareal_ptr->nhits; ++j){
  //   // first we need to populate the sktqaz_.ihacab (PMT cable numbers) and IHTFLZ (the hit flags)
  //   // the upper 16 bits are the cable numbers
  //   sktqaz_.ihacab[j] = (tqareal_ptr->cables.at(j) & 0x0000FFFF) - QB_OD_OFFSET;
  //   // the lower 16 are the hit flags
  //   sktqaz_.ihtflz[j] = tqareal_ptr->cables.at(j) >> 16;
  //   //if (!(sktqaz_.ihtflz[j] & 1)){std::cout << "reject hit j="<<j<<std::endl;continue;} // ignore hits not in 1.3us window around the primary trigger
  //   if((sktqaz_.ihtflz[j] & 0x01)==0){continue;}
  //   std::cout << "accept hit j="<<j<<std::endl;
  //   // number of OD hits:
  //   ++skqa_.nqask;
  //   // total OD charge:
  //   skqa_.qasmsk += tqareal_ptr->Q.at(j);
  //   // max charge deposited on OD PMT:
  //   if (tqareal_ptr->Q.at(j) > skqa_.qamxsk){
  //     skqa_.qamxsk = tqareal_ptr->Q.at(j);
  //     // the PMT number with that charge:
  //     skqa_.mxqask = sktqaz_.ihacab[j];
  //   }
  //   // min hit time in OD:
  //   if (skta_.tamnsk > tqareal_ptr->T.at(j)){
  //     skta_.tamnsk = tqareal_ptr->T.at(j);
  //     // the PMT number with that time:
  //     skta_.mntask = sktqaz_.ihacab[j];
  //   }
  //   // max hit time in OD:
  //   if (skta_.tamxsk < tqareal_ptr->T.at(j)){
  //     skta_.tamxsk= tqareal_ptr->T.at(j);
  //     // the PMT number with that time:
  //     skta_.mxtask = sktqaz_.ihacab[j];
  //   }

  //   //now the hits in the OD
  //   skta_.task[sktqaz_.ihacab[j] - 1] = tqareal_ptr->T.at(j);
  //   skqa_.qask[sktqaz_.ihacab[j] - 1] = tqareal_ptr->Q.at(j);

  //   if (j <= 5){
  //     std::cout << "skta_.task[sktqaz_.ihacab["<<j<<"] - 1] : " << skta_.task[sktqaz_.ihacab[j] - 1];
  //     std::cout << "skqa_.qask[sktqaz_.ihacab["<<j<<"<<] - 1] : " << skqa_.qask[sktqaz_.ihacab[j] - 1];
  //   }
    
  // }
    
  //Let's check everything's worked:
  PrintTQCommons(true, 5);
  //PrintTQCommons(false, 5);

  // if (debug == "DEBUG"){
  //   Compare_skq(skq_, skq_dupl);
  //   Compare_skt(skt_, skt_dupl);
  // }
 
  //now let's swap the pointers back:  
  // std::swap(skt_, m_data->skt_common_dupl);
  // std::swap(skq_, m_data->skq_common_dupl);
  
  return true;
}

bool CopyHits::Finalise(){

  return true;
}

void CopyHits::GetReader(){
  std::string tree_reader_str = "";
  m_variables.Get("tree_reader_str", tree_reader_str);
  if (tree_reader_str.empty() || m_data->Trees.count(tree_reader_str) == 0){
    throw std::invalid_argument("no valid treereader specified!");
  }
  tree_reader_ptr = m_data->Trees.at(tree_reader_str);
  if (tree_reader_ptr == nullptr){
    throw std::runtime_error("couldn't get treereader");
  }
  return;
}

bool CopyHits::PrintTQCommons(const bool& ID, const int& nhits){
  std::cout<<"\nPrinting Hits in "<<(ID ? "skq_, skt_" : "skta_, skqa_")<<"\n"<<std::endl;
	
  if(ID){
		
    if(skq_.nqisk==0){
      std::cout<<"skq_ had no ID hits"<<std::endl;
    } else {
			
      std::cout<<"Num ID hits in-1.3us (skq_.nqisk): "<<skq_.nqisk<<"\n"
	       <<"Total ID charge (skq_.qismsk): "<<skq_.qismsk<<" [p.e.]\n"
	       <<"ID PMT "<<skq_.mxqisk<<" saw the most charge of any ID PMT, with "<<skq_.qimxsk<<" [p.e.]\n"	
	       <<"ID PMT "<<skt_.mntisk<<" had the earliest hit in 1.3us of any ID PMT at "<<skt_.timnsk<<" [ns]\n"
	       <<"ID PMT "<<skt_.mxtisk<<" had the last hit in 1.3us of any ID PMT at "<<skt_.timxsk<<" [ns]\n";
			
      std::cout<<"first "<<nhits<<" ID hits (skq_,skt_: bad channel mask applied? within 1.3us of IT0XSK only):\n";
      for(int i=0, j=0; i<skq_.nqisk; ++i){
	int cableNumber = skchnl_.ihcab[i];
	if(cableNumber<=0 && cableNumber>MAXPM) continue;
	int iab = skchnl_.iab[cableNumber-1];
	// indexing seems like this from $SKOFL_ROOT/examples/skrd/dump_tq.F
	std::cout<<"ID Hit "<<j<<"\n"
		 <<"\tPMT: "<<cableNumber<<", iab: "<<iab<<" = "<<((iab==1) ? "A" : "B")<<"\n"
		 <<"\tT: skchnl_.itabsk: "<<skchnl_.itabsk[iab-1][cableNumber-1] // seems to be 0
		 <<", skt_.tisk: "<<skt_.tisk[cableNumber-1]<<"\n"
		 <<"\tQ: skchnl_.iqabsk: "<<skchnl_.iqabsk[iab-1][cableNumber-1] // seems to be 0
		 <<", skq_.qisk: "<<skq_.qisk[cableNumber-1]<<"\n";
	++j;
	if(j>=nhits) break;
      }
      std::cout<<std::endl;
    }
		
  } else {
		
    if(skqa_.nqask==0){
      std::cout<<"skqa had no OD hits"<<std::endl;
    } else {
			
      std::cout<<"Num OD hits in 1.3us gate (skqa_.nqask): "<<skqa_.nqask<<"\n"
	       <<"Total OD charge (skqa_.qasmsk): "<<skqa_.qasmsk<<" [p.e.]\n"
	       <<"OD PMT "<<skqa_.mxqask<<" saw the most charge of any OD PMT, with "<<skqa_.qamxsk<<" [p.e.]\n"
	       <<"OD PMT "<<skta_.mntask<<" had the earliest hit of any OD PMT at "<<skta_.tamnsk<<" [ns]\n"
	       <<"OD PMT "<<skta_.mxtask<<" had the last hit of any OD PMT at "<<skta_.tamxsk<<" [ns]\n";
			
      std::cout<<"first "<<nhits<<" OD hits (skqa_,skta_: bad channel mask applied? in 1.3us only):\n";
      for(int i=0, j=0; i<skqa_.nqask; ++i){ // FIXME are these in 1.3us or [-5,1]us?
	int cableNumber = sktqaz_.ihacab[i]; // note ihacab is in sktqaz_ NOT skchnla_ !!
	// ihacab does not contain OD PMT offset, so we need to add it,
	// BUT skta_.task and skqa_.qask are indexed WITHOUT the offset!
	int cableNumber2 = cableNumber + QB_OD_OFFSET;
	if(cableNumber2<=QB_OD_OFFSET || cableNumber2>(QB_OD_OFFSET+MAXPMA)) continue;
	std::cout<<"OD Hit "<<j<<"\n"
		 <<"\tPMT: "<<cableNumber2<<"\n"
		 <<"\tT: skta_.task: "<<skta_.task[cableNumber-1]<<"\n"  // XXX note indexing is
		 <<"\tQ: skqa_.qask: "<<skqa_.qask[cableNumber-1]<<"\n"; // WITHOUT QB_OD_OFFSET!!
	++j;
	if(j>=nhits) break;
      }
      std::cout<<std::endl;
    }
		
  }  
  return true;
}

void CopyHits::Compare_skq(const skq_common& orig, const skq_common& dupl) const {
  // nqisk:
  if (orig.nqisk != dupl.nqisk){std::cout << "skq_.nqisk does not match after copy!" << std::endl;}
  // qismsk:
  if (orig.qismsk != dupl.qismsk){std::cout << "skq_.qismsk does not match after copy!" << std::endl;}
  // qimxsk:
  if (orig.qimxsk != dupl.qimxsk){std::cout << "skq_.qimxsk does not match after copy!" << std::endl;}
  // mxqisk:
  if (orig.mxqisk != dupl.mxqisk){std::cout << "skq_.mxqisk does not match after copy!" << std::endl;}
  // qisk[MAXPM] 
  if (!CompareArray(orig.qisk, dupl.qisk, MAXPM)){std::cout << "skq_.qisk[MAXPM] does not match after copy!" << std::endl;}
  return;
}

void CopyHits::Compare_skt(const skt_common& orig, const skt_common& dupl) const {
  // timnsk:
  if (orig.timnsk != dupl.timnsk){std::cout << "skq_.timnsk does not match after copy!" << std::endl;}
  // timxsk:
  if (orig.timxsk != dupl.timxsk){std::cout << "skq_.timxsk does not match after copy!" << std::endl;}
  // mntisk:
  if (orig.mntisk != dupl.mntisk){std::cout << "skq_.mntisk does not match after copy!" << std::endl;}
  // mxtisk:
  if (orig.mxtisk != dupl.mxtisk){std::cout << "skq_.mxtisk does not match after copy!" << std::endl;}
  // tisk[MAXPM] 
  if (!CompareArray(orig.tisk, dupl.tisk, MAXPM)){std::cout << "skq_.tisk[MAXPM] does not match after copy!" << std::endl;}
  return;
}

template <typename T>
bool CopyHits::CompareArray(T* a, T* b, const int& N) const{
  bool match = true;
  for (int i = 0; i < N; ++i){match *= (a[i] == b[i]);}
  return match;
}
  


