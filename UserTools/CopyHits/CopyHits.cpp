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
     There is probably some deep clandestine reason why set_timing_gate() is removing these hits
     Rather than finding this, we'll just do things quick and dirty - copying the hits from the TQREAL branch to the required commons.
   */

  TQReal* tqreal_ptr = nullptr;
  if (!tree_reader_ptr->Get("TQREAL", tqreal_ptr) || tqreal_ptr == nullptr){
    throw std::runtime_error("CopyHits::Execute: failed to get TQREAL branch");
  }

  TQReal* tqareal_ptr = nullptr;
  if (!tree_reader_ptr->Get("TQAREAL", tqareal_ptr) || tqareal_ptr == nullptr){
    throw std::runtime_error("CopyHits::Execute: failed to get TQAREAL branch");
  }
  
  std::cout << "got TQREAL and TQAREAL branches, here are some hits:" << std::endl;
  std::cout << "TQREAL::nhits = " << tqreal_ptr->nhits << std::endl;
  for (int i = 0; i < std::min(tqreal_ptr->nhits, 5); ++i){
    std::cout << "TQREAL hit no: " << i
	      << " has time " << tqreal_ptr->T.at(i)
	      << " and charge " << tqreal_ptr->Q.at(i)
	      << std::endl;
  }
  for (int i = 0; i < std::min(tqareal_ptr->nhits, 5); ++i){
    std::cout << "TQAREAL hit no: " << i
	      << " has time " << tqareal_ptr->T.at(i)
	      << " and charge " << tqareal_ptr->Q.at(i)
	      << std::endl;
  }

  /* 
     Now let's copy the info over to the skt_ and skq_ commons
     
     TQReal::T() has indexing running from 0 to TQReal::nhits
     skt_.tisk[] is indexed with the (PMT number - 1) given by skchnl_.ihcab

     For example, we need skt_.tisk[skchnl_.ihcab[i] - 1] == TQReal::T(i) for i : 0 -> TQReal::T().size()-1

   */

  // first skq_  and skqa_ common blocks - minus the hits:
  // number of ID hits:
  skq_.nqisk = tqreal_ptr->nhits;
  // number of OD hits:
  skqa_.nqask = tqareal_ptr->nhits;
  // total charge in ID:
  skq_.qismsk = std::accumulate(tqreal_ptr->Q.begin(), tqreal_ptr->Q.end(), 0);
  // total charge in OD:
  skqa_.qasmsk = std::accumulate(tqareal_ptr->Q.begin(), tqareal_ptr->Q.end(), 0);

  const auto max_ID_charge_it = std::max_element(tqreal_ptr->Q.begin(), tqreal_ptr->Q.end());
  // charge on ID PMT with max charge:
  skq_.qimxsk = *max_ID_charge_it;
  // PMT number of that PMT:
  skq_.mxqisk = skchnl_.ihcab[std::distance(tqreal_ptr->Q.begin(), max_ID_charge_it)];
  const auto max_OD_charge_it = std::max_element(tqareal_ptr->Q.begin(), tqareal_ptr->Q.end());
  // charge on OD PMT with max charge:
  skqa_.qamxsk = *max_OD_charge_it;
  // PMT number of that PMT:
  skqa_.mxqask = skchnl_.ihcab[std::distance(tqareal_ptr->Q.begin(), max_OD_charge_it)] + QB_OD_OFFSET;

  // then skt_ and skta_ common blocks - again minus the hits:
  const auto min_ID_time_it = std::min_element(tqreal_ptr->T.begin(),tqreal_ptr->T.end());
  // min hit time in ID:
  skt_.timnsk = *min_ID_time_it;
  const auto min_OD_time_it = std::min_element(tqareal_ptr->T.begin(),tqareal_ptr->T.end());
  // min hit time in OD:
  skta_.tamnsk = *min_OD_time_it;
  const auto max_ID_time_it = std::max_element(tqreal_ptr->T.begin(),tqreal_ptr->T.end());
  // max hit time in ID:
  skt_.timxsk = *max_ID_time_it;
  const auto max_OD_time_it = std::max_element(tqareal_ptr->T.begin(),tqareal_ptr->T.end());
  // max hit time in OD:
  skta_.tamxsk = *max_OD_time_it;
  // PMT number for ID PMT that collected the hit with min time:
  skt_.mntisk = skchnl_.ihcab[std::distance(tqreal_ptr->T.begin(), min_ID_time_it)];
  // PMT number for OD PMT that collected the hit with min time:
  skta_.mntask = skchnl_.ihcab[std::distance(tqreal_ptr->T.begin(), max_ID_time_it)] + QB_OD_OFFSET;
  // PMT number for ID PMT that collected the hit with max time:
  skt_.mxtisk = skchnl_.ihcab[std::distance(tqareal_ptr->T.begin(), min_OD_time_it)];
  // PMT number for OD PMT that collected the hit with max time:
  skta_.mxtask = skchnl_.ihcab[std::distance(tqareal_ptr->T.begin(), max_OD_time_it)] + QB_OD_OFFSET;
  
  // now the hits in the ID
  for (int i = 0; i < tqreal_ptr->nhits; ++i){
    skt_.tisk[skchnl_.ihcab[i] - 1] = tqreal_ptr->T.at(i);
    skq_.qisk[skchnl_.ihcab[i] - 1] = tqreal_ptr->Q.at(i);
  }

  // and finally the hits in the OD
  for (int j = 0; j < tqareal_ptr->nhits; ++j){
    skta_.task[skchnl_.ihcab[j] - 1 + QB_OD_OFFSET] = tqareal_ptr->T.at(j);
    skqa_.qask[skchnl_.ihcab[j] - 1 + QB_OD_OFFSET] = tqareal_ptr->Q.at(j);
  }

  //Let's check everything's worked:
  PrintTQCommons(true, 5);
  PrintTQCommons(false, 5);
    
   
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
