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

  TQReal* tqreal_ptr = nullptr;
  if (!tree_reader_ptr->Get("TQREAL", tqreal_ptr) || tqreal_ptr == nullptr){
    throw std::runtime_error("CopyHits::Execute: failed to get TQREAL branch");
  }

  // TQReal* tqareal_ptr = nullptr;
  // if (!tree_reader_ptr->Get("TQAREAL", tqareal_ptr) || tqareal_ptr == nullptr){
  //   throw std::runtime_error("CopyHits::Execute: failed to get TQAREAL branch");
  // }
  
  std::cout << "got TQREAL and TQAREAL branches, here are some hits:" << std::endl;
  std::cout << "TQREAL::nhits = " << tqreal_ptr->nhits << std::endl;
  for (int i = 0; i < std::min(tqreal_ptr->nhits, 5); ++i){
    std::cout << "TQREAL hit no: " << i
	      << " has time " << tqreal_ptr->T.at(i)
	      << " and charge " << tqreal_ptr->Q.at(i)
	      << std::endl;
  }
  // for (int i = 0; i < std::min(tqareal_ptr->nhits, 5); ++i){
  //   std::cout << "TQAREAL hit no: " << i
  // 	      << " has time " << tqareal_ptr->T.at(i)
  // 	      << " and charge " << tqareal_ptr->Q.at(i)
  // 	      << std::endl;
  // }

  /* 
     Now let's copy the info over to the skt_ and skq_ commons
     
     TQReal::T() has indexing running from 0 to TQReal::nhits-1
     skt_.tisk[] is indexed with the (PMT number - 1) given by skchnl_.ihcab

     For example, we need skt_.tisk[skchnl_.ihcab[i] - 1] == TQReal::T(i) for i : 0 -> TQReal::T().size()-1

   */

  // let's zero/max things where appropriate
  skq_.nqisk = 0;
  skq_.mxqisk = 0;
  skt_.timnsk = 9999999;
  skt_.timxsk = 0;

  skqa_.nqask = 0;
  skqa_.qasmsk = 0;
  skta_.tamnsk = 9999999;
  skta_.mxtask = 0;
  
  // ID:
  for (int i = 0; i < tqreal_ptr->nhits; ++i){
    // first we need to populate the skchnl_.ichab (PMT cable numbers) and IHTIFLZ (the hit flags)
    // the upper 16 bits are the cable numbers
    skchnl_.ihcab[i] = tqreal_ptr->cables.at(i) & 0x0000FFFF;
    // the lower 16 are the hit flags
    sktqz_.ihtiflz[i] = tqreal_ptr->cables.at(i) >> 16;
    if (!(sktqz_.ihtiflz[i] & 1)){continue;} // ignore hits not in 1.3us window around the primary trigger
    // number of ID hits:
    ++skq_.nqisk;
    // total ID charge:
    skq_.qismsk += tqreal_ptr->Q.at(i);
    // max charge deposited on ID PMT:
    if (tqreal_ptr->Q.at(i) > skq_.qimxsk){
      skq_.qimxsk = tqreal_ptr->Q.at(i);
      // the PMT number with that charge:
      skq_.mxqisk = skchnl_.ihcab[i];
    }
    // min hit time in ID:
    if (skt_.timnsk > tqreal_ptr->T.at(i)){
      skt_.timnsk = tqreal_ptr->T.at(i);
      // the PMT number with that time:
      skt_.mntisk = skchnl_.ihcab[i];
    }
    // max hit time in ID:
    if (skt_.timxsk < tqreal_ptr->T.at(i)){
      skt_.timxsk = tqreal_ptr->T.at(i);
      // the PMT number with that time:
      skt_.mxtisk = skchnl_.ihcab[i];
    }
    // then the actual hits in the ID:
    skt_.tisk[skchnl_.ihcab[i] - 1] = tqreal_ptr->T.at(i);
    skq_.qisk[skchnl_.ihcab[i] - 1] = tqreal_ptr->Q.at(i);
  }

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
