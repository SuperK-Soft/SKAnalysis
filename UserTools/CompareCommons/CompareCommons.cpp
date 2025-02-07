#include "CompareCommons.h"

CompareCommons::CompareCommons():Tool(){}


bool CompareCommons::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  return true;
}


bool CompareCommons::Execute(){

  // first print out subtrigger times to compare software trigger and SLESearch tool:

  // std::vector<int> t0_sub = {};
  // if (!m_data->CStore.Get("trigger_times", t0_sub)){throw std::runtime_error("CompareCommons::Execute couldn't retrieve subtrigger times");}

  // std::vector<double> SLE_times = {};
  // if (!m_data->CStore.Get("SLE_times", SLE_times)){throw std::runtime_error("CompareCommons::Execute couldn't retrieve SLESearch subtrigger times");}

  // if (t0_sub.empty() || SLE_times.empty()){
  //   std::cout << "triggers empty" << std::endl;
  //   return true;
  // }

  // std::cout << "software trigger size: " << t0_sub.size() << std::endl;
  // std::cout << "SLE Search size: " << SLE_times.size() << std::endl;
  
  // std::cout << "software trigger:" << std::endl;
  // for (const auto& t : t0_sub){
  //   std::cout << t + skheadqb_.it0sk << std::endl;
  // }

  // std::cout << "SLESearch:" << std::endl;
  // for (const auto& t : SLE_times){
  //   std::cout << t << std::endl;
  // }

  //now compare commons and commons duplicated
  Compare_skq(skq_, skchnl_, m_data->skq_common_dupl, m_data->skchnl_common_dupl);
  Compare_skt(skt_, skchnl_, m_data->skt_common_dupl, m_data->skchnl_common_dupl);
  
  return true;
}


bool CompareCommons::Finalise(){

  return true;
}

void CompareCommons::Compare_skq(const skq_common& orig, const skchnl_common& origchnl, const skq_common& dupl, const skchnl_common& duplchnl) const {
  // nqisk:                                                                                                                                                 
  if (orig.nqisk != dupl.nqisk){
    std::cout << "skq_.nqisk does not match after copy!" << std::endl;
    std::cout << "orig.nqisk: " << orig.nqisk << ", dupl.nqisk: " << dupl.nqisk << std::endl; 
  }
  // qismsk:                                                                                                                                                
  if (orig.qismsk != dupl.qismsk){
    std::cout << "skq_.qismsk does not match after copy!" << std::endl;
    std::cout << "photoelectrons" << std::endl;
    std::cout << "orig.qismsk: " << orig.qismsk << ", dupl.qismsk: " << dupl.qismsk << std::endl; 
  }
  
  // qimxsk:
  if (orig.qimxsk != dupl.qimxsk){
    std::cout << "skq_.qimxsk does not match after copy!" << std::endl;
    std::cout << "photoelectrons" << std::endl;
    std::cout << "orig.qimxsk: " << orig.qimxsk << ", dupl.qimxsk: " << dupl.qimxsk << std::endl;
  }

  // mxqisk:
  if (orig.mxqisk != dupl.mxqisk){
    std::cout << "skq_.mxqisk does not match after copy!" << std::endl;
    std::cout << "pmt number that acquired max charge" << std::endl;
    std::cout << "orig.mxqisk: " << orig.mxqisk << ", dupl.mxqisk: " << dupl.mxqisk << std::endl; 
  
  }
  // qisk[MAXPM]
  // if (!CompareArray(orig.qisk, dupl.qisk, MAXPM)){
  //   std::cout << "skq_.qisk[MAXPM] does not match after copy!" << std::endl;
    std::cout << "some values from the arrays: " << std::endl;
    std::cout << "photoelectrons" << std::endl;
    for (int i = 0; i < 20/*std::min(std::min(20, orig.nqisk), dupl.nqisk)*/; ++i){
      int orig_cableNumber = origchnl.ihcab[i];
      int dupl_cableNumber = duplchnl.ihcab[i];
      std::cout << "orig.qisk["<<orig_cableNumber<<"-1] " <<  orig.qisk[orig_cableNumber-1] << ", dupl.qisk["<<dupl_cableNumber<<"-1] " << dupl.qisk[dupl_cableNumber-1] << std::endl;
    }
  // }
  return;
}

void CompareCommons::Compare_skt(const skt_common& orig, const skchnl_common& origchnl, const skt_common& dupl, const skchnl_common& duplchnl) const {
  // timnsk:                                                                       
  if (orig.timnsk != dupl.timnsk){
    std::cout << "skt_.timnsk does not match after copy!" << std::endl;
    std::cout << "ns:  ";
    std::cout << "orig.timnsk: " << orig.timnsk << ", dupl.timnsk: " << dupl.timnsk << std::endl;
  }
  // timxsk:                                                                       
  if (orig.timxsk != dupl.timxsk){
    std::cout << "ns:  ";
    std::cout << "skt_.timxsk does not match after copy!" << std::endl;
    std::cout << "orig.timxsk: " << orig.timxsk << ", dupl.timxsk: " << dupl.timxsk << std::endl;
  }
  // mntisk:                                                                       
  if (orig.mntisk != dupl.mntisk){
    std::cout << "ns:  ";
    std::cout << "skt_.mntisk does not match after copy!" << std::endl;
    std::cout << "orig.mntisk: " << orig.mntisk << ", dupl.mntisk: " << dupl.mntisk << std::endl;
  }
  // mxtisk:                                                                       
  if (orig.mxtisk != dupl.mxtisk){
    std::cout << "ns:  ";
    std::cout << "skt_.mxtisk does not match after copy!" << std::endl;
    std::cout << "orig.mxtisk: " << orig.mxtisk << ", dupl.mxtisk: " << dupl.mxtisk << std::endl;
  }
  // tisk[MAXPM]                                                                   
  // if (!CompareArray(orig.tisk, dupl.tisk, MAXPM)){
  //   std::cout << "skt_.tisk[MAXPM] does not match after copy!" << std::endl;
    std::cout << "some values from the arrays: " << std::endl;
    std::cout << "ns:  ";
    for (int i = 0; i < 20; i++){
      int orig_cableNumber = origchnl.ihcab[i];
      int dupl_cableNumber = duplchnl.ihcab[i];
      std::cout << "orig.tisk["<<orig_cableNumber<<"-1] " <<  orig.tisk[orig_cableNumber-1] << ", dupl.tisk["<<dupl_cableNumber<<"-1] " << dupl.tisk[dupl_cableNumber-1] << std::endl;}
  // }
  return;
}

template <typename T>
bool CompareCommons::CompareArray(T* a, T* b, const int& N) const{
  bool match = true;
  for (int i = 0; i < N; ++i){match *= (a[i] == b[i]);}
  return match;
}

