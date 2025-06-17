#ifndef CompareSubtriggers_H
#define CompareSubtriggers_H

#include "TH1D.h"

#include <string>
#include <iostream>

#include "Tool.h"

class CompareSubtriggers: public Tool {

 public:

  CompareSubtriggers();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  sktqz_common sktqz_after;
  skchnl_common skchnl_after;

  std::string trigger_time_names = "";
  
  int TreeReaderLUN;
  int neglun;
  
  int trigger_idx = 0;
  int GetReaderLUN();

  bool made_plots = false;
  TH1D hits_before_inwindow;
  TH1D hits_after_inwindow;
  TH1D hits_before_outwindow;
  TH1D hits_after_outwindow;

  bool saved = false;
  
  void CompareSKTQZ(const sktqz_common& before, const skchnl_common& before_chnl,
		    const sktqz_common& after, const skchnl_common& after_chnl);

  void Compare_skq(const skq_common& orig, const skchnl_common& origchnl,
		 const skq_common& dupl, const skchnl_common& duplchnl) const;

  void Compare_skt(const skt_common& orig, const skchnl_common& origchnl, const skt_common& dupl,
		   const skchnl_common& duplchnl) const;


  
  
  template <typename T>
  bool CompareArray(T* a, T* b, const int& N);

  template <typename T>
  void PrintPMTArray(T* a, const skchnl_common& a_chnl,
                                   T* b, const skchnl_common& b_chnl);
  
  template <typename T>
  void PrintArray(T* a, T* b);

  
};


#endif
