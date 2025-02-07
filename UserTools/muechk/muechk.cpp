#include "muechk.h"

#include "MTreeReader.h"

#include "fortran_routines.h"


#include "TH1D.h"

muechk::muechk():Tool(){}

bool muechk::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  std::string reader_name = "";
  m_variables.Get("reader_name", reader_name);
  if (m_data->Trees.count(reader_name) != 1){
    throw std::runtime_error("muechk::Initialise: couldn't get treereader");
  }
  tree_reader_ptr = m_data->Trees.at(reader_name);

  nmue_plot = TH1D("nmue_plot", "nmue_plot", 20, 0, 0);
  
  return true;
}

bool muechk::Execute(){

  LoweInfo* lowe_ptr = nullptr;
  bool ok = tree_reader_ptr->Get("LOWE", lowe_ptr);
  if (!ok || lowe_ptr == nullptr){
    throw std::runtime_error("couldn't get lowe branch");
  }

  std::cout << "muechk::Execute: bsenergy: " << lowe_ptr->bsenergy <<std::endl;
  std::cout << "muechk::Execute: bsvertex[0]: " << lowe_ptr->bsvertex[0]<<std::endl;

  //HACK FOR THE TIME BEING:
  skwaterlen_.skwaterlen = 10000.;
  
  if (skwaterlen_.skwaterlen == 0.0){
    throw std::runtime_error("skwaterlen_.skwaterlen has value 0! apfit will exit with \"ripecorr.F : APABSPT iregal 0.00000000\" ");
  }
  
  /*
    apfit: apply reconstruction algorithms to one event:
    look in /usr/local/sklib_gcc8/atmpd-trunk/src/recon/aplib/apfit.F for description of input arguments 
    - 0 refers to "apply all reconstruction algorithms for fully contained events (FC)" 

    apfit calls apreset and muechk

    apreset: reset common blocks, calls:
    1) apclrall to clear apmring.h, apfitinf.h, appatsp.h, apmue,h
    2) appclrsep to clear apringsp.h
    3) skread(-1) to clear sktq.h <-- MAKE SURE YOU DON'T NEED THIS

    muechk: counts decay-e event , fills APMUE common
    1st arg: vertex of parent event, can be bsvertex
    2nd arg: verbosity, silent if == 1
    
    yes this tool should be called "apfit" in hindsight, but in hindsight lots of things seem simple don't they... [insert joke about doing a PhD here]

  */

  std::cout << "calling runinf_(): " << std::endl;
  runinfsk_();
  
  std::cout << "calling apfit(0):" << std::endl;
  int flag = 0;
  // apclrall_();
  // apclrsep_();
  apfit_(&flag);
  // apreset_();
  // int dummy_silent = 1;
  // muechk_(lowe_ptr->bsvertex, &dummy_silent);

  int nmue = apmue_.apnmue;
  nmue_plot.Fill(nmue);
  
  std::cout << "muechk::Execute: nmue = " << nmue << std::endl;
  for (int i = 0; i < nmue; ++i){
    std::cout << "muechk::Execute: found a decay electron with time: " << apmue_.apmuetime[i] << std::endl;
  }
  
  m_data->CStore.Set("nmue", nmue);
  
  return true;
}

bool muechk::Finalise(){

  nmue_plot.SaveAs("nmue_plot.root");
  
  return true;
}

extern "C" void lfhits_(int* lfflag, float* sg, float* bg, float* tspread,
		       float* signif, int* lfhit, float* tlf, float* qlf,
		       float* xyzlf, float* ttrg, int* n200, int* nbgbef,
		       int* nbgaft, int* nnba);


extern "C" void lfhits(int* lfflag, float* sg, float* bg, float* tspread,
	    float* signif, int* lfhit, float* tlf, float* qlf,
	    float* xyzlf, float* ttrg, int* n200, int* nbgbef,
	    int* nbgaft){

  std::cout << "calling decoy lfhits" << std::endl;
  
  int my_nnba = 0;
  lfhits_(lfflag, sg, bg, tspread,
	 signif, lfhit, tlf, qlf,
	 xyzlf, ttrg, n200, nbgbef,
	 nbgaft, &my_nnba);
  return;
  
}
