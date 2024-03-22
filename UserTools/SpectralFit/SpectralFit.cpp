#include "SpectralFit.h"

SpectralFit::SpectralFit():Tool(){}

bool SpectralFit::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;
  
  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  
  GetPDFs();
  
  return true;
}

bool SpectralFit::Execute(){

  // line up all regions along the x axis to make a global function
  // normalise distributions across all regions
  // construct functional fit
  // do fit
  // split up back into region plots
  // be happy

  // for (int i = 0; i < N_distributions){
    
  // }
  
  return true;
}

bool SpectralFit::Finalise(){
  return true;
}

void SpectralFit::GetPDFs(){
  return;
}
