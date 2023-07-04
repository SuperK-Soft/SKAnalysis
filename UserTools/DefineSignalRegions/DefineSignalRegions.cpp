#include "DefineSignalRegions.h"

DefineSignalRegions::DefineSignalRegions():Tool(){}


bool DefineSignalRegions::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  return true;
}


bool DefineSignalRegions::Execute(){

  /*
    1. Get event 
    2. Get reconstructed angle
    3. Get number of neutrons in event
    4. split signal into six regions
    5. fill into 1 of 6 histograms based on that split
   */
  
  return true;
}


bool DefineSignalRegions::Finalise(){

  return true;
}
