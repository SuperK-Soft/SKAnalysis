#include "BuildHist.h"

BuildHist::BuildHist():Tool(){}


bool BuildHist::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  return true;
}


bool BuildHist::Execute(){

  return true;
}


bool BuildHist::Finalise(){

  return true;
}
