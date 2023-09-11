#include "muechk.h"

#include "MTreeReader.h"

#include "fortran_routines.h"

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
 
  return true;
}

bool muechk::Execute(){

  LoweInfo* lowe_ptr = nullptr;
  bool ok = tree_reader_ptr->Get("LOWE", lowe_ptr);
  if (!ok || lowe_ptr == nullptr){
    throw std::runtime_error("couldn't get lowe branch");
  }

  int dummy_silent = 1;

  Log("muechk::Execute: skheadg_.sk_geometry: "+std::to_string(skheadg_.sk_geometry), 0, 0);
  
  muechk_(lowe_ptr->bsvertex, &dummy_silent); //common APMUE now filled

  // char nmuestr[80];
  // char tmpstr[20];
  // int maxi = 0;
  // int i = 0;
  int nmue = apmue_.apnmue;
  // sprintf(nmuestr,"==== nmue = %d ",nmue);
  // if (nmue > 0) {
  //   sprintf(tmpstr,"  t(us) = ");
  //   strcat(nmuestr,tmpstr);
  //   if (nmue > 10)
  //     maxi = 10;
  //   else
  //     maxi = nmue;
  //   for (i=0; i < maxi; i++) {
  //     sprintf(tmpstr,"%5.2f ",apmue_.apmuetime[i]);
  //      strcat(nmuestr,tmpstr);
  //   }
  // }
  // printf("%s\n",nmuestr);

  Log("muechk::Execute: nmue = "+std::to_string(nmue),0,0);
  for (int i = 0; i < nmue; ++i){
    if (i < 10){
      Log("muechk::Execute: found a decay electron with time: "+std::to_string(apmue_.apmuetime[i]),0, 0);
    }
  }
  
  return true;
}

bool muechk::Finalise(){

  return true;
}
