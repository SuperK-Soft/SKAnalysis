#include "call_muechk.h"

#include "MTreeReader.h"

#include "fortran_routines.h"

call_muechk::call_muechk():Tool(){}

bool call_muechk::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  if (m_data->Trees.count("reader")!=1){
    throw std::runtime_error("couldn't get treereader");
  }

  tree_reader_ptr = m_data->Trees.at("reader");
 
  return true;
}

bool call_muechk::Execute(){

  LoweInfo* lowe_ptr = nullptr;
  bool ok = tree_reader_ptr->Get("LOWE", lowe_ptr);
  if (!ok || lowe_ptr == nullptr){
    throw std::runtime_error("couldn't get lowe branch");
  }

  int dummy_silent = 1;

  std::cout << "skheadg_.sk_geometry: " << skheadg_.sk_geometry << "\n";
  
  int nmue = 0;
  muechk_(lowe_ptr->bsvertex, &dummy_silent); //common APMUE now filled

  char nmuestr[80];
  char tmpstr[20];
  int maxi = 0;
  int i = 0;
  nmue = apmue_.apnmue;
  sprintf(nmuestr,"==== nmue = %d ",nmue);
  if (nmue > 0) {
    sprintf(tmpstr,"  t(us) = ");
    strcat(nmuestr,tmpstr);
    if (nmue > 10)
      maxi = 10;
    else
      maxi = nmue;
    for (i=0; i < maxi; i++) {
      sprintf(tmpstr,"%5.2f ",apmue_.apmuetime[i]);
       strcat(nmuestr,tmpstr);
    }
  }
  printf("%s\n",nmuestr);

  
  return true;
}

bool call_muechk::Finalise(){

  return true;
}
