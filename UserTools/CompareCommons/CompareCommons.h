#ifndef CompareCommons_H
#define CompareCommons_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class CompareCommons
*
* Please fill out the descripton and author information.
*
* $Author: ?.????? $
* $Date: ????/??/?? $
* $Contact: ???@km.icrr.u-tokyo.ac.jp
*/

class CompareCommons: public Tool {

 public:

  CompareCommons();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  void Compare_skq(const skq_common&, const skchnl_common&, const skq_common&,  const skchnl_common&) const;
  void Compare_skt(const skt_common&, const skchnl_common&, const skt_common&, const skchnl_common&) const;
  template <typename T>
  bool CompareArray(T*, T*, const int&) const;
  
  
  
};


#endif
