#ifndef MYTOOL_H
#define MYTOOL_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class MyTool
*
* Please fill out the descripton and author information.
*
* $Author: ?.????? $
* $Date: ????/??/?? $
* $Contact: ???@km.icrr.u-tokyo.ac.jp
*/

class MyTool: public Tool {

 public:

  MyTool();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

};


#endif
