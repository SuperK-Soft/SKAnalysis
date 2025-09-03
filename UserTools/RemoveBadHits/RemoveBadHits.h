#ifndef RemoveBadHits_H
#define RemoveBadHits_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"

/**
* \class RemoveBadHits
*
* Please fill out the descripton and author information.
*
* $Author: ?.????? $
* $Date: ????/??/?? $
* $Contact: ???@km.icrr.u-tokyo.ac.jp
*/

class RemoveBadHits: public Tool {

 public:

  RemoveBadHits();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
  MTreeReader* tree_reader;

};


#endif
