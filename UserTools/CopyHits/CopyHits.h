#ifndef CopyHits_H
#define CopyHits_H

#include <string>
#include <iostream>

#include "MTreeReader.h"

#include "Tool.h"

class CopyHits: public Tool {

 public:

  MTreeReader* tree_reader_ptr = nullptr;
  
  CopyHits();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();
  void GetReader();
  
 private:

  bool PrintTQCommons(const bool&, const int&);
  
};


#endif
