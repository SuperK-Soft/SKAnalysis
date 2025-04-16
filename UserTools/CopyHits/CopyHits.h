#ifndef CopyHits_H
#define CopyHits_H

#include <string>
#include <iostream>

#include "MTreeReader.h"
#include "fortran_routines.h"

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

  // skq_common skq_dupl;
  // skt_common skt_dupl;
  
  bool PrintTQCommons(const bool&, const int&);

  void Compare_skq(const skq_common&, const skq_common&) const;
  void Compare_skt(const skt_common&, const skt_common&) const;

  template <typename T>
  bool CompareArray(T*, T*, const int&) const;
};


#endif
