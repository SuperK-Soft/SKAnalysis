#ifndef DATAMODEL_H
#define DATAMODEL_H

#include <map>
#include <string>
#include <vector>
#include "Utilities.h"

//#include "TTree.h"

#include "Store.h"
#include "BStore.h"
#include "Logging.h"

/**
 * \class DataModel
 *
 * This class Is a transient data model class for your Tools within the ToolChain. If Tools need to comunicate they pass all data objects through the data model. There fore inter tool data objects should be deffined in this class. 
 *
 *
 * $Author: B.Richards $ 
 * $Date: 2019/05/26 18:34:00 $
 *
 */

class DataModel {


 public:
  
  DataModel(); ///< Simple constructor 

  //TTree* GetTTree(std::string name);
  //void AddTTree(std::string name,TTree *tree);
  //void DeleteTTree(std::string name,TTree *tree);

  Logging *Log; ///< Log class pointer for use in Tools, it can be used to send messages which can have multiple error levels and destination end points
  
  Store vars; ///< This Store can be used for any variables. It is an inefficent ascii based storage and command line arguments will be placed in here along with ToolChain variables
  BStore CStore; ///< This is a more efficent binary Store that can be used to store a dynamic set of inter Tool variables, very useful for constants and and flags hence the name CStore
  std::map<std::string,BStore*> Stores;  ///< This is a map of named BStore pointers which can be deffined to hold a nammed collection of any type of BStore. It is usefull to store data collections that needs subdividing into differnt stores.
  
  
  
 private:


  //std::map<std::string,TTree*> m_trees; 
  
  
  
};



#endif
