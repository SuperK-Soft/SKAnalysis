#ifndef StoreToTTree_H
#define StoreToTTree_H

#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <unistd.h>  // system
#include <cxxabi.h>  // demangle

#include "TFile.h"
#include "TTree.h"
//#include "TROOT.h"
//#include "TSystem.h"
#include "TInterpreter.h"
#include "TCint.h"
#include "TString.h"

#include "Store.h"
#include "BStore.h"
#include "Constants.h"  // TInterpreterErrors, fundamental_types, container_types

/**
* \class StoreToTTree
*
* StoreToTTree::MakeBranches takes in a Store and makes a TTree branch for each variable in the Store.
* StoreToTTree::FillBranches transfers the current contents of the BStore into a new TBranch entry.
* StoreToTTree accepts both (ASCII) Store class and (Binary) BStore classes.
*
* $Author: M.O'Flaherty $
* $Date: 2021/07/08 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/

// TODO: check BStore pointer handling is correct
// TODO: replace GetMap function once (B)Store map is made public
// TODO: check instances of TInterpreterErrors for kFatal (or kDangerous?)
//       and replace the internal TInterpreter when found, in case it is corrupt.

class StoreToTTree {
    public:
    StoreToTTree();
    ~StoreToTTree();
    
    // User functions
    void SetVerbosity(int verb){ verbosity=verb; }
    void MakeBranches(TTree* tree, Store* store);
    void MakeBranches(TTree* tree, BStore* store);
    void MakeBranches(TTree* tree, std::map<std::string, Store*> &stores);
    bool FillBranches(TTree* tree, Store* store);
    bool FillBranches(TTree* tree, BStore* store);
    bool FillBranches(TTree* tree, std::map<std::string, Store*> &stores);
    
    private:
    TCint* meInterpreter=nullptr;
    TInterpreter::EErrorCode error;
    std::map<std::string, bool> known_types;
    std::map<std::string, std::string> underlying_types;
    
    // verbosity levels: if 'verbosity' < this level, the message type will be logged.
    int verbosity=0;
    int get_ok=0;
    std::string toolName="StoreToTTree";
    
    // strongly discouraged
    bool allowRaggedBranches=false;
    
    // helper functions
    bool RegisterTypes(std::string thetype);
    bool RegisterType(std::string thetype);
    std::string RemoveQualifiers(std::string thetype, bool allpointers=false);
    template<typename T> std::map<std::string,std::string> GetMap(T* store);
    
};

// obsolete - can't handle commas in variables in Stores,
// json streamer for BStores sometimes corrupts... we have better methods now
template<typename T>
std::map<std::string,std::string> StoreToTTree::GetMap(T* store){
    // temporary stopgap while binary stores aren't available
    // we can use operator>> to dump the map into a json string
    std::map<std::string,std::string> store_as_map;
    std::string store_as_json;
    (*store) >> store_as_json;
    // now parse the json to produce a map
    // trim curly braces
    store_as_json = store_as_json.substr(1,store_as_json.length()-2);
    std::stringstream sss;
    sss << store_as_json;
    std::string pairstring;
    while(std::getline(sss,pairstring,',')){
       std::stringstream ssss;
       ssss << pairstring;
       std::string keystring;
       std::string valstring;
       std::getline(ssss,keystring,':');
       std::getline(ssss,valstring,':');
       // trim quotes
       keystring = keystring.substr(1,keystring.length()-2);
       valstring = valstring.substr(1,valstring.length()-2);
       store_as_map.emplace(keystring,valstring);
   }
   return store_as_map;
}

#endif
