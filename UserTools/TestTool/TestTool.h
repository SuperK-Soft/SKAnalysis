#ifndef TestTool_H
#define TestTool_H

#include <string>
#include <iostream>
#include <map>

#include "Tool.h"
#include "Algorithms.h" // getOutputFromFunctionCall
#include "Constants.h"  // fundamental_types, container_types, TInterpreterErrors
#include "TFile.h"
#include "TTree.h"
#include <map>
#include <cxxabi.h>

//#include "TROOT.h"
//#include "TSystem.h"
#include "TInterpreter.h"
#include "TCint.h"
#include "TString.h"
#include <unistd.h>  // system
#include <cxxabi.h>  // demangle

/**
 * \class TestTool
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class TestTool: public Tool {


 public:

  TestTool(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Executre function used to perform Tool perpose. 
  bool Finalise(); ///< Finalise funciton used to clean up resorces.


 private:
        TFile* f;
        TTree* t;
        TTree* t2;
        
        TCint* meInterpreter=nullptr;
        TInterpreter::EErrorCode error;
        std::map<std::string, bool> known_types;
        std::map<std::string, std::string> underlying_types;
        
        // verbosity levels: if 'verbosity' < this level, the message type will be logged.
        int verbosity=1;
        int v_error=0;
        int v_warning=1;
        int v_message=2;
        int v_debug=3;
        std::string logmessage="";
        int get_ok=0;
        std::string toolName="TestTool";

        void MakeBranches(TTree* tree, Store* store);
        void MakeBranches(TTree* tree, BStore* store);
        void SetBranches(TTree* tree, Store* store);  // never used, addresses need to be
        void SetBranches(TTree* tree, BStore* store); // set internally as part of FillBranches
        bool FillBranches(TTree* tree, Store* store);
        bool FillBranches(TTree* tree, BStore* store);
        // overloads to handle a map of stores (hold tool configurations)
        void MakeBranches(TTree* tree, std::map<std::string, Store*> &stores);
        void SetBranches(TTree* tree, std::map<std::string, Store*> &stores);
        bool FillBranches(TTree* tree, std::map<std::string, Store*> &stores);
        
        bool RegisterTypes(std::string thetype);
        bool RegisterType(std::string thetype);
        std::string RemoveQualifiers(std::string thetype, bool allpointers=false);
        
        template<typename T>
        std::map<std::string,std::string> GetMap(T* store){
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
               try{
               // trim quotes
               keystring = keystring.substr(1,keystring.length()-2);
               valstring = valstring.substr(1,valstring.length()-2);
               } catch(...){
               std::cout<<"keystring length was "<<keystring.length()<<std::endl;
               std::cout<<"valstring length was "<<valstring.length()<<std::endl;
               }
               store_as_map.emplace(keystring,valstring);
           }
           return store_as_map;
        }
        
};

/* whoops this is a compile-time version, we need runtime.
template<typename T>{
  T temp;
}

template<typename C, typename T = typename C::value_type>
bool RegisterType(C const& container){
  T temp;
  return RegisterType(temp);
}

template<typename U, typename V>
bool RegisterType(std::pair<U,V> const& container){
  bool ok;
  U temp;
  ok = RegisterType(temp);
  V temp2;
  ok &= RegisterType(temp2);
  return ok;
}
*/

#endif
