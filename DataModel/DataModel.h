#ifndef DATAMODEL_H
#define DATAMODEL_H

#if PYTHON==1
// Pybind11 modules
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#endif

#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <functional>

//#include "TTree.h"
#include "TApplication.h"

#include "Store.h"
#include "BStore.h"
#include "Logging.h"
#include "LoggingLevels.h"
#include "Utilities.h"
#include "StoreToTTree.h"
#include "Constants.h"

#include "EventCandidates.h"
#include "EventParticles.h"
#include "EventTrueCaptures.h"
#include "PMTHitCluster.h"

#include "ParticleCand.h"

#include "skroot_loweC.h"

class MTreeReader;
class MTreeSelection;
class TreeReader;
class ConnectionTable;

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
  ~DataModel(); ///< Simple destructor

  //TTree* GetTTree(std::string name);
  //void AddTTree(std::string name,TTree *tree);
  //void DeleteTTree(std::string name,TTree *tree);
  TApplication* GetTApp();
  ConnectionTable* GetConnectionTable(int sk_geometry=-1);

  Logging *Log; ///< Log class pointer for use in Tools, it can be used to send messages which can have multiple error levels and destination end points
  std::map<std::string,MTreeReader*> Trees; ///< A map of MTreeReader pointers, used to read ROOT trees
  std::map<std::string,MTreeSelection*> Selectors; ///< A map of MTreeSelection pointers used to read event selections
  std::unordered_map<std::string, std::function<bool()>> hasAFTs;
  std::unordered_map<std::string, std::function<bool()>> loadSHEs;
  std::unordered_map<std::string, std::function<bool()>> loadAFTs;
  std::unordered_map<std::string, std::function<bool(int)>> loadCommons;
  std::unordered_map<std::string, std::function<int(long, bool)>> getEntries;
  
  Store vars; ///< This Store can be used for any variables. It is an inefficent ascii based storage and command line arguments will be placed in here along with ToolChain variables
  BStore CStore; ///< This is a more efficent binary Store that can be used to store a dynamic set of inter Tool variables, very useful for constants and and flags hence the name CStore
  std::map<std::string,BStore*> Stores;  ///< This is a map of named BStore pointers which can be deffined to hold a nammed collection of any type of BStore. It is usefull to store data collections that needs subdividing into differnt stores.
  
  // This function is used to register a TreeReader tool's member functions with the DataModel,
  // which provides access from other Tools
  bool RegisterReader(std::string readerName, MTreeReader* reader, std::function<bool()> hasAFT, std::function<bool()> loadSHE, std::function<bool()> loadAFT, std::function<bool(int)> loadCommon, std::function<int(long,bool)> getTreeEntry);
  int getTreeEntry(std::string ReaderName="", long entrynum=0);
  // These retain function pointers to call the corresponding TreeReader functions.
  // The TreeReader instance is obtained from the name specified in their config file.
  bool HasAFT(std::string ReaderName="");
  bool LoadSHE(std::string ReaderName="");
  bool LoadAFT(std::string ReaderName="");
  bool LoadEntry(int entry_i, std::string ReaderName="");
  // tracking fortran logic unit numbers (LUNs, file handles)
  int GetNextLUN(int lun=10, std::string reader="reader");
  int GetLUN(std::string reader);
  bool FreeLUN(int lun, std::string reader="");
  // wrapper to check if we've called this yet, since we should probably only call it the once?
  void KZInit();
  bool kz_initialized=false;
  
  
  // Event vars
  BStore* eventVariables_p; // TODO replace with a pointer and update tools to use -> instead of .
  BStore &eventVariables;   // use references to preserve current behaviour...
  
  // Hits
  PMTHitCluster eventPMTHits;
  // Candidates
  EventCandidates eventCandidates;
  // Primaries (MC)
  EventParticles eventPrimaries;
  // Secondaries (MC)
  EventParticles eventSecondaries;
  // Neutron captures (MC)
  EventTrueCaptures eventTrueCaptures;
  
  std::map<std::string, Store*> tool_configs;
  StoreToTTree StoreConverter;
  
  //flag for a new muon - used for the spallation reduction
  bool newMuon = false;
  
  //flag for a new relic candidate - used for the spallation reduction
  bool newRelic = false;
  
  //deques of the struct ParticleCand to store event info of ALL muon candidates and relic candidates
  std::deque<ParticleCand> muonCandDeque;
  std::deque<ParticleCand> relicCandDeque;
  
  //deque for muons that need to be reconstructed (if they have matched to a relic they need to be reconstructed)
  std::vector<ParticleCand> muonsToRec;
  
  //vector for relic candidate events that needs to be written out - used for the spallation reduction
  std::vector<ParticleCand> writeOutRelics;
  
  //map for lowe information to store when reconstructing events mid-chain and wanting to retrieve the information later
  std::map<long, skroot_lowe_common> loweCommonBufferMap;
  
  //vector to store weightings in from various SRN models at energy intervals of 0.5 MeV.
  std::vector<std::vector<float>> SRNWeights;
  
  //vector to store modelNames for the SRN Models in
  std::vector<std::string> SRNModelNames;
  
  bool applyReweight = false;
  
 private:


  //std::map<std::string,TTree*> m_trees; 
  TApplication* rootTApp=nullptr;
  ConnectionTable* connectionTable=nullptr;
  
  
  
};



#endif
