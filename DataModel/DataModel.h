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
#include "Algorithms.h"

#include "EventCandidates.h"
#include "EventParticles.h"
#include "EventTrueCaptures.h"
#include "PMTHitCluster.h"

#include "MParticle.h"
#include "NCapture.h"
#include "NCaptCandidate.h"
#include "LoweCandidate.h"

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "MParticle.h"
#include "MVertex.h"

#include "ParticleCand.h"
#include "skroot_loweC.h"

#include "MTreeSelection.h"

class MTreeReader;
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
  
  // some of our helper classes in the DataModel dir could benefit from access
  // to the datamodel. make it a singleton (we only have one anyway)
  static DataModel* GetInstance(){ return thisptr; }
  static DataModel* thisptr;

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
  std::unordered_map<std::string, std::function<int(long)>> getEntrys;
  MTreeReader* GetTreeReader();
  
  Store vars; ///< This Store can be used for any variables. It is an inefficent ascii based storage and command line arguments will be placed in here along with ToolChain variables
  BStore CStore; ///< This is a more efficent binary Store that can be used to store a dynamic set of inter Tool variables, very useful for constants and and flags hence the name CStore
  std::map<std::string,BStore*> Stores;  ///< This is a map of named BStore pointers which can be deffined to hold a nammed collection of any type of BStore. It is usefull to store data collections that needs subdividing into differnt stores.
  
  // This function is used to register a TreeReader tool's member functions with the DataModel,
  // which provides access from other Tools
  bool RegisterReader(std::string readerName, MTreeReader* reader, std::function<bool()> hasAFT={}, std::function<bool()> loadSHE={}, std::function<bool()> loadAFT={}, std::function<bool(int)> loadCommon={}, std::function<int(long)> getTreeEntry={});
  int getTreeEntry(std::string ReaderName="", long entrynum=0, bool justdoit=false);
  // These retain function pointers to call the corresponding TreeReader functions.
  // The TreeReader instance is obtained from the name specified in their config file.
  bool HasAFT(std::string ReaderName="");
  bool LoadSHE(std::string ReaderName="");
  bool LoadAFT(std::string ReaderName="");
  bool LoadCommons(int entry_i, std::string ReaderName="");
  // tracking fortran logic unit numbers (LUNs, file handles)
  int GetNextLUN(std::string reader="reader", int lun=0);
  int GetLUN(std::string reader);
  bool FreeLUN(int lun, std::string reader="");
  // wrapper to check if we've called this yet, since we should probably only call it the once?
  void KZInit();
  bool kz_initialized=false;
  bool GeoSet(int sk_geometry_in);
  bool BonsaiInit();
  bool bonsai_initialised=false;
  
  TFile* OpenFileForWriting(std::string file, bool alreadyopenonly=false);
  TFile* OpenFileForReading(std::string file, bool fromdisk=true);
  bool CloseFile(std::string file);
  bool CloseFile(TFile* fptr);
  
  // helper functions for working with MTreeSelections
  // (basically just needed for the case of adding cuts to multiple selectors at once)
  template<typename... Args>
  bool AddCut(std::string selector, std::string cut, std::string description, Args... rest);
  template<typename... Args>
  bool ApplyCut(std::string selector, std::string cut, double val, Args... rest);
  template<typename... Args>
  bool AddPassingEvent(std::string selector, std::string cut, Args... rest);
  
  // Event vars
  BStore* eventVariables_p; // TODO replace with a pointer and update tools to use -> instead of .
  BStore &eventVariables;   // use references to preserve current behaviour...
  
  // NTag classes
  PMTHitCluster eventPMTHits;
  EventCandidates eventCandidates;
  EventParticles eventPrimaries;
  EventParticles eventSecondaries;
  EventTrueCaptures eventTrueCaptures;
  
  const TDatabasePDG* pdgdb = TDatabasePDG::Instance();
  std::vector<MParticle> eventParticles;
  std::vector<MVertex> eventVertices;
  // generalised neutron captures
  std::map<std::string,std::vector<NCaptCandidate>> NCaptureCandidates;
  std::vector<LoweCandidate> LoweCandidates;
  std::vector<NCapture> NCapturesTrue;
  
  std::map<std::string, Store*> tool_configs;
  StoreToTTree StoreConverter;
  
  // for muon-lowe matching
  // ----------------------
  bool newMuon = false;   // flag for a new muon
  bool newRelic = false;  //flag for a new relic candidate
  
  //deques of ALL muon candidates and relic candidates (before/during matching)
  std::deque<ParticleCand> muonCandDeque;
  std::deque<ParticleCand> relicCandDeque;
  
  //deque of muons that need to be reconstructed (i.e. those matched to a relic candidate)
  std::vector<ParticleCand> muonsToRec;
  
  //vector of relic candidate ready to be written out (muon scan done)
  std::vector<ParticleCand> writeOutRelics;
  
  //cached lowe common blocks, for use during matching
  std::map<long, skroot_lowe_common> loweCommonBufferMap;
  
  // -----------------------
  
 private:


  //std::map<std::string,TTree*> m_trees; 
  TApplication* rootTApp=nullptr;
  ConnectionTable* connectionTable=nullptr;
  
  // output ROOT files, for sharing between Tools
  // use OpenFile and CloseFile functions to access this.
  // key should be a full filepath (since it may be required to create the file)
  // value first entry is file pointer, second is number of active users
  std::map<std::string, std::pair<TFile*, int>> OutFiles;
  
  
  
};


template<typename... Args>
bool DataModel::AddCut(std::string selector, std::string cut, std::string description, Args... rest){
	if(selector!="all" && Selectors.count(selector)==0){
		std::cerr<<"DataModel::AddCut Error! Unrecognised selector "<<selector
		         <<" for cut "<<cut<<std::endl;
		return false;
	}
	bool ret=true;
	if(selector=="all"){
		// add this event to all selectors
		for(auto sel = Selectors.begin(); sel!=Selectors.end(); ++sel){
			ret &= sel->second->AddCut(cut, description, rest...);
		}
	} else {
		ret = Selectors.at(selector)->AddCut(cut, description, rest...);
	}
	return ret;
}

template<typename... Args>
bool DataModel::ApplyCut(std::string selector, std::string cut, double val, Args... rest){
	if(selector!="all" && Selectors.count(selector)==0){
		std::cerr<<"DataModel::ApplyCut Error! Unrecognised selector "<<selector
		         <<" for cut "<<cut<<std::endl;
		return false;
	}
	bool ret=true;
	if(selector=="all"){
		// add this event to all selectors
		for(auto sel = Selectors.begin(); sel!=Selectors.end(); ++sel){
			ret &= sel->second->ApplyCut(cut, val, rest...);
		}
	} else {
		ret = Selectors.at(selector)->ApplyCut(cut, val, rest...);
	}
	return ret;
}

template<typename... Args>
bool DataModel::AddPassingEvent(std::string selector, std::string cut, Args... rest){
	if(selector!="all" && Selectors.count(selector)==0){
		std::cerr<<"DataModel::AddPassingEvent Error! Unrecognised selector "<<selector
		         <<" for cut "<<cut<<std::endl;
		return false;
	}
	bool ret=true;
	if(selector=="all"){
		// add this event to all selectors
		for(auto sel = Selectors.begin(); sel!=Selectors.end(); ++sel){
			ret &= sel->second->AddPassingEvent(cut, rest...);
		}
	} else {
		ret = Selectors.at(selector)->AddPassingEvent(cut, rest...);
	}
	return ret;
}


#endif
