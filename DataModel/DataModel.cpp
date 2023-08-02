#include "DataModel.h"
#include "ConnectionTable.h"
#include "skheadC.h" // for skheadg_.sk_geometry needed to construct the ConnectionTable
#include "fortran_routines.h"

DataModel::DataModel() : eventVariables_p(new BStore(false,true)), eventVariables(*eventVariables_p) {
	Log=0;
	// WriteOutput tool uses StoreToTTree to write these Stores to ROOT file.
	// the current implementation of StoreToTTree only works with Trees and Stores on the heap
	// but we have NTag tools that expect the members to be objects not pointers. TODO fix this.
	
	// make a global TApplication for ROOT drawing
	rootTApp = new TApplication("rootTApp",0,0);
}

DataModel::~DataModel(){
	//if(rootTApp) delete rootTApp; // segfaults on application termination. J. Fannon 
	//if(connectionTable) delete connectionTable;   // segfaults on application termination, maybe?
}

/*
TTree* DataModel::GetTTree(std::string name){

  return m_trees[name];

}


void DataModel::AddTTree(std::string name,TTree *tree){

  m_trees[name]=tree;

}


void DataModel::DeleteTTree(std::string name,TTree *tree){

  m_trees.erase(name);

}

*/

ConnectionTable* DataModel::GetConnectionTable(int sk_geometry){
	if(connectionTable==nullptr){
		if(sk_geometry<0) sk_geometry = skheadg_.sk_geometry;
		if(sk_geometry<1) sk_geometry = 6;  // default
		connectionTable = new ConnectionTable(sk_geometry);
	}
	return connectionTable;
}

TApplication* DataModel::GetTApp(){
	if(rootTApp==nullptr){
		rootTApp = new TApplication("rootTApp",0,0);
	}
	return rootTApp;
}

bool DataModel::RegisterReader(std::string readerName, MTreeReader* reader, std::function<bool()> hasAFT, std::function<bool()> loadSHE, std::function<bool()> loadAFT, std::function<bool(int)> loadCommon, std::function<int(long, bool)> getTreeEntry){
	Trees.emplace(readerName, reader);
	hasAFTs.emplace(readerName, hasAFT);
	loadSHEs.emplace(readerName, loadSHE);
	loadAFTs.emplace(readerName, loadAFT);
	loadCommons.emplace(readerName, loadCommon);
	getEntries.emplace(readerName, getTreeEntry);
	return true;
}

int DataModel::getTreeEntry(std::string ReaderName, long entrynum){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(hasAFTs.size()) ReaderName = hasAFTs.begin()->first;
	}
	if(getEntries.count(ReaderName)){
		return getEntries.at(ReaderName)(entrynum, true);
	} else {
		std::cerr << "getTreeEntry requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return 0;
}

bool DataModel::HasAFT(std::string ReaderName){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(hasAFTs.size()) ReaderName = hasAFTs.begin()->first;
	}
	if(hasAFTs.count(ReaderName)){
		return hasAFTs.at(ReaderName)();
	} else {
		std::cerr<<"HasAFT requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadSHE(std::string ReaderName){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(hasAFTs.size()) ReaderName = hasAFTs.begin()->first;
	}
	if(loadSHEs.count(ReaderName)){
		return loadSHEs.at(ReaderName)();
	} else {
		std::cerr<<"LoadSHE requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadAFT(std::string ReaderName){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(hasAFTs.size()) ReaderName = hasAFTs.begin()->first;
	}
	if(loadAFTs.count(ReaderName)){
		return loadAFTs.at(ReaderName)();
	} else {
		std::cerr<<"LoadAFT requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadEntry(int entry_i, std::string ReaderName){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(hasAFTs.size()) ReaderName = hasAFTs.begin()->first;
	}
	if(loadCommons.count(ReaderName)){
		return loadCommons.at(ReaderName)(entry_i);
	} else {
		std::cerr<<"LoadEntry requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

void DataModel::KZInit(){
	if(!kz_initialized){
		kz_initialized=true;
		kzinit_();
	}
}

// return a new LUN. We accept a hint, but will only apply it if not already assigned.
int DataModel::GetNextLUN(int lun, std::string reader){
	// each LUN (logic unit number, a fortran file handle (ID) and/or an ID used
	// by the SuperManager to identify the TreeManager associated with a file) must be unique.
	std::map<std::string,int> lunlist;
	CStore.Get("LUNList",lunlist);
	// check if this LUN is free, otherwise print a warning and assign a different LUN
	// since we map names to LUNs, to check if a LUN is free we need to scan by value, not by key.
	// easiest way to do this while also sorting by LUN is to reverse the map.
	// sorting allows us to immediately know the next free LUN in case this one is in use.
	std::map<int,std::string> revlist;
	for(auto&& apair : lunlist){ revlist.emplace(apair.second,apair.first); }
	if(revlist.count(lun)){
		int reqLUN=lun;
		lun = revlist.rbegin()->first; // get the last assigned LUN
		++lun;                         // we'll use the next one
		std::cerr<<"DataModel::GenerateNewLUN Warning! Cannot assign LUN "<<reqLUN
			 <<" as it is already taken. Assigning "<<lun<<" instead."<<std::endl;
	}
	//assign the LUN
	lunlist.emplace(reader,lun);
	CStore.Set("LUNList",lunlist);
	return lun;
}

int DataModel::GetLUN(std::string reader){
	// get the LUN associated with a given TreeReader via its name.
	// can be useful e.g. to get the associated TreeManager
	std::map<std::string,int> lunlist;
	CStore.Get("LUNList",lunlist);
	if(lunlist.count(reader)){
		return lunlist.at(reader);
	} else {
		return -1;
	}
}

bool DataModel::FreeLUN(int lun, std::string reader){
	
	// get the map of open TreeManagers
	std::map<std::string,int> lunlist;
	CStore.Get("LUNList",lunlist);
	
	if(reader!=""){
		// reader name given, try to find the LUN associated to this reader
		if(lunlist.count(reader)==0){
			std::string errmsg = "DataModel::FreeLUN error! Could not find LUN list entry for reader ";
			errmsg += reader;
			// if we were given a LUN number, see if we recognise that instead
			// we won't close on it though, as we're not sure it's the right one....
			if(lun>0){
				// we need to build the reverse map first
				std::map<int,std::string> revlist;
				for(auto&& apair : lunlist){ revlist.emplace(apair.second,apair.first); }
				// now check if we know this LUN number
				if(revlist.count(lun)){
					errmsg += ". Specified LUN number ("+std::to_string(lun)
					       +") is associated with reader " + revlist.at(lun)
					       +". No LUNs freed.";
				} else {
					errmsg += ". Nor was there any entry matching given LUN number "
					       + std::to_string(lun);
				}
			}
			std::cerr<<errmsg<<std::endl;
			return false;
		// else we do have a reader by this name. If we were given a LUN number too, check it matches.
		} else if(lun>0 && lunlist.at(reader)!=lun){
			// uh-oh, mismatch.
			std::cerr<<"DataModel::FreeLUN error! LUN value in lunlist ("<<lunlist.at(reader)
			         <<") does not match the given reader "<<reader
			         <<"! No LUNs freed."<<std::endl;
			return false;
		} else {
			// else either no LUN number given, or it matches the reader name. Proceed.
			std::cout<<"Freeing LUN "<<lunlist.at(reader)<<" associated to reader "
			         <<reader<<std::endl;
			lunlist.erase(reader);
			return true;
		}
	} else {
		// no reader name given. to erase by LUN number, we need to build the reverse map
		std::map<int,std::string> revlist;
		for(auto&& apair : lunlist){ revlist.emplace(apair.second,apair.first); }
		// check we know this LUN number
		if(revlist.count(lun)){
			std::cout<<"Freeing LUN "<<lun<<" associated to reader "
			         <<revlist.at(lun)<<std::endl;
			lunlist.erase(revlist.at(lun));
			return true;
		} else {
			std::cerr<<"DataModel::FreeLUN error! LUN value given ("<<lun
			    <<") not found in DataModel lun list!"<<std::endl;
			return false;
		}
	}
	// dummy
	return false;
}

