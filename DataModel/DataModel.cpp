#include "DataModel.h"
#include "ConnectionTable.h"
#include "skheadC.h" // for skheadg_.sk_geometry needed to construct the ConnectionTable
#include "fortran_routines.h"

DataModel* DataModel::thisptr=0;

DataModel::DataModel() : eventVariables_p(new BStore(false,true)), eventVariables(*eventVariables_p) {
	Log=0;
	// WriteOutput tool uses StoreToTTree to write these Stores to ROOT file.
	// the current implementation of StoreToTTree only works with Trees and Stores on the heap
	// but we have NTag tools that expect the members to be objects not pointers. TODO fix this.
	
	// make a global TApplication for ROOT drawing
	rootTApp = new TApplication("rootTApp",0,0);
	thisptr = this;
}

DataModel::~DataModel(){
	//if(rootTApp) delete rootTApp;                 // segfaults on application termination, maybe?
	//if(connectionTable) delete connectionTable;   // segfaults on application termination, maybe?
	eventParticles.clear();
	eventVertices.clear();
	NCaptureCandidates.clear();
	NCapturesTrue.clear();
	// finalise bonsai if it was initialised
	if(bonsai_initialised){
		cfbsexit_();
	}
}

// open a file, making it if necessary. Useful for adding data to the same file from multiple Tools.
TFile* DataModel::OpenFileForWriting(std::string file, bool alreadyopenonly){
	if(file==""){
		std::cerr<<"DataModel::OpenFileForWriting Error! called with empty file name!"<<std::endl;
		return nullptr;
	}
	// since we use the file path as a key we need to remove any ambiguity
	std::string file1;
	int retval = SystemCall(std::string("readlink -f ")+file, file1);
	if(retval!=0 || file1.empty()){
		std::cerr<<"DataModel::OpenFileForWriting Error! readlink '"<<file<<"' resolved to empty?"<<std::endl;
		return nullptr;
	}
	file=file1;
	TFile* outfile=nullptr;
	if(OutFiles.count(file)!=0){
		outfile = OutFiles.at(file).first;
		++OutFiles.at(file).second;
	} else {
		if(alreadyopenonly){
			std::cout<<"DataModel::OpenFileForWriting file "<<file<<" not found"<<std::endl;
			return nullptr;
		}
		outfile = new TFile(file.c_str(), "RECREATE");
		if(outfile==nullptr || outfile->IsZombie()){
			std::cerr<<"DataModel::OpenFileForWriting Error! Unable to open file "<<file<<std::endl;
			if(outfile){ outfile->Close(); delete outfile; }
			return nullptr;
		}
		OutFiles.emplace(file, std::pair<TFile*, int>{outfile,1});
	}
	return outfile;
}

// retrieve handle to a file already open by another Tool
TFile* DataModel::OpenFileForReading(std::string file, bool fromdisk){
	if(file==""){
		std::cerr<<"DataModel::OpenFileForReading Error! called with empty file name!"<<std::endl;
		return nullptr;
	}
	// since we use the file path as a key we need to remove any ambiguity
	file = GetStdoutFromCommand(std::string("readlink -f ")+file);
	TFile* outfile=nullptr;
	if(OutFiles.count(file)!=0){
		outfile = OutFiles.at(file).first;
		++OutFiles.at(file).second;
		return outfile;
	}
	if(!fromdisk){
		std::cerr<<"DataModel::OpenFileForReading Error! file "<<file<<" not found in file list!"<<std::endl;
		return nullptr;
	}
	return TFile::Open(file.c_str(), "READ");
}

bool DataModel::CloseFile(std::string file){
	// since we use the file path as a key we need to remove any ambiguity
	file = GetStdoutFromCommand(std::string("readlink -f ")+file);
	if(OutFiles.count(file)==0){
		std::cerr<<"DataModel::CloseFile Error! Unknown file "<<file<<std::endl;
		return false;
	} else {
		int users = --OutFiles.at(file).second;
		if(users==0){
			TFile* f = OutFiles.at(file).first;
			f->Close();
			delete f;
			OutFiles.erase(file);
		}
	}
	return true;
}

bool DataModel::CloseFile(TFile* fptr){
	if(fptr==nullptr){
		std::cerr<<"DataModel::CloseFile Error! called with null file pointer!"<<std::endl;
		return false;
	}
	std::string key=fptr->GetName();
	return CloseFile(key);
//	// scan for this file
//	for(auto&& afile : OutFiles){
//		if(afile.second.first==fptr){
//			int users = --afile.second.second;
//			if(users==0){
//				TFile* f = afile.second.first;
//				f->Close();
//				delete f;
//				OutFiles.erase(afile.first);
//			}
//			return true;
//		}
//	}
//	std::cerr<<"DataModel::CloseFile Error! Unknown file "<<fptr->GetName()<<std::endl;
//	return false;
}

ConnectionTable* DataModel::GetConnectionTable(int sk_geometry){
	if(connectionTable==nullptr){
		if(sk_geometry<0) sk_geometry = skheadg_.sk_geometry;
		if(sk_geometry<1) sk_geometry = 6;  // default
		// as of 2023-05-09 the ConnectionTable constructor only accepts
		// a value of 5 for SK-5 geometry; anything else results in SK-4.
		// see $SKOFL_ROOT/src/ConnectionTableReader/ConnectionTable.cc
		if(sk_geometry!=4 && sk_geometry!=5){
			sk_geometry = (sk_geometry<4) ? 4 : 5;
			std::cerr<<"DataModel::GetConnectionTable: Warning! ConnectionTable only supports"
			         <<"SK-4 or SK-5; using gometry SK-"<<sk_geometry<<std::endl;
		}
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

bool DataModel::RegisterReader(std::string readerName, MTreeReader* reader, std::function<bool()> hasAFT, std::function<bool()> loadSHE, std::function<bool()> loadAFT, std::function<bool(int)> loadCommon, std::function<int(long)> getTreeEntry){
	Trees.emplace(readerName, reader);
	hasAFTs.emplace(readerName, hasAFT);
	loadSHEs.emplace(readerName, loadSHE);
	loadAFTs.emplace(readerName, loadAFT);
	loadCommons.emplace(readerName, loadCommon);
	getEntrys.emplace(readerName, getTreeEntry);
	return true;
}

MTreeReader* DataModel::GetTreeReader(){
	// most of the time there's one primary TreeReader associated with an SKROOT file which is
	// populating common blocks, and possibly further TreeReaders for auxiliary info we can identify
	// the primary TreeReader as its LUN will be set in the common block variable skheadf_.root_id
	int file_lun = skheadf_.root_id;
	std::map<std::string,int> lunlist;
	CStore.Get("LUNList",lunlist);
	for(auto&& apair : lunlist){
		if(apair.second==file_lun){
			if(Trees.count(apair.first)) return Trees.at(apair.first);
			break;
		}
	}
	std::cerr<<"DataModel::GetTreeReader no reader associated with skheadf_.root_id LUN "
	         <<file_lun<<std::endl;
	return nullptr;
}

int DataModel::getTreeEntry(std::string ReaderName, long entrynum){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(getEntrys.size()) ReaderName = getEntrys.begin()->first;
	}
	if(getEntrys.count(ReaderName)){
		// check the user registered a valid pointer to this method.
		// Typical situation where this may not be the case is if the MTreeReader is associated
		// to a Tree which is being actively written by an upstream Tool, so we can't have users
		// switching events
		if(getEntrys.at(ReaderName)){
			return getEntrys.at(ReaderName)(entrynum);
		} else {
			std::cerr<<"DataModel::getTreeEntry is not available for treeReader "<<ReaderName
			         <<" (perhaps this Tree is being generated by an upstream Tool?)"<<std::endl;
			return false;
		}
	} else {
		std::cerr << "DataModel::getTreeEntry requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return 0;
}

bool DataModel::HasAFT(std::string ReaderName){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(hasAFTs.size()) ReaderName = hasAFTs.begin()->first;
	}
	if(hasAFTs.count(ReaderName)){
		// check the user registered a valid pointer to this method.
		// Typical situation where this may not be the case is if the MTreeReader is not
		// associated with a TreeReader Tool, and is just a standalone MTreeReader class instance.
		if(hasAFTs.at(ReaderName)){
			return hasAFTs.at(ReaderName)();
		} else {
			std::cerr<<"DataModel::HasAFT is not available for treeReader "<<ReaderName<<std::endl;
			return false;
		}
	} else {
		std::cerr<<"DataModel::HasAFT requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadSHE(std::string ReaderName){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(loadSHEs.size()) ReaderName = loadSHEs.begin()->first;
	}
	if(loadSHEs.count(ReaderName)){
		// check the user registered a valid pointer to this method.
		// Typical situation where this may not be the case is if the MTreeReader is not
		// associated with a TreeReader Tool, and is just a standalone MTreeReader class instance.
		if(loadSHEs.at(ReaderName)){
			return loadSHEs.at(ReaderName)();
		} else {
			std::cerr<<"DataModel::LoadSHE is not available for treeReader "<<ReaderName<<std::endl;
			return false;
		}
	} else {
		std::cerr<<"DataModel::LoadSHE requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

bool DataModel::LoadAFT(std::string ReaderName){
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(loadAFTs.size()) ReaderName = loadAFTs.begin()->first;
	}
	if(loadAFTs.count(ReaderName)){
		// check the user registered a valid pointer to this method.
		// Typical situation where this may not be the case is if the MTreeReader is not
		// associated with a TreeReader Tool, and is just a standalone MTreeReader class instance.
		if(loadAFTs.at(ReaderName)){
			return loadAFTs.at(ReaderName)();
		} else {
			std::cerr<<"DataModel::LoadAFT is not available for treeReader "<<ReaderName<<std::endl;
			return false;
		}
	} else {
		std::cerr<<"DataModel::LoadAFT requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

// TODO... remove this? Its use is sketchy. So far nothing uses it.
bool DataModel::LoadCommons(int entry_i, std::string ReaderName){
	// to clarify what this method is for; when reading skroot files the TreeReader Tool
	// may be configured to read multiple entries at once, caching the common block contents
	// from each entry in an internal vector. Downstream Tools may then call LoadCommons(n)
	// to populate the ordinary, global common blocks using a given cached entry,
	// thereby effectively giving access to multiple events in one Execute call.
	// HOWEVER: LoadCommons ONLY fills SOME common blocks (as many as i'm aware of, basically)
	// and DOES NOT load the corresponding TTree entry, so branch access from the corresponding
	// MTreeReader will NOT be correct! XXX USE WITH CAUTION XXX
	if(ReaderName==""){
		// if no name given but we have only one TreeReader Tool, use that
		if(loadCommons.size()) ReaderName = loadCommons.begin()->first;
	}
	if(loadCommons.count(ReaderName)){
		// check the user registered a valid pointer to this method.
		// Typical situation where this may not be the case is if the MTreeReader is not
		// associated with a TreeReader Tool, and is just a standalone MTreeReader class instance.
		if(loadCommons.at(ReaderName)){
			return loadCommons.at(ReaderName)(entry_i);
		} else {
			std::cerr<<"DataModel::LoadEntry is not available for treeReader "<<ReaderName<<std::endl;
			return false;
		}
	} else {
		std::cerr<<"DataModel::LoadEntry requested for Unknown reader "<<ReaderName<<std::endl;
	}
	return false;
}

void DataModel::KZInit(){
	// wrapper to make sure this only gets invoked once
	if(!kz_initialized){
		kz_initialized=true;
		kzinit_();
	}
}

bool DataModel::GeoSet(int sk_geometry_in){
	// wrapper to make sure this only gets invoked once
	if(sk_geometry_in<=0){
		std::cerr<<"Error! DataModel::GeoSet called with invalid SK Geometry: "<<sk_geometry_in
		         <<std::endl;
		return false;
	}
	if(skheadg_.sk_geometry>0 && sk_geometry_in!=skheadg_.sk_geometry){
		std::cerr<<"*************"
		         <<"!!!WARNING!!!"
		         <<"*************"
		         <<"DataModel::GeoSet called with SK Geometry "<<sk_geometry_in
		         <<" but geometry has already been set to "
		         <<skheadg_.sk_geometry<<"\nGeometry will NOT be changed!"<<std::endl;
		//exit(-1); // should we do this? printouts are easy to miss...
		return false;
	} else {
		skheadg_.sk_geometry = sk_geometry_in;
		geoset_();
	}
	return true;
}

bool DataModel::BonsaiInit(){
	if(!bonsai_initialised){
		// initialize bonsai
		int MAXPM_var = MAXPM;
		float* xyzpm = &geopmt_.xyzpm[0][0];
		cfbsinit_(&MAXPM_var, xyzpm);
		bonsai_initialised=true;
		return true;
	}
	return false;
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

