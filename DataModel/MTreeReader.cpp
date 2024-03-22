/* vim:set noexpandtab tabstop=4 wrap */
#include "MTreeReader.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLeafElement.h"
//#include "TParameter.h"

#include "type_name_as_string.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm> // std::find

#include "Algorithms.h"  // CheckPath

bool Notifier::Notify(){
	if(verbosity) std::cout<<"Notifier for "<<treeReader->GetName()<<" loading new TTree"<<std::endl;
	//treeReader->GetTree()->Show();
	treeReader->GetTree()->GetTree()->GetEntry(0);
	//treeReader->GetTree()->Show();
	//return treeReader->UpdateBranchPointers(true);  // not sufficient!
	return treeReader->ParseBranches();               // probably overkill - TODO happy medium
}

// TODO constructor/loader for tchains or tree pointers

MTreeReader::MTreeReader(std::string iname, std::string fpath, std::string treename){
	name = iname;
	Load(fpath, treename);
}

MTreeReader::MTreeReader(std::string iname) : name(iname){};

std::string MTreeReader::GetName(){
	return name;
}

void MTreeReader::SetName(std::string iname){
	name=iname;
}

int MTreeReader::Load(std::string filename, std::string treename){
	// determine if filename is an actual file or a pattern (TChain)
	std::string pathtype;
	bool pathexists = CheckPath(filename, pathtype);
	int ok=1;
	if(pathexists && pathtype=="f"){
		// given a file: load it as normal
		ok = LoadFile(filename);
		if(not ok) return ok;
		LoadTree(treename);
		if(not ok) return ok;
		// preload first entry to get branches
		if(verbosity) std::cout<<"getting entry 0"<<std::endl;
		thetree->GetEntry(0);
		ok = ParseBranches();
		return ok;
	} else if(pathexists && pathtype=="d"){
		// given a directory: we can't work with just this
		std::cerr<<"!!! MTreeReader constructor called with a path to a directory !!!"<<std::endl
				 <<"Please pass either a path to a file or a glob pattern"<<std::endl
				 <<"The passed path was "<<filename<<std::endl;
	} else if(pathexists){
		// stat indicated that this path exists, but was neither path nor directory?
		std::cerr<<"!!! MTreeReader constructor called with a path to an unknown resource type !!!"<<std::endl
				 <<"Please pass either a path to a file or a glob pattern"<<std::endl
				 <<"The passed path was "<<filename<<", type was "<<pathtype<<std::endl;
	} else {
		// path does not exist - it could be a glob pattern
		// try to construct a TChain from it
		TChain* chain = new TChain(treename.c_str());
		if(verbosity) std::cout<<"loading TChain '"<<treename<<"' with files "<<filename<<std::endl;
		int filesadded = chain->Add(filename.c_str());
		// ↑ note this does not check the files contain the correct TTree!
		if(filesadded==0){
			std::cerr<<"!!! MTreeReader constructor called with a non-existant input path !!!"<<std::endl
					 <<"An attempt to use this as a pattern to construct a TChain found no files!"<<std::endl
					 <<"The passed path was "<<filename<<std::endl;
		} else {
			int localEntry = chain->LoadTree(0);
			if(localEntry<0){
				std::cerr<<"!!! MTreeReader constructor called with a non-existant input path !!!"<<std::endl
						 <<"An attempt to use this as a pattern to construct a TChain found "
						 <<filesadded<<" files, but loading the first TTree failed with code "
						 <<localEntry<<"!"<<std::endl
						 <<"Is the tree name \""<<treename<<"\" correct?"<<std::endl
						 <<"The passed path was "<<filename<<std::endl;
			} else {
				Load(chain);
			}
		}
	}
	return ok;
}

int MTreeReader::Load(std::vector<std::string> filelist, std::string treename){
	// Construct a TChain and add all files from the list
	if(verbosity) std::cout<<"loading TChain '"<<treename<<"' with "<<filelist.size()<<" files"<<std::endl;
	if(filelist.size()==0){
		std::cerr<<"!!! MTreeReader constructor called with empty file list !!!"<<std::endl;
		return -1;
	}
	TChain* chain = new TChain(treename.c_str());
	for(std::string& afile : filelist){
		chain->Add(afile.c_str());
		// ↑ note this does not check the files contain the correct TTree!
	}
	int localEntry = chain->LoadTree(0);
	if(localEntry<0){
		std::cerr<<"!!! MTreeReader constructor found no valid files !!!"<<std::endl
				 <<"An attempt to use construct a TChain scanned "<<filelist.size()
				 <<" files, but loading the first TTree failed with code "
				 <<localEntry<<"!"<<std::endl
				 <<"Is the tree name \""<<treename<<"\" correct?"<<std::endl
				 <<"The first path was "<<filelist.front()<<std::endl;
	} else {
		Load(chain);
	}
	return 1;
}

int MTreeReader::Load(TTree* thetreein){
	thetree = thetreein;
	
	// preload first entry to get branches
	if(verbosity) std::cout<<"getting entry 0"<<std::endl;
	thetree->GetEntry(0);
	
	int ok = ParseBranches();
	// do this after ParseBranches as TChains may return nullptr if no file has been loaded yet
	thefile = thetree->GetCurrentFile();
	notifier.SetReader(this);
	thetree->SetNotify(&notifier);
	return ok;
}

int MTreeReader::LoadFile(std::string filename){
	if(verbosity) std::cout<<"getting file "<<filename<<std::endl;
	thefile = TFile::Open(filename.c_str());
	if(not thefile){
		std::cerr<<"MTreeReader failed to load file "<<filename<<std::endl;
		return 0;
	}
	return 1;
}

int MTreeReader::LoadTree(std::string treename){
	if(verbosity) std::cout<<"getting tree "<<treename<<std::endl;
	thetree = (TTree*)thefile->Get(treename.c_str());
	if(not thetree){
		std::cerr<<"MTreeReader could not find tree "<<treename<<" in file "
				 <<thetree->GetCurrentFile()->GetName()<<std::endl;
		return 0;
	}
	// check the tree has entries
	if(verbosity) std::cout<<"getting num entries"<<std::endl;
	auto numentries = thetree->GetEntriesFast();
	if(numentries==0){
		std::cerr<<"MTreeReader tree "<<thetree->GetName()<<" in file "<<thetree->GetCurrentFile()->GetName()
				 <<" has no entries"<<std::endl;
		return 0;
	}
	if(verbosity) std::cout<<"tree has "<<numentries<<" entries"<<std::endl;
	return 1;
}

int MTreeReader::ParseBranches(){
	
	branch_pointers.clear();
	branch_istobject.clear();
	leaf_pointers.clear();
	branch_titles.clear();
	branch_value_pointers.clear();
	branch_types.clear();
	branch_isobject.clear();
	branch_isobjectptr.clear();
	branch_isarray.clear();
	branch_dimensions.clear();
	branch_dims_cache.clear();
	
	if(verbosity) std::cout<<name<<" updating branch addresses"<<std::endl;
	//for(int i=0; i<thetree->GetListOfLeaves()->GetEntriesFast(); ++i){
	for(int i=0; i<thetree->GetListOfBranches()->GetEntriesFast(); ++i){
		//TLeaf* lf=(TLeaf*)thetree->GetListOfLeaves()->At(i);
		TBranch* br=(TBranch*)thetree->GetListOfBranches()->At(i);
		// trying to work around splitting of branches here, so that we parse
		// the list of branches, but without processing branches that represent
		// members of a split class. In a split branch, the branch from
		// 'TTree->GetListOfBranches()' has itself 'TBranch->GetListOfBranches > 0'
		// whereas non-split branches have ==0. But, both seem to have
		// 'TBranch->GetListOfLeaves()==1'. .. are we safe to assume always 1?
		TLeaf* lf = (TLeaf*)br->GetListOfLeaves()->At(0);
		std::string branchname = lf->GetName();
		leaf_pointers.emplace(branchname,lf);
		branch_pointers.emplace(branchname,lf->GetBranch());
		// the following is fine for objects, primitives or containers
		// but only returns the primitive type for c-style arrays
		branch_types.emplace(branchname,lf->GetTypeName());
		// the branch title includes dimensionality for c-style arrays
		// e.g. "mybranchname     int[nparts][3]/F"
		std::string branchtitle=lf->GetBranch()->GetTitle();
		branch_titles.emplace(branchname,branchtitle);
		
		// handle object pointers
		if (lf->IsA() == TLeafElement::Class()) {
			// could be TObjects, or could be stl containers
			// is intptr_t any better than (void*)? probably not.
			// note that 'lf->GetValuePointer()' returns 0 for objects!
			TBranchElement* bev = (TBranchElement*)lf->GetBranch();
			intptr_t objpp=reinterpret_cast<intptr_t>(bev->GetObject());
			branch_value_pointers.emplace(branchname,objpp);
			branch_isobject.emplace(branchname,true);
			branch_isobjectptr.emplace(branchname,false);
			branch_isarray.emplace(branchname,false);
			
			// both classes inheriting from TObject and STL containers
			// are flagged as TLeafElements, so need further check
			TClass* ac = TClass::GetClass(lf->GetTypeName());
			if(ac!=nullptr){
				branch_istobject.emplace(branchname,ac->InheritsFrom("TObject"));
			} else {
				std::cerr<<"Unknown class for branch "<<lf->GetBranch()->GetTitle()
						 <<"! Please make a dictionary."<<std::endl;
				// what are the consequences of this? Do we need to remove it from
				// the maps? Will the TreeReader fail if we do not?
				// is it sufficient simply to instead just say it does not inherit from TObject?
				branch_istobject.emplace(branchname,false);
			}
		} else if( (strcmp(br->ClassName(),"TBranchObject")==0) ||
		           (strcmp(lf->GetTypeName(),"TClonesArray")==0) ){
			// some classes, such as TVector3, are stored as TBranchObjects
			// rather than TBranchElements. Similar to basic primitives
			// we need to get the pointer from TLeaf::GetValuePointer(),
			// but similar to TBranchElement, TBranchObject::GetAddress() returns
			// a pointer to a pointer to the object, rather than a pointer to the object itself.
			// (TBranchElement::GetObject performs one dereference, returning a
			// pointer to the object directly, but is not a method for TBranchObjects)
			// The handling of TClonesArrays is the same; although they identify as TBranchElements,
			// TBranchElement::GetObject returns nullptr
			void** ptrptr=reinterpret_cast<void**>(lf->GetValuePointer());
			//intptr_t objpp = reinterpret_cast<intptr_t>(*ptrptr);
			// don't dereference it: it seems like this can change from entry to entry!
			intptr_t objpp = reinterpret_cast<intptr_t>(ptrptr);
			// n.b. we can use TClonesArray->GetClass()->GetName() to get the class name of held elements
			branch_value_pointers.emplace(branchname,objpp);
			branch_isobject.emplace(branchname,false);
			branch_isobjectptr.emplace(branchname,true);
			branch_isarray.emplace(branchname,false);
			// not sure if we can assume that objects in a TBranchObject derive from TObject...maybe?
			TClass* ac = TClass::GetClass(lf->GetTypeName());
			if(ac!=nullptr){
				branch_istobject.emplace(branchname,ac->InheritsFrom("TObject"));
			} else {
				branch_istobject.emplace(branchname,false);
			}
		}
		// handle arrays
		//else if(lf->GetLen()>1){  // flattened length. Unsuitable when dynamic size happens to be 1!
		else if(branchtitle.find_first_of("[",0)!=std::string::npos){  // hope for no '[' in branch names
			// we'll need to parse the title to retrieve the actual dimensions
			intptr_t objpp=reinterpret_cast<intptr_t>(lf->GetValuePointer());
			branch_value_pointers.emplace(branchname,objpp);
			branch_isobject.emplace(branchname,false);
			branch_isobjectptr.emplace(branchname,false);
			branch_isarray.emplace(branchname,true);
			branch_istobject.emplace(branchname,false);
		}
		// handle basic types
		else {
			intptr_t objpp=reinterpret_cast<intptr_t>(lf->GetValuePointer());
			branch_value_pointers.emplace(branchname,objpp);
			branch_isobject.emplace(branchname,false);
			branch_isobjectptr.emplace(branchname,false);
			branch_isarray.emplace(branchname,false);
			branch_istobject.emplace(branchname,false);
		}
	}
	
	// for all array branches, parse their titles to extract information about dimensions
	for(auto&& abranch : branch_isarray){
		if(abranch.second) ParseBranchDims(abranch.first);
	}
	
	// build a string list of branch names
	branchnamestring="";
	for(auto&& abranch : branch_titles){
		if(!branchnamestring.empty()) branchnamestring += ", ";
		branchnamestring += abranch.first;
	}
	branchnamestring = "{" + branchnamestring + "}";
	
	return 1;
}

int MTreeReader::UpdateBranchPointer(std::string branchname){
	if(leaf_pointers.count(branchname)==0) return 0;
	TLeaf* lf = leaf_pointers.at(branchname);
	intptr_t objpp;
	if(branch_isobject.at(branchname)){
		// for objects TLeaf::GetValuePointer returns a pointer to a pointer
		// to the held object. We can remove the additional level of indirection
		// by calling TBranchElement::GetObject instead
		TBranchElement* bev = (TBranchElement*)lf->GetBranch();
		objpp=reinterpret_cast<intptr_t>(bev->GetObject());
	} else if(branch_isobjectptr.at(branchname)){
		// for some classes such as TClonesArrays and TVector3s,
		// TLeafElement::GetValuePointer returns a pointer to a pointer,
		// but the parent branch is a TBranchObject not TBranchElement
		// and does not support a GetObject method.
		// In such cases we need to do the additional dereference ourselves
		void** ptrptr=reinterpret_cast<void**>(lf->GetValuePointer());
		objpp = reinterpret_cast<intptr_t>(*ptrptr);
	} else {
		// for simple types TLeaf::GetValuePointer returns
		// a pointer to the object itself.
		objpp=reinterpret_cast<intptr_t>(lf->GetValuePointer());
	}
	branch_value_pointers.at(branchname)=objpp;
	return 1;
}

int MTreeReader::UpdateBranchPointers(bool all){
	// dynamic arrays and/or objects may move when switching files in a TChain
	for(auto&& abranch : branch_isarray){
		if(all || abranch.second){
			// fixed sized arrays won't need to be reallocated
			// so only update pointers if we don't have a cached (constant) size
			if(all || branch_dims_cache.count(abranch.first)==0){
				int ok = UpdateBranchPointer(abranch.first);
				if(not ok) return ok;
			}
		}
	}
	return 1;
}

int MTreeReader::ParseBranchDims(std::string branchname){
	// parse branch title for sequences of type '[X]' suggesting an array.
	// extract 'X'. Scan the list of branch names for 'X', in which case
	// this is a variable length array, otherwise it's a fixed size so use stoi.
	std::string branchtitle = branch_titles.at(branchname);
	size_t startpos = 0;
	size_t endpos = 0;
	// each subsequent entry of the array represents a dimension
	// the pair will either hold a branch name (in the first element)
	// or a numeric size (in the second element)
	std::vector<std::pair<std::string,int>> this_branch_dimensions;
	while(true){
		startpos = branchtitle.find_first_of("[", endpos+1);
		if(startpos==std::string::npos) break;
		endpos = branchtitle.find_first_of("]",startpos+1);
		if(endpos==std::string::npos){
			std::cerr<<"Trailing '[' character while parsing branch titles for branch "
					 <<branchname<<", title "<<branchtitle<<"!"<<std::endl;
			return 0;
		}
		std::string sizestring = branchtitle.substr(startpos+1,endpos-startpos-1);
		//std::cout<<"extracted array size label "<<sizestring<<" for branch "<<branchname
		//		 <<" dimension "<<this_branch_dimensions.size()<<std::endl;
		// check if this string is the name of another branch - i.e. variable size array
		if(branch_pointers.count(sizestring)){
			// it is - we'll need to retrieve the corresponding branch entry to get the size on each entry
			this_branch_dimensions.emplace_back(std::pair<std::string, int>{sizestring,0});
		} else {
			// it ought to be a static size, try to convert from string to int
			// catch the exception thrown in the event that it can't be converted
			// note this only falls over when the start of the string isn't numeric
			// it copes with negatives, but truncates to integer, and ignores anthing
			// after the first non-digit character.
			try{
				this_branch_dimensions.emplace_back(std::pair<std::string,int>{"",std::stoi(sizestring)});
			}
			catch(const std::invalid_argument& ia) {
				std::cerr<<"Failed to extract array dimension from branch title "<<branchtitle
						 <<" - not a recognised branchname and string to int conversion failed with "
						 <<ia.what()<<std::endl;
				return 0;
			}
			catch(const std::out_of_range& oor){
				std::cerr<<"Failed to extract array dimension from branch title "<<branchtitle
						 <<" - not a recognised branchname and string to int conversion failed with "
						 <<oor.what()<<std::endl;  // occurs if string exceeds range representable by int
				return 0;
			}
		}
		// loop back round for any further dimensions
	}
	branch_dimensions.emplace(branchname,this_branch_dimensions);
	if(branch_dimensions.at(branchname).size()==0){
		std::cerr<<"Failed to identify any dimensions for branch "<<branchname
				 <<" despite TLeaf::GetLength() returning >1"<<std::endl;
		return 0;
	}
	
	int vlevel=2;
	if(verbosity>vlevel) std::cout<<"end of branch title parsing, found "
								  <<branch_dimensions.at(branchname).size()<<" dimensions, [";
	// loop over the vector of dimensions
	std::string dims_string = "";
	bool allstatics=true;
	std::vector<size_t> cached_dims;
	for(auto&& dims_pair : branch_dimensions.at(branchname)){
		if(dims_pair.first==""){
			if(verbosity>vlevel) std::cout<<"(N)"<<dims_pair.second;  // static numeric size
			cached_dims.push_back(dims_pair.second);
			dims_string.append(std::string("[") + std::to_string(dims_pair.second) + std::string("]"));
		} else {
			if(verbosity>vlevel) std::cout<<"(B)"<<dims_pair.first; // name of branch that stores variable size
			dims_string.append(std::string("[") + dims_pair.first + std::string("]"));
			allstatics=false;
		}
		if((verbosity>vlevel)&&(dims_pair!=(branch_dimensions.at(branchname).back()))) std::cout<<"], [";
	}
	// if all dimensions are static we can cache the results for quicker lookup
	if(allstatics){ branch_dims_cache.emplace(branchname,cached_dims); }
	// append dimensions to the type string, since they aren't properly indicated by TLeaf::GetTypeName
	if(branch_types.count(branchname)){
		branch_types.at(branchname).append(dims_string);
	}
	
	return 1;
}

std::vector<size_t> MTreeReader::GetBranchDims(std::string branchname){
	// get the dimensions of the array for this entry
	// if all dimensions are constant we should have them cached
	if(branch_dims_cache.count(branchname)) return branch_dims_cache.at(branchname);
	
	// otherwise we must retrieve at least one size from another branch
	std::vector<size_t> dimstemp;  // dimensions for this entry
	if(branch_dimensions.count(branchname)==0){
		std::cerr<<"GetBranchDims called but no dimensions for this branch!"<<std::endl;
		return dimstemp;
	}
	// loop over dimensions
	for(auto&& adim : branch_dimensions.at(branchname)){
		if(adim.first==""){
			// this dimension is constant
			dimstemp.push_back(adim.second);
		} else {
			// this dimensions is a branch name - get the entry value of that branch
			std::string sizebranchname = adim.first;
			int lengththisentry;
			int get_ok = GetBranchValue(sizebranchname, lengththisentry);
			if(not get_ok){
				std::cerr<<"Failed to retrieve value for branch "<<sizebranchname
						 <<" required while obtaining this entry's dimensions for array in branch "
						 <<branchname<<std::endl;
				return dimstemp; // TODO throw an exception? We should do more than just print an error...
			}
			dimstemp.push_back(lengththisentry);
		}
	}
	return dimstemp;
}

int MTreeReader::Clear(){
	// loop over all branches
	for(auto&& isobject : branch_istobject){
		// skip if doesn't inherit from TObject so may not have Clear() method
		// XXX note, maybe we should check if it has a 'clear' method (stl container)
		// and invoke that if not? Should be safe even without doing that though.
		if(not isobject.second) continue;
		// get pointer to the object otherwise
		TObject* theobject = nullptr;
		if(branch_isobjectptr.at(isobject.first)){
			TObject** ptrptr = reinterpret_cast<TObject**>(branch_value_pointers.at(isobject.first));
			if(not ptrptr){
				std::cerr<<"MTreeReader AutoClear error: failure to get pointer to TObject* "
						 <<"for branch "<<isobject.first<<std::endl;
				continue;  // TODO throw suitable exception
			}
			theobject = *ptrptr;
		} else {
			theobject = reinterpret_cast<TObject*>(branch_value_pointers.at(isobject.first));
		}
		if(not theobject){
			// no object... is this an error?
			std::cerr<<"MTreeReader AutoClear error: failure to get pointer to TObject "
					 <<"for branch "<<isobject.first<<std::endl;
			continue;  // TODO throw suitable exception
		}
		if(verbosity>4) std::cout<<"Clearing "<<isobject.first<<std::endl;
		theobject->Clear();
		if(verbosity>10) std::cout<<"done"<<std::endl;
	}
	return 1; // TODO check for errs
}

int MTreeReader::GetEntry(long entry_number, bool skipTreeRead){
	// in case we've already got this entry loaded, nothing to do
	if(currentEntryNumber==entry_number) return 1;
	
	// if we've been requested to invoke Clear() on all objects before each Get, do so
	if(verbosity>3) std::cout<<"MTreeReader GetEntry "<<entry_number<<std::endl;
	if(autoclear){
		int clear_ok = Clear();
		if(not clear_ok){ return -10; }
	}
	
	// if we're processing a chain, load the tree first
	int status = thetree->LoadTree(entry_number);
	if(status<0){
		/*
		-1: The chain is empty.
		-2: The requested entry number is less than zero or too large for the chain or TTree.
		-3: The file corresponding to the entry could not be correctly opened
		-4: The TChainElement corresponding to the entry is missing or the TTree is missing from the file.
		-5: Internal error, please report the circumstance when this happen as a ROOT issue.
		-6: An error occurred within the notify callback.
		*/
		//  0 from GetEntry indicates entry does not exist     |
		// -1 from LoadTree indicates empty tchain             | these 3 are roughly equivalent
		// -2 from LoadTree indicates <0 or off end of chain   |
		
		// -1 from GetEntry indicates IO error                 | these 2 (5) are roughly equivalent
		// -3->-6 from LoadTree indicate IO errors             | 
		// first set: terminate toolchain, second set: maybe try to continue to next entry?
		if(status<-2){
			std::cerr<<"MTreeReader error loading next TTree from TChain! "
					 <<"TChain::LoadTree returned "<<status<<"\n";
		}
		int bytesread;
		if(status>-3) bytesread = 0;  // treat as "end of tchain" - merge "empty tchain" with "entry doesn't exist"
		else bytesread = status;      // treat as "IO error" - no loss of info
		return bytesread;
	}
	
	// check for tree changes
	if(currentTreeNumber!=thetree->GetTreeNumber()){
		if(verbosity) std::cout<<"MTreeReader "<<name<<" changed tree number from "<<currentTreeNumber
		         <<" to "<<thetree->GetTreeNumber()<<" while going from entry "<<currentEntryNumber
		         <<" to "<<entry_number<<std::endl;
		// new tree
		currentTreeNumber = thetree->GetTreeNumber();
		thefile = thetree->GetCurrentFile();
		// TODO maybe implement some mechanism of notifying requestors?
		// maybe build a list of function pointers to invoke?
	}
	
	int bytesread=1;
	if(!skipTreeRead){
		// load data from tree
		// The function returns the number of bytes read from the input buffer.
		// If entry does not exist the function returns 0. If an I/O error occurs, the function returns -1.
		bytesread = thetree->GetEntry(entry_number);
		if(status<0){
			std::cerr<<"MTreeReader error loading next TTree from TChain! "
					 <<"TChain::GetEntry returned "<<status<<"\n";
		}
	}
	if(bytesread>0) currentEntryNumber = entry_number;
	return bytesread;
}

long MTreeReader::GetEntriesFast(){
	return thetree->GetEntriesFast();
}

long MTreeReader::GetEntries(){
	return thetree->GetEntries();
}

MTreeReader::~MTreeReader(){
	if(iownthisfile){
		//if(thechain) thechain->ResetBranchAddresses();  // are these mutually exclusive?
		if(thetree) thetree->ResetBranchAddresses();      // 
		if(thefile){
			thefile->Close();
			delete thefile;
		}
	}
}

void MTreeReader::SetOwnsFile(bool ownsfile){
	// Set whether the MTreeReader owns the file it's reading from.
	// If so it will close and free the file in the destructor.
	// If not, the file will be assumed to be closed and freed elsewhere.
	// (e.g. when using a TreeManager, the TreeManager destructor closes the file).
	iownthisfile=ownsfile;
}

// misc operations
void MTreeReader::SetVerbosity(int verbin){
	verbosity=verbin;
	notifier.SetVerbosity(verbin);
}

void MTreeReader::SetAutoClear(bool autoclearin){
	autoclear=autoclearin;
}

// file/tree level getters
TFile* MTreeReader::GetFile(){
	return thefile;
}

// bit of a confusing set here.
TTree* MTreeReader::GetCurrentTree(){
	// return the tree containing the current element if we're processing a TChain
	// (still valid if processing a single TTree)
	return thetree->GetTree();
}

TChain* MTreeReader::GetChain(){
	// return a TChain: we can only do this if we're processing a TChain
	// otherwise return a nullptr (we shouldn't return a child class pointer to a base class object)
	TChain* c = dynamic_cast<TChain*>(thetree); // returns nullptr if this isn't actually a TChain
	if(c==nullptr){
		std::cerr<<"Warning: MTreeReader::GetChain called when processing a TTree"<<std::endl;
	}
	return c;
}

TTree* MTreeReader::GetTree(){
	// return the TChain* cast to a TTree*, if processing a TChain, or the TTree* being processed otherwise.
	//  A TTree/TChain agnostic way to get the full extend of whatever's being process
	return thetree;
}


uint64_t MTreeReader::GetEntryNumber(){
	return currentEntryNumber;
}

// branch map getters
std::map<std::string,std::string> MTreeReader::GetBranchTypes(){
	return branch_types;
}

std::map<std::string,std::string> MTreeReader::GetBranchTitles(){
	return branch_titles;
}

std::map<std::string, intptr_t> MTreeReader::GetBranchAddresses(){
	// dangerous! Some of these are pointers to pointers!
	return branch_value_pointers;
}

// specific branch getters
TBranch* MTreeReader::GetBranch(std::string branchname){
	if(branch_pointers.count(branchname)){
		return branch_pointers.at(branchname);
	} else {
		std::cerr<<"No such branch "<<branchname<<std::endl;
		return nullptr;
	}
}

std::string MTreeReader::GetBranchType(std::string branchname){
	if(branch_types.count(branchname)){
		return branch_types.at(branchname);
	} else {
		std::cerr<<"No such branch "<<branchname<<std::endl;
		return "";
	}
	return ""; // dummy to silence warning
}

int MTreeReader::DisableBranches(std::vector<std::string> branchnames){
	int success=1;
	// disable branches by name
	for(auto&& branchname : branchnames){
		if(branch_pointers.count(branchname)){
			branch_pointers.at(branchname)->SetStatus(0);
		} else {
			std::cerr<<"No such branch "<<branchname<<std::endl;
			success=0;
		}
	}
	return success;
}

int MTreeReader::EnableBranches(std::vector<std::string> branchnames){
	int success=1;
	// disable branches by name
	for(auto&& branchname : branchnames){
		if(branch_pointers.count(branchname)){
			branch_pointers.at(branchname)->SetStatus(1);
		} else {
			std::cerr<<"No such branch "<<branchname<<std::endl;
			success=0;
		}
	}
	return success;
}

int MTreeReader::OnlyDisableBranches(std::vector<std::string> branchnames){
	// enable all branches except those named
	int num_named_branches=branchnames.size();
	for(auto&& abranch : branch_pointers){
		if(std::find(branchnames.begin(),branchnames.end(),abranch.first)!=branchnames.end()){
			abranch.second->SetStatus(0);
			--num_named_branches;
		} else {
			abranch.second->SetStatus(1);
		}
	}
	// validation branch that we recognised all the branches
	for(auto&& abranch : branchnames){
		if(branch_pointers.count(abranch)==0){
			std::cerr<<"Warning! Did not recognise branch "<<abranch
			         <<" in active branches list!"<<std::endl;
		}
	}
	// return whether we found all branches in the list given
	return (num_named_branches==0);
}

int MTreeReader::OnlyEnableBranches(std::vector<std::string> branchnames){
	// disable all branches except those named
	int num_named_branches=branchnames.size();
	for(auto&& abranch : branch_pointers){
		if(std::find(branchnames.begin(),branchnames.end(),abranch.first)!=branchnames.end()){
			abranch.second->SetStatus(1);
			--num_named_branches;
		} else {
			abranch.second->SetStatus(0);
		}
	}
	// validation branch that we recognised all the branches
	for(auto&& abranch : branchnames){
		if(branch_pointers.count(abranch)==0){
			std::cerr<<"Warning! Did not recognise branch "<<abranch
			         <<" in active branches list!"<<std::endl;
		}
	}
	// return whether we found all branches in the list given
	return (num_named_branches==0);
}

// for SKROOT files this is set in TreeReader tool... is this a good idea?
void MTreeReader::SetMCFlag(bool MCin){
	isMC = MCin;
}

bool MTreeReader::GetMCFlag(){
	return isMC;
}
