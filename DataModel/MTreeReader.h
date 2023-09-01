/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef MTreeReader_H
#define MTreeReader_H

#include <string>
#include <map>
#include <utility> // pair

#include "basic_array.h"

#include "TObject.h"

class TFile;
class TChain;
class TTree;
class TBranch;
class TLeaf;
class MTreeReader;

class Notifier : public TObject {
	public:
	bool Notify();
	MTreeReader* treeReader=nullptr;
	int verbosity;
	void SetReader(MTreeReader* in){ treeReader=in; }
	void SetVerbosity(int verbin){ verbosity=verbin; }
};

class MTreeReader {
	friend class Notifier;
	public:
	
	MTreeReader(std::string iname, std::string filename, std::string treename);
	MTreeReader(std::string iname="myReader");
	~MTreeReader();
	int Load(std::string filename, std::string treename);
	int LoadFile(std::string filename);
	int LoadTree(std::string treename);
	int Load(TTree* thetreein);
	int Load(std::vector<std::string> filelist, std::string treename);
	void SetClosed();
	
	// get a pointer to an object
	template<typename T>
	int GetBranchValue(std::string branchname, const T* &pointer_in){
		if(branch_value_pointers.count(branchname)==0){
			std::cerr<<"No such branch '"<<branchname<<"'"<<std::endl;
			std::cerr<<"known branches: "<<branchnamestring<<std::endl;
			return 0;
		}
		pointer_in = reinterpret_cast<const T*>(branch_value_pointers.at(branchname));
		if(verbosity>3) std::cout<<"retrieved pointer to "<<type_name<T>()<<" at "<<pointer_in<<std::endl;
		return 1;
	}
	
	// get a pointer to an object - non-const version required by some existing algorithms
	template<typename T>
	int GetBranchValue(std::string branchname, T* &pointer_in){
		const T* tmp = const_cast<const T*>(pointer_in);
		int ret = GetBranchValue(branchname, tmp);
		if(ret) pointer_in = const_cast<T*>(tmp);
		return ret;
	}
	
	// get the object itself - valid for primitives only (could add objects if we made copies)
	template<typename T>
	int GetBranchValue(std::string branchname, T& ref_in){
		// check we know this branch
		if(branch_value_pointers.count(branchname)==0){
			std::cerr<<"No such branch '"<<branchname<<"'"<<std::endl;
			std::cerr<<"known branches: "<<branchnamestring<<std::endl;
			return 0;
		}
		// check if the branch is a primitive
		if(branch_isobject.at(branchname)||branch_isarray.at(branchname)){
			std::cerr<<"Branch "<<branchname
				 <<" is not a primitive; please pass a suitable const pointer"
				 <<" or basic_array to GetBranchValue()"<<std::endl;
			// TODO copy-construct an object, if they really want
			// requires a suitable copy constructor (or operator=) exists for the class
			return 0;
		}
		// else for primitives, de-reference the pointer to allow the user a copy
		T* objp = reinterpret_cast<T*>(branch_value_pointers.at(branchname));
		ref_in = *objp;
		return 1;
	}
	
	// specialization for array references
	
	// 1D variant
	template<typename T, std::size_t N>
	int GetBranchValue(std::string branchname, T (&ref_in)[N]){
		return GetArrayBranchValue(branchname, &ref_in[0], N);
	}
	
	// 2D variant
	template<typename T, std::size_t N, std::size_t M>
	int GetBranchValue(std::string branchname, T (&ref_in)[N][M]){
		return GetArrayBranchValue(branchname, &ref_in[0][0], N, M);
	}
	
	// 3D variant... that's as far as i'll go.
	template<typename T, std::size_t N, std::size_t M, std::size_t S>
	int GetBranchValue(std::string branchname, T (&ref_in)[N][M][S]){
		return GetArrayBranchValue(branchname, &ref_in[0][0][0], N, M, S);
	}
	
	// the above defer to this - copy the data to ther user's array
	// (the user really ought not to pass us an array reference, as it requires
	//  an unnecessary copy relative to just giving us a pointer we can direct,
	//  but not accepting array references looks like a bug to the user)
	template<typename T>
	int GetArrayBranchValue(std::string branchname, T* arr_in, std::size_t NCOL, std::size_t NROW=1, std::size_t NAISLE=1){
		// check we know this branch
		if(branch_value_pointers.count(branchname)==0){
			std::cerr<<"No such branch '"<<branchname<<"'"<<std::endl;
			std::cerr<<"known branches: "<<branchnamestring<<std::endl;
			return 0;
		}
		// check if the branch is an array - this template specialization is only for arrays
		if(not branch_isarray.at(branchname)){
			std::cerr<<"Branch "<<branchname
				 <<" is not an array; please check your datatype to GetBranchValue()"<<std::endl;
			return 0;
		}
		// check the passed array has suitable dimensions.
		// first we need to know the actual array dimensions
		std::vector<size_t> branchdims = GetBranchDims(branchname);
		// for dynamic arrays we may need to update our pointer to the stored array
		// not sure if we should bail if the user is trying to put a dynamic array
		// into a static-sized array variable.... continue for now.
		UpdateBranchPointer(branchname);
		
		// the user's array must be at least as large as required
		// first check the number of dimensions is sufficient
		int ndims = 1;
		if(NROW>1) ++ndims;
		if(NAISLE>1) ++ndims;
		if(ndims!=branchdims.size()){
			std::cerr<<"passed an array reference of dimensionality "<<ndims
					 <<" for branch "<<branchname<<" which has dimensionality "
					 <<branchdims.size()<<std::endl;
			return 0;
		}
		
		// next check each dimension has sufficient capacity
		bool cap_ok = true;
		if(NCOL<branchdims.at(0)) cap_ok = false;
		if(ndims>1 && NROW<branchdims.at(1)) cap_ok = false;
		if(ndims>2 && NAISLE<branchdims.at(2)) cap_ok = false;
		if(not cap_ok){
			std::cerr<<"passed an array with dimensions ["<<NCOL<<"]";
			if(NROW>1) std::cerr<<"["<<NROW<<"]";
			if(NAISLE>1) std::cerr<<"["<<NAISLE<<"] ";
			std::cerr<<" which is insufficient for the data dimensions ";
			for(int i=0; i<ndims; ++i) std::cerr<<"["<<branchdims.at(i)<<"]";
			std::cerr<<std::endl;
			return 0;
		}
		
		// copy the data to the user's array
		int data_cols = branchdims.at(0);
		int data_rows = (ndims>1) ? branchdims.at(1) : 1;
		int data_aisles = (ndims>2) ? branchdims.at(2) : 1;
		T* objp = reinterpret_cast<T*>(branch_value_pointers.at(branchname));
		for(int aisle=0; aisle<NAISLE; ++aisle){
			for(int row=0; row<NROW; ++row){
				for(int col=0; col<NCOL; ++col){
					std::size_t dest_flat_index = aisle*NCOL*NROW + row*NCOL + col;
					std::size_t source_flat_index = aisle*data_rows*data_cols + row*data_cols + col;
					if((col+1)>branchdims.at(0) ||
					   (ndims>1 && (row+1)>branchdims.at(1)) ||
					   (ndims>2 && (aisle+1)>branchdims.at(2))){
						arr_in[dest_flat_index] = 0;
					} else {
						arr_in[dest_flat_index] = objp[source_flat_index];
					}
				}
			}
		}
		return 1;
	}
	
	// specialization for arrays using basic_array
	template<typename T>
	int GetBranchValue(std::string branchname, basic_array<T>& ref_in){
		// check we know this branch
		if(branch_value_pointers.count(branchname)==0){
			std::cerr<<"No such branch '"<<branchname<<"'"<<std::endl;
			std::cerr<<"known branches: "<<branchnamestring<<std::endl;
			return 0;
		}
		// check if the branch is an array - this template specialization is only for arrays
		if(not branch_isarray.at(branchname)){
			std::cerr<<"Branch "<<branchname
				 <<" is not an array; please check your datatype to GetBranchValue()"<<std::endl;
			return 0;
		}
		// for dynamic arrays we may need to update our pointer to the stored array
		UpdateBranchPointer(branchname);
		// next we need to know the array dimensions, which may vary by entry
		std::vector<size_t> branchdims = GetBranchDims(branchname);
		// finally construct and return the wrapper
		ref_in = basic_array<T>(branch_value_pointers.at(branchname),branchdims);
		return 1;
	}
	
	// aliases of GetBranchValue
	template<typename T>
	int Get(std::string branchname, const T* &pointer_in){
		return GetBranchValue(branchname, pointer_in);
	}
	
	template<typename T>
	int Get(std::string branchname, T& ref_in){
		return GetBranchValue(branchname, ref_in);
	}
	
	template<typename T>
	int Get(std::string branchname, basic_array<T>& ref_in){
		return GetBranchValue(branchname, ref_in);
	}
	
	// misc operations
	void SetVerbosity(int verbin);
	
	// file/tree level getters
	std::string GetName();
	void SetName(std::string iname);
	TFile* GetFile();
	TTree* GetTree();
	TTree* GetCurrentTree();
	TChain* GetChain();
	uint64_t GetEntryNumber();
	
	// tree operations
	int Clear();
	void SetAutoClear(bool autoclearin);
	int GetEntry(long entry_number, bool skipTreeRead=false);
	long GetEntriesFast();
	long GetEntries();
	int DisableBranches(std::vector<std::string> branchnames);
	int EnableBranches(std::vector<std::string> branchnames);
	int OnlyEnableBranches(std::vector<std::string> branchnames);
	int OnlyDisableBranches(std::vector<std::string> branchnames);
	
	// maps of branch properties
	std::map<std::string,std::string> GetBranchTypes();
	std::map<std::string,intptr_t> GetBranchAddresses();
	std::map<std::string, std::string> GetBranchTitles();
	
	// specific branch properties
	TBranch* GetBranch(std::string branchname);
	std::string GetBranchType(std::string branchname);
	std::vector<size_t> GetBranchDims(std::string branchname);
	
	// random assistive functions
	void SetMCFlag(bool MCin);
	bool GetMCFlag();
	
	// functions
	int ParseBranches();
	int ParseBranchDims(std::string branchname);
	int UpdateBranchPointer(std::string branchname);
	int UpdateBranchPointers(bool all=false);
	
	private:
	
	// variables
	std::map<std::string,TBranch*> branch_pointers;  // branch name to TBranch*
	std::map<std::string,bool> branch_istobject;     // branch inherits from TObject so has Clear method
	std::map<std::string,TLeaf*> leaf_pointers;      // branch name to TLeaf*
	std::map<std::string,std::string> branch_titles; // branch name to title (including type string)
	std::map<std::string,intptr_t> branch_value_pointers; // branch name to pointer to value, cast to intptr_t
	std::map<std::string,std::string> branch_types;  // branch name to string describing type - not good for arrays
	std::map<std::string,bool> branch_isobject;      // does branch hold an object
	std::map<std::string,bool> branch_isarray;       // does branch hold a (c-style) array
	std::map<std::string,std::vector<std::pair<std::string,int>>> branch_dimensions; // dims of variable size arrays
	std::map<std::string,std::vector<size_t>> branch_dims_cache; // dims of constant sized arrays
	
	TFile* thefile=nullptr;
	TTree* thetree=nullptr;          // generic, if working with a tchain we cast it to a TTree
	bool autoclear=false;            // call 'Clear' method on all object branches before GetEntry
	int verbosity=1;                 // TODO add to constructor
	uint64_t currentEntryNumber=0;
	int currentTreeNumber=0;
	bool isMC=false;

	Notifier notifier;

	std::string name="";
	std::string branchnamestring="{}";
	
};

/*
// mechanism to check if a given type has a Clear() method.
// ROOT seems to flag both objects and stl containers as inheriting from TObject
// with InheritsFrom("TLeafElement")
// from https://stackoverflow.com/a/29772824
template<typename T>
struct has_clear {
	// NOTE: sig_matches() must come before fn_exists() as it is used for its type.
	// Also, no function bodies are needed as they are never called.
	
	// This matching sig results in a return type of true_type
	template<typename Q>
	static auto sig_matches(void(Q::*)()) -> std::true_type;
	
	// If the member function Q::Clear exists and a sig_matches() function
	// exists with the required sig, then the return type is the return type of
	// sig_matches(), otherwise this function can't exist because at least one
	// the types don't exist so match against fn_exists(...).
	template <typename Q>
	static auto fn_exists(std::nullptr_t) -> decltype(sig_matches<Q>(&Q::Clear));
	
	// Member function either doesn't exist or doesn't match against a 
	// sig_matches() function.
	template<typename Q>
	static auto fn_exists(...) -> std::false_type;
	
	// Intermediate storage of type for clarity
	typedef decltype(fn_exists<T>(nullptr)) type;
	
	// Storing the resulting value
	static int const value = type::value;
};
*/


#endif // defined MTreeReader_H
