#ifndef MTREECUT_H
#define MTREECUT_H
#include <string>
#include <vector>
#include <set>
#include <limits>

#include "TEntryList.h"
#include "TTree.h"
#include "TObjString.h"

#include "SerialisableObject.h"  // so we can put these in a BStore
#include "BinaryStream.h"        // so we can put these in a BStore

namespace {
	constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();
	constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();
}

class MTreeReader;

class MTreeCut : public SerialisableObject {
	
	public:
	// for all types
	MTreeCut(TFile* outfilein, std::string cutname, std::string description, double low=DOUBLE_MIN, double high=DOUBLE_MAX);
	MTreeCut(std::string cutname, TEntryList* inelist, TTree* intree);
	~MTreeCut();
	std::string mode="";  // can be "read" or "write". Determines whether destructor performs cleanup.
	std::string cut_name;
	std::string cut_description;   // TODO store the cut values etc in here! if mu_time<250, store the value 250!
	TFile* outfile=nullptr;
	TEntryList* ttree_entries=nullptr;
	TTree* additional_indices=nullptr;
	int type=-1;
	// -1; uninitialized,
	//  0: simple, an "event" is a TTree entry
	//  1: an "event" is a single element of an array branch in TTree entry
	//  2: an "event" is a combination of >1 elements in multiple *independent* array branches in a TTree entry
	// 
	// e.g., if a TTree has branches:
	// "num_prompt_events" (int)
	// "prompt_t" (float[num_prompt_events])
	// "prompt_e" (float[num_prompt_events])
	// "num_delayed_events" (int)
	// "delayed_t" (float[num_delayed_events])
	// "delayed_e" (float[num_delayed_events])
	// 
	// then a cut "num_delayed_events>0" would be type 0
	//      a cut "prompt_e > 8.0" would be type 1 (an "event" represents a single entry in the 'prompt' branch arrays that passes this cut)
	//      a cut "prompt_e > 8.0 && delayed_e > 5.0" would be a type 2 (an "event" represents the COMBINATION of a single entry in the 'prompt' branch
	//                                                                    arrays TOGETHER WITH a single entry in the INDEPENDENT 'delayed' branch arrays)
	// NOTE:  cut "prompt_e > 8.0 && prompt_t < 25e3" would NOT be type 2, because these arrays are NOT independent:
	// i.e. the passing index in the prompt_e array is necessarily the same as the passing index in the prompt_t array,
	// as these arrays represent different properties of the same underlying object (they sharing indexing).
	// so this is just the '&&' of two type 1 cuts, which still represents a type 1 cut.
	
	// for type 1
	std::string additional_branchname;
	std::vector<std::string> linked_branch_list;
	std::set<size_t> indexes_this_entry;               //! use a set to prevent duplicates.
	std::set<size_t>* indexes_this_entry_p=nullptr;    // required when reading, must stay in memory!
	
	// for type 2
	std::vector<std::string> additional_branchnames;
	std::vector<std::vector<std::string>> linked_branch_lists;
	std::set<std::vector<size_t>> indices_this_entry;
	std::set<std::vector<size_t>>* indices_this_entry_p=nullptr;
	
	private:
	MTreeReader* theReader=nullptr;
	Long64_t current_entry=-1;
	Long64_t tlist_entry=-1;
	Long64_t total_entries=-1;
	
	public:
	void Initialize(int type_in, MTreeReader* readerIn, TTree* distro_tree=nullptr);
	void Initialize(int type_in, std::string indexcutbranch, std::vector<std::string> linkedbranches, MTreeReader* readerIn, TTree* distro_tree=nullptr);
	void Initialize(int type_in, std::vector<std::string> indexcutbranches, std::vector<std::vector<std::string>> linkedbranches, MTreeReader* readerIn, TTree* distro_tree=nullptr);
	bool MakeDistroBranches(TTree* distros_tree);
	void SetBranchAddresses();
	bool GetMetaInfo();
	void SetDescription(std::string desc);
	
	bool Enter();
	bool Enter(size_t index);
	bool Enter(std::vector<size_t>& indices);
	
	// check if the value passes the required check and if so calls Enter
	bool Apply(double value);
	bool Apply(double value, size_t index);
	bool Apply(double value, std::vector<size_t> indices);
	
	bool Flush();
	void Write();
	Long64_t GetCurrentEntry();
	Long64_t GetNextEntry();
	std::set<size_t> GetPassingIndexes();
	std::set<std::vector<size_t>> GetPassingIndices();
	
	protected:
	// required for SerialisableObjects
	bool Serialise(BinaryStream &bs);
	bool Print();
	std::string GetVersion();
	bool got_low_thresh=false;
	double low_thresh=DOUBLE_MIN;
	bool got_high_thresh=false;
	double high_thresh=DOUBLE_MAX;
	bool got_bool_req=false;
	
	// if saving distributions
	bool save_distros=false;
	double branch_val;
	bool pass_val;
	TBranch* values_branch=nullptr;
	TBranch* pass_branch=nullptr;
	
};

#endif
