#ifndef MTREECUT_H
#define MTREECUT_H
#include <string>
#include <vector>
#include <set>

#include "TEntryList.h"
#include "TTree.h"
#include "TObjString.h"

#include "SerialisableObject.h"  // so we can put these in a BStore
#include "BinaryStream.h"        // so we can put these in a BStore

class MTreeCut : public SerialisableObject {
	
	public:
	// for all types
	MTreeCut(TFile* outfilein, std::string cutname);
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
	//  1: an "event" is described by a specific element of an array branch in TTree entry
	//  2: an "event" is described by a set of >1 elements in various array branches in a TTree entry
	// 
	// e.g., if TTree has branches:
	// "num_hits" (float)
	// "hit_charge" (float[num_hits])
	// "hit_time" (float[num_hits])
	// then a cut "num_hits>5" would be type 0
	//      a cut "hit_charge>5.f" would be type 1 (an "event" represents a single hit that passes this cut)
	//      a cut "hit_charge>5.f && hit_time<25e3" would be type 2
	
	// for type 1
	std::string additional_branchname;
	std::vector<std::string> linked_branch_list;
	std::set<size_t> indexes_this_entry;               //! use a set to prevent duplicates.
	std::set<size_t>* indexes_this_entry_p=nullptr;   // required when reading, must stay in memory!
	
	// for type 2
	std::vector<std::string> additional_branchnames;
	std::vector<std::vector<std::string>> linked_branch_lists;
	std::set<std::vector<size_t>> indices_this_entry;
	std::set<std::vector<size_t>>* indices_this_entry_p=nullptr;
	
	private:
	Long64_t current_entry=-1;
	Long64_t tlist_entry=-1;
	Long64_t total_entries=-1;
	
	public:
	void Initialize(int type_in);
	void Initialize(int type_in, std::string indexcutbranch, std::vector<std::string> linkedbranches);
	void Initialize(int type_in, std::vector<std::string> indexcutbranches, std::vector<std::vector<std::string>> linkedbranches);
	void SetEntryList(TEntryList* inelist);
	void SetTree(TTree* intree);
	void SetSaveFile(TFile* outfilein);
	void SetBranchAddresses();
	bool GetMetaInfo();
	bool Enter(Long64_t entry_number, TTree* treeptr=nullptr);
	bool Enter(Long64_t entry_number, size_t index, TTree* treeptr=nullptr);
	bool Enter(Long64_t entry_number, std::vector<size_t>& indices, TTree* treeptr=nullptr);
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
	
};

#endif
