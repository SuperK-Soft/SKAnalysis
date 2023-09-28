#ifndef MTREESELECTION_H
#define MTREESELECTION_H

#include <map>
#include <utility>
#include <string>
#include <iostream>

#include "MTreeCut.h"

#include "SerialisableObject.h"  // so we can put these in a BStore
#include "BinaryStream.h"        // so we can put these in a BStore

#include "TFile.h"
class TTree;
//class TFile;
class MTreeReader;

//template<typename T>
//std::pair<intptr_t, size_t> create_cut_pair(T branch, size_t index);

//struct pairBuilder {
//	std::pair<intptr_t, size_t> cut_pair;
//	template<typename T>
//	pairBuilder(const T& branchptr, const size_t& index) : cut_pair(create_cut_pair(value, index)) {};
//	pairBuilder(const std::string & name, std::initializer_list<pairBuilder> values) : cut_pair(create_cut_pair(name, values)) {};
//};

class MTreeSelection : public SerialisableObject {
	
	public:
	MTreeSelection();
	MTreeSelection(MTreeReader* treereader, std::string fname, std::string distrofname="");  // for writing
	MTreeSelection(std::string cutFilein);                       // for reading
	~MTreeSelection();
	bool SetTreeReader(MTreeReader* treereaderin);
	void MakeOutputFile(std::string fname, std::string distrofname="");
	bool NoteCut(std::string cutname, std::string description, double low=DOUBLE_MIN, double high=DOUBLE_MAX);
	bool AddCut(std::string cutname, std::string description, double low=DOUBLE_MIN, double high=DOUBLE_MAX);
	bool AddCut(std::string cutname, std::string description, std::string branchname, double low=DOUBLE_MIN, double high=DOUBLE_MAX); // type 1
	bool AddCut(std::string cutname, std::string description, std::vector<std::string> branchnames, double low=DOUBLE_MIN, double high=DOUBLE_MAX);  // type 1 or 2
	bool CheckCut(std::string cutname);   // just check we know this cut
	void IncrementEventCount(std::string cutname);
	// apply cut and add if it passes (new way)
	bool ApplyCut(std::string cutname, double val);
	bool ApplyCut(std::string cutname, double val, size_t index);
	bool ApplyCut(std::string cutname, double val, std::vector<size_t> indices);
	// just add directly (old way, only call for passing events)
	bool AddPassingEvent(std::string cutname);
	bool AddPassingEvent(std::string cutname, size_t index);
	bool AddPassingEvent(std::string cutname, std::vector<size_t> indices);
	
	/*
	template<typename T, class = typename std::enable_if<!std::is_same<T,TTree>::value>::type>
	bool AddPassingEvent(std::string cutname, const T* branch, size_t index);
	
	template<typename T>
	bool AddPassingEvent(std::string cutname, basic_array<T*>& branch, size_t index);
	
	bool AddPassingEvent(std::string cutname, std::vector<std::pair<intptr_t, size_t>>& indices);
	
	template <typename T, typename... Rest>
	bool AddPassingEvent(std::string cutname, T& abranch, size_t index, Rest... rest);
	
	template <typename T, typename... Rest>
	bool AddPassingEvent(std::string cutname, basic_array<T*>& abranch, size_t index, Rest... rest);
	
	template <typename T, typename... Rest>
	bool AddPassingEvent(std::string cutname, std::vector<std::pair<intptr_t, size_t>>& pairs, T abranch, size_t index, Rest... rest);
	
	template <typename T, typename... Rest>
	bool AddPassingEvent(std::string cutname, std::vector<std::pair<intptr_t, size_t>>& pairs, basic_array<T*>& abranch, size_t index, Rest... rest);
	*/
	
	std::string BranchAddressToName(intptr_t branchptr);
	void PrintCuts();
	bool Write();
	
	bool LoadCutFile(std::string cutFilein);
	Long64_t GetNextEntry(std::string cutname="");
	bool GetPassesCut(std::string cutname);
	bool GetPassesCut(std::string cutname, size_t index);
	bool GetPassesCut(std::string cutname, std::vector<size_t> indices);
	std::set<size_t> GetPassingIndexes(std::string cutname);
	std::set<std::vector<size_t>> GetPassingIndices(std::string cutname);
	Long64_t GetEntries(std::string cutname);
	bool SetEntries(Long64_t nentries);
	MTreeReader* GetTreeReader();
	std::string GetTopCut();
	
	private:
	std::vector<std::string> FindLinkedBranches(std::string cut_branch);
	
	// track num events passing cuts.
	std::vector<std::string> cut_order;
	std::map<std::string, uint64_t> cut_tracker;
	std::map<std::string, MTreeCut*> cut_pass_entries;
	
	MTreeReader* treereader=nullptr;
	std::map<intptr_t, std::string> branch_addresses;
	std::map<std::string, intptr_t> known_branches;
	
	//BoostStore* outstore=nullptr;
	TFile* outfile=nullptr;
	bool initialized=false; // written initial meta-data  - FIXME redundant/broken? (order of cuts?)
	
	// if making distributions of the variables we're cutting on
	TFile* distros_file=nullptr;
	TTree* distros_tree=nullptr;
	
	// involved in reading
	TFile* cutfile=nullptr;
	Long64_t current_entry=-1;
	std::map<std::string, bool> did_pass_cut;
	
	// required for SerialisableObjects
	protected:
	bool Serialise(BinaryStream &bs);
	bool Print();
	std::string GetVersion();
};


/*
// ↓ the additional template argument 'class' is required with the enable_if to prevent ambiguity with
// ↓ AddPassingEvent(std::string cutname, TTree* tree, Long64_t entry_number) - i.e. the below with type T=TTree
template<typename T, class>
bool MTreeSelection::AddPassingEvent(std::string cutname, const T* branch, size_t index){
	
	// check we know this cut
	if(cut_pass_entries.count(cutname)==0){
		std::cerr<<"MTreeSelection::AddPassingEvent called with unknown cut "<<cutname
				 <<"\nPlease call MTreeSelection::AddCut on all cuts in order first"<<std::endl;
		return false;
	}
	
	// a more convenient interface for the user: they give us the address
	// of the array variable, and we use it to look up the branch name.
	// it's easier because variables are at the site where this gets called,
	// whereas branch names in SK are often arbitrary and confusing
	
	// initialize the cut if not already
	if(cut_pass_entries.at(cutname)->type<0){
		std::string branchname = BranchAddressToName(reinterpret_cast<intptr_t>(branch));
		cut_pass_entries.at(cutname)->Initialize(1, branchname, FindLinkedBranches(branchname));
	}
	
	// save the TTree entry and index array.
	bool newentry = cut_pass_entries.at(cutname)->Enter(treereader->GetEntryNumber(), index, treereader->GetTree());
	
	// increment number of events passing the cut
	if(not newentry) return false;  // prevent double counting
	IncrementEventCount(cutname);
	
	return true;
}

template<typename T>
bool MTreeSelection::AddPassingEvent(std::string cutname, basic_array<T*>& branch, size_t index){
	return AddPassingEvent(cutname, branch.data(), index);
}

bool MTreeSelection::AddPassingEvent(std::string cutname, std::vector<std::pair<intptr_t, size_t>>& indices){
	
	// check we know this cut
	if(cut_pass_entries.count(cutname)==0){
		std::cerr<<"MTreeSelection::AddPassingEvent called with unknown cut "<<cutname
				 <<"\nPlease call MTreeSelection::AddCut on all cuts in order first"<<std::endl;
		return false;
	}
	
	// initialize the cut if not already
	if(cut_pass_entries.at(cutname)->type<0){
		std::vector<std::string> cutbranches;
		std::vector<std::vector<std::string>> linked_branch_lists;
		for(auto branch : indices){
			std::string branchname = BranchAddressToName(branch.first);
			cutbranches.emplace_back(branchname);
			if(branchname=="") std::cerr<<"empty branchname for "<<cutname<<std::endl;
			linked_branch_lists.emplace_back(FindLinkedBranches(branchname));
		}
		cut_pass_entries.at(cutname)->Initialize(2, cutbranches, linked_branch_lists);
	}
	
	// need to strip out just the indices... :/
	std::vector<size_t> justindices(indices.size());
	for(int i=0; i<indices.size(); ++i){
		justindices[i]=indices[i].second;
	}
	
	// save the TTree entry number
	bool newentry = cut_pass_entries.at(cutname)->Enter(treereader->GetEntryNumber(), justindices, treereader->GetTree());
	
	// increment number of events passing the cut
	if(not newentry) return false;  // prevent double counting
	IncrementEventCount(cutname);
	
	return true;
}

template <typename T, typename... Rest>
bool MTreeSelection::AddPassingEvent(std::string cutname, T& abranch, size_t index, Rest... rest) {
	std::vector<std::pair<intptr_t, size_t>> pairs{{reinterpret_cast<intptr_t>(abranch), index}};
	return AddPassingEvent(cutname, pairs, rest...);
}

template <typename T, typename... Rest>
bool MTreeSelection::AddPassingEvent(std::string cutname, basic_array<T*>& abranch, size_t index, Rest... rest) {
	std::vector<std::pair<intptr_t, size_t>> pairs{{reinterpret_cast<intptr_t>(abranch.data()), index}};
	return AddPassingEvent(cutname, pairs, rest...);
}

template <typename T, typename... Rest>
bool MTreeSelection::AddPassingEvent(std::string cutname, std::vector<std::pair<intptr_t, size_t>>& pairs, T abranch, size_t index, Rest... rest){
	pairs.emplace_back(reinterpret_cast<intptr_t>(abranch), index);
	return AddPassingEvent(cutname, pairs, rest...);
}

template <typename T, typename... Rest>
bool MTreeSelection::AddPassingEvent(std::string cutname, std::vector<std::pair<intptr_t, size_t>>& pairs, basic_array<T*>& abranch, size_t index, Rest... rest){
	pairs.emplace_back(reinterpret_cast<intptr_t>(abranch.data()), index);
	return AddPassingEvent(cutname, pairs, rest...);
}
*/

#endif
