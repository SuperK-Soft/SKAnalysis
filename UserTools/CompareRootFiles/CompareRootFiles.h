/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef CompareRootFiles_H
#define CompareRootFiles_H

#include <string>
#include <iostream>
#include <set>

#include "Tool.h"
#include "MTreeReader.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

/**
* \class CompareRootFiles
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
struct data_instance;

struct data_instance {
	std::string name;
	std::string type_as_string;
	void* address=nullptr;
	long offset=0;
	TBranch* branch_ptr=nullptr;
	
	int instance_type;
	// -1: unknown
	//  0: primitive
	//  1: statically sized array
	//  2: dynamically sized array
	//  3: stl container
	//  4: class object
	
	// for primitives we need the address and its datatype, or at least its size in bytes
	long item_size;
	
	// for arrays we also need the dimensions. these may be static...
	std::vector<int> static_dims;
	// ... and/or one dimension may also be of variable size specified by another branch
	std::string dimension_branch="";
	int* dimension_ptr=nullptr;
	
	// for stl containers, to compare two instances we need to compare the sizes,
	// and then iterate over the elements and compare those. For this we need to know
	// the container type and address (for getting the size and addresses of elements),
	// and the type of the contained element type (for comparing them)
	data_instance* contained_type = nullptr;
	
	// finally for class objects, to compare two instances we need to compare all members.
	// to do this we need a map of the members
	std::map<std::string, data_instance> members;
};

struct branch_structure {
	TBranch* branch_ptr=nullptr;
	std::string type_as_string;
	data_instance held_data;
};

struct shared_tree {
	TTree* file1_tree;
	TTree* file2_tree;
	std::map<std::string, branch_structure> file1_branches;
	std::map<std::string, branch_structure> file2_branches;
};

class CompareRootFiles: public Tool {
	
	public:
	CompareRootFiles();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	bool LoadConfig(std::string configfile);
	std::string GetBranchType(TBranch* branch);
	void ParseSharedTrees();
	void ParseTree(std::pair<const std::string, shared_tree> &tree);
	bool ParseBranch(branch_structure &branch);
	bool ObjectToDataInstance(data_instance& thedata);
	bool CompareBranchMembers(data_instance &branch1, data_instance &branch2, int* less=nullptr);
	int GetAllClassMembers(TClass* cl, std::vector<std::pair<TDataMember*, Long_t>> &members);
	int GetAllClassMembersOld(TClass* cl, Long_t offset, std::vector<std::pair<TDataMember*, Long_t>> &members);
	void TransferAddressesToOffsets();
	void TransferBranchAddresses(data_instance &thedata);
	void DoTransferAddresses(data_instance &thedata);
	void UpdateAddresses(data_instance& thedata);
	std::string GetVarName(std::string type_as_string);
	TMethodCall* GetSizeCaller(data_instance& thedata);
	TMethodCall* GetAtCaller(data_instance& thedata);
	bool GetListOfHeaders(data_instance &thedata, std::vector<std::string> &headerlist);
	const char* CppName(const char* type);
	const char* CppName(std::string type);
	bool GetNextMatchingEntries(std::pair<const std::string, shared_tree>& apair);
	
	// methods for building dictionaries on the fly
	bool LoadDictionary(data_instance& thedata);
	bool GetListOfHeaders(data_instance &thedata, std::vector<std::pair<std::string,std::string>> &headerlist);
	std::string BuildLinkDef(data_instance& thedata, std::vector<std::pair<std::string,std::string>> &headerlist);
	std::string BuildDictionary(std::string classname, std::vector<std::pair<std::string,std::string>> &headerlist);
	std::vector<std::string> GetListOfImplementationFiles(std::vector<std::pair<std::string,std::string>> &headerlist);
	std::string CompileDictionary(std::string type_as_string, std::string dictionaryfile, std::vector<std::string> implementationlist);
	
	// tool variables
	// ==============
	std::string toolName;
	std::string filename_1;
	std::string filename_2;
	TFile* file1=nullptr;
	TFile* file2=nullptr;
	std::map<std::string, shared_tree> common_trees;
	std::map<std::string, bool> active_trees; // those we're still reading entries from
	int entry_number=0;
	int max_entries=-1;
	std::map<std::string,std::string> varnames;
	std::map<std::string, TMethodCall*> sizecallers;
	std::map<std::string, TMethodCall*> atcallers;
	std::vector<std::string> source_paths;
	std::vector<std::string> include_paths;
	std::set<std::string> loaded_libraries;
	float ftolerance=0.005f;
	std::string index_name;
	data_instance* file1_index = nullptr;
	data_instance* file2_index = nullptr;
	std::map<std::string, long> entry_numbers_1; // where we are in each TTree in file 1
	std::map<std::string, long> entry_numbers_2; // where we are in each TTree in file 2
	long entry_number_1;
	long entry_number_2;
	int matching_entries = 0;
	int mismatching_entries = 0;
	std::string entry_string;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	std::stringstream smessage;
	int get_ok=0;
	
	
};


#endif
