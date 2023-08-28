/* vim:set noexpandtab tabstop=4 wrap */
#include "CompareRootFiles.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <numeric>
#include <functional>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TKey.h"
#include "TClass.h"
#include "TBaseClass.h"
#include "TDataMember.h"
#include "TRealData.h"
#include "TMethodCall.h"
#include "TClassEdit.h"
#include "TClonesArray.h"

CompareRootFiles::CompareRootFiles():Tool(){}

bool CompareRootFiles::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	LoadConfig(configfile);
	
	// Open root files
	file1 = TFile::Open(filename_1.c_str(),"READ");
	file2 = TFile::Open(filename_2.c_str(),"READ");
	if(file1==nullptr || file2==nullptr){
		if(file1==nullptr) Log(m_unique_name+" Error opening file "+filename_1,v_error,m_verbose);
		if(file2==nullptr) Log(m_unique_name+" Error opening file "+filename_2,v_error,m_verbose);
		exit(-1);
		return false;
	}
	
	// map of name to type for things in each file
	std::map<std::string, std::string> file1_contents;
	std::map<std::string, std::string> file2_contents;
	// We're primarily interested in comparing all trees we find in both files.
	std::map<std::string, TTree*> file1_trees;
	std::map<std::string, TTree*> file2_trees;
	
	// get a map of the contents of the file
	TKey *key=nullptr;
	TIter getnextkey(file1->GetListOfKeys());
	Log(m_unique_name+" looping over keys in file 1", v_debug, m_verbose);
	while ((key = (TKey*)getnextkey())){
		if(key==nullptr){
			Log(m_unique_name+" Error! null key while getting TKeys in file "+filename_2+"!",v_error,m_verbose);
			continue;
		}
		Log(m_unique_name+" processing key "+key->GetName(),v_debug,m_verbose);
		// get the type of the object
		file1_contents.emplace(key->GetName(), key->GetClassName());
		// check if it's a tree: if so we'll compare branches
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if(!cl->InheritsFrom("TTree")) continue;
		Log(m_unique_name+": it's a TTree",v_debug,m_verbose);
		file1_trees.emplace(key->GetName(),(TTree*)key->ReadObj());
	} // end loop over keys in file
	
	// do the same for our second file
	// get a map of the contents of the file
	key=nullptr;
	getnextkey = TIter(file2->GetListOfKeys());
	Log(m_unique_name+" looping over keys in file 2", v_debug, m_verbose);
	while ((key = (TKey*)getnextkey())){
		if(key==nullptr){
			Log(m_unique_name+" Error! null key while getting TKeys in file "+filename_2+"!",v_error,m_verbose);
			continue;
		}
		Log(m_unique_name+" processing key "+key->GetName(),v_debug,m_verbose);
		// get the type of the object
		file2_contents.emplace(key->GetName(), key->GetClassName());
		// check if it's a tree: if so we'll compare branches
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if(!cl->InheritsFrom("TTree")) continue;
		Log(m_unique_name+": it's a TTree",v_debug,m_verbose);
		file2_trees.emplace(key->GetName(),(TTree*)key->ReadObj());
	} // end loop over keys in file
	
	// compare the file contents:
	std::map<std::string, std::string> combined;
	combined.insert(file1_contents.begin(), file1_contents.end());
	combined.insert(file2_contents.begin(), file2_contents.end());
	bool all_keys_same=true;
	for(auto& akey : combined){
		if(file1_contents.count(akey.first) && file2_contents.count(akey.first)){
			if(file1_contents.at(akey.first)==file2_contents.at(akey.first)){
				Log(m_unique_name+": both files contain a "+akey.second+" called "+akey.first,v_debug,m_verbose);
			} else {
				Log(m_unique_name+": Type mismatch! file "+filename_1+" contains a "
					+file1_contents.at(akey.first)+" named "+akey.first+", while "
					+filename_2+" contains a "+file2_contents.at(akey.first)
					+" of the same name!",v_error,m_verbose);
				all_keys_same=false;
			}
		} else if(file1_contents.count(akey.first)){
			Log(m_unique_name+": File content mismatch! file "+filename_1+" contains a "
				+file1_contents.at(akey.first)+" named "+akey.first
				+" which is not in "+filename_2,v_error,m_verbose);
			all_keys_same=false;
		} else {
			Log(m_unique_name+": File content mismatch! file "+filename_2+" contains a "
				+file2_contents.at(akey.first)+" named "+akey.first
				+" which is not in "+filename_1,v_error,m_verbose);
			all_keys_same=false;
		}
	}
	if(all_keys_same){
		Log(m_unique_name+": Files "+filename_1+" and "+filename_2+" contain the same set of keys",v_warning,m_verbose);
	}
	
	// find any trees shared by both and we'll compare those
	std::map<std::string, TTree*> combined_trees;
	combined_trees.insert(file1_trees.begin(), file1_trees.end());
	combined_trees.insert(file2_trees.begin(), file2_trees.end());
	for(auto& akey : combined_trees){
		if(file1_trees.count(akey.first) && file2_trees.count(akey.first)){
			Log(m_unique_name+": found common tree "+akey.first,v_debug,m_verbose);
			if(file1_trees.at(akey.first)->GetEntriesFast() != file2_trees.at(akey.first)->GetEntriesFast()){
				Log(m_unique_name+": File content mismatch! file "+filename_1+" tree "+akey.first
					+" has "+toString(file1_trees.at(akey.first)->GetEntriesFast())
					+" entries, while this tree in "+filename_2+" has "
					+toString(file2_trees.at(akey.first)->GetEntriesFast())+" entries!",v_error,m_verbose);
			} else {
				Log(m_unique_name+": Both trees have the same number of entries",v_debug,m_verbose);
			}
			
			Log(m_unique_name+": checking for branches in tree "+akey.first,v_debug,m_verbose);
			// both files have a TTree by the same name, but do they share some or all branches?
			// make a list of branches in this tree in file 1....
			std::map<std::string,std::string> file1_branches;
			for(int branchi=0; branchi<file1_trees.at(akey.first)->GetListOfBranches()->GetEntries(); ++branchi){
				TBranch* thebranch = (TBranch*)file1_trees.at(akey.first)->GetListOfBranches()->At(branchi);
				file1_branches.emplace(thebranch->GetName(), GetBranchType(thebranch));
				Log(m_unique_name+" file 1 has branch "+thebranch->GetName(),v_debug,m_verbose);
			}
			// make a list of branches in this tree in file 2....
			std::map<std::string,std::string> file2_branches;
			for(int branchi=0; branchi<file2_trees.at(akey.first)->GetListOfBranches()->GetEntries(); ++branchi){
				TBranch* thebranch = (TBranch*)file2_trees.at(akey.first)->GetListOfBranches()->At(branchi);
				file2_branches.emplace(thebranch->GetName(), GetBranchType(thebranch));
				Log(m_unique_name+" file 2 has branch "+thebranch->GetName(),v_debug,m_verbose);
			}
			// merge both lists
			std::map<std::string, std::string> combined_branches;
			combined_branches.insert(file1_branches.begin(),file1_branches.end());
			combined_branches.insert(file2_branches.begin(),file2_branches.end());
			bool all_branches_same = (combined_branches.size()==file1_branches.size() &&
			                          combined_branches.size()==file2_branches.size());
			
			// scan over the combined list for those that are in both
			for(auto& bkey : combined_branches){
				if(file1_branches.count(bkey.first) && file2_branches.count(bkey.first)){
					// ensure held type is same - TODO: for classes we should also check they are the same versio
					if(file1_branches.at(bkey.first)==file2_branches.at(bkey.first)){
						// common branch!
						Log(m_unique_name+" found common branch "+bkey.first
							+" of type "+bkey.second,v_debug,m_verbose);
						if(common_trees.count(akey.first)==0){
							// first common branch of this tree - register this tree as a shared tree
							shared_tree new_tree;
							new_tree.file1_tree = file1_trees.at(akey.first);
							new_tree.file2_tree = file2_trees.at(akey.first);
							
							branch_structure file1_branch;
							file1_branch.branch_ptr = file1_trees.at(akey.first)->GetBranch(bkey.first.c_str());
							file1_branch.type_as_string = bkey.second;
							new_tree.file1_branches.emplace(bkey.first,file1_branch);
							
							branch_structure file2_branch;
							file2_branch.branch_ptr = file2_trees.at(akey.first)->GetBranch(bkey.first.c_str());
							file2_branch.type_as_string = bkey.second;
							new_tree.file2_branches.emplace(bkey.first,file2_branch);
							
							common_trees.emplace(akey.first,new_tree);
							entry_numbers_1.emplace(akey.first, -1);
							entry_numbers_2.emplace(akey.first, -1);
						} else {
							// else we know about this tree, add the branch as a common branch
							shared_tree& new_tree = common_trees.at(akey.first);
							branch_structure file1_branch;
							file1_branch.branch_ptr = file1_trees.at(akey.first)->GetBranch(bkey.first.c_str());
							file1_branch.type_as_string = bkey.second;
							new_tree.file1_branches.emplace(bkey.first,file1_branch);
							
							branch_structure file2_branch;
							file2_branch.branch_ptr = file2_trees.at(akey.first)->GetBranch(bkey.first.c_str());
							file2_branch.type_as_string = bkey.second;
							new_tree.file2_branches.emplace(bkey.first,file2_branch);
						}
					} else if(file1_branches.count(bkey.first)){
						Log(m_unique_name+": File content mismatch! file "+filename_1
							+" tree "+akey.first+" branch "+bkey.first+" is of type "
							+file1_branches.at(bkey.first)
							+" but this branch in file "+filename_2+" is of type "
							+file2_branches.at(bkey.first)
							+"!",v_error,m_verbose);
						all_branches_same=false;
					}
				} else if(file1_branches.count(bkey.first)){
					Log(m_unique_name+": File content mismatch! file "+filename_1
						+" tree "+akey.first+" contains branch "+bkey.first
						+" which is not in file "+filename_2,v_error,m_verbose);
					all_branches_same=false;
				} else {
					Log(m_unique_name+": File content mismatch! file "+filename_2
						+" tree "+akey.first+" contains branch "+bkey.first
						+" which is not in file "+filename_1,v_error,m_verbose);
					all_branches_same=false;
				}
			}
			
			Log(m_unique_name+" completed scan for common branches in tree "+akey.first,v_debug,m_verbose);
			
			// having scanned for common branches, did we find at least one?
			if(common_trees.count(akey.first)==0){
				Log(m_unique_name+": File content mismatch! Found no common branches in tree "
					+akey.first+"!",v_error,m_verbose);
			} else if(all_branches_same){
				Log(m_unique_name+": All branches in tree "+akey.first+" match in both files",v_debug,m_verbose);
			}
		} // else this tree is in one file but not the other.
		  // we already alerted the user while comparing file keys.
	}
	
	// we need to do a bit more processing to extract information about the common trees
	ParseSharedTrees();
	
	// Each Tool Execution will get the next entry from all shared trees.
	// But since some trees may have fewer entries, we may run off the end.
	// So track which trees still have entries. Initially all trees are active.
	for(auto&& atree : common_trees){
		active_trees.emplace(atree.first,true);
	}
	
	// for STL containers and all data_instances recursively within them,
	// the addresses currently stored are in fact offsets relative to the position
	// of the base STL container element, which is not yet known.
	// Transfer these addresses into offsets.
	TransferAddressesToOffsets();
	
	gInterpreter->ProcessLine("#include <iostream>");
	
	// if we have an index variable to try to match entries in shared trees,
	// it may not be a branch holding a primitive - it may be embedded in a class
	// (e.g. Header->nevsk). The easiest way to scan through all our data_instances
	// is just to try to compare an entry.
	// TODO one index per shared tree
	if(index_name!=""){
		Log(m_unique_name+" Scanning for index "+index_name,v_debug,m_verbose);
		//for(auto&& apair : common_trees){
			shared_tree& atree = common_trees.begin()->second;
			
			// get the TTree entries
			if( (atree.file1_tree->GetEntry(0)>0) &&
				(atree.file2_tree->GetEntry(0)>0) ){
				
				// scan through branches, recursing through all class members etc
				auto file1_it = atree.file1_branches.begin();
				while(file1_it!=atree.file1_branches.end()){
					Log(m_unique_name+": index scanning branch "+file1_it->first,v_debug,m_verbose);
					branch_structure& file1_branch = file1_it->second;
					branch_structure& file2_branch = atree.file2_branches.at(file1_it->first);
					int less = 0;
					// silence comparison printouts, this is only a test.
					// FIXME even this won't silence things completely, since gInterpreter printouts
					// don't respect 'm_verbose'. Also need to be able to debug print messages about the scan?
					int m_verbose_tmp = m_verbose;
					m_verbose = -1;
					Log(m_unique_name+" Scanning for index "+index_name
					         +", any following inequality printouts can be ignored",v_warning,m_verbose);
					CompareBranchMembers(file1_branch.held_data, file2_branch.held_data, &less);
					Log(m_unique_name+" Scan complete.",v_warning,m_verbose);
					m_verbose = m_verbose_tmp;
					if(file1_index!=nullptr) break; // found it
					++file1_it;
				}
				if(file1_index==nullptr || file2_index==nullptr){
					Log(m_unique_name+" Error! Failed to find index "+index_name,v_error,m_verbose);
					return false;
				}
			} else {
				// no entries in this TTree!?
				Log(m_unique_name+" Error! No entries in shared tree "+atree.file1_tree->GetName(),v_error,m_verbose);
				return false;
			}
		//}   // loop over shared trees
	}
	
	return true;
}

std::string CompareRootFiles::GetBranchType(TBranch* branch){
	TClass* cl=nullptr;
	EDataType edatatype;
	branch->GetExpectedType(cl,edatatype);
	std::string branchtype;
	if(cl==nullptr) branchtype = TDataType::GetTypeName(edatatype);   // branch holds a basic datatype
	else            branchtype = cl->GetName();                       // branch holds a complex object
	branchtype = CppName(branchtype); // un-ROOTish any typenames
	// note that arrays of primitives will also be reported via EDataType as their element type.
	// We could try to use TLeafElement::GetLen and TLeafElement::GetCount,
	// but it's actually easier just to manually parse the branch title.
	std::string title=branch->GetTitle();
	if(title.find('[')!=std::string::npos){
		std::string dims = title.substr(title.find('['),title.find('/')-title.find('['));
		branchtype=branchtype+dims;
	}
	return branchtype;
}

void CompareRootFiles::ParseSharedTrees(){
	for(auto&& atree : common_trees){
		Log(m_unique_name+": getting preliminary information about shared tree "+atree.first,v_debug,m_verbose);
		ParseTree(atree);
	}
}

void CompareRootFiles::ParseTree(std::pair<const std::string, shared_tree> &tree){
	Log(m_unique_name+": Parsing common tree "+tree.first
		+" in file "+tree.second.file1_tree->GetCurrentFile()->GetName(),v_debug,m_verbose);
	// we need to load all branch addresses, so get first entry
	if(tree.second.file1_tree->LoadTree(0)<0 || tree.second.file2_tree->LoadTree(0)<0){
		Log(m_unique_name+" tree "+tree.first+" has no entries to compare!",v_debug,m_verbose);
		return;
	}
	tree.second.file1_tree->GetEntry(0);
	tree.second.file2_tree->GetEntry(0);
	for(auto&& abranch : tree.second.file1_branches){
		// check we have a valid branch address
		if(abranch.second.branch_ptr==nullptr){
			Log(m_unique_name+" no valid address for branch "+abranch.first
				+" in call to ParseTree!",v_error,m_verbose);
			continue;
		}
		Log(m_unique_name+": getting preliminary information about branch "+abranch.first,v_debug,m_verbose);
		bool ok = ParseBranch(abranch.second);
		if(not ok){
			Log(m_unique_name+" Error parsing branch "+abranch.first+" in tree "
				+tree.second.file1_tree->GetName(),v_error,m_verbose);
		}
	}
	Log(m_unique_name+": Parsing common tree "+tree.first
		+" in file "+tree.second.file2_tree->GetCurrentFile()->GetName(),v_debug,m_verbose);
	for(auto&& abranch : tree.second.file2_branches){
		// check we have a valid branch address
		if(abranch.second.branch_ptr==nullptr){
			Log(m_unique_name+" no valid address for branch "+abranch.first
				+" in call to ParseTree!",v_error,m_verbose);
			continue;
		}
		Log(m_unique_name+": getting preliminary information about branch "+abranch.first,v_debug,m_verbose);
		bool ok = ParseBranch(abranch.second);
		if(not ok){
			Log(m_unique_name+" Error parsing branch "+abranch.first+" in tree "
				+tree.second.file2_tree->GetName(),v_error,m_verbose);
		}
	}
}

bool CompareRootFiles::ParseBranch(branch_structure &branch){
//	std::cout<<m_unique_name<<" parsebranch called with branch_ptr = "<<branch.branch_ptr
//	        <<" and type string "<<branch.type_as_string<<std::endl;
//	Log(smessage.str(),v_debug,m_verbose);
	// transfer some basic data from the branch pointer to the top level data instance
	std::string type_as_string = branch.type_as_string.substr(0,branch.type_as_string.find('['));
	Log(m_unique_name+" type stripped of arr dims = "+type_as_string,v_debug,m_verbose);
	branch.held_data.name = std::string(branch.branch_ptr->GetTree()->GetName())+"."+branch.branch_ptr->GetName();
	branch.held_data.type_as_string = type_as_string;
	branch.held_data.branch_ptr = branch.branch_ptr;
	// for primitives and arrays of primitives the branch will be a simple "TBranch",
	// while for classes and STL containers it will be a "TBranchElement" or "TBranchObject"
//	std::cout<<m_unique_name<<" getting branch type from branch pointer "<<branch.branch_ptr<<std::endl;
//	Log(smessage.str(),v_debug,m_verbose);
	std::string branchclass = branch.branch_ptr->ClassName();
	Log(m_unique_name+" branchclass is "+branchclass,v_debug,m_verbose);
	// XXX could probably refactor this into the TBranch, TBranchElement and TBranchObject if/else case below...
	// TBranch::GetAddress and TBranchObject::GetAddress return nullptr, while
	// TBranchElement::GetAddress returns a *pointer to a pointer* to the object
	void** ptrptr = (void**)branch.branch_ptr->GetAddress();
//	std::cout<<m_unique_name<<" branch data address "<<ptrptr<<std::endl;
//	Log(smessage.str(),v_debug,m_verbose);
	if(ptrptr!=nullptr){
		branch.held_data.address = *ptrptr;
	} else {
		// for these we need to go via the TLeaf
		Log(m_unique_name+" getting valuepointer",v_debug,m_verbose);
		TLeaf* lfptr = (TLeaf*)branch.branch_ptr->GetListOfLeaves()->At(0);
		if(branchclass=="TBranch"){
			// for TBranch, which holds primitives, this returns the address of the object directly
			branch.held_data.address = (void*)lfptr->GetValuePointer();
		} else {
			// for TBranchObject, it again holds a pointer to a pointer to the object
			ptrptr = (void**)lfptr->GetValuePointer();
			branch.held_data.address = *ptrptr;
		}
	}
	
	// we could just call ObjectToDataInstance here to parse the branch datatype,
	// but for top level arrays we have the additional possibility of a variable
	// dimension, where the length for each event is specified in another branch.
	// so we'll do a bit manually.
	
	if(branchclass=="TBranch"){
		Log(m_unique_name+" TBranch",v_debug,m_verbose);
		// check if the branch type is an array
		if(branch.type_as_string.find('[')!=std::string::npos){
			Log(m_unique_name+" array type",v_debug,m_verbose);
			// it's an array
			branch.held_data.instance_type = 1; // assume statically sized for now
			// parse those dimensions. could be something like "[ndim][3][2]" as a complex example.
			std::string dimstring = branch.type_as_string.substr(branch.type_as_string.find('['),std::string::npos);
			size_t startpos = 0;
			size_t endpos = 0;
			while(startpos!=std::string::npos){
				endpos=dimstring.find(']',startpos)-1;
				TString nextdim = dimstring.substr(startpos+1,endpos-startpos);
				startpos=dimstring.find('[',endpos);
				// nextdim could be a branch name or a number.
				if(nextdim.IsDigit()){
					branch.held_data.static_dims.push_back(nextdim.Atoi());
				} else {
					branch.held_data.instance_type = 2;  // actually it's dynamically sized
					branch.held_data.dimension_branch = nextdim.Data();
					// get pointer to the dimension integer
					TBranch* b1 = branch.branch_ptr->GetTree()->FindBranch(nextdim.Data());
					if(b1==nullptr){
						Log(m_unique_name+" Error! unable to find branch "+nextdim.Data()
							+" in tree "+branch.branch_ptr->GetTree()->GetName()
							+" in file "+branch.branch_ptr->GetTree()->GetCurrentFile()->GetName()
							+" which should specify dynamic size of item "+branch.held_data.name,
							v_debug,m_verbose);
						return false;
					} else {
						TLeafI* l1=(TLeafI*)b1->GetListOfLeaves()->At(0);
						branch.held_data.dimension_ptr = (int*)l1->GetValuePointer();
					}
				}
			}
			// get element size
			TDataType element_type(CppName(type_as_string));
			TClass* element_class = TClass::GetClass(type_as_string.c_str());
			if(element_type.GetType()>0){
				branch.held_data.item_size = element_type.Size();
			} else if(element_class){
				branch.held_data.item_size = element_class->Size();
			} else {
				Log(m_unique_name+" unable to get size of array elements!",v_error,m_verbose);
			}
			Log(m_unique_name+" element size is "+toString(branch.held_data.item_size),v_debug,m_verbose);
		} else {
			Log(m_unique_name+" basic primitive",v_debug,m_verbose);
			// just a primitive
			branch.held_data.instance_type = 0;
			TClass* cl;
			EDataType e;
			branch.branch_ptr->GetExpectedType(cl,e);
			std::string elementtypename = TDataType::GetTypeName(e);
			elementtypename = CppName(elementtypename);
			TDataType ttype(elementtypename.c_str());
			branch.held_data.item_size = ttype.Size();
		}
	} else if(branchclass=="TBranchElement" || branchclass=="TBranchObject"){
		Log(m_unique_name+" element or object, invoking ObjectToDataInstance",v_debug,m_verbose);
		// let ObjectToDataInstance parse the data type
		// presume that TBranchElement and TBranchObject are compatible for our uses... XXX
		// TBranchObject seems to be used for e.g. a TVector3
		ObjectToDataInstance(branch.held_data);
	} else {
		Log(m_unique_name+" Error! Unknown branch type "+branchclass+" for branch "
			+branch.branch_ptr->GetName()+" in tree "+branch.branch_ptr->GetTree()->GetName(),
			v_error,m_verbose);
		return false;
	}
	Log(m_unique_name+" branch parsed",v_debug,m_verbose);
	
	return true;
}

bool CompareRootFiles::ObjectToDataInstance(data_instance& thedata){
	Log(m_unique_name+" ObjectToDataInstance called",v_debug,m_verbose);
	// build a map of the class members as data instances.
	// get the TClass for this object from its name. Note that TClass::GetClass
	// won't work if we give it a pointer, so strip those
	std::string strippedtype = TClassEdit::ShortType(thedata.type_as_string.c_str(),1);
	Log(m_unique_name+"strippedtype is "+strippedtype,v_debug,m_verbose);
	TClass* cl = TClass::GetClass(strippedtype.c_str());
//	std::cout<<m_unique_name<<" corresponding TClass is "<<cl<<std::endl;
//	Log(smessage.str(),v_debug,m_verbose);
	if(cl==nullptr || strcmp(cl->GetName(),"string")==0){
		// check if primitive type is recognised - try to construct a TDataType from the name
		TDataType basic_type(strippedtype.c_str());
		Log(m_unique_name+" TDataType is "+basic_type.GetTypeName(),v_debug,m_verbose);
		if(basic_type.GetType()>0){
			// it's a valid primitive type
			Log(m_unique_name+" valid primitive type",v_debug,m_verbose);
			thedata.item_size = basic_type.Size();
			Log(m_unique_name+" size "+toString(thedata.item_size)+" bytes",v_debug,m_verbose);
			thedata.instance_type = 0;
			Log(m_unique_name+" instance set to 0",v_debug,m_verbose);
			Log(m_unique_name+" type as string is "+thedata.type_as_string,v_debug,m_verbose);
			// we need to check for arrays manually, they return the same TDataType
			if(thedata.type_as_string.find('[')!=std::string::npos){
				thedata.instance_type = 1;
				Log(m_unique_name+" instance type updated to 1",v_debug,m_verbose);
				std::string dimstring =
					thedata.type_as_string.substr(thedata.type_as_string.find('['),std::string::npos);
				Log(m_unique_name+" dimstring is "+dimstring,v_debug,m_verbose);
				size_t startpos = 0;
				size_t endpos = 0;
				while(startpos!=std::string::npos){
					endpos=dimstring.find(']',startpos)-1;
					std::string nextdim = dimstring.substr(startpos+1,endpos-startpos);
					startpos=dimstring.find('[',endpos);
					thedata.static_dims.push_back(atoi(nextdim.c_str()));
				}
				Log(m_unique_name+" dimensions parsed",v_debug,m_verbose);
			} // else no dimensions, not an array.
		} else if(strcmp(cl->GetName(),"string")==0){
			// string should be treated as a primitive, because the comparison method is just `==`
			thedata.item_size = sizeof(std::string);
			thedata.instance_type = 0;
		} else {
			// not a class, nor a recognised primitive type?
			Log(m_unique_name+" Error! Unknown type '"+thedata.type_as_string+"'!",v_error,m_verbose);
			thedata.instance_type = -1;
			return false;
		}
		
	} else if(IsStlContainer(strippedtype)){
		Log(m_unique_name+" stl container",v_debug,m_verbose);
		// stl containers, even of primitives, will still have a TClass
		thedata.instance_type = 3;
		// extract the type of the contained elements
		TString contained_type_string = TClassEdit::ShortType(thedata.type_as_string.c_str(),8);
		// manual version - this could choke if the allocator is explicitly given in the type, for e.g.
		//TString contained_type_string = thedata.type_as_string.substr(thedata.type_as_string.find('<')+1,
		//              thedata.type_as_string.find_last_of('>') - thedata.type_as_string.find('<') - 1);
		// trim leading and trailing whitespace
		//contained_type_string.Remove(TString::kBoth, ' ');
		Log(m_unique_name+" element type string is "+contained_type_string,v_debug,m_verbose);
		
		TDataType element_type(contained_type_string.Data());
		TClass* element_class = TClass::GetClass(contained_type_string.Data());
		if(element_type.GetType()>0){
			thedata.item_size = element_type.Size();
		} else if(element_class){
			thedata.item_size = element_class->Size();
		} else {
			Log(m_unique_name+" unable to get size of STL container elements!",v_error,m_verbose);
		}
		Log(m_unique_name+" element size is "+toString(thedata.item_size),v_debug,m_verbose);
		
		// parse the contained element types, which may themselves be complex objects
		thedata.contained_type = new data_instance;
//		std::cout<<m_unique_name<<" member data_instance at "<<thedata.contained_type<<std::endl;
//		Log(smessage.str(),v_debug,m_verbose);
		thedata.contained_type->name = thedata.name + ".at(*)";
		thedata.contained_type->type_as_string = contained_type_string;
		thedata.contained_type->branch_ptr = thedata.branch_ptr;
		// we don't know the address of elements of an stl container they're filled,
		// so we'll have to recursively update this and all members thereof at runtime.
//		std::cout<<m_unique_name<<" propagated branchptr is "<<thedata.contained_type->branch_ptr<<std::endl;
//		Log(smessage.str(),v_debug,m_verbose);
		ObjectToDataInstance(*thedata.contained_type);
		Log(m_unique_name+"returning from STL container ObjectToDataInstance",v_debug,m_verbose);
		
	} else {
		Log(m_unique_name+" class instance",v_debug,m_verbose);
		// else it's a proper class
		thedata.instance_type = 4;
		
		// TClonesArray... we'll need to handle those the same way as STL containers
		if(strippedtype=="TClonesArray"){
			// A TClonesArray stores a fixed data type, but the type isn't fixed to the class
			// so we need to get it from the instance
			TClonesArray* carr = (TClonesArray*)thedata.address;
			if(carr==nullptr){
				Log(m_unique_name+" Error! nullptr when parsing TClonesArray held type name!",v_error,m_verbose);
				return false;
			}
			// get the held type name from the clonesarray
			thedata.contained_type = new data_instance;
			thedata.contained_type->name = thedata.name + ".At(*)";
			TString strippedclonetype = TClassEdit::ShortType(carr->GetClass()->GetName(),1);
			thedata.contained_type->type_as_string = strippedclonetype;
			thedata.contained_type->branch_ptr = thedata.branch_ptr;
			
			// grab the size
			TDataType element_type(strippedclonetype.Data());
			TClass* element_class = TClass::GetClass(strippedclonetype.Data());
			if(element_type.GetType()>0){  // yes, 0 is "notype"
				thedata.item_size = element_type.Size();
			} else if(element_class){
				thedata.item_size = element_class->Size();
			} else {
				Log(m_unique_name+" unable to get size of STL container elements!",v_error,m_verbose);
			}
			
			// parse information about the cloned type
			ObjectToDataInstance(*thedata.contained_type);
			return true;
		}
		// TObjArray... well. The contained object type can vary from element to element.
		// we'll parse the class at runtime for each element... maybe inefficient.
		// we could also just compare as if it were a primitive...? ie directly with `==`
		if(strippedtype=="TObjArray"){
			thedata.contained_type = new data_instance;
			thedata.contained_type->name = thedata.name + ".At(*)";
			thedata.contained_type->branch_ptr = thedata.branch_ptr;
			// i think this is sufficient...?
		}
		
		// loop over the members and convert their types for comparison too
		// we should account for base class members too
		std::vector<std::pair<TDataMember*, Long_t>> members; // members and their base class offsets
		int total_members = GetAllClassMembers(cl, members);
		Log(m_unique_name+" processing "+toString(total_members)+" members of class "+cl->GetName(),v_debug,m_verbose);
		for(int member_i=0; member_i<total_members; ++member_i){
			TDataMember* member = members.at(member_i).first;
			Long_t member_offset = members.at(member_i).second;
			Log(m_unique_name+" member "+member->GetName()+" has offset "+toString(member_offset),v_debug,m_verbose);
			//if(member->GetOffset()==0) continue;    // skip dummy member such as fgIsA - N/A when using RealData
			
			// make a new member instance to hold info about this member variable
			data_instance member_instance;
			
			member_instance.name = thedata.name + "." + member->GetName();
			// we need to use "GetFullTypeName" not just "GetTypeName"
			// to identify members that are stored as pointers!
			member_instance.type_as_string = member->GetFullTypeName();
			Log(m_unique_name+"member type is "+member_instance.type_as_string,v_debug,m_verbose);
			// note that often, even members that are not pointers as defined by the class header,
			// are nonetheless stored as pointers!
			if(member->IsaPointer()){
				Log(m_unique_name+" pointer, double dereferencing address",v_debug,m_verbose);
				// in this case we will obtain the address of a pointer to the object
				void** ptrptr = (void**)((char*)thedata.address + member_offset);
				// we have a gotcha here: if this parent object is an element of a container
				// then its address (thedata.address) is actually an offset relative to the
				// container address, which isn't known until the container is filled by an entry.
				// in that case, we can't (yet) dereference the resulting pointer without segging.
				// so, in such cases, we'll record the address of the pointer, and note it as such.
				if(member_instance.name.find(".at(*)")!=std::string::npos ||
				   member_instance.name.find(".At(*)")!=std::string::npos){
					// it's a container element
					Log(m_unique_name+" this is a member of a container element, can't dereference yet",v_debug,m_verbose);
					// XXX FIXME this does not work!
					member_instance.address = ptrptr;
					member_instance.is_ptr = true;
				} else {
					// safe to reference
					member_instance.address = *ptrptr;
				}
			} else {
				member_instance.address = (void**)((char*)thedata.address + member_offset);
			}
			
			// we could now just call ObjectToDataInstance to convert the type,
			// but the fact that we have a TDataMember instead of a TClass
			// offers us easier access to information, so let's use it.
			// FIXME is there really any advantage here? What are we gaining?
			// this is separating out STL containers, but not e.g slonesarrays etc
			// i think in particular it's because we have GetMaxIndex for handling arrays
			// but STL containers we can bundle in with classes.
			
			// TDataMember::GetDataType returns a TDataType for basic types only
			if(member->GetDataType()!=nullptr){
				Log(m_unique_name+"basic type",v_debug,m_verbose);
				// a basic type
				member_instance.instance_type = 0;
				member_instance.item_size = member->GetUnitSize();
				
				// also need to check for arrays. Only static dimensions here, no branches within a class.
				int array_dim_i=0;
				while(member->GetMaxIndex(array_dim_i)>=0){
					Log(m_unique_name+"next array dimension "
					    +toString(member->GetMaxIndex(array_dim_i)),v_debug,m_verbose);
					member_instance.static_dims.push_back(member->GetMaxIndex(array_dim_i));
					++array_dim_i;
					member_instance.instance_type = 1;
				}
			} else {
				Log(m_unique_name+" class type",v_debug,m_verbose);
				// if it's not a primitive and it's not an STL container, it's a class.
				member_instance.instance_type = 4;
				// if it's a class, it has its own members - parse those too...
				Log(m_unique_name+" calling ObjectToDataInstance for member "+member_instance.name,v_debug,m_verbose);
				ObjectToDataInstance(member_instance);
				Log(m_unique_name+" returned from ObjectToDataInstance for member "
				    +member_instance.name,v_debug,m_verbose);
			}
			
			// add this member to the set of contents to compare
			thedata.members.emplace(member->GetName(), member_instance);
			
		} // end loop over members of this class
	} // end if this is a class.
	
	return true;
}

int CompareRootFiles::GetAllClassMembersOld(TClass* cl, Long_t offset, std::vector<std::pair<TDataMember*, Long_t>> &members){
	// we need to recurse through all base classes, and their base classes, and so on.
	// start with direct members
	int skipped=0;
	for(int member_i=0; member_i<cl->GetListOfDataMembers()->GetEntries(); ++member_i){
		TDataMember* member = (TDataMember*)cl->GetListOfDataMembers()->At(member_i);
		Log(m_unique_name+" member "+member->GetName()+" of class "+cl->GetName()
		         +" has offset "+toString(member->GetOffset()),v_debug,m_verbose);
		if(member->GetOffset()==0){
			++skipped;
			continue;    // skip dummy member such as fgIsA
		}
		members.emplace_back(std::pair<TDataMember*,Long_t>{member, offset+member->GetOffset()});
	}
	// scan base classes
	for(int base_i=0; base_i<cl->GetListOfBases()->GetEntries(); ++base_i){
		TBaseClass* bcl = (TBaseClass*)cl->GetListOfBases()->At(base_i); // get offset from parent class
		Long_t baseoffset = bcl->GetDelta();
		if(baseoffset<0){
			Log(m_unique_name+" Error! Negative base class offset for base class "+bcl->GetName()
				+" within class "+cl->GetName(),v_error,m_verbose);
			continue; // i guess we'll skip these base members....
		}
		TClass* bclcl = (TClass*)bcl->GetClassPointer();
		skipped += GetAllClassMembersOld(bclcl, baseoffset, members);
	}
	// another way might be to use GetListOfRealData instead of GetListOfDataMembers,
	// which gets us a list of TRealData*, from each of which we can get the corresponding
	// TDataMember* and the offset... supposedly this list contains all persitent members,
	// accounting for base classes and skips non-persistent data members....
	return skipped;
}

int CompareRootFiles::GetAllClassMembers(TClass* cl, std::vector<std::pair<TDataMember*, Long_t>> &members){
	// RealData lists all persistent members of a class, including those from base classes.
	cl->BuildRealData();
	for(int member_i=0; member_i<cl->GetListOfRealData()->GetEntries(); ++member_i){
		TRealData* nextrealdata = (TRealData*)cl->GetListOfRealData()->At(member_i);
		TClass* cl = nextrealdata->IsA();

		
		TDataMember* nextmember = nextrealdata->GetDataMember();
		// it seems as though the following are part of ROOT's internal streamer functionality,
		// and are very unlikely to be anything we care about being different
		// (we are of course assuming a user has not decided to use them for meaningful data,
		//  or generated a class with members that store meaningful data with these names)
		if(strcmp(nextmember->GetName(),"fBits")==0 || strcmp(nextmember->GetName(),"fUniqueID")==0) continue;
		members.emplace_back(std::pair<TDataMember*, Long_t>{nextmember, nextrealdata->GetThisOffset()});
	}
	return members.size();
}

void CompareRootFiles::TransferAddressesToOffsets(){
	for(auto&& apair : common_trees){
		shared_tree& atree = apair.second;
		auto file1_it = atree.file1_branches.begin();
		while(file1_it!=atree.file1_branches.end()){
			branch_structure& file1_branch = file1_it->second;
			branch_structure& file2_branch = atree.file2_branches.at(file1_it->first);
			TransferBranchAddresses(file1_branch.held_data);
			TransferBranchAddresses(file2_branch.held_data);
			++file1_it;
		}
	}
}

void CompareRootFiles::TransferBranchAddresses(data_instance &thedata){
	// the contained_type is a data_instance representing a container element
	// whose address is not yet known.
	// transfer its address and all its childrens addresses to offsets
	if(thedata.contained_type!=nullptr){
		DoTransferAddresses(*thedata.contained_type);
	}
}

void CompareRootFiles::DoTransferAddresses(data_instance &thedata){
	thedata.offset = (long)thedata.address;
	thedata.address = nullptr;
	for(auto&& amember : thedata.members){
		DoTransferAddresses(amember.second);
	}
	if(thedata.contained_type!=nullptr){
		DoTransferAddresses(*thedata.contained_type);
	}
}


//////////////////////////////////////////////////////////////
//                          Execute                         //
//////////////////////////////////////////////////////////////

bool CompareRootFiles::Execute(){
	
	// for all common trees, compare all common branches.
	bool more_entries=false;
	Log(m_unique_name+" doing comparison "+toString(entry_number),v_debug,m_verbose);
	bool all_equal=true;
	
	for(auto&& apair : common_trees){
		// skip this shared tree if we already hit the end of it.
		if(!active_trees.at(apair.first)) continue;
		
		shared_tree& atree = apair.second;
		
		// find the next entry to compare from each TTree
		Log(m_unique_name+" Looking for next entry with matching indices "+index_name,v_debug,m_verbose);
		get_ok = GetNextMatchingEntries(apair);
		// if there are no more entries in this tree, mark it as inactive
		if(!get_ok) active_trees.at(apair.first) = false;
		more_entries |= get_ok;
		
		entry_number_1 = entry_numbers_1.at(apair.first);
		entry_number_2 = entry_numbers_2.at(apair.first);
		entry_string = toString(entry_number_1)+"/"+toString(entry_number_2);
		
		if(get_ok){
			Log(m_unique_name+": comparing tree " + atree.file1_tree->GetName() + " entry "+toString(entry_number_1)
			    + " in file " + atree.file1_tree->GetCurrentFile()->GetName() + " with entry "
			    + toString(entry_number_2)+" in file "+atree.file2_tree->GetCurrentFile()->GetName(),
			    v_debug,m_verbose);
			
			// compare branches in this TTree entry
			auto file1_it = atree.file1_branches.begin();
			bool all_branches_equal=true;
			while(file1_it!=atree.file1_branches.end()){
				Log(m_unique_name+": comparing branch "+file1_it->first,v_debug,m_verbose);
				branch_structure& file1_branch = file1_it->second;
				branch_structure& file2_branch = atree.file2_branches.at(file1_it->first);
				all_branches_equal &= CompareBranchMembers(file1_branch.held_data, file2_branch.held_data);
				++file1_it;
			}
			all_equal &= all_branches_equal;
			if(all_branches_equal) Log(m_unique_name+": Entry "+toString(entry_number_1)+" in file "+
					      atree.file1_tree->GetCurrentFile()->GetName()+" Tree "+atree.file1_tree->GetName()+
					      " is identical to entry "+toString(entry_number_2)+" in file "+
					      atree.file2_tree->GetCurrentFile()->GetName(),
					      v_warning,m_verbose);
			
		} // no more entries in this TTree
	}
	Log(m_unique_name+": "+toString(entry_number)+" entries compared so far",v_debug,m_verbose);
	if(all_equal){
		Log(m_unique_name+": Comparison "+toString(++matching_entries)+" found matching entries so far",
		    v_message,m_verbose);
	} else {
		Log(m_unique_name+": Comparison "+toString(++mismatching_entries)+" found mismatching entries so far",
		    v_message,m_verbose);
	}
	++entry_number;
	
	// break once we have no more entries in any trees
	if(not more_entries){
		Log(m_unique_name+" Finished all trees, ending toolchain.",v_message,m_verbose);
		m_data->vars.Set("StopLoop",1);
	}
	if(entry_number>=max_entries && max_entries>0){
		Log(m_unique_name+" Processed max entries, ending toolchain.",v_message,m_verbose);
		m_data->vars.Set("StopLoop",1);
	}
	
	return true;
}


bool CompareRootFiles::CompareBranchMembers(data_instance &branch1, data_instance &branch2, int* less){
	// if we've not yet found the index branches, check if this is them, and note the instances if so.
	if(less!=nullptr && file1_index==nullptr){
		Log(m_unique_name+" Checking data_instance "+branch1.name+" against index_name "+index_name,v_debug,m_verbose);
		if(branch1.name==index_name){
			smessage <<m_unique_name<<" Match! Noting file1_index = "<<&branch1<<", file2_index = "<<&branch2;
			Log(smessage.str(),v_debug,m_verbose);
			file1_index = &branch1;
			file2_index = &branch2;
		}
	}
	
	// to compare the two branches we need to compare their held data
	// quick sanity check we're comparing things of equivalent type
	if(branch1.name!=branch2.name || branch1.type_as_string != branch2.type_as_string ||
		branch1.instance_type!=branch2.instance_type || branch1.instance_type < 0){
		Log(m_unique_name+" Error comparing branches! Comparison called between "
			+branch1.name+", of type "+branch1.type_as_string+", determined to be of instance_type "
			+toString(branch1.instance_type)+" and branch "+branch2.name+", of type "+branch2.type_as_string
			+", determined to be of instance_type "+toString(branch2.instance_type),v_error,m_verbose);
		return false;
	}
	// to compare the objects we'll also need their addresses.
	// in some cases (pointer members of container objects) we only have a pointer to their address
	// so we'll need to de-reference it now to get the underlying object
	void* branch1_add = branch1.address;
	void* branch2_add = branch2.address;
	if(branch1.is_ptr){
		//Log(m_unique_name+" comparing objects by pointers",v_debug,m_verbose);
		Log(m_unique_name+" skipping comparison of "+branch1.name
		    +" as we do not yet support pointers to types (such as char*) in containers",v_debug,m_verbose);
		// this uhhh doesn't work.
		/*
		void** branch1_ptr = (void**)branch1_add;
		void** branch2_ptr = (void**)branch2_add;
		branch1_add = *branch1_ptr;
		branch2_add = *branch2_ptr;
		std::cout<<"dereferenced addresses are "<<branch1_add<<", and "<<branch2_add<<std::endl;
		*/
		return true; // TODO FIXME FIXME
	}
	
	// ok, the way we compare will depend on the instance_type
	bool are_equal=true;
	std::string stype = branch1.type_as_string;
	int dynamic_dim;
	switch(branch1.instance_type){
		case 0: {
			// basic type
//			std::cout << m_unique_name << " basic datatype comparison of type " << stype
//			         << " at addresses " << branch1_add << " and " << branch2_add
//			         << " for branch " << branch1.name<<std::endl;
//			Log(smessage.str(),v_debug,m_verbose);
			std::string vname = GetVarName(stype);
			TString cmd = TString::Format("bool* eq_ptr = (bool*)%p; ", (void*)&are_equal);
			cmd += TString::Format("%s* %s_1 = (%s*)%p; ", 
			         stype.c_str(), vname.c_str(), stype.c_str(),(void*)branch1_add);
			cmd += TString::Format("%s* %s_2 = (%s*)%p; ",
			         stype.c_str(), vname.c_str(), stype.c_str(), (void*)branch2_add);
			if(stype!="string" && stype!= "bool"){
				// assume numeric, and that this is suitable?
				cmd += TString::Format("*eq_ptr = (std::fabs(*%s_1-*%s_2) < %f); ",
				         vname.c_str(), vname.c_str(), ftolerance); // FIXME tolerance
			} else {
				cmd += TString::Format("*eq_ptr = (*%s_1==*%s_2); ",
				         vname.c_str(), vname.c_str());
			}
			Log(m_unique_name+" Processing "+cmd,v_debug,m_verbose);
			gInterpreter->ProcessLine(cmd);
			if(not are_equal){
				// mismatching item!
				cmd  = TString::Format("std::cout<<(*%s_1)<<\" != \"<<(*%s_2)<<std::endl;",
				               vname.c_str(), vname.c_str());
				Log(m_unique_name+" Processing "+cmd,v_debug,m_verbose);
				gInterpreter->ProcessLine(cmd);  // FIXME m_verbose
				Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name,
					v_error,m_verbose);
				
				// if requested, also return which is less
				if((less!=nullptr) && (stype!="string") && (stype!="bool")){
					// try to also return which is less
					cmd = TString::Format("*eq_ptr = (*%s_1)<(*%s_2); ",
						     vname.c_str(), vname.c_str()); // FIXME tolerance
					Log(m_unique_name+" Processing "+cmd,v_debug,m_verbose);
					gInterpreter->ProcessLine(cmd);
					*less = (are_equal) ? 1 : 0;
					are_equal = false; // restore this, gets overwritten by the second ProcessLine call
				}
			}
			return are_equal;
		};
		// jump to case 2 (case 1 is next)
		case 2: {
			// dynamic array type
//			std::cout  << m_unique_name << " dynamic array comparison with element type " << stype
//			          << " with array addresses " << branch1_add << " and " << branch2_add
//			          << " for branch " << branch1.name<<std::endl;
//			Log(smessage.str(),v_debug,m_verbose);
			
			if(branch1.dimension_ptr==nullptr || branch2.dimension_ptr==nullptr){
				// if we were unable to get the pointer to the dimension branch, bail out
				return false;
			}
			int dynamic_dim_1 = *branch1.dimension_ptr;
			int dynamic_dim_2 = *branch2.dimension_ptr;
			if(dynamic_dim_1<0){
				Log(m_unique_name+" Error! Branch "+branch1.dimension_branch
					+" in tree "+branch1.branch_ptr->GetTree()->GetName()
					+" in file "+branch1.branch_ptr->GetTree()->GetCurrentFile()->GetName()
					+" which should specify dynamic size of item "+branch1.name
					+" reports a negative size for entry "+toString(entry_number_1),v_error,m_verbose);
				return false;
			} else if(dynamic_dim_2<0){
				Log(m_unique_name+" Error! Branch "+branch2.dimension_branch
					+" in tree "+branch2.branch_ptr->GetTree()->GetName()
					+" in file "+branch2.branch_ptr->GetTree()->GetCurrentFile()->GetName()
					+" which should specify dynamic size of item "+branch2.name
					+" reports a negative size for entry "+toString(entry_number_2),v_error,m_verbose);
				return false;
			} else if(dynamic_dim_1!=dynamic_dim_2){
				// dynamic dimension mismatch
//				std::cout<<dynamic_dim_1<<" != "<<dynamic_dim_2<<std::endl; // XXX m_verbose
				Log(toString(dynamic_dim_1)+" != "+toString(dynamic_dim_2),v_error,m_verbose);
				Log(m_unique_name+"Mismatch! Entry "+entry_string+" item "
				    +branch2.branch_ptr->GetTree()->GetName()+"::"+branch2.dimension_branch
				    +", which defines the length of the array in branch "+branch2.branch_ptr->GetName(),
				    v_error,m_verbose);
				dynamic_dim = std::min(dynamic_dim_1,dynamic_dim_2);
				Log(m_unique_name+" Comparing fewest ("+toString(dynamic_dim)+" entries",v_message,m_verbose);
			} else {
				// sizes are the same.
				dynamic_dim = dynamic_dim_1;
			}
			// check we have something to compare. An array length of 0 is not an error, here.
			if(dynamic_dim==0) return true;
			// otherwise, we next need to compare the static dimensions, before then comparing elements.
			// That code is the same as for the fixed-size array case, so rather than duplicating it,
			// allow this case structure to 'fall through' to the next one.
		};
		case 1: {
			// static array type
//			std::cout  << m_unique_name << " static array comparison with element type " << stype
//			          << " and array addresses " << branch1_add << " and " << branch2_add
//			          << " for branch " << branch1.name<<std::endl;
//			Log(smessage.str(),v_debug,m_verbose);
			
			// sanity check: compare the dimensionality. This should be fixed at Initialise.
			std::vector<int> static_dims = branch1.static_dims;
			if(branch1.static_dims!=branch2.static_dims){
				// if not, take the smaller of the two
				static_dims.resize(std::min(branch1.static_dims.size(),branch2.static_dims.size()));
				for(int i=0; i<static_dims.size(); ++i){
					static_dims.at(i) = std::min(branch1.static_dims.at(i),branch2.static_dims.at(i));
					if(static_dims.at(i)==0){
						if(branch1.static_dims.at(i)==0){
							Log(m_unique_name+" Warning! Static dimension "+toString(i)+" for array "+branch1.name
								+" in file "+branch1.branch_ptr->GetTree()->GetCurrentFile()->GetName()
								+" is 0!",v_error,m_verbose);
						}
						if(branch2.static_dims.at(i)==0){
							Log(m_unique_name+" Warning! Static dimension "+toString(i)+" for array "+branch2.name
								+" in file "+branch2.branch_ptr->GetTree()->GetCurrentFile()->GetName()
								+" is 0!",v_error,m_verbose);
						}
						return false;
					}
				}
				Log(m_unique_name+" Warning! Different static dimensions of item "+branch1.name
					+", will only be comparing the smallest common dimensions!",v_debug,m_verbose);
			}
			// if it's a static sized array, ensure we have at least one valid dimension
			if(branch1.instance_type==1 && static_dims.size()==0){
				// nothing to do, at least one of the branches has a 0 dimension array
				Log(m_unique_name+" branch "+branch1.name+" appears to be a static array but has 0 size?",
				    v_error,m_verbose);
				return false;
			}
			// if we're dealing with a dynamic array, combine the static and dynamic dimensions
			std::vector<int> array_dims;
			if(branch1.instance_type==2){
				array_dims.push_back(dynamic_dim); // dynamic dimension is always the first
			}
			array_dims.insert(array_dims.end(),static_dims.begin(), static_dims.end());
			
			// check we have a consistent item size, which we need to derive the correct
			// addresses of the constituent elements.
			if(branch1.item_size!=branch2.item_size || branch1.item_size<=0){
				Log(toString(branch1.item_size)+" != "+toString(branch2.item_size),v_error,m_verbose);
				Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name
				    +" item_size is different!",v_error,m_verbose);
				return false;
			}
			// flatten out the array to make the loop over elements easier
			int total_size = std::accumulate(array_dims.begin(), array_dims.end(), 1, std::multiplies<int>());
			for(int i=0; i<total_size; ++i){
				// compare element i
				void* address_1 = (void*)((char*)branch1_add + (i*branch1.item_size));
				void* address_2 = (void*)((char*)branch2_add + (i*branch1.item_size));
				bool is_equal;
				std::string vname = GetVarName(stype);
				TString cmd = TString::Format("bool* eq_ptr = (bool*)%p; ", (void*)&is_equal);
				cmd += TString::Format("%s* %s_1 = (%s*)%p; ",
				        stype.c_str(), vname.c_str(), stype.c_str(), (void*)address_1);
				cmd += TString::Format("%s* %s_2 = (%s*)%p; ",
				        stype.c_str(), vname.c_str(), stype.c_str(), (void*)address_2);
				cmd += TString::Format("*eq_ptr = (std::fabs(*%s_1-*%s_2) < %f); ",
				        vname.c_str(), vname.c_str(), ftolerance);
				Log(m_unique_name+" Processing "+cmd,v_debug,m_verbose);
				gInterpreter->ProcessLine(cmd);
				// FIXME the age old problem of floating point comparison.
				// ideally it would be relative, so maybe like 1% of the value.
				// can we: load a function into the gInterpreter to do this?
				// and/or: embed a suitable function here?
				if(not is_equal){
					// mismatching element!
					// FIXME we have to print the values from within the interpreter,
					// since only there are the types understood to achieve proper streaming,
					// but this does not (yet) take a m_verbose argument.
					cmd += TString::Format("std::cout<<(*%s_1)<<\" != \"<<(*%s_2)<<std::endl;",
						                   vname.c_str(), vname.c_str());
					Log(m_unique_name+" Processing "+cmd,v_debug,m_verbose);
					gInterpreter->ProcessLine(cmd);
					// convert index back to multidimensional array format for printing
					std::string dimstring="";
					int remdr = i;
					for(int j=1; j<array_dims.size()+1; ++j){
						int nextdimsize = 
							std::accumulate(array_dims.begin()+j, array_dims.end(), 1, std::multiplies<int>());
						dimstring.append(std::string("[")+std::to_string(remdr/nextdimsize)+"]");
						remdr = (remdr % nextdimsize);
					}
					Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name+" index "
						+dimstring,v_error,m_verbose);
					are_equal=false;  // flag that the entry is not entirely identical
				}
			} // end of loop over elements
			return are_equal;
		}  // end case static array
		case 3: {
			// stl container.
//			std::cout  << m_unique_name << " stl container comparison with type " << stype
//			          << " at addresses " << branch1_add << " and " << branch2_add
//			          << " for branch " << branch1.name<<std::endl;
//			Log(smessage.str(),v_debug,m_verbose);
			
			// First thing is we need to get the size
			size_t size1, size2;
			/*
			std::string vname = GetVarName(stype);
			TString cmd;
			cmd  = TString::Format("%s* %s_1 = (%s*)%p; ",
			                       stype.c_str(), vname.c_str(), stype.c_str(), (void*)branch1_add);
			cmd += TString::Format("%s* %s_2 = (%s*)%p; ",
			                       stype.c_str(), vname.c_str(), stype.c_str(), (void*)branch2_add);
			cmd += TString::Format("size_t* s1 = (size_t*)%p; ", (void*)&size1);
			cmd += TString::Format("size_t* s2 = (size_t*)%p; ", (void*)&size2);
			cmd += TString::Format("*s1 = %s_1->size(); *s2 = %s_2->size();",
			               vname.c_str(), vname.c_str());
			std::cout<<"Processing "<<cmd<<std::endl;
			gInterpreter->ProcessLine(cmd);
			std::cout<<"did processline for size"<<std::endl;
			*/
			// the issue with using the gInterpreter is it may not know about contained classes
			// so "std::container<MyClass>* thing = ..." will fail.
			// But instead we can get access to the size method via TMethodCall
			TMethodCall* sizecaller = GetSizeCaller(branch1);
			if(sizecaller==nullptr){
				Log(m_unique_name+" Null TMethodCall for getting size of type "+stype,v_error,m_verbose);
				return false;
			} else if(not sizecaller->IsValid()){
				Log(m_unique_name+" invalid TMethodCall for getting size of type "+stype,v_error,m_verbose);
				return false;
			}
			double retval;  // it returns a double, not an int...
			sizecaller->Execute(branch1_add, retval);
			Log(m_unique_name+" retval1 is "+toString(retval),v_debug,m_verbose);
			size1 = retval;
			sizecaller->Execute(branch2_add, retval);
			Log(m_unique_name+" retval2 is "+toString(retval),v_debug,m_verbose);
			size2 = retval;
			
			if(size1!=size2){
				printf("%lu != %lu\n",size1,size2);
				Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name
					+" vector sizes are different!",v_error,m_verbose);
			}
			
			// currently for comparing elements we use at to get the element or address of the element
			TMethodCall* atcaller = GetAtCaller(branch1);
			if(atcaller==nullptr){
				Log(m_unique_name+" Null TMethodCall for getting element from type "+stype,v_error,m_verbose);
				return false;
			} else if(not atcaller->IsValid()){
				Log(m_unique_name+" invalid TMethodCall for getting element from type "+stype,v_error,m_verbose);
				return false;
			}
			
			// we'll also need the type of the contained elements
			std::string containedtypename = branch1.contained_type->type_as_string.c_str();
			// convert to generic c++ typename if necessary
			containedtypename = CppName(containedtypename);
			TDataType basic_type(containedtypename.c_str());
			int type_num = basic_type.GetType();
			Log(m_unique_name+" element of type "+stype+" is of type "+basic_type.GetTypeName()
			      +" corresponding to EDataType "+toString(type_num),v_debug,m_verbose);
			
			// compare as many as we have
			size_t minsize = std::min(size1,size2);
			Log(m_unique_name+" comparing "+toString(minsize)+" elements",v_debug,m_verbose);
			bool is_equal;
			for(int i=0; i<minsize; ++i){
				Log(m_unique_name+" comparing next element "+toString(i),v_debug,m_verbose);
				void* add1; void* add2;
				/*
				// XXX can't do the below for types not recognised by gInterpreter
				// as the variable vname can't be made, since it is of unknown type.
				// not sure i trust our method for getting sizes, so let's do it a bit more robust
				cmd  = TString::Format("void** v1 = (void**)%p; ", (void*)&add1);
				cmd += TString::Format("void** v2 = (void**)%p; ", (void*)&add2);
				cmd += TString::Format("*v1 = (void*)&((*%s_1)[%d]); ", vname.c_str(), i);
				cmd += TString::Format("*v2 = (void*)&((*%s_2)[%d]); ", vname.c_str(), i);
				std::cout<<"Processing "<<cmd<<std::endl;
				gInterpreter->ProcessLine(cmd);
				std::cout<<"did processline to get element addresses"<<std::endl;
				*/
				// instead we'll call container::at, which returns either the value (for basics)
				// or a pointer to the value (for complex types). This will work for
				// vector, array, map, and deque, but not list, etc that do not have 'at'.
				// we could make an iterator, but then we'd:
				// 1. need to instantiate an instance, and set it to the container begin
				// 2. we can use a TMethodCall to invoke operator++, which returns
				// either the element, or a pointer to our iterator (which has advanced)
				// 3. we then need to get the address from the iterator. I don't think we can do this.
				
				// FIXME replace the use of 'at' which is specialised for some types (map, vector)
				// with the use of an iterator, which would work for all stl containers.
				// need to instantiate a TMethodCall for container<myclass>::begin()
				// use it to get a pointer to an instance of container<class>::iterator,
				// then make a TMethodCall for container<class>::iterator operator++,
				// and use the returned result to update the iterator,
				// and finally dereference the iterator..............
				// or just do it like above with ProcessLine calls.
				
				atcaller->ResetParam();    // must call this first
				Long_t index = i;          // SetParam only takes limited types, have to convert to Long_t
				atcaller->SetParam(index); // set the index to retrieve
				
				// hacky but use a temporary to hold basic return types
				Long_t tmplong_1, tmplong_2;
				Double_t tmpdouble_1, tmpdouble_2;
				Float_t tmpfloat_1, tmpfloat_2;
				if(type_num >0 && type_num!=7 && type_num!= 20){ // basic, except char* and void...
					Log(m_unique_name+" basic datatype for TMethodCall of at",v_debug,m_verbose);
					if(type_num==5  || type_num == 19){
						// Float_t, Float16_t
						Log(m_unique_name+" float-like number",v_debug,m_verbose);
						atcaller->Execute(branch1_add, tmpdouble_1);
						atcaller->Execute(branch2_add, tmpdouble_2);
//						std::cout << m_unique_name <<" tmpdouble_1 is "<<tmpdouble_1<<" at "<<&tmpdouble_1<<" "
//						         <<"tmpdouble_2 is "<<tmpdouble_2<<" at "<<&tmpdouble_2<<std::endl;
//						Log(smessage.str(),v_debug,m_verbose);
						tmpfloat_1 = static_cast<float>(tmpdouble_1);
						tmpfloat_2 = static_cast<float>(tmpdouble_2);
//						std::cout << m_unique_name << " tmpfloat_1 is "<<tmpfloat_1<<" at "<<&tmpfloat_1<<" "
//						         <<"tmpfloat_2 is "<<tmpfloat_2<<" at "<<&tmpfloat_2<<std::endl;
//						Log(smessage.str(),v_debug,m_verbose);
						add1=(void*)&tmpfloat_2;
						add2=(void*)&tmpfloat_2;
					} else if(type_num == 8 || type_num == 9){
						// Double_t, Double32_t
						Log(m_unique_name+" double-like number",v_debug,m_verbose);
						atcaller->Execute(branch1_add, tmpdouble_1);
						atcaller->Execute(branch2_add, tmpdouble_2);
//						std::cout << m_unique_name <<" tmpdouble_1 is "<<tmpdouble_1<<" at "<<&tmpdouble_1<<" "
//						         <<"tmpdouble_2 is "<<tmpdouble_2<<" at "<<&tmpdouble_2<<std::endl;
//						Log(smessage.str(),v_debug,m_verbose);
						add1=(void*)&tmpdouble_1;
						add2=(void*)&tmpdouble_2;
					} else {
						// everything else is integer-like
						atcaller->Execute(branch1_add, tmplong_1);
						atcaller->Execute(branch2_add, tmplong_2);
						add1=(void*)&tmplong_1;
						add2=(void*)&tmplong_2;
					}
				} else {
					// complex held type, return will be a pointer
					atcaller->Execute(branch1_add, tmplong_1);
					atcaller->Execute(branch2_add, tmplong_2);
					add1=(void*)tmplong_1;
					add2=(void*)tmplong_2;
				}
				
				if(add1==nullptr || add2==nullptr){
					Log(m_unique_name+" Error! nullptr returned from STL comparison for element "
						+toString(i)+" of branch "+branch1.name,v_error,m_verbose);
						return false; // assume we're not going to get good addresses for any further elements
				}
				
				// pass these addresses to the data_instance describing element types
				branch1.contained_type->address = add1;
				branch2.contained_type->address = add2;
				// we need to recursively update all member addresses too,
				// since they currently only hold offsets to the parent object
				UpdateAddresses(*branch1.contained_type);
				UpdateAddresses(*branch2.contained_type);
				
				// compare these elements
				is_equal = CompareBranchMembers(*branch1.contained_type, *branch2.contained_type, less);
				if(not is_equal){
					// mismatching element!
					Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name+" index "
						+toString(i),v_error,m_verbose);
					are_equal=false;  // flag that the entry is not entirely identical
				}
			} // move to next element
			return are_equal;
		};
		case 4: {
			// class case.
//			std::cout  << m_unique_name << " class comparison with type " << stype
//			          << " at addresses " << branch1_add << " and " << branch2_add
//			          << " for branch " << branch1.name<<std::endl;
//			Log(smessage.str(),v_debug,m_verbose);
			
			// first check for TClonesArray, basically the same as STL container version.
			// could probably refactor this code a lot better....
			if(stype=="TClonesArray" || stype=="TObjArray"){
				Log(m_unique_name+" ROOT array type comparison",v_debug,m_verbose);
				// First thing is we need to get the size
				size_t size1, size2;
				TString cmd;
				cmd  = TString::Format("TObjArray* tc1 = (TObjArray*)%p; ", (void*)branch1_add);
				cmd += TString::Format("TObjArray* tc2 = (TObjArray*)%p; ", (void*)branch2_add);
				cmd += TString::Format("size_t* s1 = (size_t*)%p; ", (void*)&size1);
				cmd += TString::Format("size_t* s2 = (size_t*)%p; ", (void*)&size2);
				cmd += TString("*s1 = tc1->GetEntriesFast(); *s2 = tc2->GetEntriesFast();");
				gInterpreter->ProcessLine(cmd);
				if(size1!=size2){
					printf("%lu != %lu\n",size1,size2);
					Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name
						+" TClonesArray sizes are different!",v_error,m_verbose);
				}
				// compare as many as we have
				size_t minsize = std::min(size1,size2);
				bool is_equal;
				for(int i=0; i<minsize; ++i){
					// get the next element address
					void* add1; void* add2;
					cmd  = TString::Format("void** v1 = (void**)%p; ", (void*)&add1);
					cmd += TString::Format("void** v2 = (void**)%p; ", (void*)&add2);
					cmd += TString::Format("*v1 = (void*)tc1->At(%d); ", i);
					cmd += TString::Format("*v2 = (void*)tc2->At(%d);", i);
					gInterpreter->ProcessLine(cmd);
					if(add1==nullptr || add2==nullptr){
						Log(m_unique_name+" Error! Entry "+entry_string+" item "+branch1.name
							+" returned nullptr for entry, despite being in range??",v_error,m_verbose);
						continue; // uhhh ?? this shouldn't happen. continue, i guess?
					}
					// pass these addresses to the data_instance describing element types
					branch1.contained_type->address = add1;
					branch2.contained_type->address = add2;
					// we need to recursively update all member addresses too,
					// since they currently only hold offsets to the parent object
					UpdateAddresses(*branch1.contained_type);
					UpdateAddresses(*branch2.contained_type);
					
					// if it's a TObjArray each element may be of a different type
					if(stype=="TObjArray"){
						// grab the types of this particular pair of elements.
						std::string type1; std::string type2;
						cmd  = TString::Format("TObject* oc1 = (TObject*)%p; ", (void*)add1);
						cmd += TString::Format("TObject* oc2 = (TObject*)%p; ", (void*)add2);
						cmd += TString::Format("string* t1 = (string*)%p; ", (void*)&type1);
						cmd += TString::Format("string* t2 = (string*)%p; ", (void*)&type2);
						cmd += TString::Format("*t1 = oc1->GetClass()->GetName(); ");
						cmd += TString::Format("*t2 = oc2->GetClass()->GetName();");
						gInterpreter->ProcessLine(cmd);
						if(type1!=type2){
							printf("%s != %s\n",type1.c_str(), type2.c_str()); // XXX m_verbose
							Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name
								+" found objects of different type in TObjArray element "+toString(i),
								v_error,m_verbose);
							are_equal=false;
							continue;
						} else {
							// pass this info to the data_instance
							branch1.contained_type->type_as_string = type1;
							branch2.contained_type->type_as_string = type2;
							// since we only just learnt about this type we need to parse it.
							bool parse_ok1 = ObjectToDataInstance(*branch1.contained_type);
							bool parse_ok2 = ObjectToDataInstance(*branch2.contained_type);
							if(!parse_ok1 || !parse_ok2){
								Log(m_unique_name+" Error! Failed to parse TObject of type "
									+type1+" during TObjArray handling!",v_error,m_verbose);
								are_equal=false;
								continue;
							}
						}
					}
					
					// compare these elements
					is_equal = CompareBranchMembers(*branch1.contained_type, *branch2.contained_type, less);
					if(not is_equal){
						// mismatching element!
						Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name
							+" index "
							+toString(i),v_error,m_verbose);
						are_equal=false;  // flag that the entry is not entirely identical
					}
				} // move to next element
				return are_equal;
			}
			
			// ok so it's a class. loop over the members and compare those.
			Log(m_unique_name+" comparing "+toString(branch1.members.size())+" members of "+stype,v_debug,m_verbose);
			for(std::pair<const std::string, data_instance>& nextmember : branch1.members){
				Log(m_unique_name+" comparing member "+nextmember.first,v_debug,m_verbose);
				if(branch2.members.count(nextmember.first)==0){
					Log(m_unique_name+" Error! CompareBranchMembers called with two data instances "
						+"containing classes with different members! instance "
						+branch1.name+" has member "+nextmember.first
						+" which is not present in second file data_instance!",v_error,m_verbose);
					continue;
				}
				data_instance& next_branch2_member = branch2.members.at(nextmember.first);
				Log(m_unique_name+" invoking member comparison",v_debug,m_verbose);
				are_equal = CompareBranchMembers(nextmember.second, next_branch2_member, less);
				Log(m_unique_name+" member "+nextmember.first+" comparison was "+toString(are_equal),v_debug,m_verbose);
				if(not are_equal){
					// mismatching element!
					Log(m_unique_name+" Mismatch! Entry "+entry_string+" item "+branch1.name,
						v_error,m_verbose);
				}
			} // move to next member
			// should have printed any info about discrepancies as we went
			return are_equal;
		};
		default: {
			Log(m_unique_name+" unknown instance_type "+toString(branch1.instance_type)
				+"for item "+branch1.name,v_error,m_verbose);
			return false;
		};
	} // end switch statement
	// shouldn't get to here.
	Log(m_unique_name+"Uncaught exit in CompareBranchMembers for instance type "
		+toString(branch1.instance_type),v_error,m_verbose);
	return false;
}

void CompareRootFiles::UpdateAddresses(data_instance& thedata){
	// if this data_instance is a class, recursively update all its member addresses
	// based on a newly set address for this class instance
	for(auto&& amember : thedata.members){
		void* memberadd = (void*)((char*)thedata.address + amember.second.offset); // calculate absolute address
		amember.second.address = memberadd;  // use it to set the member's address
		UpdateAddresses(amember.second);
	}
}

std::string CompareRootFiles::GetVarName(std::string type_as_string){
	if(varnames.count(type_as_string)==0){
		varnames.emplace(type_as_string,std::string("var_")+std::to_string(varnames.size()));
	}
	return varnames.at(type_as_string);
}

TMethodCall* CompareRootFiles::GetSizeCaller(data_instance& thedata){
	std::string type_as_string = thedata.type_as_string;
	if(sizecallers.count(type_as_string)==0){
		TClass* cl = TClass::GetClass(type_as_string.c_str());
		TMethodCall* m = new TMethodCall(cl,"size","");
		if(not m->IsValid()){
			// we don't seem to know about this class' size() method.
			if(loaded_libraries.count(type_as_string)==0){
				// Try to build a dictionary.
				Log(m_unique_name+" invalid TMethodCall for getting size of type "+type_as_string
					+", trying to build and load library",v_warning,m_verbose);
				bool loaded_dict = LoadDictionary(thedata);
				// see if we succeeded
				if(not loaded_dict){
					Log(m_unique_name+" failed to build dictionary to access size method of class "
						+type_as_string,v_error,m_verbose);
					sizecallers.emplace(type_as_string,(TMethodCall*)nullptr);
					return nullptr;
				}
				// succeeded in loading a dictionary, now try again to get size method
				gInterpreter->SetClassInfo(cl,true);   // this reloads the TClass
				delete m;   // presumably this is now invalidated
				m = new TMethodCall(cl,"size","");
				if(not m->IsValid()){
					Log(m_unique_name+" invalid TMethodCall for getting size of type "+type_as_string
						+", despite attempt to load library",v_error,m_verbose);
					delete m;
					sizecallers.emplace(type_as_string,(TMethodCall*)nullptr);
					return nullptr;
				} else {
					// nice! register our now functional method
					sizecallers.emplace(type_as_string,m);
				}
			} else {
				// else dictionary was already generated
				Log(m_unique_name+" invalid TMethodCall for getting size of type "+type_as_string
					+", despite previous attempt to load library",v_error,m_verbose);
				sizecallers.emplace(type_as_string,(TMethodCall*)nullptr);
			}
		} else {
			// no dictionary generation needed, register the method
			sizecallers.emplace(type_as_string,m);
		}
	} // else we've seen this class before
	
	// verbose printout....too verbose?
	if(sizecallers.at(type_as_string)==nullptr){
		Log(m_unique_name+" invalid TMethodCall for getting size of type "+type_as_string
		    +", despite previous attempt to load library",v_debug,m_verbose);
	}
	
	return sizecallers.at(type_as_string);
}

TMethodCall* CompareRootFiles::GetAtCaller(data_instance& thedata){
	std::string type_as_string = thedata.type_as_string;
	if(atcallers.count(type_as_string)==0){
		TClass* cl = TClass::GetClass(type_as_string.c_str());
		// I do not know why this two step process is needed to make it recognise it needs an int argument
		TMethodCall* mc = new TMethodCall(cl,"at","");
		mc->InitWithPrototype(cl,"at","int");
		if(not mc->IsValid()){
			if(loaded_libraries.count(type_as_string)==0){
				Log(m_unique_name+" Error building TMethodCall for "+type_as_string+".at()"
					+", trying to build and load dictionary",v_warning,m_verbose);
				bool loaded_dict = LoadDictionary(thedata);
				// see if we succeeded
				if(not loaded_dict){
					Log(m_unique_name+" failed to build dictionary to access TMethodCall for "
						+type_as_string+".at()",v_error,m_verbose);
					atcallers.emplace(type_as_string,(TMethodCall*)nullptr);
					return nullptr;
				}
				// succeeded in loading a dictionary, now try again to get size method
				gInterpreter->SetClassInfo(cl,true);   // this reloads the TClass
				delete mc;   // presumably this is now invalidated
				mc = new TMethodCall(cl,"at","");
				mc->InitWithPrototype(cl,"at","int");
				if(not mc->IsValid()){
					Log(m_unique_name+" invalid TMethodCall for "+type_as_string+".at(), "
						+"despite attempt to load library",v_error,m_verbose);
					delete mc;
					atcallers.emplace(type_as_string,(TMethodCall*)nullptr);
					return nullptr;
				} else {
					// nice! register our now functional method
					atcallers.emplace(type_as_string,mc);
				}
			} else {
				Log(m_unique_name+" invalid TMethodCall for "+type_as_string+".at(), "
				    +"despite previous attempt to load library",v_error,m_verbose);
				atcallers.emplace(type_as_string,(TMethodCall*)nullptr);
				return nullptr;
			}
		} else {
			// no dictionary generation needed, register the method
			atcallers.emplace(type_as_string,mc);
		}
	} // else we've seen this class before
	
	// verbose warning... too verbose?
	if(atcallers.at(type_as_string)==nullptr){
		Log(m_unique_name+" invalid TMethodCall for "+type_as_string+".at(), "
		    +"despite previous attempt to load library",v_debug,m_verbose);
	}
	return atcallers.at(type_as_string);
}

bool CompareRootFiles::Finalise(){
	
	return true;
}

bool CompareRootFiles::LoadConfig(std::string configfile){
	
	Log(m_unique_name+" reading configuration file "+configfile,v_debug,m_verbose);
	// read the config file
	std::ifstream fin (configfile.c_str());
	if(not fin.is_open()){
		Log(m_unique_name+" failed to read configuration file "+configfile,v_error,m_verbose);
		return false;
	}
	
	std::string Line;
	bool settingIncludePaths=false;
	bool settingSourcePathsList=false;
	
	// scan over lines in the config file
	while (getline(fin, Line)){
		Log(m_unique_name+" parsing config line \""+Line+"\"",v_debug,m_verbose);
		// skip empty lines
		if (Line.empty()) continue;
		std::string LineCopy = Line; // make a copy so we can print it in case of parsing error
		// trim preceding whitespace
		Line.erase(0,Line.find_first_not_of(" \t\015"));
		// skip comment lines
		if(Line[0] == '#') continue;
		// trim line end comments (everything after and including a '#' character)
		if(Line.find('#')!=std::string::npos) Line.erase(Line.find_first_of('#'),std::string::npos);
		// trim trailing whitespace
		if(Line.find_last_not_of(" \t\n\015\014\013")!=std::string::npos)
			Line.erase(Line.find_last_not_of(" \t\n\015\014\013")+1,std::string::npos);
		
		// split apart the key and value
		size_t splitpos = Line.find_first_of(" \t\n\015\014\013");
		std::string thekey   = Line.substr(0,splitpos);
		std::string thevalue = "";
		if(splitpos!=std::string::npos){
			thevalue = Line.substr(splitpos+1,std::string::npos);
		}
		bool push_variable=true;  // whether to record this in m_variables (skip paths and flag lines)
		
		// first check if we're parsing a list of paths for building dictionaries
		if (thekey=="StartIncludeList"){
			settingIncludePaths = true;
			push_variable=false;
		}
		else if(thekey=="EndIncludeList"){
			settingIncludePaths = false;
			push_variable=false;
		}
		else if(settingIncludePaths){
			include_paths.push_back(Line);
			push_variable=false;
		}
		else if (thekey=="StartSourcePathsList"){
			settingSourcePathsList = true;
			push_variable=false;
		}
		else if(thekey=="EndSourcePathsList"){
			settingSourcePathsList = false;
			push_variable=false;
		}
		else if(settingSourcePathsList){
			source_paths.push_back(Line);
			push_variable=false;
		}
		else if(thekey=="verbosity") m_verbose = stoi(thevalue);
		else if(thekey=="filename_1") filename_1 = thevalue;
		else if(thekey=="filename_2") filename_2 = thevalue;
		else if(thekey=="max_entries") max_entries = stoi(thevalue);
		else if(thekey=="ftolerance") ftolerance = stof(thevalue);
		else if(thekey=="index_name") index_name = thevalue;
		else {
			Log(m_unique_name+" unrecognised option in config file line: \""+LineCopy,v_error,m_verbose);
		}
		if(push_variable){ m_variables.Set(thekey,thevalue); }
	}
	// done parsing, close config file
	fin.close();
	
	// append a set of standard include and sourcefile paths
	char* empty = const_cast<char*>("");
	char* SKOFL_ROOT = getenv("SKOFL_ROOT");
	if(SKOFL_ROOT==nullptr) SKOFL_ROOT = empty;
	char* ATMPD_ROOT = getenv("ATMPD_ROOT");
	if(ATMPD_ROOT==nullptr) ATMPD_ROOT = empty;
	
	// XXX why is there both $SKOFL_ROOT/inc and $SKOFL_ROOT/include?
	// the contents are not the same, and there are some files in both that are not in the other....
	// same with ATMPD_ROOT...
	include_paths.push_back("./include");
	include_paths.push_back(std::string(SKOFL_ROOT) + "/inc");
	include_paths.push_back(std::string(SKOFL_ROOT) + "/inc/lowe");
	include_paths.push_back(std::string(SKOFL_ROOT) + "/include");
	include_paths.push_back(std::string(SKOFL_ROOT) + "/include/skonl"); 
	include_paths.push_back(std::string(SKOFL_ROOT) + "/include/lowe"); 
	include_paths.push_back(std::string(ATMPD_ROOT) + "/inc");
	include_paths.push_back(std::string(ATMPD_ROOT) + "/include");
	// TODO we could try scanning the Dependencies folder for any folders named "inc" or "include"
	
	// now sourcefiles paths
	source_paths.push_back("./DataModel");
	source_paths.push_back(std::string(SKOFL_ROOT) + "/src");
//	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/skroot");
//	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/softtrg");
//	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/ConnectionTableReader");
//	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/skrd");
////	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/iolib");
////	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/astro");
////	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/library");
////	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/monlib");
////	source_paths.push_back(std::string(SKOFL_ROOT) + "/src/mufitpe");
	source_paths.push_back(std::string(SKOFL_ROOT) + "/src");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/analysis");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/calib");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/atmflux");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/neut");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/programs");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/recon");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/reduc");
//	source_paths.push_back(std::string(ATMPD_ROOT) + "/src/usmc");
	
	return true;
}
//////////////////////////////////////////////////////////////////
//             Compiling Dictionaries on Request                //
//////////////////////////////////////////////////////////////////

bool CompareRootFiles::LoadDictionary(data_instance& thedata){
	loaded_libraries.emplace(thedata.type_as_string);
	// there are 6 steps to making CINT recognise a new type.
	
	// STEP 1. get the header for this class, along with a list of all
	//         class member types and their respective headers as well
	std::vector<std::pair<std::string,std::string>> headerlist;
	GetListOfHeaders(thedata, headerlist);
	if(headerlist.size()==0){
		Log(m_unique_name+" Error! Found no headers when building dictionary for class "+thedata.type_as_string,
		    v_error,m_verbose);
		return false;
	}
	
	// STEP 2. use these to create a LinkDef file.
	std::string linkdef_filename = BuildLinkDef(thedata, headerlist);
	if(linkdef_filename==""){
		Log(m_unique_name+" Error! Failed to build linkdef file for class "+thedata.type_as_string,v_error,m_verbose);
		return false;
	}
	
	// STEP 3. use the linkdef file to build the dictionary sourcefile
	std::string dictionary_filename = BuildDictionary(thedata.type_as_string, headerlist);
	if(dictionary_filename==""){
		Log(m_unique_name+" Error! failed to build dictionary for type "+thedata.type_as_string,v_error,m_verbose);
		return false;
	}
	
	// STEP 4. get the implementation files for these classes.
	std::vector<std::string> implementationlist = GetListOfImplementationFiles(headerlist);
	if(implementationlist.size()==0){
		Log(m_unique_name+" Warning! No implementation files found when building dictionary file class "
		    +thedata.type_as_string,v_warning,m_verbose);
	}
	
	// STEP 5. use these to compile the dictionary and sourcefiles into a library
	std::string library_file = CompileDictionary(thedata.type_as_string, dictionary_filename, implementationlist);
	if(library_file==""){
		Log(m_unique_name+" Error! Failed to build dictionary library for class "+thedata.type_as_string,
		    v_error,m_verbose);
		return false;
	}
	
	// STEP 6. load the dictionary.
	std::string cmd = "gSystem->Load(\""+library_file+"\");";
	int ret = gInterpreter->ProcessLine(cmd.c_str());
	if(ret!=0){
		Log(m_unique_name+" Error! Failed to load dictionary file for class "+thedata.type_as_string,v_error,m_verbose);
		return false;
	}
	
	return true;
}

/*
// TODO add support for containers that take two arguments e.g. maps.
{
	std::string type_as_string = thedata.type_as_string;
	// if it's a container, strip off the contained type
	if(type_as_string.find('<')!=std::string::npos){
		type_as_string = type_as_string.substr(0,type_as_string.find('<'));
		// you'd think we could use: 
		// int kind = TClassEdit::STLKind(type_as_string);
		// int nargs = TClassEdit::STLArgs(kind);
		// but this returns all kinds of garbage... we'll just have to make our own map.
		if(containers_of_two_types.count(type_as_string){
			// 2 templated types.
			// e.g. map, unordered map, multimap... 
		} else {
			// 1 templated type.
			// e.g. vector, array, set, multiset, forward_list, deque...
		}
	}
}
*/

bool CompareRootFiles::GetListOfHeaders(data_instance &thedata, std::vector<std::pair<std::string,std::string>> &headerlist){
	// note: handle members FIRST
	// this is because, when a class contains a member that is a container of a class,
	// the Linkdef file we build must contain embedded types BEFORE those that contain them.
	// e.g. if EventTrueCaptures contains a vector<TrueCapture> member, then the linkdef
	// must have `#pragma link vector<TrueCapture>` BEFORE `#pragma link EventTrueCaptures`
	for(std::pair<const std::string, data_instance>& amember : thedata.members){
		GetListOfHeaders(amember.second, headerlist);
	}
	
	std::string stype = thedata.type_as_string;
	while(IsStlContainer(stype)){
		// no point trying to get the headers from this, won't be found.
		// extract the contained type manually
		stype = TClassEdit::ShortType(thedata.type_as_string.c_str(),8);
	}
	TClass* cl = TClass::GetClass(stype.c_str());;
	if(cl){
		std::string headername = cl->GetDeclFileName();
		// strip off any preceding path stubs
		if(headername.find('/')!=std::string::npos){
			headername=headername.substr(headername.find_last_of('/')+1,std::string::npos);
		}
		// include the class name even if no header.
		// this will ensure e.g. #pragma link lines for vector<MyClass>
		// which has no header but is required!
		headerlist.emplace_back(cl->GetName(),headername);
	}
	// add the container specialization so it'll get added to the pragma lines
	if(stype!=thedata.type_as_string){
		headerlist.emplace_back(thedata.type_as_string,"");
	}
	return true;
}

std::string CompareRootFiles::BuildLinkDef(data_instance& thedata, std::vector<std::pair<std::string,std::string>> &headerlist){
	// to build a dictionary requires a Linkdef file.
	std::string linkdef_filename = thedata.type_as_string+"_LinkDef.h";
	std::replace(linkdef_filename.begin(),linkdef_filename.end(),'<','_');
	std::replace(linkdef_filename.begin(),linkdef_filename.end(),'>','_');
	std::ofstream linkdef;
	linkdef.open(linkdef_filename.c_str(), std::ios::out);
	if(not linkdef.is_open()) return "";
	
	// include all the header files
	for(auto&& aclass : headerlist){
		// if we have a header for this class
		if(aclass.second!=""){
			linkdef << "#include \""<<aclass.second<<"\"\n";
		}
	}
	linkdef << "#ifdef __CINT__\n";
	// now we need to include a pragma line for every class that ROOT needs to know about.
	for(auto&& aclass : headerlist){
		// should we use 'class+protected' instead of class+?
		// see https://root.cern.ch/root/htmldoc/guides/users-guide/AddingaClass.html which says:
		/* With '#pragma link c++ class MyClass' the Dictionary of all public members of MyClass
		will be generated. If the 'class+protected' flag is used, the dictionary for protected members
		will also be generated. This 'class+protected' flag will help you only for plain protected
		member access, but not for virtual function resolution. */
		linkdef << "#pragma link C++ class+protected "<<aclass.first<<"+;\n";
		// if the type is a container, we also need to specifically request iterators,
		// and for some reason the != operator of the iterator!!
		if(IsStlContainer(aclass.first)){
			linkdef << "#pragma link C++ class "<<aclass.first<<"::iterator;\n";  // no trailing +!!!
			linkdef << "#pragma link C++ function operator != ( "
			        <<aclass.first<<"::iterator, "<<aclass.first<<"::iterator );\n";
			// other iterator types
			linkdef << "#pragma link C++ class "<<aclass.first<<"::const_iterator;\n";  // no trailing +!!!
			linkdef << "#pragma link C++ function operator != ( "
			        <<aclass.first<<"::const_iterator, "<<aclass.first<<"::const_iterator );\n";
			// XXX does everything that provides an interator provide a reverse iterator?
			linkdef << "#pragma link C++ class "<<aclass.first<<"::reverse_iterator;\n";  // no trailing +!!!
			linkdef << "#pragma link C++ function operator != ( "
			        <<aclass.first<<"::reverse_iterator, "<<aclass.first<<"::reverse_iterator );\n";
		}
	}
	// i guess turn on nested classes and typedefs too... Not sure how this works.
	linkdef << "#pragma link C++ nestedclass;\n";
	linkdef << "#pragma link C++ nestedtypedef;\n";
	linkdef << "#endif\n";
	linkdef.close();
	
	headerlist.emplace_back("",linkdef_filename);
	
	return linkdef_filename;
}

std::vector<std::string> CompareRootFiles::GetListOfImplementationFiles(std::vector<std::pair<std::string,std::string>> &headerlist){
	std::vector<std::string> implfilelist;
	// it would be lovely if we could similarly use TClass::GetImplFileName()
	// but this only seems to be populated for ROOT classes. So we'll do a manual search.
	
	// loop over our classes
	for(auto&& aclass : headerlist){
		if(aclass.first=="" || aclass.second=="") continue;
		// this is trickier as we can't get them from the TClass,
		// so we'll look in some standard locations, and anywhere else
		// the user tells us within the config options for this tool.
		// we'll also assume the implementation filename is the same
		// as the headerfile name, with an extension of the form '[cC]*'
		std::string filename = aclass.second;
		filename = filename.substr(0,filename.find_last_of('.'));
		int matches=0;
		for(auto&& searchpath : source_paths){
			std::string cmd = std::string("find ")+searchpath
			    +" -regextype egrep -regex '.*?/?" + filename + "\\.[cC].*'";
			// try to find the file
			std::string return_string = getOutputFromFunctionCall(safeSystemCall, cmd);
			//std::string return_string = GetStdoutFromCommand(cmd);
			if(return_string!=""){
				// check we don't have more than one...there's probably a more efficient way to do this...
				std::stringstream ssl;
				ssl << return_string;
				std::string nextfilestring;
				while(getline(ssl,nextfilestring)){
					implfilelist.push_back(nextfilestring); // i mean we could just add all of them....
					++matches;
				}
				if(matches>1){
					Log(m_unique_name+" Error! Found more than one possible implementation"
					     +" file candidate for header "+aclass.second,v_error,m_verbose);
					auto it=implfilelist.rbegin();
					for(int i=0; i<matches; ++i){
						Log(m_unique_name+" Candidate: "+(*it),v_debug,m_verbose);
						++it;
					}
				}
				break; // presumably not more than one directory will have a match...
			}
			// else no match, check next path
		}
		if(matches==0){
			Log(m_unique_name+" Error! Could not find any implementation file matching header "
			   +aclass.second,v_error,m_verbose);
		}
	}
	return implfilelist;
}

std::string CompareRootFiles::BuildDictionary(std::string type_as_string, std::vector<std::pair<std::string,std::string>> &headerlist){
	// convert vector of headers into a list
	std::string headerstring;
	for(auto&& aheader : headerlist){
		if(aheader.second!=""){ headerstring.append(aheader.second + " "); }
	}
	// convert the list of header file search paths to a string too
	std::string includestring;
	for(auto&& apath : include_paths){
		includestring.append(std::string(" -I") + apath);
	}
	
	// combine things into the command to build the dictionary
	std::string dictfilename = type_as_string+"Dict.cxx";
	std::replace(dictfilename.begin(),dictfilename.end(),'<','_');
	std::replace(dictfilename.begin(),dictfilename.end(),'>','_');
	std::string cmd = "rootcint -f "+dictfilename+" -c -p -fPIC "+includestring
	                 +" `root-config --cflags` "+headerstring;
	// e.g. rootcint -f EventTrueCapturesDict.cxx -c -p -fPIC -I/HOME/ntag/NTag_ToolFramework/include
	//                 `root-config --cflags` EventTrueCaptures.h EventTrueCaptures_LinkDef.h
	Log(m_unique_name+" generating dictionary source file with command\n"+cmd,v_debug,m_verbose);
	int ret = safeSystemCall(cmd);
	if(ret==0) return dictfilename;  // 0 for success FIXME WEXITSTATUS isn't correct in Algorithms.cpp...
	// backup method: check we made the file
	//std::cout<<"safeSystemCall returned "<<ret<<std::endl;
	std::string detected_type="?";
	bool exists = CheckPath(dictfilename, detected_type);
	//std::cout<<"checkpath says exists="<<exists<<" and type is "<<detected_type<<std::endl;
	if(exists && detected_type=="f") return dictfilename;
	return "";
}


std::string CompareRootFiles::CompileDictionary(std::string type_as_string, std::string dictionaryfile, std::vector<std::string> implementationlist){
	// convert the list of header file search paths to a string too
	std::string includestring;
	for(auto&& apath : include_paths){
		includestring.append(std::string(" -I") + apath);
	}
	
	// build a list of sourcefiles
	std::string sourcesstring;
	for(auto&& afile : implementationlist){
		sourcesstring.append(afile+" ");
	}
	
	// name of the library
	std::string library_file = "lib" + type_as_string + "Dict.so";
	std::replace(library_file.begin(),library_file.end(),'<','_');
	std::replace(library_file.begin(),library_file.end(),'>','_');
	
	// put it all together
	std::string cmd = "g++ -shared -fPIC -std=c++11 `root-config --cflags --libs` "
	                + dictionaryfile + " " + sourcesstring + includestring + " -o "+library_file;
	Log(m_unique_name+" compiling dictionary with command:\n"+cmd,v_debug,m_verbose);
	
	// compile!
	int ret = safeSystemCall(cmd);
	if(ret==0) return library_file; // FIXME this is returning 1 despite success
	// backup method: check we made the file
	//std::cout<<"safeSystemCall returned "<<ret<<std::endl;
	std::string detected_type="?";
	bool exists = CheckPath(library_file, detected_type);
	//std::cout<<"checkpath says exists="<<exists<<" and type is "<<detected_type<<std::endl;
	if(exists && detected_type=="f") return library_file;
	return "";
}

// when trying to get the size of datatypes, we need something that returns the size of the datatype
// based on a string of the name. For classes we can construct a TClass from the class name and use that.
// For primitives we can use TDataType, which may be constructed from the type name, e.g. TDataType("int").
// The trouble is various methods (such as TBranch::GetExpectedType) return either a ROOTish name ("Int_t")
// or an EDataType, for which we can use the static method TDataType::GetTypeName(EDataType) to again get
// the ROOTish name.
// But, for some stupid reason, TDataType does not recognise ROOTish names in it's named constructor
// - TDataType("Int_t") is not valid! So we need to Un-ROOTish the names to construct the TDataType.
const char* CompareRootFiles::CppName(const char* type){
	if(strcmp(type,"Char_t")==0) return "char";
	if(strcmp(type,"Short_t")==0) return "short";
	if(strcmp(type,"Int_t")==0) return "int";
	if(strcmp(type,"Long_t")==0) return "long";
	if(strcmp(type,"Float_t")==0) return "float";
	if(strcmp(type,"Int_t")==0) return "int";
	if(strcmp(type,"char*")==0) return "char*";
	if(strcmp(type,"Double_t")==0) return "double";
	if(strcmp(type,"Double32_t")==0) return "double";
	if(strcmp(type,"UChar_t")==0) return "unsigned char";
	if(strcmp(type,"UShort_t")==0) return "unsigned short";
	if(strcmp(type,"UInt_t")==0) return "unsigned int";
	if(strcmp(type,"ULong_t")==0) return "unsigned long";
	if(strcmp(type,"UInt_t")==0) return "unsigned int";
	if(strcmp(type,"Long64_t")==0) return "long long";
	if(strcmp(type,"ULong64_t")==0) return "unsigned long long";
	if(strcmp(type,"Bool_t")==0) return "bool";
	if(strcmp(type,"Float16_t")==0) return "float";
	return type;
}

const char* CompareRootFiles::CppName(std::string type){
	return CppName(type.c_str());
}

bool CompareRootFiles::GetNextMatchingEntries(std::pair<const std::string, shared_tree>& apair){
	// scan until we find an entry in each tree where the index variable matches
	
	shared_tree& atree = apair.second;
	
	// advance both entries by one
	long entry_number_1 = entry_numbers_1.at(apair.first)+1;
	long entry_number_2 = entry_numbers_2.at(apair.first)+1;
	bool load_next_tree1_entry = true;
	bool load_next_tree2_entry = true;
	
	while(true){
		get_ok = 1;
		// check we have entries available
		if(load_next_tree1_entry) get_ok &= ( atree.file1_tree->LoadTree(entry_number_1)>=0 );
		if(load_next_tree2_entry) get_ok &= ( atree.file2_tree->LoadTree(entry_number_2)>=0 );
		if(get_ok){
			
			// get the TTree entries
			if(load_next_tree1_entry) atree.file1_tree->GetEntry(entry_number_1);
			if(load_next_tree2_entry) atree.file2_tree->GetEntry(entry_number_2);
			entry_string = toString(entry_number_1)+"/"+toString(entry_number_2);
			
			// TODO FIXME right now we just have one 'index variable', but we'll need one
			// for every TTree being compared (or at least to map those required to
			// their corresponding TTree)
			if(index_name==""){
				// no index variable - just do 1:1 entry comparison.
				break;
			} else if(file1_index==nullptr || file2_index==nullptr){
				// unknown index variable
				Log(m_unique_name+" Index variable "+index_name+" not found",v_error,m_verbose);
				return false;
			}
			
			// compare the indexes to see if these represent a comparable event
			int less = -1;
			get_ok = CompareBranchMembers(*file1_index, *file2_index, &less);
			Log(m_unique_name+" Index variable comparison returned "+toString(get_ok)
			    +", with less value "+toString(less),v_debug,m_verbose);
			
			if(get_ok){
				// same index - break and compare these two entries
				break;
			} else if(less==1){
				// index_1 < index_2 -- get next entry from TTree in file 1
				load_next_tree1_entry = true;
				++entry_number_1;
				load_next_tree2_entry = false;
			} else if(get_ok==0){
				// index_1 > index_2 -- get next entry from TTree in file 2
				load_next_tree1_entry = false;
				load_next_tree2_entry = true;
				++entry_number_2;
			} else {
				// non-numeric data_instance type - cannot compare
				Log(m_unique_name+" Error! Index variable "+index_name+
				    " does not appear to be of numeric type",v_error,m_verbose);
				return false;
			}
			
		} else {
			// ran off end of TTree
			logmessage = m_unique_name+" Reached off end of Tree "+atree.file1_tree->GetName();
			Log(logmessage,v_warning,m_verbose);
			return false;
		}
	}
	entry_numbers_1.at(apair.first) = entry_number_1;
	entry_numbers_2.at(apair.first) = entry_number_2;
	return true;
}
