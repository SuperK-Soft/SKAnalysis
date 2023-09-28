#include "MTreeCut.h"

#include "TROOT.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TObjectTable.h"

#include <iostream>
#include <sstream> // std::stringstream

////// debug
//#include "Algorithms.h"         // for getOutputFromFunctionCall

#include "MTreeReader.h"

MTreeCut::MTreeCut(TFile* outfilein, std::string cutname, std::string description, double low, double high)
   : cut_name(cutname), cut_description(description), outfile(outfilein) {
	mode="write";
	if(low!=DOUBLE_MIN || high!=DOUBLE_MAX){
		low_thresh = low;
		high_thresh = high;
		got_low_thresh=true;
		got_high_thresh=true;
	} else {
		got_bool_req=true;
	}
}

MTreeCut::MTreeCut(std::string cutname, TEntryList* inelist, TTree* intree) : additional_indices(intree), ttree_entries(inelist) {
	mode="read";
	total_entries = ttree_entries->GetN();
	GetMetaInfo();
	SetBranchAddresses();
	GetNextEntry(); // load first entry
}

Long64_t MTreeCut::GetEntries(){
	//return ttree_entries->GetN();
	return total_entries;
}

bool MTreeCut::Serialise(BinaryStream &bs){
	if(!(bs & mode)) return false;
	if(!(bs & cut_name)) return false;
	if(!(bs & cut_description)) return false;
//	if(!(bs & ttree_entries)) return false;
//	if(!(bs & additional_indices)) return false;
	if(!(bs & current_entry)) return false;
	if(!(bs & tlist_entry)) return false;
	if(!(bs & type)) return false;
	// for type 1
	if(!(bs & additional_branchname)) return false;
	if(!(bs & linked_branch_list)) return false;
	if(!(bs & indexes_this_entry)) return false;
	// for type 2
	if(!(bs & additional_branchnames)) return false;
	if(!(bs & linked_branch_lists)) return false;
	if(!(bs & indices_this_entry)) return false;
	return true;
}

bool MTreeCut::Print(){ return true; } // dummy, required for SerialisableObjects
std::string MTreeCut::GetVersion(){ return "0.0"; }

void MTreeCut::SetBranchAddresses(){
	if(type>0){
		additional_indices->SetBranchAddress("TreeEntry",&current_entry);
	}
	if(type==1){
		indexes_this_entry_p = &indexes_this_entry;
		additional_indices->SetBranchAddress("AdditionalIndices",&indexes_this_entry_p);
	} else if(type==2){
		indices_this_entry_p = &indices_this_entry;
		additional_indices->SetBranchAddress("AdditionalIndices",&indices_this_entry_p);
	}
}

MTreeCut::~MTreeCut(){
	if(mode=="write"){
		if(additional_indices){ delete additional_indices; additional_indices=nullptr; }
		if(ttree_entries){ delete ttree_entries; ttree_entries = nullptr; }
		if(values_branch){ values_branch->ResetAddress(); values_branch=nullptr; }
		if(pass_branch){ pass_branch->ResetAddress(); pass_branch=nullptr; }
	}
}

// for type 0
void MTreeCut::Initialize(int type_in, MTreeReader* readerIn, TTree* distro_tree){
	auto currdir = gDirectory;
	outfile->cd();
	
	type=type_in;
	ttree_entries = new TEntryList(cut_name.c_str(),cut_name.c_str());
	ttree_entries->SetName(TString::Format("TEntryList_%s",cut_name.c_str()));
	
	// even though it'll have no branches, we still need a TTree for the meta info
	additional_indices = new TTree(cut_name.c_str(),cut_description.c_str());
	
	// store meta info
	TNamed* thecutname = new TNamed("cut_name", cut_name.c_str());
	additional_indices->GetUserInfo()->Add(thecutname);
	
	TNamed* thecutdescription = new TNamed("cut_description",cut_description.c_str());
	additional_indices->GetUserInfo()->Add(thecutdescription);
	
	TParameter<Int_t>* thecuttype = new TParameter<Int_t>("cut_type",type);
	additional_indices->GetUserInfo()->Add(thecuttype);
	
	theReader = readerIn;
	if(distro_tree) MakeDistroBranches(distro_tree);
	
	currdir->cd();
}

//for type 1
void MTreeCut::Initialize(int type_in, std::string indexcutbranch, std::vector<std::string> linkedbranches, MTreeReader* readerIn, TTree* distro_tree){
	auto currdir = gDirectory;
	outfile->cd();
	
	type=type_in;
	ttree_entries = new TEntryList(cut_name.c_str(),cut_name.c_str());
	ttree_entries->SetName(TString::Format("TEntryList_%s",cut_name.c_str()));
	additional_branchname = indexcutbranch;
	linked_branch_list = linkedbranches;
	
	// the TTree will store the additional indices required to specify array indices within this TTree entry
	additional_indices = new TTree(cut_name.c_str(),cut_description.c_str());
	additional_indices->Branch("TreeEntry",&current_entry);
	additional_indices->Branch("AdditionalIndices",&indexes_this_entry);
	
	// store meta info
	TNamed* thecutname = new TNamed("cut_name", cut_name.c_str());
	additional_indices->GetUserInfo()->Add(thecutname);
	
	TNamed* thecutdescription = new TNamed("cut_description",cut_description.c_str());
	additional_indices->GetUserInfo()->Add(thecutdescription);
	
	TParameter<Int_t>* thecuttype = new TParameter<Int_t>("cut_type",type);
	additional_indices->GetUserInfo()->Add(thecuttype);
	
	// As in the example above we often have several branches storing arrays that share an index;
	// the entries 'hit_charge[0]' and 'hit_time[0]' both represent the same hit.
	// So, a cut on 'hit_charge' should be reflected in which entries are valid in 'hit_time'.
	// In this method the 'hit_charge' branch might be the 'indexcutbranch', and the 'hit_time'
	// would be one of potentially many 'linked branches'.
	// Of course, there's some degeneracy here in what's the 'cut branch' and what's a 'linked branch':
	// the underlying point is "hit 0 passed the cut", so it doesn't really matter which is the
	// cut branch and which the linked branch: in either case hit_charge[0] and hit_time[0]
	// are both pieces of information about a passing entry. Even the idea that 'hit_charge'
	// was the branch that was cut on is a restrictive notion: the cut may have been complex,
	// operating on several branches and other derived variables. Point is, as long as the given
	// 'indexcutbranch' one of the branches with the right indexing, it'll work out.
	// Why store the special "cut branch" at all? .... just in case. ;)
	
	TNamed* thecutbranch = new TNamed("cut_branch",indexcutbranch.c_str());
	additional_indices->GetUserInfo()->Add(thecutbranch);
	
	TObjArray* linkedbrancharr = new TObjArray();
	linkedbrancharr->SetOwner(true);  // own its contents, handle cleanup
	linkedbrancharr->SetName("linked_branch_list");
	for(auto&& alinkbr : linkedbranches){
		linkedbrancharr->Add((TObject*)(new TObjString(alinkbr.c_str())));
	}
	additional_indices->GetUserInfo()->Add(linkedbrancharr);
	
	theReader = readerIn;
	if(distro_tree) MakeDistroBranches(distro_tree);
	
	currdir->cd();
}

// for type 2
void MTreeCut::Initialize(int type_in, std::vector<std::string> indexcutbranches, std::vector<std::vector<std::string>> linkedbranches, MTreeReader* readerIn, TTree* distro_tree){
	auto currdir = gDirectory;
	outfile->cd();
	
	type=type_in;
	ttree_entries = new TEntryList(cut_name.c_str(),cut_name.c_str());
	ttree_entries->SetName(TString::Format("TEntryList_%s",cut_name.c_str()));
	additional_branchnames = indexcutbranches;
	linked_branch_lists = linkedbranches;
	
	// TTree to store additional indices
	additional_indices = new TTree(cut_name.c_str(),cut_description.c_str());
	additional_indices->Branch("TreeEntry",&current_entry);
	additional_indices->Branch("AdditionalIndices",&indices_this_entry);
	
	// store meta info
	TNamed* thecutname = new TNamed("cut_name", cut_name.c_str());
	additional_indices->GetUserInfo()->Add(thecutname);
	
	TNamed* thecutdescription = new TNamed("cut_description",cut_description.c_str());
	additional_indices->GetUserInfo()->Add(thecutdescription);
	
	TParameter<Int_t>* thecuttype = new TParameter<Int_t>("cut_type",type);
	additional_indices->GetUserInfo()->Add(thecuttype);
	
	// setup branches to store additional indices
	TObjArray* indexcutbranchesarr = new TObjArray();
	indexcutbranchesarr->SetOwner(true);  // own its contents, handle cleanup
	indexcutbranchesarr->SetName("cut_branches");
	for(auto&& aindexbr : indexcutbranches){
		indexcutbranchesarr->Add((TObject*)(new TObjString(aindexbr.c_str())));
	}
	additional_indices->GetUserInfo()->Add(indexcutbranchesarr);
	
	TObjArray* linkarra = new TObjArray();  // will be an array of arrays
	linkarra->SetOwner(true);
	linkarra->SetName("linked_branch_lists");
	for(int cutbranchi=0; cutbranchi<linkedbranches.size(); ++cutbranchi){
		TObjArray* temparr = new TObjArray();
		temparr->SetOwner(true);
		// set name to the cut branch these are linked to
		temparr->SetName(indexcutbranches.at(cutbranchi).c_str());
		for(auto&& alink : linkedbranches.at(cutbranchi)){
			temparr->Add((TObject*)(new TObjString(alink.c_str())));
		}
		linkarra->Add((TObject*)temparr);
	}
	additional_indices->GetUserInfo()->Add(linkarra);
	
	theReader = readerIn;
	if(distro_tree) MakeDistroBranches(distro_tree);
	
	currdir->cd();
	
}

bool MTreeCut::MakeDistroBranches(TTree* distro_tree){
	values_branch=distro_tree->Branch(cut_name.c_str(), &branch_val);
	pass_branch=distro_tree->Branch((cut_name+"_pass").c_str(), &pass_val);
	save_distros = (values_branch!=nullptr && pass_branch!=nullptr);
	if(!save_distros){
		std::cerr<<"MTreeCut::MakeDistroBranches failed to make branches: "<<values_branch<<", "<<pass_branch<<std::endl;
	}
	return save_distros;
}

bool MTreeCut::Apply(double value){
	pass_val = true;
	if(got_low_thresh  && value < low_thresh)   pass_val=false;
	if(got_high_thresh && value > high_thresh)  pass_val=false;
	if(got_bool_req    && value == 0.0)       pass_val=false;
	if(save_distros){
		branch_val = value;
		values_branch->Fill();
		pass_branch->Fill();
	}
	if(!pass_val) return false;
	return Enter();
}

bool MTreeCut::Apply(double value, size_t index){
	pass_val = true;
	if(got_low_thresh  && value < low_thresh)   pass_val=false;
	if(got_high_thresh && value > high_thresh)  pass_val=false;
	if(got_bool_req    && value == 0.0)       pass_val=false;
	if(save_distros){
		branch_val = value;
		values_branch->Fill();
		pass_branch->Fill();
	}
	if(!pass_val) return false;
	return Enter(index);
}

bool MTreeCut::Apply(double value, std::vector<size_t> indices){
	pass_val = true;
	if(got_low_thresh  && value < low_thresh)   pass_val=false;
	if(got_high_thresh && value > high_thresh)  pass_val=false;
	if(got_bool_req    && value == 0.0)       pass_val=false;
	if(save_distros){
		branch_val = value;
		values_branch->Fill();
		pass_branch->Fill();
	}
	if(!pass_val) return false;
	return Enter(indices);
}

// type 0
bool MTreeCut::Enter(){
	if(type!=0){
		std::cerr<<"MTreeCut::Enter() called with no additional indices on MTreeCut "<<cut_name
		         <<" but its type is "<<type<<"!"<<std::endl;
		return false;
	}
	Long64_t entry_num = theReader->GetEntryNumber();
	TTree* t = theReader->GetTree();
	if(!t){ std::cerr<<"MTreeCut::Enter "<<cut_name<<" TREE IS NULL!"<<std::endl; return false; }
	bool newentry = ttree_entries->Enter(entry_num, t);
	//return ttree_entries->Enter(theReader->GetEntryNumber(), theReader->GetTree());
	if(newentry) ++total_entries;
	return newentry;
}

// type 1
bool MTreeCut::Enter(size_t index){
	if(type!=1){
		std::cerr<<"MTreeCut::Enter() called with a single index on MTreeCut "<<cut_name
		         <<" but its type is "<<type<<"!"<<std::endl;
		return false;
	}
	Long64_t entry_number = theReader->GetEntryNumber();
	// if starting a new TTree entry, write out all passing indices for the last entry
	if((current_entry!=entry_number)&&(indexes_this_entry.size()!=0)){
		additional_indices->Fill();
		indexes_this_entry.clear();
	}
	// add the entry_number to the TEntryList, if it isn't already
	bool newtreeentry = ttree_entries->Enter(entry_number, theReader->GetTree());
	// sanity check
	if((current_entry!=entry_number)&&(newtreeentry==false)){
		std::cerr<<"Out of order call to MTreeCut::Enter! All passing sub-indices for a given "
		         <<"TTree entry must be given in sequence!"<<std::endl;
		// we can't do: Enter(entry_number=0, index=1), Enter(entry_number=1, index=5),
		// then Enter(entry_number=0, index=2); i.e. entry number 0, then 1, then 0 again.
		// To accommodate that we'd need to go back to a previous additional_indices entry
		// and update with more indices - TTrees don't support that.
	}
	current_entry=entry_number;
	// add the passing index to the list of indices
	auto ret = indexes_this_entry.emplace(index);
	// std::set::emplace returns a std::pair of an iterator to element and a boolean of whether it's new
	// return whether this combination was a unique new combination
	return (newtreeentry || ret.second);
}

// type 2
bool MTreeCut::Enter(std::vector<size_t>& indices){
	if(type!=2){
		std::cerr<<"MTreeCut::Enter() called with a vector of indices on MTreeCut "<<cut_name
		         <<" but its type is "<<type<<"!"<<std::endl;
		return false;
	}
	Long64_t entry_number = theReader->GetEntryNumber();
	// if starting a new TTree entry, write out all passing indices for the last entry
	if((current_entry!=entry_number)&&(indices_this_entry.size()!=0)){
		additional_indices->Fill();
		indices_this_entry.clear();
	}
	// add the entry_number to the TEntryList, if it isn't already
	bool newtreeentry = ttree_entries->Enter(entry_number, theReader->GetTree());
	// sanity check
	if((current_entry!=entry_number)&&(newtreeentry==false)){
		std::cerr<<"Out of order call to MTreeCut::Enter! All passing sub-indices for a given "
		         <<"TTree entry must be given in sequence!"<<std::endl;
	}
	current_entry=entry_number;
	auto ret = indices_this_entry.emplace(indices);
	// return whether this combination was a unique new combination
	return (newtreeentry || ret.second);
}

bool MTreeCut::Flush(){
	// call before saving to write out the last entry
	if(additional_indices&&current_entry>=0){
		additional_indices->Fill();
		if(type==1) indexes_this_entry.clear();
		if(type==2) indices_this_entry.clear();
		current_entry=-1;
		return true;
	}
	return false;
}

void MTreeCut::Write(){
	auto currdir = gDirectory;
	outfile->cd();
	if(type<0){
		// in the event that we've not yet added any entries, let's just put a flag in the TFile.
		// leave the 'cut_name' TFile key free, so that we can use it as normal if we get a later Write call
		// once we have some entries.
		std::string flagstring = cut_name+"_empty_write";
		TNamed flag(flagstring.c_str(), flagstring.c_str());
		flag.Write(flagstring.c_str(), TObject::kOverwrite);
	} else {
		ttree_entries->Write("",TObject::kOverwrite);
		Flush();
		additional_indices->Write("",TObject::kOverwrite);
	}
	currdir->cd();
}


////////////////////////////////////////////////////////////////////////
// ↑↑ Methods for Writing MTreeCuts
// ================================
// ↓↓ Methods for Reading MTreeCuts
////////////////////////////////////////////////////////////////////////

bool MTreeCut::GetMetaInfo(){
	// retrieve meta info from additional_indices TTree
	// our tree has:
	// a branch TreeEntry storing the next passing entry
	// a branch AdditionalIndices storing a vector of additional indices, if applicable
	// a TList of meta info in the UserInfo
	
	TList* next_meta_info = additional_indices->GetUserInfo();
	// we need to retrieve:
	// a TNamed with name = "cut_description" and title storing a description of the cut (e.g. threshold)
	// a TParameter<Int_t> with name = "cut_type" and value of the cut type (0-2)
	// for type 1 cuts we also have:
	// a TNamed with name = "cut_branch" and title giving the name of the branch the index relates to
	// a TObjArray with name = "linked_branch_list" that stores TObjStrings with the branches
	// 'linked' to this branch - i.e. arrays that share an indexing with the cut branch
	// for type 2 cuts we also have:
	// a TObjArray with name = "cut_branches" that stores TObjStrings with the branches the indices relate to
	// a TObjArray with name "linked_branch_lists" that stores an array of TObjArrays, each of which
	// is an array of linked branches for the corresponding branch, as per type 1.
	for(int obj_i=0; obj_i<next_meta_info->GetEntries(); ++obj_i){
		std::string obj_name = next_meta_info->At(obj_i)->GetName();
		if(obj_name=="cut_name"){
			TNamed* thecutname = (TNamed*)next_meta_info->At(obj_i);
			cut_name = thecutname->GetTitle();
		} else if(obj_name=="cut_description"){
			TNamed* thecutdescription = (TNamed*)next_meta_info->At(obj_i);
			cut_description = thecutdescription->GetTitle();
		} else if(obj_name=="cut_type"){
			TParameter<Int_t>* thecuttype = (TParameter<Int_t>*)next_meta_info->At(obj_i);
			type = thecuttype->GetVal();
		} else if(obj_name=="cut_branch"){
			TNamed* this_cut_branch = (TNamed*)next_meta_info->At(obj_i);
			additional_branchname = this_cut_branch->GetTitle();
		} else if(obj_name=="linked_branch_list"){
			TObjArray* linked_branch_arr = (TObjArray*)next_meta_info->At(obj_i);
			for(int link_i=0; link_i<linked_branch_arr->GetEntries(); ++link_i){
				linked_branch_list.push_back(linked_branch_arr->At(link_i)->GetName());
			}
		} else if(obj_name=="cut_branches"){
			TObjArray* branchnames_arr = (TObjArray*)next_meta_info->At(obj_i);
			for(int br_i=0; br_i<branchnames_arr->GetEntries(); ++br_i){
				additional_branchnames.push_back(branchnames_arr->At(br_i)->GetName());
			}
		} else if(obj_name=="linked_branch_lists"){
			TObjArray* linked_list_arr = (TObjArray*)next_meta_info->At(obj_i);
			for(int branch_i=0; branch_i<linked_list_arr->GetEntries(); ++branch_i){
				std::vector<std::string> temp_vec;
				TObjArray* linked_branch_arr = (TObjArray*)linked_list_arr->At(branch_i);
				for(int link_i=0; link_i<linked_branch_arr->GetEntries(); ++link_i){
					temp_vec.push_back(linked_branch_arr->At(link_i)->GetName());
				}
				linked_branch_lists.push_back(temp_vec);
			}
		}
		// else an object in the UserInfo we don't know?
	}
	return true;
}

Long64_t MTreeCut::GetNextEntry(){
	if(mode=="read"){
		if(tlist_entry>total_entries) return -1; // end of TEntryList
		++tlist_entry;
		current_entry = ttree_entries->GetEntry(tlist_entry);
		if(type>0){
			additional_indices->GetEntry(tlist_entry);
		}
	}
	// bypass when writing
	return current_entry;
}

Long64_t MTreeCut::GetCurrentEntry(){
	return current_entry;
}

std::set<size_t> MTreeCut::GetPassingIndexes(){
	return indexes_this_entry;
}

std::set<std::vector<size_t>> MTreeCut::GetPassingIndices(){
	return indices_this_entry;
}
