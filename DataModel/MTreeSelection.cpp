/* vim:set noexpandtab tabstop=4 wrap */
#include "MTreeCut.h"
#include "MTreeReader.h"
#include "MTreeSelection.h"
#include "Constants.h"

//#include "BoostStore.h"
#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"
#include "TParameter.h"
#include "TKey.h"

#include <limits>

MTreeSelection::MTreeSelection(){}  // required to declare them in headers

MTreeSelection::MTreeSelection(MTreeReader* treereaderin, std::string fname, std::string distrofname){
	if(treereaderin) SetTreeReader(treereaderin);  // for writing a cut file
	MakeOutputFile(fname, distrofname);
}

MTreeSelection::MTreeSelection(std::string cutFilein){
	LoadCutFile(cutFilein);
}

MTreeSelection::~MTreeSelection(){
	if(outfile){
		outfile->Close();
//		for(auto&& acut : cut_pass_entries){   // gets cleaned up byTFile::Close
//			delete acut.second;
//			acut.second=nullptr;
//		}
		delete outfile;
	}
	if(cutfile){
		cutfile->Close();
	}
	if(distros_file){
		distros_file->Close();
		delete distros_file;
	}
}

// required for SerialisableObjects
bool MTreeSelection::Serialise(BinaryStream &bs){
	if(!(bs & cut_order)) return false;
	if(!(bs & cut_tracker)) return false;
	if(!(bs & cut_pass_entries)) return false;
	if(!(bs & branch_addresses)) return false;
	if(!(bs & known_branches)) return false;
//	if(!(bs & outfile)) return false;
//	if(!(bs & cutfile)) return false;
	if(!(bs & current_entry)) return false;
	if(!(bs & did_pass_cut)) return false;
	return true;
}

bool MTreeSelection::Print(){ PrintCuts(); return true; } // required for SerialisableObjects
std::string MTreeSelection::GetVersion(){ return "0.0"; } // required for SerialisableObjects

void MTreeSelection::MakeOutputFile(std::string fname, std::string distrofname){
	if(outfile==nullptr) outfile = new TFile(fname.c_str(),"RECREATE");
	if(distros_file==nullptr && distrofname!=""){
		distros_file = new TFile(distrofname.c_str(), "RECREATE");
		distros_tree = new TTree("vals", "Distributions of variables used for selection");
	}
}

bool MTreeSelection::SetTreeReader(MTreeReader* treereaderin){
	if(treereaderin==nullptr){
		std::cerr<<"MTreeSelection::SetTreeReader called with nullptr!"<<std::endl;
		return false;
	}
	treereader=treereaderin;
	return true;
}

bool MTreeSelection::NoteCut(std::string cutname, std::string description, double low, double high){
	auto it = cut_tracker.emplace(cutname,0);
	if(!it.second){
		std::cerr<<"MTreeSelection::NoteCut - cut "<<cutname<<" already exists!"<<std::endl;
		return false;
	}
	cut_order.push_back(cutname);
	did_pass_cut.emplace(cutname,false);
	// N.B. we must have an output file before we make the TTrees in the MTreeCut.
	cut_pass_entries.emplace(cutname, new MTreeCut(outfile, cutname, description, low, high));
	return true;
}

bool MTreeSelection::AddCut(std::string cutname, std::string description, bool savedist, double low, double high){
	bool ok = NoteCut(cutname, description, low, high);
	if(not ok){
		std::cerr<<"Failed to make cut "<<cutname<<", is cut name unique?"<<std::endl;
		return false;
	}
	TTree* dtree = (savedist ? distros_tree : nullptr);
	cut_pass_entries.at(cutname)->Initialize(0, treereader, dtree);
	return true;
}

bool MTreeSelection::AddCut(std::string cutname, std::string description, bool savedist, std::string branchname, double low, double high){
	bool ok = NoteCut(cutname, description, low, high);
	if(not ok){
		std::cerr<<"Failed to make cut "<<cutname<<", is cut name unique?"<<std::endl;
		return false;
	}
	TTree* dtree = (savedist ? distros_tree : nullptr);
	if(branchname=="") std::cerr<<"empty branchname passed for cut "<<cutname<<std::endl; // XXX
	cut_pass_entries.at(cutname)->Initialize(1, branchname, FindLinkedBranches(branchname), treereader, dtree);
	return true;
}

bool MTreeSelection::AddCut(std::string cutname, std::string description, bool savedist, std::vector<std::string> branchnames, double low, double high){
	bool ok = NoteCut(cutname, description, low, high);
	if(not ok){
		std::cerr<<"Failed to make cut "<<cutname<<", is cut name unique?"<<std::endl;
		return false;
	}
	TTree* dtree = (savedist ? distros_tree : nullptr);
	if(branchnames.size()==0){
		cut_pass_entries.at(cutname)->Initialize(0, treereader, dtree);
	} else if(branchnames.size()==1){
		cut_pass_entries.at(cutname)->Initialize(1, branchnames[0], FindLinkedBranches(branchnames[0]), treereader, dtree);
	} else {
		std::vector<std::vector<std::string>> linked_branch_lists;
		for(auto&& branchname : branchnames){
			if(branchname=="") std::cerr<<"empty string in AddPassingEvent for cut "<<cutname; // XXX
			linked_branch_lists.emplace_back(FindLinkedBranches(branchname));
		}
		cut_pass_entries.at(cutname)->Initialize(2, branchnames, linked_branch_lists, treereader, dtree);
	}
	return true;
}

void MTreeSelection::PrintCuts(){
	for(int i=0; i<cut_order.size(); ++i){
		std::cout<<((i==0) ? "\n" : "")<<"cut "<<i<<": "<<cut_order.at(i)
				 <<" => "<<cut_tracker.at(cut_order.at(i))<<"\n";
	}
}

void MTreeSelection::IncrementEventCount(std::string cutname){
	// track remaining numbers of events after each cut
	if(cut_tracker.count(cutname)==0){
		cut_tracker.emplace(cutname,1);
	} else {
		++cut_tracker.at(cutname);
	}
}

bool MTreeSelection::CheckCut(std::string cutname){
	// check we know this cut
	if(cut_pass_entries.count(cutname)==0){
		std::cerr<<"MTreeSelection::AddPassingEvent called with unknown cut "<<cutname
				 <<"\nPlease call MTreeSelection::AddCut first"<<std::endl;
		return false;
	}
	return true;
}

bool MTreeSelection::ApplyCut(std::string cutname, double val){
	bool newentry = cut_pass_entries[cutname]->Apply(val);
	if(not newentry) return false;
	IncrementEventCount(cutname);
	return true;
}

bool MTreeSelection::ApplyCut(std::string cutname, double val, size_t index){
	bool newentry = cut_pass_entries[cutname]->Apply(val, index);
	if(not newentry) return false;
	IncrementEventCount(cutname);
	return true;
}

bool MTreeSelection::ApplyCut(std::string cutname, double val, std::vector<size_t> indices){
	bool newentry = cut_pass_entries[cutname]->Apply(val, indices);
	if(not newentry) return false;
	IncrementEventCount(cutname);
	return true;
}

bool MTreeSelection::AddPassingEvent(std::string cutname){
	
	// save the TTree entry number of the passing event
	bool newentry = cut_pass_entries.at(cutname)->Enter();
	
	// increment number of events passing the cut
	if(not newentry) return false;   // prevent double counting
	IncrementEventCount(cutname);
	
	return true;
}

// provided Initialize is called first we can skip most of the arguments
bool MTreeSelection::AddPassingEvent(std::string cutname, size_t index){
	
	// save the TTree entry number and the array index.
	// A simple TEventList allows one to record which TTree entries passed a cut.
	// But we have additional arguments this time, because here a 'passing event'
	// is not just a TTree entry, but more specifically *one element of an array*
	// within a single TTree entry. We need to record the TTree entry number,
	// the branch name holding the array, and the array index.
	bool newentry = cut_pass_entries.at(cutname)->Enter(index);
	
	// increment number of events passing the cut
	if(not newentry) return false;  // prevent double counting
	IncrementEventCount(cutname);
	
	return true;
}

// but what if the "event" is specified by not just one, but several array indices?
// this time we pass a vector of pairs, each pair being the branch pointer and the index number
// call this by the syntax `AddPassingEvent("mycut",{{branchptr1,index1},{branchptr2,index2},...});
bool MTreeSelection::AddPassingEvent(std::string cutname, std::vector<size_t> indices){
	
	// save the TTree entry number and the array indices.
	bool newentry = cut_pass_entries.at(cutname)->Enter(indices);
	
	// increment number of events passing the cut
	if(not newentry) return false;  // prevent double counting
	IncrementEventCount(cutname);
	
	return true;
}

//bool MTreeSelection::Write(std::string outfilename){
//	if(outstore==nullptr){
//		outstore = new BoostStore(true,BSTORE_BINARY_FORMAT);   // typechecking enabled, single-entry binary
//	}
//	for(auto&& acut : cut_pass_entries){
//		acut.Flush();  // must call before write
//	}
//	outstore->Set("cut_order",cut_order);
//	outstore->Set("cut_tracker",cut_tracker);
//	outstore->Set("cut_pass_entries",cut_pass_entries);
//	outstore->Save(outfile.c_str());
//	outstore->Close(); // imperative!!!
//	return true;
//}

bool MTreeSelection::Write(){
	TDirectory* currdir = gDirectory;  // so we can reset it
	outfile->cd();
	if(not initialized){
		// add the meta information recording cuts in this MTreeSelection,
		// their order and how many events passed each cut
		/*
		TObjArray cut_order_obj;
		cut_order_obj.SetOwner(true);
		cut_order_obj.SetName("cut_order");
		*/
		TObjArray cut_tracker_obj;
		cut_tracker_obj.SetOwner(true);
		cut_tracker_obj.SetName("cut_tracker");
		for(auto&& acutname : cut_order){
			//cut_order_obj.Add((TObject*)(new TObjString(acutname.c_str())));
			cut_tracker_obj.Add((TObject*)(new TParameter<Long64_t>(acutname.c_str(), cut_tracker.at(acutname.c_str()))));
		}
		//cut_order_obj.Write("cut_order", TObject::kSingleKey);      // need to add TObject::kSingleKey to prevent
		cut_tracker_obj.Write("cut_tracker", TObject::kSingleKey);  // each element being written separately!
		initialized=true;
	} else {
		// update the num events passing each cut
		TObjArray* cut_tracker_obj_p = (TObjArray*)outfile->Get("cut_tracker");
		for(int cut_i=0; cut_i<cut_tracker.size(); ++cut_i){
			TParameter<Long64_t>* current_cut = (TParameter<Long64_t>*)cut_tracker_obj_p->At(cut_i);
			current_cut->SetVal(cut_tracker.at(current_cut->GetName()));
		}
		cut_tracker_obj_p->Write("cut_tracker",TObject::kOverwrite&&TObject::kSingleKey);
	}
	// write out the MTreeCuts saving which entries passed each cut
	for(auto&& acut : cut_pass_entries){
		acut.second->Write();
	}
	
	// write distributions if we're keeping them
	if(distros_tree){
		// update number of entries in the tree
		Long64_t nentries=0;
		for(int i=0; i<distros_tree->GetListOfBranches()->GetEntriesFast(); ++i){
			TBranch* br=(TBranch*)distros_tree->GetListOfBranches()->At(i);
			int thisentries = br->GetEntries();
			if(thisentries>nentries) nentries=thisentries;
		}
		distros_tree->SetEntries(nentries);
		
		distros_file->cd();
		distros_tree->Write("",TObject::kOverwrite&&TObject::kSingleKey);
	}
	
	// reset ROOT directory
	currdir->cd();
	return true;
}

std::vector<std::string> MTreeSelection::FindLinkedBranches(std::string cut_branch){
	
	// linked branches are a concept only relevant for branches holding arrays
	if(!treereader->GetBranchIsArray(cut_branch)) return std::vector<std::string>{};
	
	// when we're cutting out a particular index from an array in a branch,
	// there may be other arrays in other branches that reference the same object.
	// e.g. we may have a TTree entry containing an array of muons, which may have
	// arrays muon_times[num_muons], muon_charges[num_muons], muon_goodness[num_muons]...
	// if we're selecting by the cut "muon_times[muon_i]>250", and find that muon_i
	// passes the cut, then not only is the element "muon_times[muon_i]" relevant,
	// but so is muon_charges[muon_i], and muno_goodness[muon_i] ... etc.
	// i.e., all these branches are 'linked' - they reference the same object.
	// Let's try to be clever and identify all linked branches for a cut,
	// based on the fact that all those branches are arrays of size 'num_muons'.
	// sadly, this will only work with c-style arrays in branches. :(
	std::map<std::string,std::string> branch_titles = treereader->GetBranchTitles();
	if(branch_titles.count(cut_branch)==0){
		std::cerr<<"MTreeSelection::FindLinkedBranches could not find the title for branch "
				 <<cut_branch<<"! Link branches manually!"<<std::endl;
		return std::vector<std::string>{};
	}
	
	std::string cut_branch_title = branch_titles.at(cut_branch);
	// search for the array size specifier in the title
	size_t startpos = cut_branch_title.find_first_of("[");
	size_t endpos = cut_branch_title.find_first_of("]",startpos+1);
	if((startpos==std::string::npos)||(endpos==std::string::npos)){
		std::cerr<<"MTreeSelection::FindLinkedBranches could not extract array size from branch "
				 <<cut_branch<<" with title "<<cut_branch_title<<std::endl;
		return std::vector<std::string>{};
	}
	std::string sizestring = cut_branch_title.substr(startpos+1,endpos-startpos-1);
	// check if this string is the name of another branch - i.e. variable size array
	// (we need this to be the case to sensibly link branches, really)
	if(branch_titles.count(sizestring)==0){
		// it's not. :( Damn, must be a fixed size (e.g. muon_times[50])
		// we could check... (only good for positive integers, which is enough)
		bool not_num = (sizestring.empty()) ||
					   (sizestring.find_first_not_of("0123456789") != std::string::npos);
		if(not_num){
			std::cerr<<"MTreeSelection::FindLinkedBranches got bad static array size?"
					 <<"cut branch "<<cut_branch<<" has title "<<cut_branch_title
					 <<" with extracted size "<<sizestring<<", but no such branch could"
					 <<" be found, nor does this seem to be a number???"<<std::endl;
		}
		return std::vector<std::string>{};
	}
	// scan other branch titles for the signature that they're arrays of this size
	std::vector<std::string> linked_branchnames;
	for(auto&& abranch : branch_titles){
		if(abranch.first==cut_branch) continue; // we already know about this one
		startpos = abranch.second.find_first_of("[");
		endpos = abranch.second.find_first_of("]",startpos+1);
		if((startpos==std::string::npos)||(endpos==std::string::npos)) continue; // not an array
		std::string asizestring = abranch.second.substr(startpos+1,endpos-startpos-1);
		if(asizestring==sizestring){
			linked_branchnames.push_back(abranch.first); // a match! Shared dimension with cut_branch
		}
		// i think it's safe to assume we only need to check the first dimension for multi-dimensional branches?
	}
	return linked_branchnames;
}

// ↑↑ Methods for Writing a Cut File ↑↑
////////////////////////////////////////////////////////////////////////
// =====================================================================
////////////////////////////////////////////////////////////////////////
// ↓↓ Methods for Reading A Cut File ↓↓

bool MTreeSelection::LoadCutFile(std::string cutfilename){
	TDirectory* currdir = gDirectory;  // so we can reset it
	
	// open the input file
	cutfile = TFile::Open(cutfilename.c_str(), "READ");
	if(cutfile==nullptr){
		std::cerr<<"MTreeSelection::LoadCutFile could not find "<<cutfilename<<"!"<<std::endl;
		return false;
	}
	
	// The file contains two TObjArrays recording the cuts in this MTreeSelection,
	// their order of application and how many events passed each cut
	
	// a TObjArray with key "cut_order" stores TObjStrings with the cut names in order of application
	//TObjArray* cut_order_obj=nullptr;
	
	// a TObjArray with key "cut_tracker" stores TParameter<Long64_t> with name=<cut_name>, value=# passing events
	TObjArray* cut_tracker_obj=nullptr;
	
	// for each cut we also have a TEntryList and a TTree
	// the TEntryList key is "TEntryList_<cutname>" and has title <cutname>
	std::map<std::string, TEntryList*> cut_entrylists;
	// the TTree key is <cutname> and stores meta info and subindices for passing events, if applicable
	std::map<std::string, TTree*> cut_trees;
	
	// loop over the keys in the TFile and retrieve all this stuff
	TKey *key=nullptr;
	TIter getnextkey(cutfile->GetListOfKeys());
	// loop over objects in file
	while ((key = (TKey*)getnextkey())){
		if(key==nullptr){
			std::cerr<<"null key!"<<std::endl;
			continue;
		}
		std::string keyname=key->GetName();
		//if(keyname=="cut_order"){
		//	cut_order_obj = (TObjArray*)key->ReadObj();
		//	continue;
		//}
		if(keyname=="cut_tracker"){
			cut_tracker_obj = (TObjArray*)key->ReadObj();
			continue;
		}
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if(cl->InheritsFrom("TEntryList")){
			std::string cutname = keyname.substr(keyname.find_first_of('_')+1,std::string::npos);
			cut_entrylists.emplace(cutname,(TEntryList*)key->ReadObj());
			continue;
		}
		if(cl->InheritsFrom("TTree")){
			cut_trees[keyname]=(TTree*)key->ReadObj();  // override a nullptr from empty_write if necessary
			continue;
		}
		// we may instead have a TNamed with key "<cutname>_empty_write" if no events passed the cut
		// which indicates a Write() call was made to a cut that had no entries.... redundant?
		// it does mean that we don't know the cut details (e.g. type, branch, etc.)
		// TODO FIXME - just have a method to configure the cut branch(es) and then just have a
		// AddPassingEvent method that takes a name and just Long64_t entry number and size_t ints
		// (might need to disable the variadic expansion version for this...?)
		// (or use enable_if to prevent it being used when the third parameter (first branch) is of type size_t)
		size_t pos=keyname.find("_empty_write");
		if(pos!=std::string::npos){
			cut_trees.emplace(keyname.substr(0,pos),(TTree*)nullptr);  // will not override a valid entry if it exists
			continue;
		}
		// otherwise... dunno what this key is.
		std::cout<<"MTreeSelection ignoring unknown key in input file "<<keyname<<std::endl;
	} // end loop over keys in file
	
	// ok, now need to parse the retrieved objects
	/*
	if(cut_order_obj==nullptr){
		std::cerr<<"Did not find cut_order in input file! Is this a valid MTreeSelection output?"<<std::endl;
		return false;
	} else {
		for(int cut_i=0; cut_i<cut_order_obj->GetEntries(); ++cut_i){
			TObjString* current_cut = (TObjString*)cut_order_obj->At(cut_i);
			cut_order.push_back(current_cut->GetName());
		}
	}
	*/
	if(cut_tracker_obj==nullptr){
		std::cerr<<"Did not find cut_order in input file! Is this a valid MTreeSelection output?"<<std::endl;
		return false;
	} else {
		for(int cut_i=0; cut_i<cut_tracker_obj->GetEntries(); ++cut_i){
			TParameter<Long64_t>* current_cut = (TParameter<Long64_t>*)cut_tracker_obj->At(cut_i);
			cut_order.push_back(current_cut->GetName());
			cut_tracker.emplace(current_cut->GetName(), current_cut->GetVal());
			// TODO sanity check - check that the name is aligned with cut_order (making it redundant)?
			// update: it does seem to be redundant; removed
		}
	}
	
	// build the MTreeCut objects
	for(int cut_i=0; cut_i<cut_order.size(); ++cut_i){
		std::string next_cut_name = cut_order.at(cut_i);
		TEntryList* next_elist = cut_entrylists.at(next_cut_name);
		TTree* next_tree = cut_trees.at(next_cut_name);
		cut_pass_entries.emplace(next_cut_name, new MTreeCut(next_cut_name, next_elist, next_tree));
	}
	
	// reset ROOT directory
	currdir->cd();
	return true;
}

Long64_t MTreeSelection::GetNextEntry(std::string cutname){
	if(treereader!=nullptr){
		/*
		when we are also filling the MTreeSelection, we aren't controlling the reading of the TTree.
		Instead GetNextEntry simply marks each cut as passing or failing, 
		based on whether it's 'currententry' member is the same as the MTreeReader current entry.
		(if it doesn't pass the cut, AddPassingEvent won't have been called, and it's currentEntryNumber
		will be behind the MTreeReader).
		Note we may have calls where all MTreeCuts are marked as failing - processing tools
		should be able to cope with this.
		*/
		// loop over all cuts and see which have a current entry matching this entry
		// this tells us whether this entry passed the cut or not.
		for(auto&& acut : cut_pass_entries){
			did_pass_cut[acut.first] = (acut.second->GetCurrentEntry()==current_entry);
		}
	} else if(cutname==""){
		/*
		advance to the next entry for which at least one of the active MTreeCuts has that entry in the TEntryList.
		To do this, call MTreeCut->GetEntryList->GetNextEntry() on each active cut passing the last processed event.
		so: the first time all cuts have nextentry=-1, so call getnextentry on all.
		if returned numbers are 'cut_a:0, cut_b:0, cut_c:10, cut_d:120'
		then the next active cuts are a and b, for both of which the next entry to process is entry 0.
		Load the TTree entry 0, mark cuts a and b as 'passing' and the others as 'failing'. return.
		On the next GetNextEntry, we call MTreeCut::GetNextEntry on all cuts whose status is 'passing'
		(cuts a and b), since these have now been processed. then proceed as above -
		if the returned numbers are cut_a:5, cut_b:10, then the next active entry is 5: load TTree entry 5,
		mark cut a as passing, all others as failing.
		on the subsequent GetNextEntry, advance only cut a - let's say the returned number is 10.
		Now, cuts 0, 1 and 2 are all passing. Load TTree entry 10. And so on.
		Other tools then process all indices of the passing cuts.
		*/
		
		// no specific cut: load the next entry that passes any of the cuts
		std::map<Long64_t, std::vector<std::string>> current_entry_numbers;
		for(auto&& acut : cut_pass_entries){
			Long64_t the_entry_number = acut.second->GetCurrentEntry();
			if(current_entry_numbers.count(the_entry_number)==0){
				current_entry_numbers.emplace(the_entry_number, std::vector<std::string>{acut.first});
			} else {
				current_entry_numbers.at(the_entry_number).push_back(acut.first);
			}
		}
		// using the fact that maps are inherently ordered, the first entry is the lowest entry number.
		// That'll be the one we just processed. Advance all the MTreeCuts that were holding that entry.
		std::vector<std::string> cuts_to_advance = current_entry_numbers.begin()->second;
		// remove the entry numbers for all cuts we're advancing
		current_entry_numbers.erase(current_entry_numbers.begin());
		for(std::string acut : cuts_to_advance){
			// update as we go
			Long64_t the_entry_number = cut_pass_entries[acut]->GetNextEntry();
			if(current_entry_numbers.count(the_entry_number)==0){
				current_entry_numbers.emplace(the_entry_number, std::vector<std::string>{acut});
			} else {
				current_entry_numbers.at(the_entry_number).push_back(acut);
			}
		}
		// the first entry is the lowest passing entry number, which will be the entry we'll now process
		current_entry = current_entry_numbers.begin()->first;
		// the value of the first entry is the name of all the cuts for which this entry passed.
		// reinit the list of cuts that passed
		for(auto&& acut : did_pass_cut) acut.second=false;
		// and set those that did
		for(auto&& apassedcut : current_entry_numbers.begin()->second) did_pass_cut[apassedcut]=true;
		
	} else {
		/*
		if we're interested in a specific cut, we should always call GetNextEntry on that cut.
		then we need to advance any cuts that are now behind, before checking which are passing.
		*/
		
		// load the next passing entry of the requested cut
		current_entry = cut_pass_entries[cutname]->GetNextEntry();
		if(current_entry<0) return -1; // end of TEntryList
		// see which other cuts pass this. To do this we'll first need to find any cuts
		// whose current entry is lower than the one we're now processing and advance them.
		for(auto&& acut : cut_pass_entries){
			if(acut.first==cutname) continue; // already did this one
			Long64_t next_passing_entry = acut.second->GetCurrentEntry();
			while((next_passing_entry<current_entry) && (next_passing_entry>=0)){
				next_passing_entry = acut.second->GetNextEntry();
			};
		}
		// now loop over all cuts and see which have a current entry matching this entry
		// this tells us whether this entry passed the cut or not.
		for(auto&& acut : cut_pass_entries){
			did_pass_cut[acut.first] = (acut.second->GetCurrentEntry()==current_entry);
		}
	}
	return current_entry;
}

bool MTreeSelection::GetPassesCut(std::string cutname){
	if(did_pass_cut.count(cutname)==0){
		std::cerr<<"MTreeSelection::GetPassesCut called with unknown cut "<<cutname<<std::endl;
		return false;
	}
	return did_pass_cut[cutname];
}

// TODO something to check that the indices are for the right branches, and in the right order?
// maybe accept a branch pointer and check?
// or in the case of multiple indices, re-order if necessary?
bool MTreeSelection::GetPassesCut(std::string cutname, size_t index){
	if(did_pass_cut.count(cutname)==0){
		std::cerr<<"MTreeSelection::GetPassesCut called with unknown cut "<<cutname<<std::endl;
		return false;
	}
	if(not did_pass_cut[cutname]) return false;
	// otherwise check this index
	return cut_pass_entries[cutname]->GetPassingIndexes().count(index);
}

bool MTreeSelection::GetPassesCut(std::string cutname, std::vector<size_t> indices){
	if(did_pass_cut.count(cutname)==0){
		std::cerr<<"MTreeSelection::GetPassesCut called with unknown cut "<<cutname<<std::endl;
		return false;
	}
	if(not did_pass_cut[cutname]) return false;
	// otherwise check indices
	return cut_pass_entries[cutname]->GetPassingIndices().count(indices);
}

std::set<size_t> MTreeSelection::GetPassingIndexes(std::string cutname){
	if(did_pass_cut.count(cutname)==0){
		std::cerr<<"MTreeSelection::GetPassingIndexes called with unknown cut "<<cutname<<std::endl;
		return std::set<size_t>{};
	}
	if(did_pass_cut[cutname]){
		return cut_pass_entries[cutname]->GetPassingIndexes();
	} else {
		// if it didn't pass the cut, it has no indices this entry
		return std::set<size_t>{};
	}
}

std::set<std::vector<size_t>> MTreeSelection::GetPassingIndices(std::string cutname){
	if(did_pass_cut.count(cutname)==0){
		std::cerr<<"MTreeSelection::GetPassingIndices called with unknown cut "<<cutname<<std::endl;
		return std::set<std::vector<size_t>>{};
	}
	if(did_pass_cut[cutname]){
		return cut_pass_entries[cutname]->GetPassingIndices();
	} else {
		// if it didn't pass the cut, it has no indices this entry
		return std::set<std::vector<size_t>>{};
	}
}

MTreeReader* MTreeSelection::GetTreeReader(){
	return treereader;
}

std::string MTreeSelection::GetTopCut(){
	return ((cut_order.size()==0) ? "" : cut_order.front());
}

std::string MTreeSelection::BranchAddressToName(intptr_t branchptr){
	if(branch_addresses.count(branchptr)==0){
		std::string indexbranch="";
		// we haven't yet seen this branch, or the branch address has changed
		// scan the current branch addresses to find which one fills this variable
		std::map<std::string, intptr_t> branchadds = treereader->GetBranchAddresses();
		for(auto&& abranch : branchadds){
			if(abranch.second==branchptr){
				indexbranch=abranch.first;
				break;
			}
		}
		if(indexbranch==""){
			std::cerr<<"MTreeSelection::BranchAddressToName failed to find name for branch!"<<std::endl;
			return "";
		}
		branch_addresses.emplace(branchptr, indexbranch);
		// it shouldn't happen often, but each time a branch address changes,
		// it'll orphan an entry in our branch address to name map.
		// Keep a reverse map, and use it to remove the old entry if one exists.
		if(known_branches.count(indexbranch)){
			branch_addresses.erase(known_branches.at(indexbranch));
			known_branches.at(indexbranch) = branchptr;
		} else {
			known_branches.emplace(indexbranch,branchptr);
		}
		return indexbranch;
	} else {
		return branch_addresses.at(branchptr);
	}
}

Long64_t MTreeSelection::GetEntries(std::string cutname){
	if(cutfile==nullptr){
		return cut_tracker.at(cutname);
	}
	return cut_pass_entries.at(cutname)->GetEntries();
}
