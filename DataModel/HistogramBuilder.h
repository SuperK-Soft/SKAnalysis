#ifndef HISTOGRAM_BUILDER_H
#define HISTOGRAM_BUILDER_H
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <type_traits>
#include <map>
#include <iostream>
#include <cxxabi.h>

class HistogramBuilder {
	
	public:
	HistogramBuilder();
	~HistogramBuilder();
	TFile* MakeFile(std::string filename, std::string treename="tree", bool notree_in=false);
	TTree* MakeTree(std::string treename="tree");
	void SaveHists(bool dosavehists);
	void SetVerbosity(int verb);
	
	bool Save();
	bool Close();
	TFile* GetFile();
	TTree* GetTree();
	int SetTreeEntries();
	size_t GetNHists();
	std::map<std::string,TH1*> GetHists();
	TH1* GetHist(std::string hname, std::string cut, int unbinned=-1);  // user should pass unbinned == 0 or 1
	TH1* GetHist(size_t i);
	
	template<typename T>
	TH1* AddHist(std::string histname, T type);
	template<typename T>
	TH1* AddHist(std::string histname, T type1, T type2);
	template<typename T>
	TH1* AddHist(std::string histname, T type1, T type2, T type3);
	
	template<typename T>
	TBranch* AddBranch(std::string branchname, T type, typename std::enable_if<std::is_arithmetic<T>::value, bool>::type potato=true);
	template<typename T>
	TH1* AddHist(std::string histname, T type, std::array<double,3> binning);
	template<typename T>
	TH2* AddHist(std::string histname, T type, std::array<double,6> binning);
	template<typename T>
	TH3* AddHist(std::string histname, T type, std::array<double,9> binning);
	
	private:
	template<typename T>
	TH1* AddHist(std::string histname, std::vector<T> type);
	std::string AddHist(std::string histname);
	
	template<typename T>
	TH1D* NewHist(std::string name, std::string unique_name, T dummy, std::array<double,3> binning, typename std::enable_if<std::is_floating_point<T>::value, bool>::type potato= true);
	template<typename T>
	TH1I* NewHist(std::string name, std::string unique_name, T dummy, std::array<double,3> binning, typename std::enable_if<std::is_integral<T>::value, bool>::type potato= true);
	template<typename T>
	TH2D* NewHist(std::string name, std::string unique_name, T dummy, std::array<double,6> binning, typename std::enable_if<std::is_floating_point<T>::value, bool>::type potato= true);
	template<typename T>
	TH2I* NewHist(std::string name, std::string unique_name, T dummy, std::array<double,6> binning, typename std::enable_if<std::is_integral<T>::value, bool>::type potato= true);
	template<typename T>
	TH3D* NewHist(std::string name, std::string unique_name, T dummy, std::array<double,9> binning, typename std::enable_if<std::is_floating_point<T>::value, bool>::type potato= true);
	template<typename T>
	TH3I* NewHist(std::string name, std::string unique_name, T dummy, std::array<double,9> binning, typename std::enable_if<std::is_integral<T>::value, bool>::type potato= true);
	
	public:
	template <typename T>
	bool Fill(std::string name, T val);
	template <typename T>
	bool Fill(std::string name, T val1, T val2);
	template <typename T>
	bool Fill(std::string name, T val1, T val2, T val3);
	template <typename T>
	bool Fill(std::string name, std::vector<T> vals);
	private:
	template <typename T>
	bool BranchFill(std::string uniquebranchname, T val);
	
	private:
	TFile* ofile=nullptr;
	std::string ofilename="";
	bool fileistmp=false;
	bool savehists=true;
	TTree* tree=nullptr;
	bool notree=false; // whether to just make histograms, not a TTree
	int m_verbose=0;
	
	static size_t hbuildercounter;
	std::map<std::string, std::string> histnames;
	std::map<std::string, TH1*> hists;
	std::map<std::string, TBranch*> branches;
	std::map<std::string, void*> branchvals;
	
	// char, double, float, int, short
	enum class branchType : char { C, D, F, I, S };
	static const std::map<std::string, branchType> typechars;
	std::map<std::string, branchType> branchtypes;
	
	
};

//=========================================================================================//

template<typename T>
TBranch* HistogramBuilder::AddBranch(std::string branchname, T type, typename std::enable_if<std::is_arithmetic<T>::value, bool>::type potato){
	
	// sanity check
	if(notree){
		std::cerr<<"HistogramBuilder::AddBranch called but 'notree' is set!"<<std::endl;
		return nullptr;
	}
	
	// make a tree if we haven't got one
	if(tree==nullptr){
		// make a default tree
		std::cerr<<"HistogramBuilder::AddBranch making tree with name 'data'"<<std::endl;
		if(MakeTree("data")==nullptr) return nullptr;
	}
	
	// check that this branch doesn't exist already
	if(tree->GetBranch(branchname.c_str())!=nullptr){
		std::cerr<<"HistogramBuilder::AddBranch "<<branchname<<" already exists!"<<std::endl;
		return nullptr;
	}
	
	// we need to be consistent in using a variable of the same type for TTree::Branch
	// and subsequent Fill calls. SetBranchAddress will fail if you try to update it
	// to point to a new variable of different type.
	std::string thetype = abi::__cxa_demangle(typeid(T).name(),nullptr, nullptr, nullptr);
	if(typechars.count(thetype)==0){
		std::cout<<"HistogramBuilder::AddBranch type '"<<thetype<<"' is not supported"<<std::endl;
		return nullptr;
	}
	branchType typechar = typechars.at(thetype);
	branchtypes.emplace(branchname, typechar);
	T* val = new T;
	branchvals[branchname] = (void*)(val);
	branches.emplace(branchname, tree->Branch(branchname.c_str(), val));
	
	return branches.at(branchname);
	
}

template<typename T>
TH1* HistogramBuilder::AddHist(std::string histname, T type){
	return AddHist(histname, type, std::array<double,3>{});
}

template<typename T>
TH1* HistogramBuilder::AddHist(std::string histname, T type1, T type2){
	return AddHist(histname, type1, std::array<double,6>{});
}

template<typename T>
TH1* HistogramBuilder::AddHist(std::string histname, T type1, T type2, T type3){
	return AddHist(histname, type1, std::array<double,9>{});
}


template<typename T>
TH1* HistogramBuilder::AddHist(std::string histname, T type, std::array<double,3> binning){
	
	// get a unique name
	std::string unique_name = AddHist(histname);
	if(unique_name.empty()) return nullptr;
	
	// make and build the histogram
	hists.emplace(histname, NewHist(histname, unique_name, type, binning));
	
	return hists.at(histname);
}

template<typename T>
TH2* HistogramBuilder::AddHist(std::string histname, T type, std::array<double,6> binning){
	
	// get a unique name
	std::string unique_name = AddHist(histname);
	if(unique_name.empty()) return nullptr;
	
	// make and build the histogram
	hists.emplace(histname, NewHist(histname, unique_name, type, binning));
	
	return dynamic_cast<TH2*>(hists.at(histname));
}

template<typename T>
TH3* HistogramBuilder::AddHist(std::string histname, T type, std::array<double,9> binning){
	
	// get a unique name
	std::string unique_name = AddHist(histname);
	if(unique_name.empty()) return nullptr;
	
	// make and build the histogram
	hists.emplace(histname, NewHist(histname, unique_name, type, binning));
	
	return dynamic_cast<TH3*>(hists.at(histname));
}

template<typename T>
TH1* HistogramBuilder::AddHist(std::string histname, std::vector<T> type){
	if(type.size()==1) return AddHist(histname, type[0], std::array<double,3>{});
	if(type.size()==2) return AddHist(histname, type[0], std::array<double,6>{});
	if(type.size()==3) return AddHist(histname, type[0], std::array<double,9>{});
	std::cerr<<"HistogramBuilder::AddHist called with unsupported N dims: "<<type.size()<<std::endl;
	return nullptr;
}

std::string HistogramBuilder::AddHist(std::string histname){
	
	// sanity check that this histogram doesn't exist already
	if(hists.count(histname)){
		std::cerr<<"HistogramBuilder::AddHist "<<histname<<" already exists!"<<std::endl;
		return "";
	}
	
	// this is not required, just tidies up some names
	bool trim=false;
	size_t pos = histname.find_last_of('_');
	if(pos!=std::string::npos){
		std::string tail = histname.substr(pos,std::string::npos);
		int suffix;
		int matched=sscanf(tail.c_str(),"_%d",&suffix);
		if(matched){
			trim=true;
		}
	}
	
	// make a unique name for this histogram
	std::string unique_name="";
	if(trim && histnames.count(histname)==0 && gROOT->FindObject(histname.c_str())==nullptr){
		// for the first case we can just use the user name as the unique name
		// ... not sure i like the asymmetry.
		unique_name = histname;
	} else {
		// after that we start appending random numbers
		do {
			unique_name = histname+"_"+std::to_string(hbuildercounter++);
		} while(gROOT->FindObject(unique_name.c_str())!=nullptr);
	}
	std::cout<<"HistogramBuilder::AddHist making new histogram with name '"<<histname<<"'"<<std::endl;
	histnames.emplace(histname,unique_name);
	
	return unique_name;
}
//=========================================================================================//

// Histogram factory. You know, we probably don't need all the types.
// In fact, honestly, TH*D would probably work for anything....

// float, double, long double -> use TH1D
template<typename T>
TH1D* HistogramBuilder::NewHist(std::string name, std::string unique_name, T dummy, std::array<double,3> binning, typename std::enable_if<std::is_floating_point<T>::value, bool>::type potato){
	if(binning[0]==0){
		binning[0]=200;
	}
	return new TH1D(unique_name.c_str(), name.c_str(), binning[0], binning[1], binning[2]);
}
// bool, char, short, int, long, long long, and similar -> use TH1I
template<typename T>
TH1I* HistogramBuilder::NewHist(std::string name, std::string unique_name, T dummy, std::array<double,3> binning, typename std::enable_if<std::is_integral<T>::value, bool>::type potato){
	if(binning[0]==0){
		binning[0]=200;
	}
	return new TH1I(unique_name.c_str(), name.c_str(), binning[0], binning[1], binning[2]);
}

// float, double, long double -> use TH2D
template<typename T>
TH2D* HistogramBuilder::NewHist(std::string name, std::string unique_name, T dummy, std::array<double,6> binning, typename std::enable_if<std::is_floating_point<T>::value,bool>::type potato){
	if(binning[0]==0){
		binning[0]=200;
		binning[3]=200;
	}
	return new TH2D(unique_name.c_str(), name.c_str(), binning[0], binning[1], binning[2], binning[3], binning[4], binning[5]);
}
// bool, char, short, int, long, long long, and similar -> use TH2I
template<typename T>
TH2I* HistogramBuilder::NewHist(std::string name, std::string unique_name, T dummy, std::array<double,6> binning, typename std::enable_if<std::is_integral<T>::value, bool>::type potato){
	if(binning[0]==0){
		binning[0]=200;
		binning[3]=200;
	}
	return new TH2I(unique_name.c_str(), name.c_str(), binning[0], binning[1], binning[2], binning[3], binning[4], binning[5]);
}

// float, double, long double -> use TH3D
template<typename T>
TH3D* HistogramBuilder::NewHist(std::string name, std::string unique_name, T dummy, std::array<double,9> binning, typename std::enable_if<std::is_floating_point<T>::value, bool>::type potato){
	if(binning[0]==0){
		binning[0]=200;
		binning[3]=200;
		binning[6]=200;
	}
	return new TH3D(unique_name.c_str(), name.c_str(), binning[0], binning[1], binning[2], binning[3], binning[4], binning[5], binning[6], binning[7], binning[8]);
}
// bool, char, short, int, long, long long, and similar -> use TH3I
template<typename T>
TH3I* HistogramBuilder::NewHist(std::string name, std::string unique_name, T dummy, std::array<double,9> binning, typename std::enable_if<std::is_integral<T>::value, bool>::type potato){
	if(binning[0]==0){
		binning[0]=200;
		binning[3]=200;
		binning[6]=200;
	}
	return new TH3I(unique_name.c_str(), name.c_str(), binning[0], binning[1], binning[2], binning[3], binning[4], binning[5], binning[6], binning[7], binning[8]);
}

//=========================================================================================//

// 1D
template <typename T>
bool HistogramBuilder::Fill(std::string name, T val){
	
	// ensure we have something to Fill
	if(!notree && branches.count(name)==0){
		std::cout<<"HistogramBuilder::Fill making new branch with name '"<<name<<"'"<<std::endl;
		if(!AddBranch(name, val)) return false;
	} else if(notree && hists.count(name)==0){
		if(!AddHist(name, val)) return false;
	}
	
	// fill tree if we have one
	if(branches.count(name)){
		BranchFill(name, val);
	}
	
	// fill hist if we have one
	if(hists.count(name)){
		hists.at(name)->Fill(val);
	}
	
	return true;
}

// 2D
template <typename T>
bool HistogramBuilder::Fill(std::string name, T val1, T val2){
	return Fill(name, std::vector<T>{val1, val2});
}

// 3D
template <typename T>
bool HistogramBuilder::Fill(std::string name, T val1, T val2, T val3){
	return Fill(name, std::vector<T>{val1, val2, val3});
}

// ND - maybe this should be private...?
template <typename T>
bool HistogramBuilder::Fill(std::string name, std::vector<T> vals){
	
	// sanity check
	if(vals.size()==0){
		std::cerr<<"HistogramBuilder::Fill called with name "<<name<<" but no values!"<<std::endl;
		return false;
	}
	
	// ensure we have something to Fill
	if(!notree && branches.count(name+"_0")==0){
		std::cout<<"HistogramBuilder::Fill making new branches for name '"<<name<<"'"<<std::endl;
		// two values so we'll need two branches
		for(int i=0; i<vals.size(); ++i){
			std::string uniquebranchname = name+"_"+std::to_string(i);
			if(!AddBranch(uniquebranchname, vals[0])) return false;
		}
	} else if(notree && hists.count(name)==0){
		// only supported for up to 3D
		if(!AddHist(name, vals)) return false;
	}
	
	// fill tree if we have one
	if(!notree){
		for(int i=0; i<vals.size(); ++i){
			std::string uniquebranchname = name+"_"+std::to_string(i);
			BranchFill(uniquebranchname, vals[i]);
		}
	}
	
	// fill hist if we have one
	if(hists.count(name)){
		TH1* hist = hists.at(name);
		if(hist->GetDimension()!=vals.size()){
			std::cerr<<"HistogramBuilder::Fill called for "<<name<<", but "<<vals.size()
			         <<" values given do not match dimensionality of histogram of "<<hist->GetDimension()
			         <<std::endl;
			return false;
		}
		     if(vals.size()==1){ hist->Fill(vals[0]); }
		else if(vals.size()==2){ TH2* hist2 = (TH2*)(hist); hist2->Fill(vals[0],vals[1]); }
		else if(vals.size()==3){ TH3* hist3 = (TH3*)(hist); hist3->Fill(vals[0],vals[1],vals[2]); }
	}
	
	return true;
}

template <typename T>
bool HistogramBuilder::BranchFill(std::string uniquebranchname, T val){
	
	TBranch* branchp = branches.at(uniquebranchname);
	void* branchval = branchvals.at(uniquebranchname);
	branchType branchtype = branchtypes.at(uniquebranchname);
	
	int nbytes=0;
	
	// copy passed value to our held branch variable and fill
	switch (branchtype) {
		case branchType::C: {
			char* valp = static_cast<char*>(branchval);
			*valp = val;
			nbytes = branchp->Fill();
			break;
		}
		case branchType::S: {
			short* valp = static_cast<short*>(branchval);
			*valp = val;
			nbytes = branchp->Fill();
			break;
		}
		case branchType::I: {
			int* valp = static_cast<int*>(branchval);
			*valp = val;
			nbytes = branchp->Fill();
			break;
		}
		case branchType::F: {
			float* valp = static_cast<float*>(branchval);
			*valp = val;
			nbytes = branchp->Fill();
			break;
		}
		case branchType::D: {
			double* valp = static_cast<double*>(branchval);
			*valp = val;
			nbytes = branchp->Fill();
			break;
		}
	}
	return nbytes;
}


#endif
