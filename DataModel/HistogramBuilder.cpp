#include "HistogramBuilder.h"

size_t HistogramBuilder::hbuildercounter = 0;

const std::map<std::string, HistogramBuilder::branchType> HistogramBuilder::typechars{
	{"char", HistogramBuilder::branchType::C},
	{"unsigned char",HistogramBuilder::branchType::C},
	{"short",HistogramBuilder::branchType::S},
	{"unsigned short",HistogramBuilder::branchType::S},
	{"int",HistogramBuilder::branchType::I},
	{"unsigned int",HistogramBuilder::branchType::I},
	{"long",HistogramBuilder::branchType::I},
	{"unsigned long",HistogramBuilder::branchType::I},
	{"long long",HistogramBuilder::branchType::D},
	{"unsigned long long",HistogramBuilder::branchType::D},
	{"float",HistogramBuilder::branchType::F},
	{"double",HistogramBuilder::branchType::D}
};

HistogramBuilder::HistogramBuilder(){}

TFile* HistogramBuilder::MakeFile(std::string filename, std::string treename, bool notree_in){
	// check if we already have a file open
	if(ofile!=nullptr){
		std::cerr<<"HistogramBuilder::MakeFile called to make file '"<<filename<<"'"
		         <<" but this HistogramBuilder is already associated with file '"
		         <<ofile->GetName()<<"'"<<std::endl;
		return nullptr;
	}
	// attempt to make the requested file
	if(filename!=""){
		// if filename given
		ofilename = filename;
		ofile = new TFile(ofilename.c_str(),"RECREATE");
		if(ofile==nullptr || ofile->IsZombie()){
			std::cerr<<"HistogramBuilder: Error making new file '"<<filename<<"'"<<std::endl;
			if(ofile) ofile->Close();
			delete ofile;
			ofile=nullptr;
			return nullptr;
		}
	} else {
		// no name given: make some temporary output file...?
		ofilename = std::string(std::tmpnam(nullptr));
		ofile = new TFile(ofilename.c_str(),"CREATE");
		if(ofile==nullptr || ofile->IsZombie()){
			std::cerr<<"HistogramBuilder::MakeFile error attempting to make temporary file '"
			         <<ofilename<<"'"<<std::endl;
			if(ofile) ofile->Close();
			delete ofile;
			ofile=nullptr;
			return nullptr;
		}
		std::cerr<<"HistogramBuilder::MakeFile created temporary file '"<<ofilename<<"'"<<std::endl;
		fileistmp=true;
	}
	notree = notree_in;
	if(!notree){
		MakeTree(treename);
	}
	
	return ofile;
}

TTree* HistogramBuilder::MakeTree(std::string treename){
	notree=false;
	
	// ensure we have a file open
	if(ofile==nullptr){
		std::cerr<<"HistogramBuilder::NameTree needs an output file! Call MakeFile first"<<std::endl;
		return nullptr;
	}
	
	if(tree==nullptr){
		ofile->cd();
		tree = new TTree(treename.c_str(),treename.c_str());
	} else {
		tree->SetName(treename.c_str());
	}
	
	return tree;
}

int HistogramBuilder::SetTreeEntries(){
	// returns num entries
	if(tree) return tree->SetEntries(-1);  // update num entries in tree based on num entries in branches
	return 0;
}

bool HistogramBuilder::Save(){
	
	if(ofile==nullptr){
		std::cerr<<"HistogramBuilder::Save called but no output file!"<<std::endl;
		return false;
	}
	ofile->cd();
	if(tree){
		SetTreeEntries();
		tree->Write("",TObject::kOverwrite);
	}
	if(savehists){
		for(auto&& ahist : hists){
			ahist.second->Write("",TObject::kOverwrite);
		}
	}
	return true;
}

bool HistogramBuilder::Close(){
	if(tree || ofile) Save();
	if(tree) tree->ResetBranchAddresses();
	if(ofile) ofile->Close();
	// these are owned by the file...
	/*
	if(m_verbose) std::cout<<"HistogramBuilder: deleting histograms"<<std::endl;
	for(auto&& ahist : hists){
		delete ahist.second;
	}
	std::cout<<"HistogramBuilder: histograms deleted"<<std::endl;
	*/
	hists.clear();
	delete ofile;
	ofile=nullptr;
	tree=nullptr;
	if(m_verbose) std::cout<<"HistogramBuilder: Close done"<<std::endl;
	return true;
}

HistogramBuilder::~HistogramBuilder(){
	
	Close();
	
	if(fileistmp && !ofilename.empty()){
		if(m_verbose) std::cerr<<"HistogramBuilder: removing temporary output file "<<ofilename<<std::endl;
		std::string cmd="rm '"+ofilename+"'";
		system(cmd.c_str());
		std::cout<<"HistogramBuilder: tempfile removed"<<std::endl;
	}
	
}

void HistogramBuilder::SetVerbosity(int verb){
	m_verbose = verb;
	return;
}

void HistogramBuilder::SaveHists(bool dosavehists){
	savehists=dosavehists;
	return;
}

TFile* HistogramBuilder::GetFile(){
	return ofile;
}

TTree* HistogramBuilder::GetTree(){
	// user needs to bear in mind that, depending on when they call this,
	// since we fill branches independently, the number of entries in the tree
	// may not be reported correctly!
	SetTreeEntries();
	return tree;
}

size_t HistogramBuilder::GetNHists(){
	return hists.size();
}

std::map<std::string,TH1*> HistogramBuilder::GetHists(){
	return hists;
}

TH1* HistogramBuilder::GetHist(std::string name, std::string cut, int unbinned){
	
	//std::cout<<"getting hist '"<<name<<"' with cut '"<<cut<<"', unbinned: "<<unbinned<<std::endl;
	
	std::string histname = name;
	std::string branchname = name;
	
	// see if the user specified a branchname and histogram name
	size_t pos = name.find(">>");
	if(pos!=std::string::npos){
		branchname = name.substr(0,pos);
		histname = name.substr(pos+2,std::string::npos);
	}
	
	// if this histogram already exists...
	if(hists.count(histname)>0){
		//std::cout<<"known hist"<<std::endl;
		// did the user pass options suggesting they wanted us to *make* this histogram?
		if(!cut.empty() || unbinned!=-1){
			
			/*
			//std::cerr<<"HistogramBuilder::GetHist called to create histogram '"<<histname<<"' "
			//	     <<" with options, but it already exists!"<<std::endl;
			
			// XXX what do we do in this situation? Do we...
			// 1. apply the options the user has specified, clobbering the current histogram?
			// 2. return the one that exists, ignoring the options
			// 3. return nullptr, so that the user does not accidentally get something they don't expect
			// 4. make a new histogram with a unique name and the options given <<< choose this
			
			// option 1: return the one that exists
			return hists.at(histname);
			
			// option 2: reset and refill (effectively replacing) it?
			TH1* thehist = hists->at(histname);
			thehist->Reset();
			thehist->Clear();
			int dims = thehist->GetDimension();
			switch (dims) {
				case 1: thehist->SetBins(200,0,0); break;
				case 2: thehist->SetBins(200,0,0,200,0,0); break;
				case 3: thehist->SetBins(200,0,0,200,0,0,200,0,0); break;
			}
			
			// option 3: return null
			return nullptr;
			*/
			
			// we choose option 4: make a new histogram with the given options
			// we will need a unique *user* name to use as a key for the histnames map.
			int i=0;
			do {
				histname = histname+"_"+std::to_string(hbuildercounter+i);
				++i;
			} while(gROOT->FindObject(histname.c_str())!=nullptr);
			
		} else {
			// otherwise no options given, the histogram already exists, return it
			return hists.at(histname);
		}
	}
	//std::cout<<"unknown hist"<<std::endl;
	
	// histogram doesn't exist: time to make it
	if(branches.count(branchname)>0){
		
		//std::cout<<"known branch"<<std::endl;
		// make 1D histogram from a known branch
		// get branch type
		if(branchtypes.count(branchname)==0){ // should not happen
			std::cerr<<"don't know the type of branch '"<<branchname<<"'"<<std::endl;
			return nullptr;
		}
		branchType datatype = branchtypes.at(branchname);
		
		// make the histogram (this gives it a unique name)
		if((datatype==branchType::F) || (datatype==branchType::D)){
			if(!AddHist(histname, double(0))) return nullptr;
		} else {
			if(!AddHist(histname, int(0))) return nullptr;
		}
		
		// fill it
		std::string drawcmd = branchname+">>"+histnames.at(histname)+"(200)";
		//std::cout<<"drawcmd: '"<<drawcmd<<"'"<<std::endl;
		SetTreeEntries(); // ensure our tree knows how many entries it has
		int ret = tree->Draw(drawcmd.c_str(), cut.c_str());
		//std::cout<<"draw returned "<<ret<<std::endl;
		
		// error checks
		if(ret<0){
			std::cerr<<"HistogramBuilder::GetHist error "<<ret<<" drawing '"<<drawcmd<<"'"<<std::endl;
		} else if(ret==0){
			// n.b. if no events are found, no histogram will be made!
			// this is partly why we make the histogram ourselves first. (also to give it a unique name).
			std::cerr<<"HistogramBuilder::GetHist found no entries from command '"<<drawcmd<<"'"<<std::endl;
		}
		
		return hists.at(histname);
		
	} else if(branches.count(branchname+"_0")){
		//std::cout<<"found "<<branchname<<"_0 branch"<<std::endl;
		
		// we have a branch with a name that suggests this is a multi-dimensional parameter,
		// split across multiple branches. scan for other branches to determine dimensionality
		// and make a suitable dimensional histogram.
		if(branches.count(branchname+"_1")==0){
			// if it's multidimensional we ought to have a second dimension at least...?
			std::cerr<<"HistogramBuilder::GetHist requesting histogram from branch '"<<branchname
			         <<"'; found branch "+branchname+"_0, but not "+branchname+"_1?"<<std::endl;
			return nullptr;
		} else if(branches.count(branchname+"_2")==0){
			// make 2D histogram
			
			// get branch type (assume they are all the same)
			if(branchtypes.count(branchname+"_0")==0){ // should not happen
				std::cerr<<"don't know the type of branch '"<<branchname<<"_0'"<<std::endl;
				return nullptr;
			}
			branchType datatype = branchtypes.at(branchname+"_0");
			
			// make the histogram (this gives it a unique name)
			if((datatype==branchType::F) || (datatype==branchType::D)){
				if(!AddHist(histname, double(0), double(0))) return nullptr;
			} else {
				if(!AddHist(histname, int(0), int(0))) return nullptr;
			}
			
			// fill it
			std::string drawcmd = branchname+"_0:"+branchname+"_1>>"+histnames.at(histname);
			if(unbinned==0) drawcmd += "(200,0,0,200,0,0)";
			SetTreeEntries(); // ensure our tree knows how many entries it has
			int ret = tree->Draw(drawcmd.c_str(),cut.c_str());
			
			// error checks
			if(ret<0){
				std::cerr<<"HistogramBuilder::GetHist error "<<ret<<" drawing '"<<drawcmd<<"'"<<std::endl;
			} else if(ret==0){
				// n.b. if no events are found, no histogram will be made!
				// this is partly why we make the histogram ourselves first. (also to give it a unique name).
				std::cerr<<"HistogramBuilder::GetHist found no entries from command '"<<drawcmd<<"'"<<std::endl;
			}
			
			return hists.at(histname);
			
		} else if(branches.count(branchname+"_3")==0){
			
			// get branch type (assume they are all the same)
			if(branchtypes.count(branchname+"_0")==0){ // should not happen
				std::cerr<<"don't know the type of branch '"<<branchname<<"_0'"<<std::endl;
				return nullptr;
			}
			branchType datatype = branchtypes.at(branchname+"_0");
			
			// make the histogram (this gives it a unique name)
			if((datatype==branchType::F) || (datatype==branchType::D)){
				if(!AddHist(histname, double(0), double(0), double(0))) return nullptr;
			} else {
				if(!AddHist(histname, int(0), int(0), int(0))) return nullptr;
			}
			
			// fill it
			std::string drawcmd = branchname+"_0:"+branchname+"_1:"+branchname+"_2>>"+histnames.at(histname);
			if(unbinned==0) drawcmd += "(200,0,0,200,0,0,200,0,0)";
			SetTreeEntries(); // ensure our tree knows how many entries it has
			int ret = tree->Draw(drawcmd.c_str(),cut.c_str());
			
			// error checks
			if(ret<0){
				std::cerr<<"HistogramBuilder::GetHist error "<<ret<<" drawing '"<<drawcmd<<"'"<<std::endl;
			} else if(ret==0){
				// n.b. if no events are found, no histogram will be made!
				// this is partly why we make the histogram ourselves first. (also to give it a unique name).
				std::cerr<<"HistogramBuilder::GetHist found no entries from command '"<<drawcmd<<"'"<<std::endl;
			}
			
			return hists.at(histname);
		} else {
			std::cerr<<"HistogramBuilder::GetHist found too many branches associated with name '"
			         <<branchname<<"' - 4D histograms are not supported"<<std::endl;
			return nullptr;
		}
	} else {
		// else no histogram of this name, no branch of this name, nor any branch with a name
		// suggesting it's a mutli-dimensional parameter. No idea what they're asking for.
		std::cerr<<"HistogramBuilder::GetHist Error! Request for unknown variable "<<name<<std::endl;
		return nullptr;
	}
	
	return nullptr;
}

TH1* HistogramBuilder::GetHist(size_t i){
	if(i<hists.size()) return std::next(hists.begin(),i)->second;
	std::cerr<<"HistogramBuilder::GetHist Error! Request for histogram "<<i<<" / "
	         <<hists.size()<<std::endl;
	return nullptr;
}


