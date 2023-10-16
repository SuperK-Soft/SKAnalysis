#include "AddTree.h"

#include "MTreeReader.h"
#include "TTree.h"

AddTree::AddTree():Tool(){}


bool AddTree::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	std::string renameDataTo, newTreeName, treeReaderName;
	m_variables.Get("renameDataTo", renameDataTo);
	m_variables.Get("newTreeName", newTreeName);
	m_variables.Get("treeReaderName", treeReaderName);
	int lun2=0;
	m_variables.Get("LUN",lun2);
	// give it a name to put it in the datamodel map for downstream Tools
	std::string treeWriterName;
	m_variables.Get("treeWriterName", treeWriterName);
	
	// set the associated file as active
	int lun1 = m_data->GetLUN(treeReaderName);
	TreeManager* mgr1 = skroot_get_mgr(&lun1);
	if(mgr1==nullptr){
		Log(m_unique_name+" Error! Couldn't find manager associated with existing tree '"+treeReaderName+"'",
		    v_error,m_verbose);
		return false;
	}
	TTree* otree1 = mgr1->GetOTree();
	if(otree1==nullptr){
		Log(m_unique_name+" Error! Null output tree for existing tree '"+treeReaderName+"'",
		    v_error,m_verbose);
		return false;
	}
	// to prevent clash we may wish to rename the existing hard-coded 'data' tree
	if(renameDataTo!=""){
		otree1->SetName(renameDataTo.c_str());
	}
	ofile = otree1->GetCurrentFile();
	if(ofile==nullptr){
		Log(m_unique_name+" Error! No output file associated to existing output tree?!",v_error,m_verbose);
		return false;
	}
	ofile->cd();
	
	// make a new TreeManager in Write mode, but give it an empty file name.
	// this causes the internal call to 'new TFile(filename, "RECREATE")'
	// to fail, but this is not checked, so the code subsequently carries on,
	// making a new managed Tree, associated with the current ROOT directory
	// (i.e. our existing file)
	lun2 = m_data->GetNextLUN(treeWriterName, lun2);
	Log(m_unique_name+": The following error '<TFile::TFile>: file name is not specified'"
	                  " may be safely ignored",v_warning,m_verbose);
	skroot_open_write_(&lun2, "", 0);
	TreeManager* mgr2 = skroot_get_mgr(&lun2);
	if(mgr2==nullptr){
		Log(m_unique_name+" Error! TreeManager not found!",v_error,m_verbose);
		return false;
	}
	thistree = mgr2->GetOTree();
	if(thistree==nullptr){
		Log(m_unique_name+" Error! New output tree is null!",v_error,m_verbose);
		return false;
	}
	thistree->SetName(newTreeName.c_str());
	// don't think this is actually needed
	//thistree->SetDirectory(ofile);
	
	return true;
}


bool AddTree::Execute(){
	
	return true;
}


bool AddTree::Finalise(){
	
	return true;
}

AddTree::~AddTree(){
	// the last thing is we need to make sure that this new tree gets written
	// to disk before the parent file is closed. The parent file is managed by
	// another TreeManager, which in turn is owned by the SuperManager singleton
	// (which creates it during skroot_open*).
	// When the application terminates the SuperManager's destructor will call the
	// destructor of all the TreeManagers, causing each to Write and Close their files.
	// So, it should be fine if we call Write on this Tree in Finalise, because the
	// parent file will not be Closed and deleted until the application terminates.
	// ...
	// unless we close it manually! Whoops! CloseLUN now commented out in TreeReader::Finalise...
	// FIXME find a better solution?
	std::string newTreeName;
	m_variables.Get("newTreeName", newTreeName);
	//std::cout<<m_unique_name<<" Destructor: writing out tree "<<newTreeName<<" which has "<<thistree->GetEntries()
	//         <<" entries to file "<<ofile->GetName()<<" which is zombie: "<<ofile->IsZombie()<<std::endl;
	
	ofile->cd();
	int nbytes = thistree->Write("",TObject::kOverwrite);
	//std::cout<<"wrote "<<nbytes<<" bytes"<<std::endl;
}
