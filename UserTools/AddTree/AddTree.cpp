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
	
	if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+" failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	TTree* existingTree = myTreeReader->GetTree();
	
	// to prevent clash we may wish to rename the existing hard-coded 'data' tree
	if(renameDataTo!="") existingTree->SetName(renameDataTo.c_str());
	
	// set the associated file as active
	myTreeReader->GetFile()->cd();
	
	// make a new TreeManager in Write mode, but give it an empty file name.
	// this causes the internal call to 'new TFile(filename, "RECREATE")'
	// to fail, but this is not checked, so the code subsequently carries on,
	// making a new managed Tree, associated with the current ROOT directory
	// (i.e. our existing file)
	int LUN = m_data->GetNextLUN(LUN, treeReaderName);
	skroot_open_write_(&LUN, "", 0);
	TreeManager* mgr = skroot_get_mgr(&LUN);
	thistree = mgr->GetTree();
	thistree->SetName(newTreeName.c_str());
	// don't think this is actually needed
	//thistree->SetDirectory(myTreeReader->GetFile());
	
	return true;
}


bool AddTree::Execute(){
	
	return true;
}


bool AddTree::Finalise(){
	// the last thing is we need to make sure that this new tree gets written
	// to disk before the parent file is closed. The parent file is managed by
	// another TreeManager, which in turn is owned by the SuperManager singleton
	// (which creates it during skroot_open*).
	// When the application terminates the SuperManager's destructor will call the
	// destructor of all the TreeManagers, causing each to Write and Close their files.
	// So, it should be fine if we call Write on this Tree in Finalise, because the
	// parent file will not be Closed and deleted until the application terminates.
	myTreeReader->GetFile()->cd();
	thistree->Write("",TObject::kOverwrite);
	
	return true;
}
