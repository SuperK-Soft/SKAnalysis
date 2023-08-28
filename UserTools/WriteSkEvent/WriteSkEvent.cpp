#include "WriteSkEvent.h"

WriteSkEvent::WriteSkEvent():Tool(){}


bool WriteSkEvent::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	std::string treeReaderName;
	get_ok = m_variables.Get("treeReaderName", treeReaderName);
	m_variables.Get("delete_outside_hits", delete_outside_hits);
	
	// we don't actually need the reader itself, just its logic unit number (LUN)
	LUN = m_data->GetLUN(treeReaderName);
	if(LUN<0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1); // fatal error
		return false;
	}
	
	// sanity check
	TreeManager* mgr = skroot_get_mgr(&LUN);
	if(mgr->GetMode()==2){
		// file is opened in READ mode! no output file to write!
		Log(m_unique_name+" Error! TreeReader "+treeReaderName+" is opened in SKROOT_READ mode,"
		    " so cannot write outputs to it!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1); // fatal error
		return false;
	}
	
	return true;
}


bool WriteSkEvent::Execute(){
	
	Log(m_unique_name+" Executing...",v_debug,m_verbose);
	
	// remove hits outside 1.3 microsec window from primary trigger
	if(delete_outside_hits) delete_outside_hits_();
	
	// pass information from common blocks to output TTree
	skroot_set_tree_(&LUN);
	
	// skroot_set_tree.F is a wrapper that calls:
	// * skroot_set_header
	// * skroot_set_tqreal
	// * skroot_set_tqareal
	
	// others include:
	// - skroot_set_lowe_
	// - skroot_set_mu_
	// - skroot_set_upmu_
	// - skroot_set_atmpd_
	// - skroot_set_sle_
	// - skroot_set_prevt0_
	// - skroot_set_mc_
	// - skroot_set_secondary_
	// - skroot_set_softtrglist_
	// we leave these to other Tools to call when appropriate.
	
	// FIXME should we also leave skroot_set_tree to another Tool?
	// on the one hand, this is general information: which Tool should do it?
	// on the other, these routines pull information from the global common blocks,
	// which are used by and manipulated by all manner of arbitrary fortan routines.
	// So results from one Tool may, intentionally or not, be overwritten by another downstream Tool.
	// Does requiring Tools to call `skroot_set_*` give better control over what is going into the file?
	
	// XXX note that skroot_set_tqreal/tqareal will need to be called again after delete_outside_hits
	
	// invokes TTree::Fill
	skroot_fill_tree_(&LUN);
	
	return true;
}


bool WriteSkEvent::Finalise(){
	
	return true;
}
