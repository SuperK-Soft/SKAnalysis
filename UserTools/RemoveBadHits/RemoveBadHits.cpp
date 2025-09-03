#include "RemoveBadHits.h"

RemoveBadHits::RemoveBadHits():Tool(){}


bool RemoveBadHits::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	
	std::string tree_reader_str = "";
	m_variables.Get("reader", tree_reader_str);
	if(m_data->Trees.count(tree_reader_str) == 0){
		throw std::runtime_error(m_unique_name+" Failed to get treereader "+tree_reader_str+"!");
	}
	tree_reader = m_data->Trees.at(tree_reader_str);
	
	
	return true;
}


bool RemoveBadHits::Execute(){
	
	int hit_idx;
	for (hit_idx = 0; hit_idx < sktqz_.nqiskz; ++hit_idx){
		if(sktqz_.tiskz[hit_idx]>600E3) break;
	}
	Log(m_unique_name+": truncating nqiskz from "+toString(sktqz_.nqiskz)+" to "+toString(hit_idx-1),v_debug,m_verbose);
	sktqz_.nqiskz=hit_idx-1;
	rawtqinfo_.nqisk_raw=hit_idx-1;
	
	TQReal* TQREAL = nullptr;
	tree_reader->Get("TQREAL", TQREAL);
	TQREAL->nhits=sktqz_.nqiskz;
	TQREAL->T.resize(sktqz_.nqiskz);
	TQREAL->Q.resize(sktqz_.nqiskz);
	TQREAL->cables.resize(sktqz_.nqiskz);
	
	return true;
}


bool RemoveBadHits::Finalise(){
	
	return true;
}
