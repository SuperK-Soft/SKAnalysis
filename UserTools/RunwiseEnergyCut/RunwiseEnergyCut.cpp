#include "RunwiseEnergyCut.h"
#include "fortran_routines.h"
#include "Constants.h"

RunwiseEnergyCut::RunwiseEnergyCut():Tool(){}


bool RunwiseEnergyCut::Initialise(std::string configfile, DataModel &data){
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	
	if(!ParseOptions(configfile)){
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// add a cut to the selector if being used
	get_ok = m_variables.Get("selectorName", selectorName);
	if(get_ok){
		std::stringstream description{"runwise energy cuts:\n"};
		for(auto&& acut : cuts){
			description<<"{"<<acut.first.first<<" < Run < "<<acut.first.second
			         <<" : "<<acut.second.first<<" < bsenergy < "<<acut.second.second<<"}\n";
		}
		m_data->AddCut(selectorName, m_unique_name, description.str());
	}
	
	return true;
}


bool RunwiseEnergyCut::Execute(){
	
	bool muon = false;
	m_data->vars.Get("newMuon", muon);
	if(muon){
		return true;
	}
	
	float reconEnergy = skroot_lowe_.bsenergy;
	
	bool rejected=false;
	for(auto&& acut : cuts){
		int runmin = acut.first.first;
		int runmax = acut.first.second;
		if(runmin>0 && runmin <= skhead_.nrunsk) continue;
		if(runmax>0 && runmax >= skhead_.nrunsk) continue;
		if(reconEnergy < acut.second.first){  rejected=true; break; }
		if(reconEnergy > acut.second.second){ rejected=true; break; }
	}
	
	if(!rejected){
		Nskipped++;
		m_data->vars.Set("Skip", true);
		return true;
	}
	
	if(!selectorName.empty() && !rejected) m_data->ApplyCut(selectorName, m_unique_name, reconEnergy);
	
	Log(m_unique_name+" Event passed with energy: "+toString(skroot_lowe_.bsenergy),v_debug,m_verbose);
	
	return true;
}


bool RunwiseEnergyCut::Finalise(){
	
	Log(m_unique_name+": Number of events skipped due to RunEnergy cut: "+toString(Nskipped),v_debug,m_verbose);
	
	return true;
}

bool RunwiseEnergyCut::ParseOptions(std::string& configfile){
	
	std::string line;
	std::ifstream infile(configfile);
	if(!infile.is_open()){
		Log(m_unique_name+" Error opening config file "+configfile,v_error,m_verbose);
		return false;
	}
	std::stringstream ss;
	std::string key;
	int startrun, endrun;
	float minE, maxE;
	while(getline(infile, line)){
		if(line.empty()) continue;
		ss.str(line);
		if(!(ss >> key)) break; // ?? what happened?
		if(key[0]=='#') continue;
		if(key!="cut") continue;
		if(!(ss >> startrun >> endrun >> minE >> maxE)) continue;
		cuts.emplace_back(std::pair<int,int>{startrun, endrun}, std::pair<int,int>{minE, maxE});
	}
	
	if(m_verbose >= v_debug){
		std::cout<<m_unique_name<<" got cuts:\n";
		for(auto&& acut : cuts){
			std::cout<<acut.first.first<<" < Run < "<<acut.first.second
			         <<" : Require "<<acut.second.first<<" < bsenergy < "<<acut.second.second<<"\n";
		}
		std::cout<<std::flush;
	}
	
	if(cuts.empty()){
		Log(m_unique_name+" Error! Found no valid cut specifications in config file!",v_error,m_verbose);
		return false;
	}
	
	return true;
}
