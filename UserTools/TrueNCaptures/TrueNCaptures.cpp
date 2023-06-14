#include "TrueNCaptures.h"

TrueNCaptures::TrueNCaptures():Tool(){
	// get the name of the tool from its class name
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}


bool TrueNCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	m_variables.Get("verbosity",verbosity);
	
	 return true;
}


bool TrueNCaptures::Execute(){
	
	m_data->NCapturesTrue.clear();
	
	// loop over true particles and pull out neutron capture events
	//std::cout<<"Scanning "<<m_data->eventParticles.size()<<" for neutrons"<<std::endl;
	for(int i=0; i<m_data->eventParticles.size(); ++i){
		MParticle& aparticle = m_data->eventParticles.at(i);
		//std::cout<<"particle "<<i<<" has pdg "<<aparticle.pdg<<std::endl;
		// check if it's a neutron
		if(aparticle.pdg==2112){
			//std::cout<<"neutron at index "<<i<<std::endl;
			bool captured=false;
			// get termination processes
			if(aparticle.GetEndProcesses()){
				// check if it captured
				for(int j=0; j<aparticle.GetEndProcesses()->size(); ++j){
					int process_code = aparticle.GetEndProcesses()->at(j);
					if(process_code==18){
						captured=true;
						break;
					}
				}
			} else {
				// uhh, we have no process codes for its termination.
				// no way to know if it captured or not.... assume yes? XXX
				Log(m_unique_name+" found neutron but no termination processes recorded!"
				   " Tentatively assuming it underwent capture...",v_warning,verbosity);
				captured=true;
			}
			if(!captured){
				Log(m_unique_name+" found true neutron, but it did not capture",v_debug,verbosity);
				continue;
			}
			Log(m_unique_name+" found true neutron capture",v_debug,verbosity);
			NCapture acapture;
			acapture.SetNeutronIndex(i);
			m_data->NCapturesTrue.push_back(acapture);
		}
	}
	
	PrintCaptures();
	
	return true;
}

bool TrueNCaptures::Finalise(){
	
	 return true;
}

bool TrueNCaptures::PrintCaptures(){
	std::cout<<"This event contained "<<m_data->NCapturesTrue.size()
	         <<" recorded true neutron captures"
	         <<"\n==========================================\n";
	for(int i=0; i<m_data->NCapturesTrue.size(); ++i){
		if(i>0) std::cout<<"------------------------------------------\n";
		std::cout<<"Capture "<<i<<"\n";
		NCapture& acapture = m_data->NCapturesTrue.at(i);
		acapture.Print();
	}
	std::cout<<"=========================================="<<std::endl;
	return true;
}

bool TrueNCaptures::MakePlots(){
	// TODO
	// pi chart of capture nuclide
	// following: stack by nuclide, gd by isotopes then H
	// distribution of neutron travel time
	// distribution of neutron travel distance
	// distribution of number of daughter gammas
	// distribution of total gamma enery
	// distribution of individual gamma energy
	// distribution of num conversion electrons
	// distribution of total electron energy
	// distribution of individual electron energy
	
	return true;
}

