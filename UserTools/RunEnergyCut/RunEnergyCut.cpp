#include "RunEnergyCut.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "skroot_loweC.h"

RunEnergyCut::RunEnergyCut():Tool(){}


bool RunEnergyCut::Initialise(std::string configfile, DataModel &data){
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	
	m_variables.Get("treeReaderName", treeReaderName);
	
	if(m_data->Trees.count(treeReaderName)==0){
	Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
	return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	return true;
}


bool RunEnergyCut::Execute(){
	
	bool muon = false;
	m_data->vars.Get("newMuon", muon);
	if(muon){
		return true;
	}
	
	myTreeReader->Get("HEADER", myHeader);
	myTreeReader->Get("LOWE", myLowe);
	
	float reconEnergy = skroot_lowe_.bsenergy;
	
	//	float reconEnergy = myLowe->bsenergy;
	
	// set upper limit at 100 MeV
	if(reconEnergy > 100.){
		Nskipped++;
		m_data->vars.Set("Skip", true);
		return true;
	}
	// check for run number. At run 68671 the SHE threshold was decreased from 70 hits to 68 hits, so the lower
	// limit for the search is reduced by 2 MeV.
	if(myHeader->nrunsk < 68671){
		if(reconEnergy < 10.){
			//std::cout << "Event " << myHeader->nevsk << " skipped by energy: " << reconEnergy << std::endl;
			Nskipped++;
			m_data->vars.Set("Skip", true);
			return true;
		}
	}else{
		if(reconEnergy < 8.){
			//std::cout << "Event " << myHeader->nevsk << " skipped by energy: " << reconEnergy << std::endl;
			Nskipped++;
			m_data->vars.Set("Skip", true);
			return true;
		}
	}
	
	//	std::cout << "Passed run energy: " << skroot_lowe_.bsenergy << std::endl;
	
	return true;
}


bool RunEnergyCut::Finalise(){
	
	std::cout << "Number of events skipped due to RunEnergy cut: " << Nskipped << std::endl;
	
	return true;
}
