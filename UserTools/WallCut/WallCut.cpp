#include "WallCut.h"
#include "geotnkC.h" // for HIINTK, DIINTK
#include "fortran_routines.h"
#include "Constants.h"
#include "skroot_loweC.h"

#include <cmath>

WallCut::WallCut():Tool(){}


bool WallCut::Initialise(std::string configfile, DataModel &data){
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("cutThreshold",cutThreshold);
	
	// add a cut to the selector if being used
	get_ok = m_variables.Get("selectorName", selectorName);
	if(get_ok){
		std::string description = "bswall < " + toString(cutThreshold);
		m_data->AddCut(selectorName, m_unique_name, description,true);
	}
	
	return true;
}


bool WallCut::Execute(){
	
	EventType eventType;
	m_data->vars.Get("eventType", eventType);
	if(eventType!=EventType::LowE) return true;
	
	float wallDistance = 0;
	
	// distance to closest wall
	//wallDistance = skroot_lowe_.linfo[8];    // skroot_lowe_clwallsk
	wallDistance = skroot_lowe_.linfo[9];      // skroot_lowe_bswallsk
	// or if it hasn't been calculated yet
	//wallDistance = bswallsk = wallsk_(&skroot_lowe_.bsvertex[0]);
	
	// back-projected distance to wall
	//wallDistance = skroot_lowe_.linfo[5];     // skroot_lowe_bseffwal
	//wallDistance = skroot_lowe_.linfo[3];     // skroot_lowe_cleffwal
	
	// deprecated manual version
	//wallDistance = DistanceToWall();
	
	bool rejected=false;
	if(wallDistance < cutThreshold){
		rejected=true;
		Nskipped++;
		m_data->vars.Set("Skip",true);
		//std::cout << "Skipped due to wall" << std::endl;
		return true;
	}
	
	if(!selectorName.empty() && !rejected) m_data->ApplyCut(selectorName, m_unique_name, wallDistance);
	
	return true;
}

float WallCut::DistanceToWall(){
	
	float* reconVertex = &skroot_lowe_.bsvertex[0];
	
	static const float IVRadius = DIINTK/2.f;   // inner volume radius, 1690 cm
	static const float posHeight = HIINTK/2.f;  // inner volume half-height, 1810 cm
	
	float xyDistance = IVRadius - sqrt(pow(reconVertex[0], 2.f) + pow(reconVertex[1], 2.f));
	float zDistance = posHeight - abs(reconVertex[2]);
	
	return std::min(xyDistance, zDistance);
	
}


bool WallCut::Finalise(){
	
	Log(m_unique_name+": Number of events rejected: "+toString(Nskipped),v_debug,m_verbose);
	
	return true;
}
