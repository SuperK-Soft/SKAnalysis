#include "MVertex.h"
#include "MParticle.h"
#include "DataModel.h"

MVertex::MVertex(){
	m_data = DataModel::GetInstance();
	extraInfo = new BStore{true,constants::BSTORE_BINARY_FORMAT};
}

MVertex::~MVertex(){
	if(extraInfo) delete extraInfo;
	extraInfo=nullptr;
}

MVertex::MVertex(MVertex&& rhs){
	
	m_data = rhs.m_data;
	pos = rhs.pos;
	time = rhs.time;
	type = rhs.type;
	target_pdg = rhs.target_pdg;
	processes = rhs.processes;
	incident_particle_idx = rhs.incident_particle_idx;
	direct_parent = rhs.direct_parent;
	incident_particle_pdg = rhs.incident_particle_pdg;
	incident_particle_mom = rhs.incident_particle_mom;
	
	// reason for manually doing all this:
	extraInfo = rhs.extraInfo;
	rhs.extraInfo = nullptr;
	
}

void MVertex::SetIncidentParticle(int idx){
	// if -1, unrecorded.
	// if 0+, index of the direct parent.
	// If <0, absolute value - 2 is index of the nearest recorded parent.
	if(idx<-1){
		idx = std::abs(idx) - 2;
		direct_parent = false;
	} else {
		direct_parent = true;
	}
	incident_particle_idx = idx;
}

MParticle* MVertex::GetIncidentParticle(){
	if(incident_particle_idx>=0 && incident_particle_idx<m_data->eventParticles.size()){
		return &m_data->eventParticles.at(incident_particle_idx);
	}
	// else no associated particle
	return nullptr;
}

bool MVertex::IsParentDirect(){
	return direct_parent;
}

std::string MVertex::PrintProcesses(){
	if(processes.size()==0){
		return "{?}";
	}
	std::string processnames = "{";
	for(int i=0; i<processes.size(); ++i){
		if(i>0) processnames+= ", ";
		int aprocid = processes.at(i);
		if(aprocid<0){
			// G4 process ID
			processnames += G4_process_code_to_string(std::abs(aprocid));
		} else {
			processnames += G3_process_code_to_string(aprocid);
		}
	}
	processnames += "}";
	return processnames;
}

void MVertex::Print(bool verbose){
	if(verbose){
		std::cout<<"eventVertices index: ";
		bool found=false;
		for(int i=0; i<m_data->eventVertices.size(); ++i){
			if(this==&m_data->eventVertices.at(i)){
				std::cout<<i<<"\n";
				found = true;
				break;
			}
		}
		if(!found) std::cout<<"?\n";
	}
	std::cout<<"\tPosition [cm]: ("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")\n"
	         <<"\tTime [ns]: "<<time<<"\n"
	         <<"\tType: "<<type<<"\n"
	         <<"\tProcsses: "<<PrintProcesses()<<"\n";
	MParticle* ipart = GetIncidentParticle();
	if(ipart==nullptr){
		std::cout<<"\tIncident particle: Not Set\n";
	} else {
		std::cout<<"\tIncident particle idx: "<<incident_particle_idx<<"\n"
		         <<"\tDirect parent? "<<direct_parent<<"\n";
	}
	std::cout<<"\tTrue incident particle pdg: "<<incident_particle_pdg<<"\n"
	         <<"\tTrue incident particle momentum [MeV/c]: "<<toString(incident_particle_mom)<<"\n";
	//std::cout<<"\tExtra Info: "; extraInfo->Print(false);
}

