#include "MVertex.h"
#include "MParticle.h"
#include "DataModel.h"

MVertex::MVertex(){
	m_data = DataModel::GetInstance();
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

// TODO add Print method
