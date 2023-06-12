#ifndef MVERTEX_H
#define MVERTEX_H
#include "TVector3.h"
#include "BStore.h"
#include <vector>

class MParticle;
class DataModel;

class MVertex {
	public:
	MVertex();
	// FIXME should probably make these arguments of the constructor so that they're mandatory
	TVector3 pos{999,999,999};  // [cm]
	double time=-999;  // [ns]
	// should provide setters/getters for any optional components of these (i.e., probably all of them)
	int type=-1;                   // what is this; primary/secondary/???
	int target_pdg=-1;
	std::vector<int> processes; // FIXME g3 process code? g4? string? needs to match MParticle
	
	private:
	int incident_particle_idx=-1;
	bool direct_parent;
	
	public:
	void SetIncidentParticle(int idx);
	
	// some particles may store an indirect parent index, in which case the parent pdg
	// may not be the pdg of the particle that actually created the particle.
	// in this case the below will be the actual pdg code of the direct parent
	int incident_particle_pdg;
	TVector3 incident_particle_mom{999,999,999};  // at the interaction point
	
	// contents may vary by source
	BStore extraInfo; // e.g. medium_id, kcase_code
	
	MParticle* GetIncidentParticle();
	bool IsParentDirect();
	
	DataModel* m_data=nullptr;
	
};

#endif

