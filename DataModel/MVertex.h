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
	// TODO make these arguments of the constructor so that they're mandatory
	// (since we don't have getters and can't really identify if they're valid based on default values)
	TVector3 pos{999,999,999};  // [cm]
	double time=-999;  // [ns]
	
	// provide setters/getters for optional components (i.e. probably all the rest)?
	// or are the default sufficient to identify that they haven't been set...maybe
	int type=-1; // see ReadMCParticles; TODO make this an enum class and make a map to string in Constants
	int target_pdg=-1;
	std::vector<int> processes; // g3 process code if >0, g4 process code if <0.
	                            // for now G4 codes are only used when there is no G3 equivalent (in SKG4)
	
	private:
	int incident_particle_idx=-1;
	bool direct_parent;
	
	public:
	void SetIncidentParticle(int idx);
	
	// some particles may store an indirect parent index, in which case the parent pdg
	// may not be the pdg of the particle that actually created the particle.
	// in this case the below will be the actual pdg code of the direct parent
	int incident_particle_pdg=-1;
	TVector3 incident_particle_mom{999,999,999};  // at the interaction point
	
	// contents may vary by source
	BStore extraInfo; // e.g. medium_id, kcase_code
	
	MParticle* GetIncidentParticle();
	bool IsParentDirect();
	
	DataModel* m_data=nullptr;
	void Print(bool verbose=false);
	std::string PrintProcesses();
	
};

#endif

