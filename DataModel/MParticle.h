#ifndef MPARTICLE_H
#define MPARTICLE_H
#include "TVector3.h"
#include "BStore.h"
#include <vector>

class MVertex;
class DataModel;

class MParticle {
	public:
	MParticle();
	~MParticle();
	MParticle(MParticle&) = delete; // must manually declare copy construction that deals with BStore
	MParticle(MParticle&&);
	int pdg=-1; // pdg code
	int start_vtx_idx=-1;
	int end_vtx_idx=-1;
	
	private:
	int parent_idx=-1; // index in eventParticles vector, kinda.
	// if -1, unrecorded.
	// if 0+, index of the direct parent.
	// If <0, absolute value - 2 is index of the nearest recorded parent.
	bool direct_parent;   // whether the parent index represents an indirect grandparent (or above)
	public:
	void SetParentIndex(int idx);
	std::vector<int> daughters;
	// contents may vary by source
	BStore* extraInfo=nullptr; // e.g. ceren_flag, end_flag
	
	// don't like this solution, but to know whether a particle has its
	// start or end momentum set (as not all sources record them)
	// require the use of setters, and provide getters
	// FIXME perhaps we could change this to use std::optional
	private:
	TVector3 start_mom{0,0,0};  // [MeV/c]
	TVector3 end_mom{0,0,0};    // [MeV/c]
	double startE;              // [MeV]
	double endE;                // [MeV]
	public:
	void SetStartMom(double* nums);
	void SetEndMom(double* nums);
	void SetStartMom(float* nums);
	void SetEndMom(float* nums);
	void SetStartMom(const TVector3& nums);
	void SetEndMom(const TVector3& nums);
	void SetStartMom(double x, double y, double z);
	void SetEndMom(double x, double y, double z);
	
	double* GetStartTime();
	TVector3* GetStartPos();
	TVector3* GetStartMom();
	double* GetStartE();
	double* GetEndTime();
	TVector3* GetEndPos();
	TVector3* GetEndMom();
	double* GetEndE();
	std::vector<int>* GetStartProcesses();
	std::vector<int>* GetEndProcesses();
	bool IsParentDirect();
	
	MVertex* GetStartVertex();
	MVertex* GetEndVertex();
	MParticle* GetParent();
	int GetNearestParentIndex();  // should give a valid index or -1
	
	std::string PrintStartPos();
	std::string PrintStartTime();
	std::string PrintStartMom();
	std::string PrintStartE();
	std::string PrintStartProcesses();
	std::string PrintEndPos();
	std::string PrintEndTime();
	std::string PrintEndMom();
	std::string PrintEndE();
	std::string PrintEndProcesses();
	void Print(bool verbose=false);
	
	DataModel* m_data=nullptr;
	
};
#endif

