#ifndef NCAPTURE_H
#define NCAPTURE_H
#include <vector>
#include "TVector3.h"

class DataModel;
class MParticle;

class NCapture {
	public:
	NCapture();
	
	bool SetNeutronIndex(int neutron_index); // use setter to validate index
	
	//=========//
	// MC info //
	//=========//
	
	// following are calculated on demand from event particles vector
	double* GetTime();  // [us]
	TVector3* GetPos(); // [cm]
	bool NGammas(int& ngammas);
	bool NConversiones(int& nconve);
	bool SumGammaE(double& sumgammae); // [MeV]
	bool SumConversioneE(double& sumconvee);  // [MeV]
	bool NeutronTravelDist(double& ntraveldist); // [cm]
	bool NeutronTravelTime(double& ntravelt); // [us]
	bool GetDaughters(); // sets internal variable only
	bool GetDaughters(std::vector<int>& daughters);
	MParticle* GetNeutron();
	MParticle* GetDaughterNuclide();
	void Print(bool verbose=false);
	
	private:
	bool got_capture_t=false;
	double capture_time=999;
	// index in event particles vector
	int neutron_trackid=-1;
	int daughter_nuclide_idx=-1;
	// indices in event particles vector
	bool got_daughters = false;
	bool got_daughter_nuclide = false;
	std::vector<int> daughters;
	DataModel* m_data=nullptr;
};
#endif
