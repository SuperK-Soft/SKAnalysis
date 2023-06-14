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
	double* GetTime();
	TVector3* GetPos();
	bool NGammas(int& ngammas);
	bool NConversiones(int& nconve);
	bool SumGammaE(double& sumgammae);
	bool SumConversioneE(double& sumconvee);
	bool NeutronTravelDist(double& ntraveldist);
	bool NeutronTravelTime(double& ntravelt);
	bool GetDaughters(); // sets internal variable only
	bool GetDaughters(std::vector<int>& daughters);
	MParticle* GetTrueNeutron();
	MParticle* GetDaughterNuclide();
	void Print();
	
	private:
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
