#include "NCapture.h"
#include "DataModel.h"
#include "MParticle.h"

NCapture::NCapture(){
	m_data = DataModel::GetInstance();
}

bool NCapture::SetNeutronIndex(int neutron_index){
	if(neutron_index>=0 && neutron_index<m_data->eventParticles.size()
	  && m_data->eventParticles.at(neutron_index).pdg==2112){
		neutron_idx = neutron_index;
	} else {
		// out of bounds, or not a neutron
		m_data->Log->Log("NCapture::SetNeutronIndex Error! eventParticles index given "
		    "is not that of a neutron!",0,0);
		return false;
	}
	
	got_capture_t=false;
	got_daughters = false;
	got_daughter_nuclide = false;
	return true;
}

bool NCapture::SetIBDPositronIndex(int positron_index){
	if(positron_index>=0 && positron_index<m_data->eventParticles.size()
	  && m_data->eventParticles.at(positron_index).pdg==-11){
		positron_idx = positron_index;
	} else {
		// out of bounds, or not a positron
		m_data->Log->Log("NCapture::SetIBDPositronIndex Error! eventParticles index given "
		    "is not that of a positron!",0,0);
		return false;
	}
	return true;
}

MParticle* NCapture::GetNeutron(){
	if(neutron_idx>=0 && neutron_idx<m_data->eventParticles.size()){
		return &m_data->eventParticles.at(neutron_idx);
	}
	return nullptr;
}

MParticle* NCapture::GetIBDPositron(){
	if(positron_idx>=0 && positron_idx<m_data->eventParticles.size()){
		return &m_data->eventParticles.at(positron_idx);
	}
	return nullptr;
}

MParticle* NCapture::GetDaughterNuclide(){
	// if we've already found it, return it
	if(got_daughter_nuclide){
		if(daughter_nuclide_idx>=0) return &m_data->eventParticles.at(daughter_nuclide_idx);
		else return nullptr;
	}
	daughter_nuclide_idx = -1; // shouldn't be necessary
	got_daughter_nuclide=true;
	
	// if we have no daughters, we can't find the nuclide
	if(!GetDaughters()) return nullptr;
	
	// scan daughters for nuclide
	for(int i=0; i<daughters.size(); ++i){
		int daughter_idx = daughters.at(i);
		// pdg nuclear codes are 10-digit numbers +-10LZZZAAAI
		// BUT skdetsim returns some "custom" pdg codes for nuclides in the range 100,000+
		if(daughter_idx>=0 && daughter_idx < m_data->eventParticles.size() // do we need this?
		   && m_data->eventParticles.at(daughter_idx).pdg>1E5){
			if(daughter_nuclide_idx>=0){
				/*
				m_data->Log->Log("NCapture::GetDaughterNuclide Error! Found multiple daughter nuclei "
				    "from neutron capture!",0,0);
				std::cerr<<"this neutron had "<<daughters.size()<<" daughters, of types:\n";
				for(int j=0; j<daughters.size(); ++j){
					std::cout<<"\tdaughter "<<j<<", index: "<<daughters.at(j)
					         <<", type: "<<PdgToString(m_data->eventParticles.at(daughters.at(j)).pdg)<<"\n";
				}
				std::cout<<std::endl;
				*/
				
				// we do sometimes have events that have multiple nuclide daughters from the neutron:
				// e.g. deuterium + C12, or deuterium + O16
				// how can we pick the right nuclide? well, first let's prioritise the presence of D2 or Gd
				int this_daughter_pdg = m_data->eventParticles.at(daughter_idx).pdg;
				int old_daughter_pdg = m_data->eventParticles.at(daughter_nuclide_idx).pdg;
				bool this_daughter_likely=false;
				bool old_daughter_likely=false;
				if(this_daughter_pdg==100045 || this_daughter_pdg==1000641560 || this_daughter_pdg==1000641580){
					this_daughter_likely=true;
				}
				if(old_daughter_pdg==100045 || old_daughter_pdg==1000641560 || old_daughter_pdg==1000641580){
					old_daughter_likely=true;
				}
				if(this_daughter_likely != old_daughter_likely){  // one and only one is likely
					if(this_daughter_likely) daughter_nuclide_idx = daughter_idx;
					continue;
				}
				else if(this_daughter_likely){ // both are likely... an unlikely situation
					// turns out we do seem to have this  - occasionally we have events with more than one Gd nuclei!
					// guess we'll just have to accept that. Is it a bug in SKG4? Some unaccounted mechanism? who knows.
					m_data->Log->Log("NCapture::GetDaughterNuclide Error! Found multiple daughter nuclei "
					                 "from neutron capture!",0,0);
					std::cerr<<"this neutron had "<<daughters.size()<<" daughters, of types:\n";
					for(int j=0; j<daughters.size(); ++j){
						std::cerr<<"\tdaughter "<<j<<", index: "<<daughters.at(j)
						         <<", type: "<<PdgToString(m_data->eventParticles.at(daughters.at(j)).pdg)<<"\n";
					}
					std::cerr<<std::endl;
					// i guess we just abort the search and take the first one we found...?
					// alternatively return neither?
					//return nullpr;
					break;
				}
				// else neither are likely. i guess retain the first one we found, but keep looking
				continue;
			}
			daughter_nuclide_idx = daughter_idx;
		}
	}
	if(daughter_nuclide_idx<0) return nullptr;
	return &m_data->eventParticles.at(daughter_nuclide_idx);
}

bool NCapture::GetDaughters(){
	// if we've already done this, no need to do it again
	if(got_daughters) return true;
	got_daughters = true;
	// check we have corresponding MC info
	if(m_data->eventParticles.size()==0) return false;
	// check we have a valid index
	if(neutron_idx<0 || neutron_idx>=m_data->eventParticles.size()) return false;
	daughters = m_data->eventParticles.at(neutron_idx).daughters;
	return true;
}

bool NCapture::GetDaughters(std::vector<int>& daughtersin){
	daughtersin.clear();
	if(!GetDaughters()) return false;
	daughtersin=daughters;
	return true;
}


bool NCapture::NGammas(int& ngammas){
	ngammas=0;
	if(!GetDaughters()) return false;
	for(int i=0; i<daughters.size(); ++i){
		int daughter_idx = daughters.at(i);
		if(daughter_idx<0 || daughter_idx>=m_data->eventParticles.size()) continue;
		if(m_data->eventParticles.at(daughter_idx).pdg==22){
			++ngammas;
		}
	}
	return true;
}

bool NCapture::NConversiones(int& nconversiones){
	nconversiones=0;
	if(!GetDaughters()) return false;
	for(int i=0; i<daughters.size(); ++i){
		int daughter_idx = daughters.at(i);
		if(daughter_idx<0 || daughter_idx>=m_data->eventParticles.size()) continue;
		if(m_data->eventParticles.at(daughter_idx).pdg==11){
			++nconversiones;
		}
	}
	return true;
}

// this may miss energy from gammas if their start momentum (energy) is not recorded
bool NCapture::SumGammaE(double& sumgammae){
	sumgammae=0;
	if(!GetDaughters()) return false;
	for(int i=0; i<daughters.size(); ++i){
		int daughter_idx = daughters.at(i);
		if(daughter_idx<0 || daughter_idx>=m_data->eventParticles.size()) continue;
		if(m_data->eventParticles.at(daughter_idx).pdg==22){
			if(m_data->eventParticles.at(daughter_idx).GetStartMom()){
				sumgammae+=m_data->eventParticles.at(daughter_idx).GetStartMom()->Mag();
			}
		}
	}
	return true;
}

// this may miss energy from electrons if their start momentum (energy) is not recorded
bool NCapture::SumConversioneE(double& sumconvee){
	sumconvee=0;
	if(!GetDaughters()) return false;
	for(int i=0; i<daughters.size(); ++i){
		int daughter_idx = daughters.at(i);
		if(daughter_idx<0 || daughter_idx>=m_data->eventParticles.size()) continue;
		if(m_data->eventParticles.at(daughter_idx).pdg==11){
			TVector3* electronmom = m_data->eventParticles.at(daughter_idx).GetStartMom();
			if(electronmom){
				double electronmass = m_data->pdgdb->GetParticle(11)->Mass();
				double electronE = sqrt(electronmom->Mag2() + pow(electronmass,2.)) - electronmass;
				sumconvee += electronE;
			}
		}
	}
	return true;
}

double* NCapture::GetTime(){
	if(got_capture_t) return &capture_time;
	got_capture_t = true;
	MParticle* neutron = GetNeutron();
	if(neutron==nullptr) return nullptr;
	capture_time = *neutron->GetEndTime();
	return &capture_time;
}

TVector3* NCapture::GetPos(){
	MParticle* neutron = GetNeutron();
	if(neutron==nullptr) return nullptr;
	return neutron->GetEndPos();
}

bool NCapture::NeutronTravelDist(double& ntraveldist){
	ntraveldist=0;
	MParticle* neutron = GetNeutron();
	if(neutron==nullptr) return false;
	TVector3* startpos = neutron->GetStartPos();
	TVector3* endpos = neutron->GetEndPos();
	if(startpos!=nullptr && endpos!=nullptr){
		ntraveldist = (*startpos-*endpos).Mag();
	} else {
		return false;
	}
	return true;
}

bool NCapture::NeutronTravelTime(double& ntravelt){
	ntravelt=0;
	MParticle* neutron = GetNeutron();
	if(neutron==nullptr) return false;
	double* startt = neutron->GetStartTime();
	double* endt = neutron->GetEndTime();
	if(startt!=nullptr && endt!=nullptr){
		ntravelt = (*endt-*startt);
	} else {
		return false;
	}
	return true;
}

void NCapture::Print(bool verbose){
	std::cout<<"\tcapture time [us]: "<<(GetTime() ? toString(*GetTime()/1000.) : "?")<<"\n"
	         <<"\tcapture position [cm]: "<<(GetPos() ? toString(*GetPos()) : "(?,?,?)")<<"\n"
	         <<"\tdaughter nuclide pdg: "<<(GetDaughterNuclide() ? toString(GetDaughterNuclide()->pdg) : "?")<<"\n";
	int ng, ne;
	double sge, see, ntd, ntt;
	std::cout<<"\tnum gammas: "<<(NGammas(ng) ? toString(ng) : "?")<<"\n"
	         <<"\ttotal gamma E [MeV]: "<<(SumGammaE(sge) ? toString(sge) : "?")<<"\n"
	         <<"\tnum conversion electrons: "<<(NConversiones(ne) ? toString(ne) : "?")<<"\n"
	         <<"\ttotal conversion election E [MeV]: "<<(SumConversioneE(see) ? toString(see) : "?")<<"\n"
	         <<"\tneutron travel distance [cm]: "<<(NeutronTravelDist(ntd) ? toString(ntd) : "?")<<"\n"
	         <<"\tneutron travel time [us]: "<<(NeutronTravelTime(ntt) ? toString(ntt/1000.) : "?")<<"\n";
	if(verbose){
		std::cout<<"Neutron MParticle: ";
		MParticle* neutron = GetNeutron();
		if(neutron) neutron->Print(true);
		else std::cout<<"NOT FOUND\n";
	}
}
