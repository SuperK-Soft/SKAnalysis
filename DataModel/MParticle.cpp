#include "MParticle.h"
#include "MVertex.h"
#include "DataModel.h"
#include "Constants.h"

#include <cstdlib> // abs

namespace {
	// use TObject::fBits bit 15 to flag that the value has been set.
	const uint32_t initbit{15};
}

MParticle::MParticle(){
	m_data = DataModel::GetInstance();
	extraInfo = new BStore{true,constants::BSTORE_BINARY_FORMAT};
}

MParticle::~MParticle(){
	if(extraInfo) delete extraInfo;
	extraInfo=nullptr;
}

MParticle::MParticle(MParticle&& rhs){
	
	m_data = rhs.m_data;
	pdg = rhs.pdg;
	start_vtx_idx = rhs.start_vtx_idx;
	end_vtx_idx = rhs.end_vtx_idx;
	parent_idx = rhs.parent_idx;
	direct_parent = rhs.direct_parent;
	daughters = rhs.daughters;
	
	start_mom = rhs.start_mom;
	end_mom = rhs.end_mom;
	startE = rhs.startE;
	endE = rhs.endE;
	
	extraInfo = rhs.extraInfo;
	rhs.extraInfo = nullptr;
	
}

void MParticle::SetParentIndex(int idx){
	// if -1, unrecorded.
	// if 0+, index of the direct parent.
	// If <0, absolute value - 2 is index of the nearest recorded parent.
	if(idx<-1){
		idx = std::abs(idx) - 2;
		direct_parent = false;
	} else {
		direct_parent = true;
	}
	parent_idx = idx;
}

bool MParticle::IsParentDirect(){
	return direct_parent;
}

void MParticle::SetStartMom(double* nums){
	if(nums==nullptr){
		std::cerr<<"MParticle::SetStartMom called with nullptr!"<<std::endl;
		return;
	}
	start_mom = TVector3{nums};
	start_mom.SetBit(initbit);
}

void MParticle::SetStartMom(float* nums){
	if(nums==nullptr){
		std::cerr<<"MParticle::SetStartMom called with nullptr!"<<std::endl;
		return;
	}
	start_mom = TVector3{nums};
	start_mom.SetBit(initbit);
}

void MParticle::SetStartMom(const TVector3& nums){
	start_mom = TVector3{nums};
	start_mom.SetBit(initbit);
}

void MParticle::SetStartMom(double x, double y, double z){
	start_mom = TVector3{x, y, z};
	start_mom.SetBit(initbit);
}

void MParticle::SetEndMom(double* nums){
	if(nums==nullptr){
		std::cerr<<"MParticle::SetEndMom called with nullptr!"<<std::endl;
		return;
	}
	end_mom = TVector3{nums};
	end_mom.SetBit(initbit);
}

void MParticle::SetEndMom(float* nums){
	if(nums==nullptr){
		std::cerr<<"MParticle::SetEndMom called with nullptr!"<<std::endl;
		return;
	}
	end_mom = TVector3{nums};
	end_mom.SetBit(initbit);
}

void MParticle::SetEndMom(const TVector3& nums){
	end_mom = TVector3{nums};
	end_mom.SetBit(initbit);
}

void MParticle::SetEndMom(double x, double y, double z){
	end_mom = TVector3{x, y, z};
	end_mom.SetBit(initbit);
}

// Getters ==========

MVertex* MParticle::GetStartVertex(){
	if(start_vtx_idx>=0 && start_vtx_idx < m_data->eventVertices.size()){
		return &m_data->eventVertices.at(start_vtx_idx);
	}
	// else no associated vertex
	return nullptr;
}

TVector3* MParticle::GetStartPos(){
	if(start_vtx_idx>=0 && start_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(start_vtx_idx).pos;
	}
	// else no associated vertex
	return nullptr;
}

double* MParticle::GetStartTime(){
	if(start_vtx_idx>=0 && start_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(start_vtx_idx).time;
	}
	// else no associated vertex
	return nullptr;
}

TVector3* MParticle::GetStartMom(){
	if(start_mom.TestBit(initbit)){
		return &start_mom;
	}
	// else no associated vertex
	return nullptr;
}

double* MParticle::GetStartE(){
	if(start_mom.TestBit(initbit)){
		double mass = PdgToMass(pdg);
		if(mass<0) return nullptr; // e.g. unrecognised pdg
		startE = sqrt(start_mom.Mag2() + pow(mass,2.)) - mass;
		//startE = start_mom.Mag2() / (2.*mass);
		return &startE;
	}
	// else no associated vertex
	return nullptr;
}

std::vector<int>* MParticle::GetStartProcesses(){
	if(start_vtx_idx>=0 && start_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(start_vtx_idx).processes;
	}
	// else no associated vertex
	return nullptr;
}

MVertex* MParticle::GetEndVertex(){
	if(end_vtx_idx>=0 && end_vtx_idx < m_data->eventVertices.size()){
		return &m_data->eventVertices.at(end_vtx_idx);
	}
	// else no associated vertex
	return nullptr;
}

TVector3* MParticle::GetEndPos(){
	if(end_vtx_idx>=0 && end_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(end_vtx_idx).pos;
	}
	// else no associated vertex
	return nullptr;
}

double* MParticle::GetEndTime(){
	if(end_vtx_idx>=0 && end_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(end_vtx_idx).time;
	}
	// else no associated vertex
	return nullptr;
}

TVector3* MParticle::GetEndMom(){
	if(end_mom.TestBit(initbit)){
		return &end_mom;
	}
	// else no associated vertex
	return nullptr;
}

double* MParticle::GetEndE(){
	if(end_mom.TestBit(initbit)){
		double mass = PdgToMass(pdg);
		if(mass<0) return nullptr; // e.g. unrecognised pdg
		endE = sqrt(end_mom.Mag2() + pow(mass,2.)) - mass;
		//endE = start_mom.Mag2() / (2.*mass);
		return &endE;
	}
	// else no associated vertex
	return nullptr;
}

std::vector<int>* MParticle::GetEndProcesses(){
	if(end_vtx_idx>=0 && end_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(end_vtx_idx).processes;
	}
	// else no associated vertex
	return nullptr;
}

MParticle* MParticle::GetParent(){
	// note: should also check `direct_parent` member to see if this is direct or indirect.
	if(parent_idx>=0 && parent_idx<m_data->eventParticles.size()){
		return &m_data->eventParticles.at(parent_idx);
	}
	// else no associated particle
	return nullptr;
}

int MParticle::GetNearestParentIndex(){
	// -1: unrecorded
	// 0+: index of parent. Use IsParentDirect() to see if it's a direct parent.
	// should not have any other negative numbers.
	if(parent_idx<-1){
		std::cerr<<"Error! MParticle with parent index "<<parent_idx<<": should not be <-1!"<<std::endl;
	}
	return parent_idx;
}

// ===================
// printers, for convenience

std::string MParticle::PrintStartPos(){
	std::string ret;
	TVector3* sp = GetStartPos();
	if(sp!=nullptr){
		ret = "(" + toString(sp->X())+", "+toString(sp->Y())
			             +", "+ toString(sp->Z())+")";
	} else {
		ret="(?,?,?)";
	}
	return ret;
}

std::string MParticle::PrintStartTime(){
	std::string ret;
	double* starttp = GetStartTime();
	if(starttp){
		ret=toString(*starttp);
	} else {
		ret="?";
	}
	return ret;
}

std::string MParticle::PrintStartE(){
	std::string ret;
	double* startep = GetStartE();
	if(startep){
		ret=toString(*startep);
	} else {
		ret="?";
	}
	return ret;
}

std::string MParticle::PrintStartMom(){
	std::string ret;
	if(start_mom.TestBit(initbit)){
		ret = "(" + toString(start_mom.X())+", "+toString(start_mom.Y())
			             +", "+ toString(start_mom.Z())+")";
	} else {
		// not set
		ret="(?,?,?)";
	}
	return ret;
}

std::string MParticle::PrintStartProcesses(){
	// if >0, g3 process codes
	// if <0, g4 process codes
	// TODO we have both maps - make a Constants function to turn them into strings
	std::string proc_string="{";
	if(start_vtx_idx>=0 && start_vtx_idx<m_data->eventVertices.size()){
		std::vector<int>& processes = m_data->eventVertices.at(start_vtx_idx).processes;
		for(int i=0; i<processes.size(); ++i){
			if(i>0) proc_string += ", ";
			proc_string+= toString(processes.at(i));
		}
	} else {
		// else no associated vertex
		proc_string+="?";
	}
	proc_string+="}";
	return proc_string;
}

std::string MParticle::PrintEndPos(){
	std::string ret;
	TVector3* ep = GetEndPos();
	if(ep!=nullptr){
		ret = "(" + toString(ep->X())+", "+toString(ep->Y())
			             +", "+ toString(ep->Z())+")";
	} else {
		ret="(?,?,?)";
	}
	return ret;
}

std::string MParticle::PrintEndTime(){
	std::string ret;
	double* endtp = GetEndTime();
	if(endtp){
		ret=toString(*endtp);
	} else {
		ret="?";
	}
	return ret;
}

std::string MParticle::PrintEndE(){
	std::string ret;
	double* endep = GetEndE();
	if(endep){
		ret=toString(*endep);
	} else {
		ret="?";
	}
	return ret;
}

std::string MParticle::PrintEndMom(){
	std::string ret;
	if(end_mom.TestBit(initbit)){
		ret = "(" + toString(end_mom.X())+", "+toString(end_mom.Y())
			             +", "+ toString(end_mom.Z())+")";
	} else {
		// not set
		ret="(?,?,?)";
	}
	return ret;
}

std::string MParticle::PrintEndProcesses(){
	// if >0, g3 process codes
	// if <0, g4 process codes
	// TODO we have both maps - make a Constants function to turn them into strings
	std::string proc_string="{";
	if(end_vtx_idx>=0 && end_vtx_idx<m_data->eventVertices.size()){
		std::vector<int>& processes = m_data->eventVertices.at(end_vtx_idx).processes;
		for(int i=0; i<processes.size(); ++i){
			if(i>0) proc_string += ", ";
			proc_string+= toString(processes.at(i));
		}
	} else {
		// else no associated vertex
		proc_string+="?";
	}
	proc_string+="}";
	return proc_string;
}

void MParticle::Print(bool verbose){
	if(verbose){
		std::cout<<"eventParticles index: ";
		bool found=false;
		for(int i=0; i<m_data->eventParticles.size(); ++i){
			if(this==&m_data->eventParticles.at(i)){
				std::cout<<i<<"\n";
				found = true;
				break;
			}
		}
		if(!found) std::cout<<"?\n";
	}
	std::cout<<"\tPDG code: "<<pdg<<"\n"
	         <<"\tcreation vertex index: "<<start_vtx_idx<<"\n"
	         <<"\tcreation time [ns]: "<<PrintStartTime()<<"\n"
	         <<"\tcreation position [cm]: "<<PrintStartPos()<<"\n"
	         <<"\ttermination vertex index: "<<end_vtx_idx<<"\n"
	         <<"\ttermination time [ns]: "<<PrintEndTime()<<"\n"
	         <<"\ttermination pos [cm]: "<<PrintEndPos()<<"\n"
	         <<"\tinitial momentum [MeV/c]: "<<PrintStartMom()<<"\n"
	         <<"\tinitial energy [MeV]: "<<PrintStartE()<<"\n"
	         <<"\tfinal momentum [MeV/c]: "<<PrintEndMom()<<"\n"
	         <<"\tfinal energy [MeV]: "<<PrintEndE()<<"\n"
	         <<"\tparent particle index in this array: "<<GetNearestParentIndex();
	if(GetParent()==nullptr){
		std::cout<<" (This was a primary particle)\n";
	} else {
		if(!IsParentDirect()) std::cout<<" (This is an indirect parent)";
		std::cout<<"\n\tcreated by an interaction of type(s) " << PrintStartProcesses()
		         <<" between an incident particle of type ";
		if(GetStartVertex()!=nullptr){
			int incident_pdg_i = GetStartVertex()->incident_particle_pdg;
			std::string incident_pdg = (incident_pdg_i<0) ? "?" : toString(incident_pdg_i);
			std::cout<<incident_pdg<<" and kinetic energy "
			         <<GetStartVertex()->incident_particle_mom.Mag()<<" MeV ";
		} else {
			std::cout<<"? and kinetic energy ? MeV ";
		}
		std::cout<<"with a target particle of type ";
		if(GetStartVertex()!=nullptr){
			int target_pdg_i = GetStartVertex()->target_pdg;
			std::string target_pdg_s = (target_pdg_i<0) ? "?" : toString(target_pdg_i);
			std::cout<<target_pdg_s;
		} else {
			std::cout<<"?";
		}
		std::cout<<" within a medium of type ";
		if((GetStartVertex()->extraInfo) && GetStartVertex()->extraInfo->Has("medium_id")){
			int medium_id=-1;
			GetStartVertex()->extraInfo->Get("medium_id",medium_id);
			std::cout<<((medium_id<0) ? "?" : toString(medium_id))<<"\n";
		} else {
			std::cout<<"?\n";
		}
	}
	std::cout<<"\ttermination process code list: " << PrintEndProcesses()<<"\n";
	if(verbose){
		std::cout<<"Start MVertex: ";
		MVertex* vtx = GetStartVertex();
		if(vtx!=nullptr) vtx->Print(true);
		else std::cout<<"NOT FOUND\n";
		std::cout<<"End MVertex: ";
		vtx = GetEndVertex();
		if(vtx!=nullptr) vtx->Print(true);
		else std::cout<<"NOT FOUND\n";
	}
}
