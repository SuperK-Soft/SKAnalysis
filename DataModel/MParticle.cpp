#include "MParticle.h"
#include "MVertex.h"
#include "DataModel.h"

#include <cstdlib> // abs

MParticle::MParticle(){
	m_data = DataModel::GetInstance();
}

namespace {
	// use TObject::fBits bit 15 to flag that the value has been set.
	const uint32_t initbit{15};
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
	start_mom = TVector3{nums};
	start_mom.SetBit(initbit);
}

void MParticle::SetEndMom(double* nums){
	end_mom = TVector3{nums};
	end_mom.SetBit(initbit);
}

void MParticle::SetStartMom(float* nums){
	start_mom = TVector3{nums};
	start_mom.SetBit(initbit);
}

void MParticle::SetEndMom(float* nums){
	end_mom = TVector3{nums};
	end_mom.SetBit(initbit);
}

void MParticle::SetStartMom(const TVector3& nums){
	start_mom = TVector3{nums};
	start_mom.SetBit(initbit);
}

void MParticle::SetEndMom(const TVector3& nums){
	end_mom = TVector3{nums};
	end_mom.SetBit(initbit);
}

void MParticle::SetStartMom(double x, double y, double z){
	start_mom = TVector3{x, y, z};
	start_mom.SetBit(initbit);
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
	return nullptr; //TVector3{999,999,999};
}

double* MParticle::GetStartTime(){
	if(start_vtx_idx>=0 && start_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(start_vtx_idx).time;
	}
	// else no associated vertex
	return nullptr; //-999;
}

TVector3* MParticle::GetStartMom(){
	if(start_mom.TestBit(initbit)){
		return &start_mom;
	}
	// else no associated vertex
	return nullptr;
}

std::vector<int>* MParticle::GetStartProcesses(){
	if(start_vtx_idx>=0 && start_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(start_vtx_idx).processes;
	}
	// else no associated vertex
	return nullptr; //std::vector<int>{};
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
	return nullptr; //TVector3{999,999,999};
}

double* MParticle::GetEndTime(){
	if(end_vtx_idx>=0 && end_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(end_vtx_idx).time;
	}
	// else no associated vertex
	return nullptr; //-999;
}

TVector3* MParticle::GetEndMom(){
	if(end_mom.TestBit(initbit)){
		return &end_mom;
	}
	// else no associated vertex
	return nullptr;
}

std::vector<int>* MParticle::GetEndProcesses(){
	if(end_vtx_idx>=0 && end_vtx_idx<m_data->eventVertices.size()){
		return &m_data->eventVertices.at(end_vtx_idx).processes;
	}
	// else no associated vertex
	return nullptr; //std::vector<int>{};
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
	return parent_idx;
}

// ===================
// printers, for convenience

std::string MParticle::PrintStartPos(){
	std::string ret;
	TVector3* sp = GetStartPos();
	//if(start_vtx_idx!=TVector3{999,999,999}){
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
		ret+="?";
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
	//if(end_vtx_idx!=TVector3{999,999,999}){
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
		ret+="?";
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

void MParticle::Print(){
	MParticle& parti = *this;
	std::cout<<"\tPDG code: "<<parti.pdg<<"\n"
	         <<"\tcreation vertex index: "<<parti.start_vtx_idx<<"\n"
	         <<"\tcreation time [ns]: "<<parti.PrintStartTime()<<"\n"
	         <<"\tcreation position [cm]: "<<parti.PrintStartPos()<<"\n"
	         <<"\ttermination vertex index: "<<parti.end_vtx_idx<<"\n"
	         <<"\ttermination time [ns]: "<<parti.PrintEndTime()<<"\n"
	         <<"\ttermination pos [cm]: "<<parti.PrintEndPos()<<"\n"
	         <<"\tinitial momentum [GeV/c]: "<<parti.PrintStartMom()<<"\n"
	         <<"\tfinal momentum [GeV/c]: "<<parti.PrintEndMom()<<"\n"
	         <<"\tparent particle index in this array: "<<parti.GetNearestParentIndex();
	if(parti.GetParent()==nullptr){
		std::cout<<" (This was a primary particle)\n";
	} else {
		if(!parti.IsParentDirect()) std::cout<<" (This is an indirect parent)";
		std::cout<<"\n\tcreated by an interaction of type(s) " << parti.PrintStartProcesses()
		         <<" between an incident particle of type ";
		if(parti.GetStartVertex()!=nullptr){
			int incident_pdg_i = parti.GetStartVertex()->incident_particle_pdg;
			std::string incident_pdg = (incident_pdg_i<0) ? "?" : toString(incident_pdg_i);
			std::cout<<incident_pdg<<" and kinetic energy "
			         <<parti.GetStartVertex()->incident_particle_mom.Mag()<<" GeV ";
		} else {
			std::cout<<"? and kinetic energy ? GeV ";
		}
		std::cout<<"with a target particle of type ";
		if(parti.GetStartVertex()!=nullptr){
			int target_pdg_i = parti.GetStartVertex()->target_pdg;
			std::string target_pdg_s = (target_pdg_i<0) ? "?" : toString(target_pdg_i);
			std::cout<<target_pdg_s;
		} else {
			std::cout<<"?";
		}
		std::cout<<" within a medium of type ";
		if(parti.GetStartVertex()->extraInfo.Has("medium_id")){
			int medium_id=-1;
			parti.GetStartVertex()->extraInfo.Get("medium_id",medium_id);
			std::cout<<((medium_id<0) ? "?" : toString(medium_id))<<"\n";
		} else {
			std::cout<<"?\n";
		}
	}
	std::cout<<"\ttermination process code list: " << parti.PrintEndProcesses()<<"\n";
}
