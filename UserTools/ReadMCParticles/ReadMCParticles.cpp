/* vim:set noexpandtab tabstop=4 wrap */
#include "ReadMCParticles.h"

#include "MTreeReader.h"

ReadMCParticles::ReadMCParticles():Tool(){
	// get the name of the tool from its class name
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}

namespace {
	// used for removing redundant particles/vertices
	const double TIME_TOLERANCE=1;    // 1ns?
	const double POS_TOLERANCE=0.1;   // 1mm?
	const double MC_TIME_OFFSET=1000; // SKG4 adds 1,000 to all hit times, so to put MC particles in the same time frame, add this
}

bool ReadMCParticles::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	
	Log(m_unique_name+": Initializing",v_debug,verbosity);
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);
	std::string treeReaderName;
	m_variables.Get("TreeReaderName",treeReaderName);
	m_variables.Get("dataSrc",dataSrc);
	m_variables.Get("debugEntryNum",debugEntryNum);
	
	// Get the TreeReader
	// ------------------
	if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+": Failed to find TreeReader "+treeReaderName+" in DataModel!",0,0);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	// note in CStore where the MCTruth info came from
	std::string mcparticlesfile = myTreeReader->GetFile()->GetName();
	m_data->CStore.Set("mcparticlesfile",mcparticlesfile);
	
	return true;
}


bool ReadMCParticles::Execute(){
	
	Log(m_unique_name+": Executing event "+toString(myTreeReader->GetEntryNumber()),v_debug,verbosity);
	
	m_data->eventVertices.clear();
	m_data->eventParticles.clear();
	
	if(debugEntryNum>=0 && myTreeReader->GetEntryNumber()==debugEntryNum) verbosity=99;
	
	// N.B. times printed by the Print functions here are before adding MC_TIME_OFFSET
	if(dataSrc==0){
		if(verbosity>5) PrintSecondaryInfo();
		GetSecondaryInfo();
	} else if(dataSrc==1){
		if(verbosity>5) PrintSecondaryVectors(true);
		GetSecondaryVectors();
	}
	
	if(verbosity>2) PrintEvent();
	
	// MC adds 1,000 to all hit times, which means hits and reconstructed variables are all offset by 1,000
	// wrt to MC truth variables (i.e. particle times). Add that offset now so we align them.
	for(auto&& avertex : m_data->eventVertices){
		avertex.time += MC_TIME_OFFSET;
	}
	
	//if(debugEntryNum>=0 && myTreeReader->GetEntryNumber()==debugEntryNum) m_data->vars.Set("StopLoop",1);
	
	return true;
}

bool ReadMCParticles::PrintEvent(){
	
	std::cout<<"This event contained "<<m_data->eventParticles.size()<<" recorded particles and "
	          <<m_data->eventVertices.size()<<" recorded vertices:\n"
	          <<"==========================================\n";
	for(int pi=0; pi<m_data->eventParticles.size(); ++pi){
		if(pi!=0) std::cout<<"------------------------------------------\n";
		std::cout<<"Particle "<<pi<<"\n";
		MParticle& parti = m_data->eventParticles.at(pi);
		parti.Print();
	}
	std::cout<<"=========================================="<<std::endl;
	
	return true;
}

bool ReadMCParticles::Finalise(){
	
	return true;
}

// ==================================================
// Getters for various sources
// ==================================================

bool ReadMCParticles::PrintSecondaryVectors(bool checkconsistency){
	
	sec_info = nullptr;
	get_ok = (myTreeReader->GetBranchValue("SECONDARY",sec_info));
	if(!get_ok){
		Log(m_unique_name+": Error getting Secondary branch from tree!",v_error,verbosity);
		return false;
	}
	// Print method should have been const-qualified. Hack around it.
	//SecondaryInfo* sci = const_cast<SecondaryInfo*>(sec_info);
	//sci->Print();
	std::cout<<"\n==========================================\n";
	
	std::map<std::string, size_t> vec_sizes;
	bool vertexes_consistent=true;
	vec_sizes.emplace("vertex_pos",sec_info->vertex_pos.size());
	vec_sizes.emplace("vertex_time",sec_info->vertex_time.size());
	vec_sizes.emplace("vertex_incident_particle",sec_info->vertex_incident_particle.size());
	vec_sizes.emplace("vertex_incident_particle_pdg_code",sec_info->vertex_incident_particle_pdg_code.size());
	vec_sizes.emplace("vertex_incident_particle_momentum",sec_info->vertex_incident_particle_momentum.size());
	vec_sizes.emplace("vertex_target_pdg_code",sec_info->vertex_target_pdg_code.size());
	vec_sizes.emplace("vertex_medium_id",sec_info->vertex_medium_id.size());
	vec_sizes.emplace("vertex_process_codes",sec_info->vertex_process_codes.size());
	vec_sizes.emplace("vertex_kcase_code",sec_info->vertex_kcase_code.size());
	int asize=-1;
	for(std::pair<const std::string, size_t>& avec : vec_sizes){
		if(asize<0) asize = avec.second;
		else if(asize!=avec.second) vertexes_consistent=false;
	}
	if(vertexes_consistent){
		std::cout<<"had "<<asize<<" vertices\n";
	} else {
		Log(m_unique_name+" Error! Secondary vertex vectors are not of consistent size!",v_error,verbosity);
		for(std::pair<const std::string, size_t>& avec : vec_sizes){
			std::cout<<"size of "<<avec.first<<" = "<<avec.second<<std::endl;
		}
	}
	
	vec_sizes.clear();
	bool tracks_consistent=true;
	vec_sizes.emplace("track_creation_vtx",sec_info->track_creation_vtx.size());
	vec_sizes.emplace("track_pdg_code",sec_info->track_pdg_code.size());
	vec_sizes.emplace("track_ini_momentum",sec_info->track_ini_momentum.size());
	vec_sizes.emplace("track_parent",sec_info->track_parent.size());
	vec_sizes.emplace("track_creation_toffset",sec_info->track_creation_toffset.size());
	vec_sizes.emplace("track_termination_vtx",sec_info->track_termination_vtx.size());
	asize=-1;
	for(std::pair<const std::string, size_t>& avec : vec_sizes){
		if(asize<0) asize = avec.second;
		else if(asize!=avec.second) tracks_consistent=false;
	}
	if(tracks_consistent){
		std::cout<<"had "<<asize<<" tracks\n";
	} else {
		Log(m_unique_name+" Error! Secondary track vectors are not of consistent size!",v_error,verbosity);
		for(std::pair<const std::string, size_t>& avec : vec_sizes){
			std::cout<<"size of "<<avec.first<<" = "<<avec.second<<std::endl;
		}
	}
	
	// direct dump of vectors
	try {
		for(int i=0; i<sec_info->vertex_pos.size(); ++i){
			std::cout<<"Vertex "<<i<<"\n"
			         <<"\tpos = (";
			for(int j=0; j<sec_info->vertex_pos.at(i).size(); ++j){
				if(j>0) std::cout<<", ";
				std::cout<<sec_info->vertex_pos.at(i).at(j);
			}
			std::cout<<")\n"
			         <<"\ttime = "<<sec_info->vertex_time.at(i)<<"\n"
			         <<"\tincident_particle (index) = "<<sec_info->vertex_incident_particle.at(i)<<"\n"
			         <<"\tincident_particle_pdg_code = "
			         <<sec_info->vertex_incident_particle_pdg_code.at(i)<<"\n"
			         <<"\tincident_particle_momentum = (";
			for(int j=0; j<sec_info->vertex_incident_particle_momentum.at(i).size(); ++j){
				if(j>0) std::cout<<", ";
				std::cout<<sec_info->vertex_incident_particle_momentum.at(i).at(j);
			}
			std::cout<<")\n"
			         <<"\ttarget_pdg_code = "<<sec_info->vertex_target_pdg_code.at(i)<<"\n"
			         <<"\tmedium_id = "<<sec_info->vertex_medium_id.at(i)<<"\n"
			         <<"\tkcase_code = "<<sec_info->vertex_kcase_code.at(i)<<"\n"
			         <<"\tnum process_codes = "<<sec_info->vertex_process_codes.at(i).size()<<"\n"
			         <<"\tprocess_codes = {";
			for(int j=0; j<sec_info->vertex_process_codes.at(i).size(); ++j){
				if(j>0) std::cout<<", ";
				std::cout<<sec_info->vertex_process_codes.at(i).at(j);
			}
			std::cout<<"}\n";
		}
		
		std::cout<<"n tracks: "<<sec_info->track_creation_vtx.size()<<"\n";
		for(int i=0; i<sec_info->track_creation_vtx.size(); ++i){
			std::cout<<"Track "<<i<<"\n"
			         <<"\tpdg_code = "<<sec_info->track_pdg_code.at(i)<<"\n"
			         <<"\tini_momentum = (";
			for(int j=0; j<sec_info->track_ini_momentum.at(i).size(); ++j){
				if(j>0) std::cout<<", ";
				std::cout<<sec_info->track_ini_momentum.at(i).at(j);
			}
			std::cout<<")\n"
			         <<"\tcreation_toffset = "<<sec_info->track_creation_toffset.at(i)<<"\n"
			         <<"\tparent (index) = "<<sec_info->track_parent.at(i)<<"\n"
			         <<"\tcreation_vtx (index) = "<<sec_info->track_creation_vtx.at(i)<<"\n"
			         <<"\ttermination_vtx (index) = "<<sec_info->track_termination_vtx.at(i)<<"\n";
		}
	} catch(std::exception& e){
		Log(m_unique_name+"::PrintSecondaryVectors caught "+e.what(),v_error,verbosity);
	}
	std::cout<<std::endl;
	
	if(!checkconsistency) return true;
	bool consistent = vertexes_consistent && tracks_consistent;
	
	// consistency checks
	try {
		for(int tracki=0; tracki<sec_info->track_creation_vtx.size(); ++tracki){
			// check start vertex
			int start_vtx_idx = sec_info->track_creation_vtx.at(tracki);
			if(start_vtx_idx<0){
				// shouldn't happen: we should always have a start vertex
				Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
					+" has negative start vtx idx ("+toString(start_vtx_idx)+")!",v_error,verbosity);
					consistent=false;
			} else if(start_vtx_idx>=sec_info->vertex_time.size()){
				// references index out of bounds! uh-oh
				Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
					+" has creation vtx idx out of bounds ("+toString(start_vtx_idx)
					+"/"+toString(sec_info->vertex_time.size())+")!",v_error,verbosity);
					consistent=false;
			} else {
				// valid creation vertex.
				// check that the parent of this track matches the incident particle of the creation vertex
				if(std::abs(sec_info->track_parent.at(tracki))==
				   std::abs(sec_info->vertex_incident_particle.at(start_vtx_idx))){
					Log(m_unique_name+": secondary vectors track "+toString(tracki)+" parent track "
					   +toString(sec_info->track_parent.at(tracki))+" is consitent with incident track "
						" of creation vertex",v_debug,verbosity);
				} else {
					Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
					   +" has parent track index "+toString(sec_info->track_parent.at(tracki))
					   +" but that vertex has incident track index "
					   +toString(sec_info->vertex_incident_particle.at(start_vtx_idx)),v_error,verbosity);
					consistent=false;
				}
				// check that if this is marked primary, that the creation vertex incident_pdg_code
				// also suggests no particle (-1)
				if(sec_info->track_parent.at(tracki)==-1){
					if(sec_info->vertex_incident_particle_pdg_code.at(start_vtx_idx)!=-1){
						Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
							+" has no track parent (primary particle), but start vtx (idx "
							+toString(start_vtx_idx)+") has incident pdg code !=-1 ("
							+toString(sec_info->vertex_incident_particle_pdg_code.at(start_vtx_idx))
							+")!",v_error,verbosity);
						consistent=false;
					}
				}
			}
			
			// check end vertex
			int end_vtx_idx = sec_info->track_termination_vtx.at(tracki);
			if(end_vtx_idx<0){
				// shouldn't happen: we should always have an end vertex
				Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
					+" has negative end vtx idx ("+toString(end_vtx_idx)+")!",v_error,verbosity);
					consistent=false;
			} else if(end_vtx_idx >= sec_info->vertex_time.size()){
				// references index out of bounds! uh-oh
				Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
					+" has termination vtx idx out of bounds ("+toString(end_vtx_idx)
					+"/"+toString(sec_info->vertex_time.size())+")!",v_error,verbosity);
					consistent=false;
			} else {
				// check that the termination vertex incident particle is this one
				if(sec_info->vertex_incident_particle.at(end_vtx_idx)==tracki){
					Log(m_unique_name+": secondary vectors track "+toString(tracki)+" end vertex "
						"incident particle index is consitent with this track ",v_debug,verbosity);
				} else {
					Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
					   +" termination vertex has incident particle index "
					   +toString(sec_info->vertex_incident_particle.at(end_vtx_idx))
					   +" but this track is index "+toString(tracki),v_error,verbosity);
					consistent=false;
				}
			}
			
			// check direct parent type matches the type of incident particle at creation
			int parent_index = sec_info->track_parent.at(tracki);
			if(parent_index>0){
				// direct parent
				int parent_pdg = sec_info->track_pdg_code.at(parent_index);
				int vtx_parent_pdg = sec_info->vertex_incident_particle_pdg_code.at(start_vtx_idx);
				if(parent_pdg!=vtx_parent_pdg){
					Log(m_unique_name+": Error! secondary vectors track "+toString(tracki)
					   +" has parent type "+toString(parent_pdg)+" (from parent track "+toString(parent_index)
					   +" but creation vertex "+toString(start_vtx_idx)+" says incident particle type was "
					   +toString(vtx_parent_pdg)+"!",v_error,verbosity);
					consistent=false;
				}
			}
			
			// check that no track is its own parent
			if(sec_info->track_parent.at(tracki)==tracki){
				Log(m_unique_name+": Error; secondary vectors track "+toString(tracki)
				   +" is its own parent!",v_error,verbosity);
				consistent=false;
			}
		}
	} catch(std::exception& e){
		Log(m_unique_name+"::PrintSecondaryVectors caught "+e.what(),v_error,verbosity);
	}
	
	return consistent;
}

// Secondary vectors; populated by skdetsim
// ========================================
bool ReadMCParticles::GetSecondaryVectors(){
	// the easy one; i added this to skdetsim. maybe i'm biased.
	Log(m_unique_name+" Getting MC Particles from SECONDARY vectors",v_debug,verbosity);
	
	// SecondaryInfo vectors contain both primaries and secondaries
	sec_info = nullptr;
	get_ok = (myTreeReader->GetBranchValue("SECONDARY",sec_info));
	if(!get_ok){
		Log(m_unique_name+": Error getting Secondary branch from tree!",v_error,verbosity);
		return false;
	}
	
	// print event; should give ~same output as ReadMCParticles::PrintEvent once we're done
	if(verbosity>3){
		PrintSecondaryVectors();
	}
	
	// loop over vertices.
	for(int i=0; i<sec_info->vertex_time.size(); ++i){
		// all information about a vertex relates to the interaction point and is always valid
		// (unless N/A: e.g. incident particle, in the case of a primary vertex)
		MVertex avertex;
		avertex.pos = TVector3{sec_info->vertex_pos.at(i).data()};
		avertex.time = sec_info->vertex_time.at(i);
		avertex.type = (sec_info->vertex_incident_particle.at(i) == -1) ? 1 : 2;   // 1=primary, 2=secondary
		avertex.SetIncidentParticle(sec_info->vertex_incident_particle.at(i));
		avertex.incident_particle_pdg = sec_info->vertex_incident_particle_pdg_code.at(i);
		avertex.incident_particle_pdg = avertex.incident_particle_pdg;
		avertex.incident_particle_mom = TVector3{sec_info->vertex_incident_particle_momentum.at(i).data()};
		avertex.processes = sec_info->vertex_process_codes.at(i);
		
		// the interaction target pdg code gets set in skdetsim's MICAP algorithm, but this is (was)
		// only used for neutron interactions < 20MeV; so originally this would only have been valid
		// for low e neutrons. However, since skdetsim-gd replaced MICAP with Geant4 algorithms,
		// MICAP is no longer used at all, and i know of no suitable alternative.
		//avertex.target_pdg = sec_info->vertex_target_pdg_code.at(i);
		
		// Instead this variable got repurposed to hold a volume id, which can be mapped to
		// something like 'Inner Water volume' instead.
		// See gumed.h for numbers and sggeom_sk1.F for the corresponding names.
		avertex.extraInfo.Set("volume_id",sec_info->vertex_target_pdg_code.at(i));
		
		// we also have medium id. this can map to e.g. 'WATER' or 'ACRYLIC'
		// (i'm not sure where the corresponding mapping can be found, exactly...)
		avertex.extraInfo.Set("medium_id",sec_info->vertex_medium_id.at(i));
		// the kcase code is apparently another number describing the interaction process,
		// but again i cannot find a map of codes to their meaning in the G3 manual.
		avertex.extraInfo.Set("kcase_code",sec_info->vertex_kcase_code.at(i));
		
		m_data->eventVertices.push_back(avertex);
	}
	
	// loop over particles
	for(int i=0; i<sec_info->track_pdg_code.size(); ++i){
		MParticle aparticle;
		aparticle.pdg = sec_info->track_pdg_code.at(i);
		aparticle.start_vtx_idx = sec_info->track_creation_vtx.at(i);
		aparticle.end_vtx_idx = sec_info->track_termination_vtx.at(i);
		aparticle.SetStartMom(sec_info->track_ini_momentum.at(i).data());
		aparticle.SetEndMom(m_data->eventVertices.at(aparticle.end_vtx_idx).incident_particle_mom);
		// the track_parent index may refer to a direct parent, or if the parent wasn't stored,
		// the nearest recorded indirect parent. This is indicated by the sign of the index here,
		// which gets interpreted by the MParticle class; use MParticle::IsParentDirect() to check.
		aparticle.SetParentIndex(sec_info->track_parent.at(i));
		
		// daughter creation vertices represent the parent particle interaction that
		// resulted in their generation. But for decay products, they may be emitted some time
		// after their corresponding isotope was created (i.e. some time after the corresponding
		// parent particle interaction). The following offset accounts for any delay between
		// parent interaction and this daughter particle's emission.
		aparticle.extraInfo.Set("time_offset",sec_info->track_creation_toffset.at(i));
		
		m_data->eventParticles.push_back(aparticle);
		
		// add this as a daughter of its parent.
		// i believe we never have a situation where particle at index A has parent
		// at index B, with B > A. If that happens the parent MParticle won't yet exist.
		if(i<aparticle.GetNearestParentIndex()){
			Log(m_unique_name+" Error! Particle at index "+toString(i)+" has parent at index "
			    +toString(sec_info->track_parent.at(i))+"; unable to set this as its daughter!",
			    v_error,verbosity);
			// this is easy to fix: just add another loop after this loop creating particles
			// just for setting daughters. For now though, don't think it's necessary.
			// FIXME maybe remove this when we're more confident....
			std::cerr<<"\nDumping event "<<myTreeReader->GetEntryNumber()<<" info for debug\n";
			PrintSecondaryInfo();
			PrintSecondaryVectors(true);
			m_data->vars.Set("StopLoop",1);
			return false;
		} else {
			MParticle* parentp = aparticle.GetParent();
			if(parentp!=nullptr){
				parentp->daughters.push_back(i);
			}
		}
	}
	
	return true;
}

// Secondary and MCInfo c-style arrays
// populated by skdetsim and SKG4
// ==================================================
bool ReadMCParticles::GetSecondaryInfo(){
	Log(m_unique_name+" Getting MC Particles from SECONDARY c-style arrays",v_debug,verbosity);
	
	// MCInfo contains primary particles
	// SecondaryInfo contains secondary particles
	get_ok = (myTreeReader->GetBranchValue("SECONDARY",sec_info)) &&
	         (myTreeReader->GetBranchValue("MC",mc_info));
	if(!get_ok){
		Log(m_unique_name+": Error getting Secondary and MCInfo branches from tree!",v_error,verbosity);
		return false;
	}
	
	// debug; print MCInfo
	if(verbosity>3){
		// Print method should have been const-qualified. Hack around it.
		MCInfo* mci = const_cast<MCInfo*>(mc_info);
	}
	
	// primary particles
	int n_outgoing_primaries = mc_info->nvc;
	basic_array<int*> primary_PDG_code{intptr_t(mc_info->ipvc),n_outgoing_primaries};                // PDG code
	basic_array<float(*)[3]> primary_start_mom{intptr_t(mc_info->pvc),n_outgoing_primaries};         // [MeV?] dir & magnitude
	//basic_array<float*> primary_start_e{intptr_t(mc_info->energy),n_outgoing_primaries}            // [MeV] i believe this is incorrect, use pvc
	basic_array<int*> primary_start_vtx_id{intptr_t(mc_info->ivtivc), n_outgoing_primaries};         // id in primary vertex array, 1-based
	basic_array<int*> primary_end_vtx_id{intptr_t(mc_info->ivtfvc), n_outgoing_primaries};           // id in primary vertex array, 1-based
	// extra info?
	basic_array<float(*)[3]> primary_start_vertex{intptr_t(mc_info->pvtxvc),n_outgoing_primaries};   // [cm] - isn't this redundant with start_vtx_id?
	basic_array<int*> primary_parent_idx{intptr_t(mc_info->iorgvc),n_outgoing_primaries};            // isn't this always N/A by defn of being primary? seems to be 0 for IBD sims
	
	// sec_info also contains these arrays, which may be about primaries?
	//npvcscnd ?  number of outgoing primaries again?
	//nchildsvc[100] nchilds for primaries?
	//ichildidxv[100] ichildidx for primaries?
	
	// dunno what these mean
	basic_array<int*> primary_ceren_flag{intptr_t(mc_info->icrnvc), n_outgoing_primaries};           // ??
	basic_array<int*> primary_end_flag{intptr_t(mc_info->iflgvc), n_outgoing_primaries};             // ??
	
	// primary vertices
	int n_primary_vertices = mc_info->nvtxvc;
	basic_array<float(*)[3]> primary_vtx_pos{intptr_t(mc_info->pvtxvc),n_primary_vertices};
	basic_array<float*> primary_vtx_time{intptr_t(mc_info->timvvc),n_primary_vertices};
	basic_array<int*> primary_vtx_type{intptr_t(mc_info->iflvvc),n_primary_vertices};
	basic_array<int*> primary_vtx_parent{intptr_t(mc_info->iparvc),n_primary_vertices};
	// seems like generally primary vertex times are 0, positions are same as those in primary particles array,
	// parent indices are 0, types are 1. Basically, not much useful info here, at least for IBD sims...
	
	// secondaries info
	int n_secondaries = sec_info->nscndprt;
	basic_array<int*> secondary_PDG_code{intptr_t(sec_info->iprtscnd),n_secondaries};             // PDG code
	basic_array<float(*)[3]> secondary_start_vertex{intptr_t(sec_info->vtxscnd),n_secondaries};   // [cm?]
	basic_array<float*> secondary_start_time{intptr_t(sec_info->tscnd),n_secondaries};            // [ns]
	basic_array<float(*)[3]> secondary_start_mom{intptr_t(sec_info->pscnd),n_secondaries};        // [MeV?]
	basic_array<int*> secondary_gen_process{intptr_t(sec_info->lmecscnd),n_secondaries};          // 
	basic_array<int*> secondary_n_childs{intptr_t(sec_info->nchilds),n_secondaries};              // 
	// ^ possibly number of daughters produced in the step that produced this secondary!
	// NOT number of daughters OF the secondary itself?!?
	basic_array<int*> secondary_first_daughter_index{intptr_t(sec_info->ichildidx),n_secondaries}; // if >0, 1-based index in the secondaries array. 0 if no recorded daughters.
	basic_array<int*> secondary_parent_index{intptr_t(sec_info->iprntidx),n_secondaries};         // if >0, 1-based index in the secondaries array. If <0, 1-based index in the primaries array. if 0, unknown.
	basic_array<int*> secondary_flag{intptr_t(sec_info->iflgscnd),n_secondaries};                 // for pions; see apscndry.h for decryption
	
	// direct parent information
	basic_array<float(*)[3]> secondary_parent_mom_at_sec_creation{intptr_t(sec_info->pprnt),n_secondaries}; // 
	// the following i believe are redundant
	basic_array<int*> secondary_parent_PDG_code{intptr_t(sec_info->iprntprt),n_secondaries};                // PDG code
	basic_array<float(*)[3]> secondary_parent_init_pos{intptr_t(sec_info->vtxprnt),n_secondaries};          // [cm?] position of parent @ its birth
	basic_array<float(*)[3]> secondary_parent_init_mom{intptr_t(sec_info->pprntinit),n_secondaries};        // [MeV?] momentum of parent @ its birth
	
	// the many ways to map a secondary to its parent; who knows what they all mean?
	basic_array<int*> secondary_parent_primary_idx{intptr_t(sec_info->iprnttrk),n_secondaries};   // index of primary parent in primaries array (not populated by SKG4)
	basic_array<int*> secondary_parent_primary_pdg{intptr_t(sec_info->iorgprt),n_secondaries};    // pdg code of parent primary particle
	basic_array<int*> secondary_parent_G3_trackid{intptr_t(sec_info->itrkscnd),n_secondaries};    // parent G3 track number; any use?
	basic_array<int*> secondary_parent_G3_stackid{intptr_t(sec_info->istakscnd),n_secondaries};   // parent G3 stack number; any use?
	
	// =============
	
	// transfer information to DataModel
	
	// first primary vertices
	std::map<int,int> primary_vertex_map;   // for removing redundant copies
	Log(m_unique_name+" looping over "+toString(n_primary_vertices)+" primary vertices",v_debug,verbosity);
	for(int i=0; i<n_primary_vertices; ++i){
		// secondary information stores one vertex for each primary, but often they are all identical!
		// so first see if this there's a matching vertex already.
		Log(m_unique_name+" check for existing vtx",v_debug,verbosity);
		int vertex_idx=-1;
		for(int j=0; j<m_data->eventVertices.size(); ++j){
			if((m_data->eventVertices.at(j).time-primary_vtx_time.at(i))>TIME_TOLERANCE) continue;
			TVector3 start_pos{primary_vtx_pos.at(i).data()};
			if((m_data->eventVertices.at(j).pos-start_pos).Mag()>POS_TOLERANCE) continue;
			vertex_idx = j;
			primary_vertex_map.emplace(i,j);  // primary vertex i is now at index j
			break;
		}
		if(vertex_idx<0){
			Log(m_unique_name+" no matching index found, making a new one",v_debug,verbosity);
			MVertex avertex;
			avertex.pos = TVector3{primary_vtx_pos.at(i).data()};
			avertex.time = primary_vtx_time.at(i);
			avertex.type = primary_vtx_type.at(i); // seems to be 1 for primaries? XXX clarify!
			//std::cout<<"primary vertex "<<i<<" has type "<<avertex.type<<std::endl;
			// skdetsim/skg4 i believe set the parent to itself
			//avertex.SetIncidentParticle(primary_vtx_parent.at(i)-1); // fortran indexing
			//std::cout<<"primary vtx "<<i<<" has parent "<<primary_vtx_parent.at(i)<<std::endl;
			// instead set the parent to -1 since we use that to indicate primaries
			avertex.SetIncidentParticle(-1);
			
			// not meaningful for primaries
			//avertex.target_pdg = -1;
			//avertex.incident_particle_mom = TVector3{};
			//avertex.incident_particle_pdg = -1;
			//avertex.processes = std::vector<int>{};
			//avertex.extraInfo.Set("kcase_code",);
			//avertex.extraInfo.Set("medium_id",); 
			
			m_data->eventVertices.push_back(avertex);
		} else {
			Log(m_unique_name+" this primary vertex is the same as vertex "
			    +toString(vertex_idx)+", skipping",v_debug,verbosity);
		}
	}
	
	// then primary particles
	Log(m_unique_name+" looping over "+toString(n_outgoing_primaries)+" primary particles",v_debug,verbosity);
	for(int i=0; i<n_outgoing_primaries; ++i){
		MParticle aparticle;
		aparticle.pdg = primary_PDG_code.at(i);
		aparticle.SetStartMom(primary_start_mom.at(i).data());
		aparticle.SetParentIndex(primary_parent_idx.at(i)-1); // fortran indexing. always 0 for primaries.
		aparticle.start_vtx_idx = primary_start_vtx_id.at(i)-1; // fortran indexing
		// see if we eliminated this vertex as a redundant one
		if(primary_vertex_map.count(aparticle.start_vtx_idx)){
			aparticle.start_vtx_idx = primary_vertex_map.at(aparticle.start_vtx_idx);
		}
		aparticle.end_vtx_idx = primary_end_vtx_id.at(i)-1; // fortran indexing
		//aparticle.daughters.push_back(); // populated while scanning secondaries
		
		// not available from this source
		//aparticle.end_mom = TVector3{}; // not available
		
		// additional information
		aparticle.extraInfo.Set("ceren_flag",primary_ceren_flag.at(i));
		aparticle.extraInfo.Set("end_flag",primary_end_flag.at(i));
		
		// consistency checks
		Log(m_unique_name+" primary particle "+toString(i)+" (pdg "+toString(aparticle.pdg)+")"
		         " has primary parent idx "+toString(primary_parent_idx.at(i)-1)+
		         " and vtx id "+toString(aparticle.start_vtx_idx)+", which has parent particle "
		         +toString(primary_vtx_parent.at(aparticle.start_vtx_idx)),v_debug,verbosity);
		
		// check these vertices agree with the direct information in the primary particle arrays
		TVector3 astart{primary_start_vertex.at(i).data()};
		TVector3 otherstart{primary_vtx_pos.at(aparticle.start_vtx_idx).data()};
		if(astart!=otherstart){
			Log(m_unique_name+" Error! conflicting start vertices: primary "+toString(i)
			     +" has primary_start_vtx "+toString(astart)+", but primary vtx id "
			     +toString(aparticle.start_vtx_idx)+" which has position "+toString(otherstart)
			     ,v_error,verbosity);
		} else {
			Log(m_unique_name+" redundant primary vertex: primary_start_vertex["+toString(i)
			     +"] == primary_vtx_pos["+toString(aparticle.start_vtx_idx)+"] = "
			     +toString(astart),v_debug,verbosity);
		}
		
		m_data->eventParticles.push_back(aparticle);
	}
	
	// secondary particles
	Log(m_unique_name+" looping over "+toString(n_secondaries)+" secondaries",v_debug,verbosity);
	for(int i=0; i<n_secondaries; ++i){
		MParticle aparticle;
		//std::cout<<"init"<<std::endl;
		aparticle.pdg = secondary_PDG_code.at(i);
		aparticle.SetStartMom(secondary_start_mom.at(i).data());
		int aparent_index = secondary_parent_index.at(i)-1; // fortran indexing
		// sign indicates which set of vectors (primary or secondary) it's in.
		if(aparent_index>=(n_outgoing_primaries+n_secondaries)){
			// secondary index out of bounds!
			Log(m_unique_name+" Error! secondary parent index "+toString(aparent_index)+"/"
			   +toString(n_outgoing_primaries+n_secondaries)+" out of bounds!", v_error,verbosity);
			aparent_index=-1;
		} else if(aparent_index>=0){
			// if >= 0 this is an index in the secondaries array, but since we merged with the
			// primaries array, we need to offset by the number of primaries
			aparent_index += n_outgoing_primaries;
		} else if((std::abs(aparent_index)-1)>n_outgoing_primaries){
			// primary index out of bounds!
			Log(m_unique_name+" Error! primary parent index "+toString(std::abs(aparent_index)-2)+"/"
			    +toString(n_outgoing_primaries)+" out of bounds!", v_error,verbosity);
			aparent_index=-1;
		} else if(aparent_index<-1){
			// atmpd arrays use negative indices to indicate primary particles,
			// but we merge all particles and then use negative indices to indicate indirect parents.
			// so first make positive to indicate it's a direct parent
			aparent_index = std::abs(aparent_index);
			// and then correct for the fact that our decrement-by-one to convert fortran indices
			// to c++ indices was of the wrong sign (so would have added one, not decremented by one)
			aparent_index -= 2;
		} else if(aparent_index==-1){
			// index 0 (made -1 by our fortran indexing correction) means no recorded parent
			// maybe we have an indirect parent primary?
			Log(m_unique_name+" using primary parent index as no valid secondary parent index",
			      v_debug,verbosity);
			aparent_index = secondary_parent_primary_idx.at(i) - 1; // fortran indexing
			// if valid, indicate to MParticle class that this is an indirect parent
			if(aparent_index!=-1) aparent_index = -(aparent_index+2);
		}
		aparticle.SetParentIndex(aparent_index);
		
		// in the event that particle at index A has parent index B, where A<B,
		// our parent won't yet exist, which messes up some later code.
		// validate that this doesn't happen (pretty sure it doesn't).
		if(aparticle.GetParent()==nullptr){
			Log(m_unique_name+" Error! Parent index "+toString(aparent_index)+" not yet "
			     "in event particles!",v_error,verbosity);
			// FIXME remove if we're confident this isn't our problem
			std::cerr<<"\nDumping event "<<myTreeReader->GetEntryNumber()<<" info for debug\n";
			PrintSecondaryInfo();
			PrintSecondaryVectors(true);
			m_data->vars.Set("StopLoop",1);
			return false;
		}
		
		// not available from this source
		//aparticle.end_mom = TVector3{};
		
		// secondaries store their own vertices, so we need to make a vertex from the start/stop info.
		// But, if there are multiple recorded secondaries from one interaction, we will get duplicate
		// vertices. so first see if this there's a matching vertex already.
		Log(m_unique_name+" check for existing vtx",v_debug,verbosity);
		int vertex_idx=-1;
		for(int j=0; j<m_data->eventVertices.size(); ++j){
			if((m_data->eventVertices.at(j).time-secondary_start_time.at(i))>TIME_TOLERANCE) continue;
			TVector3 start_pos{secondary_start_vertex.at(i).data()};
			if((m_data->eventVertices.at(j).pos-start_pos).Mag()>POS_TOLERANCE) continue;
			// warning
			if(vertex_idx>0){
				Log(m_unique_name+" secondary "+toString(i)+" start vertex matches more than"
				    " one vertex in eventParticles! Reduce matching tolerance!",v_warning,verbosity);
			}
			vertex_idx = j;
			// if we don't care about the warning, we can just break after finding the first one
			break;
		}
		Log(m_unique_name+" check done, vertex_idx "+toString(vertex_idx),v_debug,verbosity);
		if(vertex_idx>0){
			// we found an existing vertex; use this one
			aparticle.start_vtx_idx = vertex_idx;
			// debug print
			if(verbosity>v_debug){
				std::cout<<"secondary "<<i
					     <<" has a start vertex at the same time and position as one already recorded: "
					     <<" this start_vtx ("<<secondary_start_time.at(i)<<", "
					     <<secondary_start_vertex.at(i).at(0)<<", "
					     <<secondary_start_vertex.at(i).at(1)<<", "
					     <<secondary_start_vertex.at(i).at(2)<<"), matched to ("
					     <<m_data->eventVertices.at(vertex_idx).time<<", "
					     <<m_data->eventVertices.at(vertex_idx).pos.X()<<", "
					     <<m_data->eventVertices.at(vertex_idx).pos.Y()<<", "
					     <<m_data->eventVertices.at(vertex_idx).pos.Z()<<")"<<std::endl;
			}
		} else {
			// no matching vertex; make one
			Log(m_unique_name+" no matching vtx, making one",v_debug,verbosity);
			MVertex start_vtx;
			start_vtx.type = 2;        // XXX clarify!
			start_vtx.pos = TVector3{secondary_start_vertex.at(i).data()};
			start_vtx.time = secondary_start_time.at(i);
			start_vtx.SetIncidentParticle(aparent_index);
			start_vtx.processes = std::vector<int>{secondary_gen_process.at(i)};
			start_vtx.incident_particle_mom = TVector3{secondary_parent_mom_at_sec_creation.at(i).data()};
			start_vtx.incident_particle_pdg = secondary_parent_PDG_code.at(i);
			
			// not available for this source
			//start_vxt.target_pdg = -1;
			
			// no extra info
			//start_vtx.extraInfo.Set();
			
			m_data->eventVertices.push_back(start_vtx);
			
			// set this new vertex as the start vertex for our current particle
			aparticle.start_vtx_idx = m_data->eventVertices.size()-1;
		}
		
		// end vertex is not stored, either as an index nor in secondary particle information.
		aparticle.end_vtx_idx = -1;  // though this is default anyway
		
		// add this particle as a daughter of its parent
		Log(m_unique_name+" daughter scan, parent_idx "+toString(aparent_index),v_debug,verbosity);
		if(aparticle.GetParent()!=nullptr){
			// parent particle is in the event particles vector
			MParticle* parent = aparticle.GetParent();
			
			// add this particle as a daughter of the parent
			parent->daughters.push_back(m_data->eventParticles.size());
			
			// santiy checks: check pdg consistency, IF the parent is direct
			if(aparticle.IsParentDirect() && parent->pdg!=secondary_parent_PDG_code.at(i)){
				Log(m_unique_name+" Warning! parent_PDG_code "+toString(secondary_parent_PDG_code.at(i))
				         +" for secondary "+toString(i)+" does not match event particle at parent index "
				         +toString(aparent_index)+" which has pdg "+toString(parent->pdg),v_debug,verbosity);
				// best i can tell this happens, and isn't our problem? seems to be indirect parents
				// e.g. reports parent pdg '22' for secondary proton indirectly from primary positron annihilation
				/*
				std::cerr<<"\nDumping event "<<myTreeReader->GetEntryNumber()<<" info for debug\n";
				PrintSecondaryInfo();
				PrintSecondaryVectors(true);
				m_data->vars.Set("StopLoop",1);
				return false;
				*/
			} else if(aparticle.IsParentDirect()) {
				Log(m_unique_name+" sanity check passed; parent_PDG_code "
				     +toString(secondary_parent_PDG_code.at(i))+" for secondary "+toString(i)
				     +" matches event particle at parent index "+toString(aparent_index),v_debug,verbosity);
			} // else not a direct parent, can't do the comparison
			
			// another: check start position consistency
			TVector3& pstart1 = m_data->eventVertices.at(parent->start_vtx_idx).pos;
			TVector3 pstart2{secondary_parent_init_pos.at(i).data()};
			if(aparticle.IsParentDirect() && (pstart1 - pstart2).Mag()>POS_TOLERANCE){
				Log(m_unique_name+" Warning! parent_init_pos for secondary "+toString(i)
				      +": "+toString(pstart2)+" does not match that of event particle at parent index "
				      +toString(aparent_index)+" which has start pos "+toString(pstart1),v_debug,verbosity);
				/*
				// best i can tell this happens... and isn't our problem? not sure why.
				std::cerr<<"\nDumping event "<<myTreeReader->GetEntryNumber()<<" info for debug\n";
				PrintSecondaryInfo();
				PrintSecondaryVectors(true);
				m_data->vars.Set("StopLoop",1);
				return false;
				*/
			} else if(aparticle.IsParentDirect()){
				Log(m_unique_name+" sanity check passed; parent_init_pos for secondary "+toString(i)
				    +" matches the start pos of event particle at parent index "
				    +toString(aparent_index),v_debug,verbosity);
			}
			// could also compare parent start momentum using secondary_parent_init_mom
			
			// one thing we could do is use a daughter's start vertex as the parent's end vertex.
			// this ASSUMES the parent terminates at the daughter start location, which may not be true!!!
			// (although it will be for neutron capture, at least...)
			if(parent->end_vtx_idx<0 && aparticle.IsParentDirect()){
				parent->end_vtx_idx = aparticle.start_vtx_idx;
				parent->SetEndMom(secondary_parent_mom_at_sec_creation.at(i).data());
			} else if(aparticle.IsParentDirect()){
				// this may happen if the parent has multiple recorded daughters
				Log(m_unique_name+" secondary "+toString(i)+" parent index "+toString(aparent_index)
				      +" already has an end vtx id "+toString(parent->end_vtx_idx)
				      +" (perhaps from another daughter?)",v_debug,verbosity);
			} // else this isn't a direct parent, so don't do this
			
		} else {
			Log(m_unique_name+" secondary "+toString(i)+" parent index out of bounds ("
			 +toString(aparent_index)+"/"+toString(m_data->eventParticles.size())+")",v_debug,verbosity);
		}
		
		// extra info
		//std::cout<<"extra info"<<std::endl;
		aparticle.extraInfo.Set("nchilds",secondary_n_childs.at(i));
		aparticle.extraInfo.Set("ichildidx",secondary_first_daughter_index.at(i)-1+n_outgoing_primaries); // correct fortran indexing and for preceding primary particles
		aparticle.extraInfo.Set("parent_primary_pdg",secondary_parent_primary_pdg.at(i));
		aparticle.extraInfo.Set("parent_primary_idx",secondary_parent_primary_idx.at(i)-1); // correct fortran indexing
		aparticle.extraInfo.Set("iflgscnd",secondary_flag.at(i));
		
		// and of course we have 3 other numbers for the parent index.
		// not even gonna store these... don't think they're useful??
		//std::cout<<"debug print"<<std::endl;
		Log(m_unique_name+" secondary "+toString(i)+" has parent idx "+toString(secondary_parent_index.at(i))
		         +", primary parent index "+toString(secondary_parent_primary_idx.at(i))
		         +", parent G3 trackid "+toString(secondary_parent_G3_trackid.at(i))
		         +", parent G3 stack "+toString(secondary_parent_G3_stackid.at(i)),v_debug,verbosity);
		
		m_data->eventParticles.push_back(aparticle);
		
	}
	
	
	return get_ok;
}

// Direct print of information in secondary c-style arrays
// ========================================================
bool ReadMCParticles::PrintSecondaryInfo(){
	
	// MCInfo contains primary particles
	// SecondaryInfo contains secondary particles
	get_ok = (myTreeReader->GetBranchValue("SECONDARY",sec_info)) &&
	         (myTreeReader->GetBranchValue("MC",mc_info));
	if(!get_ok){
		Log(m_unique_name+": Error getting Secondary and MCInfo branches from tree!",v_error,verbosity);
		return false;
	}
	
	//std::cout<<"dumping"<<std::endl;
	//mc_info->Dump();
	//sec_info->Dump();
	
	for(int i=0; i<mc_info->nvc; ++i){
		std::cout<<"primary "<<i<<"\n\t";
		std::cout<<"pdg (ipvc): "<<mc_info->ipvc[i]<<"\n\t";
		std::cout<<"vertex (pvtxvc): ("<<mc_info->pvtxvc[i][0]<<","
		                               <<mc_info->pvtxvc[i][1]<<","
		                               <<mc_info->pvtxvc[i][2]<<")\n\t";
		std::cout<<"mom (pvc): ("<<mc_info->pvc[i][0]<<","
		                         <<mc_info->pvc[i][1]<<","
		                         <<mc_info->pvc[i][2]<<")\n\t";
		std::cout<<"parent index (iorgvc): "<<mc_info->iorgvc[i]<<"\n\t";
		std::cout<<"start vtx idx (ivtivc): "<<mc_info->ivtivc[i]<<"\n\t";
		if((mc_info->ivtivc[i]-1)>=0 && (mc_info->ivtivc[i]-1)<mc_info->nvtxvc){
			std::cout<<"\tvtx: ("<<mc_info->pvtxvc[mc_info->ivtivc[i]-1][0]<<","
			                     <<mc_info->pvtxvc[mc_info->ivtivc[i]-1][1]<<","
			                     <<mc_info->pvtxvc[mc_info->ivtivc[i]-1][2]<<")\n\t";
			std::cout<<"\time: timvvc: "<<mc_info->timvvc[mc_info->ivtivc[i]-1]<<"\n\t";
			std::cout<<"\tparent (iparvc): "<<mc_info->iparvc[mc_info->ivtivc[i]-1]<<"\n\t";
			std::cout<<"\ttype (iflvvc): "<<mc_info->iflvvc[mc_info->ivtivc[i]-1]<<"\n\t";
		}
		std::cout<<"end vtx idx (ivtfvc): "<<mc_info->ivtfvc[i]<<"\n";
		if((mc_info->ivtfvc[i]-1)>=0 && (mc_info->ivtfvc[i]-1)<mc_info->nvtxvc){
			std::cout<<"\t\tvtx: ("<<mc_info->pvtxvc[mc_info->ivtfvc[i]-1][0]<<","
			                       <<mc_info->pvtxvc[mc_info->ivtfvc[i]-1][1]<<","
			                       <<mc_info->pvtxvc[mc_info->ivtfvc[i]-1][2]<<")\n\t";
			std::cout<<"\time: timvvc: "<<mc_info->timvvc[mc_info->ivtfvc[i]-1]<<"\n\t";
			std::cout<<"\tparent (iparvc): "<<mc_info->iparvc[mc_info->ivtfvc[i]-1]<<"\n\t";
			std::cout<<"\ttype (iflvvc): "<<mc_info->iflvvc[mc_info->ivtfvc[i]-1]<<"\n";
		}
	}
	std::cout<<std::endl;
	
	for(int i=0; i<sec_info->nscndprt; ++i){
		std::cout<<"secondary "<<i<<"\n\t";
		std::cout<<"pdg (iprtscnd): "<<sec_info->iprtscnd[i]<<"\n\t";
		std::cout<<"vertex (vtxscnd): ("<<sec_info->vtxscnd[i][0]<<","
		                                <<sec_info->vtxscnd[i][1]<<","
		                                <<sec_info->vtxscnd[i][2]<<")\n\t";
		std::cout<<"time (tscnd): "<<sec_info->tscnd[i]<<"\n\t";
		std::cout<<"mom (pscnd): ("<<sec_info->pscnd[i][0]<<","
		                           <<sec_info->pscnd[i][1]<<","
		                           <<sec_info->pscnd[i][2]<<")\n\t";
		std::cout<<"genproc (lmecscnd): "<<sec_info->lmecscnd[i]<<"\n\t";
		std::cout<<"parent index (iprntidx): "<<sec_info->iprntidx[i]<<"\n\t";
		std::cout<<"parent pdg (iprntprt): "<<sec_info->iprntprt[i]<<"\n\t";
		std::cout<<"parent mom at sec creation (pprnt): ("<<sec_info->pprnt[i][0]<<","
		                                                  <<sec_info->pprnt[i][1]<<","
		                                                  <<sec_info->pprnt[i][2]<<")\n\t";
		std::cout<<"parent mom at creation (pprntinit): ("<<sec_info->pprntinit[i][0]<<","
		                                                  <<sec_info->pprntinit[i][1]<<","
		                                                  <<sec_info->pprntinit[i][2]<<")\n\t";
		std::cout<<"parent pos at creation (vtxprnt): ("<<sec_info->vtxprnt[i][0]<<","
		                                                <<sec_info->vtxprnt[i][1]<<","
		                                                <<sec_info->vtxprnt[i][2]<<")\n\t";
		std::cout<<"primary parent index (iprnttrk): "<<sec_info->iprnttrk[i]<<"\n\t";
		std::cout<<"primary parent pdg (iorgprt): "<<sec_info->iorgprt[i]<<"\n\t";
		std::cout<<"parent G3 track num (itrkscnd): "<<sec_info->itrkscnd[i]<<"\n\t";
		std::cout<<"parent G3 stack num (istakscnd): "<<sec_info->istakscnd[i]<<"\n\t";
		std::cout<<"num siblings? (nchilds): "<<sec_info->nchilds[i]<<"\n\t";
		std::cout<<"daughter index (ichildidx): "<<sec_info->ichildidx[i]<<"\n\t";
		std::cout<<"pion flag (iflgscnd): "<<sec_info->iflgscnd[i]<<"\n";
	}
	return true;
}

// ATMPD c-style arrays
// =====================
bool ReadMCParticles::GetAtmpdInfo(){
	/*
	this information is in the 'h1' tree generated by converting skdetsim zbs files
	to root vfiles via:
	* 1. run skdetsim or skdetsim-gd with output option 'SKCNTL-OUTPUTTYPE 2' to generate zbs output file
	* 2. convert the zbs file to an hbk file with `fillnt_simple.sh -o (output hbook file) (input zbs file)`
	* 3. convert the hbk file to a root file with `h2root (input hbook file) (output root file)`
	The resulting files have a huge number of variables, but those relating to secondaries
	are predominantly the same as the ones in the SecondaryInfo class
	*/
	
	// numnu is 0 even when npar is >3...
	bool success = true;
	/*
//	// neutrino interaction info - first primaries array includes neutrino and target (index 0 and 1)
//	(myTreeReader.GetBranchValue("mode",nu_intx_mode))                  &&  // see neut_mode_to_string(mode)
//	(myTreeReader.GetBranchValue("numnu",tot_n_primaries))              &&  // both ingoing and outgoing
//	
//	// following are arrays of size numnu
//	(myTreeReader.GetBranchValue("ipnu",primary_pdg))                   &&  // see constants::numnu_code_to_string
//	(myTreeReader.GetBranchValue("pnu",primary_momentum))               &&  // [GeV/c]
	
	// primary event - second primaries array includes more info
	(myTreeReader.GetBranchValue("posv",primary_event_vertex))          &&  // [cm]
//	(myTreeReader.GetBranchValue("wallv",primary_event_dist_from_wall)) &&  // [cm]
	(myTreeReader.GetBranchValue("npar",n_outgoing_primaries))          &&  // should be (tot_n_primaries - 2)?
	
	// following are arrays of size npar
	(myTreeReader.GetBranchValue("ipv",primary_G3_code))                &&  // see constants::g3_to_pdg
	(myTreeReader.GetBranchValue("dirv",primary_start_mom_dir))         &&  // 
	(myTreeReader.GetBranchValue("pmomv",primary_start_mom))            &&  // [units?]
	
	// secondaries - second secondaries array...
	(myTreeReader.GetBranchValue("nscndprt",n_secondaries_1))           &&
	
	// following are arrays of size nscndprt
	(myTreeReader.GetBranchValue("iprtscnd",secondary_PDG_code_1))      &&  //
	(myTreeReader.GetBranchValue("vtxscnd",secondary_start_vertex_1))   &&  // [units?]
	(myTreeReader.GetBranchValue("tscnd",secondary_start_time_1))       &&  // [ns]? relative to event start?
	(myTreeReader.GetBranchValue("pscnd",secondary_start_mom_1))        &&  // [units?]
	(myTreeReader.GetBranchValue("lmecscnd",secondary_gen_process))     &&  // constants::G3_process_code_to_string
	(myTreeReader.GetBranchValue("nchilds",secondary_n_childs))      &&  // 
	(myTreeReader.GetBranchValue("iprntidx",parent_index))              &&  // if >0, 1-based index in this array
//	(myTreeReader.GetBranchValue("ichildidx",secondary_first_daughter_index)) &&  // if >0, 1-based index in this
	
	// further parentage information - still arrays of size nscndprt. Useful?
//	(myTreeReader.GetBranchValue("iprntprt",parent_G3_code))            &&  // or is it a PDG code?
	(myTreeReader.GetBranchValue("pprnt",parent_mom_at_sec_creation))   &&  // use w/γ to get n energy @ capture
	(myTreeReader.GetBranchValue("vtxprnt",parent_init_pos))            &&  // [cm?] parent pos @ birth
	(myTreeReader.GetBranchValue("pprntinit",parent_init_mom))          &&  // [MeV?] parent mom @ birth
//	(myTreeReader.GetBranchValue("itrkscnd",parent_G3_trackid))         &&  // how do we use this?
//	(myTreeReader.GetBranchValue("istakscnd",parent_G3_stack_trackid))  &&  // how do we use this?
	(myTreeReader.GetBranchValue("iprnttrk",parent_trackid));         //&&  // relates secondaries to primaries
	// NOTE this is carried over to daughters of secondaries, so only use as parent if iprntidx==0
//	(myTreeReader.GetBranchValue("iorgprt",parent_track_pid_code))      &&  // i'm so confused
	*/
	
	/*
	// secondaries alternative; these always seems to be empty, though.
	int n_secondaries_2;
	// the following are arrays of size npar2
	basic_array<int*> secondary_G3_code_2;
	basic_array<float(*)[3]> secondary_start_vertex_2;
	basic_array<float*> secondary_start_dist_from_wall_2;
	basic_array<float(*)[3]> secondary_start_mom_2;
	basic_array<int*> secondary_origin_2;
	(myTreeReader.GetBranchValue("npar2",n_secondaries_2))                   &&
	(myTreeReader.GetBranchValue("ipv2",secondary_G3_code_2))                &&  // 
	(myTreeReader.GetBranchValue("posv2",secondary_start_vertex_2))          &&  // [cm?] what about time?
	(myTreeReader.GetBranchValue("wallv2",secondary_start_dist_from_wall_2)) &&  // [cm?]
	(myTreeReader.GetBranchValue("pmomv2",secondary_start_mom_2))            &&  // [units?]
	(myTreeReader.GetBranchValue("iorg",secondary_origin_2))                 &&  // what is "origin"?
	*/
	
	return success;
}
