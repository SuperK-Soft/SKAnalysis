#include "TrueNCaptures.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "MTreeReader.h"

TrueNCaptures::TrueNCaptures():Tool(){
	// get the name of the tool from its class name
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}


bool TrueNCaptures::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("plotsFile",plotsfile);
	
	if(plotsfile!="") MakePlots(0);
	
	 return true;
}


bool TrueNCaptures::Execute(){
	
	Log(m_unique_name+" executing...",v_debug,verbosity);
	m_data->NCapturesTrue.clear();
	
	// loop over true particles and pull out neutron capture events
	Log(m_unique_name+" scanning "+toString(m_data->eventParticles.size())
	    +" recorded particles for neutrons",v_debug,verbosity);
	for(int i=0; i<m_data->eventParticles.size(); ++i){
		MParticle& aparticle = m_data->eventParticles.at(i);
		//std::cout<<"particle "<<i<<" has pdg "<<aparticle.pdg<<std::endl;
		// check if it's a neutron
		if(aparticle.pdg==2112){
			//std::cout<<"neutron at index "<<i<<std::endl;
			bool captured=false;
			// get termination processes
			if(aparticle.GetEndProcesses()){
				// check if it captured
				for(int j=0; j<aparticle.GetEndProcesses()->size(); ++j){
					int process_code = aparticle.GetEndProcesses()->at(j);
					if(process_code==18){
						captured=true;
						break;
					}
				}
			} else {
				// uhh, we have no process codes for its termination.
				// this can happen e.g. with SECONDARY c-style arrays where the event
				// had no recorded secondaries. In this case the neutron probably
				// did not capture, hence nothing recorded; such instances have been confirmed
				// by Secondary vectors, where the neutron termination process was 'transportation'
				Log(m_unique_name+" Warning! Found neutron but no termination processes recorded!"
				   " Assuming it did not capture...",v_debug,verbosity);
				/*
				// n.b. the event number printed by this might be off by 1 because we may not
				// be using the right TreeReader (might be one that comes after this Tool)
				std::cerr<<"Dumping info for event "
				         <<m_data->Trees.begin()->second->GetEntryNumber()<<std::endl;
				for(int i=0; i<m_data->eventParticles.size(); ++i){
					m_data->eventParticles.at(i).Print(true);
				}
				m_data->vars.Set("StopLoop",1);
				*/
			}
			if(!captured){
				Log(m_unique_name+" found true neutron, but it did not capture",v_debug,verbosity);
				continue;
			}
			Log(m_unique_name+" found true neutron capture at index "+toString(i),v_debug,verbosity);
			NCapture acapture;
			get_ok = acapture.SetNeutronIndex(i);
			if(!get_ok){
				Log(m_unique_name+" Error! Captured neutron at index "+toString(i)
				     +" rejected by SetNeutronIndex?!",v_error,verbosity);
				return false;
				// this should never happen
			}
			
			// try to find the associated true IBD positron, if applicable.
			// we want to find a positron that shares the same parent index as our neutron
			// (which may be -1 if it was a simulated IBD event and both are primaries)
			if(aparticle.GetNearestParentIndex()==-1 || aparticle.IsParentDirect()){
				// primary, or valid recorded parent
				int NeutronParentIdx = aparticle.GetNearestParentIndex();
				int PositronIndex =-1;
				// search for a positron with the same parent index
				for(int j=0; j<m_data->eventParticles.size(); ++j){
					MParticle& bparticle = m_data->eventParticles.at(j);
					if(bparticle.pdg==-11 && (bparticle.IsParentDirect() || NeutronParentIdx==-1)
					    && bparticle.GetNearestParentIndex()==NeutronParentIdx){
						if(PositronIndex<0){
							PositronIndex=j;
							// if we don't care about the possibility of finding >1 matching positron:
							//break;
						} else {
							Log(m_unique_name+" Found multiple positrons from the same parent index"
							   " as the neutron ("+toString(NeutronParentIdx)+")",v_error,verbosity);
							// we can't uniquely identify the correct positron... should be unlikely...
						}
					}
				}
				acapture.SetIBDPositronIndex(PositronIndex);
			} // else this neutron has no direct recorded parent, so no way to match.
			
			m_data->NCapturesTrue.push_back(acapture);
		}
	}
	
	if(verbosity>2) PrintCaptures();
	if(plotsfile!="") MakePlots(1);
	
	return true;
}

bool TrueNCaptures::Finalise(){
	
	if(plotsfile!="") MakePlots(2);
	
	 return true;
}

bool TrueNCaptures::PrintCaptures(){
	std::cout<<"This event contained "<<m_data->NCapturesTrue.size()
	         <<" recorded true neutron captures"
	         <<"\n==========================================\n";
	for(int i=0; i<m_data->NCapturesTrue.size(); ++i){
		if(i>0) std::cout<<"------------------------------------------\n";
		std::cout<<"Capture "<<i<<"\n";
		NCapture& acapture = m_data->NCapturesTrue.at(i);
		acapture.Print((verbosity>5));
	}
	std::cout<<"=========================================="<<std::endl;
	return true;
}

bool TrueNCaptures::MakePlots(int step){
	// pi chart of capture nuclide
	// following: stack by nuclide, gd by isotopes then H
	// distribution of neutron travel time
	// distribution of neutron travel distance
	// distribution of number of daughter gammas
	// distribution of total gamma enery
	// distribution of individual gamma energy
	// distribution of num conversion electrons
	// distribution of total electron energy
	// distribution of individual electron energy
	// remember, always make a TTree not Histograms!
	
	if(step==0){
		
		// initialisation; get or make file & tree
		fplots = m_data->OpenFileForWriting(plotsfile);
		if(fplots==nullptr || fplots->IsZombie()) return false;
		fplots->cd();
		tplots = new TTree("eventtree","True Neutron Capture Variables");
		std::vector<std::string> dvars{"prompt_x","prompt_y","prompt_z","prompt_t",
		                               "capt_x","capt_y","capt_z","capt_t",
		                               /*"prompt_dwall", "prompt_deffwall", ... other? same for capt..*/
		                               "neutron_travel_time","neutron_travel_dist",
		                               "xtravel","ytravel","ztravel","neutron_start_energy",
		                               "neutron_tot_gammaE","neutron_tot_electronE","neutron_tot_daughterE"};
		// not sure whether the entries in a std::map may be moved around when new entries are added
		// to be sure, make all the entries first, then set the addresses after
		for(auto&& avar : dvars){
			dbranchvars[avar] = 0;
		}
		for(auto&& avar : dvars){
			tplots->Branch(avar.c_str(), &dbranchvars[avar]);
		}
		std::vector<std::string> ivars{"nuclide_daughter_pdg",
		                               "neutron_n_gammas","neutron_n_electrons","neutron_n_daughters"};
		for(auto&& avar : ivars){
			ibranchvars[avar] = 0;
		}
		for(auto&& avar : ivars){
			tplots->Branch(avar.c_str(), &ibranchvars[avar]);
		}
		std::vector<std::string> vvars{"gamma_energy", "electron_energy","gamma_time","electron_time"};
		for(auto&& avar : vvars){
			vbranchvars[avar] = std::vector<double>{};
		}
		for(auto&& avar : vvars){
			tplots->Branch(avar.c_str(), &vbranchvars[avar]);
		}
		
	} else if(step==1){
		
		Log(m_unique_name+" filling plots tree with "+toString(m_data->NCapturesTrue.size())+" true captures",v_debug,verbosity);
		// execution - fill tree
		for(NCapture& acap : m_data->NCapturesTrue){
			double tmpd=0;
			ibranchvars.at("nuclide_daughter_pdg") = (acap.GetDaughterNuclide()) ? acap.GetDaughterNuclide()->pdg : 0;
			dbranchvars.at("neutron_travel_time") = (acap.NeutronTravelTime(tmpd)) ? tmpd : -1.;
			dbranchvars.at("neutron_travel_dist") = (acap.NeutronTravelDist(tmpd)) ? tmpd : -1.;
			double* startE = acap.GetNeutron() ? acap.GetNeutron()->GetStartE() : nullptr;
			dbranchvars.at("neutron_start_energy") = (startE) ? *startE : 0.;
			TVector3* cappos = acap.GetPos();
			if(cappos){
				dbranchvars.at("capt_x") = cappos->X();
				dbranchvars.at("capt_y") = cappos->Y();
				dbranchvars.at("capt_z") = cappos->Z();
			} else {
				dbranchvars.at("capt_x") = 9999;
				dbranchvars.at("capt_y") = 9999;
				dbranchvars.at("capt_z") = 9999;
			}
			double* cap_t = acap.GetTime();
			dbranchvars.at("capt_t") = (cap_t ? *cap_t : 9999);
			TVector3* startpos = (acap.GetNeutron()) ? acap.GetNeutron()->GetStartPos() : nullptr;
			if(startpos){
				dbranchvars.at("prompt_x") = startpos->X();
				dbranchvars.at("prompt_y") = startpos->Y();
				dbranchvars.at("prompt_z") = startpos->Z();
			} else {
				dbranchvars.at("prompt_x") = 9999;
				dbranchvars.at("prompt_y") = 9999;
				dbranchvars.at("prompt_z") = 9999;
			}
			double* prompt_t = (acap.GetNeutron()) ? acap.GetNeutron()->GetStartTime() : nullptr;
			dbranchvars.at("prompt_t") = (prompt_t ? *prompt_t : 9999);
			double x=0,y=0,z=0;
			if(cappos && startpos){
				x = startpos->X() - cappos->X();
				y = startpos->Y() - cappos->Y();
				z = startpos->Z() - cappos->Z();
			}
			dbranchvars.at("xtravel") = x;
			dbranchvars.at("ytravel") = y;
			dbranchvars.at("ztravel") = z;
			int tmpi=0;
			ibranchvars.at("neutron_n_gammas") = (acap.NGammas(tmpi) ? tmpi : -1);
			ibranchvars.at("neutron_n_daughters") = tmpi;
			tmpi=0;
			ibranchvars.at("neutron_n_electrons") = (acap.NConversiones(tmpi) ? tmpi : -1);
			ibranchvars.at("neutron_n_daughters") += tmpi;
			dbranchvars.at("neutron_tot_daughterE") = 0;
			double sumgammae=0, sumconvee=0;
			dbranchvars.at("neutron_tot_gammaE") = acap.SumGammaE(sumgammae) ? sumgammae : 0;
			dbranchvars.at("neutron_tot_electronE") = acap.SumConversioneE(sumconvee) ? sumconvee : 0;
			dbranchvars.at("neutron_tot_daughterE") = sumgammae + sumconvee;
			vbranchvars.at("gamma_energy").clear();
			vbranchvars.at("gamma_time").clear();
			vbranchvars.at("electron_energy").clear();
			vbranchvars.at("electron_time").clear();
			sumgammae=0; // debug...
			std::vector<int> daughters;
			if(acap.GetDaughters(daughters)){
				for(int daughteridx : daughters){
					if(daughteridx<0 || daughteridx>=m_data->eventParticles.size()){
						Log(m_unique_name+" Error! bad daughter idx: "+toString(daughteridx),v_error,verbosity);
						// this event number here may be off by one, because the first TreeReader
						// might be one relating to a later Tool. Meh, close enough.
						std::cerr<<"Dumping event "
						         <<m_data->Trees.begin()->second->GetEntryNumber()
						         <<" info"<<std::endl;
						for(int i=0; i<m_data->eventParticles.size(); ++i){
							m_data->eventParticles.at(i).Print(true);
						}
						m_data->vars.Set("StopLoop",1);
						break;
						continue;
					}
					MParticle* adaughter = &m_data->eventParticles.at(daughteridx);
					startE = adaughter->GetStartE();
					double* startT = adaughter->GetStartTime();
					if(startE && startT){
						if(adaughter->pdg==22){
							if(*startE>10){
								Log(m_unique_name+" Error! Gamma from ncapture with energy "
								    +toString(*startE)+"> 10MeV?!",v_error,verbosity);
								// this event number here may be off by one, because the first TreeReader
								// might be one relating to a later Tool. Meh, close enough.
								std::cerr<<"Dumping event "
										 <<m_data->Trees.begin()->second->GetEntryNumber()
										 <<" info"<<std::endl;
								for(int i=0; i<m_data->eventParticles.size(); ++i){
									m_data->eventParticles.at(i).Print(true);
								}
								m_data->vars.Set("StopLoop",1);
								break;
							}
							vbranchvars.at("gamma_energy").push_back(*startE);
							vbranchvars.at("gamma_time").push_back(*startT);
							sumgammae += *startE;
						} else if(adaughter->pdg==11){
							vbranchvars.at("electron_energy").push_back(*startE);
							vbranchvars.at("electron_time").push_back(*startT);
						}
					}
				}
				if(sumgammae!=dbranchvars.at("neutron_tot_gammaE")){
					std::cerr<<"sumgammae ("<<sumgammae<<") != that from neutron ("
					         <<dbranchvars.at("neutron_tot_gammaE")<<")"<<std::endl;
				}
			}
			tplots->Fill();
			//tplots->Show(tplots->GetEntries()-1);
		}
		Log(m_unique_name+": done filling tree",v_debug,verbosity);
		
	} else if(step==2){
		
		// finalise - write out file
		fplots->Write(nullptr,TObject::kOverwrite);
		if(tplots) tplots->ResetBranchAddresses();
		m_data->CloseFile(plotsfile);
		
	}
	
	return true;
}

