#include "SolarPreSelection.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "skheadC.h"
#include "ParticleCand.h"
#include <inttypes.h>
#include <algorithm>
#include <iomanip>
#include <bitset>

SolarPreSelection::SolarPreSelection():Tool(){}

bool SolarPreSelection::Initialise(std::string configfile, DataModel &data){
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("match_window", match_window);  // [s]
	match_window *= 1E9; // convert to [ns]
	match_window_ticks = match_window * COUNT_PER_NSEC;
	
	solar_nqisk_precut_thresh = 1000;
	m_variables.Get("solar_nqisk_precut_thresh", solar_nqisk_precut_thresh);  // [hits]
	
	// get RFM input reader
	std::string rfmReaderName;
	m_variables.Get("rfmReaderName", rfmReaderName);
	if(m_data->Trees.count(rfmReaderName)==0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+rfmReaderName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	rfmReader = m_data->Trees.at(rfmReaderName);
	m_data->vars.Set("SolarRfmReader",rfmReaderName);
	
	// get relics we'll be matching against
	// there's not many per run, and we need to know if a solar event will be matched to relics that may come after it,
	// so just read in all the relics up front.
	
	//m_variables.Get("relic_fname",relic_fname);
	// Not a huge fan of it, but i think the best way to do this is generate the relic filename based on the run number.
	// The problem with this is we're hiding an assumption on input file name and location in hard-code...
	// The alternatives are specifying the filename as a config constant (correct way? but a bit of a pain)
	// or coding means to derive it from the available info (e.g. input file/run number - which is tricky to implement generally)
	char relic_fname[255];
	snprintf(relic_fname,256,"/disk03/lowe12/warwick/sk6/RelicMuonMatching/test_jacks_%06d.root",skhead_.nrunsk);
	std::cout<<"run "<<skhead_.nrunsk<<" lead to file "<<relic_fname<<std::endl;
	TFile f_relic(relic_fname,"READ");
	if(f_relic.IsZombie()){
		Log(m_unique_name+" Error! failed to open relic file '"+relic_fname+"'",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	TTree* t_relic=(TTree*)f_relic.Get("relic");
	if(t_relic==nullptr || t_relic->GetEntries()==0){
		Log(m_unique_name+" Error! relic file '"+relic_fname+"' has no 'relic' tree or relic tree has no entries!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	Header* header_relic = new Header{};
	LoweInfo* lowe_relic = new LoweInfo{};
	t_relic->SetBranchAddress("HEADER",&header_relic);
	t_relic->SetBranchAddress("LOWE",&lowe_relic);
	for(int i=0; i<t_relic->GetEntries(); ++i){
		t_relic->GetEntry(i);
		// struct to encapsulate the key info
		SolarRelic a_relic;
		a_relic.out_entry_num=i;
		a_relic.nevsk=header_relic->nevsk;
		
		// combined timestamp of this event
		a_relic.ticks = (header_relic->counter_32 & ~0x1FFFF);
		a_relic.ticks = a_relic.ticks << 15;
		uint64_t iticks = *reinterpret_cast<uint32_t*>(&header_relic->t0);
		a_relic.ticks += iticks;
		
		// reconstructed vertex
		a_relic.vtx = {lowe_relic->bsvertex[0],lowe_relic->bsvertex[1],lowe_relic->bsvertex[2]};
		
		// record this relic in the map
		relics_this_run.push_back(a_relic);
	}
	std::cout<<"we had "<<relics_this_run.size()<<" relics this run"<<std::endl;
	t_relic->ResetBranchAddresses();
	f_relic.Close();
	delete header_relic;
	delete lowe_relic;
	
	// we'll do time matching of solar events to relics here, but distance matching will come later
	// (after we reconstruct the solar event), so we can't make the final matching decision here alone.
	// so put the relics in the DataModel to pass down to the SolarPostSelection tool for that.
	intptr_t solar_relics_ptr = reinterpret_cast<intptr_t>(&relics_this_run);
	m_data->CStore.Set("solar_relics", solar_relics_ptr);
	
	// add some extra branches to solar event tree
	int lun = m_data->GetLUN(rfmReaderName);
	TreeManager* mgr = skroot_get_mgr(&lun);
	solarTree = mgr->GetOTree();
	// rename tree to solar as we'll add a relic tree in a second
	solarTree->SetName("solar");
	solarTree->Branch("HwClockTicks",&thiseventticks);
	// note the set of relics (time) matched to a solar is fixed here - if this solar is later rejected, it just won't be recorded
	solarTree->Branch("MatchedEvNums",&matched_ev_nums);
	solarTree->Branch("MatchedOutEntryNums",&matched_out_entry_nums);
	solarTree->Branch("MatchedTimeDiffs",&matched_tdiffs);
	
	// also add an output relic tree with 1:1 matching to input relic tree
	solarTree->GetCurrentFile()->cd();
	relic_tree = new TTree("relic","relic");
	relic_tree->Branch("nevsk",&out_relic.nevsk);
	relic_tree->Branch("MatchedEvNums",&out_relic.matched_ev_nums);
	relic_tree->Branch("MatchedInEntryNums",&out_relic.matched_in_entry_nums);
	relic_tree->Branch("MatchedOutEntryNums",&out_relic.matched_out_entry_nums);
	relic_tree->Branch("MatchedTimeDiffs",&out_relic.matched_tdiffs);
	relic_tree->Branch("MatchedDists",&out_relic.matched_dists);
	relic_tree->Branch("RejectedBy",&out_relic.rejected_by);
	relic_tree->Branch("Rejected", &out_relic.rejected);
	
	// see if recording this as a cut, and if so make the cut recorder
	get_ok = m_variables.Get("solarSelectorName",solarSelectorName);
	if(get_ok){
		m_data->AddCut(solarSelectorName, "solar_nqisk", "solar id charge in 1.3us",true,0,solar_nqisk_precut_thresh);
		m_data->AddCut(solarSelectorName, "relic_solar_tdiff", "times to solar events within 60s",true,-60,60);
	}
	
	return true;
}


bool SolarPreSelection::Execute(){
	
	++solarcount;
	// 1. apply any pre-selection requirements of the solar event
	
	// since we're going to apply a bse < 25 MeV cut in SolarPostSelection,
	// but bonsai takes ages to run especially on high energy events,
	// apply a very loose pre-cut based on ID charge
	if(!solarSelectorName.empty()){
		m_data->ApplyCut(solarSelectorName, "solar_nqisk", skq_.nqisk);
	}
	if(skq_.nqisk > solar_nqisk_precut_thresh){
		Log(m_unique_name+" solar event over ID Q threshold, skip",v_debug,m_verbose);
		m_data->vars.Set("Skip",1);
		return true;
	}
	
	// 2. check if the current event is within 60s of any relics.
	
	// get timestamp of this event
	thiseventticks = (skheadqb_.nevhwsk & ~0x1FFFF);
	thiseventticks = thiseventticks << 15;
	uint64_t iticks = *reinterpret_cast<uint32_t*>(&skheadqb_.it0sk);
	thiseventticks += iticks;
	
	matched_ev_nums.clear();
	matched_out_entry_nums.clear();
	matched_tdiffs.clear();
	matched_indices.clear();
	
	
	// scan relics for events within 60s
	int relici=0;
	int matched_relic_count=0;
	int to_erase_count=0;
	for(std::deque<SolarRelic>::iterator it=relics_this_run.begin(); it!=relics_this_run.end(); ++it){
		SolarRelic& a_relic = *it;
		
		if(a_relic.nevsk==skhead_.nevsk) continue; // don't match it to itself
		
		uint64_t& thisrelicsticks = a_relic.ticks;
		int64_t ticksDiff = (a_relic.nevsk<skhead_.nevsk) ? (thiseventticks - thisrelicsticks) : (thisrelicsticks - thiseventticks);
		// ticksdiff should be positive (since sign is made positive by nevsk comp)- if not, rollover has happened
		if(ticksDiff<0) ticksDiff += (int64_t(1) << 47);
		double t_diff_sign = (a_relic.nevsk<skhead_.nevsk) ? -1. : 1.; // sign relative to solar (as this is solar tree)
		
		//std::cout<<"solar event "<<skhead_.nevsk<<", relic "<<relici++<<" event "<<a_relic.nevsk
		//         <<" tdiff = "/*<<ticksDiff<<" / "*/<<(ticksDiff/(COUNT_PER_NSEC*1E9))<<"s"<<std::endl;
		
		if(!solarSelectorName.empty()){
			m_data->ApplyCut(solarSelectorName, "relic_solar_tdiff", t_diff_sign*(ticksDiff/COUNT_PER_NSEC)/1.E9);
		}
		
		if(ticksDiff > match_window_ticks){
			if(a_relic.nevsk < skhead_.nevsk){
				// if this relic is >60s before current event,
				// it won't be matched to any future solars, so can be written out
				std::cout<<"relic "<<relici<<" ready for writeout"<<std::endl;
				FillRelic(a_relic);
				++to_erase_count;
				if(&a_relic==&relics_this_run.back()){
					// if we've written out all the relics, no further solars will be of interest
					m_data->vars.Set("StopLoop",1);
				}
				// don't handle this relic as a match
				continue;
			} else {
				// if this relic is >60s after the current event,
				// this and no later relics will match to the current solar event
				break;
			}
		}
		
		// if we got here, this solar event is within 60s of this relic!
		Log(m_unique_name+" matched solar with tdiff "+toString(t_diff_sign*(ticksDiff/COUNT_PER_NSEC)/1.E9)+"s",v_debug,m_verbose);
		
		// add the match in the solar tree branch variables
		matched_ev_nums.push_back(a_relic.nevsk);
		matched_indices.push_back(matched_relic_count++);
		matched_out_entry_nums.push_back(erased_count + to_erase_count + matched_relic_count); // index in full relics_this_run
		matched_tdiffs.push_back(t_diff_sign*(ticksDiff/COUNT_PER_NSEC)/1.E9);
		
	}
	
	// handle erasure of relics written out
	if(to_erase_count) relics_this_run.erase(relics_this_run.begin(),relics_this_run.begin()+to_erase_count);
	erased_count += to_erase_count;
	
	if(matched_ev_nums.empty()){
		// this solar event was not within 60s of any relics, no further interest
		Log(m_unique_name+" solar event not within 60s of any relics, skip",v_debug,m_verbose);
		m_data->vars.Set("Skip",1);
	} else {
		Log(m_unique_name+" solar event in proxomity to relic, move to reconstruct...",v_debug,m_verbose);
		// for solar events that pass post-selection, matches and time diffs
		// will be recorded in the output relic tree, so pass down to that Tool
		m_data->CStore.Set("solar_relic_matches",matched_indices);
		m_data->CStore.Set("solar_relic_tdiffs",matched_tdiffs);
		++matchedsolarcount;
	}
	
	return true;
}


bool SolarPreSelection::Finalise(){
	
	// write out any remaining relics
	for(SolarRelic& a_relic : relics_this_run){
		FillRelic(a_relic);
	}
	relics_this_run.clear();
	relic_tree->ResetBranchAddresses();
	
	std::cout<<"checked "<<solarcount<<" solar events and "<<relics_this_run.size()<<" relics"<<std::endl;
	std::cout<<"found "<<matchedsolarcount<<" solar candidates that were within 60s of a relic (before resconstruction and solar selection cuts)"<<std::endl;
	
	return true;
}

void SolarPreSelection::FillRelic(SolarRelic& a_relic){
	out_relic = a_relic;
	relic_tree->Fill();
	return;
}
