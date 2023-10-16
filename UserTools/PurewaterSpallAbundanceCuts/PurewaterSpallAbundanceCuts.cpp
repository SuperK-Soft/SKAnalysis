/* vim:set noexpandtab tabstop=4 wrap */
#include "PurewaterSpallAbundanceCuts.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <chrono>      // std::chrono::seconds
#include <cmath>       // fabs
#include <algorithm>   // std::max_element

#include "Constants.h" // muboy_class

PurewaterSpallAbundanceCuts::PurewaterSpallAbundanceCuts():Tool(){}

const double li9_endpoint = 14.5; // MeV
// various dt cuts to assess systematics on the dlt cut efficiency
// these are assorted lifetimes of N16, Li8, Li9, B12, B13 respectively
const std::vector<float> spall_lifetimes{7.13, 0.838, 0.178, 0.0202, 0.0174};

bool PurewaterSpallAbundanceCuts::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);            // how verbose to be
	m_variables.Get("treeReaderName",treeReaderName);  // reader name for input
	m_variables.Get("outputFile",outputFile);          // output file to write
	
	// Get cut thresholds
	// ------------------
	m_variables.Get("max_closest_muon_dt",max_closest_muon_dt); // cut lowe events too close to a mu (pre OR post!)
	m_variables.Get("max_closest_lowe_dx",max_closest_lowe_dx); // cut lowe events too close to another within 60s
	m_variables.Get("ntag_FOM_threshold",ntag_FOM_threshold);
	m_variables.Get("run_min",run_min);
	m_variables.Get("run_max",run_max);
	
	// get the reader for accessing input file branches
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	// Set up the tree selector to operate on entries in this tree
	myTreeSelections.SetTreeReader(myTreeReader);
	myTreeSelections.MakeOutputFile(outputFile);
	// Set the tree selector so that downstream tools can check which events pass which cuts
	m_data->Selectors.emplace(treeReaderName, &myTreeSelections);
	
	// Pre-populate event count tracker with all cuts to get the ordering of reduction right.
	// FIXME find a better way. Maybe we should split lowe events // mu-lowe pairs // mu-lowe-ntag triplets,
	// as otherwise it need not be monotonic. Easiest way for now is:
	// `grep myTreeSelections.IncrementEventCount UserTools/PurewaterSpallAbundanceCuts/PurewaterSpallAbundanceCuts.cpp`
	// TODO build cut names to reflect actual cuts - e.g. "SNR > "+ std::to_string(snr_cut)
	// HAVE config file specify "active" for all cuts
	std::vector<std::pair<std::string, std::vector<std::string>>> cut_names{
		// cut name						// list of branches whose indices are required to identify this event
		{"all",							{}},   // order of the branches specified here must match
		{"61525<run<73031",				{}},   // the order of indices given to AddPassingEvent!
		//{"SNR>0.5",					{}},   // already applied as part of online first reduction*
		{"dwall>200cm",					{}},   // already applied as part of online first reduction*
		{"dt_mu_lowe>50us",				{}},   // already applied as part of online first reduction*
		//{"Q50/N50>0.5",				{}},   // *per Ashida-san, although dt_mu_lowe>50us does cut some events?
		//{"thirdred",					{}},
		//{"max_hits_200ns_AFT<50",		{}},
		//{"min(dt_mu_lowe)>1ms",		{}},
		//{"nearest_other_lowe>490cm",	{}},
		{"lowe_energy>6MeV",			{}},
		{"pre_muon_muboy_i==0",			{"mubitrack"}},
		{"pre_mu_dt_cut_0",				{"spadt"}},      // # of pre-muons passing dlt cut after different dt cuts
		{"pre_mu_dt_cut_1",				{"spadt"}},
		{"pre_mu_dt_cut_2",				{"spadt"}},
		{"pre_mu_dt_cut_3",				{"spadt"}},
		{"pre_mu_dt_cut_4",				{"spadt"}},
		{"post_muon_muboy_i==0",		{"mubitrack"}},
		{"post_mu_dt_cut_0",			{"spadt"}},      // # of post-muons passing dlt cut after different dt cuts
		{"post_mu_dt_cut_1",			{"spadt"}},
		{"post_mu_dt_cut_2",			{"spadt"}},
		{"post_mu_dt_cut_3",			{"spadt"}},
		{"post_mu_dt_cut_4",			{"spadt"}},
		{"mu_lowe_pairs",				{"spadt"}},
		{"pre-mu_dt<0",					{"spadt"}},
		{"muboy_index==0",				{"mubitrack"}},
		{"dlt_mu_lowe>200cm",			{"spadlt"}},
		{"lowe_energy_in_li9_range",	{}},
		{"dt_mu_lowe_in_li9_range",		{"spadt"}},
		{"closest_other_mu_dt>1ms",		{}},
		{"ntag_FOM>0.995",				{}},              // no ntag index; just require event had 1+ passing ntag
		{"mu_lowe_ntag_triplets",		{"spadt","dt"}}
	};
	for(auto&& acut : cut_names) myTreeSelections.AddCut(acut.first, acut.first, false, acut.second);
	
	// pass the TreeSelections to downstream tools so they can see which events passed which selections
	intptr_t myTreeSelectionsPtr = reinterpret_cast<intptr_t>(&myTreeSelections);
	m_data->CStore.Set("SpallAbundanceSelection",myTreeSelectionsPtr);
	
	return true;
}

bool PurewaterSpallAbundanceCuts::Execute(){
	
	// retrieve branch variables
	GetBranchValues();
	
	// Apply event selections
	Log(m_unique_name+" applying cuts",v_debug,m_verbose);
	Analyse();
	
	return true;
}

bool PurewaterSpallAbundanceCuts::Finalise(){
	
	Log(m_unique_name+" event counts trace: ",v_warning,m_verbose);
	myTreeSelections.PrintCuts();
	
	// write out the event numbers that passed each cut
	myTreeSelections.Write();
	
	// add the last run duration to the livetime and note for downstream tools
	AddLastRunTime();
	m_data->CStore.Set("livetime", livetime);
	
	return true;
}

bool PurewaterSpallAbundanceCuts::GetBranchValues(){
	
	// retrieve variables from TTree to member variables
	int success = 
	(myTreeReader->Get("HEADER", HEADER)) &&
	(myTreeReader->Get("LOWE", LOWE)) &&
	(myTreeReader->Get("ThirdRed", thirdredvars)) &&
	(myTreeReader->Get("np", num_neutron_candidates)) &&
	(myTreeReader->Get("N200M", max_hits_200ns_AFT)) &&
	(myTreeReader->Get("neutron5", ntag_FOM)) &&
	(myTreeReader->Get("nmusave_pre", num_pre_muons)) &&
	(myTreeReader->Get("nmusave_post", num_post_muons)) &&
	(myTreeReader->Get("mubstatus", mu_class)) &&
	(myTreeReader->Get("mubitrack", mu_index)) &&
	(myTreeReader->Get("mubgood", mu_fit_goodness)) &&
	(myTreeReader->Get("spadt", dt_mu_lowe)) &&
	(myTreeReader->Get("spadlt", dlt_mu_lowe)) &&
	(myTreeReader->Get("multispa_dist", closest_lowe_60s));
	
	return success;
}

// #####################################################################

// main body of the tool
bool PurewaterSpallAbundanceCuts::Analyse(){
	// This gets called for each Execute iteration, to process one lowe event
	// Apply all cuts in sequence, noting which events pass each cut
	
	myTreeSelections.AddPassingEvent("all");
	entry_number = myTreeReader->GetEntryNumber();
	Log(m_unique_name+" entry "+toString(entry_number)+" run "+toString(HEADER->nrunsk),v_debug,m_verbose);
	
	Log(m_unique_name+" checking run cut",v_debug+1,m_verbose);
	if(HEADER->nrunsk < run_min) return false;    // start of SK-IV, Oct 2008
	if(HEADER->nrunsk > run_max){                 // Yang Zhang's time range, Oct 2014
		m_data->vars.Set("StopLoop",1);
		Log(m_unique_name+" entry "+toString(entry_number)+" run "+toString(HEADER->nrunsk)
					+" beyond final run number "+toString(run_max)+", stopping ToolChain",v_warning,m_verbose);
		return false;
	}
	//if (HEADER->nrunsk > 74781) return false;  // WIT started after this run. What's the significance of this?
	myTreeSelections.AddPassingEvent("61525<run<73031");
	
	IncrementLivetime();
	
	// find lowe events ✅
	// for reference, paper says 54,963 beta events... though not clear after which cuts
	
//	// n_hits_with_Q_lt_0.5pe / n_hits_total > 0.55
//	Log(m_unique_name+" checking SNR cut",v_debug+1,m_verbose);
//	// TODO how is this implemented...? (supposedly already applied as part of first reduction)
//	myTreeSelections.AddPassingEvent("SNR>0.5");
	
//	// Q50/N50 cut is something different-  not mentioned by paper, is this another of sonia's additional cuts?
//	if((double)thirdredvars->q50 / (double)(LOWE->bsn50) > 2.) return false;
//	myTreeSelections.AddPassingEvent("Q50/N50>0.5");
	
	// dwall > 2m
	Log(m_unique_name+" checking dwall cut",v_debug+1,m_verbose);
	if(thirdredvars->dwall < 200.) return false;
	myTreeSelections.AddPassingEvent("dwall>200cm");
	
	// dt_muon_lowe > 50us
	// check the closest preceding muon.
	// n.b. if muboy found multiple, they all have the same time, so we don't neeed to scan
	Log(m_unique_name+" checking afterpulsing cut",v_debug+1,m_verbose);
	if(fabs(dt_mu_lowe[num_pre_muons-1]) < 50e-6) return false;
	myTreeSelections.AddPassingEvent("dt_mu_lowe>50us");
	
//	// new set of cuts from atmospheric analysis TODO enable?
//	Log(m_unique_name+" checking third reduction cuts",v_debug+1,m_verbose);
//	if (apply_third_reduction(third, LOWE)) return false;
//	myTreeSelections.AddPassingEvent("thirdred");
	
//	// new cut of events with a low energy muon in the ncapture window TODO enable?
//	Log(m_unique_name+" checking for muons in AFT trigger",v_debug+1,m_verbose);
//	if (max_hits_200ns_AFT > 50) return false;
//	myTreeSelections.AddPassingEvent("max_hits_200ns_AFT<50");
	
//	// new cut of events within 1ms of nearest muon (pre- or post) cut TODO enable?
//	// how would this work? We're looking for spallation: we want events in proximity to a muon...
//	Log(m_unique_name+" checking for closest pre-muon >1ms",v_debug+1,m_verbose);
//	bool fail_1ms = false;
//	for ( int i = 0; i < num_pre_muons; i++ ) {
//		if ( fabs( dt_mu_lowe[ i ] ) < 0.001 ){
//			fail_1ms = true;
//			break;
//		}
//	}
//	if(fail_1ms) return false;
//	myTreeSelections.AddPassingEvent("min(dt_mu_lowe)>1ms");
	
//	// new cut to remove all lowe events in close proximity to another lowe event TODO enable?
//	Log(m_unique_name+" checking closest lowe event within 60s",v_debug+1,m_verbose);
//	if (closest_lowe_60s[0] < max_closest_lowe_dx) return false;    // not applied by Zhang
//	myTreeSelections.AddPassingEvent("nearest_other_lowe>490cm");
	
	// cut all lowe events with energy < 6 MeV (mostly non-spallation)
	Log(m_unique_name+" checking lowe energy > 6 MeV",v_debug+1,m_verbose);
	if( LOWE->bsenergy < 6.f ) return false;
	myTreeSelections.AddPassingEvent("lowe_energy>6MeV");
	
	// find muons within 30 prior to lowe events ( muons identified by: >1000pe in ID ) ✅
	// find muons within 30s post lowe events    ( muons identified by: >1000pe in ID ) ✅
	// pair all lowe events with all pre-muons ✅
	// pair all lowe events with all post-muons ✅
	
	// plot lt, dt distributions as a function of muon class, for both pre- and post-muons
	Log(m_unique_name+" looping over "+toString(num_pre_muons)+" preceding muons and "+toString(num_post_muons)
				+" following muons to fill spallation and control dl, dt histograms",v_debug,m_verbose);
	// pre muons
	for(size_t mu_i=0; mu_i<num_pre_muons; ++mu_i){
		// only consider first muboy muon (only for multi-mu events?)
		Log(m_unique_name+" checking muboy index==0",v_debug+2,m_verbose);
		if(mu_index[mu_i] > 0) continue;
		myTreeSelections.AddPassingEvent("pre_muon_muboy_i==0",mu_i);
		
		// to evaluate systematic on lt cut, apply various dt cuts and see how the lt cut efficiency varies
		// since we're interested in the effect on the spallation sample, which is given by
		// the total - post-muon sample, record both pre- and post- muon samples with various dt cuts
		for(int dt_cut_i=0; dt_cut_i<spall_lifetimes.size(); ++dt_cut_i){
			Log(m_unique_name+" checking nominal dlt cut",v_debug+2,m_verbose);
			if(fabs(dt_mu_lowe[mu_i]) < spall_lifetimes.at(dt_cut_i)){
				myTreeSelections.AddPassingEvent(std::string("pre_mu_dt_cut_")+toString(dt_cut_i),mu_i);
			}
		}
	}
	// post muons
	for(size_t mu_i=num_pre_muons; mu_i<(num_pre_muons+num_post_muons); ++mu_i){
		Log(m_unique_name+" checking muboy index==0",v_debug+2,m_verbose);
		if (mu_index[mu_i] > 0) continue;
		myTreeSelections.AddPassingEvent("post_muon_muboy_i==0",mu_i);
		
		for(int dt_cut_i=0; dt_cut_i<spall_lifetimes.size(); ++dt_cut_i){
			Log(m_unique_name+" checking nominal dlt cut",v_debug+2,m_verbose);
			if(dt_mu_lowe[mu_i] < spall_lifetimes.at(dt_cut_i)){
				myTreeSelections.AddPassingEvent(std::string("post_mu_dt_cut_")+toString(dt_cut_i),mu_i);
			}
		}
	}
	// in Finalise we'll substract the two to get dt and dlt distributions for spallation only.
	// we'll also compare across various dt cuts to get the systematic error on the spallation dlt cut.
	
	// the following cuts are based on muon-lowe pair variables, so loop over muon-lowe pairs
	Log(m_unique_name+" Looping over "+toString(num_pre_muons)
				+" preceding muons to look for spallation events",v_debug,m_verbose);
	for(size_t mu_i=0; mu_i<num_pre_muons; ++mu_i){
		myTreeSelections.AddPassingEvent("mu_lowe_pairs", mu_i);
		
		// safety check: should not consider muons after lowe event
		Log(m_unique_name+" checking dt of pre-muon <0",v_debug+2,m_verbose);
		if (dt_mu_lowe[mu_i] >= 0) break;
		myTreeSelections.AddPassingEvent("pre-mu_dt<0", mu_i);
		
		Log(m_unique_name+" checking muboy index==0",v_debug+2,m_verbose);
		if (mu_index[mu_i] > 0) continue;
		myTreeSelections.AddPassingEvent("muboy_index==0", mu_i);
		
		// Apply nominal lt cut, unless muon type was misfit or a poorly fit single muon
		Log(m_unique_name+" checking nominal dlt cut",v_debug+2,m_verbose);
		if(not (dlt_mu_lowe[mu_i] < 200 || mu_class[mu_i] == int(muboy_class::misfit) || 
			   (mu_class[mu_i] == int(muboy_class::single_thru_going) && mu_fit_goodness[mu_i] < 0.4)))
			    continue;
		myTreeSelections.AddPassingEvent("dlt_mu_lowe>200cm", mu_i);
		
		// That's all for assessing the amount of general spallation isotopes
		// ------------------------------------------------------------------
		// for Li9 we have a number of additional cuts
		
		// apply Li9 energy range cut
		Log(m_unique_name+" checking li9 energy range cut",v_debug+2,m_verbose);
		if ( LOWE->bsenergy <= 7.5 || LOWE->bsenergy >= li9_endpoint ) continue;
		myTreeSelections.AddPassingEvent("lowe_energy_in_li9_range");
		
		// apply Li9 lifetime cut
		Log(m_unique_name+" checking li9 lifetime cut",v_debug+2,m_verbose);
		if( dt_mu_lowe[mu_i] > -0.05 || dt_mu_lowe[mu_i] < -0.5 ) continue;
		myTreeSelections.AddPassingEvent("dt_mu_lowe_in_li9_range", mu_i);
		
		// no other mu within 1ms of this lowe event
		Log(m_unique_name+" checking for another muon within 1ms",v_debug+2,m_verbose);
		bool other_muon_within_1ms = false;  // XXX check we're interpreting this cut right
		for(int othermu_i=0; othermu_i<(num_pre_muons+num_post_muons); ++othermu_i){
			if(othermu_i==mu_i) continue; // looking for muons other than the current one
			if(fabs(dt_mu_lowe[othermu_i])<max_closest_muon_dt){
				other_muon_within_1ms = true;
				break;
			}
		}
		if(other_muon_within_1ms) continue;
		myTreeSelections.AddPassingEvent("closest_other_mu_dt>1ms");
		
		// search for ncapture candidates: >7 hits within 10ns T-TOF in 50ns-535us after lowe events ✅
		
		// calculate neutron FOM and cut failing ones
		Log(m_unique_name+" checking for a neutron passing BDT cut",v_debug+2,m_verbose);
		if( (num_neutron_candidates==0) || 
			(*std::max_element(ntag_FOM.begin(), ntag_FOM.end())<ntag_FOM_threshold)) continue;
		// XXX as reference, we should have 116 remaining candidate events here
		myTreeSelections.AddPassingEvent("ntag_FOM>0.995");
		
		// apparently in Zhang study no events had multiple ntag candidates
		if(num_neutron_candidates>1){
			Log(m_unique_name+" event with "+toString(num_neutron_candidates)
				+" neutron candidates!",v_debug,m_verbose);
		}
		
		// plot distribution of beta->ntag dt from passing triplets, compare to fig 5
		// Zhang had no events with >1 ntag candidate: should we only take the first? XXX
		for(size_t neutron_i=0; neutron_i<num_neutron_candidates; ++neutron_i){
			myTreeSelections.AddPassingEvent("mu_lowe_ntag_triplets", {mu_i, neutron_i});
		}
		
	} // end loop over muons
	
	return true;
}

// #####################################################################

void PurewaterSpallAbundanceCuts::IncrementLivetime(){
	// if this is the first event of a new run, add the run duration to the livetime
	if(HEADER->nrunsk > current_run){
		current_run = HEADER->nrunsk;
		// set the time from the offsets since Jan 1st 1900?
		//update the livetime with the run duration
		livetime += difftime(mktime(&runstart), mktime(&runend));
		
		// set the new run start to this event's timestamp
		runstart.tm_year = HEADER->ndaysk[0];
		runstart.tm_mon = HEADER->ndaysk[1] - 1;
		runstart.tm_mday = HEADER->ndaysk[2];
		runstart.tm_hour = HEADER->ntimsk[0];
		runstart.tm_min = HEADER->ntimsk[1];
		runstart.tm_sec = HEADER->ntimsk[2];
	}
	
	// every time update the end of run timestamp to this event's timestamp.
	runend.tm_year = HEADER->ndaysk[0];
	runend.tm_mon = HEADER->ndaysk[1] - 1;
	runend.tm_mday = HEADER->ndaysk[2];
	runend.tm_hour = HEADER->ntimsk[0];
	runend.tm_min = HEADER->ntimsk[1];
	runend.tm_sec = HEADER->ntimsk[2];
}

void PurewaterSpallAbundanceCuts::AddLastRunTime(){
	// if this is the last event we'll process, add the run livetime
	livetime += difftime(mktime(&runstart), mktime(&runend));
}

bool PurewaterSpallAbundanceCuts::apply_third_reduction(const ThirdRed *th, const LoweInfo *LOWE){
	Log(m_unique_name+" applying third reduction",v_debug,m_verbose);
	if (th->maxpre >= 12)             return true;                // remove events with high pre-activity
	if (th->nmue > 0)                 return true;                // remove events with a muon decay electron?
	if (th->dwall < 200.)             return true;                // remove events within 2m from wall
	if (th->effwall < 500.)           return true;                // remove events with back projection <5m to wall
	if (th->pilike >= 0.36)           return true;                // remove pion-like events
	if ((double)th->q50 / (double)(LOWE->bsn50) > 2) return true; // remove events with low signal-to-noise
	if (th->angle<38 || th->angle>50) return true;                // remove events with Cherenkov angle != 42°
	if (LOWE->bsenergy<8)             return true;                // remove events with energy < 8 MeV
	if ( LOWE->bsgood[1] < 0.5 )      return true;                // remove events with poor vertex fit
	if (th->ovaq < 0.25 )             return true;                // remove events poor event quality
	return false;
}
