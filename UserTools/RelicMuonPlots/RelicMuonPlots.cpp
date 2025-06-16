#include "RelicMuonPlots.h"
#include "MTreeReader.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include <cmath>
#include <bitset>
#include "geotnkC.h"  // for SK tank geometric constants
#include "type_name_as_string.h"

RelicMuonPlots::RelicMuonPlots():Tool(){}

// TODO move to Constants
namespace {
	// from $RELIC_WORK_DIR/data_reduc/spallation/spallation_observables/mu_info.C
	// charge per cm of a minimum-ionizing-particle (MIP), in hits per cm, presumably
	static const float pe_per_cm = 26.78;
	// trying to see where this comes from.... muon rate of energy loss in water:
	// 1.992 Mev cm^2/g dE/dx min   (from: Muon Stopping Power and Range Tables, D.E.Groom et al, LBNL-44742,
	// Atomic Data and Nuclear Data Tables, Vol. 76, No. 2, July 2001 - ~/Downloads/adndt.pdf)
	// at 13.5C density of water 0.999315 g/cm^3
	// (from: https://www.internetchemistry.com/chemical-data/water-density-table.php )
	// so 1.992*0.999315 = 1.99063548 MeV/cm (makes sense with commonly known value of muon loss ~2MeV/cm)
	// but how to translate this MeV to pe? Need hits per MeV; petable? surely depends on water transparency?
}

bool RelicMuonPlots::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("relicReaderName",relicReaderName);
	m_variables.Get("muReaderName",muReaderName);
	
	// get intput treereaders
	if(m_data->Trees.count(relicReaderName)==0){
		Log(m_unique_name+" Error! No TreeReader '"+relicReaderName+"' in DataModel!",v_error,m_verbose);
		return false;
	}
	relicReader = m_data->Trees.at(relicReaderName);
	
	if(m_data->Trees.count(muReaderName)==0){
		Log(m_unique_name+" Error! No TreeReader '"+muReaderName+"' in DataModel!",v_error,m_verbose);
		return false;
	}
	muReader = m_data->Trees.at(muReaderName);
	
	// make output files
	m_variables.Get("outputFile", outputFile);  // pair (spall) variables
	m_variables.Get("relicFile", relicFile);    // relic variables
	m_variables.Get("muFile", muFile);          // muon variables
	m_variables.Get("pe_table_file", pe_table_filename);
	if(pe_table_filename==""){
		char* skofl_dir = getenv("SKOFL_ROOT");
		std::string empty = "";
		if(skofl_dir==nullptr) skofl_dir = const_cast<char*>(empty.c_str());
		pe_table_filename = std::string(skofl_dir) + "/const/lowe/petable.dat";
	}
	get_ok = GetPeTable();
	if(!get_ok){
		Log(m_unique_name+" Error reading petable file '"+pe_table_filename+"'",
		    v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	MakeHists(0);
	
	return true;
}


bool RelicMuonPlots::Execute(){
	
	// get next relic candidate
	get_ok = GetRelicEvt();
	if(!get_ok){
		Log(m_unique_name+" Error getting relic entry "+toString(relicEntryNum),
		    v_error,m_verbose);
		return false;
	}
	
	// make sure it was reconstructed
	if(relicLowe->bsenergy == 0){             // reconstruction not done
		Log(m_unique_name+" relic bsenergy == 0, was lowe reconstruction done?",v_error,m_verbose);
		return false;
	} else if(relicLowe->bsenergy == 9999){   // reconstruction failed
		Log(m_unique_name+" relic bsenergy == 9999, reconstruction failed, skipping",v_message,m_verbose);
		return true;
	}
	
	// reset last muon event number so that mu_to_mu time matching is only recorded
	// for muons associated to the same relic
	lastmu_nevsk=0;
	
	// loop over muons matched to this relic
	Log(m_unique_name+"looping over "+std::to_string(relicMatchedEntryNums->size())+" muons for this relic",v_debug,m_verbose);
	//std::cout<<"relicTimeDiffs is of size "<<relicTimeDiffs->size()<<" vs "
	//         <<relicMatchedEntryNums->size()<<" matches"<<std::endl;
	
	// scan for dups
	std::set<int> unique_muon_entry_nums;
	std::vector<int> unique_entry_posns;
	for(size_t i=0; i<relicMatchedEntryNums->size(); ++i){
		int next_mu_ent_num = relicMatchedEntryNums->at(i);
		if(unique_muon_entry_nums.emplace(next_mu_ent_num).second==false){
			Log(m_unique_name+" removing duplicate muon "+toString(next_mu_ent_num)
			    +" in matches for relic "+toString(relicEntryNum),v_error,m_verbose);
		} else {
			unique_entry_posns.push_back(i);
		}
	}
	
	//for(size_t i=0; i<relicMatchedEntryNums->size(); ++i){
	for(size_t u=0; u<unique_entry_posns.size(); ++u){
		
		int i=unique_entry_posns.at(u);
		
		// get the next matched muon entry number in muon tree
		muEntryNum = relicMatchedEntryNums->at(i);
		
		// sanity check that this pair is unique
		std::pair<int,int> pair2{relicEntryNum,muEntryNum};
		if(pairs_compared2.count(pair2)!=0){
			Log(m_unique_name+" Error! Duplicate comparison between relic entry "+toString(relicEntryNum)
			    +" and muon entry "+toString(muEntryNum),v_error,m_verbose);
			return false;
		}
		pairs_compared2.emplace(pair2);
		
		// get the muon entry
		get_ok = GetMuonEvt();
		if(!get_ok){
			Log(m_unique_name+" Error getting muon entry "+toString(muEntryNum),v_error,m_verbose);
			continue;
		}
		
		// sanity check that this pair is unique
		std::pair<int,int> pair{relicEvNum,muEvNum};
		if(evmap.count(pair)!=0){
			Log(m_unique_name+" Error! Duplicate comparison between relic "+toString(relicEvNum)
			   +" and muon "+toString(muEvNum)+" for entries "+toString(relicEntryNum)
			   +" and "+toString(muEntryNum),v_error,m_verbose);
			std::cerr<<"previously seen with entries "<<evmap.at(pair).first
			         <<" and "<<evmap.at(pair).second<<std::endl;
			return false;
		}
		evmap.emplace(std::pair<std::pair<int,int>,std::pair<int,int>>{pair,pair2});
		
		dt = relicTimeDiffs->at(i);
		
		get_ok = MakePairVariables();
		if(get_ok) MakeHists(1);
		
	}
	
	// write to file every so often so we don't lose everything if toolchain crashes
	if((++loopnum%100)==0){
		hb.Save();
		mu_hb.Save();
		relic_hb.Save();
	}
	
	return true;
}


bool RelicMuonPlots::Finalise(){
	
	MakeHists(2);
	
	std::cout<<"found "<<spall_count<<" spall candidates (muon before relic) "
	         <<"and "<<rand_count<<" random candidates (muon after relic)"<<std::endl;
	
	return true;
}


bool RelicMuonPlots::GetRelicEvt(){
	
	// get next relic
	get_ok  = relicReader->Get("HEADER", relicHeader);
	get_ok &= relicReader->Get("LOWE", relicLowe);
	get_ok &= relicReader->Get("HwClockTicks", relicClockTicks);
	get_ok &= relicReader->Get("NumRollovers", relicRollovers);
	get_ok &= relicReader->Get("MatchedOutEntryNums", relicMatchedEntryNums);
	get_ok &= relicReader->Get("MatchedTimeDiff", relicTimeDiffs);
	get_ok &= relicReader->Get("TQREAL", relicTQReal);
	get_ok &= relicReader->Get("TQAREAL", relicTQAReal);
	
	if(!get_ok){
		Log(m_unique_name+" Error getting relic branches!",v_error,m_verbose);
		return false;
	}
	
	// get relic entry & event number
	relicEntryNum = relicReader->GetEntryNumber();
	relicEvNum = relicHeader->nevsk;
	
	if(relicEvNum==0){
		Log(m_unique_name+" 0 relic nevsk for entry "+toString(relicEntryNum),
		    v_error,m_verbose);
		return false;
	}
	
	// check if we've seen this relic before
	if(relic_nevsks.count(relicEvNum)!=0){
		Log(m_unique_name+" ERROR! Repeated nvesk in relic tree "+toString(relicEntryNum)
		     +" == "+toString(relic_nevsks.at(relicEvNum)),v_error,m_verbose);
		return false;
	}
	relic_nevsks.emplace(relicEvNum,relicEntryNum);
	
	if(relicMatchedEntryNums->size()!=relicTimeDiffs->size()){
		Log(m_unique_name+ "ERROR! relicMatchedEntryNums and relicTimeDiffs should be same size but aren't! ("
		     + std::to_string(relicMatchedEntryNums->size())+" vs "
		     + std::to_string(relicTimeDiffs->size()) +")",v_error,m_verbose);
		
		std::vector<int>* MatchedEvNums_vec_p = nullptr;
		std::vector<int>* MatchedInEntryNums_vec_p = nullptr;
		std::vector<int>* MatchedOutEntryNums_vec_p = nullptr;
		std::vector<bool>* MatchedEntryHasAFT_vec_p = nullptr;
		std::vector<float>* MatchedTimeDiff_vec_p = nullptr;
		std::vector<float>* MatchedParticleE_vec_p = nullptr;
		
		if(!relicReader->Get("MatchedEvNums", MatchedEvNums_vec_p)){
			std::cerr<<"error getting MatchedEvNums"<<std::endl;
		}
		if(!relicReader->Get("MatchedInEntryNums", MatchedInEntryNums_vec_p)){
			std::cerr<<"error getting MatchedInEntryNums"<<std::endl;
		}
		if(!relicReader->Get("MatchedOutEntryNums", MatchedOutEntryNums_vec_p)){
			std::cerr<<"error getting MatchedOutEntryNums"<<std::endl;
		}
		if(!relicReader->Get("MatchedEntryHasAFT", MatchedEntryHasAFT_vec_p)){
			std::cerr<<"error getting MatchedEntryHasAFT"<<std::endl;
		}
		if(!relicReader->Get("MatchedTimeDiff", MatchedTimeDiff_vec_p)){
			std::cerr<<"error getting MatchedTimeDiff"<<std::endl;
		}
		if(!relicReader->Get("MatchedParticleE", MatchedParticleE_vec_p)){
			std::cerr<<"error getting MatchedParticleE"<<std::endl;
		}
		
		std::cout<<"MatchedEvNums.size() = "<<MatchedEvNums_vec_p->size()<<std::endl;
		std::cout<<"MatchedInEntryNums.size() = "<<MatchedInEntryNums_vec_p->size()<<std::endl;
		std::cout<<"MatchedOutEntryNums.size() = "<<MatchedOutEntryNums_vec_p->size()<<std::endl;
		std::cout<<"MatchedEntryHasAFT.size() = "<<MatchedEntryHasAFT_vec_p->size()<<std::endl;
		std::cout<<"MatchedTimeDiff.size() = "<<MatchedTimeDiff_vec_p->size()<<std::endl;
		std::cout<<"MatchedParticleE.size() = "<<MatchedParticleE_vec_p->size()<<std::endl;
		
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	basic_array<float> bsvertex(relicLowe->bsvertex);
	basic_array<float> bsdir(relicLowe->bsdir);
	float bsdwall = *reinterpret_cast<float*>(&relicLowe->linfo[9]);
	float bsgood = relicLowe->bsgood[1];
	float bonsai_e = relicLowe->bsenergy;
	
	// Fill debug plots
	relic_hb.Fill("relic_vertex", bsvertex[0], bsvertex[1], bsvertex[2]);
	relic_hb.Fill("relic_dir",bsdir[0],bsdir[1],bsdir[2]);
	relic_hb.Fill("bonsai_goodness",bsgood);
	relic_hb.Fill("dist_to_wall",bsdwall);
	relic_hb.Fill("bonsai_e",bonsai_e);
	
	relic_hb.Fill("nevsk",relicEvNum);
	relic_hb.Fill("NumIDhits", relicTQReal->nhits);
	relic_hb.Fill("NumODhits", relicTQAReal->nhits);
	relic_hb.Fill("qismsk", skq_.qismsk);  // FIXME this isn't going to be populated, see below
	// ^ this won't be correct if we're not using a treereader for the relic tree,
	// but there's no equivalent member of TQREAL directly.
	
	int64_t thiseventticks = (relicHeader->counter_32 & ~0x1FFFF);
	thiseventticks = (thiseventticks << 15);
	int64_t it0sk = (int64_t(relicHeader->t0) & 0xFFFFFFFF);
	thiseventticks += it0sk;
	Log(m_unique_name+"this event at time "+std::to_string(double(thiseventticks/COUNT_PER_NSEC))
	      +" ns since last rollover",v_debug,m_verbose);
	if(thiseventticks==0){
		Log(m_unique_name+" RELIC EVENT AT 0 TICKS! counter_32: "
		   +std::to_string(relicHeader->counter_32)+"; header dump: ",v_error,m_verbose);
		relicHeader->Dump();
	}
	
	int64_t rolloverticks = relicRollovers * (int64_t(1) << 47);
	Log(m_unique_name+"adding "+std::to_string(double(rolloverticks/COUNT_PER_NSEC))+" ns for "
	   +std::to_string(relicRollovers)+" rollovers",v_debug,m_verbose);
	double totalticks = thiseventticks + rolloverticks;
	Log(m_unique_name+"total ticks: "+std::to_string(totalticks)+" or "+std::to_string(double(totalticks/COUNT_PER_NSEC))
	   +" ns \n c.f. EventTicks: "+std::to_string(relicClockTicks)+" or "+std::to_string(double(relicClockTicks/COUNT_PER_NSEC))
	   +" ns",v_debug,m_verbose);
	thiseventticks += relicRollovers * (int64_t(1) << 47);
	relic_hb.Fill("event_time_secs",double(thiseventticks/COUNT_PER_NSEC)/1.E9);
	
	if(lastrelicticks!=0){
		int64_t ticksDiff = thiseventticks - lastrelicticks;
		if(ticksDiff<0) ticksDiff += (int64_t(1) << 47);
		if(ticksDiff==0){
			std::cerr<<"relic entry "<<relicEntryNum<<" ticksDiff is 0 from last!"<<std::endl;
			assert(false);
			exit(-1);
			return false;
		}
		double ns_since_last = double(ticksDiff)/COUNT_PER_NSEC;
		if((ns_since_last/1E9)>70E3){
			std::cerr<<"TOO MUCH!!"<<std::endl;
			assert(false);
			exit(-1);
			return false;
		}
		relic_hb.Fill("relic_to_relic_secs", ns_since_last/1E9);
		lastrelicticks = thiseventticks;
	} else {
		lastrelicticks = thiseventticks;
	}
	
	/* TODO debug PLOTS
	energy
	energy vs nqisk
	energy vs qismsk
	n od hits (sanity check: mu should have some, relic should have ... fewer (quantify))
	event num nevsk (check flat distribution)
	event time (run start + clock ticks + rollovers) (check flat distribution)
	event time vs nevsk - check straight line
	event rate (time since last particle of that type - e.g. time from muon to muon, relic to relic)
	*/
	
	return get_ok;
}

bool RelicMuonPlots::GetMuonEvt(){
	
	// read the requested muon entry from file using upstream TreeReader Tool
	// we use this in SK mode to read the hits into the common blocks
	// so we can (re-)do muon dedx calculations because they didn't work..??
	get_ok = m_data->getTreeEntry(muReaderName, muEntryNum);
	if(!get_ok){
		Log(m_unique_name+" Error reading Muon entry " +toString(muEntryNum),v_error,m_verbose);
		return false;
	}
	Log(m_unique_name+" Muon event is entry "+toString(muEntryNum)+" of file "
	    +muReader->GetFile()->GetName(),v_debug,m_verbose);
	
	// otherwise read the entry ok; get branches
	get_ok  = muReader->Get("HEADER", muHeader);
	get_ok &= muReader->Get("MU", muMu);
	get_ok &= muReader->Get("HwClockTicks", muClockTicks);
	get_ok &= muReader->Get("NumRollovers", muRollovers);
	get_ok &= muReader->Get("TQREAL", muTQReal);
	get_ok &= muReader->Get("TQAREAL", muTQAReal);
	//get_ok &= muReader->Get("MatchedTimeDiff", muTimeDiffs);
	
	if(!get_ok){
		Log(m_unique_name+" Error getting muon branches!",v_error,m_verbose);
		return false;
	}
	
	muEvNum = muHeader->nevsk;
	
	if(muEvNum==0){
		Log(m_unique_name+" 0 muon nevsk for entry "+toString(muEntryNum),
		    v_error,m_verbose);
		return false;
	}
	
	if(muon_nevsks.count(muEvNum)!=0 && muon_nevsks.at(muEvNum)!=muEntryNum){
		Log(m_unique_name+" ERROR! Repeated nvesk in muon tree "+toString(muEntryNum)
		     +" == "+toString(muon_nevsks.at(muEvNum)),v_error,m_verbose);
		return false;
	}
	muon_nevsks.emplace(muEvNum,muEntryNum);
	
	basic_array<float> muff_dir(muMu->mudir);                     // direction from mfmuselect
	basic_array<float> muboy_dir(muMu->muboy_dir);                // direction from muboy
	basic_array<float> bff_dir(muMu->mubff_dir);                  // direction from bff
	
	basic_array<float> muff_entrypoint(muMu->muentpoint);         // entry point from mfmuselect
	basic_array<float[10][4]> muboy_entrypoint(muMu->muboy_entpos);  // entry points from muboy
	basic_array<float> bff_entrypoint(muMu->mubff_entpos);        // entry point from bff
	
	float mu_qismsk = muMu->muinfo[0];            // qismsk before PMT saturation correction
	float mu_qismsk_corr = muMu->muqismsk;        // qismsk after PMT saturation correction
	float muboy_goodness = muMu->mugoodness;
	muboy_class muboy_status(muboy_class{muMu->muboy_status});
	float muboy_tracklen = muMu->muboy_length;    // of primary track? or sum?
	int muboy_ntracks = muMu->muboy_ntrack;
	int muboy_index = muMu->muinfo[7];
	bool didbff = muMu->muinfo[6];
	float bff_goodness = muMu->mubff_goodness;
	if(muboy_index>=muboy_ntracks){
		Log(m_unique_name+" Error! Muboy index out of range: "+toString(muboy_index)+"/"
		    +toString(muboy_ntracks),v_error,m_verbose);
		return false;
	}
	
	// we have several sources of reconstructed variables: mfmuselect, muboy, BFF,
	// as well as both kirk and scott's dE/dx arrays. Pick one set.
	if(didbff){
		// use BFF reconstruction variables. Note that afaik BFF is only invoked
		// when BFF reports a single muon, and BFF only reconstructs a single muon
		muon_entrypoint = const_cast<float*>(bff_entrypoint.data());
		muon_direction = const_cast<float*>(bff_dir.data());
		
		// BFF doesn't give a track length; it is designed for through-going muons
		// so calculate the track length based on entry point and direction
		muon_tracklen = CalculateTrackLen(muon_entrypoint, muon_direction);
		
	} else {
		// use MuBoy reco variables. MuBoy may reconstruct multiple muons,
		// each of which has a separate entry point and dE/dx.
		// for now, although it's inefficient on storage, it's simpler for processing
		// to save them all as separate entries. We just need to pull the right index.
		muon_entrypoint = const_cast<float*>(muboy_entrypoint[muboy_index].data());
		muon_direction = const_cast<float*>(muboy_dir.data());
		
		// muboy only returns one track length, but has multiple entry points
		// so i guess we need to calculate this for secondary muons as well...?
		if(muboy_index==0){
			muon_tracklen = muboy_tracklen;
		} else {
			muon_tracklen = CalculateTrackLen(muon_entrypoint, muon_direction);
		}
	}
	
	// Fill muon plots, but only the first time we read this muon
	// (XXX don't treat muboy tracks as unique for this purpose)
	if(muon_plotted.count(muHeader->nevsk)==0){
		muon_plotted.emplace(muHeader->nevsk,1);
	} else {
		return true;
	}
	mu_hb.Fill("mfmuselect_dir",muff_dir[0],muff_dir[1],muff_dir[2]);
	mu_hb.Fill("mfmuselect_entrypos", muff_entrypoint[0], muff_entrypoint[1], muff_entrypoint[2]);
	
	mu_hb.Fill("muboy_dir",muboy_dir[0],muboy_dir[1],muboy_dir[2]);
	mu_hb.Fill("muboy_entrypos", muboy_entrypoint[muboy_index][0],
	                             muboy_entrypoint[muboy_index][1],
	                             muboy_entrypoint[muboy_index][2]);
	mu_hb.Fill("muboy_goodness", muboy_goodness);
	mu_hb.Fill("muboy_class", int(muboy_status));
	mu_hb.Fill("muboy_ntracks", muboy_ntracks);
	mu_hb.Fill("muboy_tracklen", muboy_tracklen);
	
	if(didbff){
		mu_hb.Fill("bff_dir", bff_dir[0], bff_dir[1], bff_dir[2]);
		mu_hb.Fill("bff_entrypos", bff_entrypoint[0],
	                               bff_entrypoint[1],
	                               bff_entrypoint[2]);
		mu_hb.Fill("bff_dir", bff_goodness);
	}
	
	mu_hb.Fill("nevsk",muHeader->nevsk);
	mu_hb.Fill("qismsk", mu_qismsk);
	mu_hb.Fill("qismsk_corr", mu_qismsk_corr);
	mu_hb.Fill("NumIDhits", muTQReal->nhits);
	mu_hb.Fill("NumODhits", muTQAReal->nhits);
	mu_hb.Fill("muon_tracklen", muon_tracklen);
	
	int64_t thiseventticks = (muHeader->counter_32 & ~0x1FFFF);
	thiseventticks = (thiseventticks << 15);
	int64_t it0sk = (int64_t(muHeader->t0) & 0xFFFFFFFF);
	thiseventticks += it0sk;
	Log(m_unique_name+" this event at time "+std::to_string(double(thiseventticks/COUNT_PER_NSEC))+" ns since last rollover",v_debug,m_verbose);
	if(thiseventticks==0){
		std::cerr<<"MU EVENT AT 0 TICKS"<<std::endl
		         <<"counter_32: "<<muHeader->counter_32<<"; header dump: ";
		muHeader->Dump();
	}
	
	int64_t rolloverticks = muRollovers * (int64_t(1) << 47);
	Log(m_unique_name+" adding "+std::to_string(double(rolloverticks/COUNT_PER_NSEC))+" ns for "+std::to_string(muRollovers)+" rollovers",v_debug,m_verbose);
	double totalticks = thiseventticks + rolloverticks;
	Log(m_unique_name+" total ticks: "+std::to_string(totalticks)+" or "+std::to_string(double(totalticks/COUNT_PER_NSEC))+" ns",v_debug,m_verbose);
	Log(m_unique_name+" c.f. EventTicks: "+std::to_string(muClockTicks)+" or "+std::to_string(double(muClockTicks/COUNT_PER_NSEC))+" ns",v_debug,m_verbose);
	thiseventticks += muRollovers * (int64_t(1) << 47);
	mu_hb.Fill("event_time_secs",double(thiseventticks/COUNT_PER_NSEC)/1.E9);
	
	// record time between distinct muons
	// (not multiple muboy tracks or multiple muon subtriggers in the same readout)
	if(lastmu_nevsk==0){
		lastmu_nevsk = muHeader->nevsk;
		lastmuticks = thiseventticks;
	} else if(muHeader->nevsk != lastmu_nevsk){
		int64_t ticksDiff = thiseventticks - lastmuticks;
		if(ticksDiff<0) ticksDiff += (int64_t(1) << 47);
		if(ticksDiff==0){
				std::cerr<<"mu entry "<<muEntryNum<<" ticksDiff is 0 from last!"<<std::endl;
				assert(false);
				exit(-1);
				return false;
		}
		double ns_since_last = double(ticksDiff)/COUNT_PER_NSEC;
		if(ns_since_last<120E9){
			mu_hb.Fill("mu_to_mu_secs", ns_since_last/1E9);
			Log(m_unique_name+"mu_to_mu nanoseconds: "+std::to_string(ns_since_last)+" from nevsk "+std::to_string(muHeader->nevsk)
			   +" to "+std::to_string(lastmu_nevsk),v_debug,m_verbose);
			if((ns_since_last/1E9)>70E3){
				std::cerr<<"TOO MUCH!!"<<std::endl;
				assert(false);
				exit(-1);
			}
		} else {
			// we shouldn't get this, if our muon-relic matching window is +-60s.
			Log(m_unique_name+" Error! Muon to muon time of "+toString(ns_since_last/1E9)
			   +">120 seconds!",v_error,m_verbose);
		}
		lastmu_nevsk = muHeader->nevsk;
		lastmuticks = thiseventticks;
	}
	
	return get_ok;
}

double RelicMuonPlots::CalculateTrackLen(float* muon_entrypoint, float* muon_direction, double* exitpt){
	
	// for reference HITKTK is the water volume height and DITKTK is its diameter,
	// these are #defined constants in geotnkC.h
	
	// sanity check muon is inward going or it's not going to go through the tank
	if( ((-muon_entrypoint[0]*muon_direction[0] + -muon_entrypoint[1]*muon_direction[1])<0) &&
	    (muon_entrypoint[1]>=DITKTK/2.) ){
		// not radially inwards
		Log(m_unique_name+": Muon trajectory is not into tank!",v_error,m_verbose);
		return 0;
	}
	if( (muon_entrypoint[2] >= ( HITKTK/2.) && muon_direction[2]>0) ||
	    (muon_entrypoint[2] <= (-HITKTK/2.) && muon_direction[2]<0) ){
		Log(m_unique_name+": Muon track points out of endcaps!",v_error,m_verbose);
		// pointing out of barrel
		return 0;
	}
	
	// calculate track length under assumption of a through-going muon,
	// based on its entry point and direction
	// first check for muons directed in the x-y plane
	if(std::abs(muon_direction[2]) > 0.1){
		double dist_to_endcap = (HITKTK/2.) - std::abs(muon_entrypoint[2]);
		if(std::signbit(muon_entrypoint[2]) != std::signbit(muon_direction[2])){
			dist_to_endcap = (HITKTK - dist_to_endcap);
		}
		// start with the simple case; project to the plane of the appropriate endcap
		double dxdz = muon_direction[0]/std::abs(muon_direction[2]);
		double dydz = muon_direction[1]/std::abs(muon_direction[2]);
		double proj_x = muon_entrypoint[0] + dxdz*dist_to_endcap;
		double proj_y = muon_entrypoint[1] + dydz*dist_to_endcap;
		if(((proj_x*proj_x)+(proj_y*proj_y)) < std::pow(DITKTK/2.,2.)){
			// if the projected point is within the tank radius, this is where it exits the tank
			double xtravel = proj_x-muon_entrypoint[0];
			double ytravel = proj_y-muon_entrypoint[1];
			double tracklen = std::sqrt(std::pow(xtravel,2)+
				                        std::pow(ytravel,2.)+
				                        std::pow(dist_to_endcap,2.));
			
			if(exitpt){
				exitpt[0] = proj_x;
				exitpt[1] = proj_y;
				exitpt[2] = muon_entrypoint[2] + dist_to_endcap*(muon_direction[2]>0 ? 1 : -1);
			}
			return tracklen;
		}
	}
	// if the muon is in the x-y plane, or the projected point is outside the tank,
	// then the track left through the barrel
	
	// get the angle of the track in the x-y plane
	double trackangle = std::atan2(muon_direction[1],muon_direction[0]);
	
	// get angle the track makes from entry point to the centre of the tank
	double entrypointpolarangle = std::atan2(-muon_entrypoint[1],-muon_entrypoint[0]);
	
	double chordlen;
	if(std::sqrt(std::pow(muon_entrypoint[0],2.)+std::pow(muon_entrypoint[1],2.))==DITKTK/2.){
		// with these we can calculate the angle subtended by the chord the muon track makes with:
		double angle_subtended = M_PI - 2.*(trackangle - entrypointpolarangle);
		// and from this we can work out the length of the chord:
		chordlen = std::abs(DITKTK * std::sin(angle_subtended/2.));
	} else {
		// ah, but if the muon entered through an endcap, the chord is truncated
		double a = std::sqrt(std::pow(muon_entrypoint[0],2.)+std::pow(muon_entrypoint[1],2.));
		// a is the distance from entry point to centre of endcap
		double C = entrypointpolarangle - trackangle;
		// C is the angle from the trajectory to tank centre
		// from a/SinA = c/SinC; -> sinA = a/c*SinC
		double A = std::asin(a/(DITKTK/2.) * std::sin(C));
		// A is the opening angle of the chord
		chordlen = a*std::cos(C) + (DITKTK/2.)*std::cos(A);
	}
	
	// amount of z travel is then chord length times z gradient
	double dzdr = muon_direction[2]/std::sqrt(std::pow(muon_direction[0],2.)+
	                                          std::pow(muon_direction[1],2.));
	double zdist = dzdr * chordlen;
	// sum is track length
	double tracklen = std::sqrt(std::pow(zdist,2.)+std::pow(chordlen,2.));
	
	if(exitpt){
		exitpt[0] = muon_entrypoint[0] + (chordlen * std::cos(trackangle));
		exitpt[1] = muon_entrypoint[1] + (chordlen * std::sin(trackangle));
		exitpt[2] = muon_entrypoint[2] + zdist;
	}
	
	return tracklen;
	
}

bool RelicMuonPlots::MakePairVariables(){
	
	Log(m_unique_name+"making pair variables",v_debug,m_verbose);
	
	std::pair<int,int> pair{relicEvNum,muEvNum};
	if(pairs_compared.count(pair)!=0){
		Log(m_unique_name+" Error! Duplicate comparison between relic "+toString(relicEvNum)
		   +" and muon "+toString(muEvNum)+" for entries "+toString(relicEntryNum)
		   +" and "+toString(muEntryNum),v_error,m_verbose);
		return false;
	}
	pairs_compared.emplace(pair);
	
	mu_before_relic = (muEvNum < relicEvNum);
	mu_relic_evtnum_diff = relicEvNum - muEvNum;
	if(std::abs(mu_relic_evtnum_diff)>1E6){
		Log(m_unique_name+" Warning! Large event num difference of "+toString(mu_relic_evtnum_diff)
		   +" for relic "+toString(relicEvNum)+" and muon "+toString(muEvNum)+" in entries "
		   +toString(relicEntryNum)+" and "+toString(muEntryNum),v_error,m_verbose);
		//return false;
	}
	
	Log(m_unique_name+"relic nevsk: "+std::to_string(relicEvNum)+", muon nevsk: "+std::to_string(muEvNum)
	         +"; spall cand? "+(mu_before_relic ? "Y" : "N"),v_debug,m_verbose);
	
	basic_array<float> scott_dedx(muMu->muboy_dedx);
	//float (*kirk_dedx_arr)[200] = (float(*)[200])(muMu->muinfo+10);
	float (&kirk_dedx_arr)[200] = *(float(*)[200])(&muMu->muinfo[10]);
	basic_array<float> kirk_dedx(kirk_dedx_arr);
	
	
	// debug check
	/*
	std::cout<<"scott's dedx: size "<<scott_dedx.size()<<" at "<<scott_dedx.data()<<std::endl;
	for(int i=0; i<std::min(200,scott_dedx.size()); ++i){
		if(i>0) std::cout<<", ";
		std::cout<<scott_dedx.at(i);
	}
	std::cout<<"}"<<std::endl;
	*/
	/*
	std::cout<<"kirk's dedx: size "<<kirk_dedx.size()<<" at "<<kirk_dedx.data()<<std::endl;
	for(int i=0; i<std::min(200,kirk_dedx.size()); ++i){
		if(i>0) std::cout<<", ";
		std::cout<<kirk_dedx.at(i);
	}
	std::cout<<"}"<<std::endl;
	*/
	
	// seems like scotts is the latest
	muon_dedx = const_cast<float*>(scott_dedx.data());
	//muon_dedx = const_cast<float*>(kirk_dedx.data());
	
	// find position of max dedx
	// defined as bin where a sliding window of 4.5m (9x 50cm bins) has maximum sum
	double max_edep = 0;
	int max_edep_bin=0;
	for(int i=0;i<111;i++){   // from mu_info.C - XXX why 111? muboy_dedx array is 200 bins in length...?
		// calculate sum of dedx in 9 bins from bin i to i+9
		// (use 9 as an odd number so we have one central bin)
		double e_dep_in_window = 0 ;
		for(int j=0;j<9;j++){
			e_dep_in_window = e_dep_in_window + muon_dedx[i+j]; // will be qpeak
		}
		if(e_dep_in_window > max_edep){
			max_edep_bin = i+4;  // +4 to get centre bin
			max_edep = e_dep_in_window;
		}
	}
	double max_edep_pos = 50.*max_edep_bin;
	Log(m_unique_name+"muon track max dE/dx: "+toString(max_edep)+" at "+toString(max_edep_pos)+" cm "
	         +"(bin "+toString(max_edep_bin)+"/115) along the track",v_debug,m_verbose);
	
	// re-check the coulomb-to-photoelectron conversion factor (petable entry)
	double pe_per_coulomb = 0;
	// scan through pe table for this run, and get corresponding pe_per_coulomb
	Log(m_unique_name+"scanning petable for pe_to_coulombs conversion",v_debug,m_verbose);
	while(pe_per_coulomb == 0 && pe_table_index < petable_startrun.size()){
		Log(m_unique_name+"entry "+std::to_string(pe_table_index),v_debug,m_verbose);
		if(petable_startrun[pe_table_index] <= muHeader->nrunsk && petable_endrun[pe_table_index] >= muHeader->nrunsk){
		    pe_per_coulomb = pe_to_coulombs[pe_table_index][1];
		    break;
		}
		pe_table_index++;
	}
	if(pe_table_index==petable_startrun.size()){
		Log(m_unique_name+" Error! Did not find run "+toString(muHeader->nrunsk)+" in petable!",v_error,m_verbose);
		// use a nominal value?
		pe_per_coulomb = 30;
	}
	
	// XXX not sure why pe_per_cm is in this formula?
	Log(m_unique_name+" calculating observed pe from this muon",v_debug,m_verbose);
	double pe_from_muon = muMu->muqismsk * (pe_per_cm / pe_per_coulomb);  // pe*cm^-1 / pe*C^-1 = C/cm
	// mu_info.C calculates this based on tracklen from muboy (resq_tmp), tracklen calculated
	// based on muboy secondary track entry point and muboy dir (resq_sprt_tmp, for mutiple muons)
	// and calculated tracklen BFF track (resq_sprt_tmp_bff)
	double pe_from_MIP = muon_tracklen * pe_per_cm;
	double residual_q_over_MIP = pe_from_muon - pe_from_MIP;
	Log(m_unique_name+" pe from MIP: "+std::to_string(pe_from_MIP)+", pe from muon: "+std::to_string(pe_from_muon)
	    +", residual: "+std::to_string(residual_q_over_MIP),v_debug,m_verbose);
	
	// calculate footpoint and transverse distance
	// find transverse distance from relic to muon
	float* relic_vertex = relicLowe->bsvertex;
	float appr;  // output: XXX is this "foot point"? distance along track where dlt is defined relative to?
	
	/*
	std::cout<<"calling getdl_ with:\n"
	         <<"\trelic vertex: ("<<relic_vertex[0]<<", "<<relic_vertex[1]<<", "<<relic_vertex[2]<<")\n"
	         <<"\tmuon entry point: ("<<muon_entrypoint[0]<<", "<<muon_entrypoint[1]<<", "<<muon_entrypoint[2]<<")\n"
	         <<"\tmuon entry dir: ("<<muon_direction[0]<<", "<<muon_direction[1]<<", "<<muon_direction[2]<<")"<<std::endl;
	*/
	
	getdl_(muon_direction, &relic_vertex[0], &relic_vertex[1], &relic_vertex[2], muon_entrypoint, &dlt, &appr);
	
	Log(m_unique_name+" transverse distance: "+std::to_string(dlt)+", with foot point at "+std::to_string(appr)
	    +" cm along muon track",v_debug,m_verbose);
	// XXX sometimes, but not always, appr is negative. Might be fine. Might want to check.
	
	// XXX warning: see Shinoki-san's talk from lowe meeting ~06/05/21:
	// lt distriubtion for multiple muons with low goodness (<~0.4) is quite large,
	// but cutting on muboy goodness can bias muon energy and neutron multiplicity (by affecting search window)
	// he evaluates systematics of this cut; we should look into this.
	// see his slides also for example distributions
	
	// calculate longitudinal distance from point of peak energy deposition
	dll = max_edep_pos - appr;
	Log(m_unique_name+" longitudinal distance: "+std::to_string(dll),v_debug,m_verbose);
	
	
	return true;
}


bool RelicMuonPlots::GetPeTable(){
	// read pe table file which describes conversion from coulombs to photo-electrons
	Log(m_unique_name+" Filling pe table from "+pe_table_filename,v_debug,m_verbose);
	std::ifstream pe_table_file(pe_table_filename);
	
	// first count the lines in the file so we can pre-allocate memory
	size_t numlines = std::count(std::istreambuf_iterator<char>(pe_table_file), 
	                             std::istreambuf_iterator<char>(), '\n');
	
	// Petable lookup variables
	std::string nextline;
	petable_startrun.resize(numlines);
	petable_endrun.resize(numlines);
	pe_to_coulombs.resize(numlines);
	
	// scan from the start of the file until we find the line for the current run number
	int i = 0;
	while(getline(pe_table_file, nextline)){
		std::stringstream ss(nextline);
		ss >> petable_startrun[i]
		   >> petable_endrun[i]
		   >> pe_to_coulombs[i][0]
		   >> pe_to_coulombs[i][1]
		   >> pe_to_coulombs[i][2];
		i++;
	}
	
	pe_table_file.close();
	
	// each entry covers a given span of run numbers. As we change runs we may need to
	// move to the next petable entry. Keep track of our position in the table.
	pe_table_index = 0;
	
	return true;
}

bool RelicMuonPlots::MakeHists(int step){
	
	get_ok = true;
	
	// step 0: Initialise
	// ==================
	if(step==0){
		hb.MakeFile(outputFile);
		hb.SaveHists(false);
		mu_hb.MakeFile(muFile);
		mu_hb.SaveHists(false);
		relic_hb.MakeFile(relicFile);
		relic_hb.SaveHists(false);
		
		// sanity check?
		if(hb.GetFile()==nullptr){
			Log(m_unique_name+" histogram builder 'hb' has no output file!", v_error,m_verbose);
			return false;
		}
		
	// step 1: Fill
	// ============
	} else if(step==1){
		// fill tree
		hb.Fill("mu_before_relic", mu_before_relic);
		hb.Fill("evtNumDiff",mu_relic_evtnum_diff);
		hb.Fill("dt", dt);                   // time diff
		hb.Fill("dll", dll);                 // longitudinal distance from relic to point of muon max dedx
		hb.Fill("dlt", dlt);                 // transverse distance from relic to point of muon max dedx
		
		if(mu_before_relic==0){
			spall_count++;
		} else {
			rand_count++;
		}
		
	// step 2: Finalise
	// ================
	} else {
		// names should be branch name + '_spall' or '_relic'
		std::vector<TH1*> spall_dists;
		std::vector<TH1*> random_dists;
		
		// make/get histograms
		spall_dists.push_back(hb.GetHist("dt","mu_before_relic==1"));
		random_dists.push_back(hb.GetHist("dt","mu_before_relic==0"));
		
		spall_dists.push_back(hb.GetHist("dll","mu_before_relic==1"));
		random_dists.push_back(hb.GetHist("dll","mu_before_relic==0"));
		
		spall_dists.push_back(hb.GetHist("dlt","mu_before_relic==1"));
		random_dists.push_back(hb.GetHist("dlt","mu_before_relic==0"));
		
		spall_dists.push_back(hb.GetHist("evtNumDiff","mu_before_relic==1"));
		random_dists.push_back(hb.GetHist("evtNumDiff","mu_before_relic==0"));
		
		// draw and save distributions
		hb.GetFile()->cd();
		TCanvas* c_spall = new TCanvas("c_spall","c_spall",1280,1024);
		c_spall->cd();
		for(int i=0; i<spall_dists.size(); ++i){
			TH1* spall_dist = spall_dists.at(i);
			TH1* rand_dist = random_dists.at(i);
			std::string paramname = spall_dist->GetName();
			paramname = paramname.substr(0,paramname.find('_'));
			spall_dist->SetLineColor(kRed);
			rand_dist->SetLineColor(kBlack);
			
			spall_dist->Write(std::string(paramname+"_spall").c_str());
			rand_dist->Write(std::string(paramname+"_rand").c_str());
			
			c_spall->Clear();
			// use a stack so axes ranges don't get clipped
			THStack hs(paramname.c_str(), paramname.c_str());
			// clone because the HistogramBuilder owns its histograms...
			// (technically the file it creates takes ownership of them...)
			hs.Add((TH1*)(spall_dist->Clone()));
			hs.Add((TH1*)(rand_dist->Clone()));
			hs.Draw("nostack");
			c_spall->BuildLegend();
			c_spall->Write(paramname.c_str());
		}
		
	}
	
	return get_ok;
}

