#include "PlotHitTimes.h"
#include "MTreeReader.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "ColourWheel.h"
#include "Constants.h"

PlotHitTimes::PlotHitTimes():Tool(){}

bool PlotHitTimes::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	std::string foutname = ""; //"subtriggers.root";
	m_variables.Get("outputfile",foutname);
	std::string treeReaderName;
	m_variables.Get("treeReaderName",treeReaderName);
	m_variables.Get("onlyWriteSubtriggers",onlyWriteSubtriggers);
	
	// get the reader for inputs
	if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	if(!foutname.empty()){
		fout = new TFile(foutname.c_str(), "RECREATE");
		c_subtriggers = new TCanvas("c_subtriggers","c_subtriggers",1024,800);
		gDirectory->cd();
	}

        useSLESearchTool=false;
        m_variables.Get("useSLESearchTool",useSLESearchTool);
	
	// how to find subtriggers - SLESearch Tool or get_sub_triggers
	useSLESearchTool=false;
	m_variables.Get("useSLESearchTool",useSLESearchTool);
	
	return true;
}


bool PlotHitTimes::Execute(){
	
	Log(m_unique_name+" Executing...",v_debug,m_verbose);
	
	if(sktqz_.nqiskz==0){
		Log(m_unique_name+": sktqz had no ID hits",v_message,m_verbose);
		return true;
	}
	
	std::cout<<"Primary trigger was of type "<<GetTriggerNames(skhead_.idtgsk)<<std::endl;
	
	// use sktqz as it seems to be populated either from TQReal or TQList,
	// covering frm file data, processed data, and MC files
	// debug print that we do indeed have hits
	std::cout<<"nqiskz: "<<sktqz_.nqiskz<<std::endl;
	int pnhits=std::min(5,sktqz_.nqiskz);
	for(int i=0; i<pnhits; ++i) std::cout<<"\thit at "<<sktqz_.tiskz[i]<<std::endl;
	
	// since set_timing_gate (called by GetSubtriggerFlags) updates the times in sktqz_
	// to be relative to the new t0 for each subtrigger, we'll need to make a note of those
	// times now (or we could call set_timing_gate_(t0sk) again before filling the histograms) << we do this, verfieid its ok
	///std::vector<float> times_wrt_primary_trigger(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz);
	
	// use a bitset for each hit to record which subtrigger 1.3us gate (if any) that hit is in
	// a 32-bit bitset allows us to note up to 31 subtriggers - should be plenty
	std::vector<std::bitset<32>> in_subtrigger_flags(sktqz_.nqiskz);
	
	// first note which hits are in the primary trigger 1.3us window
	int n_hits_in_prim_13=0;
	for(int i=0; i<sktqz_.nqiskz; ++i){
		std::unordered_map<std::string, int> hit_flags = GetHitFlagNames(sktqz_.ihtiflz[i]);
		if(hit_flags.at("in 1.3us")){
			in_subtrigger_flags[i].set(0);
			++n_hits_in_prim_13;
		}
	}
	Log(m_unique_name+": had "+toString(n_hits_in_prim_13)+" hits in primary 1.3us window",v_debug,m_verbose);
	
	// scan for subtriggers of given type(s) and note which hits are in each
	int n_triggers = 1;
	std::vector<int> t0_sub(1,0);  // times of subtriggers
	std::vector<TriggerType> triggers_of_interest{TriggerType::SLE}; //,LE,HE,SHE,AFT,OD};
	for(TriggerType atrigtype : triggers_of_interest){
		std::cout<<"\tSearching for " << TriggerIDToName(atrigtype)<<" triggers\n"
		         << "\t\tthreshold from runinfo: " << skruninf_.softtrg_thr[atrigtype]<<"\n";
		// check we have a threshold set for this trigger type
		if(skruninf_.softtrg_thr[atrigtype]==0){
			int default_thresh = GetTriggerThreshold(atrigtype);
			std::cout<<"\t\t->overriding with default "<<default_thresh<<"\n";
			skruninf_.softtrg_thr[atrigtype] = default_thresh;
		}
		
		// find the subtriggers
		n_triggers += GetSubtriggers(atrigtype, t0_sub);
	}
	
	// note what hits are in that subtrigger
	GetSubtriggerFlags(t0_sub, in_subtrigger_flags);
	
//	// XXX DEBUG XXX
//	if(n_triggers>1){
//		m_data->vars.Set("StopLoop",1);
//	}
	
	// we're gonna make a bunch of temporary histograms on each Execute, re-using the same names.
	// tell ROOT not to assume ownership (i.e. not to add them to its list of histograms),
	// otherwise ROOT "remembers" that a histogram of that name existed (even after it goes out of scope)
	// and spits out a warning that it is being replaced and there may be a memory leak (which there isn't).
	//TH1::AddDirectory(kFALSE);
	// this is a blunt tool as it's a global directive; instead manually remove them from gDirectory later
	
	// subsequent calls will re-read from sktqz_.itisk_raw and tbuf_raw, which are unaffected by set_timing_gate
	// so we can repeatedly call it to scan over each subtrigger
	
	// make a plot with each subtrigger selected as primary trigger
	// these should all look the same except with shifted time axes
	// ~and the primary trigger colour always at ~t=0~ << we should do this, but we don't yet
	for(int m=0; m<n_triggers; ++m){
		
		// adjust times relative to this subtrigger
		int it0xsk = skheadqb_.it0sk + t0_sub.at(m);
		if(myTreeReader->GetTree()->GetBranch("TQREAL")==nullptr){
			set_timing_gate_(&it0xsk);
		} else {
			set_timing_gate_m_(&it0xsk);
		}
		/*
		// sanity check that loading subtriggers then reloading the primary does not affect times
		// passed!
		if(m==0){
			bool allsame=true;
			for(size_t j=0; j<sktqz_.nqiskz; ++j){
				if(times_wrt_primary_trigger[j]!=sktqz_.tiskz[j]){
					std::cerr<<"reloading primary trigger changed hit "<<j<<" time from "
					         <<times_wrt_primary_trigger[j]<<" to "<<sktqz_.tiskz[j]<<std::endl;
					allsame=false;
				}
			}
			if(allsame){
				std::cout<<"Reloading primary trigger times OK!"<<std::endl;
			}
		}
		*/
		
		float earliest_hit_time = *std::min_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz);
		float latest_hit_time = *std::max_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz);
		float padding = (latest_hit_time - earliest_hit_time)*0.05; // 5% padding on range
		Log(m_unique_name+" hits span range "+toString(earliest_hit_time)+" to "
		    +toString(latest_hit_time),v_debug,m_verbose);
		
		//int nbins=1000;
		// use number of bins such that each bin represents 200ns, the window size used for finding subtriggers
		int nbins=padding*0.105;
		
		// histograms for out-of-gate, in-gate, and each 1.3us window for the (current) primary and all subtriggers
		// times will be plotted in the frame of reference of the primary trigger
		THStack h_hits("h_hits","Hit Times");
		h_hits.Add(new TH1D("h_in_gate","In-Gate Hit Times", nbins, earliest_hit_time-padding, latest_hit_time+padding));
		h_hits.Add(new TH1D("h_out_gate","Out-of-Gate Hit Times", nbins, earliest_hit_time-padding, latest_hit_time+padding));
		h_hits.Add(new TH1D("h_in_13us_0","Hit Times in 1.3us gate 0", nbins, earliest_hit_time-padding, latest_hit_time+padding));
		for(int i=0; i<n_triggers; ++i){
			std::string h_name = "h_in_13us_"+std::to_string(i+1);
			std::string h_title = "Hit Times in 1.3us gate "+std::to_string(i+1);
			h_hits.Add(new TH1D(h_name.c_str(), h_title.c_str(), nbins, earliest_hit_time-padding, latest_hit_time+padding));
		}
		
		// just for giggles lets also plot the distribution of charges for these populations
		// not sure of the histogram range, we have the max charge of any hit in the 1.3us window in skq_.qimxsk
		// which is probably the max charge of any hit in the readout, but let's double it as a precaution...
		TH1D h_in_gate_q("h_in_gate_q","In-Gate Hit Charges", nbins, 0, 2.0*skq_.qimxsk);
		TH1D h_out_gate_q("h_out_gate_q","Out-of-Gate Hit Charges", nbins, 0, 2.0*skq_.qimxsk);
		TH1D h_in_primary_1p3us_q("h_in_primary_13us_q","Hit Charges In Primary 1.3us", nbins, 0, 2.0*skq_.qimxsk);
		
		int multi_subtrigger_hits=4; // debug
		std::vector<int> hitcounts_in_1p3us(n_triggers); // debug
		TH1D* ahist = nullptr; // placeholder
		
		for(int i=0; i<sktqz_.nqiskz; ++i){
			std::unordered_map<std::string, int> hit_flags = GetHitFlagNames(sktqz_.ihtiflz[i]);
			
			// see if it's in one of the 1.3us windows
			if(in_subtrigger_flags.at(i).any()){
				// it's in one of the subtriggers
				
				// sanity check to see if this hit is in more than one subtrigger: they are not supposed to overlap.
				// (except of course with the primary trigger, which is independent from the subtrigger search)
				int in_prim = (in_subtrigger_flags.at(i).test(0)) ? 1 : 0;
				if(multi_subtrigger_hits>0 && in_subtrigger_flags.at(i).count()>(in_prim+1)){
					std::cerr<<"Hit "<<i<<" is in "<<(in_subtrigger_flags.at(i).count()-in_prim)
					         <<" subtriggers!"<<std::endl;
					if(multi_subtrigger_hits==0) std::cerr<<"suppressing further warnings..."<<std::endl;
					--multi_subtrigger_hits;
				}
				// add it to the appropriate histogram. prioritise the primary trigger
				// otherwise take the first subtrigger if multiple (it doesn't seem to happen)
				for(int j=0; j<32; ++j){
					if(in_subtrigger_flags.at(i).test(j)){
						ahist = (TH1D*)h_hits.GetHists()->At(2+j);
						ahist->Fill(sktqz_.tiskz[i]);
						++hitcounts_in_1p3us.at(j);
						if(j==0) h_in_primary_1p3us_q.Fill(sktqz_.qiskz[i]);
						break;
					}
				}
			// else see if it's in the primary trigger gate (full primary trigger not just 1.3us)
			} else if(hit_flags.at("in gate")){
				ahist = (TH1D*)h_hits.GetHists()->At(0);
				ahist->Fill(sktqz_.tiskz[i]);
				h_in_gate_q.Fill(sktqz_.qiskz[i]);
			// else it's an out-of-gate hit
			} else {
				ahist = (TH1D*)h_hits.GetHists()->At(1);
				ahist->Fill(sktqz_.tiskz[i]);
				h_out_gate_q.Fill(sktqz_.qiskz[i]);
			}
		}
		
		// debug print
		std::cout<<"1.3us window hit counts:\n";
		for(int i=0; i<hitcounts_in_1p3us.size(); ++i){
			std::cout<<"\tSubtrigger "<<i<<": "<<hitcounts_in_1p3us.at(i);
			ahist=(TH1D*)h_hits.GetHists()->At(i+2);
			std::cout<<" ("<<ahist->GetEntries()<<")\n";
		}
		if(hitcounts_in_1p3us.size()==0){
			std::cout<<"\tNone"<<std::endl;
		}
		
		// FIXME right now a subtrigger may wholly overlap with the primary trigger,
		// which will result in an empty subtrigger. Maybe we can identify and prune such cases.
		
		// we combine all our histograms into one plot to produce one histogram of all hits,
		// but since we have each group as a separate TH1 we can give them each different fill colours
		
		// in-gate but not in 1.3 window hits just fill with white (i.e. leave as default)
		// out-of-gate hits fill with black
		ahist = (TH1D*)h_hits.GetHists()->At(1);
		ahist->SetFillStyle(1001);
		ahist->SetFillColor(1);
		
		// subtriggers we'll shade with varying colours
		colourwheel.Reset();
		for(int i=2; i<h_hits.GetHists()->GetEntries(); ++i){
			ahist = (TH1D*)h_hits.GetHists()->At(i);
			ahist->SetFillStyle(3002);
			ahist->SetFillColor(colourwheel.GetNextColour());
			//std::cout<<"setting trigger "<<i<<" to colour "<<colourwheel.GetCurrentColour()<<std::endl;
		}
		
		if(fout && (!onlyWriteSubtriggers || n_triggers>1)){
			fout->cd();
			c_subtriggers->cd();
			h_hits.Draw();
			std::string canv_name = "event_"+std::to_string(myTreeReader->GetEntryNumber())+"_trigger_"+std::to_string(m);
			std::cout<<"writing event "<<canv_name<<std::endl;
			c_subtriggers->Write(canv_name.c_str());
			gDirectory->cd();
		}
		
		// NO, ROOT, YOU DO NOT OWN THESE. THEY DON'T EVEN EXIST ANY MORE. STOP IT.
		for(int i=0; i<h_hits.GetHists()->GetEntries(); ++i){
			gDirectory->GetList()->Remove(h_hits.GetHists()->At(i));
		}
		c_subtriggers->Clear();
		
	}
	
	return true;
}

int PlotHitTimes::GetSubtriggers(int subtrigtype, std::vector<int>& subtrigger_times){
	
	// scan for subtriggers of a given type
	
	// arbitrary limit of max num subtriggers. 32 should be more than sufficient
	// note this has to be equal to or smaller than the size of the bitset we use
	int MAX_SUBTRIGS=32;
	
	// trigger time of the primary trigger
	std::cout<<"primary trigger is at "<<skheadqb_.it0sk
	         <<" or "<<(double(skheadqb_.it0sk)/COUNT_PER_NSEC)<<" ns"<<std::endl;
	std::cout<<"t0 offset is "<<skruninf_.softtrg_t0_offset[TriggerType::SLE_hitsum]<<" or "
	         <<double(skruninf_.softtrg_t0_offset[TriggerType::SLE_hitsum])/COUNT_PER_NSEC<<" ns"<<std::endl;
	
	int ntrigsfound=0;
	std::vector<int> t0_sub(MAX_SUBTRIGS,-1);  // relative time of subtrigger to IT0SK
	
	// run subtrigger algorithm to search for subtriggers of this type
	if(useSLESearchTool){
		std::vector<double> sle_times;
		get_ok = m_data->CStore.Get("SLE_times",sle_times); // these are in ns from primary trigger
		if(!get_ok){
			Log(m_unique_name+" no SLE_times in CStore!",v_error,m_verbose);
			return false;
		}
		// convert nanoseconds from primary trigger to clock ticks
		for(int i=0; i<sle_times.size(); ++i){
			t0_sub.at(i) = sle_times.at(i)*COUNT_PER_NSEC;
			std::cout<<"SLESearch time "<<i<<" at "<<sle_times.at(i)
			         << " ns or " << t0_sub.at(i) << " ticks from primary trigger" << std::endl;
		}
		ntrigsfound=sle_times.size();
	} else {
		// run subtrigger algorithm to search for subtriggers of this type
		get_sub_triggers_(&subtrigtype, &ntrigsfound, t0_sub.data(), &MAX_SUBTRIGS);
	}
	
	Log(m_unique_name+" found "+toString(ntrigsfound)+" subtriggers of type "
	    +TriggerIDToName(subtrigtype),v_message,m_verbose);
	
	subtrigger_times.insert(subtrigger_times.end(), t0_sub.begin(), t0_sub.begin()+ntrigsfound);
	
	return ntrigsfound;
	
}

int PlotHitTimes::GetSubtriggerFlags(const std::vector<int>& t0_sub, std::vector<std::bitset<32>>& subtrigger_flags){
	
	// we need the lun to call skcread
	// this is needed to recalculate some variables, e.g. nqisk
	// note this should not be needed to reload hits - times and flags are updated by set_timing_gate
	int lun = m_data->GetLUN(myTreeReader->GetName());
	lun = -std::abs(lun);  // IMPORTANT: make it negative to ensure skcread does not advance TTree
	
	// for each subtrigger, scan over every hit and note which subtriggers are in its 1.3us gate
	int nhiterrs=5;
	for(int i=0; i<t0_sub.size(); ++i){
		
		// trigger time of the subtrigger
		int it0xsk = skheadqb_.it0sk + t0_sub.at(i);
		
         // FIXME a few of the following assume SLE_hitsum triggertype
         // would need to pass the trigger types vector if we need this
         // for now though nothing actually critical, only printouts and debug stuff uses it.
		std::cout<<"subtrigger "<<i<<" is at it0xsk "<<uint64_t(it0xsk)
		         <<" or "<<(t0_sub.at(i)/COUNT_PER_NSEC)<<" ns from the primary trigger"
		         <<" so peak should be at "<<(double(t0_sub.at(i))/COUNT_PER_NSEC)
		         +(double(skruninf_.softtrg_t0_offset[TriggerType::SLE_hitsum])/COUNT_PER_NSEC)<<" ns"<<std::endl;
		
		// set IT0XSK to the position of the next subtrigger
		// this updates the hit flags ihtiflz and icabiz, anc times in itiskz, tiskz
		// based on the values in commons itiskz_raw, tbuf_raw, icabbf_raw etc.
		if(myTreeReader->GetTree()->GetBranch("TQREAL")==nullptr){
			set_timing_gate_(&it0xsk);
		} else {
			// the in-1.3us flag is set based on whether sktqz_.itisk[] is within the gate
			// but itiskz is not populated for processed files (tbuf_raw and tiskz are, but not itiskz)
			// so we have a modified version that instead uses these instead
			Log(m_unique_name+" using set_timing_gate_m",v_debug,m_verbose);
			set_timing_gate_m_(&it0xsk);
		}
		
		/*
		// get in-gate times (for check?)
		int pre_t0 = skruninf_.softtrg_pre_t0[TriggerType::SLE_hitsum];  // FIXME is this the right index to use?
		int post_t0 = skruninf_.softtrg_post_t0[TriggerType::SLE_hitsum]; // FIXME is this the right index to use?
		// n.b. the signal is expected to be at it0xsk + skruninf_.sofftrg_t0_offset, NOT at it0xsk itself.
		float subtr_start_t = float(it0xsk - pre_t0) / COUNT_PER_NSEC;
		float subtr_end_t = float(it0xsk + post_t0) / COUNT_PER_NSEC;
		*/
		
		int n_in_gate_hits_this_subtrigger=0;
		
		// scan hits, check whether they're in the new 1.3us gate
		std::cout<<"looping over "<<subtrigger_flags.size()<<" hits"<<std::endl;
		for(int k=0; k<subtrigger_flags.size(); ++k){
			std::unordered_map<std::string, int> hit_flags = GetHitFlagNames(sktqz_.ihtiflz[k]);
			
			// sanity check
			/* FIXME enable when we know what index to use for pre- and post- t0's
			// fixed; pretty sure it's runinfo_.t0_offset[i] where i is the index for the type of trigger
			// as in GetTriggerName(i)
			float new_time = sktqz_.tiskz[i] - (it0xsk - skheadqb_.it0sk)/COUNT_PER_NSEC;
			if(nhiterrs>0 && hit_flags.at("in 1.3us") != (new_time > subtr_start_t && new_time < subtr_end_t)){
				Log(m_unique_name+" hit "+toString(k)+" in subtrigger "+toString(j)+" has time "
				   +toString(new_time)+" relative to subtrigger window ["+toString(subtr_start_t)
				   +","+toString(subtr_end_t)+"], but in-gate flag "+toString(hit_flags.at("1.3us")),
				   v_error,m_verbose);
				--nhiterrs;  // limit to the first N occurrences so we don't swamp the user
			}
			*/
			
			if(hit_flags.at("in 1.3us")){
				subtrigger_flags.at(k).set(i);
				++n_in_gate_hits_this_subtrigger;
			} else {
				subtrigger_flags.at(k).reset(i);
			}
		}
		
		std::cout<<"Found "<<n_in_gate_hits_this_subtrigger
		         <<" hits with in-1.3us flag set for subtrigger "<<i<<std::endl;
		
	}
	
	return 1;
}


bool PlotHitTimes::Finalise(){
	
	if(fout){
		fout->Write("",TObject::kOverwrite);
		fout->Close();
		delete fout;
	}
	delete c_subtriggers;
	return true;
}
