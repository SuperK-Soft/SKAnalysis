#include "PlotHitTimes.h"
#include "MTreeReader.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "ColourWheel.h"
#include "Constants.h"

PlotHitTimes::PlotHitTimes():Tool(){
	// get the name of the tool from its class name
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}

bool PlotHitTimes::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	std::string foutname = "subtriggers.root";
	m_variables.Get("outputfile",foutname);
	
	// get the reader for inputs
	std::string treeReaderName;
	m_variables.Get("treeReaderName",treeReaderName);
	if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	fout = new TFile(foutname.c_str(), "RECREATE");
	c_subtriggers = new TCanvas("c_subtriggers","c_subtriggers",1024,800);
	
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
	float earliest_hit_time = *std::min_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz);
	float latest_hit_time = *std::max_element(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz);
	float padding = (latest_hit_time - earliest_hit_time)*0.05; // 5% padding on range
	
	// since set_timing_gate (called by GetSubtriggerFlags) updates the times in sktqz_
	// to be relative to the new t0 for each subtrigger, we'll need to make a note of those
	// times now (or we could call set_timing_gate_(t0sk) again before filling the histograms)
	std::vector<float> times_wrt_primary_trigger(sktqz_.tiskz, sktqz_.tiskz+sktqz_.nqiskz);
	
	// use a bitset for each hit to note which 1.3us windows its in
	// a 32-bit bitset allows us to note up to 31 subtriggers - should be plenty
	std::vector<std::bitset<32>> in_subtrigger_flags(sktqz_.nqiskz);
	
	// first note which hits are in the primary trigger 1.3us window
	int n_triggers = 1;
	int n_hits_in_prim_13=0;
	for(int i=0; i<sktqz_.nqiskz; ++i){
		std::unordered_map<std::string, int> hit_flags = GetHitFlagNames(sktqz_.ihtiflz[i]);
		if(hit_flags.at("in 1.3us")){
			in_subtrigger_flags.at(i).set(0);
			++n_hits_in_prim_13;
		}
	}
	//std::cout<<"had "<<n_hits_in_prim_13<<" hits in primary 1.3us window"<<std::endl;
	
	// then scan for subtriggers of given type(s) and note which hits are in those
	std::vector<TriggerType> triggers_of_interest{TriggerType::SLE}; //,LE,HE,SHE,AFT,OD};
	for(TriggerType atrigtype : triggers_of_interest){
		int thetrigtype = (int)(atrigtype);
		n_triggers += GetSubtriggerFlags(thetrigtype, in_subtrigger_flags, n_triggers);
	}
	
	// histograms for out-of-gate, in-gate, and each 1.3us window for the primary and all subtriggers
	// times will be plotted in the frame of reference of the primary trigger
	THStack h_hits("h_hits","Hit Times");
	int nbins=1000;
	h_hits.Add(new TH1D("h_in_gate","In-Gate Hit Times", nbins, earliest_hit_time-padding, latest_hit_time+padding));
	h_hits.Add(new TH1D("h_out_gate","Out-of-Gate Hit Times", nbins, earliest_hit_time-padding, latest_hit_time+padding));
	h_hits.Add(new TH1D("h_in_13us_0","Hit Times in 1.3us gate 0", nbins, earliest_hit_time-padding, latest_hit_time+padding));
	for(int i=0; i<n_triggers; ++i){
		std::string h_name = "h_in_13us_"+std::to_string(i+1);
		std::string h_title = "Hit Times in 1.3us gate "+std::to_string(i+1);
		h_hits.Add(new TH1D(h_name.c_str(), h_title.c_str(), nbins, earliest_hit_time-padding, latest_hit_time+padding));
	}
	TH1D* ahist = nullptr; // placeholder
	
	// just for giggles lets also plot the distribution of charges for these populations
	// not sure of the histogram range, we have the max charge of any hit in the 1.3us window in skq_.qimxsk
	// which is probably the max charge of any hit in the readout, but let's double it as a precaution...
	TH1D h_in_gate_q("h_in_gate_q","In-Gate Hit Charges", nbins, 0, 2.0*skq_.qimxsk);
	TH1D h_out_gate_q("h_out_gate_q","Out-of-Gate Hit Charges", nbins, 0, 2.0*skq_.qimxsk);
	TH1D h_in_primary_1p3us_q("h_in_primary_13us_q","Hit Charges In Primary 1.3us", nbins, 0, 2.0*skq_.qimxsk);
	
	int multi_subtrigger_hits=4; // debug
	std::vector<int> hitcounts_in_1p3us(n_triggers); // debug
	
	for(int i=0; i<sktqz_.nqiskz; ++i){
		std::unordered_map<std::string, int> hit_flags = GetHitFlagNames(sktqz_.ihtiflz[i]);
		
		// see if it's in one of the 1.3us windows
		if(in_subtrigger_flags.at(i).any()){
			// it's in one of the subtriggers; add it to the appropriate histogram
			// sanity check if it's in more than one: I believe they are not supposed to overlap.
			// (this doesn't seem to happen, so that seems correct)
			if(multi_subtrigger_hits>0 && in_subtrigger_flags.at(i).count()>1){
				std::cerr<<"Hit "<<i<<" is in "<<in_subtrigger_flags.at(i).count()<<" subtriggers!"<<std::endl;
				--multi_subtrigger_hits;
			}
			// i guess take the first trigger we encounter if this happens (it doesn't seem to anyway)
			for(int j=0; j<32; ++j){
				if(in_subtrigger_flags.at(i).test(j)){
					ahist = (TH1D*)h_hits.GetHists()->At(2+j);
					ahist->Fill(times_wrt_primary_trigger[i]);
					if(j==0) h_in_primary_1p3us_q.Fill(sktqz_.qiskz[i]);
					++hitcounts_in_1p3us.at(j);
					break;
				}
			}
		// else see if it's in the primary trigger gate
		} else if(hit_flags.at("in gate")){
			ahist = (TH1D*)h_hits.GetHists()->At(0);
			ahist->Fill(times_wrt_primary_trigger[i]);
			h_in_gate_q.Fill(sktqz_.qiskz[i]);
		// else it's an out-of-gate hit
		} else {
			ahist = (TH1D*)h_hits.GetHists()->At(1);
			ahist->Fill(times_wrt_primary_trigger[i]);
			h_out_gate_q.Fill(sktqz_.qiskz[i]);
		}
	}
	
	// debug print
	std::cout<<"1.3us window hit counts:\n";
	for(int i=0; i<hitcounts_in_1p3us.size(); ++i){
		std::cout<<"Subtrigger "<<i<<": "<<hitcounts_in_1p3us.at(i);
		ahist=(TH1D*)h_hits.GetHists()->At(i+2);
		std::cout<<" ("<<ahist->GetEntries()<<")\n";
	}
	
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
	
	fout->cd();
	c_subtriggers->Clear();
	c_subtriggers->cd();
	h_hits.Draw();
	std::string canv_name = "event_"+std::to_string(myTreeReader->GetEntryNumber());
	c_subtriggers->Write(canv_name.c_str());
	
	return true;
}

int PlotHitTimes::GetSubtriggerFlags(int subtrigtype, std::vector<std::bitset<32>>& subtrigger_flags, int offset){
	// scan for subtriggers, and for every hit, note which subtriggers its in the gate of
	
	// arbitrary limit of max num subtriggers. 32 should be more than sufficient
	// note this has to be equal to or smaller than the size of the bitset we use
	int MAX_SUBTRIGS=32;
	
	// we need the lun to call skcread
	int lun = m_data->GetLUN(myTreeReader->GetName());
	lun = -std::abs(lun);  // IMPORTANT: make it negative to ensure skcread does not advance TTree
	
	// trigger time of the primary trigger
	int it0sk = skheadqb_.it0sk;
	std::cout<<"primary trigger is at "<<it0sk<<std::endl;
	
	int ntrigsfound=0;
	std::vector<int> t0_sub(MAX_SUBTRIGS,-1);  // relative time of subtrigger to IT0SK
	
	// run subtrigger algorithm to search for subtriggers of this type
	get_sub_triggers_(&subtrigtype, &ntrigsfound, t0_sub.data(), &MAX_SUBTRIGS);
	
	Log(m_unique_name+" found "+toString(ntrigsfound)+" subtriggers of type "
	      +TriggerIDToName(subtrigtype),v_message,m_verbose);
	
	// process the subtriggers
	int nhiterrs=5;
	for(int i=0; i<ntrigsfound; ++i){
		
		// trigger time of the subtrigger
		int it0xsk = it0sk + t0_sub.at(i);
		std::cout<<"subtrigger "<<i<<" is at "<<it0xsk<<std::endl;
		
		// get in-gate times (for check?)
		float pre_t0 = 0;  //skruninf_.softtrg_pre_t0[0];  FIXME not sure what index to use
		float post_t0 = 0; //skruninf_.softtrg_post_t0[0]; FIXME not sure what index to use
		float subtr_start_t = it0xsk - pre_t0;
		float subtr_end_t = it0xsk + post_t0;
		
		// set IT0XSK to the position of the next subtrigger
		set_timing_gate_(&it0xsk);
		
		int n_in_gate_hits_this_subtrigger=0;
		
		// scan hits, check whether they're in the new 1.3us gate
		for(int k=0; k<sktqz_.nqiskz; ++k){
			std::unordered_map<std::string, int> hit_flags = GetHitFlagNames(sktqz_.ihtiflz[k]);
			
			// sanity check
			/* FIXME enable when we know what index to use for pre- and post- t0's
			float new_time = sktqz_.tiskz[i] - (it0xsk - it0sk)/COUNT_PER_NSEC;
			if(nhiterrs>0 && hit_flags.at("in 1.3us") != (new_time > subtr_start_t && new_time < subtr_end_t)){
				Log(m_unique_name+" hit "+toString(k)+" in subtrigger "+toString(j)+" has time "
				   +toString(new_time)+" relative to subtrigger window ["+toString(subtr_start_t)
				   +","+toString(subtr_end_t)+"], but in-gate flag "+toString(hit_flags.at("1.3us")),
				   v_error,m_verbose);
				--nhiterrs;  // limit to the first N occurrences so we don't swamp the user
			}
			*/
			
			if(hit_flags.at("in 1.3us")){
				subtrigger_flags.at(k).set(i+offset);
				++n_in_gate_hits_this_subtrigger;
			} else {
				subtrigger_flags.at(k).reset(i+offset);
			}
		}
		
		//std::cout<<"Found "<<n_in_gate_hits_this_subtrigger
		//         <<" hits with in-1.3us flag set for subtrigger "<<i<<std::endl;
		
	}
	
	return ntrigsfound;
}


bool PlotHitTimes::Finalise(){
	
	fout->Write("",TObject::kOverwrite);
	fout->Close();
	delete fout;
	delete c_subtriggers;
	return true;
}
