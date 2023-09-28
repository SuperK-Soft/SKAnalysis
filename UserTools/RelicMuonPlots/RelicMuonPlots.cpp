#include "RelicMuonPlots.h"
#include "MTreeReader.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include <cmath>

RelicMuonPlots::RelicMuonPlots():Tool(){}


bool RelicMuonPlots::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("outputFile", outputFile);
	m_variables.Get("relicReaderName",relicReaderName);
	
	// get intput treereaders
	if(m_data->Trees.count(relicReaderName)==0){
		Log(m_unique_name+" Error! No TreeReader '"+relicReaderName+"' in DataModel!",v_error,m_verbose);
		return false;
	}
	relicReader = m_data->Trees.at(relicReaderName);
	
	// since we'll need to loop over muon entries for each relic,
	// we won't be reading one entry per loop, so a TreeReader Tool isn't really suitable.
	// instead just use the underlying MTreeReader class.
	std::string inputFile = relicReader->GetFile()->GetName();
	muReader.Load(inputFile, "mu");
	
	// make output file
	MakeHists(0);
	
	return true;
}


bool RelicMuonPlots::Execute(){
	
	// get next relic candidate
	get_ok = GetRelicEvt();
	if(!get_ok){
		Log(m_unique_name+" Error getting relic entry "+toString(relicReader->GetEntryNumber()),
		    v_error,m_verbose);
	}
	
	// make sure it was reconstructed
	if(relicLowe->bsenergy != 0){             // reconstruction not done
		Log(m_unique_name+" relic bsenergy == 0, was lowe reconstruction done?",v_error,m_verbose);
		return false;
	} else if(relicLowe->bsenergy == 9999){   // reconstruction failed
		Log(m_unique_name+" relic bsenergy == 9999, reconstruction failed, skipping",v_message,m_verbose);
		return true;
	}
	
	
	// loop over muons matched to this relic
	std::cout<<"looping over "<<muEvNums->size()<<" muons for this relic"<<std::endl;
	for(int& amuevnum : *muEvNums){
		std::cout<<"next muon entry: "<<amuevnum<<std::endl;
		get_ok = GetMuonEvt(amuevnum);
		if(!get_ok){
			Log(m_unique_name+" Error getting muon entry "+toString(amuevnum),v_error,m_verbose);
			continue;
		}
		
		MakePairVariables();
		
		// fill distributions
		MakeHists(1);
		
	}
	
	return true;
}


bool RelicMuonPlots::Finalise(){
	
	MakeHists(2);
	
	return true;
}


bool RelicMuonPlots::GetRelicEvt(){
	
	// get next relic
	get_ok  = relicReader->Get("HEADER", relicHeader);
	get_ok &= relicReader->Get("LOWE", relicLowe);
	get_ok &= relicReader->Get("MatchedEvNums", muEvNums);
	
	return get_ok;
}

bool RelicMuonPlots::GetMuonEvt(int entrynum){
	
	// load next entry data from TTree
	int bytesread = muReader.GetEntry(entrynum);
	
	// stop loop if we ran off the end of the tree
	if(bytesread<1&&bytesread>-3){
		Log(m_unique_name+" Ran off end of muon file, stopping loop",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	// stop loop if we had an error of some kind
	else if(bytesread<0){
		std::string err;
		if(bytesread==-1) err=" IO error";
		if(bytesread==-10) err=" AutoClear error";
		if(bytesread <-2) err=" Unknown error";
		m_data->vars.Set("StopLoop",1);
		Log(m_unique_name+err+" loading muon entry "+toString(entrynum),v_error,m_verbose);
		return false;
	}
	
	// otherwise read the entry ok; get branches
	get_ok  = muReader.Get("HEADER", muHeader);
	get_ok &= muReader.Get("MU", muMuInfo);
	
	return get_ok;
}

bool RelicMuonPlots::MakePairVariables(){
	
	std::cout<<"making pair variables"<<std::endl;
	
	// plot time difference
	int64_t relic_ticks = (relicHeader->counter_32*32768) + (relicHeader->t0 & 32767);
	float relic_time = relic_ticks / COUNT_PER_NSEC;
	int64_t mu_ticks = (muHeader->counter_32*32768) + (muHeader->t0 & 32767);
	float mu_time = mu_ticks / COUNT_PER_NSEC;
	dt = relic_time - mu_time;
	std::cout<<"dt = "<<(dt/1E9)<<" seconds"<<std::endl;
	// XXX compare to MatchedTimeDiff branch to ensure they're the same.
	
	// find position of max dedx
	// defined as bin where a sliding window of 4.5m (9x 50cm bins) has maximum sum
	double max_edep = 0;
	int max_edep_bin=0;
	for(int i=0;i<111;i++){   // from mu_info.C - XXX why 111? muboy_dedx array is 200 bins in length...?
		// calculate sum of dedx in 9 bins from bin i to i+9
		// (use 9 as an odd number so we have one central bin)
		double e_dep_in_window = 0 ;
		for(int j=0;j<9;j++){
			e_dep_in_window = e_dep_in_window + muMuInfo->muboy_dedx[i+j]; // will be qpeak
		}
		if(e_dep_in_window > max_edep){
			max_edep_bin = i+4;  // +4 to get centre bin
			max_edep = e_dep_in_window;
		}
	}
	double max_edep_pos = 50.*max_edep_bin;
	std::cout<<"muon track max dE/dx: "<<max_edep<<" at "<<max_edep_pos<<" cm along the track"<<std::endl;
	// FIXME getting 0 from this?
	
	// calculate footpoint and transverse distance
	// inputs:
	float (*mudir)[3]=nullptr;
	float (*muentry)[3]=nullptr;
	// choose appropriate source of muon info: BFF if it happened, muboy otherwise
	if(muMuInfo->mubff_goodness > 0.3){ // should be a sufficient check that we did BFF and it was good?
		muentry = &muMuInfo->mubff_entpos;
		mudir = &muMuInfo->mubff_dir;
	} else {
		muentry = &muMuInfo->muentpoint;
		mudir = &muMuInfo->muboy_dir;
	}
	
	// find transverse distance from relic to muon
	float* vertex_relic = relicLowe->bsvertex;
	float appr;  // output: XXX is this "foot point"? distance along track where dlt is defined relative to?
	std::cout<<"calling getdl_"<<std::endl;
	std::cout<<"relic vertex: ("<<vertex_relic[0]<<", "<<vertex_relic[1]<<", "<<vertex_relic[2]<<")\n"
	         <<"muon entry point: ("<<*muentry[0]<<", "<<*muentry[1]<<", "<<*muentry[2]<<")\n"
	         <<"muon entry dir: ("<<*mudir[0]<<", "<<*mudir[1]<<", "<<*mudir[2]<<")"<<std::endl;
	getdl_(*mudir, &vertex_relic[0], &vertex_relic[1], &vertex_relic[2], *muentry, &dlt, &appr);
	std::cout<<"transverse distance: "<<dlt<<", with foot point at "<<appr<<" cm along muon track"<<std::endl;
	// FIXME getting 0 from this
	// XXX warning: see Shinoki-san's talk from lowe meeting ~06/05/21:
	// lt distriubtion for multiple muons with low goodness (<~0.4) is quite large,
	// but cutting on muboy goodness can bias muon energy and neutron multiplicity (by affecting search window)
	// he evaluates systematics of this cut; we should look into this.
	// see his slides also for example distributions
	
	// calculate longitudinal distance from point of peak energy deposition
	dll = max_edep_pos - appr;
	std::cout<<"longitudinal distance: "<<dll<<std::endl;
	// FIXME getting 0 from this
	
	// charge per cm of a minimum-ionizing-particle (MIP)
	float peprcm = 26.78;
	// 1.992 Mev cm^2/g dE/dx min   (from: Muon Stopping Power and Range Tables, D.E.Groom et al, LBNL-44742,
	// Atomic Data and Nuclear Data Tables, Vol. 76, No. 2, July 2001 - ~/Downloads/adndt.pdf)
	// at 13.5C density of water 0.999315 g/cm^3
	// (from: https://www.internetchemistry.com/chemical-data/water-density-table.php )
	// so 1.992*0.999315 = 1.99063548 MeV/cm (makes sense with commonly known value of muon loss ~2MeV/cm)
	
	return true;
}

bool RelicMuonPlots::MakeHists(int step){
	
	get_ok = true;
	
	// step 0: Initialise
	// ==================
	if(step==0){
		hb.MakeFile(outputFile);
		hb.SaveHists(false);
		
	// step 1: Fill
	// ============
	} else if(step==1){
		std::cout<<"Filling tree"<<std::endl;
		// fill tree
		hb.Fill("type", int(std::signbit(dt))); // true if negative (muon after relic)
		hb.Fill("dt", std::abs(dt));       // time diff
		hb.Fill("dll", dll);               // longitudinal distance from relic to point of muon max dedx
		hb.Fill("dlt", dlt);               // transverse distance from relic to point of muon max dedx
		
	// step 2: Finalise
	// ================
	} else {
		// names should be branch name + '_spall' or '_relic'
		std::vector<TH1*> spall_dists;
		std::vector<TH1*> random_dists;
		
		// make/get histograms
		std::cout<<"making dt hists"<<std::endl;
		spall_dists.push_back(hb.GetHist("dt","type==0"));
		random_dists.push_back(hb.GetHist("dt","type==1"));
		
		std::cout<<"making dll hists"<<std::endl;
		spall_dists.push_back(hb.GetHist("dll","type==0"));
		random_dists.push_back(hb.GetHist("dll","type==1"));
		
		std::cout<<"making dlt hists"<<std::endl;
		spall_dists.push_back(hb.GetHist("dlt","type==0"));
		random_dists.push_back(hb.GetHist("dlt","type==1"));
		
		// draw and save distributions
		std::cout<<"writing file"<<std::endl;
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
