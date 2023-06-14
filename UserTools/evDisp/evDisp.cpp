#include "evDisp.h"

#include <utility>
#include <cmath>
#include <math.h>
#include <bitset>

#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TColor.h"
#include "TView.h"

#include "fortran_routines.h"
#include "Constants.h"

#include "TableReader.h"
#include "TableEntry.h"

evDisp::evDisp():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

namespace{
	const double caplimit = 1700;
	const double barrellimitY = 2000;
}

bool evDisp::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("treeReaderName",treeReaderName);
	m_variables.Get("plotVar",plotVar);
	m_variables.Get("dataSrc",dataSrc);
	m_variables.Get("inGateOnly",inGateOnly);
	m_variables.Get("evtSrc",evtSrc);
	m_variables.Get("plotStyle",plotStyle);
	
	if(dataSrc == 1 || (dataSrc == -1 && treeReaderName!="")){
		// if getting data from TTree, check the TreeReader
		 if(m_data->Trees.count(treeReaderName)==0){
			Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
			if(dataSrc==1) return false; // we may not need it...
		} else {
			myTreeReader = m_data->Trees.at(treeReaderName);
		}
	}
	
	// class to get PMT positions from cable number
	myConnectionTable = m_data->GetConnectionTable();
	
	gROOT->SetStyle("Plain");
	if(plotStyle==0){
		// polymarker versions
		Log(toolName+" using plotStyle 0: TGraphs",v_debug,verbosity);
		topCapHitMap = new TGraph2D();
		topCapHitMap->SetMarkerStyle(7);
		bottomCapHitMap = new TGraph2D();
		bottomCapHitMap->SetMarkerStyle(7);
		barrelHitMap = new TGraph2D();
		barrelHitMap->SetMarkerStyle(7);
	} else if(plotStyle==1){
		// histgrams for plotting
		Log(toolName+" using plotStyle 1: Histograms",v_debug,verbosity);
		topCapHeatMap = new TH2D("topCapHeatMap", "topCap", 50, -2000, 2000, 50, -2000, 2000);
		bottomCapHeatMap = new TH2D("bottomCapHeatMap", "bottomCap", 50, -2000, 2000, 50, -2000, 2000);
		barrelHeatMap = new TH2D("barrelSideHeatMap", "barrelSide", 150, -3.14, 3.14, 50, -2000, 2000);
		
	}
	
	if(plotStyle==0 || plotStyle==1){
		// Setting up the histograms for time and charge distributions, ranges might want
		// looking at in more detail.
		hitsVTimeHist = new TH1F("hitsVTimeHist", "# of hits vs time", 200, 950, 1400);
		hitsVTimeHist->SetFillColor(38);
		hitsVTimeHist->SetLineColor(38);
		hitsVTimeHist->GetXaxis()->SetTitle("Time (ns)");
		
		hitsVChargeHist = new TH1F("hitsVChargeHist", "Charge deposited per PMT", 100, 0, 10);
		hitsVChargeHist->SetFillColor(46);
		hitsVChargeHist->SetLineColor(46);
		hitsVChargeHist->GetXaxis()->SetTitle("Charge (pe)");
	}
	
	// canvas for plotting
  if(plotStyle==0 || plotStyle==1){
		displayCanvas = new TCanvas();  // if we name it, name must be unique - prevents duplicate Tools.
		((TRootCanvas*)displayCanvas->GetCanvasImp())->Resize(1024,700); // FIXME make unique names automatically
		
		/* divide up the canvas:
		 __________
		| ———  ——— |
		||TOP||BOT||
		| ———  ——— |
		| ———————— |
		|| BARREL ||
		| ———————— |
		| ———  ——— |
		||HST||HST||
		| ———  ——— |
		 ‾‾‾‾‾‾‾‾‾‾
		*/
		displayCanvas->Divide(1,3);
		displayPad = displayCanvas->cd(1);
		displayPad->Divide(2,1);
		displayCanvas->cd(1);
		displayPad->cd(1);
		gPad->SetFrameFillColor(1); // black fill
		displayPad->cd(2);
		gPad->SetFrameFillColor(1); // black fill
		displayCanvas->cd(2);
		gPad->SetFrameFillColor(1); // black fill
		botDisplayPad = displayCanvas->cd(3);
		botDisplayPad->Divide(2,1);
 	 gStyle->SetOptStat(0);   // disable stats box
	}
	
//	gStyle->SetPalette(57);  // kBird doesn't seem to exist in ROOT 5, do it ourselves
//	Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000 };

	// SuperScan-like colour pallette to replicate the colours used by the official SK plots.
	// p.e ranges that each stop represents are also the same as SuperScan plots that can be found
	// on the SK website.
	Double_t stops[17] = {0.0000, 0.0083, 0.0166, 0.0291, 0.0458, 0.0667, 0.0908, 0.1183, 0.1500, 0.1850, 0.2242, 0.2667, 0.3125, 0.3625, 0.3950, 0.4150, 1.0000};  // *see note below
	Double_t r[17]     = {0.4745, 0.3294, 0.1137, 0.3451, 0.3686, 0.3922, 0.5725, 0.6471, 0.8078, 0.9608, 0.9255, 0.8745, 0.8627, 0.8235, 0.8588, 0.9686, 1.0000};
	Double_t g[17]     = {0.2078, 0.1686, 0.1255, 0.4196, 0.5373, 0.5490, 0.6980, 0.7924, 0.8863, 0.9882, 0.9412, 0.7333, 0.7255, 0.5804, 0.4275, 0.1255, 1.0000};
	Double_t b[17]     = {0.8118, 0.8196, 0.8039, 0.8431, 0.4588, 0.1961, 0.2039, 0.2039, 0.2627, 0.3098, 0.2824, 0.2392, 0.2157, 0.1804, 0.1569, 0.1294, 1.0000};
	
	int Idx = TColor::CreateGradientColorTable(17, stops, r, g, b, 255);
	int nlevels=100;
	Int_t MyPalette[nlevels];
	for (int i=0;i<nlevels;i++) MyPalette[i] = Idx+i;
	gStyle->SetPalette(nlevels, MyPalette);
	// *note: ROOT 5 doesn't seem to stretch the colour range properly to span the full z axis range.
	// So while kBird ought to run from dark blue to yellow, the "default" (copied from ROOT 6)
	// colour palette will be almost entirely dark blue. We can mitigate this somewhat by plotting
	// on a logarithmic scale, and we can also do it by distorting the 'stops' as shown above.
	// Best results also seem to depend on the range of z values in the event...
	// Based on superscan, we clip the max charge range at 26.7, and combined with the above stops
	// we seem to get reasonable displays, if a little lacking in colour.
	
	if(plotStyle==0){
		// to get the axes ranges to work properly we need some dummy draws
		// note that for TGraph2D::Draw not to have issues it seems to need
		// at least points without degeneracies, and sufficient separation
		topCapHitMap->SetPoint(0,-10,10,26.7);
		topCapHitMap->SetPoint(1,1,1,1);
		topCapHitMap->SetPoint(2,2,2,2);
		topCapHitMap->SetMargin(0.03);
		topCapHitMap->GetHistogram()->SetBins(50,-caplimit,caplimit, 50,-caplimit, caplimit);
		bottomCapHitMap->SetPoint(0,-10,10,26.7);
		bottomCapHitMap->SetPoint(1,1,1,1);
		bottomCapHitMap->SetPoint(2,2,2,2);
		bottomCapHitMap->SetMargin(0.03);
		bottomCapHitMap->GetHistogram()->SetBins(50,-caplimit,caplimit, 50,-caplimit, caplimit);
		barrelHitMap->SetPoint(0,-10,10,26.7);
		barrelHitMap->SetPoint(1,1,1,1);
		barrelHitMap->SetPoint(2,2,2,2);
		barrelHitMap->SetMargin(0.03);
		barrelHitMap->GetHistogram()->SetBins(50,-3.2,3.2, 50,-barrellimitY, barrellimitY);
	}
	
	// if dataSrc is set to a negative value then call FindHits to look for a data source that is populated.
	// XXX check that the appropriate data structures are populated by the time we reach this
	if(dataSrc == -1){
		get_ok = FindHits();
		if(get_ok==0){
			Log(toolName+" ERROR! Couldn't find hits anywhere!", v_error, verbosity);
			return false;
		}
	}
	
	// print whether SECONDARY and LOWE branches are present and populated
	if(dataSrc == 1){
		BranchCheck();
	}
	
	return true;
}


bool evDisp::Execute(){
	
	// not sure where we want to get our hits:
	// sktqz_, rawtqinfo_, TQREAL, or somewhere else?
	
	/* rfm data files don't have TQREAL populated, even for events with a valid HEADER.
	   Instead they have a valid TQLIST branch, which is a TClonesArray of TQ class objects.
	   This gets read out by skrawread which calls TQRAWSK, which calls TQSKZ, which calls skroot_get_idtq
	   to read the TQLIST branch. The returned data gets passed into to TQSKZ, which puts it into the
	   common block array IQISKZ. After returning to TQRAWSK, this calls SKTQCONV_QB to do charge conversion
	   from TDC/ADC counts to ns/p.e., and puts the result into TBUF_RAW (or QBUF_RAW, there are basically
	   arrays of T, Q, and cable number here).
	   Finally, skroot_set_tree calls skroot_set_tqreal to transfer the converted hits
	   from QBUF_RAW etc. into the TQREAL branch.
	   The upshot is TQLIST is used to populate TQREAL after conversion from TDC/ADC counts to ns/p.e.
	   Once populated TQREAL has members `int nhits`, and a `std::vector<float> T, Q`
	   for hit times and charges and a `std::vector<int> cables` for cable numbers of the hits.
	   following from $SKOFL_ROOT/inc/sktq.h:
	   
	*     ---------
	*     Contents of /SKHEADQB/:
	*     ---------
	*     IT0SK   ; Original T0 of the event
	*     IT0XSK  ; T0 of the event for ITISKZ,ITASKZ,ITABSK,TISK,IHCAB,TASK,IHACAB,...
	*     NUMHWSK ; Number of HW triggers
	*     HWSK(I) ; TRG EVENT COUNTER list (I=1,NUMHWSK)
	*     NTRIGSK ; sub-trigger # (=(it0xsk-it0sk)/count_per_nsec/10)
	*              
	*     ---------
	*     Contents of /SKQ/:
	*     ---------
	*          QISK    ; Q of individual PMT's
	*          NQISK   ; # of hit PMT's
	*          QISMSK  ; sum of Q
	*     ---------
	*     Contents of /SKCHNL/:
	*     ---------
	*          IHCAB(I);  Hit cable number (I=1,NQISK)
	*     ---------
	*     Contents of /SKT/:
	*     ---------
	*          TISK    ; T of individual PMT's
	*     ---------
	*     Contents of /SKTQZ/:
	*     ---------
	*          NQISKZ    ; Number of ALL ID hits     (not just near sw trigger)
	*          TISKZ(I)  ; T (ns) for ALL ID hits    (I = 1, NQISKZ)
	*          QISKZ(I)  ; Q (pe) for ALL ID hits    (I = 1, NQISKZ)
	*          IHTIFLZ(I); Hit flags for ALL ID hits (I = 1, NQISKZ)
	*                        11-6   (# of TRG EVENT COUNTER - 1) * 64 (0-63)
	*                         5-4   charge range (0:Small, 1:Medium, 2:Large)
	*                         3-2   trig ID (0: Narrow,   1: Wide
	*                                        2: Pedestal, 3: Not used)
	*                         1bit  In gate (1=in gate, 0=not in gate)
	*                         0bit  In 1.3usec (1=in, 0=out)
	*          ICABIZ(I) ; Cable number for ALL ID hits (I = 1, NQISKZ)
	*/
	
	if(myTreeReader){
		// if we have both SHE and AFT available, load the appropriate common block data
		if(evtSrc==0) m_data->LoadSHE(treeReaderName);
		else          m_data->LoadAFT(treeReaderName);
	}
	if(dataSrc==1){
		GetData();  // get data from TreeReader if not using SK common blocks
	}
	
	
	long it0sk;
	long it0xsk;
	switch (dataSrc){
		case 0: {
			// sktqz_ common block
			totalPMTsActivated = sktqz_.nqiskz;
			triggerID = skhead_.idtgsk;
			break;
		}
		case 1: {
			// TQReal branch
			totalPMTsActivated = myTQReal->cables.size();
			it0sk = myHeader->t0;
			it0xsk = myTQReal->it0xsk;
			triggerID = myHeader->idtgsk;
			break;
		}
		case 2: {
			// skt_, skq_, skchnl_ commons
			totalPMTsActivated = skq_.nqisk;
			triggerID = skhead_.idtgsk;
			break;
		}
		default: {
			// unknown
			Log(toolName+" unknown dataSrc: "+std::to_string(dataSrc),v_error,verbosity);
			// TODO add support for rawtqinfo_ (see SimplifyTree Tool)
			totalPMTsActivated = 0;
			break;
		}
	}
	Log(toolName+" This event had "+std::to_string(totalPMTsActivated)+" PMTs hit",v_debug,verbosity);
	//if(totalPMTsActivated<10) return true;
	
	if(plotStyle==0){
		// reset the graphs? Things are acting strange... this seems to work
		topCapHitMap->Clear();
		bottomCapHitMap->Clear();
		barrelHitMap->Clear();
	}
	
	int top_cap_pmts_hit=0, bottom_cap_pmts_hit=0, barrel_pmts_hit=0;
	std::pair<float, float> top_cap_range{0, 1};
	std::pair<float, float> bot_cap_range{0, 1};
	std::pair<float, float> barrel_range{0, 1};
	float varMin=1E9, varMax=-1E9;
	
//	std::cout<<"********* SETTING POINTS ***********"<<std::endl;
	
	if(totalPMTsActivated == 0){
		// HEADER branch
		// extract the run number, subrun number, number of events and the trigger type
		// from the HEADER branch
		runNum = myHeader->nrunsk;
		subrunNum = myHeader->nsubsk;
		eventNum = myHeader->nevsk;
		
		readoutT0 = myHeader->t0;
		gateWidth = myHeader->gate_width;
		if(myLowe->bsenergy != 0){
			recVertex = myLowe->bsvertex;
			recEnergy = myLowe->bsenergy;
			recDirection = myLowe->bsdir;
		}
	}
	
	for (int pmtNumber = 0; pmtNumber < totalPMTsActivated; ++pmtNumber){
		switch (dataSrc){
			case 0: {
				// sktqz_ common block
				cableNumber = sktqz_.icabiz[pmtNumber];
				charge = sktqz_.qiskz[pmtNumber];             // FIXME which, if any of these 3, should use
				time = sktqz_.tiskz[pmtNumber];               // `cableNumber-1` instead of pmtNumber???
				in_gate = sktqz_.ihtiflz[pmtNumber] & 0x02;   // sktqaz_.ihtflz[pmtNumber]  for OD
				// skhead_ common block
				runNum = skhead_.nrunsk;
				subrunNum = skhead_.nsubsk;
				eventNum = skhead_.nevsk;
				readoutT0 = skheadqb_.it0sk;
				gateWidth = skheadqb_.gatewsk;
				recVertex = skroot_lowe_.bsvertex;
				break;
			}
			case 1: {
				// TQReal branch
				cableNumber = myTQReal->cables.at(pmtNumber);
				charge = myTQReal->Q.at(pmtNumber);
				time = myTQReal->T.at(pmtNumber);
				// upper 16 bits encode IHTIFLZ (tqrealsk.F::117)
				in_gate = cableNumber >> 16;
				// for PMT number extract lower 16 bits (tqrealsk.F::121)
				cableNumber = cableNumber & 0x0000FFFF;
				// convert hit time within readout window to hit time within subtrigger window
				// by subtracting the time from 40us buffer readout start to subtrigger start
				// then convert TDC ticks to nanoseconds
				// COUNT_PER_NSEC is a #defined constant in '$SKOFL_ROOT/inc/skheadC.h'
				// describing conversion from TDC ticks to nanoseconds
				time = time -(it0xsk-it0sk)/COUNT_PER_NSEC;
				// HEADER branch
				// extract the run number, subrun number, number of events and the trigger type
				// from the HEADER branch
				runNum = myHeader->nrunsk;
				subrunNum = myHeader->nsubsk;
				eventNum = myHeader->nevsk;
				readoutT0 = myHeader->t0;
				gateWidth = myHeader->gate_width;
				if(myLowe->bsenergy != 0){
					recVertex = myLowe->bsvertex;
					recEnergy = myLowe->bsenergy;
					recDirection = myLowe->bsdir;
				}
				break;
			}
			case 2: {
				//skt_, skq_, skchnl_ common blocks
				//note index for skq and skt has to be cableNumber - 1 to deal with the conversion from
				//fortran indexing to c indexing
				cableNumber = skchnl_.ihcab[pmtNumber];
				charge = skq_.qisk[cableNumber-1];
				time = skt_.tisk[cableNumber-1];
				in_gate = true;  // TODO where is ihtiflz in this case...?
				// not sure what to do about event info, only seems to be in skhead
				// for now will just use skhead
				runNum = skhead_.nrunsk;
				subrunNum = skhead_.nsubsk;
				eventNum = skhead_.nevsk;
				gateWidth = skheadqb_.gatewsk;
				readoutT0 = skheadqb_.it0sk;
				break;
			}
			default: {
				// unknown
				Log(toolName+" unknown dataSrc: "+std::to_string(dataSrc),v_error,verbosity);
				totalPMTsActivated = 0;
				break;
			}
		}
		
		//std::cout << "cable number is: " << cableNumber << std::endl;
		/*
		 # from $SKOFL_ROOT/const/connection.super.sk-4.dat
		 # -- for inner-PMTs (1-11146)
		 # -- for muon VETO(11151,11152,11153,11154)
		 # -- for calibration ID (11155-?, see skveto.h for details)
		 # -- for muon chamber(only hut3 and hut4)
		 # -- for trigger ID QB (15001-15240, see skhead.h for details)
		 # -- for anti-PMT(20001-21885)
		*/
		// MAXPM is a #defined constant representing max cable number of ID PMTs, starting from 1.
		// MAXPMA likewise represents the number of OD PMTs,
		// can't seem to see a #defined constant for OD tubeID starts... it's 20001 from above
		if(cableNumber==0 || cableNumber>MAXPM) continue;
		// TODO add OD drawing.
		
		// see above IHTIFLZ definition for in_gate test
		if(inGateOnly && (std::bitset<8*sizeof(int)>(in_gate).test(1)==0) ) continue;
		
		// get tube position
		myConnectionTable->GetTubePosition(cableNumber, tubePosition);
		tubeRadialCoordinate = sqrt(pow(tubePosition[0], 2.f) + pow(tubePosition[1],2.f));
		tubeAngularCoordinate = acos(tubePosition[0] / tubeRadialCoordinate);
		if(tubePosition[1] > 0) tubeAngularCoordinate = -tubeAngularCoordinate;
		// get location (barrel, top/bottom cap, ID/OD...)
		// enum Locations{ kIDTop, kIDWall, kIDBot, kODTop, kODWall, kODBot }; (from ConnectionTable.cc)
		int loc = myConnectionTable->GetLocation(tubePosition[0],tubePosition[1],tubePosition[2]);
		
		if(verbosity>10){
			std::cout<<"hit on PMT "<<cableNumber<<std::endl;
			std::cout<<"charge is "<<charge<<", time is "<<time
				 <<", PMT position is "<<tubePosition[0]<<", "<<tubePosition[1]
				 <<", "<<tubePosition[2]<<std::endl;
		}
		
		// sanity check limits
		if(abs(tubeAngularCoordinate)>M_PI){
			std::cerr<<"angle outside pi range: "<<tubeAngularCoordinate<<std::endl;
		} else if(abs(tubePosition[2])>barrellimitY){
			std::cerr<<"y position out of barrel range: "<<tubePosition[2]<<std::endl;
		} else if(abs(tubePosition[0])>caplimit){
			std::cerr<<"x position out of range: "<<tubePosition[0]<<std::endl;
		} else if(abs(tubePosition[1])>caplimit){
			std::cout<<"z position out of range: "<<tubePosition[1]<<std::endl;
		}
		
		// try to match colour scale of superscan.
		if(charge>26.7) charge=26.7;
//		charge /= 26.7;
		
		// for now we only plot one thing
		var = (plotVar) ? time : charge; //log(charge+1.0f);
		
//		std::cout<<"hit var is "<<var<<std::endl;
		if(var<varMin) varMin=var;
		if(var>varMax) varMax=var;
		
		// choose plot from location
		if (loc==0){
//			std::cout<<"setting top cap point {"<<tubePosition[0]
//					 <<", "<<tubePosition[1]<<", "<<var<<"}"<<std::endl;
			if(topCapHeatMap) topCapHeatMap->Fill(tubePosition[0], tubePosition[1], var);
			if(topCapHitMap)  topCapHitMap->SetPoint(top_cap_pmts_hit, tubePosition[0],
			                                         tubePosition[1], var);
			if(var<top_cap_range.first) top_cap_range.first = var;
			if(var>top_cap_range.second) top_cap_range.second = var;
			++top_cap_pmts_hit;
		} else if (loc==2){
//			std::cout<<"setting bottom cap point {"<<tubePosition[0]
//					 <<", "<<tubePosition[1]<<", "<<var<<"}"<<std::endl;
			if(bottomCapHeatMap) bottomCapHeatMap->Fill(tubePosition[0], tubePosition[1], var);
			if(bottomCapHitMap)  bottomCapHitMap->SetPoint(bottom_cap_pmts_hit,
			                                               tubePosition[0], tubePosition[1], var);
			if(var<bot_cap_range.first) bot_cap_range.first = var;
			if(var>bot_cap_range.second) bot_cap_range.second = var;
			++bottom_cap_pmts_hit;
		} else if(loc==1){
//			std::cout<<"setting barrel point {"<<tubeAngularCoordinate
//					 <<", "<<tubePosition[2]<<", "<<var<<"}"<<std::endl;
			if(barrelHeatMap) barrelHeatMap->Fill( tubeAngularCoordinate, tubePosition[2], var);
			if(barrelHitMap)  barrelHitMap->SetPoint(barrel_pmts_hit, tubeAngularCoordinate,
			                                         tubePosition[2], var);
			if(var<barrel_range.first) barrel_range.first = var;
			if(var>barrel_range.second) barrel_range.second = var;
			++barrel_pmts_hit;
		}
		
		// filling histograms of charge and time
		if(hitsVTimeHist) hitsVTimeHist->Fill(time);
		if(hitsVChargeHist) hitsVChargeHist->Fill(charge);
	}
	
//	std::cout<<"********* SETTING CONFIGS ***********"<<std::endl;
	if(plotStyle==0){
		// shrink internal arrays of TGraph2Ds so they don't plot old points
//		std::cout<<"Z val range is "<<varMin<<" to "<<varMax<<std::endl;
		
		// set a dummy point so we still get the histogram axes, rather than empty space
		if(top_cap_pmts_hit>3){
//			std::cout<<"setting N top cap hits to "<<top_cap_pmts_hit<<std::endl;
			topCapHitMap->Set(top_cap_pmts_hit);
		} else {
//			std::cout<<"setting dummy top cap hit"<<std::endl;
			topCapHitMap->SetPoint(0,-10,10,26.7);
			topCapHitMap->SetPoint(1,1,1,1);
			topCapHitMap->SetPoint(2,2,2,2);
			topCapHitMap->Set(3);
		}
		if(bottom_cap_pmts_hit>3){
//			std::cout<<"setting N bottm cap hits to "<<bottom_cap_pmts_hit<<std::endl;
			bottomCapHitMap->Set(bottom_cap_pmts_hit);
		} else {
//			std::cout<<"setting dummy bottom cap hit"<<std::endl;
			bottomCapHitMap->SetPoint(0,-10,10,26.7);
			bottomCapHitMap->SetPoint(1,1,1,1);
			bottomCapHitMap->SetPoint(2,2,2,2);
			bottomCapHitMap->Set(3);
		}
		if(barrel_pmts_hit>3){
//			std::cout<<"setting N barrel hits to "<<barrel_pmts_hit<<std::endl;
			barrelHitMap->Set(barrel_pmts_hit);
		} else {
//			std::cout<<"setting dummy barrel hit"<<std::endl;
			barrelHitMap->SetPoint(0,-10,10,26.7);
			barrelHitMap->SetPoint(1,1,1,1);
			barrelHitMap->SetPoint(2,2,2,2);
			barrelHitMap->Set(3);
		}
		
//		int npts= barrelHitMap->GetN();
//		double* xs=barrelHitMap->GetX();
//		double* ys=barrelHitMap->GetY();
//		double* zs=barrelHitMap->GetZ();
//		std::cout<<"barrel hits:\n";
//		for(int i=0; i<20; ++i){
//			std::cout<<"{"<<xs[i]<<", "<<ys[i]<<", "<<zs[i]<<"}\n";
//		}
//		std::cout<<std::endl;
		
//		std::cout<<"********* DRAWING ***********"<<std::endl;
		
//		std::cout<<"drawing top cap"<<std::endl;
		displayCanvas->cd(1);
		displayPad->cd(1);
		topCapHitMap->Draw("PCOL");
		displayCanvas->Modified(); displayCanvas->Update(); gSystem->ProcessEvents();
//		std::cout<<"gPad->GetView is "<<gPad->GetView()<<std::endl;
		if(gPad->GetView()) gPad->GetView()->TopView();
		
//		std::cout<<"drawing bottom cap"<<std::endl;
		displayCanvas->cd(1);
		displayPad->cd(2);
		bottomCapHitMap->Draw("PCOL");
		displayCanvas->Modified(); displayCanvas->Update(); gSystem->ProcessEvents();
//		std::cout<<"gPad->GetView is "<<gPad->GetView()<<std::endl;
		if(gPad->GetView()) gPad->GetView()->TopView();
		
//		std::cout<<"drawing barrel"<<std::endl;
		displayCanvas->cd(2);
		barrelHitMap->Draw("PCOL");
		displayCanvas->Modified(); displayCanvas->Update(); gSystem->ProcessEvents();
//		std::cout<<"gPad->GetView is "<<gPad->GetView()<<std::endl;
		if(gPad->GetView()) gPad->GetView()->TopView();
		
		displayCanvas->Modified(); displayCanvas->Update(); gSystem->ProcessEvents();
		
		if(histSet==0){
//			std::cout<<"drawing hits vs time"<<std::endl;
			displayCanvas->cd(3);
			botDisplayPad->cd(1);
			hitsVTimeHist->Draw();
			displayCanvas->Modified(); displayCanvas->Update(); gSystem->ProcessEvents();
			
//			std::cout<<"drawing hits vs charge"<<std::endl;
			botDisplayPad->cd(2);
			hitsVChargeHist->Draw();
			displayCanvas->Modified(); displayCanvas->Update(); gSystem->ProcessEvents();
		}
	} else if(plotStyle==1){
		// histograms
		displayCanvas->cd(1);
		displayPad->cd(1);
		topCapHeatMap->Draw("COL");
		displayPad->cd(2);
		bottomCapHeatMap->Draw("COL");
		displayCanvas->cd(2);
		barrelHeatMap->Draw("COL");
	}
	
	// Print the event number, run number, sub run number and gate width for each event
	std::cout << std::endl << "Event " << eventNum << " was from run " << runNum << ", subrun " << 
	subrunNum << "." << std::endl;
	
	// event time and date
	if(dataSrc == 1){
		tm rundate = {0};
		rundate.tm_year = myHeader->ndaysk[0];
		rundate.tm_mon = myHeader->ndaysk[1] - 1;
		rundate.tm_mday = myHeader->ndaysk[2];
		rundate.tm_hour = myHeader->ntimsk[0];
		rundate.tm_min = myHeader->ntimsk[1];
		rundate.tm_sec = myHeader->ntimsk[2];
		time_t runtime = mktime(&rundate);         // need to use mktime to derive day of week
		std::string timestring = ctime(&runtime);  // format timestamp into a string
		timestring.pop_back();                     // drop trailing newline
		std::cout << "The event occurred at " << timestring << std::endl;
	}
	
	// print the gate width
	std::cout << "Gate width/recorded readout window: " << gateWidth << " ns" << std::endl;
	
	// print out the value for T0SK. Units need to be added.
	std::cout << "Readout IT0SK: " << readoutT0/(COUNT_PER_NSEC * 1000000) << " ms" << std::endl;
	
	// check for dataSrc being the TQReal branch before printing reconstructed info. There does seem to be a
	// common block for reconstructed info (skroot_lowe) but its not populated TODO < under what circumstances?
	if(dataSrc==1){
		// if the LOWE branch has been found to be present and populated then printout the 4-vector reconstructed
		// vertex
		if(recEnergy){
			// vertex and energy are reconstructed independently where vertex recon can pass but energy
			// recon can fail so have to check separately
			// if reconstruction does not go ahead a default value of 9999 passes for each variable
			// if one value in vertex is still set at 9999 all four are generally 9999 so only check one
			// value
			if(recVertex[0] == 9999){
				std::cout << "The vertex for this event could not be reconstructed." << std::endl;
			} else {
				// if vertex passes then go on to print out
				std::cout << "Reconstructed (x, y, z, t) vertex: (";
				// loop over the bsvertex with limit on the size of the vertex
				for(int i = 0; i < recVertex.size(); i++){
					std::cout << recVertex[i];
					// print commas in between vertex elements
					if(i < recVertex.size() - 1){
						std::cout << ", ";
					}
				}
				// close vertex parentheses
				std::cout << ")" << std::endl;
			}
			// similarly as for the vertex printout check one element for non-reconstruction
			if(recDirection[0] == 9999){
				std::cout << "The direction for this event could not be reconstructed." << std::endl;
			} else{
				// printout direction if the first element is not 9999
				std::cout << "Reconstructed direction: (";
				for(int i = 0; i < recDirection.size(); i++){
					std::cout << recDirection[i];
					if(i < recDirection.size() - 1){
						std::cout << ", ";
					}
				}
				std::cout << ")" << std::endl;
			}
			// check energy for a failed/not tried reconstruction and if it isn't 9999 then print
			if(recEnergy != 9999){
				std::cout << "Reconstructed energy: " << recEnergy << " MeV" << std::endl;
			} else std::cout << "The energy for this event could not be reconstructed." << std::endl;
		}
	}
	
	// initlialise string Trigs to store trigger strings in
	std::string Trigs;
	// loop over the 31 trigger IDs and check if the trigger bit for each ID, i, is 1
	// pull trigger strings from Trigger_ID_To_Trigger using Constants::TriggerIDToTrigger 
	// and append them to trigs
	for(int i=0; i<=31; i++){
		if(triggerID.test(i)) Trigs += " " + TriggerIDToTrigger(i) + ",";
	}
	// check that Trigs is not empty, or pop_back makes evDisp fail the execute
	// if Trigs is not empty then print the trigger flags that were set for the event
	// if not print a warning for the event that no trigger falgs were found
	if(!Trigs.empty()){
		Trigs.pop_back();
		std::cout << "Trigger flags set:" << Trigs << "." << std::endl;
	} else{
		std::cout << "WARNING - There were no trigger flags set/found for this event." << std::endl;
	}
	
	if(plotStyle==0 || plotStyle==1) gPad->WaitPrimitive();
	
	// reset the charge and time histograms, this is required or else the bins never empty and the histograms
	// are filled with information from every event cummulatively
	if(hitsVChargeHist) hitsVChargeHist->Reset("ICESM");
	if (hitsVTimeHist) hitsVTimeHist->Reset("ICESM");
	
	// clear variables
	runNum = 0;
	subrunNum = 0;
	eventNum = 0;
	triggerID = std::bitset<32>{};
	readoutT0 = 0;
	gateWidth = 0;
	if(dataSrc==1){
		if(myLowe->bsenergy != 0){
			recVertex = basic_array<float>();
			recEnergy = 0.0;
			recDirection = basic_array<float>();
		}
	}
	
	
	//displayCanvas->SaveAs("HeatMap.png");
	
	return true;
}


bool evDisp::Finalise(){
	
	Log(toolName+" performing cleanup",v_debug,verbosity);
	// graphs
	if(topCapHitMap) delete topCapHitMap;
	if(bottomCapHitMap) delete bottomCapHitMap;
	if(barrelHitMap) delete barrelHitMap;
	// histograms
	if(topCapHeatMap) delete topCapHeatMap;
	if(bottomCapHeatMap) delete bottomCapHeatMap;
	if(barrelHeatMap) delete barrelHeatMap;
	if(hitsVChargeHist) delete hitsVChargeHist;
	if (hitsVTimeHist) delete hitsVTimeHist;
	// canvas
	if(displayCanvas) delete displayCanvas;
	
	return true;
}

bool evDisp::GetData(){
	myTreeReader->Get("TQREAL", myTQReal);
	myTreeReader->Get("HEADER", myHeader);
	return true;
}

int evDisp::FindHits(){
	Log(toolName+": Looking for data sources... ",v_debug,verbosity);
	// first check for TQReal
	logmsg = toolName+": TQREAL branch is ";
	if(myTreeReader->GetTree()->FindBranch("TQREAL")){
		// TQREAL branch exists, now check if it is populated with hits
		myTreeReader->Get("TQREAL", myTQReal);
		if(myTQReal->nhits > 0){
			logmsg += "present and populated, setting dataSrc = 1";
			Log(logmsg,v_debug,verbosity);
			dataSrc = 1;
			return 1;
		} else logmsg += "present but empty";
	} else logmsg += " absent";
	Log(logmsg,v_debug,verbosity);
	
	// no TQReal hits so check the sktqz_ common block
	if(sktqz_.nqiskz > 0){
		Log(toolName + ": sktqz_ block is populated, setting dataSrc = 0",v_debug,verbosity);
		dataSrc = 0;
		return 1;
	}
	Log(toolName+": sktqz_ block is empty",v_debug,verbosity);
	
	if(skt_.tisk && skq_.qisk){
		Log(toolName+" skt_ and skq_ blocks are populated, setting dataSrc = 2",v_debug,verbosity);
		dataSrc = 2;
		return 1;
	}
	Log(toolName+": skt_ and skq_ blocks are empty",v_debug,verbosity);
	// if no hits found anywhere return 0
	return 0;
}

bool evDisp::BranchCheck(){
	// check the SECONDARY branch exists
	logmsg = toolName+": SECONDARY branch is";
	if(myTreeReader->GetTree()->FindBranch("SECONDARY")){
		logmsg += " present!";
		// it exists! now check if it has any secondary particle information in it
		// first atmpd c-style arrays
		myTreeReader->Get("SECONDARY", mySecondary);
		if(mySecondary->nscndprt > 0){
			logmsg += " atmpd arrays are populated";
		} else {
			logmsg += " atmpd arrays do not seem populated";
		}
		// also check the SECONDARY vectors
		if(mySecondary->vertex_time.size() > 0){
			logmsg += ", secondary vectors are populated";
		} else {
			logmsg += ", secondary vectors do not seem populated.";
		}
	} else {
		logmsg += " absent";
	}
	Log(logmsg,v_debug,verbosity);
	
	
	//check the LOWE branch exists
	logmsg = toolName+": LOWE branch is";
	if(myTreeReader->GetTree()->FindBranch("LOWE")){
		logmsg += " present";
		// it exists! now check if it is populated. Do this by checking if bsenergy has a value.
		// FIXME perhaps there is a better way...
		myTreeReader->Get("LOWE", myLowe);
		if(myLowe->bsenergy != 0){
			// lowe branch is populated!
			logmsg+= " and appears populated!";
		// lowe branch is unpopulated
		} else {
			logmsg += " but appears unpopulated.";
		}
	} else {
		logmsg += " absent";
	}
	Log(logmsg,v_debug,verbosity);
	
	return true;
}
