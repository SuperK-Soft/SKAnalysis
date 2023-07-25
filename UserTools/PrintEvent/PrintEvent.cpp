#include "PrintEvent.h"
#include "MTreeReader.h"
#include "Constants.h"
#include "SuperWrapper.h"
#include <bitset>

PrintEvent::PrintEvent():Tool(){
	// get the name of the tool from its class name
	m_unique_name=type_name<decltype(this)>(); m_unique_name.pop_back();
}


bool PrintEvent::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	
	// get the reader for inputs
	std::string treeReaderName;
	m_variables.Get("treeReaderName",treeReaderName);
	if(m_data->Trees.count(treeReaderName)==0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	return true;
}




bool PrintEvent::Execute(){
	
	
	int file_format = skheadf_.sk_file_format;  // 0=zbs, 1=root
	std::cout<<"\n\nReading "<<(file_format ? "ROOT" : "ZBS")<<" file "
	         << myTreeReader->GetFile()->GetName()<<std::endl;
	
	PrintHeaderInfo();
	
	PrintRunInfo();
	
	PrintTriggerInfo(); // from commons... not sure this is all available from a single class.
	
	PrintHWTriggerInfo();
	
	PrintSWTriggerInfo();
	
	if(myTreeReader->GetMCFlag()) PrintMCInfo();
	
	PrintBadChannels();
	
	PrintDarkInfo();
	
	// from Reformatter.h, apparently there is a branch:
	// TClonesArray*      TRGList;   // ??? where is this? how relates to SWTRGLIST and HWTRGLIST?
	
	PrintPrevT0(); // what is this?
	
	PrintTQRealInfo();
	
	PrintHits();   // compare various sources: TQReal, commons, TODO TQLIST, TQRAWSK?
	
	//TClonesArray*       TQList;   // TODO
	//TClonesArray*     ODTQList;   // TODO
	
	std::cout<<"\n ==== COMMONS ====\n"<<std::endl;
	
	std::cout<<"SKHEADQB:\n"
	         <<"\tnevhwsk: "<<skheadqb_.nevhwsk<<"\n"    // "TRG EVENT COUNTER (where T0 have)" ..... VERY large?
	         <<"\tnevswsk: "<<skheadqb_.nevswsk<<"\n"    // software trigger id - small, != nevsk; is it a triggerID?
	         <<"\tntrigsk: "<<skheadqb_.ntrigsk<<"\n"    // sub-trigger # (=(it0xsk-it0sk)/count_per_nsec/10).
	         <<"\t\t(check: should == "<<(((skheadqb_.it0xsk-skheadqb_.it0sk)/COUNT_PER_NSEC)/10)<<")\n"
	         <<"\tnumhwsk: "<<skheadqb_.numhwsk<<"\n";   // == SoftwareTrigger::nhwtrgs == HWTRGLIST->GetEntries()?
	for(int i=0; i<skheadqb_.numhwsk; ++i){
		std::cout<<"\t\thwsk["<<i<<"]: "<<skheadqb_.hwsk[i]<<"\n"; // TRG EVENT COUNTER
	}
	
	std::cout<<"SKHEADC:\n"
	         <<"\tsk_condition: "<<skheadc_.sk_condition<<"\n"   // see sk_condition.dat
	         <<"\tsk_gd_concentration: "<<skheadc_.sk_gd_concentration<<"\n";  // Gd2(SO4)3 8H2O by weight
	
	std::cout<<"SKWATERLEN:\n"
	         <<"\tskwaterlen: "<<skwaterlen_.skwaterlen<<"\n"      // attenuation length from water job
	         <<"\tskwatergain: "<<skwaterlen_.skwatergain<<"\n";  // relative global PMT gain from water job
	
	/*
	std::cout<<"odmaskflag:\n"  // what are these? masked (bad) OD PMTs?
	         <<"\tod_mask_ver: "<<odmaskflag_.od_mask_ver<<"\n";
	         <<"\tod_mask_nhits: "<<odmaskflag_.od_mask_nhits<<"\n"
	for(int i=0; i<skqa_.odmaskflag_.od_mask_nhits; ++i){
	std::cout<<"OD PMT "<<i<<"\n:"
	         <<"\tod_mask_cable_askz: "<<odmaskflag_.od_mask_cable_askz[i]<<"\n"
	         <<"\tod_mask_nhits: "<<odmaskflag_.od_mask_nhits<<"\n"
	         <<"\tod_mask_nhits: "<<odmaskflag_.od_mask_nhits<<"\n";
	}
	*/
	
	return true;
}

bool PrintEvent::PrintTQRealInfo(bool verbose){
	// ID hit info
	TQReal* myTQReal=nullptr;
	get_ok = myTreeReader->Get("TQREAL", myTQReal);
	if(!get_ok || myTQReal==nullptr){
		Log(m_unique_name+" failed to get TQREAL from Tree",v_error,verbosity);
		return false;
	}
	// OD hit info
	TQReal* myTQAReal=nullptr;
	get_ok = myTreeReader->Get("TQAREAL", myTQAReal);
	if(!get_ok || myTQAReal==nullptr){
		Log(m_unique_name+" failed to get TQAREAL from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting TQREAL & TQAREAL Branches\n";
	
	/* members of TQReal, from tqrealroot.h
	Int_t nhits;
	Float_t pc2pe;
	Int_t tqreal_version;
	Int_t qbconst_version;
	Int_t tqmap_version;
	Int_t pgain_version;
	Int_t it0xsk;
	std::vector<int> cables;
	std::vector<float> T;
	std::vector<float> Q;
	*/
	
	std::cout<<"TQReal (ID)\n"
	         <<"\tversion: "<<myTQReal->tqreal_version<<"\n"
	         <<"\tQB const version: "<<myTQReal->qbconst_version<<"\n"
	         <<"\tTQ map version: "<<myTQReal->tqmap_version<<"\n"
	         <<"\tpgain version: "<<myTQReal->pgain_version<<"\n"   // what is pgain? PMT gain map?
	         <<"\tpC->p.e. conversion factor: "<<myTQReal->pc2pe<<"\n" // possibly obsolete? sktqC.h
	         <<"\tit0xsk [TDC count]: "<<myTQReal->it0xsk   // skheadqb_.it0xsk
	         <<" = " << (myTQReal->it0xsk)/(COUNT_PER_NSEC * 1000) << " us\n"
	         <<"\tNum hits: "<<myTQReal->nhits<<"\n"
	         << std::endl;
	
	// XXX i would expect most of these to be the same, or perhaps not relevant, for OD?
	std::cout<<"TQAReal (OD):\n"
	         <<"\tversion: "<<myTQAReal->tqreal_version<<"\n"
	         <<"\tQB const version: "<<myTQAReal->qbconst_version<<"\n"
	         <<"\tTQ map version: "<<myTQAReal->tqmap_version<<"\n"
	         <<"\tpgain version: "<<myTQAReal->pgain_version<<"\n"
	         <<"\tpC->p.e. conversion factor: "<<myTQAReal->pc2pe<<"\n"
	         <<"\tit0xsk [TDC count]: "<<myTQAReal->it0xsk
	         <<" = " << (myTQAReal->it0xsk)/(COUNT_PER_NSEC * 1000) << " us\n"
	         <<"\tNum hits: "<<myTQAReal->nhits<<"\n"
	         << std::endl;
	
	// The trigger window is defined by the time range [IT0XSK-pre_t0]->[IT0XSK+post_t0];
	// we expect a cluster of hits meeting the trigger condition to be at [IT0XSK-t0_offset]
	
	//    +-----------+--------+---------------------------------+
	//    |           |        |                                 |
	//  gate start  T(signal)  T0(IT0XSK)                       gate end
	//  =T0-pre_t0             =T(Signal)+t0_offset             =T0+post_t0
	
	// (from lowe school 01-skdetector.pdf:30)
	// The SKIV+ DAQ reads out a continuous stream of 17us blocks, multiple of which are then
	// merged and scanned for triggers online. When a trigger is found, ALL hits
	// from the relevant data chunks are written out - including some not in this trigger gate
	
	//                           trigger gate                            
	//                 /‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\                 
	//    +-----------+----------+---------+------------+---------------+
	//    ‖           |          |         ‖            |               ‖
	// block 1   gate start   trigger   block 2      gate end         block2
	//  start    (T0-pre_t0)    T0       start     (T0+post_t0)        end
	//     \____________________________________________________________/
	//                           readout gate                            
	
	// the 'in_gate' flag denotes whether the hits are in the trigger gate;
	// out of gate hits should be ignored (though they may be in-gate for a sub-trigger)
	// since the position of the trigger is random wrt to the block start/end,
	// a trigger may be wholly contained in one, two, or more blocks,
	// so the first and last hit times are an aribitrary distance from the trigger time.
	
	// the '1.3us' flag is a similar throwback to SKI-III triggers, which all had 1.3us
	// readout windows (SK-IV+ uses trigger window lengths that depend on the trigger type).
	// To emulate SKI-III style triggers, select only hits with the '1.3us' flag set.
	// (the `delete_outside_hits_` routine removes all hits with `1.3us` flag == 0)
	
	/*   Cable Numbers below from $SKOFL_ROOT/const/connection.super.sk-4.dat
	# -- for ID PMTs (1-11146) last ID PMT # (11146) #defined as MAXPM in tqrealroot.h
	# -- for muon VETO(11151,11152,11153,11154) 
	# -- for calibration ID (11155-11186, see skveto.h for details)
	# -- for trigger ID QB (15001-15024, see skhead.h IDTGSK cable # for mapping)
	# -- for OD PMTs (20001-21885)
	#       OD cable num offset (20000) is #defined as QB_OD_OFFSET in sktqC.h
	#       last OD PMT # (1885; excludes offset) #defined as MAXPMA in tqrealroot.h
	# -- for muon chamber(only hut3 and hut4)  .... no cable numbers? can't find any info about this
	*/
	
	// count the number of hits from various sources, and in or out of various gates
	std::map<std::string, int> hit_counts{ {"ID",0},
	                                       {"Veto",0},  // scintillation counters over cable holes in tank top
	                                       {"Calibration",0},
	                                       {"Trigger",0},
	                                       {"OD",0},
	                                       {"Unknown",0},
	                                       {"ID_in_gate",0},
	                                       {"OD_in_gate",0},
	                                       {"ID_in_13",0},
	                                       {"OD_in_13",0} };
	
	for (int i = 0; i < myTQReal->nhits; ++i){
		int cableNumber = myTQReal->cables.at(i);
		float charge = myTQReal->Q.at(i);
		float htime = myTQReal->T.at(i);
		
		// TQReal::cableNumber is a combination of hit flags and cable number;
		// upper 16 bits are flags, lower 16 bits are cable number
		int hitflags = cableNumber >> 16;       // upper 16 bits encode IHTIFLZ (tqrealsk.F::117)
		cableNumber = cableNumber & 0x0000FFFF; // lower 16 bits encode PMTNum  (tqrealsk.F::121)
		// bonsai demo uses 0x3fff = lower 14 bits for ihtiflz mask... only 11 seem significant anyway.
		
		if(verbose){
			// print hits if verbose
			std::string hitflagstring;
			GetHitFlagNames(hitflags, &hitflagstring);
			std::cout<<"\tHit "<<i<<", PMT: "<<cableNumber<<", Q: "<<charge<<" [units?]" // XXX
			         <<", T: "<<htime<<" [ns], flags: "<<hitflagstring<<"\n";
		}
		
		std::bitset<32> flags{hitflags};
		bool in_13us = flags.test(0);
		bool in_gate = flags.test(1);
		
		if(cableNumber>0 && cableNumber<MAXPM){
			++hit_counts.at("ID");
			if(in_13us) ++hit_counts.at("ID_in_13");
			if(in_gate) ++hit_counts.at("ID_in_gate");
		} else if(cableNumber>11150 && cableNumber<11155){
			++hit_counts.at("Veto");
		} else if(cableNumber>11154 && cableNumber<11186){
			++hit_counts.at("Calibration");
		} else if(cableNumber>15000 && cableNumber<15024){
			++hit_counts.at("Trigger");
		} else if(cableNumber>QB_OD_OFFSET && cableNumber<(QB_OD_OFFSET+MAXPMA)){
			++hit_counts.at("OD");
			if(in_13us) ++hit_counts.at("OD_in_13");
			if(in_gate) ++hit_counts.at("OD_in_gate");
		} else {
			++hit_counts.at("Unknown");
		}
	}
	
	std::cout<<"TQReal Hit counts:\n"
	         <<"\tID Hits: "<<hit_counts.at("ID")<<" of which "<<hit_counts.at("ID_in_gate")
	                        <<" were in-gate and "<<hit_counts.at("ID_in_13")<<" were in a 1.3us window\n"
	         <<"\tOD Hits: "<<hit_counts.at("OD")<<" of which "<<hit_counts.at("OD_in_gate")
	                        <<" were in-gate and "<<hit_counts.at("OD_in_13")<<" were in a 1.3us window\n"
	         <<"\tVeto Hits: "<<hit_counts.at("Veto")<<"\n"
	         <<"\tCalibration Hits: "<<hit_counts.at("Veto")<<"\n"
	         <<"\tTrigger Hits: "<<hit_counts.at("Veto")<<"\n"
	         <<"\tUnknown Hits: "<<hit_counts.at("Veto")<<std::endl;
	
	// =========================================================
	
	// same scan on TQAReal
	for(auto&& acount : hit_counts) acount.second = 0;
	for (int i = 0; i < myTQAReal->nhits; ++i){
		int cableNumber = myTQAReal->cables.at(i);
		float charge = myTQAReal->Q.at(i);
		float htime = myTQAReal->T.at(i);
		int hitflags = cableNumber >> 16;
		cableNumber = cableNumber & 0x0000FFFF;
		
		if(verbose){
			std::string hitflagstring;
			GetHitFlagNames(hitflags, &hitflagstring);
			std::cout<<"\tHit "<<i<<", PMT: "<<cableNumber<<", Q: "<<charge<<" [units?]" // XXX
			                   <<", T: "<<htime<<" [ns], flags: "<<hitflagstring<<"\n";
		}
		
		std::bitset<32> flags{hitflags};
		bool in_13us = flags.test(0);
		bool in_gate = flags.test(1);
		
		if(cableNumber>0 && cableNumber<MAXPM){
			++hit_counts.at("ID");
			if(in_13us) ++hit_counts.at("ID_in_13");
			if(in_gate) ++hit_counts.at("ID_in_gate");
		} else if(cableNumber>11150 && cableNumber<11155){
			++hit_counts.at("Veto");
		} else if(cableNumber>11154 && cableNumber<11186){
			++hit_counts.at("Calibration");
		} else if(cableNumber>15000 && cableNumber<15024){
			++hit_counts.at("Trigger");
		} else if(cableNumber>QB_OD_OFFSET && cableNumber<(QB_OD_OFFSET+MAXPMA)){
			++hit_counts.at("OD");
			if(in_13us) ++hit_counts.at("OD_in_13");
			if(in_gate) ++hit_counts.at("OD_in_gate");
		} else {
			++hit_counts.at("Unknown");
		}
	}
	
	std::cout<<"TQAReal Hit counts:\n"
	         <<"\tID Hits: "<<hit_counts.at("ID")<<" of which "<<hit_counts.at("ID_in_gate")
	                        <<" were in-gate and "<<hit_counts.at("ID_in_13")<<" were in a 1.3us window\n"
	         <<"\tOD Hits: "<<hit_counts.at("OD")<<" of which "<<hit_counts.at("OD_in_gate")
	                        <<" were in-gate and "<<hit_counts.at("OD_in_13")<<" were in a 1.3us window\n"
	         <<"\tVeto Hits: "<<hit_counts.at("Veto")<<"\n"
	         <<"\tCalibration Hits: "<<hit_counts.at("Veto")<<"\n"
	         <<"\tTrigger Hits: "<<hit_counts.at("Veto")<<"\n"
	         <<"\tUknonwn Hits: "<<hit_counts.at("Veto")<<std::endl;
	
	return true;
}

// TODO add Print methods for other ROOT branches: TQLIST, TQRAWSK,... any others?
// TODO also add to the PrintHits method to compare all sources

bool PrintEvent::PrintHits(){
	// we have many sources of hits: let's compare
	
	TQReal* myTQReal=nullptr;
	get_ok = myTreeReader->Get("TQREAL", myTQReal);
	if(!get_ok || myTQReal==nullptr){
		Log(m_unique_name+" failed to get TQREAL from Tree",v_error,verbosity);
		return false;
	}
	TQReal* myTQAReal=nullptr;
	get_ok = myTreeReader->Get("TQAREAL", myTQAReal);
	if(!get_ok || myTQAReal==nullptr){
		Log(m_unique_name+" failed to get TQAREAL from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting Hits (various sources)\n"<<std::endl;
	
	/* COMMON BLOCKS:
	skq_                       skqa_         : units = p.e. or pC? XXX
	{   XXX NOTE: index of qisk/qask are (PMT number-1) from skchnl_.ihcab
		int    nqisk;          nqask         : # hit ID PMTs ‖ # OD hits [-5,+1]us from T0 (see notes below)
		float  qismsk;         qasmsk        : total ID charge ‖ total charge of OD hits
		float  qimxsk;         qamxsk        : max charge an ID PMT ‖ max charge on an OD PMT
		int    mxqisk;         mxqask        : PMT# for max charge PMT ‖ PMT# of max charge PMT
		float  qisk[nqisk];    qask[nqask]   : charges of hits ‖ charges of hits
	}
	// sktqC.h says nqask is "num OD hits within [-5,1]us from trigger T0",
	// while nhitaz is "all OD hits (not just near TRG trigger)".
	// however, it does not say anything about nqisk being within a limited time range,
	// but then says nqiskz is "# all ID hits (not just near sw trigger",
	// so maybe nqisk is also within [-5,1]us from T0?? FIXME figure this out!
	
	skt_                      skta_          : units = [ns]? XXX, T0 ref: IT0XSK
	{   XXX NOTE: index of tisk/task are (PMT number-1) from skchnl_.ihcab
		float  timnsk;        tamnsk         : min hit time
		float  timxsk;        tamxsk         : max hit time
		int    mntisk;        mntask         : PMT# for min hit time PMT
		int    mxtisk;        mxtask         : PMT# for max hit time PMT
		float  tisk[nqisk];   task[nqask]    : times of hits
	}
	
	skchnl_                                  : units = TDC/ADC counts, T0 ref: IT0XSK
	{   FIXME are array indexes here (PMT number-1) from ihcab??
		int    iab[MAXPM];                   : channel number (1=chA, 2=chB) ... maybe ATM DAQ related?
		int    iqabsk[2][MAXPM];             : raw Q counts of A & B channels (array size? nqisk?)
		int    itabsk[2][MAXPM];             : raw T counts of A & B channels
		int    ihcab[nqisk];                 : PMT number
	}
	
	sktqz_                         sktqaz_          : units = , T0 ref: IT0XSK
	{   XXX NOTE: index of all arrays in sktqz are linear, NOT (PMT number-1)
		int    nqiskz;             nhitaz           : # ID hits ‖ # OD hits ("all, not just near TRG trigger")
		int    ihtiflz[nqiskz];    ihtflz[nhitaz]   : hit flags*
		int    icabiz[nqiskz];     icabaz[nhitaz]   : PMT#
		int    itiskz[nqiskz];     itaskz[nhitaz]   : raw T count
		int    iqiskz[nqiskz];     iqaskz[nhitaz]   : raw Q count (bits 0-10), hit flags (bits 11-15)**
		float  tiskz[nqiskz];      taskz[nhitaz]    : hit T [ns]
		float  qiskz[nqiskz];      qaskz[nhitaz]    : hit charge [pe]
		[no ID variable]           ihacab[nqask]    : OD PMT# without 20k offset (SKI-III).
		// note ihacab in sktqz is presumably OD version of ihcab in skchnl_ as there is no skchnla_!
		// its size is different than the rest of the arrays in sktqz_!
	}
	// * for ihtiflz/ihtflz see Constants.cpp::GetHitFlagNames
	// ** for iqiskz/iqaskz see Constants.cpp::GetHitChargeAndFlags
	
	// "information before bad channel mask": XXX
	// seem to be the same variables as in other commons just with extra entries for bad PMTs...?
	rawtqinfo_                           : units = mixed, T0 ref: ???
	{   XXX NOTE: index of all arrays in sktqz are linear, NOT (PMT number-1)
		int    nqisk_raw;                : # ID hits
		int    nhitaz_raw;               : # OD hits
		float  qbuf_raw[30*MAXPM];       : ID hit charges [units: p.e. or pC? XXX]
		float  qaskz_raw[30*MAXPMA];     : OD hit charges [units: p.e. or pC? XXX]
		float  tbuf_raw[30*MAXPM];       : ID hit times [ns]
		float  taskz_raw[30*MAXPMA];     : OD hit times [ns]
		int    icabbf_raw[30*MAXPM];     : ID PMT# (lower 16 bits) + hit flags (upper 16 bits) (XXX ==ihtiflz?)
		int    icabaz_raw[30*MAXPMA];    : OD PMT# (lower 16 bits) + hit flags (upper 16 bits) (XXX ==ihtflz?)
		float  pc2pe_raw;                : conversion factor pC to p.e. (possibly "obsolete since SK2"?)
		int    itiskz_raw[30*MAXPM];     : ID hits raw T counts
		int    itaskz_raw[30*MAXPMA];    : OD hits raw T counts
		int    iqiskz_raw[30*MAXPM];     : ID hits raw Q & hit flags
		int    iqaskz_raw[30*MAXPMA];    : OD hits raw Q & hit flags
	}
	*/
	
	// presumably these common blocks are all filled by SKRAWREAD/SKREAD....? FIXME verify?
	
	// TODO add notes explaining the differences
	std::cout<<"Num ID hits:\n"
	         <<"\tTQReal.nhits: "<<myTQReal->nhits<<"\n" // TODO == nqisk or nqiskz?
	         <<"\tskq_.nqisk (in-gate only): "<<skq_.nqisk<<"\n" // XXX check these are all and only in-gate hits
	         <<"\tsktqz_.nqiskz (all hits): "<<sktqz_.nqiskz<<"\n" // == skq_.nqisk
	         <<"\trawtqinfo_.nqisk_raw (all hits): "<<rawtqinfo_.nqisk_raw<<"\n"
	         <<"Total ID charge (skq_.qismsk): "<<skq_.qismsk<<" [UNITS?]\n" // all these are floats...
	         <<"ID PMT "<<skq_.mxqisk<<" saw the most charge of any ID PMT, with "<<skq_.qimxsk<<" [UNITS?].\n"
	         <<"ID PMT "<<skt_.mntisk<<" had the earliest hit of any ID PMT at "<<skt_.timnsk<<" [UNITS?]\n"
	         <<"ID PMT "<<skt_.mxtisk<<" had the last hit of any ID PMT at "<<skt_.timxsk<<" [UNITS?]\n";
	         // note these 'max/min, earliest/latest' variables are only for in-gate hits
	         // so e.g. skt_.timnsk == skt_.tisk[skchnl_.ihcab[0]]
	
	if(myTQReal->nhits==0){
		std::cout<<"TQReal had no ID hits"<<std::endl;
	} else {
		std::cout<<"first 5 ID hits (TQReal: bad channel mask applied, not just near sw trigger?):\n";
		for (int i = 0, j=0; i<myTQReal->nhits; ++i){
			int cableNumber = myTQReal->cables.at(i);
			int hitflags = cableNumber >> 16;
			cableNumber = cableNumber & 0xFFFF;
			if(cableNumber<=0 && cableNumber>MAXPM) continue;
			std::string hitflagstring;
			GetHitFlagNames(hitflags, &hitflagstring);
			std::cout<<"ID Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<"\n"
				     <<"\tT: TQReal.t: "<<myTQReal->T.at(j)<<" [ns]\n"
				     <<"\tQ: TQReal.q: "<<myTQReal->Q.at(j)<<" [units?]\n"
				     <<"\fhit flags: "<<hitflagstring<<"\n";
			++j;
			if(j>5) break;
		}
	}
	
	if(skq_.nqisk==0){
		std::cout<<"skq_ had no ID hits"<<std::endl;
	} else {
		std::cout<<"first 5 ID hits (skq_,skt_: bad channel mask applied? near sw trigger only?):\n";
		for(int i=0, j=0; i<skq_.nqisk; ++i){
			int cableNumber = skchnl_.ihcab[i];
			if(cableNumber<=0 && cableNumber>MAXPM) continue;
			int iab = skchnl_.iab[cableNumber-1];
			// indexing seems like this from $SKOFL_ROOT/examples/skrd/dump_tq.F
			std::cout<<"ID Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<", iab: "<<iab<<" = "<<((iab==1) ? "A" : "B")
				     <<"\tT: skchnl_.itabsk: "<<skchnl_.itabsk[iab-1][cableNumber-1] // seems to be 0
				     <<", skt_.tisk: "<<skt_.tisk[cableNumber-1]<<"\n"
				     <<"\tQ: skchnl_.iqabsk: "<<skchnl_.iqabsk[iab-1][cableNumber-1] // seems to be 0
				     <<", skq_.qisk: "<<skq_.qisk[cableNumber-1]<<"\n";
			++j;
			if(j>5) break;
		}
	}
	
	if(sktqz_.nqiskz==0){
		std::cout<<"sktqz had no ID hits"<<std::endl;
	} else {
		std::cout<<"first 5 ID hits (sktqz_: bad channel mask applied? not just near sw trigger):\n";
		for(int i=0, j=0; i<sktqz_.nqiskz; ++i){
			int cableNumber = sktqz_.icabiz[i];
			if(cableNumber<=0 && cableNumber>MAXPM) continue;
			// XXX note indices are [i] NOT [cableNumber-1]!
			std::string ihtiflz_flagstring="";
			GetHitFlagNames(sktqz_.ihtiflz[i], &ihtiflz_flagstring);
			std::string iqiskz_flagstring="";
			int iqiskz_q_counts, iqiskz_flags;
			GetHitChargeAndFlags(sktqz_.iqiskz[i], iqiskz_q_counts, iqiskz_flags, &iqiskz_flagstring);
			std::cout<<"ID Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<"\n"
				     <<"\tT: sktqz_.itiskz: "<<sktqz_.itiskz[i]
				     <<" [ticks], sktqz_.tiskz: "<<sktqz_.tiskz[i]<<" [ns]\n"
				     <<"\tQ: sktqz_.iqiskz: "<<iqiskz_q_counts<<
				     " [cts], sktqz_.qiskz: "<<sktqz_.qiskz[i]<<" [UNITS?]\n" // XXX units?
				     <<"\thit flags: sktqz_.iqiskz: "<<iqiskz_flagstring<<"\n"     // < FIXME these variables
				     <<"\thit flags: sktqz_.ihtiflz: "<<ihtiflz_flagstring<<"\n";  // < do not agree!
			++j;
			if(j>5) break;
		}
	}
	
	if(rawtqinfo_.nqisk_raw==0){
		std::cout<<"rawtqinfo had no ID hits"<<std::endl;
	} else {
		std::cout<<"first 5 ID hits (rawtqinfo_: before bad channel mask, not just near sw trigger?):\n";
		for(int i=0, j=0; i<rawtqinfo_.nqisk_raw; ++i){
			int cableNumber = rawtqinfo_.icabbf_raw[i] & 0xFFFF; // FIXME Constants function for icabbf/icabaz?
			if(cableNumber<=0 || cableNumber>MAXPM) continue;
			// XXX note indices are [i] NOT [cableNumber-1]!
			std::string icabbf_flagstring="";
			GetHitFlagNames((rawtqinfo_.icabbf_raw[i] >> 16), &icabbf_flagstring);
			std::string iqiskz_flagstring="";
			int iqiskz_q_counts, iqiskz_flags;
			GetHitChargeAndFlags(rawtqinfo_.iqiskz_raw[i], iqiskz_q_counts, iqiskz_flags, &iqiskz_flagstring);
			std::cout<<"ID Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<"\n"
				     <<"\tT: rawtqinfo_.itiskz_raw: "<<rawtqinfo_.itiskz_raw[i] // == sktqz.itiskz
				     <<", rawtqinfo_.tbuf_raw: "<<rawtqinfo_.tbuf_raw[i]<<"\n"  // == sktqz.tiskz
				     <<"\tQ: rawtqinfo_.iqiskz_raw: "<<(rawtqinfo_.iqiskz_raw[i] & 0x07FF) // bits 0-10
				     <<", rawtqinfo_.qbuf_raw: "<<rawtqinfo_.qbuf_raw[i]<<"\n"
				     <<"\thit flags: rawtqinfo_.iqiskz_raw: "<<iqiskz_flagstring<<"\n"
				     <<"\thit flags: rawtqinfo_.icabbf_raw: "<<icabbf_flagstring<<"\n";
			++j;
			if(j>5) break;
		}
	}
	
	//  =======================
	
	std::cout<<"Num OD hits:\n"
	         <<"\tTQAReal.nhits: "<<myTQAReal->nhits<<"\n"
	         <<"\tskqa_.nqask: "<<skqa_.nqask<<"\n"
	         <<"\tsktqaz_.nhitaz: "<<sktqaz_.nhitaz<<"\n"
	         <<"\trawtqinfo_.nhitaz_raw: "<<rawtqinfo_.nhitaz_raw<<"\n" // XXX != sktqaz_.nhitaz unlike ID hits!
	         <<"Total OD charge (skqa_.qasmsk): "<<skqa_.qasmsk<<" [UNITS?]\n"
	         <<"OD PMT "<<skqa_.mxqask<<" saw the most charge of any OD PMT, with "<<skqa_.qamxsk<<" [UNITS?]\n"
	         <<"OD PMT "<<skta_.mntask<<" had the earliest hit of any OD PMT at "<<skta_.tamnsk<<" [ns?]\n"
	         <<"OD PMT "<<skta_.mxtask<<" had the last hit of any OD PMT at "<<skta_.tamxsk<<" [UNITS?]\n";
	
	if(myTQAReal->nhits==0){
		std::cout<<"TQAReal had no OD hits"<<std::endl;
	} else {
		std::cout<<"first 5 ID hits (TQAReal: bad channel mask applied, not just near sw trigger?):\n";
		for (int i = 0, j=0; i<myTQAReal->nhits; ++i){
			int cableNumber = myTQAReal->cables.at(i);
			int hitflags = cableNumber >> 16;
			cableNumber = cableNumber & 0xFFFF;
			if(cableNumber<=QB_OD_OFFSET || cableNumber>(QB_OD_OFFSET+MAXPMA)) continue;
			std::string hitflagstring;
			GetHitFlagNames(hitflags, &hitflagstring);
			std::cout<<"OD Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<"\n"
				     <<"\tT: TQAReal.t: "<<myTQAReal->T.at(j)<<" [ns]\n"
				     <<"\tQ: TQAReal.q: "<<myTQAReal->Q.at(j)<<" [units?]\n"
				     <<"\fhit flags: "<<hitflagstring<<"\n";
			++j;
			if(j>5) break;
		}
	}
	
	if(skqa_.nqask==0){
		std::cout<<"skqa had no OD hits"<<std::endl;
	} else {
		std::cout<<"first 5 OD hits (skqa_,skta_: bad channel mask applied? near sw trigger only?):\n";
		for(int i=0, j=0; i<skqa_.nqask; ++i){
			int cableNumber = sktqaz_.ihacab[i]; // note ihacab is in sktqaz_ not skchnla_, which does not exist!
			// note also that ihacab does not contain OD PMT offset, so we need to add it
			cableNumber += QB_OD_OFFSET;
			if(cableNumber<=QB_OD_OFFSET || cableNumber>(QB_OD_OFFSET+MAXPMA)) continue;
			// sktqC.h says:
			// NQASK  ; # of hits (>= # of PMTs)
			// QASK   ; Q of individual PMT's (earliest hit for each PMT)
			// suggesting that NQASK is larger than the number of elements in QASK.
			// so how many elements does QASK contain? ... 
			// $SKOFL_ROOT/examples/skrd/dump_tq.F dumps NQASK elements of task, qask, as does the mc ver
			// $SKOFL_ROOT/examples/summer-test/read_rfm.F does the same
			// $SKOFL_ROOT/src/skrd/tqask.F sets NQASK to size of QASK while filling QASK/TASK
			// $SKOFL_ROOT/lowe/sklowe/lfnhita.F counts anti-hits in some timing window
			// and actually scans MAXPMA elements for task values in the timing window, suggesting either
			// MAXPMA elements, or at least suitable default values (e.g. 0) for unused elements
			// $SKOFL_ROOT/examples/skrd/sample_read.F loops MAXPMA times, printing elements where qisk[i]!=0
			// (though it also prints task[i] and qask[i], not task[ihacab[i]] and qask[ihacab[i]])
			std::cout<<"OD Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<"\n"
				     <<"\tT: skta_.task: "<<skta_.task[cableNumber-1]<<"\n"  // XXX garbage numbers? (NaN!)
				     <<"\tQ: skqa_.qask: "<<skqa_.qask[cableNumber-1]<<"\n"; // XXX ==0?!!
			++j;
			if(j>5) break;
		}
	}
	
	if(sktqaz_.nhitaz==0){
		std::cout<<"sktqaz had no OD hits"<<std::endl;
	} else {
		std::cout<<"first 5 OD hits (sktqaz_: bad channel mask applied? not just near sw trigger):\n";
		for(int i=0, j=0; i<sktqaz_.nhitaz; ++i){
			int cableNumber = sktqaz_.icabaz[i];
			if(cableNumber<=QB_OD_OFFSET || cableNumber>(QB_OD_OFFSET+MAXPMA)) continue;
			// XXX note indices are [i] NOT [cableNumber-1]!
			std::string ihtflz_flagstring="";
			GetHitFlagNames(sktqaz_.ihtflz[i], &ihtflz_flagstring);
			std::string iqaskz_flagstring="";
			int iqaskz_q_counts, iqaskz_flags;
			GetHitChargeAndFlags(sktqaz_.iqaskz[i], iqaskz_q_counts, iqaskz_flags, &iqaskz_flagstring);
			std::cout<<"OD Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<"\n"                             // XXX also all seem to be 0.
				     <<"\tT: sktqaz_.itaskz: "<<sktqaz_.itaskz[i]               // have we set an SKREAD option
				     <<" [ticks], sktqaz_.taskz: "<<sktqaz_.taskz[i]<<" [ns]\n" // which is ignoring OD hits?
				     <<"\tQ: sktqaz_.iqaskz: "<<iqaskz_q_counts
				     <<" [cts], sktqaz_.qaskz: "<<sktqaz_.qaskz[i]<<"\n"
				     <<"\thit flags: sktqaz_.iqaskz: "<<iqaskz_flagstring<<"\n"
				     <<"\thit flags: sktqaz_.ihtflz: "<<ihtflz_flagstring<<"\n";
			++j;
			if(j>5) break;
		}
	}
	
	if(rawtqinfo_.nhitaz_raw==0){
		std::cout<<"rawtqinfo had no OD hits"<<std::endl;
	} else {
		std::cout<<"first 5 OD hits (before bad ch masking):\n";
		for(int i=0, j=0; i<rawtqinfo_.nhitaz_raw; ++i){
			int cableNumber = rawtqinfo_.icabbf_raw[i] & 0xFFFF;
			if(cableNumber<=QB_OD_OFFSET || cableNumber>(QB_OD_OFFSET+MAXPMA)) continue;
			std::cout<<"OD Hit "<<j<<"\n"
				     <<"\tPMT: "<<cableNumber<<"\n"
				     <<"\tT: rawtqinfo_.itaskz_raw: "<<rawtqinfo_.itaskz_raw[i]
				     <<", rawtqinfo_.taskz_raw: "<<rawtqinfo_.taskz_raw[i]<<"\n"
				     <<"\tQ: rawtqinfo_.iqaskz_raw: "<<(rawtqinfo_.iqaskz_raw[i] & 0x07FF)
				     <<", rawtqinfo_.qaskz_raw: "<<rawtqinfo_.qaskz_raw[i]<<"\n";
			++j;
			if(j>5) break;
		}
	}
	
	// =========
	
	// probably don't need to keep this code once we verify the indexing
	if(sktqz_.nqiskz==0){
		std::cout<<"can't do check of sktqz vs skt because no ID hits"<<std::endl;
	} else {
		std::cout<<"checking first 10 in-gate hits ID in sktqz_[i] vs skt_[ihcab[i]]:\n";
		for(int i=0, j=0; i<sktqz_.nqiskz; ++i){
			if((sktqz_.ihtiflz[i] & 0x01)==0) continue; // not in-gate
			int cableNumber = sktqz_.icabiz[i];
			int cableNumber2 = skchnl_.ihcab[j];
			if(cableNumber!=cableNumber2){
				std::cout<<"in-gate ID Hit "<<j<<"\n"
					     <<"\tPMT: "<<cableNumber<<" != "<<cableNumber2<<"\n"; // FIXME still getting this
			} else {
				std::string ihtiflz_flagstring="";
				GetHitFlagNames(sktqz_.ihtiflz[i], &ihtiflz_flagstring);
				std::string iqiskz_flagstring="";
				int iqiskz_q_counts, iqiskz_flags;
				GetHitChargeAndFlags(sktqz_.iqiskz[i], iqiskz_q_counts, iqiskz_flags, &iqiskz_flagstring);
				int iab = skchnl_.iab[cableNumber2-1];
				std::cout<<"in-gate ID Hit "<<j<<"\n"
				         <<"\tPMT: "<<cableNumber<<"\n"
				         <<"\tT: sktqz_.itiskz: "<<sktqz_.itiskz[i]
				         <<" [ticks], sktqz_.tiskz: "<<sktqz_.tiskz[i]<<" [ns]\n"
				         <<"\tQ: sktqz_.iqiskz: "<<iqiskz_q_counts
				         <<" [cts], sktqz_.qiskz: "<<sktqz_.qiskz[i]<<"\n"
				         <<"\thit flags: sktqz_.iqiskz: "<<iqiskz_flagstring<<", "
				         <<"\thit flags: sktqz_.ihtiflz: "<<ihtiflz_flagstring<<"\n"
				         <<"vs\n"
				         <<"\tPMT: sktqz_.ihacab: "<<cableNumber2<<"\n"
				         <<"\tT: skchnl_.itabsk: "<<skchnl_.itabsk[iab-1][cableNumber2-1]
				         <<", skt_.tisk: "<<skt_.tisk[cableNumber2-1]<<"\n"
				         <<"\tQ: skchnl_.iqabsk: "<<skchnl_.iqabsk[iab-1][cableNumber2-1]
				         <<", skq_.qisk: "<<skq_.qisk[cableNumber2-1]<<"\n";
			}
			++j;
			if(j>10) break;
		}
	}
	
	if(sktqaz_.nhitaz==0){
		std::cout<<"can't do check of sktqaz vs skta because no OD hits"<<std::endl;
	} else {
		std::cout<<"checking first 10 in-gate OD hits in sktqaz_[i] vs skta_[ihacab[i]]:\n";
		for(int i=0, j=0; i<sktqaz_.nhitaz; ++i){
			if((sktqaz_.ihtflz[i] & 0x02)==0) continue; // not in-gate
			int cableNumber = sktqaz_.icabaz[i];
			int cableNumber2 = sktqaz_.ihacab[j];
			cableNumber2 += QB_OD_OFFSET;
			if(cableNumber!=cableNumber2){
				std::cout<<"in-gate OD Hit "<<j<<"\n"
					     <<"\tPMT: "<<cableNumber<<" != "<<cableNumber2<<"\n";
			} else {
				std::string ihtflz_flagstring="";
				GetHitFlagNames(sktqaz_.ihtflz[i], &ihtflz_flagstring);
				std::string iqaskz_flagstring="";
				int iqaskz_q_counts, iqaskz_flags;
				GetHitChargeAndFlags(sktqaz_.iqaskz[i], iqaskz_q_counts, iqaskz_flags, &iqaskz_flagstring);
				std::cout<<"in-gate OD Hit "<<j<<"\n"
					     <<"\tPMT: sktqaz_.icabaz: "<<cableNumber<<"\n"
					     <<"\tT: sktqaz_.itaskz: "<<sktqaz_.itaskz[i]
					     <<" [ticks], sktqaz_.taskz: "<<sktqaz_.taskz[i]<<" [ns]\n"
					     <<"\tQ: sktqaz_.iqaskz: "<<iqaskz_q_counts
					     <<"[cts], sktqaz_.qaskz: "<<sktqaz_.qaskz[i]<<"\n"
					     <<"\thit flags: sktqaz_.iqaskz: "<<iqaskz_flagstring<<", "
					     <<"\thit flags: sktqaz_.ihtflz: "<<ihtflz_flagstring<<"\n"
					     <<"\n"
					     <<"\tPMT: sktqaz_.ihacab: "<<cableNumber2<<"\n"
					     <<"\tT: skta_.task: "<<skta_.task[cableNumber2-1]<<"\n"
					     <<"\tQ: skqa_.qask: "<<skqa_.qask[cableNumber2-1]<<"\n";
			}
			++j;
			if(j>10) break;
		}
	}
	
	// checking number of entries in skqa_/skta_ to see if it's really < nqask
	if(skqa_.nqask==0){
		std::cout<<"can't do check of number of elements in skqa/skta because no OD hits"<<std::endl;
	} else {
		int n_set_OD_hits=0;
		int lastCableNumber=sktqaz_.ihacab[0];
		int print_multiple_hits=0;
		for(int i=0; i<MAXPMA; ++i){
			int cableNumber = sktqaz_.ihacab[i];
			if(cableNumber==0) continue; // important
			cableNumber += QB_OD_OFFSET;
			if(skta_.task[cableNumber-1]==0 || skqa_.qask[cableNumber-1]==0) continue;
			++n_set_OD_hits;
			if(cableNumber==lastCableNumber && print_multiple_hits<5){
				std::cout<<"found >1 hit on OD PMT "<<cableNumber<<std::endl;
				++print_multiple_hits;
			}
			lastCableNumber = cableNumber;
		}
		std::cout<<"counted "<<n_set_OD_hits<<" entries in skta_.task/skqa_.qask, vs skqa_.nqask = "
			     <<skqa_.nqask<<std::endl;
		// "counted 12 entries in skta_.task/skqa_.qask, vs skqa_.nqask = 17" ... XXX
	}
	
	//==========
	
	return true;
}

bool PrintEvent::PrintTriggerInfo(){
	
	std::cout<<"\nPrinting Trigger Info (commons)\n"<<std::endl;
	
	std::cout<<"skheadqb_.it0sk: "<<skheadqb_.it0sk<<"\n"
	         <<"skheadqb_.nevswsk: "<<skheadqb_.nevswsk<<"\n"
	         <<"skruninf_.softtrg_pre_t0: "<<skruninf_.softtrg_pre_t0[skheadqb_.nevswsk]<<"\n"
	         <<"skruninf_.softtrg_post_t0: "<<skruninf_.softtrg_post_t0[skheadqb_.nevswsk]<<"\n"
	         <<"skruninf_.softtrg_t0_offset: "<<skruninf_.softtrg_t0_offset[skheadqb_.nevswsk]<<"\n";
	         // seems these skruninf_ variables are not set in rfm files: all 0
	
	std::cout<<"Trigger gate spanned: "
	         <<((skheadqb_.it0sk-skruninf_.softtrg_pre_t0[skheadqb_.nevswsk])*COUNT_PER_NSEC)
	         <<"ns to "
	         <<((skheadqb_.it0sk+skruninf_.softtrg_post_t0[skheadqb_.nevswsk])*COUNT_PER_NSEC)
	         <<" with trigger signal at "
	         <<((skheadqb_.it0sk-skruninf_.softtrg_t0_offset[skheadqb_.nevswsk])*COUNT_PER_NSEC)<<"\n";
	
	// TODO add more info such as trigger type? hmm... don't want to overlap the other printers....
	return true;
}

bool PrintEvent::PrintSubTriggers(bool verbose){
	
	Header* myHeader=nullptr;
	get_ok = myTreeReader->Get("HEADER", myHeader);
	if(!get_ok || myHeader==nullptr){
		Log(m_unique_name+" failed to get HEADER from Tree",v_error,verbosity);
		return false;
	}
	
	return true; // TODO
	
	std::cout<<"\nPrinting subtriggers\n"<<std::endl;
	
	long it0sk = myHeader->t0;        // skheadqb_.it0sk;
	// IT0SK   ; Original T0 of the event
	// IT0XSK  ; T0 of the event for ITISKZ,ITASKZ,ITABSK,TISK,IHCAB,TASK,IHACAB,...
	
	// IT0SK marks the time of the trigger that resulted in a given readout, but there may have been
	// other trigger conditions met within that readout window (e.g. the SLE NHITS
	// condition is met at some time within the 40us readout window of an LE primary trigger)
	// such embedded subtriggers may be later extracted and saved as unique events, but these events
	// will share the same IT0SK. To identify their unique trigger times they will be assigned unique
	// values of IT0XSK - the time of the delayed trigger. (for the primary readout, IT0SK == IT0XSK).
	
	// we can retrieve the times of "subtriggers" by calling `get_sub_triggers`
	// and then call `set_timing_gate(it0sk+t0_sub(i))` (where t0_sub[i] is the time of subtrigger i)
	// followed by `skread(-LUN)` to reclculate the following variables:
	// it0xsk (new it0sk), tisk, task, qisk, qask, and hit flags.
	// note it0sk, tiskz, taskz are not updated (i think?)
	
	// so to process all subtriggers the following steps are required:
	// 1. runinf_ to retrieve trigger settings for the given run
	// 2. softtrg_set_cond_ to configure these settings (1&2 may only be required if you want to change trigger settings?)
	// 3. skrawread, to get next event
	// 4. get_sub_triggers(trigid, ntrig, t0_sub, MAX_TRIG) to get subtriggers of given type
	//    trigid: input int, format of idtgsk, which trigger to look for (one bit set only? TODO check)
	//    ntrig: output - # subtriggers found
	//    t0_sub: output array of subtrigger times
	//    MAX_TRIG: max # to search for
	//    (interally calls sofftrg_inittrgtbl to apply software trigger)
	// 5. set_timing_gate(it0sk+t0_sub(i)) to set next t0
	// 6. skread(-LUN) to re-load event and recalculate times and other variables
	// 7. optionally update PREVT0 and tdiff branches.... ?? unclear what this is being set to...*
	// 8. optionally use `delete outside_hits` to remove hits outside a 1.3us window around IT0XSK
	//    to emulate an SKI-III readout; if you don't call this and save subtriggers as new events,
	//    you'll end up with many events that essentially hold the same hits, so its wasteful....
	//    although mue_decay.F example does exactly this!
	
	// *see $SKOFL_ROOT/examples/lowe/mue_decay.F for example that selects parent mu & decay-e,
	// applies lfmufit to primary SHE triggers and saves event, then applies lfallfit_sk4_data_
	// to subtriggers & saves each subtrigger as a new event.
	// This also updates PREVT0 and ltdiff, but i don't follow the process....
	// apparently lomufit_sle also does this (10a-slereduction_lomufit_sle.pdf:10)
	// TODO
	
	// TODO  print first/last hit times in subtrigger windows or sth.
	/*
		// convert hit times in primary trigger window to hit time within a subtrigger
		// by subtracting [the time from 40us buffer readout start (IT0SK) to subtrigger start (IT0XSK)],
		// converted from TDC ticks to nanoseconds:
		hit_time -= (it0xsk-it0sk)/COUNT_PER_NSEC;
	*/
	
	return true;
}

bool PrintEvent::PrintLowEInfo(){
	LoweInfo* lowe=nullptr;
	get_ok = myTreeReader->Get("LOWE", lowe);
	if(!get_ok || lowe==nullptr){
		Log(m_unique_name+" failed to get LOWE from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting LOWE Branch\n"<<std::endl;
	
	if(lowe->bsenergy != 0){ // reconstruction not done
		std::cout<<"LOWE::bsenergy==0, appears LOWE branch is not populated"<<std::endl;
	}
	else if(lowe->bsenergy == 9999){
		std::cout<<"LOWE::bsenergy==9999, appears bonsai energy reconstruction failed"<<std::endl;
	}
	if(lowe->bsenergy!=0 && lowe->bsvertex[0] == 9999){
		// vertex and energy are reconstructed independently, both can succeed/fail independently
		std::cout<<"LOWE::bsvertex[0]==9999, appears bonsai position reconstruction failed"<<std::endl;
	}
	
	// TODO
	lowe->Dump();
	/*  from loweroot.h, see also skroot_loweC.h for common block equivalents
		Float_t  bsvertex[4];  // bonsai vertex:  x, y, z, t
		Float_t  bsresult[6];  // bonsai results: theta, phi, alpha, cos_theta, epsilon, like_t
		Float_t  bsdir[3];     // bonsai direction(x,y,z) based on bonsai vertex
		Float_t  bsgood[3];    // bonsai goodness: likelihood, goodness, ?
		Float_t  bsdirks;      // bonsai direction KS 
		Float_t  bseffhit[12]; // bonsai effective hits at fixed water transparencies
		Float_t  bsenergy;     // bonsai energy
		Int_t    bsn50;        // bonsai # of hits in 50 nsec after TOF subtraction
		Float_t  bscossun;     // bonsai cossun
		Float_t  clvertex[4];  // clusfit vertex: x, y, z, t
		Float_t  clresult[4];  // clusfit results: theta, phi, cos_theta, good_t
		Float_t  cldir[3];     // clusfit direction(x,y,z) based on clusfit vertex
		Float_t  clgoodness;   // clusfit goodness 
		Float_t  cldirks;      // clusfit direction KS
		Float_t  cleffhit[12]; // clusfit effective hits at fixed water transparencies
		Float_t  clenergy;     // clusfit energy
		Int_t    cln50;        // clusfit # of hits in 50 nsec after TOF subtraction
		Float_t  clcossun;     // clusfit cossun
		Int_t    latmnum;      // ATM number
		Int_t    latmh;        // ATM hit
		Int_t    lmx24;        // max 24
		Double_t ltimediff;    // time to the previous LE event (in raw data)
		Float_t  lnsratio;     // Noise-Signal ratio
		Float_t  lsdir[3];     // solar direction at the time (x,y,z)
		Int_t    spaevnum;     // event number of the parent muon 
		Float_t  spaloglike;   // spallation log likelihood
		Float_t  sparesq;      // spallation residual charge
		Float_t  spadt;        // spallation delta-T between parent muon
		Float_t  spadll;       // longitudinal distance
		Float_t  spadlt;       // traversal distance 
		Int_t    spamuyn;      // spallation muyn
		Float_t  spamugdn;     // spallation muon goodness
		Float_t  posmc[3];     // MC true vertex position
		Float_t  dirmc[2][3];  // MC true direction (1st and 2nd particles)
		Float_t  pabsmc[2];    // MC absolute momentum (1st & 2nd)
		Float_t  energymc[2];  // MC generated energy(s) (1st & 2nd)
		Float_t  darkmc;       // MC dark rate for generation
		Int_t    islekeep;     // SLE keep flag
		// below is added on 30-OCT-2007 y.t. (still version 1)
		Float_t  bspatlik;     // bonsai pattern likelihood for being single electron-like event
		Float_t  clpatlik;     // bonsai pattern likelihood for being single electron-like event
		Float_t  lwatert;      // water transparency value (at reconstruction?)
		Int_t    lninfo;       // # of extra low-e inforamation
		Int_t    linfo[255];   // extra low-e information (see skroot_lowe.h)
		// some of interest:
		linfo[9,10] dwall for clusfit, bonsai
		linfo[5,6] effwall for clusfit, bonsai
		linfo[7,8] effective # hits for clusfit, bonsai ("at given water transparency"?)
		linfo[25] bonsai clik (cluster likelihood to reject noise events near wall)
		linfo[26] bonsai ovaq (bonsai goodness^2 - dirks^2) - combined 1D vertex & angle goodness metric
		linfo[27-33] ariadne values (bsn20, hits, adir, amsg, aratio, anscat, acosscat)? FIXME
		// n.b. from skroot_lowe.h, some aliases are given to some linfo elements:
		linfo=61-65: software trigger
		integer swtrig, swtrig_thr(0:3)
		EQUIVALENCE (linfo(61), swtrig)
		EQUIVALENCE (linfo(62), swtrig_thr(0)) ! threshold of trigid = 0
	*/
	return true;
}

bool PrintEvent::PrintATMPDInfo(){
	AtmpdInfo* atmpd=nullptr;
	get_ok = myTreeReader->Get("ATMPD", atmpd);
	if(!get_ok || atmpd==nullptr){
		Log(m_unique_name+" failed to get ATMPD from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting ATMPD Branch\n"<<std::endl;
	
	// TODO
	atmpd->Dump();
	return true;
}

bool PrintEvent::PrintUpMuInfo(){
	UpmuInfo* upmu=nullptr;
	get_ok = myTreeReader->Get("UPMU", upmu);
	if(!get_ok || upmu==nullptr){
		Log(m_unique_name+" failed to get UPMU from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting UPMU Branch\n"<<std::endl;
	
	// TODO
	upmu->Dump();
	return true;
}

bool PrintEvent::PrintMuInfo(){
	MuInfo* muinfo=nullptr;
	get_ok = myTreeReader->Get("MU", muinfo);
	if(!get_ok || muinfo==nullptr){
		Log(m_unique_name+" failed to get UPMU from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting MU Branch\n"<<std::endl;
	
	// TODO
	muinfo->Dump();
	/* from loweroot.h
		Float_t  muentpoint[3];
		Float_t  mudir[3];
		Double_t mutimediff;
		Float_t  mugoodness;
		Float_t  muqismsk;
		Int_t    muyn;
		Int_t    mufast_flag;
		Int_t    muboy_status;         // muboy status
		Int_t    muboy_ntrack;         // number of tracks
		Float_t  muboy_entpos[10][4];  // up to 10 tracks
		Float_t  muboy_dir[3];         // common direction
		Float_t  muboy_goodness;       // goodness
		Float_t  muboy_length;         // track length
		Float_t  muboy_dedx[200];      // dE/dX histogram
		Float_t  mubff_entpos[3];      // bff entpos
		Float_t  mubff_dir[3];         // bff direction
		Float_t  mubff_goodness;       // bff goodness
		Int_t    muninfo;              // number of additional data in muinfo
		Int_t    muinfo[255];          // additional data
	*/
	return true;
}

bool PrintEvent::PrintSLEInfo(){
	SLEInfo* sleinfo=nullptr;
	get_ok = myTreeReader->Get("SLE", sleinfo);
	if(!get_ok || sleinfo==nullptr){
		Log(m_unique_name+" failed to get SLE from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting SLE Branch\n"<<std::endl;
	
	// TODO
	sleinfo->Dump();
	/* from loweroot.h
		Float_t  wallcut;
		Int_t    nsel;
		Float_t  itbsvertex[4];
		Float_t  itbsresult[6];
		Float_t  itbsgood[3];
		Int_t    nbonsai;
		Float_t  itcfvertex[4];
		Float_t  itcfresult[4];
		Float_t  itcfgoodness;
		Int_t    nclusfit;
	*/
	return true;
}

bool PrintEvent::PrintPrevT0(){
	
	PrevT0* prevt0=nullptr;
	get_ok = myTreeReader->Get("PREVT0", prevt0);
	if(!get_ok || prevt0==nullptr){
		Log(m_unique_name+" failed to get PREVT0 from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting PREVT0 Branch\n"<<std::endl;
	
	std::cout<<"Time from this trigger to last trigger of each given type (PrevT0):\n"
	         <<"\tsw_event_no: "<<prevt0->sw_event_no<<"\n"; // skheadqb_.nevsk (XXX NOT skheadqb_.nevswsk!), == SoftwareTrigger::sw_event_no
	// XXX seems these will read 0 from the start of a run until at least there is a valid past trigger of that type?
	for(int i=0; i<NTRGTYPES; ++i){
		std::cout<<"\tfor trigger "<<TriggerIDToTrigger(i)<<"\n"
		         <<"\t\tprev_t0: "<<prevt0->prev_t0[i]<<"\n"
		         <<"\t\tprev_counter_32: "<<prevt0->prev_counter_32[i]<<"\n"
		         <<"\t\tprev_clock_ticks (clk48): ["<<prevt0->clk48[i][0]
		                                      <<", "<<prevt0->clk48[i][1]
		                                      <<", "<<prevt0->clk48[i][2]<<"]"<<std::endl;
	}
	
	// see $SKOFL_ROOT/lowe/sklowe/tdiff.F; tdiff(nt48sk, ltimediff) routine
	// (which for SKIV+ actually ignores its nt48sk argument and calls skroot_get_prevt0)
	// note from source code: "timediff of after trigger is not valid because it0xsk of AFT seems to be zero"!
	// also from lowe school slides: for SKIV+ ltimediff will roll over if time since last event is >2s!
	
	/* see also skprevt0C.h? how is this info related?
	/*   === previous timing information === */
	/*       iprev_ctr32(0:31)  ; previous HW TRG counter(32bits)*/
	/*       iprev_tdc  (0:31)  ; previous TDC count     (15bits)*/
	/*       iprev_clk48(3,0:31); previous TRG 48bit clk.(16bits x 3)*/
	/*   === timing difference information === */
	/*       idiff_ctr32(0:31)  ; hw trg counter difference*/
	/*       idiff_tdc  (0:31)  ; tdc count difference*/
	/*       rtdiff     (0:31)  ; timing difference      (micro sec.)*/
	/*                 rtdiff = idiff_ctr32 * (256./15.) + idiff_tdc / 1920.*/
	/*       iovlperr           ; should be 0. If not 0, reject this event.*/
	/*                            if this event was contained by the previous*/
	/*                            event, integer number (trig. id) will be set.*/
	/* *** for backward compatibility ; requested from Shiozawa-san,*/
	/*       mintrgid           ; trigger ID for minium T0*/
	/*       tdiffmin           ; minimum Tdiff (in micro sec.)*/
	/*                          ;   To define `minimum', only the following*/
	/*                          ;       trigger IDs are considered.*/
	/*                          ;                 */
	/*                          ;   LE,HE,OD,Periodic,After/CAL, */
	/*                          ;   VETO start, VETO stop,*/
	/*                          ;   SHE, Relic (after) trigger*/
	
	return true;
}


bool PrintEvent::PrintHWTriggerInfo(){
	
	TClonesArray* hwtriggers=nullptr;
	get_ok = myTreeReader->Get("HWTRGLIST", hwtriggers);
	if(!get_ok || hwtriggers==nullptr){
		Log(m_unique_name+" failed to get HWTRGLIST from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting HWTRGLIST Branch\n"<<std::endl;
	
	std::cout<<"This event had "<<hwtriggers->GetEntries()<<" HW triggers"<<std::endl;
	for(int i=0; i<hwtriggers->GetEntries(); ++i){  // skheadqb_.numhwsk? == SoftwareTrigger::nhwtrgs?
		
		HardwareTrigger* trg = dynamic_cast<HardwareTrigger*>((*hwtriggers)[i]);
		
		// make trigger timestamp
		tm trigdate = {0};
		int time_ymd = trg->time_ymd;   // year*10000+(mon+1)*100+day
		trigdate.tm_year = time_ymd / 10000;
		time_ymd = (time_ymd % 10000);
		trigdate.tm_mon = (time_ymd / 100) -1;
		time_ymd = (time_ymd%100);
		trigdate.tm_mday = time_ymd;
		
		int time_hms = trg->time_hms;   // hour*10000+min*100+sec
		trigdate.tm_hour = (time_hms / 10000);
		time_hms = (time_hms % 10000);
		trigdate.tm_min = (time_hms / 100);
		time_hms = (time_hms % 100);
		trigdate.tm_sec = time_hms;
		
		std::string trgtimestring;
		char trgbuffer[20];
		strftime(trgbuffer,20,"%F %T",&trigdate);
		trgtimestring = trgbuffer;
		
		std::cout<<"HW Trigger "<<i<<"\n"
		         <<"\tstatus: "<<trg->status<<"\n"
		         <<"\toccurred at: "<<trgtimestring<<"\n"
		         <<"\ttime: "<<trg->time<<"\n" // XXX what is this? some large integer?
		         <<"\tcounter32: "<<trg->counter_32<<"\n"
		         <<"\tclock ticks (clk48): ["<<trg->clk48[0]<<", "<<trg->clk48[1]<<", "<<trg->clk48[2]<<"]\n"
		         <<"\thw counter: "<<trg->hw_trg_ctr<<"\n"
		         <<"\tID: "<<trg->hw_trg_id<<"\n"
		         <<"\tas type:"<<TriggerIDToTrigger(trg->hw_trg_id)<<std::endl; // seems valid, usually 'Periodic'
		
	}
	
	/*
	from sktrg.h:
	sktrg_. common {
	  int    nsktrghw;                  // Total number of TRG segments
	  int    nruntrghw[nsktrghw];       // Run number
	  int    nevtrghw[nsktrghw];        // TRG 32bit counter for this TRG segment
	  int    idtrghw[nsktrghw];         // Trigger ID from TRG module
	  int    nclk48trghw[3][nsktrghw];  // 48bit clk
	  int    stattrghw[nsktrghw];       // Status
	}; */
	
	return true;
}

bool PrintEvent::PrintSWTriggerInfo(){
	
	SoftwareTrigger* swtrg=nullptr;
	get_ok = myTreeReader->Get("SOFTWARETRG",swtrg);
	if(!get_ok  || swtrg==nullptr){
		Log(m_unique_name+" failed to get SOFTWARETRG from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting SOFTWARETRG Branch\n"<<std::endl;
	
	std::cout<<"SW Trigger info:\n";
	std::cout<<"\tRun mode: "<<RunModeToName(swtrg->run_mode)<<std::endl;    // skhead_.mdrnsk? == RunInfo::mdrnsk?
	std::cout<<"\tsw_event_no: "<<swtrg->sw_event_no<<std::endl;             // skheadqb_.nevswsk - seems to be the same as nevsk...
	std::cout<<"\tsoftware trigger version: "<<swtrg->sw_trg_ver<<"\n"       // i guess this is a version number???
	         <<"\tnum hw triggers: "<<swtrg->nhwtrgs<<"\n"                   // skheadqb_.numhwsk? == HWTRGLIST->GetEntries(). don't ask me why this is in the SW trigger.
	         <<"\ttrg ID: "<<swtrg->sw_trg_id<<"\n"
	         <<"\tas type: "<<TriggerIDToTrigger(swtrg->sw_trg_id)<<"\n"     // seems valid
	         <<"\tMasked triggers: "<<swtrg->sw_trg_mask<<"\n"
	         <<"\tas types: "<<GetTriggerNames(swtrg->sw_trg_mask)<<"\n" // skheadqb_.trgmask? == Header::sw_trg_mask? is this a bitmask of ignored triggers??? FIXME don't think this is right.
	         <<"\tclock ticks (clk48): ["<<swtrg->clk48[0]<<", "<<swtrg->clk48[1]<<", "<<swtrg->clk48[2]<<"]\n"
	         <<"\tcounter32: "<<swtrg->counter_32<<"\n"
	         <<"\tt0: "<<swtrg->t0<<"\n"   // skheadqb_.it0sk, == Header::t0
	         <<"\tgate width: "<<swtrg->gate_width<<"\n"  // skheadqb_.gatewsk == Header::gate_width
	         <<"\tnhits: "<<swtrg->nhits<<std::endl;  // what kind? TQ? TQReal? ID? OD? ????
	
	// see also softtrg info in RunInfo class or skruninfo_ common block.
	// how are these different???
	
	// from sktrighitC.h:
	// is swtrg->nhits related to sktrighit_.ntgsk?
	// latter is commented as "number of hits within 40us window between 15000 and 15024"....?
	// sktrighit_ also contains 'icabtg', 'qtgsk' and 'ttgsk': cable #, Q and T of hits in trigger
	// (time comment = 'apply T0 correction'???)
	
	// how does this relate to the SWTRGLIST (SoftTrgList*) branch? not in rfm files...
	// defined in softtrgroot.h w/ members:
	/*
		std::vector<int> swtrg_Overlap;
		std::vector<int> swtrg_Version;
		std::vector<int> swtrg_ID;
		std::vector<int> swtrg_T0;
		std::vector<int> swtrg_GateStart;
		std::vector<int> swtrg_GateStop;
	*/
	// seems to be generated by TreeManager when converting from zbs, maybe...
	// equivalent common block: swtrglist_ (from softtrg_listC.h):
	/*
	  int    nswtrgs;
	  int    swtrgovlp[SWTRG_LISTSIZE];
	  int    swtrgver[SWTRG_LISTSIZE];
	  int    swtrgid[SWTRG_LISTSIZE];
	  int    swtrgt0[SWTRG_LISTSIZE];
	  int    swtrggs[SWTRG_LISTSIZE];
	  int    swtrgge[SWTRG_LISTSIZE];
	*/
	
	// for MC we have various variables in MCInfo class: labelled
	// 'DSTRGOFF: made for software trigger inside the detector simulation' .... whatever that means.
	/*
	Int_t    it0_offset;              // IT0_OFFSET of event: random clock position
	Int_t    it0sk_geantt0;           // geant0 in clock counts
	Int_t    numdsswtrgs;             // number of software triggers (1st : primary, 2nd- : delayed)
	Int_t    trigbit[MAXSWTRGS];      // sw trig id:(-1=NONE),(Bit flag:1=LE,2=HE,3=OD)
	Int_t    it0sk_temp[MAXSWTRGS];   // it0sk of trigger
	Float_t  prim_pret0[MAXSWTRGS];   // geant_t0 relative to software trigger
	Int_t    prim_trg[MAXSWTRGS];     // trigger candidate number
	*/
	
	return true;
}

bool PrintEvent::PrintHeaderInfo(){
	
	Header* myHeader=nullptr;
	get_ok = myTreeReader->Get("HEADER",myHeader);
	if(!get_ok  || myHeader==nullptr){
		Log(m_unique_name+" failed to get HEADER from Tree",v_error,verbosity);
		return false;
	}
	std::cout<<"\nPrinting HEADER Branch\n"<<std::endl;
	
	std::cout<<"SK Geom: "<<myHeader->sk_geometry  // this is 0 for rfm files!!
	         <<", (c.f. skheadg_.sk_geometry: "<<skheadg_.sk_geometry<<")"<<std::endl;
	
	//                   // from the HEADER branch        // from the common blocks
	std::cout<<    "Run: " <<myHeader->nrunsk             // skhead_.nrunsk;
	         <<", Subrun: "<<myHeader->nsubsk             // skhead_.nsubsk;
	         <<", Event: " <<myHeader->nevsk<<std::endl;  // skhead_.nevsk;
	
	// get run start - FIXME actually maybe this is event time... it is later than the time in RunInfo class
	tm rundate = {0};
	rundate.tm_year = myHeader->ndaysk[0];
	rundate.tm_mon = myHeader->ndaysk[1] - 1;
	rundate.tm_mday = myHeader->ndaysk[2];
	rundate.tm_hour = myHeader->ntimsk[0];
	rundate.tm_min = myHeader->ntimsk[1];
	rundate.tm_sec = myHeader->ntimsk[2];
	// we also have: myHeader->ntimsk[3] which holds a count of (1/100)th secs
	
	std::string timestring;
	/*
	// ctime formats as 'Thu Jul 6 19:42:45 2023'
	time_t runtime = mktime(&rundate);
	timestring = ctime(&runtime);   // format timestamp into a string
	timestring.pop_back();          // drop trailing newline
	*/
	// We can specify a custom format, such as '2023-07-06 19:42:45' with strftime
	char buffer[20];
	strftime(buffer,20,"%F %T",&rundate);
	std::string nsstring = std::to_string(myHeader->ntimsk[3]/100.); // trim off the leading '0.'
	timestring = std::string(buffer) + nsstring.substr(nsstring.find('.'),std::string::npos);
	
	std::cout << "The event occurred at " << timestring << std::endl;
	
	// "48-bit clock" from skheadC.h
	std::cout<<"clock ticks (nt48sk): ["<<myHeader->nt48sk[0]<<", "
	                                    <<myHeader->nt48sk[1]<<", "
	                                    <<myHeader->nt48sk[2]<<"]"<<std::endl;
	/*
	from lowe school: 48bit clock was used in SKI-III to get time difference between two events.
	(excludes SLE events). 1 tick = 20ns, tick count reset at the start of each run.
	NT48SK was unstable at start of SKIV, so instead a PREVT0 branch is added.
	*/
	
	int runMode = myHeader->mdrnsk;           // skhead_.mdrnsk == RunInfo::run_mode
	std::cout<<"This run was a "<<RunModeToName(runMode)<<" run"<<std::endl;
	
	std::cout<<"Triggers bits in this event: "<<GetTriggerNames(myHeader->idtgsk)<<"\n";  // skhead_.idtgsk
	
	std::cout<<"Event flags in this event: "<<GetEventFlagNames(myHeader->ifevsk)<<"\n";  // skhead_.ifevsk
	
	std::cout<<"Event header info:\n"
	         <<"\tit0sk [ticks]: "<<myHeader->t0<<"\n" // skheadqb_.it0sk XXX rfm files both same but negative??
	         <<"\tcounter_32: "<<myHeader->counter_32<<"\n"
	         <<"\tnum hardware triggers: "<<myHeader->ntrg<<"\n"
	         <<"\tsoftware trigger id: "<<myHeader->swtrg_id<<"\n"  // == SoftwareTrigger::sw_trg_id
	         <<"\tas type: "<<TriggerIDToTrigger(myHeader->swtrg_id)<<"\n" // FIXME not sure this is correct
	         <<"\tsoftware trigger mask: "<<myHeader->sw_trg_mask<<"\n"
	         <<"\tas types: "<<GetTriggerNames(myHeader->sw_trg_mask)<<"\n" // skheadqb_.trgmask? == SoftwareTrigger::sw_trg_mask == RunInfo::trigger_mask, what is this? are these the enabled triggers??
	         // FIXME if so, 'Normal Run' (rfm_run086730.000001.root) includes LINAC, LINAC MW, but not SHE, AFT?
	         <<"\ttrigger gate width: "<<myHeader->gate_width<<"\n";   // skheadqb_.gatewsk == SoftwareTrigger::gate_width
	         //<<"\t\"contents\": "<<myHeader->contents<<std::endl;     // skheadqb_.contsk? seriously wtf is this
	
	/* further members:
	Int_t ltcgps;   //      skheada_.ltcgps ; Local time clock at last GPS time
	Int_t nsgps;    //      skheada_.nsgps  ; GPS time (sec)
	Int_t nusgps;   //      skheada_.nusgps ; GPS time (usec)
	Int_t ltctrg;   //      skheada_.ltctrg ; Local time clock at TRG
	Int_t ltcbip;   //      skheada_.ltcbip ; Local time clock at end of BIP
	Int_t itdct0[4];//      skheada_.itdct0 ; TDC T0 (TRG) time for hut I (I = 1,4)
	Int_t iffscc;   //      skheada_.iffscc ; FSCC busy flags
	Int_t icalva;   //      skheada_.icalva ; Calibration constant version
	*/
	
	// note: data may also have the following branches, but they contain nothing interesting, i think
	//EventHeader* Head;
	//EventTrailer* Trail;
	
	return true;
}

bool PrintEvent::PrintMCInfo(){
	
	// for MC we have some variables in MCInfo that may be similar/relevant?
	MCInfo* mcinfo=nullptr;
	get_ok = myTreeReader->Get("MC",mcinfo);
	if(!get_ok  || mcinfo==nullptr){
		Log(m_unique_name+" failed to get HEADER from Tree",v_error,verbosity);
		return false;
	}
	
	std::cout<<"\nPrinting MC Branch\n"<<std::endl;
	
	// see mcinfo.h
	std::cout<<"SK Geom: (MCInfo): "<<mcinfo->ivmcp<<"\n" // MC version? apparently 1000+SK_geometry? unclear
	         <<"Run: "<<mcinfo->mcrun<<"\n"
	         <<"it0_offset: "<<mcinfo->it0_offset<<"\n"         // IT0_OFFSET of event (random)
	         <<"it0sk_geantt0: "<<mcinfo->it0sk_geantt0<<"\n"   // geant0 in clock counts
	         <<"dark rate [Hz]: "<<mcinfo->darkds<<"\n"
	         <<"trigger threshold: "<<mcinfo->trigds<<"\n"  // for which type???
	         <<"numdsswtrgs: "<<mcinfo->numdsswtrgs<<"\n";
	for(int i=0; i<mcinfo->numdsswtrgs; ++i){
		std::string triggernames;
		if(mcinfo->trigbit[i]==-1){
			triggernames="None";
		} else {
			if(mcinfo->trigbit[i] & 0x1) triggernames+=", LE";
			if(mcinfo->trigbit[i] & 0x2) triggernames+=", HE";
			if(mcinfo->trigbit[i] & 0x3) triggernames+=", OD";
			if(triggernames.size()){
				triggernames = triggernames.substr(1,std::string::npos);
			}
			triggernames = "[" + triggernames+" ]";
		}
		std::cout<<"mc sw trig "<<i<<":\n"
		         <<"\ttrigger candidate (prim_trg): "<<mcinfo->prim_trg[i]<<"\n" // meaning?
		         <<"\ttrigger it0sk [ticks]: "<<mcinfo->it0sk_temp[i]<<"\n"
		         <<"\ttrigger bits: "<<mcinfo->trigbit[i]<<triggernames<<"\n"
		         <<"\tgeant_t0 relative time: "<<mcinfo->prim_pret0[i]<<"\n";
	}
	
	/*
	// other variables "corresponding to MCPAR MCONV in skdetsim/dsprsv.F"
	Float_t  tcnvsk;  // skparm_.tcnvsk,  channel connversion coefficient for TDC count
	Float_t  qcnvsk;  // skparm_.qcnvsk,  channel connversion coefficient for ADC count
	Float_t  dthrsk;  // skparm_.dthrsk,  threshold count of discriminator
	Float_t  darkds;  // dark noise rate of PMT's (Hz)
	Float_t  tresds;  // timing resolution (nsec)
	Float_t  qresds;  // pulse height resolution (fraction of sigma)
	Float_t  twinds;  // hitsum pulse window (nsec)
	Float_t  trigds;  // trigger threshold (hit)
	Float_t  gateds;  // dynamic range (nsec)
	Float_t  beftds;  // ???
	Float_t  deadds;  // dead time for reject reflection (nsec)
	Float_t  sigwds;  // charge integration time (nsec)
	
	Int_t    mcninfo;    // where is this mapped? anything of interest in here?
	Int_t    mcinfo[mcninfo];
	*/
	
	return true;
}

bool PrintEvent::PrintRunInfo(){
	
	// RunInfo is stored in the data tree userinfo (data only)
	RunInfo* myRunInfo = nullptr;
	// we can retrieve it via the TreeManager...
	std::string treeReaderName = myTreeReader->GetName();
	int LUN = m_data->GetLUN(treeReaderName);
	if(LUN<0){
		Log(m_unique_name+" failed to find TreeReader "+treeReaderName+" in DataModel",v_error,verbosity);
		return false;
	}
	TreeManager* manager = GetTreeManager(LUN);
	myRunInfo = manager->GetRINFO();  // returns null if not present
	
	/*
	// or we get can get it from the TTree ourselves
	// (note: TreeManager uses iterators but has no loop, so only seems to check first entry)
	TList* RareList = (TList*) tree->GetUserInfo()->FindObject("Rare");
	if(RareList && RareList->GetSize()){
		TObject* rare =  RareList->At(0);
		if (strcmp(rare->GetName().Data(),"RUNINFO")==0){
			myRunInfo = (RunInfo*)rare;
		}
	}
	// other info in the UserInfo, which probably may be retreived from TreeManager in same way (check this):
	// Pedestal* PEDESTAL;
	// SlowControl* SLWCTRL;
	*/
	
	if(myRunInfo==nullptr){
		std::cout<<" No RunInfo"<<std::endl;
		return true;
	}
	
	std::cout<<"\nPrinting RunInfo\n"<<std::endl;
	
	// get run start timestamp
	tm runstartdate{*localtime(&(myRunInfo->run_start_time_sec))};
	
	/*
	// runinf_ common instead stores 'start_ymd_sk' and 'start_hms_sk' variables.
	// this is how the TreeManager generates those variables:
	int start_ymd_sk = runstartdate.tm_year * 10000 + (runstartdate.tm_mon + 1) * 100 + runstartdate.tm_mday;
	int start_hms_sk = runstartdate.tm_hour * 10000 +  runstartdate.tm_min      * 100 + runstartdate.tm_sec;
	
	// and code to split these back into components:
	int start_ymd_sk = skruninf_.start_ymd_sk;   // year*10000+(mon+1)*100+day
	runstartdate.tm_year = start_ymd_sk / 10000;
	start_ymd_sk = (start_ymd_sk % 10000);
	runstartdate.tm_mon = (start_ymd_sk / 100) -1;
	start_ymd_sk = (start_ymd_sk%100);
	runstartdate.tm_mday = start_ymd_sk;
	
	int start_hms_sk = skruninf_.start_hms_sk;   // hour*10000+min*100+sec
	runstartdate.tm_hour = (start_hms_sk / 10000);
	start_hms_sk = (start_hms_sk % 10000);
	runstartdate.tm_min = (start_hms_sk / 100);
	start_hms_sk = (start_hms_sk % 100);
	runstartdate.tm_sec = start_hms_sk;
	*/
	
	std::string runstarttimestring;
	/*
	// ctime formats as 'Thu Jul 6 19:42:45 2023'
	time_t runstarttime = mktime(&runstartdate);
	runstarttimestring = ctime(&runstarttime);   // format timestamp into a string
	runstarttimestring.pop_back();               // drop trailing newline
	*/
	// We can specify a custom format, such as '2023-07-06 19:42:45' with strftime
	char startbuffer[20];
	strftime(startbuffer,20,"%F %T",&runstartdate);
	char nsbuffer[20];
	snprintf(nsbuffer, 20, "%09ld", myRunInfo->run_start_time_nsec);
	runstarttimestring = std::string(startbuffer) + "." + nsbuffer;
	
	std::cout<<"Run start: "<<runstarttimestring<<std::endl;
	
	// repeat for run end...
	std::string runendtimestring="[not set]";
	if(myRunInfo->run_end_time_sec!=0){
		// from several data files it seems run_end_time_* variables are always 0?
		tm runenddate{*localtime(&(myRunInfo->run_end_time_sec))};   // skruninf_.end_time_sec_sk
		memset(startbuffer,0,sizeof(startbuffer));
		strftime(startbuffer,20,"%F %T",&runenddate);
		memset(nsbuffer,0,sizeof(nsbuffer));
		snprintf(nsbuffer, 20, "%09ld", myRunInfo->run_end_time_nsec);
		runendtimestring = std::string(startbuffer) + "." + nsbuffer;
	}
	std::cout<<"Run end: "<<runendtimestring<<std::endl;
	
	if(myRunInfo->run_end_time_sec!=0){
		int run_duration_sec_sk = myRunInfo->run_end_time_sec - myRunInfo->run_start_time_sec;
		int run_hours = (run_duration_sec_sk / (3600));
		run_duration_sec_sk = (run_duration_sec_sk % (3600));
		int run_mins = (run_duration_sec_sk / 60);
		int run_secs = (run_duration_sec_sk % 60);
		printf("Run duration: %d:%d:%d\n",run_hours, run_mins, run_secs);
	} else {
		std::cout<<"Run duration: ?"<<std::endl;
	}
	
	std::cout<<"Run Num: "<<myRunInfo->run_no<<"\n"
	         <<"Run title: "<<myRunInfo->run_title<<"\n" // skruninf_.run_title_sk (same as run type name)
	         <<"Run type: "<<RunModeToName(myRunInfo->run_mode)<<"\n" // skhead_.mdrnsk == Header::mdrnsk?
	         <<"Shift leader: "<<myRunInfo->shift_leader<<"\n" // skruninf_.shift_leader_sk
	         <<"Shift member: "<<myRunInfo->shift_members<<"\n" // skruninf_.shift_member_sk
	         <<"Run end message: "<<myRunInfo->run_end_message<<"\n"
	         <<"Run counter_32: "<<myRunInfo->counter_32<<"\n"; // XXX seems not set (in rfm files?)
	
	// from $SKOFL_ROOT/include/skruninf.h
	std::bitset<32> mask(myRunInfo->trigger_mask);  // skruninf_.softtrg_mask (!=idtgsk)
	for(int i=0; i<32; ++i){
		// TODO skip triggers for which trigger bit i is not a meaningful one
		std::cout<<"Trigger "<<TriggerIDToTrigger(i)<<"\n"
			     //<<(myRunInfo->detector[i] ? " (OD trigger) " : " (ID trigger) ")<<"\n" // skruninf_.softtrg_detector - literally only OD for the DO trigger.
			     <<"\tenabled? "<<(mask.test(i) ? "Y" : "N")<<"\n"  // skruninf_.softtrg_mask
			     <<"\tthreshold: "<<myRunInfo->thr[i]<<" [hits]\n"  // skruninf_.softtrg_thr
			     // FIXME print these in counts as well, to make it clear thats what units they're in.
			     <<"\tpre-trigger window: "<<(myRunInfo->pre_t0[i]*COUNT_PER_NSEC)<<" [ns]\n" // skruninf_.softtrg_pre_t0
			     <<"\tpost-trigger window: "<<(myRunInfo->post_t0[i]*COUNT_PER_NSEC)<<" [ns]\n" // skruninf_.softtrg_post_t0
			     <<"\tt0 offset: "<<(myRunInfo->t0_offset[i]*COUNT_PER_NSEC)<<" [ns]\n"; // skruninf_.softtrg_t0_offset
		// n.b. COUNT_PER_NSEC = 1.92(count/nsec), #defined constant in skheadC.h
	}
	
	//UInt_t dummy[512-1-5*32]; // what does this hold? anything of interest?
	
	return true;
}

bool PrintEvent::PrintBadChannels(){
	
	std::cout<<"\nPrinting Bad Channels\n"<<std::endl;
	
	// TODO
	/*     from 'skbadc.h' or skbadc0C.h */
	/*     COMMON block COMBAD: infomation of bad channels*/
	/*    ===== inner detector =====*/
	/*       NBAD           : number of badch*/
	/*       IBAD(1:MAXPM)  : ch is bad if IBAD(ch).ne.0*/
	/*       ISQBAD(1:NBAD) : cable number of badch*/
	/*       NBAD0           : number of badch              imaskbad.eq.0*/
	/*       IBAD0(1:MAXPM)  : ch is bad if IBAD(ch).ne.0   imaskbad.eq.0*/
	/*       ISQBAD0(1:NBAD) : cable number of badch        imaskbad.eq.0*/
	/*    ===== outer detector =====*/
	/*       NBADA            : number of badch*/
	/*       IBADA(1:MAXPMA)  : ch is bad if IBADA(ch).ne.0*/
	/*       ISQBADA(1:NBADA) : cable number of badch*/
	/*    ===== common option =====*/
	/*       IMASKBAD       : kind of bad channels for current common data*/
	/*                             0: all kinds of bad channels are selected (normal)*/
	/*                            -1: mask only "badch.00*" and 'criterion 5'*/
	/*                          2**0: mask "badch.dat"*/
	/*                          2**1: mask criterion 1  inner dead 1*/
	/*                          2**2: mask criterion 2  inner dead 2*/
	/*                          2**3: mask criterion 3  inner noisy*/
	/*                          2**4: mask criterion 4  outer */
	/*                          2**5: mask criterion 5  inner HK PMT*/
	/*       IMASKBADOPT    : kind of bad channels for next skbadch() action*/
	/*       log_level_skbadch : log level of findconst() for skbadch*/
	/*                 1) Lots of output*/
	/*                 2) only prints filenames*/
	/*                 3) only prints when not found*/
	/*                 4) do not print*/
	/*                 others) only prints when found*/
	/*
	combad_: ID bads common block
		int    nbad;
		int    ibad[MAXPM];
		int    isqbad[MAXNBD];  // MAXNBD = #defined 11146
		int    imaskbad;
		int    imaskbadopt;
		int    log_level_skbadch;
	combada_: OD bads common block
		int    nbada;
		int    ibada[MAXPMA];
		int    isqbada[MAXNBDA]; // MAXNBDA = #defined 1885
	combad0_:  how is this different to combad_ ??
		int    nbad0;              // number of badch
		int    ibad0[1:MAXPM];     // ch is bad if IBAD(ch).ne.0
		int    isqbad0[1:NBAD];    // cable number of badch
		int    imaskbad0;          // kind of bad channels for current common data:
		                              0: all kinds of bad channels are selected (normal)
		                              "force to be ZERO" ???
		int    imaskbadopt0;       // kind of bad channels for next skbadch0() action
	combad00_: ID bads not masked?
		int    nbad0;
		int    ibad0[MAXPM];
		int    isqbad0[MAXNBD];
	*/
	
	return true;
}

bool PrintEvent::PrintDarkInfo(){
	
	std::cout<<"\nPrinting Dark Info\n"<<std::endl;
	// TODO
	/*     from 'skbadc.h'*/
	/*     COMMON block COMDARK: infomation of dark rate (filled by src/skrd/skdark.F)*/
	/*    ===== inner detector =====*/
	/*       nrun_dark          : run number of the current darkrate()*/
	/*       nsub_dark          : # of used subruns to make the table*/
	/*       nfill_dark         : # of used events to make the table*/
	/*       dark_ave           : averaged dark rate of the run*/
	/*       dark_rate(1:MAXPM) : dark rate of each channel of the run*/
	/*       log_level_skdark : log level of findconst() for skdark*/
	/*                 1) Lots of output*/
	/*                 2) only prints filenames*/
	/*                 3) only prints when not found*/
	/*                 4) do not print*/
	/*                 others) only prints when found*/
	/*
	comdark_: dark common
	  int    nrun_dark;
	  int    nsub_dark;
	  int    nfill_dark;
	  float  dark_ave;
	  float  dark_rate[MAXPM];
	  float  dark_rate_od[MAXPMA];
	  float  dark_rate_od_subped[MAXPMA];
	  int    log_level_skdark;
	*/
	return true;
}


bool PrintEvent::Finalise(){

  return true;
}
