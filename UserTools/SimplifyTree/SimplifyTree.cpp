/* vim:set noexpandtab tabstop=4 wrap */
#include "SimplifyTree.h"

#include "ConnectionTable.h"
#include "fortran_routines.h"
#include <bitset>
#include <unistd.h>

SimplifyTree::SimplifyTree():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}

bool SimplifyTree::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);            // how verbose to be
	std::string treeReaderName;
	m_variables.Get("treeReaderName",treeReaderName);
	m_variables.Get("OutputDir",outputdir);
	
	// check for input TreeReader
	 if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		// retrieve the TreeReader
		iTreeReader = m_data->Trees.at(treeReaderName);
		
		// only used by debug code
		LUN = m_data->GetLUN(treeReaderName);
		
		// old method: use TreeReader Tool in root2root mode, accessing hits from the TQReal branch
		// of the output file as it is populated.Produces a redundant output file, so we no longer use it.
		/*
		// get the LUN associated with it
		LUN = m_data->GetLUN(treeReaderName);
		
		// use this to get the corresponding TreeManager
		TreeManager* mgr = skroot_get_mgr(&LUN);
		
		// and use the TreeManager to get access to the output TTree as it is being written
		// (saves us converting the file and then having to run a second toolchain to read it)
		TTree* otree = mgr->GetOTree();
		
		// note that this is equivalent to reading the RAWTQINFO common block (being used to populate it),
		// so accesses hit information without bad channel masking being applied by skread.
		// it may therefore give different hit info than if the TQREAL branch of the output file
		// was read back in with bad channel masking being applied.
		get_ok = oTreeReader.Load(otree);
		*/
	}
	
	// class to get PMT positions from cable number
	myConnectionTable = m_data->GetConnectionTable();
	
	return true;
}


bool SimplifyTree::Execute(){
	
	Log(toolName+": Executing",v_debug,verbosity);
	if(iexecute%1000==0) Log(toolName+": "+std::to_string(iexecute),v_message,verbosity);
	++iexecute;
	
	// check if processing a new file
	if(infilename!=iTreeReader->GetFile()->GetName()){
		
		// new input file -> new output file
		Log(toolName+": New input file, new output file",v_debug,verbosity);
		
		// write out and close the current output file, if there is one
		if(fout){
			fout->Write("",TObject::kOverwrite);
			outtree->ResetBranchAddresses();
			fout->Close();
			delete fout; fout=nullptr;
		}
		
		// make an output file named after the input file
		infilename = iTreeReader->GetFile()->GetName();
		// strip path, put in requested (default=current) directory
		std::string outfname = basename(infilename.c_str());
		outfname = outfname.substr(0, outfname.size()-5);
		outfname += "_simple.root";
		if(outputdir!="") outfname=outputdir+"/"+outfname;
		Log(toolName+": New output simple file: "+outfname,v_debug,verbosity);
		fout = new TFile(outfname.c_str(), "recreate");
		if(!fout || fout->IsZombie()){
			Log(toolName+": Error making output file "+outfname,v_error,verbosity);
			return false;
		}
		fout->cd();
		
		// make the output TTree
		outtree = new TTree("data","SK Hits");
		outtree->Branch("sk_phase",&sk_phase);
		outtree->Branch("run_num",&run_num);
		outtree->Branch("subrun_num",&subrun_num);
		outtree->Branch("event_num",&event_num);
		outtree->Branch("event_timestamp",&timestring);
		outtree->Branch("trigger_flags",&trigger_flags);
		outtree->Branch("event_flags",&event_flags);
		outtree->Branch("readout_t0",&readout_t0);
		outtree->Branch("trigger_t0",&trigger_t0);
		outtree->Branch("readout_len",&readout_len);
		outtree->Branch("hit_id",&hit_id);
		outtree->Branch("hit_q",&hit_q);
		outtree->Branch("hit_t",&hit_t);
		outtree->Branch("hit_x",&hit_x);
		outtree->Branch("hit_y",&hit_y);
		outtree->Branch("hit_z",&hit_z);
		outtree->Branch("hit_theta",&hit_theta);
		outtree->Branch("hit_loc",&hit_loc);
		outtree->Branch("hit_flags",&hit_flags);
	}
	
	// clear the hit vectors
	hit_id.clear();       // PMT number
	hit_t.clear();        // hit time
	hit_q.clear();        // hit charge
	hit_x.clear();        // PMT X
	hit_y.clear();        // PMT Y
	hit_z.clear();        // PMT Z
	hit_theta.clear();    // PMT azimuthal angle in barrel
	hit_loc.clear();      // PMT location identifier: ID / OD / top cap/ bottom cap / barrel
	hit_flags.clear();   // was hit in tight timing window, and other flags
	
	// basic event info
	sk_phase = skheadg_.sk_geometry;
	run_num = skhead_.nrunsk;
	subrun_num = skhead_.nsubsk;
	event_num = skhead_.nevsk;
	trigger_flags = skhead_.idtgsk;
	event_flags = skhead_.ifevsk;
	readout_t0 = skheadqb_.it0sk;
	trigger_t0 = skheadqb_.it0xsk;
	readout_len = skheadqb_.gatewsk;
	
	// event time and date
	tm rundate = {0};
	rundate.tm_year = skhead_.ndaysk[0];
	rundate.tm_mon = skhead_.ndaysk[1] - 1;
	rundate.tm_mday = skhead_.ndaysk[2];
	rundate.tm_hour = skhead_.ntimsk[0];
	rundate.tm_min = skhead_.ntimsk[1];
	rundate.tm_sec = skhead_.ntimsk[2];
	time_t runtime = mktime(&rundate);         // need to use mktime to derive day of week
	timestring = ctime(&runtime);              // format timestamp into a string
	timestring.pop_back();                     // drop trailing newline
	//std::cout<<"The event occurred at "<<timestring<<std::endl;
	
	/*
	// for files with TQReal populated, we could read that branch
	// while rfm files do not have TQREAL populated, if we run the TreeReader in root2root mode,
	// we can use the following to populate the TQREAL branch of an output file
	skroot_set_tree_(&LUN);
	
	// we could then get the TQREAL branch from the output file, to access the converted data
	get_ok = GetBranchValues();
	*/
	
	int thisevttrigflg=-99; // debug variable
	
	// get event time and subevent time
	it0sk = skheadqb_.it0sk;
	it0xsk = skheadqb_.it0xsk; // same as it0sk for primary trigger (first in readout)
	// there may be other triggers (e.g. SLE NHITS trigger) captured within the readout
	// window associated by this primary trigger (e.g. the SLE NHITS condition is met at some
	// later time within the 40us LE primary trigger readout window)
	// we do not extract such "subtriggers", which would have different IT0XSK values
	
	// loop twice; once over ID, once over OD hits
	for(int id_od=0; id_od<2; ++id_od){
		
		std::string id_or_od = (id_od) ? "ID" : "OD";
		
		int n_hits;
		if(id_or_od=="ID"){
			// ID
			//n_hits = myTQReal->cables.size();
			n_hits = rawtqinfo_.nqisk_raw;
		} else {
			// OD
			//n_hits = myTQAReal->cables.size();
			n_hits = rawtqinfo_.nhitaz_raw;
		}
		
		Log(toolName+": had "+std::to_string(n_hits)+" "+id_or_od+" hits",v_debug,verbosity);
		
		int n_ingate_hits=0;
		
		// loop over hits
		for (int ihit = 0; ihit < n_hits; ++ihit){
			
			// get hit details
			
			/*
			// from TQREAL/TQAREAL object
			TQReal* tq_real = (id_or_od=="ID") ? myTQReal : myTQAREAL;
			int cableNumber = tq_real->cables.at(ihit);
			float charge = tq_real->Q.at(ihit);
			float time = tq_real->T.at(ihit);
			long it0xsk = tq_real->it0xsk;
			long it0sk = myHeader->t0;
			*/
			
			// straight from common blocks used to fill TQREAL object
			
			if(id_or_od=="ID"){
				cableNumber = rawtqinfo_.icabbf_raw[ihit];
				charge = rawtqinfo_.qbuf_raw[ihit];
				time = rawtqinfo_.tbuf_raw[ihit];
			} else {
				cableNumber = rawtqinfo_.icabaz_raw[ihit];
				charge = rawtqinfo_.qaskz_raw[ihit];
				time = rawtqinfo_.taskz_raw[ihit];
			}
			
			// the upper 16 bits are a set of hit flags
			// do this first as the shift does not modify cableNumber
			int hitflags = cableNumber >> 16;
			// mask off the lower 16 bits for the PMT number
			cableNumber = cableNumber & 0x0000FFFF;
			
			std::bitset<8*sizeof(int)> hitflagsbits(hitflags);
			if(hitflagsbits.test(1)){
				//std::cout<<"in gate flag from icabbf_raw for hit "<<ihit<<std::endl;
				++n_ingate_hits;
			}
			
			/*
			if(hitflags!=0) std::cout<<""<<hitflags<<std::endl;
			int hitflags2 = sktqz_.ihtiflz[ihit];
			std::bitset<8*sizeof(int)> hitflags2bits(hitflags2);
			if(hitflags2bits.test(1)) std::cout<<"in gate flag from ihtiflz for hit "<<ihit<<std::endl;
			int hitflags3 = sktqz_.iqiskz[ihit];
			std::bitset<8*sizeof(int)> hitflags3bits(hitflags3);
			if(hitflags3bits.test(11)) std::cout<<"in gate flag from iqiskz for hit "<<ihit<<std::endl;
			*/
			
			/*
			for reference:
			 
			 from $SKOFL_ROOT/const/connection.super.sk-4.dat
				# -- for inner-PMTs (1-11146)
				# -- for muon VETO(11151,11152,11153,11154)
				# -- for calibration ID (11155-?, see skveto.h for details)
				# -- for muon chamber(only hut3 and hut4)
				# -- for trigger ID QB (15001-15240, see skhead.h for details)
				# -- for anti-PMT(20001-21885)
			
			 from $SKOFL_ROOT/inc/sktq.h, meaning of bits of hitflags (IHTIFLZ):
				#    11-6   (# of TRG EVENT COUNTER - 1) * 64 (0-63)
				#    5-4    charge range (0:Small, 1:Medium, 2:Large)
				#    3-2    trig ID (0: Narrow, 1: Wide, 2: Pedestal, 3: Not used)
				#    1:     In gate (1=in gate, 0=not in gate)
				#    0:     In 1.3usec (1=in, 0=out)
				
			 from $SKOFL_ROOT/skhead.h, meaning of bits in trigger_flags (IDTGSK):
				 bit   Mask         Trigger (SK-IV)                       Trigger (SKI-III)
				 ------------------------------------------------------------------------
				  0  0x00000001     LE  (by software trig.)               LE
				  1  0x00000002     HE  (by software trig.)               HE
				  2  0x00000004     SLE (by software trig.)               SLE
				  3  0x00000008     OD  (by software trig.)               OD (normal run) or
				                                                          Fission trig (Ni run)
				  4  0x00000010     (same as SK-I,II,III)                 Periodical firing of:
				                                                          a. nothing (null trigger)
				                                                          b. TQ map laser
				                                                          c. water attenuation meas. laser
				                                                          d. Xe ball
				  5  0x00000020     (same as SK-I,II,III)                 After trigger (normal run) or
				                                                          Periodical firing of (Calib run):
				                                                          a. Laser trigger
				                                                          b. Xe trigger
				                                                          c. Ni trigger
				                                                          d. Linac trigger
				  6  0x00000040     (same as SK-I,II,III)                 VETO START
				  7  0x00000080     (same as SK-I,II,III)                 VETO STOP
				  8  0x00000100     N/A                                   N/A
				  9  0x00000200     N/A                                   N/A
				 10  0x00000400     N/A                                   N/A
				 11  0x00000800     Random Wide Trigger                   N/A
				 12  0x00001000     Laser (ID, Usho Laser)                N/A
				 13  0x00002000     LED                                   N/A
				 14  0x00004000     Ni                                    N/A
				 15  0x00008000     Laser (OD, AutoTQlaser)               N/A
				 16  0x00010000     LE  (by hitsum signal)                N/A
				 17  0x00020000     HE  (by hitsum signal)                N/A
				 18  0x00040000     SLE (by hitsum signal)                N/A
				 19  0x00080000     OD  (by hitsum signal)                N/A
				 20  0x00100000     N/A                                   N/A
				 21  0x00200000     N/A                                   N/A
				 22  0x00400000     SN Burst                              N/A
				 23  0x00800000     mu->e decay                           N/A
				 24  0x01000000     LINAC                                 N/A
				 25  0x02000000     LINAC microwave                       N/A
				 26  0x04000000     N/A                                   N/A
				 27  0x08000000     Periodic (simple periodic trigger)    N/A
				 28  0x10000000     SHE trigger    (by software trig.)    N/A
				 29  0x20000000     After trigger  (by software trig.)    N/A
				 30  0x40000000     pedestal event (by software trig.)    N/A
				 31  0x80000000     T2K event      (by software trig.)    N/A
				
			 from $SKOFL_ROOT/skhead.h, meaning of bits in event_flags (IFEVSK):
				 bit   Mask     SK-I,II,III                               SK-IV
				 ----------------------------------------------------------------------
				 0  0x00000001  ATM                                       QBEE TQ
				 1  0x00000002  TRG                                       HARD TRG
				 2  0x00000004  SMP REGISTER                              QBEE STAT
				 3  0x00000008  SCALER                                    DB_STAT_BLOCK
				 4  0x00000010  PEDESTAL START                            CORRUPTED_CHECKSUM
				 5  0x00000020  PEDESTAL DATA(ATM)                        MISSING SPACER
				 6  0x00000040  PEDESTAL HISTOGRAM                        PED_HIST_BLOCK
				 7  0x00000080  PEDESTAL END                              N/A
				 8  0x00000100  END OF RUN                                N/A
				 9  0x00000200  PEDESTAL(ON)                              PEDESTAL ON
				 0  0x00000400                                            RAW_AMT_BLOCK
				 1  0x00000800  GPS DATA                                  GPS DATA
				 2  0x00001000  CAMAC ADC                                 PEDESTAL_CHECK
				 3  0x00002000  ANTI DATA                                 SEND_BLOCK
				 4  0x00004000  INNER SLOW DATA                           INNER SLOW DATA
				 5  0x00008000  RUN INFORMATION                           RUN INFORMATION
				 6  0x00010000  ERROR (TKO-PS)                            PREV T0 BLOCK
				 7  0x00020000  ERROR (HV-PS)                             N/A
				 8  0x00040000  ERROR (TEMPERARTURE)                      FE_TRL_BLOCK
				 9  0x00080000                                            SPACER_BLOCK
				 0  0x00100000  UNCOMPLETED ATM DATA                      INCOMPLETE TQ
				 1  0x00200000  INVALID     ATM DATA                      CORRUPT TQ BLOCK
				 2  0x00400000                                            TRG MISMATCH TQ
				 3  0x00800000                                            QBEE ERROR
				 4  0x01000000  ERROR (DATA)                              SORT_BLOCK
				 5  0x02000000  UNREASONABLE DATA                         CORRUPTED_BLOCK
				 6  0x04000000  LED BURST ON                              LED BURST ON
				 7  0x08000000                                            EVNT TRAILER
				 8  0x10000000  INNER DETECTOR OFF                        INNER DETECTOR OFF
				 9  0x20000000  ANTI  DETECTOR OFF                        ANTI  DETECTOR OFF
				 0  0x40000000                                            T2K GPS
				 1  0x80000000  TRG IS AVAILABLE                          (EVNT HDR)&(SOFTWARE TRG)
			
			*/
			
			// skip channels that aren't relevant PMTs
			if( id_or_od=="ID" && (cableNumber==0 || cableNumber>MAXPM) ) continue;
			if( id_or_od=="OD" && (cableNumber<20001 || cableNumber>(20000+MAXPMA)) ) continue;
			
			// use cable number to get PMT position
			myConnectionTable->GetTubePosition(cableNumber, tubePosition);
			
			// get location (barrel, top/bottom cap, ID/OD...)
			// enum Locations{ kIDTop, kIDWall, kIDBot, kODTop, kODWall, kODBot }; (from ConnectionTable.cc)
			int loc = myConnectionTable->GetLocation(tubePosition[0],tubePosition[1],tubePosition[2]);
			
			// calculate angle in barrel
			double tubeR = sqrt(pow(tubePosition[0], 2.f) + pow(tubePosition[1],2.f));
			if(tubeR>0){
				tubetheta = acos(tubePosition[0] / tubeR);
			} else {
				tubetheta = 0;
			}
			if(tubePosition[1] > 0) tubetheta = -tubetheta;
			
			// convert hit time within readout window to hit time within subtrigger window
			// by subtracting the time from 40us buffer readout start to subtrigger start
			// then convert TDC ticks to nanoseconds
			// COUNT_PER_NSEC is a #defined constant in '$SKOFL_ROOT/inc/skheadC.h'
			// describing conversion from TDC ticks to nanoseconds
			time = time -(it0xsk-it0sk)/COUNT_PER_NSEC;
			
			// transfer data to output variables
			hit_id.push_back(cableNumber);
			hit_q.push_back(charge);
			hit_t.push_back(time);
			hit_flags.push_back(hitflags);
			hit_x.push_back(tubePosition[0]);
			hit_y.push_back(tubePosition[1]);
			hit_z.push_back(tubePosition[2]);
			hit_theta.push_back(tubetheta);
			hit_loc.push_back(loc);
			
			// consistency check - IHTIFLZ bits 2-3 are "trigger ID", so should be the same for all hits, right?
			// seems to pass this test, so we may only need to check the first one if we're interested.
			// in fact, it seems to always be 0.
			/*
			if(id_or_od=="ID"){
				int trigid = (hitflags >> 2) & 0x00000003;
				if(trigid>3){
					std::cerr<<"trigid>3; hitflags is: "
					         <<std::bitset<8*sizeof(int)>(hitflags).to_string()<<std::endl;
				}
				if(ihit==0){
					thisevttrigflg = trigid;
				} else {
					if(thisevttrigflg!=trigid){
						std::cout<<"inconsistent trig flag in hit "<<ihit<<" in event "<<event_num<<std::endl;
						std::cout<<"previous trigger ID was "<<thisevttrigflg<<", this hit has trigger ID "
							     <<trigid<<std::endl;
						thisevttrigflg=trigid;
					}
				}
			}
			*/
			
		}  // end loop over hits
		
		/*
		if(id_or_od=="ID"){
			std::cout<<"event had "<<n_ingate_hits<<" hits in gate"<<std::endl;
			std::bitset<8*sizeof(int)> trigger_bits(trigger_flags);
			if(trigger_bits.test(28)){ std::cout<<"SHE trigger event "<<event_num<<std::endl; }
			else if(trigger_bits.test(1)){ std::cout<<"HE trigger event "<<event_num<<std::endl; }
			else if(trigger_bits.test(0)){ std::cout<<"LE trigger event "<<event_num<<std::endl; }
		}
		*/
		
		if(verbosity>3 && id_or_od=="ID"){
			// some trigger flag checks just for info
			// hardware NHITS triggers do not seem to be set, at least for rfm_run080008.000008.root
			std::bitset<8*sizeof(int)> trigger_bits(trigger_flags);
			std::cout<<"SLE HWTRG="<<trigger_bits.test(18)<<", SLE SWTRG="<<trigger_bits.test(2)<<"\n"
				     <<"LE HWTRG="<<trigger_bits.test(16)<<", LE SWTRG="<<trigger_bits.test(0)<<"\n"
				     <<"HE HWTRG="<<trigger_bits.test(17)<<", HE SWTRG="<<trigger_bits.test(1)<<"\n"
				     <<"SHE SWTRG="<<trigger_bits.test(28)<<"\n"
				     <<"AFT HWTRG="<<trigger_bits.test(5)<<", AFT SWTRG="<<trigger_bits.test(29)<<"\n"
				     <<"PED SWTRG="<<trigger_bits.test(30)<<"\n";
			// some event flag checks to try to make sense of them: seems like the following are always 1
			// at least in events normally returned by the TreeReader Tool
			std::bitset<8*sizeof(int)> event_bits(event_flags);
			std::cout<<"qbee_tq="<<event_bits.test(0)<<", hw_trg="<<event_bits.test(1)<<"\n";
			// some hit flag checks to try to make sense of them: these seem to be always 0,
			// at least in events normally returned by the TreeReader Tool
			std::cout<<"first hit IHTIFLZ trigger ID: ";
			if(hit_t.size()) std::cout<<thisevttrigflg<<"\n";
			else std::cout<<"N/A (no hits)\n\n";
		}
		
	}  // end loop over ID/OD
	
	/*
	// a readout may contain multiple NHITS software triggers...
	// do we need to scan over these to split this up? maybe...!
	int num_sw_trgs = skheadqb_.ntrigsk;
	std::cout<<"main trigger: it0sk is "<<skheadqb_.it0sk<<", it0xsk is "<<skheadqb_.it0xsk<<std::endl;
	
	// note - not all trigger types supported, just the NHITS ones, PERIODIC and PEDESTAL.
	std::vector<std::pair<std::string,int>> trgtypes{{"SLE",2},{"LE",0},{"HE",1},{"SHE",28}};
	for(auto&& trgtype : trgtypes){
		int trgtypenum = trgtype.second;
		int ntrgsoftype=-99;
		std::vector<int> subtrg_t0s(10);
		int max_trgs = 10;
		// check how many we had and get their respective t0s
		get_sub_triggers_(&trgtypenum, &ntrgsoftype, subtrg_t0s.data(), &max_trgs);
		std::cout<<"we had "<<ntrgsoftype<<" triggers of type "<<trgtype.first<<" at:\n";
		for(int x=0; x<ntrgsoftype; ++x){
			std::cout<<"\t"<<subtrg_t0s.at(x)<<"\n";
		}
		// shift the timing gate and count the number of hits in each subtrigger
		for(int x=0; x<ntrgsoftype; ++x){
			int newt0 = readout_t0 + subtrg_t0s.at(x);
			set_timing_gate_(&newt0);
			// newt0 and it0xsk are the same, it0sk is constant.
			std::cout<<"subtrigger "<<x<<": it0sk is "<<skheadqb_.it0sk
			         <<", newt0 is "<<newt0<<", it0xsk is "<<skheadqb_.it0xsk<<std::endl;
			// recalculate hit info (includes hit flags ihtiflz i think)
			int LUN2=-LUN;
			skcread_(&LUN2, &get_ok);
			// nqisk_raw is always the same but nqisk does change.
			std::cout<<"nqisk in this subtrigger is "<<skq_.nqisk
			         <<", nqisk_raw is "<<rawtqinfo_.nqisk_raw<<std::endl;
			//delete_outside_hits_();
			// calling delete_outside_hits_ results in nqisk_raw matching nqisk for the first subtrigger,
			// but subsequently all other subtriggers have both nqisk and nqisk_raw 0!
			// seems like skcread does not re-retrieve the deleted hits?!?
		}
	}
	*/
	
	outtree->Fill();
	
	// invoke TTree::Fill on output file in root2root mode. We don't need the resulting file, but why not?
	// XXX XXX BEWARE: THIS CLEARS TQREAL!!! (and what else?) ONLY DO IT UNTIL AFTER WE'RE DONE WITH THEM!
	// skroot_fill_tree_(&LUN);
	
	return true;
}


bool SimplifyTree::Finalise(){
	
	// write last file, if there is one
	if(fout){
		fout->Write("",TObject::kOverwrite);
		outtree->ResetBranchAddresses();
		fout->Close();
		delete fout;
	}
	
	// tell the oTreeReader not to try to delete the file it's reading,
	// that will already be done by the TreeManager
	//oTreeReader.SetClosed();
	
	return true;
}

bool SimplifyTree::GetBranchValues(){
	get_ok = true;
	get_ok &= oTreeReader.Get("TQREAL", myTQReal);
	get_ok &= oTreeReader.Get("TQAREAL", myTQAReal);
	get_ok &= oTreeReader.Get("HEADER", myHeader);
	return get_ok;
}
