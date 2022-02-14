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
	
	// check for input TreeReader
	 if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		// retrieve the TreeReader
		iTreeReader = m_data->Trees.at(treeReaderName);
		
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
		
		// write out and close the current output file, if there is one
		if(fout){
			fout->Write("",TObject::kOverwrite);
			outtree->ResetBranchAddresses();
			fout->Close();
			delete fout; fout=nullptr;
		}
		
		// make an output file named after the input file
		infilename = iTreeReader->GetFile()->GetName();
		// strip path, put in current directory, not alongside input file
		std::string outfname = basename(infilename.c_str());
		outfname = outfname.substr(0, outfname.size()-5);
		outfname += "_simple.root";
		Log(toolName+": New output simple file: "+outfname,v_debug,verbosity);
		fout = new TFile(outfname.c_str(), "recreate");
		
		// make the output TTree
		outtree = new TTree("data","SK Hits");
		outtree->Branch("hit_id",&hit_id);
		outtree->Branch("hit_q",&hit_q);
		outtree->Branch("hit_t",&hit_t);
		outtree->Branch("hit_x",&hit_x);
		outtree->Branch("hit_y",&hit_y);
		outtree->Branch("hit_z",&hit_z);
		outtree->Branch("hit_theta",&hit_theta);
		outtree->Branch("hit_loc",&hit_loc);
		outtree->Branch("hit_ingate",&hit_ingate);
	}
	
	/*
	// for files with TQReal populated, we could read that branch
	// this would populate TQReal branch for raw data input files, if TreeReader is in root2root mode.
	skroot_set_tree_(&LUN);
	
	// get next TQREAL branch, which accesses the populated data
	get_ok = GetBranchValues();
	*/
	
	// clear the hit vectors
	hit_id.clear();       // PMT number
	hit_t.clear();        // hit time
	hit_q.clear();        // hit charge
	hit_x.clear();        // PMT X
	hit_y.clear();        // PMT Y
	hit_z.clear();        // PMT Z
	hit_theta.clear();    // PMT azimuthal angle in barrel
	hit_loc.clear();      // PMT location identifier: ID / OD / top cap/ bottom cap / barrel
	hit_ingate.clear();   // was hit in tight timing window, and other flags
	
	// loop twice; once over ID, once over ID
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
				it0sk = skheadqb_.it0sk;
				it0xsk = skheadqb_.it0xsk;
			} else {
				cableNumber = rawtqinfo_.icabaz_raw[ihit];
				charge = rawtqinfo_.qaskz_raw[ihit];
				time = rawtqinfo_.taskz_raw[ihit];
				it0sk = skheadqb_.it0sk;
				it0xsk = skheadqb_.it0xsk;
			}
			
			// we only need the lower 16 bits of the PMT number
			int in_gate = cableNumber >> 16;
			cableNumber = cableNumber & 0x0000FFFF;
			
			/*
			 for reference:
			 from $SKOFL_ROOT/const/connection.super.sk-4.dat
				# -- for inner-PMTs (1-11146)
				# -- for muon VETO(11151,11152,11153,11154)
				# -- for calibration ID (11155-?, see skveto.h for details)
				# -- for muon chamber(only hut3 and hut4)
				# -- for trigger ID QB (15001-15240, see skhead.h for details)
				# -- for anti-PMT(20001-21885)
			 from $SKOFL_ROOT/inc/sktq.h, meaning of bits of in_gate (IHTIFLZ):
				#    11-6   (# of TRG EVENT COUNTER - 1) * 64 (0-63)
				#    5-4    charge range (0:Small, 1:Medium, 2:Large)
				#    3-2    trig ID (0: Narrow, 1: Wide, 2: Pedestal, 3: Not used)
				#    1:     In gate (1=in gate, 0=not in gate)
				#    0:     In 1.3usec (1=in, 0=out)
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
			if(loc==1 || loc==4){
				double tubeR = sqrt(pow(tubePosition[0], 2.f) + pow(tubePosition[1],2.f));
				tubetheta = acos(tubePosition[0] / tubeR);
				if(tubePosition[1] > 0) tubetheta = -tubetheta;
			} else {
				tubetheta=0;
			}
			
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
			hit_ingate.push_back(in_gate);
			hit_x.push_back(tubePosition[0]);
			hit_y.push_back(tubePosition[1]);
			hit_z.push_back(tubePosition[2]);
			hit_theta.push_back(tubetheta);
			hit_loc.push_back(loc);
			
		}
		
	}
	
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
