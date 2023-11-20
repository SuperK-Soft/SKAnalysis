/*
* Outputs to LOWE branch but the MC variables may differ from the fortran code.
* Assuming IBD:
* posmc[3]      - mc true position of the positron
* pabsmc[2]     - mc absolute momentum of the positron and neutron
* energymc[2]   - mc energy of the positron and neutron
* dirmc[2][3]   - mc true direction of the positron and neutron
*/

#include "CombinedFitter.h"

#include "Algorithms.h"
#include "Constants.h"

#include <string>
#include <iostream>
#include <bitset>

// including bonsai in skofl rather than putting it all in the DataModel
#include "searchgrid.h"
#include "bscalls.h"
#include "pairlikelihood.h"
#include "fourhitgrid.h"
#include "combinedgrid.h"
#include "goodness.h"

#include <TCanvas.h>
#include <TRandom.h>

CombinedFitter::CombinedFitter():Tool(){}


bool CombinedFitter::Initialise(std::string configfile, DataModel &data){
	
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",m_verbose);    // how verbose to be
	m_variables.Get("readerName",readerName);  // name given to the TreeReader used for file handling
	m_variables.Get("dataSrc",dataSrc);        // where to get the data from (common blocks/tqreal)
	m_variables.Get("bonsaiSrc",bonsaiSrc);    // which bonsai to use (skofl or local)
	m_variables.Get("addNoise",addNoise);      // whether to add noise to AFT
	m_variables.Get("outputFile",outputFile);  // name of file to save ntuples to
	
	// use the readerName to find the LUN associated with this file
	std::map<std::string,int> lunlist;
	m_data->CStore.Get("LUNList",lunlist);
	if(lunlist.count(readerName)==0){
		Log(m_unique_name+" error! No LUN associated with readerName "+readerName,v_error,m_verbose);
		return false;
	}
	lun = lunlist.at(readerName);
	
	// determine if this file is MC or not, so we know whether or not to get
	// the MC information
	t = skroot_get_tree(&lun);
	if(t) MC = (t->FindBranch("MC")!=nullptr);
	
	// make output file and tree
	if(outputFile == nullptr) outputFile = "CombinedFitter_output.root";
	ifstream f(outputFile.c_str());
	if (f.good()) fout = new TFile(outputFile.c_str(), "UPDATE");
	else fout = new TFile(outputFile.c_str(),"CREATE");
	outputTree = new TTree("corefit", "Coincidence reconstruction results");
	
	// Set the branches we need to fill in the ntuple
	outputTree->Branch("x",&x,"x/F");
	outputTree->Branch("y",&y,"y/F");
	outputTree->Branch("z",&z,"z/F");
	outputTree->Branch("x_prev",&x_prev,"x_prev/F");
	outputTree->Branch("y_prev",&y_prev,"y_prev/F");
	outputTree->Branch("z_prev",&z_prev,"z_prev/F");
	outputTree->Branch("x_combined",&x_combined,"x_combined/F");
	outputTree->Branch("y_combined",&y_combined,"y_combined/F");
	outputTree->Branch("z_combined",&z_combined,"z_combined/F");
	outputTree->Branch("mcx",&mcx,"mcx/F");
	outputTree->Branch("mcy",&mcy,"mcy/F");
	outputTree->Branch("mcz",&mcz,"mcz/F");
	outputTree->Branch("mcx_prev",&mcx_prev,"mcx_prev/F");
	outputTree->Branch("mcy_prev",&mcy_prev,"mcy_prev/F");
	outputTree->Branch("mcz_prev",&mcz_prev,"mcz_prev/F");
	outputTree->Branch("mcx_ncapture",&mcx_ncapture,"mcx_ncapture/F");
	outputTree->Branch("mcy_ncapture",&mcy_ncapture,"mcy_ncapture/F");
	outputTree->Branch("mcz_ncapture",&mcz_ncapture,"mcz_ncapture/F");
	outputTree->Branch("mct_ncapture",&mct_ncapture,"mct_ncapture/F");
	outputTree->Branch("pdg",&pdg,"pdg/I");
	outputTree->Branch("pdg_prev",&pdg_prev,"pdg_prev/I");
	outputTree->Branch("mc_energy",&mc_energy,"mc_energy/F");
	outputTree->Branch("mc_energy_prev",&mc_energy_prev,"mc_energy_prev/F");
	outputTree->Branch("nhitsAFT_raw",&nhitsAFT_raw,"nhitsAFT_raw/I");
	outputTree->Branch("tgood",&tgood,"tgood/F");
	outputTree->Branch("tgood_prev",&tgood_prev,"tgood_prev/F");
	outputTree->Branch("tgood_combined",&tgood_combined,"tgood_combined/F");
	outputTree->Branch("tgood_combined_prev",&tgood_combined_prev,"tgood_combined_prev/F");
	outputTree->Branch("bsn50",&bsn50,"bsn50/I");
	outputTree->Branch("bsn10",&bsn10,"bsn10/I");
	outputTree->Branch("bsn50AFT",&bsn50AFT,"bsn50AFT/I");
	outputTree->Branch("bsn10AFT",&bsn10AFT,"bsn10AFT/I");
	outputTree->Branch("n50_prev",&n50_prev,"n50_prev/I");
	outputTree->Branch("n10_prev",&n10_prev,"n10_prev/I");
	outputTree->Branch("n50_aft",&n50_aft,"n50_aft/I");
	outputTree->Branch("n10_aft",&n10_aft,"n10_aft/I");
	
	// initialize water transparency table
	// (this will be for energy reconstruction I presume)
	skrunday_();
	skwt_gain_corr_();
	
	// get detector parameters from common blocks (geopmtC.h)
	numPMTs = MAXPM; //number of PMTs in the detector
	float* xyzpm = &geopmt_.xyzpm[0][0];//lf_allfit_new.cc
	//nrunsk_last = skhead_.nrunsk;
	//nsubsk_last = skhead_.nsubsk;
	
	// Initialize BONSAI - create BONSAI objects from the PMT position array
	//------------------
	
	// Initialise bonsai for coincidence reconstruction
	bsgeom = new pmt_geometry(numPMTs,xyzpm);
	bspairlike = new pairlikelihood(bsgeom->cylinder_radius(),bsgeom->cylinder_height());
	bspairfit = new bonsaifit(bspairlike);
	
	// Use built-in BONSAI functions (bonsai/bscalls.cc::)
	m_data->BonsaiInit();
	
	// check where we are getting the hit info from
	if (dataSrc==0) 
		Log(m_unique_name+": Getting hit info from common blocks skt_/skq_",v_message,m_verbose);
	else if (dataSrc==1) {
		Log(m_unique_name+": Getting raw hit info from common block sktqz_",v_message,m_verbose);
		Log(m_unique_name+": Warning: this method is not yet functional (need to remove bad channels, etc)",v_warning,m_verbose);
	}
	else {
		Log(m_unique_name+": Getting raw hit info from tqreal branch",v_message,m_verbose);
		Log(m_unique_name+": Warning: this method may not give exactly the same result as skt_/skq_ common blocks",v_warning,m_verbose);
	}
	
	return true;
}

bool CombinedFitter::Execute(){
	
	
	// Get the branches we need from the tree
	// TODO we don't need to call this for each entry
	SuperManager* Smgr = SuperManager::GetManager();
	TreeManager* mgr = Smgr->GetTreeManager(lun);
	TQReal *tqreal = mgr->GetTQREALINFO();
	MCInfo *mc = mgr->GetMC();
	int evid = 0;
	TBranch *ev_id = t->Branch("evid",&evid,"evid/D");
	
	if((nread%1)==0){
		Log(m_unique_name+" recostructing event "+toString(nread),v_message,m_verbose);
	}
	++nread;
	
	// Set the run number if run number is not real/applicable
	if(MC && skhead_.nrunsk==999999){
		//Log(m_unique_name+" warning: no run number!!",v_warning,m_verbose);
		skhead_.nrunsk = 75000;
	}
	// Putting this here temporarily for SKG4 sim
	if(MC && skhead_.nrunsk==85000){
		//Log(m_unique_name+" warning: no run number!!",v_warning,m_verbose);
		skhead_.nrunsk = 75000;
	}
	// once per run update the water transparency and SLE threshold
	if(skhead_.nrunsk!=nrunsk_last){
		// Get the SLE trigger threshold for this run
		int idetector[32], ithr[32], it0_offset[32], ipret0[32] ,ipostt0[32];
		softtrg_get_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);
		SLE_threshold = ithr[2];
		Log(m_unique_name+": SLE threshold"+toString(SLE_threshold),v_debug,m_verbose);
		
		// Update the water transparency
		int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];
		lfwater_(&days_to_run_start, &watert);
		Log(m_unique_name+" loaded new water transparency value "+toString(watert)
			+" for run "+toString(skhead_.nrunsk),v_debug,m_verbose);
		nrunsk_last = skhead_.nrunsk;
	}
	
	if(MC && (skhead_.nrunsk!=nrunsk_last || skhead_.nsubsk!=nsubsk_last)){
		int ierr;
		skbadch_(&skhead_.nrunsk,&skhead_.nsubsk,&ierr);
		nrunsk_last = skhead_.nrunsk;
		nsubsk_last = skhead_.nsubsk;
		if(skhead_.nrunsk!=nrunsk_last) darklf_(&skhead_.nrunsk);
		// Get the dark rate
		// TODO sort this out for later geometries
		int sk_geometry = skheadg_.sk_geometry;
		if (sk_geometry>=4) darkmc = mc->darkds*0.7880; //(1/1.269)
		if (darkmc==0) darkmc = 5000;
	}
	
	// Start by clearing all variables
	
	// set up variables for various fits
	float bswallsk=0;	// distance to wall
	float bseffh=0;		// effective hits
	float dist=0;		// distance to Sun in AU
	float bseffwal=0;	// effective distance to poswal
	int bsn20raw=0;		// raw hits in 20 ns window around peak of light
	int nflf = 0;		// number of flashers?
	float bsovaq = 0;	// fit quality derived from bsgood and bsdirks
	int bsn20 = 0;		// no. hits within 20 ns around peak of light
	float bsenergy_lfdir2=0;// energy from lfdir2
	float amsg = -1; 	// if the direction fit was successful
	float aratio = 0;	// ratio of used directions
	float acosscat = 0; 	// cos of max scattering angle
	int anscat = 0;		// number of possible directions
	
	float bsdir_lfdir2[3] = {0};// direction from first fit
	float poswal[3] = {0};	// extrapolated position at ID or OD wall
	float adir[3] = {0};  	// direction fit
	
	int lnqisk = skq_.nqisk;// total number of hits in event
	
	// mc vertex, direction, momentum, energy variables
	float *posmc = skroot_lowe_.posmc; //MC vtx (mm) 1 particle
	float posmcAFT[3]; 
	float (*dirmc)[3] = skroot_lowe_.dirmc; //MC dir 2 particles
	float *pabsmc = skroot_lowe_.pabsmc; 	//MC absolute momentum 1st and 2nd particles
	float *energymc = skroot_lowe_.energymc; // MC energy (MeV) 1st and 2nd particles
	
	float bsvertexAFT[4];
	float bsresultAFT[6];
	float bsgoodAFT[3];
	float bsvertexCombined[4];
	float *bsvertex = skroot_lowe_.bsvertex;
	float *bsresult = skroot_lowe_.bsresult;
	float *bsgood = skroot_lowe_.bsgood;
	
	// clear variables for low & mu fitters
	// also clears mc variables
	// sets bs variables to 9999
	lfclear_all_();
	for (int i=0;i<6;i++) bsresultAFT[i]=9999;
	for (int i=0;i<4;i++){
		bsvertexAFT[i]=9999;
		bsvertexCombined[i]=9999;
	}
	for (int i=0;i<3;i++){ 
		bsgoodAFT[i]=9999;
		posmcAFT[i] = 9999;
	}
	// clear variables to be saved to the ntuples
	x=9999; y=9999; z = 9999;
	x_prev=9999; y_prev=9999; z_prev = 9999;
	x_combined=9999; y_combined=9999; z_combined = 9999;
	mcx=9999; mcy=9999; mcz = 9999;
	mcx_prev=9999; mcy_prev=9999; mcz_prev = 9999;
	pdg=9999; pdg_prev = 9999;
	mcx_ncapture=9999; mcy_ncapture=9999; mcz_ncapture = 9999;
	mc_energy=9999; mc_energy_prev = 9999;
	mct_ncapture=9999;
	tgood=9999; tgood_prev=9999; tgood_combined=9999; tgood_combined_prev=9999;
	
	int max_nqisk_for_clusfit = 800; //(MC) ? 800 : 1000;  //  max no. hits to reconstruct
	
	//--------------Get the MC information we want----------------------//
	// Assumes IBD events - TODO does this need to work for other types?
	// NB in SKG4 the mc info is stored in reverse order:
	// neutron, positron, antineutrino for IBD events
	// This saves the info for the positron, and neutron where applicable
	if (MC) {
		
		float p_x = 9999;
		float p_y = 9999;
		float p_z = 9999;
		float p_abs = 9999;
		
		double restmass;
		
		// Get the MC info for the prompt event
		// This part is the same for SKG4 and SKDetsim
		
		// Get the mc x, y and z co-ordinates of the prompt
		posmc[0] = mc->pvtxvc[1][0];
		posmc[1] = mc->pvtxvc[1][1];
		posmc[2] = mc->pvtxvc[1][2];
		
		// Calculate the absolute momentum pabs of the prompt
		p_x = mc->pvc[1][0];//x momentum
		p_y = mc->pvc[1][1];//y momentum
		p_z = mc->pvc[1][2];//z momentum
		p_abs = sqrt(pow(p_x,2)+pow(p_y,2)+pow(p_z,2));
		pabsmc[0] = p_abs;
		
		//calculate the mc direction
		if (p_abs>0) {
			dirmc[0][0] = p_x/p_abs;
			dirmc[0][1] = p_y/p_abs;
			dirmc[0][2] = p_z/p_abs;
		}
		
		// Accessing the energy is different in skdetsim and SKG4
		if (mc->energy[0]!=0) { 
			//In SKG4, the energy is already calculated
			energymc[0] = mc->energy[1]; // positron energy
			pdg_prev = mc->ipvc[1];
		} else {
			// In SKDetsim, mc energy is 0 - we have to calculate it
			// Get the particle rest mass and calculate energy
			pdg_prev = mc->ipvc[1]; // positron
			restmass = PdgToMass(pdg_prev);
			energymc[0] = sqrt(pow(restmass,2)+pow(p_abs,2));
		}
		
		// Get the neutron MC info
		// Here, the whole process is different for SKG4 and SKDetim
		if (mc->energy[0]!=0) {
			// SKG4 (energy is zero in SKDetsim)
			// In SKG4, vtx and momentum are in secondaries branch
			SecondaryInfo* secondaries = mgr->GetSECONDARY();
			
			// Get the mc vertex of the neutron-capture event
			// For this, I'm going to use the creation vertex of
			// the gammas produced in the neutron capture.
			// In SKG4, the neutron is the parent in the Secdata tree 
			// and the gammas, conversion electrons, deuteron, Gd-158, etc are secondaries.
			int nsecondaries = secondaries->nscndprt; // number of secondaries (usually gammas)
			for (int isecondary = 0; isecondary<nsecondaries; isecondary++){
				pdg = secondaries->iprntprt[isecondary]; // parent of the secondary particle
				int pdg_secondary = secondaries->iprtscnd[isecondary]; // pdg of the secondary
				bool from_ncapture = secondaries->lmecscnd[isecondary]==18; // production process
				if (pdg_secondary == 22 && from_ncapture && pdg==2112){
					// save the creation vertex of the gamma from the neutron capture
					mcx_ncapture = secondaries->vtxscnd[isecondary][0];
					mcy_ncapture = secondaries->vtxscnd[isecondary][1];
					mcz_ncapture = secondaries->vtxscnd[isecondary][2];
					posmcAFT[0] = secondaries->vtxprnt[isecondary][0];
					posmcAFT[1] = secondaries->vtxprnt[isecondary][1];
					posmcAFT[2] = secondaries->vtxprnt[isecondary][2];
					mct_ncapture = secondaries->tscnd[isecondary];
					Log(m_unique_name+" ncapture time "+toString(mct_ncapture),v_debug,m_verbose);
					break; // stop looking once we have found one
				}
			}
			
			// Calculate the absolute initial momentum of the delayed
			p_x = mc->pvc[0][0];//x momentum
			p_y = mc->pvc[0][1];//y momentum
			p_z = mc->pvc[0][2];//z momentum
			p_abs = sqrt(pow(p_x,2)+pow(p_y,2)+pow(p_z,2));
			pabsmc[1] = p_abs;
			
			// Calculate the initial mc direction of the delayed
			if (p_abs>0) {
				dirmc[1][0] = p_x/p_abs;
				dirmc[1][1] = p_y/p_abs;
				dirmc[1][2] = p_z/p_abs;
			}
			
			// Then get the initial energy of the neutron from the MC branch
			energymc[1] = mc->energy[0];
		} else {
			// SKDetsim. NB outputs the initial vertex, not the neutron-capture vertex
			// Get the mc x, y and z co-ordinates of the delayed
			posmcAFT[0] = mc->pvtxvc[2][0];
			posmcAFT[1] = mc->pvtxvc[2][1];
			posmcAFT[2] = mc->pvtxvc[2][2];
			
			//calculate the absolute momentum of the delayed
			p_x = mc->pvc[2][0];//x momentum
			p_y = mc->pvc[2][1];//y momentum
			p_z = mc->pvc[2][2];//z momentum
			p_abs = sqrt(pow(p_x,2)+pow(p_y,2)+pow(p_z,2));
			pabsmc[1] = p_abs;
			
			//calculate the mc direction
			if (p_abs>0){
				dirmc[1][0] = p_x/p_abs;
				dirmc[1][1] = p_y/p_abs;
				dirmc[1][2] = p_z/p_abs;
			}
			
			// Calculate the energy
			pdg = mc->ipvc[2];
			restmass = PdgToMass(pdg);
			energymc[1] = sqrt(pow(restmass,2)+pow(p_abs,2));
		}
	} // endif MC
	mcx_prev = posmc[0];
	mcy_prev = posmc[1];
	mcz_prev = posmc[2];
	mcx = posmcAFT[0];
	mcy = posmcAFT[1];
	mcz = posmcAFT[2];
	mc_energy = energymc[1];
	mc_energy_prev = energymc[0];
	
	// Initialise with zeros the vectors which will be passed 
	// to the fitter as arrays 
	std::vector<float> charges;  // final array of hit charges
	std::vector<float> times;    // final array of hit times
	std::vector<int> cableIDs;   // final array of PMT numbers
	
//	lbfset0(cableIDs,&nqisk); // not needed if we're using vectors
//	lbfset0(charges,&nqisk);
//	lbfset0(times,&nqisk);
	int nhit=0;
	
	// Check if the first trigger is SHE
	std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
	int SHE = trigger_bits.test(28);
	int SLE = trigger_bits.test(2);
	int OD = trigger_bits.test(3);
	if (SLE && !OD) {
		/*************************************************************/
		// Reconstruct the prompt first
		
		// Fill the arrays of hit info which will be passed to fitter
		// Get the hits stored in the common blocks (skt_/skq_).
		// These are in-gate hits with bad/missing channels removed.
		SetPromptHits(charges,times,cableIDs,nhit);
		// Get the last hit in the SHE event.
		// This will be used to start the search for an AFT.
		float first_prompt = *min_element(begin(times),end(times));
		float last_prompt = *max_element(begin(times),end(times));	
		last_prompt = mct_ncapture; // using the mc time for now, to check the fit
		
		// MCInfo.prim_pret0 saves the time difference between the trigger time and geant_t0
		Log(m_unique_name+": first prompt hit time "+toString(first_prompt)+", pret0 "+toString(mc->prim_pret0[0]),v_debug,m_verbose);
		float prompt_trigger_t0 = mc->prim_pret0[0]+500.;
		
		// Do the single-event BONSAI fit
		int nsel;
		float bscharges[nhit];
		float bstimes[nhit];
		int bscableIDs[nhit];
		std::copy(charges.begin(),charges.end(),bscharges);
		std::copy(times.begin(),times.end(),bstimes);
		std::copy(cableIDs.begin(),cableIDs.end(),bscableIDs);
		if (nhit>9 && nhit<=max_nqisk_for_clusfit) {
			int nbf = bonsaifit_(&bsvertex[0],&bsresult[0],&bsgood[0],&nsel,&nhit,&bscableIDs[0],&bstimes[0],&bscharges[0]);
		}
		tgood_prev = bsgood[1];
		
		/*************************************************************/
		// Now reconstruct the AFT, if there is one
		// Get the raw hits from the tqreal branch in MC
		// Create vectors to store the info about all hits
		// 'raw' in the following means pre-skread i.e. all hits
		int nhitsAll = tqreal->nhits;
		std::vector<int> cableIDsRaw;
		std::vector<float> chargesRaw;
		std::vector<float> timesRaw;
		int nhitsRaw = 0;
		for (int ihit = 0;ihit<nhitsAll;ihit++){
			
			// Get the raw values from the tqreal branch
			// NB these are raw in the sense that they are for all 
			// of the hits, although in data some pre-processing
			// has been done by this stage e.g. to convert charge
			// to pe.(by skrawread)
			if (tqreal->T[ihit] >= last_prompt){
				chargesRaw.push_back(tqreal->Q[ihit]);
				cableIDsRaw.push_back(tqreal->cables[ihit] & ( (1<<16)-1));//cable ID (i.e. PMT number) is stored in 1st 16 bits
				timesRaw.push_back(tqreal->T[ihit]);
				nhitsRaw++;
			}
		}
		nhitsAFT_raw = timesRaw.size();
		Log(m_unique_name+": number of hits in the AFT before dark rate - "+toString(timesRaw.size()),v_debug,m_verbose);
		// Get the hits for the peak number of events in 200 ns
		std::vector<int> cableIDsAFT;
		std::vector<float> chargesAFT;
		std::vector<float> timesAFT;
		int nhitsAFT = 0;
		if (nhitsRaw>0){//SLE_threshold){
			nhitsAFT = SetAftHits(prompt_trigger_t0,SLE_threshold,addNoise,numPMTs,darkmc,last_prompt,chargesRaw,timesRaw,cableIDsRaw,nhitsRaw,chargesAFT,timesAFT,cableIDsAFT);
		}
		
		Log(m_unique_name+": number of in-gate hits in the AFT trigger - "+toString(nhitsAFT),v_debug,m_verbose);
//		for (int hit=0;hit<nhitsAFT;hit++){
//			Log(m_unique_name+"AFT hit time,charge,cable: "+toString(timesAFT[hit])+", "+toString(chargesAFT[hit])+", "+toString(cableIDsAFT[hit]),v_debug,m_verbose);
//			timesAFT[hit]-=(timesAFT[0]-bstimes[0]);
//		}
		// Do the single fit
		float bschargesAFT[nhitsAFT];
		float bstimesAFT[nhitsAFT];
		int bscableIDsAFT[nhitsAFT];
		std::copy(chargesAFT.begin(),chargesAFT.end(),bschargesAFT);
		std::copy(timesAFT.begin(),timesAFT.end(),bstimesAFT);
		std::copy(cableIDsAFT.begin(),cableIDsAFT.end(),bscableIDsAFT);
		nsel=0;
		if (bsvertex[0]<9999 && nhitsAFT>9 && nhitsAFT<=max_nqisk_for_clusfit) {
			int nbf = bonsaifit_(&bsvertexAFT[0],&bsresultAFT[0],&bsgoodAFT[0],&nsel,&nhitsAFT,&bscableIDsAFT[0],&bstimesAFT[0],&bschargesAFT[0]);
			tgood = bsgoodAFT[1];
			
			/*********************************************************/
			// Now that we have done the single fit for each of the SHE and AFT,
			// do the combined fit
			goodness* bsgdnSHE = new goodness(bspairlike->sets(),bspairlike->chargebins(),bsgeom,nhit,bscableIDs,bstimes,bscharges);
			fourhitgrid *bsgridSHE = new fourhitgrid(bsgeom->cylinder_radius(),bsgeom->cylinder_height(),bsgdnSHE);
			goodness *bsgdnAFT = new goodness(bspairlike->sets(),bspairlike->chargebins(),bsgeom,nhitsAFT,bscableIDsAFT,bstimesAFT,bschargesAFT);
			fourhitgrid *bsgridAFT = new fourhitgrid(bsgeom->cylinder_radius(),bsgeom->cylinder_height(),bsgdnAFT);
			combinedgrid *bspairgrid;
			bspairgrid = new combinedgrid(bsgeom->cylinder_radius(),
			                              bsgeom->cylinder_height(),
			                              bsgridSHE,bsgridAFT);
			int nselSHE = bsgdnSHE->nselected();
			int nselAFT = bsgdnAFT->nselected();
			bspairlike->set_hits(bsgdnSHE,bsgdnAFT);
			int useAngle = 1;
			bspairlike->maximize(bspairfit, bspairgrid,useAngle);
			bsvertexCombined[0] = bspairfit->xfit();
			bsvertexCombined[1] = bspairfit->yfit();
			bsvertexCombined[2] = bspairfit->zfit();
			bspairlike->ntgood(0,bsvertexCombined,darkmc,tgood_combined_prev);
			bspairlike->ntgood(1,bsvertexCombined,darkmc,tgood_combined);
		}
		
		// Now get the n50 and n10 values from the combined fit
		if (bsvertex[0]<9999){
			std::vector<int> cableIDs_bsn50;
			bsn50 = CalculateNX(50,bsvertex,cableIDs,times,cableIDs_bsn50);
		
			std::vector<int> cableIDs_bsn10;
			bsn10 = CalculateNX(10,bsvertex,cableIDs,times,cableIDs_bsn10);
		}
		if (bsvertexAFT[0]<9999){
			std::vector<int> cableIDs_bsn50AFT;
			bsn50AFT = CalculateNX(50,bsvertexAFT,cableIDsAFT,timesAFT,cableIDs_bsn50AFT);
			std::vector<int> cableIDs_bsn10AFT;
			bsn10AFT = CalculateNX(10,bsvertexAFT,cableIDsAFT,timesAFT,cableIDs_bsn10AFT);
		}
		if (bsvertexCombined[0]<9999){
			std::vector<int> cableIDs_n50_prev;
			n50_prev = CalculateNX(50,bsvertexCombined,cableIDs,times,cableIDs_n50_prev);
			std::vector<int> cableIDs_n10_prev;
			n10_prev = CalculateNX(10,bsvertexCombined,cableIDs,times,cableIDs_n10_prev);
			std::vector<int> cableIDs_n50_aft;
			n50_aft = CalculateNX(50,bsvertexCombined,cableIDsAFT,timesAFT,cableIDs_n50_aft);
			std::vector<int> cableIDs_n10_aft;
			n10_aft = CalculateNX(10,bsvertexCombined,cableIDsAFT,timesAFT,cableIDs_n10_aft);
		}
		
	} // SHE
	
	// Save the results of the fit
	x_prev = bsvertex[0];
	y_prev = bsvertex[1];
	z_prev = bsvertex[2];
	x = bsvertexAFT[0];
	y = bsvertexAFT[1];
	z = bsvertexAFT[2];
	x_combined = bsvertexCombined[0];
	y_combined = bsvertexCombined[1];
	z_combined = bsvertexCombined[2];
	
	outputTree->Fill();
	
	// TODO move the rest of the fitting to loop over both triggers
	
	// First, we get the atm hit parameters
	int nvalid,listatm[3],nhitatm[3];
	int atmchnum[3],atmratio[3];
	int num = 3;
	skatmmap2_(&num,&nvalid,&listatm[0],&nhitatm[0],&atmratio[0],&atmchnum[0]);
	skroot_lowe_.latmnum = atmchnum[0];
	skroot_lowe_.latmh = nhitatm[0];
	
	// Get the number of hits around a specified cable mx24
	// mxqisk # of max Q pmt
	int mxhit8;
	lfneihit_(&skq_.mxqisk,&mxhit8,&skroot_lowe_.lmx24);
//	tdiff(nt48sk,ltimediff);
	
	// NS ratio
	// calculate the ratio of hits in q bins
	int nqrng[3]; // nhit[3] number of hit divided region of q in 600ns & 1800ns
	// nqrng[0] <= -1 pC, -1 pC <= nqrng[1] <= 1 pC, nqrng[2] > 1 pC
	float rqrng[3]; // ratio of hits for each q bin
	lfqhit_(&nqrng[0],&rqrng[0]);
	if ((nqrng[0]+nqrng[1]+nqrng[2])>0)
	{
		skroot_lowe_.lnsratio = float(nqrng[1])/float(nqrng[0]+nqrng[1]+nqrng[2]);
	}
	
	if (nhit>9 && nhit<=max_nqisk_for_clusfit) {
	
		//------------------------------------------------------------------------
		// Now get the other values based on the reconstructed vertex
		//------------------------------------------------------------------------
		
		// Calculate n50 (hits in (-10,40) ns around the fitted time of emission
		// lfallfit:283 lfnhit(50.,bsvertex,bsn50,bshit) 
		// lfnhit(twindow,vertex,nhit,ihitcab) 
		// input: twindow (width of timing window), vertex[3] (x,y,z vertex)
		// output: nhit (maximal number of hits in timing window), ihitcab[nhit] (array with cable numbers of these hits)
		if (bsvertex[0]<9999)
		{
			std::vector<int> cableIDs_n50;
			Log(m_unique_name+" calculating NX",v_debug,m_verbose);
			skroot_lowe_.bsn50 = CalculateNX(50,bsvertex,cableIDs,times,cableIDs_n50);
			
			// TODO can we avoid using more of the fortran routines?
			
			//--------------------------------------------------------------------------
			// These parameters are saved to LOWE linfo array
			//--------------------------------------------------------------------------
			// There is an issue saving some of the parameters to linfo.
			// linfo is an int array but some of the variables stored in linfo are
			// floats or doubles. While some of the floats seem to be stored fine in
			// linfo, others are not. For now, I have just put a note by the ones that 
			// are not.
			
			// Distance to wall (needed for energy calculation)
			bswallsk = wallsk_(&bsvertex[0]);
			if (skroot_lowe_.lninfo < 10 )
			{
				skroot_lowe_.lninfo = 10;
			}
			
			if (bswallsk>0.0)
			{
				// Direction from lfdir2: "direction as 0th version"
				// Calculates the direction (costheta) from the hit positions and bsvertex
				// Also calculates dirks - essentially (there are corrections) sum over all hits of:
				// alog10(exp(linear function of costheta)/2000.)*pmt angular-dependent efficiency*charge
				lfdir2_(&bsvertex[0],&bsdir_lfdir2[0],&skroot_lowe_.bsdirks);
				
				// Energy based on lfdir2
				skroot_lowe_.lwatert = watert;
				// TODO use the correct versions for the geometry
				int flag_log = 0;
				int flag_skip = 0;
				if (MC)
				{
					lfneweff_sk3_final_(&bsvertex[0],&bsdir_lfdir2[0],&skroot_lowe_.bseffhit[0]);// but for sk6 and higher?
					lfeffwt_sk3_(&skroot_lowe_.bseffhit[0],&skroot_lowe_.lwatert,&bseffh);// but for sk6 and higher?
				}
				else
				{
					lfneweff_sk4_final_qe43_(&bsvertex[0], &bsdir_lfdir2[0], &skroot_lowe_.bseffhit[0], &skroot_lowe_.lwatert, &flag_log, &flag_skip);
					bseffh = skroot_lowe_.bseffhit[0];
				}
				int enelf_version = 4;
				bsenergy_lfdir2 = enelf_sk4_(&bseffh,&enelf_version);
				if (skroot_lowe_.lninfo < 14)
				{
					skroot_lowe_.lninfo = 14;
				}
				
				// Direction as a function of reconstructed vertex, and energy from first direction fit: 
				// dir[0]=cosphi, dir[1]=sinphi, dir[2]=costheta
				lfdir4_(&bsvertex[0],&bsenergy_lfdir2,&skroot_lowe_.bsdir[0],&skroot_lowe_.bsdirks);
				// Calculate direction KS: measure of the spread of hits around the fitted direction
				if (skroot_lowe_.bsdirks != -1) // if initial direction fit was successful
				{
					// final calculation of bsdirks
					lfdirks_(&bsvertex[0],&skroot_lowe_.bsdir[0],&skroot_lowe_.bsdirks);
				}
				else
				{
					skroot_lowe_.bsdirks = 2.0; // why???
				}
				
				// Cossun
				slsundir3_(&skhead_.ndaysk[0], &skhead_.ntimsk[0],&skroot_lowe_.bsdir[0],&skroot_lowe_.lsdir[0],&skroot_lowe_.bscossun,&dist);
				
				// Energy from Bonsai vertex (direction passed to function but is unused)
				// First calculate effective hits as a function of vertex, dark rate/current, pmt occupancy, effective photocoverage
				if (MC)
				{
					lfneweff_sk3_final_(&bsvertex[0],&skroot_lowe_.bsdir[0],&skroot_lowe_.bseffhit[0]); // but for sk6 and higher?
					// Then adjust effective hits for water transparency
					lfeffwt_sk3_(&skroot_lowe_.bseffhit[0],&skroot_lowe_.lwatert,&bseffh); // but for sk6 and higher?
				}
				else 
				{
					lfneweff_sk4_final_qe43_(&bsvertex[0], &bsdir_lfdir2[0], &skroot_lowe_.bseffhit[0], &skroot_lowe_.lwatert, &flag_log, &flag_skip);
				}
				bseffh = skroot_lowe_.bseffhit[0];
				skroot_lowe_.bsenergy = enelf_sk4_(&bseffh,&enelf_version);
				
				// patlik: likelihood of n50 hit pattern given the reconstructed vertex, direction and energy
				int ierr;
				skroot_lowe_.bspatlik = patliklf_sk4_(&bsvertex[0],&skroot_lowe_.bsdir[0],&skroot_lowe_.bsenergy,&ierr);
				
				lfflasher_(&skroot_lowe_.lmx24,&skq_.qimxsk,&nflf);
				if (skroot_lowe_.lninfo < 3)
				{
					skroot_lowe_.lninfo = 3;
				}
				
				// bseffwall effective distance from all
				int ID = 1; // to calculate distance from inner-detector wall (black sheet)
				bseffwal = effwallf_(&ID,&bsvertex[0],&skroot_lowe_.bsdir[0],&poswal[0]);
				if (skroot_lowe_.lninfo < 36)
				{
					skroot_lowe_.lninfo = 36;
				}
				
				// N20 raw hits at bonsai vertex
				// TODO get this with CalculateNX
				int bshit[numPMTs];
				int twindow = 20;
				int hitCables[numPMTs];
				lfnhit_tisk_(&twindow, &bsvertex[0], &bsn20raw, &hitCables[0], &bshit[0]);
				
				/*
				// Commented out for the moment due to crash when calling lfnhit2_ 
				//bsr02
				float bsr02;
				int bsn20=0;
				lfnhit2_(&twindow,&bsvertex[0], &bsn20, &hitCables[0], &bshit[0]);
				if (bsn20 == 0)
				{
					bsr02 = -1.0;
				}
				else
				{
					bsr02 = lfcal_r02_(&bsn20,&hitCables[0],&bsvertex[0]);
				}
				
				// bonsai clik
				float bsclik = (float)bsn20raw / bseffh * bsr02;
				*/
				
				// bonsai ovaq (from timing and direction ks goodness)
				bsovaq = pow(skroot_lowe_.bsgood[1],2)-pow(skroot_lowe_.bsdirks,2);
				if (skroot_lowe_.lninfo < 26) 
				{
					skroot_lowe_.lninfo = 26;
				}
				
				// Direction fit with ariadne.cc
				// linfo=27-33 Ariadne
				if (skroot_lowe_.lninfo < 33)
				{
					skroot_lowe_.lninfo = 33;
				}
				
				// Find the peak number of hits in 20 ns
				std::vector<int> cableIDs_n20;
				int bsn20 = CalculateNX(20,bsvertex,cableIDs,times,cableIDs_n20);
				//lfnhit_(&twindow,&bsvertex[0],&bsn20,&bshit[0]); // TODO can I use the c++ one?
				
				// Do the ariadne direction fit
				// Get the positions of the n20 hits
				float hits[3*bsn20];
				if (bsn20 <= 390)
				{
					for (int hit=0; hit<bsn20; hit++)
					{
						//TODO -1 in lfallfit but do we need it here?
						hits[3*hit]=geopmt_.xyzpm[0][cableIDs_n20[hit]];
						hits[3*hit+1]=geopmt_.xyzpm[1][cableIDs_n20[hit]];
						hits[3*hit+2]=geopmt_.xyzpm[2][cableIDs_n20[hit]];
					} // end loop over hits
					float vertex[3];
					vertex[0] = bsvertex[0];
					vertex[1] = bsvertex[1];
					vertex[2] = bsvertex[2];
					lfariadne_(&vertex[0],&bsn20,(&hits)[0],&adir[0],&amsg,&aratio,&anscat,&acosscat);
					
//					Attempt at C++ version of ariadne - need to get ariadne to compile properly
//					float costheta = 0.719;// TODO check this in sk code
//					ariadne *ari = new ariadne(vertex,bsn20,hits,costheta);
//					ari->fit();
//					amsg = ari->dir_goodness();
//					if (amsg>=0)
//					{
//						ari->dir(adir);
//						acosscat = ari->cos_scat_angle();
//						aratio = ari->dir_quality();
//						anscat = ari->nr_scat();
//					}
					
				} // endif bsn20 <= 390
			} // endif bswallsk > 0.0
		} // endif bsvertex[0] < 9999
	} // endif nhit > 9 && nhit < 800
	// Fill linfo array
	skroot_lowe_.linfo[0]  = nflf;
	skroot_lowe_.linfo[2]  = lnqisk;
	skroot_lowe_.linfo[4]  = dist;
	skroot_lowe_.linfo[5]  = bseffwal;
	skroot_lowe_.linfo[7]  = bseffh;
	skroot_lowe_.linfo[9]  = bswallsk;
	skroot_lowe_.linfo[10] = bsdir_lfdir2[0]; // TODO bsdir_lfdir2 not saving
	skroot_lowe_.linfo[11] = bsdir_lfdir2[1]; // TODO not saving
	skroot_lowe_.linfo[12] = bsdir_lfdir2[2]; // TODO not saving
	skroot_lowe_.linfo[13] = bsenergy_lfdir2;
	skroot_lowe_.linfo[22] = bsn20raw;
//	skroot_lowe_.linfo[23] = bsr02;
//	skroot_lowe_.linfo[24] = bsclik;
	skroot_lowe_.linfo[25] = bsovaq;          // TODO bsovaq not saving to linfo array
	skroot_lowe_.linfo[26] = adir[0];         // TODO not saving
	skroot_lowe_.linfo[27] = adir[1];         // TODO not saving
	skroot_lowe_.linfo[28] = adir[2];         // TODO not saving
	skroot_lowe_.linfo[29] = amsg;
	skroot_lowe_.linfo[30] = aratio;          // TODO not saving
	skroot_lowe_.linfo[31] = anscat;          // TODO not saving
	skroot_lowe_.linfo[32] = acosscat;        // TODO not saving
	skroot_lowe_.linfo[33] = poswal[0];
	skroot_lowe_.linfo[34] = poswal[1];
	skroot_lowe_.linfo[35] = poswal[2];
	
	// TODO do we need the clusfit reconstruction here?
	//clusfit
	//float *clvertex = skroot_lowe_.clvertex;
	//float *clresult = skroot_lowe_.clresult;
	//float *cldir = skroot_lowe_.cldir;
	//float  clgoodness = skroot_lowe_.clgoodness;
	//float  cldirks = skroot_lowe_.cldirks;
	//float  *cleffhit = skroot_lowe_.cleffhit;
	//float  clenergy = skroot_lowe_.clenergy;
	//int    cln50 = skroot_lowe_.cln50;
	
	// pass reconstructed variables from skroot_lowe_ common block to skroot file
	// TODO or do this without the common block
	skroot_set_lowe_(&lun,                      &skroot_lowe_.bsvertex[0], &skroot_lowe_.bsresult[0],
	                 &skroot_lowe_.bsdir[0],    &skroot_lowe_.bsgood[0],   &skroot_lowe_.bsdirks,
	                 &skroot_lowe_.bseffhit[0], &skroot_lowe_.bsenergy,    &skroot_lowe_.bsn50,
	                 &skroot_lowe_.bscossun,    &skroot_lowe_.clvertex[0], &skroot_lowe_.clresult[0],
	                 &skroot_lowe_.cldir[0],    &skroot_lowe_.clgoodness,  &skroot_lowe_.cldirks,
	                 &skroot_lowe_.cleffhit[0], &skroot_lowe_.clenergy,    &skroot_lowe_.cln50,
	                 &skroot_lowe_.clcossun,    &skroot_lowe_.latmnum,     &skroot_lowe_.latmh,
	                 &skroot_lowe_.lmx24,       &skroot_lowe_.ltimediff,   &skroot_lowe_.lnsratio,
	                 &skroot_lowe_.lsdir[0],    &skroot_lowe_.spaevnum,    &skroot_lowe_.spaloglike,
	                 &skroot_lowe_.sparesq,     &skroot_lowe_.spadt,       &skroot_lowe_.spadll,
	                 &skroot_lowe_.spadlt,      &skroot_lowe_.spamuyn,     &skroot_lowe_.spamugdn,
	                 &skroot_lowe_.posmc[0],    &skroot_lowe_.dirmc[0],    &skroot_lowe_.pabsmc[0],
	                 &skroot_lowe_.energymc[0], &skroot_lowe_.darkmc,      &skroot_lowe_.islekeep,
	                 &skroot_lowe_.bspatlik,    &skroot_lowe_.clpatlik,    &skroot_lowe_.lwatert,
	                 &skroot_lowe_.lninfo,      &skroot_lowe_.linfo[0]);
	
	
	// TODO???? remove hits outside 1.3 microsec
	delete_outside_hits_();
	
	// store header & TQ info.
	// skroot_set_tree.F calls skroot_set_header, skroot_set_tqreal and skroot_set_tqareal
	// pulls variables from fortran common blocks and passes them to the TreeManager
	skroot_set_tree_(&lun);
	
	// invokes TTree::Fill. Only use it in SKROOT mode WRITE or COPY!
	skroot_fill_tree_(&lun);
	//ev_id->Fill();
	//mgr->fill_tree(); // crashing for some reason
	mgr->Clear(); // zero out structures for next entry ! don't forget !
	
	//prev_t_nsec = t_nsec; // TODO reinstate when time is sorted
	
	return true;
}



bool CombinedFitter::Finalise()
{
	
	// plots for sanity check?
	// write the output ntuples to file
	fout->Write();
	delete outputTree;
	fout->Close();
	delete fout;
	outputTree = 0;
	fout = 0;
	
	return true;
}


void CombinedFitter::lbfset0(void *dat, int *numdat){
	
	size_t nsize;
	int iset;
	nsize = ((size_t)(*numdat)) * (size_t)4;
	iset = 0;
	memset(dat,iset,(size_t)nsize);
}

int CombinedFitter::CalculateNX(int timewindow, float* vertex, std::vector<int> cableIDs, std::vector<float> times, std::vector<int>& cableIDs_twindow) {
	if (times.size() <= 0)
	{
		return 0;
	}
	
	// Based on Fortran routine lfnhit ($skoflroot/lowe/sklowe/lfnhit.F)
	float cns2cm =21.58333; // speed of light in medium
	// Find tof subtracted times for all hits at reconstructed vertex
	std::vector<float> tof;
	for (int hit=0; hit<times.size(); hit++) {
		tof.push_back(times[hit]-sqrt(pow((vertex[0]-geopmt_.xyzpm[cableIDs[hit]-1][0]),2)+pow((vertex[1]-geopmt_.xyzpm[cableIDs[hit]-1][1]),2)+pow((vertex[2]-geopmt_.xyzpm[cableIDs[hit]-1][2]),2))/cns2cm);
	}
	
	// Sort in ascending order
	auto tof_sorted = tof;
	sort(tof_sorted.begin(),tof_sorted.end());
	
	// Find the centre of the distribution
	
	int bsnwindow = 1;
	int hstart_test = 0;
	int hstart = 0; 
	int hstop = 0;
	while (hstart_test < times.size()-bsnwindow) {
		hstop = hstart_test+bsnwindow;
		while((hstop<times.size()) && (tof_sorted[hstop]-tof_sorted[hstart_test]<=timewindow))
		{
			hstart = hstart_test;
			bsnwindow++;
			hstop++;
		}
		hstart_test++;
	}
	
	hstop=hstart+bsnwindow-1;
	// Make a list of cable IDs for the hits in the time window
	int bsnwin = 0;
	for(int hit=0; hit<times.size(); hit++) {
		if (tof[hit]<tof_sorted[hstart]) continue;
		if (tof[hit]>tof_sorted[hstop]) continue;
		cableIDs_twindow.push_back(cableIDs[hit]);
		bsnwin++;
	}
	if (bsnwin!=bsnwindow){
	   Log(m_unique_name+": bsnwin error, "+toString(bsnwin)+" != "+toString(bsnwindow),v_error,m_verbose);
	}
	return(bsnwindow);
}



void CombinedFitter::SetPromptHits(std::vector<float>& charges, std::vector<float>& times, std::vector<int> &cableIDs, int &nhit)
{
	
	// Get the hits stored in the common blocks (skt_/skq_).
	// These are in-gate hits with bad/missing channels removed.
	nhit = skq_.nqisk;
	for (int i=0;i<skq_.nqisk; i++) {
		// skt_, skq_, skchnl_ common blocks
		cableIDs.push_back(skchnl_.ihcab[i]);
		charges.push_back(skq_.qisk[skchnl_.ihcab[i]-1]);
		times.push_back(skt_.tisk[skchnl_.ihcab[i]-1]);
	}
	
}

void CombinedFitter::SingleEventFit(std::vector<float> charges, std::vector<float> times, std::vector<int> cableIDs, int nhit, float *bsvertex, float *bsresult, float *bsgood)
{
	// Do the BONSAI fit
	// set the relevant variables for the reconstruction results
	// TODO decide which variables we need
	int nsel;
	
	// Bonsai
	float bscharges[nhit];
	float bstimes[nhit];
	int bscableIDs[nhit];
	std::copy(charges.begin(),charges.end(),bscharges);
	std::copy(times.begin(),times.end(),bstimes);
	std::copy(cableIDs.begin(),cableIDs.end(),bscableIDs);
	int nselected=0;
	int nbf = bonsaifit_(&bsvertex[0],&bsresult[0],&bsgood[0],&nselected,&nhit,&bscableIDs[0],&bstimes[0],&bscharges[0]);
	/*else {
		goodness *bshits = new goodness(bslike->sets(),bslike->chargebins(),bsgeom,nhit,bscableIDs,bstimes,bscharges);
		int nsel = bshits->nselected();
		if (bshits->nselected()<4) { 
			return false;
		}
		fourhitgrid* bsgrid = new fourhitgrid(bsgeom->cylinder_radius(),bsgeom->cylinder_height(),bshits);
		
		bslike->set_hits(bshits);
		bslike->maximize(bsfit,bsgrid);
		int nbf = bslike->nfit();
		if (nbf >0) {
			// get best reconstructed vertex and time
			*bsvertex =bsfit->xfit(); 								//reco x
			bsvertex[1] = bsfit->yfit(); 							//reco y
			bsvertex[2] = bsfit->zfit();							//reco z
			bsvertex[3] = bslike->get_zero();//is there an offset?  //reco t0
			// get maximal log likelihood
			float gdn[bslike->sets()];
			bsgood[2] = bslike->goodness(*bsgood,bsvertex,gdn);		//not sure what all the bsgood vals are
			bslike->tgood(bsvertex,0,bsgood[1]);					//position (timing) goodness
			//bsnsel[1] = bslike->nwind(bsvertex,-6,12);
			*skroot_lowe_.bsgood = bsfit->maxq();								//max charge???
			bsfit->fitresult();										
			bslike->get_dir(skroot_lowe_.bsresult);								//direction
			skroot_lowe_.bsresult[5] = bslike->get_ll0();						//maximal log likelihood
			bslike->set_hits(NULL);				
		}
		
	}*/
	
}

int CombinedFitter::SetAftHits(float prompt_trigger_t0, int SLE_threshold, int addNoise, int numPMTs, float darkmc, float last_prompt, std::vector<float> chargesRaw, std::vector<float> timesRaw, std::vector<int> cableIDsRaw, int nhitsRaw, std::vector<float>& chargesAFT, std::vector<float>& timesAFT, std::vector<int>& cableIDsAFT)
{
	// Remove bad/missing channels and hits that were included in the prompt
//	std::vector<float> charges_tmp;
//	std::vector<float> times_tmp;
//	std::vector<int> cableIDs_tmp;
	std::vector<HitInfo> hits_tmp;
	
	// Alternative method to access the raw hits
	// Slightly different to getting from tqreal branch
	// (a few hits are missing this way - is it doing some of the below work for us?)
//	std::vector<float> chargesRawz;
//	std::vector<float> timesRawz;
//	std::vector<int> cableIDsRawz;
//	int nhit = sktqz_.nqiskz;
//	for (int i=0;i<nhit; i++) {
//		cableIDsRawz.push_back(sktqz_.icabiz[i]);
//		chargesRawz.push_back(sktqz_.qiskz[i]);
//		timesRawz.push_back(sktqz_.tiskz[i]);
//	}
	
	// Create a vector of IDs for removing bad channels
	std::vector<int> badIDs;
	for (int ibadch=0; ibadch< combad_.nbad; ibadch++){
		badIDs.push_back(combad_.isqbad[ibadch]);
	}
	
	// Loop over all raw hits after the SHE and remove OD hits, bad channels, etc
	for (int ihit = 0; ihit<nhitsRaw; ihit++){
		// only inner hits
		if (cableIDsRaw[ihit]<0 || cableIDsRaw[ihit]>=numPMTs)
			continue;
		// ignore hits included in the prompt
		if (timesRaw[ihit]<=last_prompt)
			continue;
		
		// Remove bad channels TODO missing channels
		bool bad = (find(badIDs.begin(),badIDs.end(),cableIDsRaw[ihit]) != badIDs.end());
		if (bad){
			Log(m_unique_name+" removing bad channel "+toString(cableIDsRaw[ihit]),v_debug,m_verbose);
			continue; 
		}
		
		// TODO remove missing channels
		//int *missingIDs = ?_.isqmis;
		//bool missing = (find(missingIDs.begin(),missingIDs.end(),cableIDsRaw[ihit]) != missingIDs.end());
		
		// save the selected hits into arrays
		hits_tmp.push_back({chargesRaw[ihit], timesRaw[ihit], cableIDsRaw[ihit]});
		
	}// nhitsRaw
	
	// Sort the hits in order of ascending time
	sort(hits_tmp.begin(),hits_tmp.end(),
		[](HitInfo const& i, HitInfo const& j) {return i.time < j.time;});
	
	if (addNoise){
		Log(m_unique_name+" Warning, adding dark noise",v_debug,m_verbose);
		TRandom rnd;
		// Get a random dark rate for this event
		float darkRate = numPMTs*darkmc;// total expected dark rate in Hz
		float ndark = rnd.Poisson(darkRate); // dark rate
		ndark *= 1e-9;
		Log(m_unique_name+" Dark hits per ns "+toString(ndark),v_debug,m_verbose);
		float tstartNoise = hits_tmp[0].time-300;
		float tendNoise = hits_tmp[hits_tmp.size()-1].time+100;
		float noiseWindow = tendNoise - tstartNoise; // total in the AFT
		ndark *= noiseWindow;// total 
		Log(m_unique_name+" Dark rate "+toString(darkRate),v_debug,m_verbose);
		Log(m_unique_name+" Noise window "+toString(noiseWindow),v_debug,m_verbose);
		Log(m_unique_name+" Total dark hits "+toString(ndark),v_debug,m_verbose);
		// loop over randomly generated dark hits and assign random dark
		// rate where event rate is below dark rate for hits in trigger
		for (int darkhit = 0; darkhit<ndark; darkhit++){
			
			// pick a random cable ID (PMT) and random time
			int cableDark = (int)(numPMTs*rnd.Rndm())+1;
			bool bad = (find(badIDs.begin(),badIDs.end(),cableDark) != badIDs.end());
			// Remove bad channels and hits not in the inner detector
			if (bad) 
				continue;
			if (cableDark<0 || cableDark>=numPMTs)
				continue;
			// TODO remove missing channels
			
			// pick a random time between the end of the SHE - 1000 ns (the neutron
			// capture at the moment, for testing purposes) and 1000 ns after 
			// the last neutron-capture hit
			float timeDark = tstartNoise + noiseWindow*rnd.Rndm();
			
			// Randomly generate a cable ID and charge and add the hit to 
			// our list of raw hits
			hits_tmp.push_back({abs(rnd.Gaus(1,0.7)), timeDark, cableDark});
			nhitsRaw++;
		} // end of loop over dark hits
		
		// Re-sort the hits in order of ascending time
		sort(hits_tmp.begin(),hits_tmp.end(),
			[](HitInfo const& i, HitInfo const& j) {return i.time < j.time;});
	} // endif
	
	
	// Find the 200 ns window with the maximum number of hits (n200Max)
	int n200max = 0.;
	float t_trigger = hits_tmp[0].time;
	float t_gate_start = t_trigger;
	float t_gate_end = t_trigger+1300;
	int triggerwindow = 200; //ns
	
	// loop over all of the hits after the prompt event to see if there is 
	// a SLE trigger
	int nhitsAFT = hits_tmp.size();
	int iend_tmp=0;
	
	for (int istart=0; istart< nhitsAFT-1; istart++){
		int n200=0;
		for (int iend=istart+1; iend<nhitsAFT; iend++){
			if ( (hits_tmp[iend].time - hits_tmp[istart].time) > triggerwindow )
				break;
			else{
				n200++;
				iend_tmp = iend;
			}
		}
		
		if (n200>n200max){
			n200max = n200;
			// Set the gate width (1.3 usec gate):
			// 	 |---------|---------------------------------------------|
			//-0.3usec  t_trigger                                   +1usec
			t_trigger = hits_tmp[iend_tmp].time;
			t_gate_start = t_trigger-300;
			t_gate_end = t_trigger + 1000;
		}
	}
	
	if (addNoise && n200max<SLE_threshold) return(0);
	if (!addNoise) n200max = hits_tmp.size();
	
	// Remove hits outside the 1.3 usec gate
	hits_tmp.erase(remove_if(begin(hits_tmp),end(hits_tmp),
		[t_gate_start,t_gate_end](HitInfo const& hit){
			return (hit.time<t_gate_start || hit.time>t_gate_end);
		}
	),hits_tmp.end());
	
	// Save the in-gate AFT hits 
//	float offset = t_trigger+prompt_trigger_t0+500;
	for ( auto ihit : hits_tmp) {
		cableIDsAFT.push_back(ihit.cableID);
		chargesAFT.push_back(ihit.charge);
		timesAFT.push_back(ihit.time-t_trigger + prompt_trigger_t0); // 500 time offset and make prompt and delayed trigger times the same
	}
	Log(m_unique_name+": prompt trigger "+toString(prompt_trigger_t0)+", delayed trigger time "+toString(t_trigger)+", time of first hit in delayed trigger "+toString(timesAFT[0]),v_debug,m_verbose);
	
	return(hits_tmp.size());
}
