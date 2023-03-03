/*
* Outputs to LOWE branch but the MC variables may differ from the fortran code.
* Assuming IBD:
* posmc[3]		- mc true position of the positron
* pabsmc[2]		- mc absolute momentum of the positron and neutron
* energymc[2]		- mc energy of the positron and neutron
* dirmc[2][3]		- mc true direction of the positron and neutron
*/

#include "CombinedFitter.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <string>
#include <iostream>
#include <bitset>
#include "SuperManager.h"

#include "Bonsai/searchgrid.h"
#include "Bonsai/bscalls.h"
#include "Bonsai/ariadne.h"
#include "TCanvas.h"

// declarations and #includes for SK fortran routines
#include "fortran_routines.h"
#include "SK_helper_functions.h"

using namespace std;

CombinedFitter::CombinedFitter():Tool(){

	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();

}


bool CombinedFitter::Initialise(std::string configfile, DataModel &data){


	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();

	m_data= &data;

	Log(toolName+": Initializing",v_debug,verbosity);

	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity",verbosity);    // how verbose to be
	m_variables.Get("readerName",readerName);  // name given to the TreeReader used for file handling
	m_variables.Get("dataSrc",dataSrc);   	   // where to get the data from (common blocks/tqreal)
	m_variables.Get("bonsaiSrc",bonsaiSrc);    // which bonsai to use (skofl or local)

	// use the readerName to find the LUN associated with this file
	std::map<std::string,int> lunlist;
	m_data->CStore.Get("LUNList",lunlist);
	if(lunlist.count(readerName)==0){
		Log(toolName+" error! No LUN associated with readerName "+readerName,v_error,verbosity);
		return false;
	}
	lun = lunlist.at(readerName);

	// determine if this file is MC or not, so we know whether or not to get
	// the MC information
	t = skroot_get_tree(&lun);
	if(t) MC = (t->FindBranch("MC")!=nullptr);

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
	
	// Use built-in BONSAI functions (bonsai/bscalls.cc::)
	if (bonsaiSrc==0) {
		std::cout << "Initialising BONSAI using built-in function cfbsinit_." << std::endl;
		cfbsinit_(&numPMTs, xyzpm);
	}

	// Use BONSAI direct (currently local)
	else{
		std::cout << "Initialising BONSAI direct." << std::endl;
		bsgeom = new pmt_geometry(numPMTs,xyzpm);
		bslike = new likelihood(bsgeom->cylinder_radius(),bsgeom->cylinder_height());
		bsfit = new bonsaifit(bslike);
	}
	std::cout << "BONSAI initialised" << std::endl;

	// check where we are getting the hit info from
	if (dataSrc==0) 
		Log(toolName+": Getting hit info from common blocks skt_/skq_",v_message,verbosity);
	else if (dataSrc==1) {
		Log(toolName+": Getting raw hit info from common block sktqz_",v_message,verbosity);
		Log(toolName+": Warning: this method is not yet functional (need to remove bad channels, etc)",v_warning,verbosity);
	}
	else {
		Log(toolName+": Getting raw hit info from tqreal branch",v_message,verbosity);
		Log(toolName+": Warning: this method may not give exactly the same result as skt_/skq_ common blocks",v_warning,verbosity);
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

	if((nread%100)==0){
		Log(toolName+" read loop "+toString(nread)+", current run "+toString(skhead_.nrunsk),v_message,verbosity);
	}
	++nread;

	// Set the run number if run number is not real/applicable
	if(MC && skhead_.nrunsk==999999){
		Log(toolName+" warning: no run number!!",v_warning,verbosity);
		skhead_.nrunsk = 75000;
	}
	if(MC && skhead_.nrunsk==85000){
		Log(toolName+" warning: no run number!!",v_warning,verbosity);
		skhead_.nrunsk = 75000;
	}

	// once per run update the water transparency
	if(skhead_.nrunsk!=nrunsk_last){
		int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];
		lfwater_(&days_to_run_start, &watert);
		Log(toolName+" loaded new water transparency value "+toString(watert)
	    	+" for run "+toString(skhead_.nrunsk),v_debug,verbosity);
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
		float darkmc;
		int sk_geometry = skheadg_.sk_geometry;
//		if (sk_geometry==4) darkmc = mc->darkds*0.7880; //(1/1.269)
//		else if (sk_geometry>=5) darkmc = mc->darkds*0.7880; //(1/1.269)
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

	// clear variables for low & mu fitters
	// also clears mc variables
	// sets bs variables to 99999
	lfclear_all_();

	int NHITCUT = 800; //(MC) ? 800 : 1000;  //  max no. hits to reconstruct
	
	//--------------Get the MC information we want----------------------//
	// Assumes IBD events - TODO does this need to work for other types?
	// NB in SKG4 the mc info is stored in reverse order:
	// neutron, positron, antineutrino for IBD events
	// This saves the info for the positron, and neutron where applicable
	if (MC) {
		// mc vertex, direction, momentum, energy variables
		float *posmc = skroot_lowe_.posmc; //MC vtx (mm) 1 particle
		float (*dirmc)[3] = skroot_lowe_.dirmc; //MC dir 2 particles
		float *pabsmc = skroot_lowe_.pabsmc; 	//MC absolute momentum \
							1st and 2nd particles
		float *energymc = skroot_lowe_.energymc; // MC energy (MeV) \
							1st and 2nd particles

		float p_x = 9999;
		float p_y = 9999;
		float p_z = 9999;
		float p_abs = 9999;
		double restmass;
		int pdg;
		
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
		} else {
			// In SKDetsim, mc energy is 0 - we have to calculate it
			// Get the particle rest mass and calculate energy
			pdg = mc->ipvc[1]; // positron
			restmass = PdgToMass(pdg);
			energymc[0] = sqrt(pow(restmass,2)+pow(p_abs,2));
		}

		// Get the neutron MC info
		// Here, the whole process is different for SKG4 and SKDetim
		if (mc->energy[0]!=0) {
			// SKG4 (energy is zero in SKDetsim)
			// In SKG4, vtx and momentum are in secondaries branch
			SecondaryInfo* secondaries = mgr->GetSECONDARY();
			
			// Get the mc vertex of the neutron event
			// TODO - how shall we save this info?
			posmc_d[0] = secondaries->vtxscnd[0][0];
			posmc_d[1] = secondaries->vtxscnd[0][1];
			posmc_d[2] = secondaries->vtxscnd[0][2];
			

			// Calculate the absolute momentum of the delayed
			p_x = mc->pvc[0][0];//x momentum
			p_y = mc->pvc[0][1];//y momentum
			p_z = mc->pvc[0][2];//z momentum
			p_abs = sqrt(pow(p_x,2)+pow(p_y,2)+pow(p_z,2));
			pabsmc[1] = p_abs;
			
			// Calculate the mc direction
			if (p_abs>0) {
				dirmc[1][0] = p_x/p_abs;
				dirmc[1][1] = p_y/p_abs;
				dirmc[1][2] = p_z/p_abs;
			}
		
			// Then get the energy of the neutron from the MC branch
			energymc[1] = mc->energy[0];
		} else {
			// SKDetsim
			// Get the mc x, y and z co-ordinates of the delayed
			// TODO - how shall we save this info?
			posmc_d[0] = mc->pvtxvc[2][0];
			posmc_d[1] = mc->pvtxvc[2][1];
			posmc_d[2] = mc->pvtxvc[2][2];

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

	// Initialise with zeros the vectors which will be passed 
	// to the fitter as arrays 
	vector<float> charges;  	// final array of hit charges
	vector<float> times;	// final array of hit times
	vector<int> cableIDs;	// final array of PMT numbers
	float *bsvertex = skroot_lowe_.bsvertex;
	float *bsresult = skroot_lowe_.bsresult;
	float *bsgood = skroot_lowe_.bsgood;

//	lbfset0(cableIDs,&nqisk); // not needed if we're using vectors
//	lbfset0(charges,&nqisk);
//	lbfset0(times,&nqisk);
	int nhit=0;

	// Check if the first trigger is SHE
	std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
	int SHE = trigger_bits.test(28);
	if (SHE) {
		/*************************************************************/
		// Reconstruct the prompt first

		// We shouldn't have explicitly get the SHE but may need it for timing
		// in data so leaving it here for now commented out. Seems to change
		// the hits so care must be taken if reinstating.
		
		/*int in_gate_only = 0;
		int hw_ctr = 0; // time of hw trigger (also SHE in MC) 
		int max_trig = 500;
		int num_trig = softtrg_inittrgtbl_(&skhead_.nrunsk, &hw_ctr, &in_gate_only, &max_trig); 
		double t0_SHE=0;
		// Get the time of the first trigger
		for (int trig = 0; trig<num_trig; trig++) {
			if (swtrgtbl_.swtrgtype[trig]==28){
				t0_SHE = swtrgtbl_.swtrgt0ctr[trig];
			}
		}
		
		Shouldn't need to explicitly get the hits for the SHE as they
		should already be loaded. Commenting out for now. Might be needed
		if reinstating the above
		Set the hits for the SHE trigger
		int new_gate = skheadqb_.it0sk + t0_SHE;
		set_timing_gate_(&new_gate);
		do skread and get the sktq_
		skread(-lun,*1003,*1003,*1003,*1003)
		skroot_set_tree(&lun);
		skroot_fill_tree(&lun);*/

		Log(toolName+"it0sk = "+toString(skheadqb_.it0sk/1.92/1000.),v_debug,verbosity);
		
		// Fill the arrays of hit info which will be passed to fitter
		// Get the hits stored in the common blocks (skt_/skq_).
		// These are in-gate hits with bad/missing channels removed.
		SetPromptHits(charges,times,cableIDs,nhit);
		float endSHE = *max_element(begin(times),end(times));	
		
		// Do the single-event BONSAI fit
		int nsel;
		float bscharges[nhit];
		float bstimes[nhit];
		int bscableIDs[nhit];
		std::copy(charges.begin(),charges.end(),bscharges);
		std::copy(times.begin(),times.end(),bstimes);
		std::copy(cableIDs.begin(),cableIDs.end(),bscableIDs);
		if (nhit<=NHITCUT) {
			int nbf = bonsaifit_(&bsvertex[0],&bsresult[0],&bsgood[0],&nsel,&nhit,&bscableIDs[0],&bstimes[0],&bscharges[0]);
		}


		/*************************************************************/
		// Now reconstruct the AFT, if there is one

		// Get the raw hits from the tqreal branch
		// Create vectors to store the info about all hits
		// 'raw' in the following means pre-skread i.e. all hits
		int nhitsRaw = tqreal->nhits;
		vector<int> cableIDsRaw;
		vector<float> chargesRaw;
		vector<float> timesRaw;

		for (int ihit = 0;ihit<nhitsRaw;ihit++){

			// Get the raw values from the tqreal branch
			// NB these are raw in the sense that they are for all 
			// of the hits, although in data some pre-processing
			// has been done by this stage e.g. to convert charge
			// to pe.(by skrawread)
			chargesRaw.push_back(tqreal->Q[ihit]);
			cableIDsRaw.push_back(tqreal->cables[ihit] & ( (1<<16)-1));//cable ID (i.e. PMT number) is stored in 1st 16 bits
			timesRaw.push_back(tqreal->T[ihit]);
		}

		// Get the hits for the peak number of events in 200 ns
		int nhits_AFT = 0;
		vector<int> cableIDsAFT;
		vector<float> chargesAFT;
		vector<float> timesAFT;
		int nhitsAFT = SetAftHits(endSHE,chargesRaw,timesRaw,cableIDsRaw,nhitsRaw,chargesAFT,timesAFT,cableIDsAFT);

		// Do the single fit
		float bschargesAFT[nhitsAFT];
		float bstimesAFT[nhitsAFT];
		int bscableIDsAFT[nhitsAFT];
		std::copy(chargesAFT.begin(),chargesAFT.end(),bschargesAFT);
		std::copy(timesAFT.begin(),timesAFT.end(),bstimesAFT);
		std::copy(cableIDsAFT.begin(),cableIDsAFT.end(),bscableIDsAFT);
		nsel=0;
		//TODO write a version of lf_clear_all_ to set these to 9999
		float bsvertexAFT[4]={9999,9999,9999,9999};
		float bsresultAFT[6]={9999,9999,9999,9999,9999,9999};
		float bsgoodAFT[3]={9999,9999,9999};
		if (nhitsAFT<=NHITCUT) {
			int nbf = bonsaifit_(&bsvertexAFT[0],&bsresultAFT[0],&bsgoodAFT[0],&nsel,&nhitsAFT,&bscableIDsAFT[0],&bstimesAFT[0],&bschargesAFT[0]);
		}
		
		/*********************************************************/
		// Now that we have done the single fit for each of the SHE and AFT,
		// do the combined fit
		// TODO
	
	} // SHE

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
//			tdiff(nt48sk,ltimediff);

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

	if (nhit>9 && nhit<=NHITCUT) {




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
			vector<int> cableIDs_n50;
			Log(toolName+" calculating NX",v_debug,verbosity);
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
				vector<int> cableIDs_n20;
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
			//		Attempt at C++ version of ariadne - need to get ariadne to compile properly
			//		float costheta = 0.719;// TODO check this in sk code
			//		ariadne *ari = new ariadne(vertex,bsn20,hits,costheta);
			//		ari->fit();
			//		amsg = ari->dir_goodness();
			//		if (amsg>=0)
			//		{
			//			ari->dir(adir);
			//			acosscat = ari->cos_scat_angle();
			//			aratio = ari->dir_quality();
			//			anscat = ari->nr_scat();
			//		}

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
//  skroot_lowe_.linfo[23] = bsr02;
//	skroot_lowe_.linfo[24] = bsclik;
	skroot_lowe_.linfo[25] = bsovaq;  	// TODO bsovaq not saving to linfo array
	skroot_lowe_.linfo[26] = adir[0]; 	// TODO not saving
	skroot_lowe_.linfo[27] = adir[1]; 	// TODO not saving
	skroot_lowe_.linfo[28] = adir[2]; 	// TODO not saving
	skroot_lowe_.linfo[29] = amsg;
	skroot_lowe_.linfo[30] = aratio;  	// TODO not saving
	skroot_lowe_.linfo[31] = anscat;  	// TODO not saving
	skroot_lowe_.linfo[32] = acosscat;	// TODO not saving
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
    std::cout << "Finished reconstructing event " << ev << std::endl;

    // terminate bonsai?
	// plots for sanity check?
    //TTree* t = skroot_get_tree(&lun);
	//TH1D *h1 = new TH1D();
	//t->Draw("sqrt(pow(LOWE.posmc[0]-LOWE.bsvertex[0],2)+pow(LOWE.posmc[1]-LOWE.bsvertex[1],2)+pow(LOWE.posmc[2]-LOWE.bsvertex[2],2))");
	ht->Draw();
  	return true;
}


void CombinedFitter::lbfset0(void *dat, int *numdat){

    size_t nsize;
    int iset;
    nsize = ((size_t)(*numdat)) * (size_t)4;
    iset = 0;
    memset(dat,iset,(size_t)nsize);
}

int CombinedFitter::CalculateNX(int timewindow, float* vertex, vector<int> cableIDs, vector<float> times, vector<int>& cableIDs_twindow) {
	if (skq_.nqisk <= 0)
	{
		return 0;
	}

	// Based on Fortran routine lfnhit ($skoflroot/lowe/sklowe/lfnhit.F)
	float cns2cm =21.58333; // speed of light in medium
	// Find tof subtracted times for all hits at reconstructed vertex
	vector<float> tof;
	Log(toolName+" getting tof subtracted times for all hits",v_debug,verbosity);
	for (int hit=0; hit<skq_.nqisk; hit++) {
		   tof.push_back(times[hit]-sqrt(pow((vertex[0]-geopmt_.xyzpm[cableIDs[hit]-1][0]),2)+pow((vertex[1]-geopmt_.xyzpm[cableIDs[hit]-1][1]),2)+pow((vertex[2]-geopmt_.xyzpm[cableIDs[hit]-1][2]),2))/cns2cm);
	}

	Log(toolName+" sorting tof subtracted times",v_debug,verbosity);
	// Sort in ascending order
	auto tof_sorted = tof;
	sort(tof_sorted.begin(),tof_sorted.end());
		   
	// Find the centre of the distribution
	Log(toolName+" finding the centre of the tof subtracted times distribution",v_debug,verbosity);

	int bsnwindow = 1;
	int hstart_test = 0 ;
	int hstart = 0; 
	int hstop = 0;
	while (hstart_test < skq_.nqisk-bsnwindow) {
		hstop = hstart_test+bsnwindow;
		while((hstop<skq_.nqisk) && (tof_sorted[hstop]-tof_sorted[hstart_test]<=timewindow))
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
    for(int hit=0; hit<skq_.nqisk; hit++) {
        if (tof[hit]<tof_sorted[hstart]) continue;
        if (tof[hit]>tof_sorted[hstop]) continue;
        cableIDs_twindow.push_back(cableIDs[hit]);
		bsnwin++;
    }
    if (bsnwin!=bsnwindow) printf("bsnwin error %d!=%d\n",bsnwin,bsnwindow);
    return(bsnwindow);
}



void CombinedFitter::SetPromptHits(vector<float>& charges, vector<float>& times, vector<int> &cableIDs, int &nhit)
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

void CombinedFitter::SingleEventFit(vector<float> charges, vector<float> times, vector<int> cableIDs, int nhit, float *bsvertex, float *bsresult, float *bsgood)
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
	printf("bsvertex[0] %f \n",bsvertex[0]);	
	/*else {
		goodness *bshits = new goodness(bslike->sets(),bslike->chargebins(),bsgeom,nhit,bscableIDs,bstimes,bscharges);
		int nsel = bshits->nselected();
		if (bshits->nselected()<4) { 
			std::cout << "Event " << ev << ", " << bshits->nselected() << " selected hits not enough." << std::endl;
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

int CombinedFitter::SetAftHits(float endSHE, vector<float> chargesRaw, vector<float> timesRaw, vector<int> cableIDsRaw, int nhitsRaw, vector<float>& chargesAFT, vector<float>& timesAFT, vector<int>& cableIDsAFT)
{
	// Remove bad/missing channels and hits that were included in the prompt
	int nhitsTmp = 0;
//	vector<float> charges_tmp;
//	vector<float> times_tmp;
//	vector<int> cableIDs_tmp;
	vector<HitInfo> hits_tmp;
	
	// Alternative method to access the raw hits
	// Slightly different to getting from tqreal branch
	// (a few hits are missing this way)
//	vector<float> chargesRawz;
//	vector<float> timesRawz;
//	vector<int> cableIDsRawz;
//	int nhit = sktqz_.nqiskz;
//	for (int i=0;i<nhit; i++) {
//		cableIDsRawz.push_back(sktqz_.icabiz[i]);
//		chargesRawz.push_back(sktqz_.qiskz[i]);
//		timesRawz.push_back(sktqz_.tiskz[i]);
//	}

	for (int ihit = 0; ihit<nhitsRaw; ihit++){
		// ignore hits with no real pmt id
		if (cableIDsRaw[ihit]<0 || cableIDsRaw[ihit]>=MAXPM)
			continue;
		// ignore hits included in the prompt
		if (timesRaw[ihit]<=endSHE)
			continue;
		// Remove bad channels TODO missing channels
		vector<int> badIDs;
		for (int ibadch=0; ibadch< combad_.nbad; ibadch++){
			badIDs.push_back(combad_.isqbad[ibadch]);
		}
		//int *missingIDs = ?_.isqmis;
		bool bad = (find(badIDs.begin(),badIDs.end(),cableIDsRaw[ihit]) != badIDs.end());
		//bool missing = (find(missingIDs.begin(),missingIDs.end(),cableIDsRaw[ihit]) != missingIDs.end());
		if (bad)
			continue; 

		// save the selected hits into arrays
		hits_tmp.push_back({chargesRaw[ihit], timesRaw[ihit], cableIDsRaw[ihit]});	
		nhitsTmp++;
	}// nhitsRaw
	
	// Sort the hits in order of ascending time
	sort(hits_tmp.begin(),hits_tmp.end(),
		[](HitInfo const& i, HitInfo const& j) {return i.time < j.time;});
		

	// TODO Find the 200 ns window with the maximum number of hits (n200Max)

	// TODO Save the hits in the peak 200 ns window

	// No dark noise in the MC so for now we are just going to save the AFT hits
	for ( auto ihit : hits_tmp) {
		Log(toolName+"AFT hit time: "+toString(ihit.time),v_debug,verbosity);
		cableIDsAFT.push_back(ihit.cableID);
		chargesAFT.push_back(ihit.charge);
		timesAFT.push_back(ihit.time);
	}	
	return(hits_tmp.size());
}
