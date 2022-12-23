#include "CombinedFitter.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"

#include <string>
#include <vector>
#include <iostream>
#include <bitset>

#include "SuperManager.h"

#include "Bonsai/searchgrid.h"
#include "Bonsai/bscalls.h"
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
    m_variables.Get("verbosity",verbosity);            // how verbose to be
    m_variables.Get("readerName",readerName);          // name given to the TreeReader used for file handling
	m_variables.Get("dataSrc",dataSrc);   				// where to get the data from (common blocks/tqreal)
	m_variables.Get("bonsaiSrc",bonsaiSrc);   			// which bonsai to use (skofl or local)

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
  	
    // TODO Set the output tree
	// Set up a tree to save the results of the reconstruction
	// and other relevant data as ntuples.
    /*
    if ( fname ) {
        fout = new TFile(fname, "RECREATE");
    }
    else fout = new TFile("SK2p2MeV_output.root", "RECREATE");

    theOTree = new TTree("sk2p2", "SK 2.2 MeV");
    // 2.2 MeV info
    theOTree->Branch("np",    &res.np,     "np/I");
    theOTree->Branch("N200M", &res.N200M,  "N200M/I");
    theOTree->Branch("T200M", &res.T200M,  "T200M/F");
    theOTree->Branch("N10",    res.N10,    "N10[np]/I");
    theOTree->Branch("N200",    res.N200,    "N200[np]/I");
    theOTree->Branch("N10d",    res.N10d,    "N10d[np]/I");
    theOTree->Branch("Nc",     res.Nc,     "Nc[np]/I");
    theOTree->Branch("Nback",   res.Nback,   "Nback[np]/I");
    theOTree->Branch("N300",   res.N300,   "N300[np]/I");
    theOTree->Branch("trms",   res.trms,   "trms[np]/F");
    theOTree->Branch("trmsdiff",   res.trmsdiff,   "trmsdiff[np]/F");
    theOTree->Branch("fpdist",   res.fpdist,   "fpdist[np]/F");
    theOTree->Branch("bpdist",   res.bpdist,   "bpdist[np]/F");
    theOTree->Branch("fwall",   res.fwall,   "fwall[np]/F");
    theOTree->Branch("bwall",   res.bwall,   "bwall[np]/F");
    theOTree->Branch("pvx",   res.pvx,   "pvx[np]/F");
    theOTree->Branch("pvy",   res.pvy,   "pvy[np]/F");
    theOTree->Branch("pvz",   res.pvz,   "pvz[np]/F");
    theOTree->Branch("bse",   res.bse,   "bse[np]/F");
    theOTree->Branch("mintrms_3",   res.mintrms_3,   "mintrms_3[np]/F");
    theOTree->Branch("mintrms_6",   res.mintrms_6,   "mintrms_6[np]/F");
    //*********new added by Yang Zhang****
    theOTree->Branch("Q10",   res.Q10,   "Q10[np]/F");
    theOTree->Branch("Qrms",   res.Qrms,   "Qrms[np]/F");
    theOTree->Branch("Qmean",   res.Qmean,   "Qmean[np]/F");
    theOTree->Branch("thetarms",   res.thetarms,   "thetarms[np]/F");
    theOTree->Branch("NLowtheta",   res.NLowtheta,   "NLowtheta[np]/I");
    //************************************
    theOTree->Branch("phirms", res.phirms, "phirms[np]/F");
    theOTree->Branch("bsdirks", res.bsdirks, "bsdirks[np]/F");
    theOTree->Branch("thetam", res.thetam, "thetam[np]/F");
    theOTree->Branch("dt",     res.dt,     "dt[np]/F");
    theOTree->Branch("dtn",     res.dtn,     "dtn[np]/F");
    theOTree->Branch("nvx",     res.nvx,     "nvx[np]/F");
    theOTree->Branch("nvy",     res.nvy,     "nvy[np]/F");
    theOTree->Branch("nvz",     res.nvz,     "nvz[np]/F");
    theOTree->Branch("tindex",     res.tindex,     "tindex[np]/I");
    theOTree->Branch("n40index",     res.n40index,     "n40index[np]/I");

    theOTree->Branch("Neff",   res.Neff,   "Neff[np]/I");
    theOTree->Branch("ratio",  res.ratio,  "ratio[np]/F");
    theOTree->Branch("Nc1",    res.Nc1,    "Nc1[np]/I");
    theOTree->Branch("NhighQ", res.NhighQ, "NhighQ[np]/I");
    theOTree->Branch("NlowQ", res.NlowQ, "NlowQ[np]/I");

    theOTree->Branch("Nlow1",   res.Nlow1,   "Nlow1[np]/I");
    theOTree->Branch("Nlow2",   res.Nlow2,   "Nlow2[np]/I");
    theOTree->Branch("Nlow3",   res.Nlow3,   "Nlow3[np]/I");
    theOTree->Branch("Nlow4",   res.Nlow4,   "Nlow4[np]/I");
    theOTree->Branch("Nlow5",   res.Nlow5,   "Nlow5[np]/I");
    theOTree->Branch("Nlow6",   res.Nlow6,   "Nlow6[np]/I");
    theOTree->Branch("Nlow7",   res.Nlow7,   "Nlow7[np]/I");
    theOTree->Branch("Nlow8",   res.Nlow8,   "Nlow8[np]/I");
    theOTree->Branch("Nlow9",   res.Nlow9,   "Nlow9[np]/I");
    theOTree->Branch("ncomb3", &res.ncomb3,  "ncomb3/I");
    //theOTree->Branch("theta", res.theta, "theta[ncomb3]/F");
    theOTree->Branch("cable", res.cable, "cable[ncomb3]/I");
    /////////////////////////////////////////////////////////////////////
    //     theOTree->Branch("ncomb", &res.ncomb,  "ncomb/I");
    //     theOTree->Branch("dalpha", res.dalpha, "dalpha[ncomb]/F");
    //     theOTree->Branch("ncomb1", &res.ncomb1,  "ncomb1/I");
    //     theOTree->Branch("dalpha1", res.dalpha1, "dalpha1[ncomb1]/F");
    //     theOTree->Branch("ncomb2", &res.ncomb2,  "ncomb2/I");
    //     theOTree->Branch("dalpha2", res.dalpha2, "dalpha2[ncomb2]/F");
    /////////////////////////////////////////////////////////////////////
*/


    // initialize water transparency table
	// (this will be for energy reconstruction I presume)
    skrunday_();
    skwt_gain_corr_();
   
    // get detector parameters from common blocks (geopmtC.h)
    int numPMTs = MAXPM; //number of PMTs in the detector
    float* xyzpm = &geopmt_.xyzpm[0][0];//lf_allfit_new.cc
//  nrunsk_last = skhead_.nrunsk;
//	nsubsk_last = skhead_.nsubsk;

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
	else if (dataSrc==1)
	{
		Log(toolName+": Getting raw hit info from common block sktqz_",v_message,verbosity);
		Log(toolName+": Warning: this method is not yet functional (need to remove bad channels, etc)",v_warning,verbosity);
	}
	else 
	{
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

    if(MC && skhead_.nrunsk==999999){
        Log(toolName+" warning: no run number!!",v_warning,verbosity);
        skhead_.nrunsk = 75000;
    }

    // once per run update the water transparency
    if(skhead_.nrunsk!=nrunsk_last){
        int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];  // defined in skdayC.h
        lfwater_(&days_to_run_start, &watert);
        Log(toolName+" loaded new water transparency value "+toString(watert)
            +" for run "+toString(skhead_.nrunsk),v_debug,verbosity);
        nrunsk_last = skhead_.nrunsk;
    }

    if(MC){
        if(skhead_.nrunsk!=nrunsk_last || skhead_.nsubsk!=nsubsk_last){
            int ierr;
            skbadch_(&skhead_.nrunsk,&skhead_.nsubsk,&ierr);
            nrunsk_last = skhead_.nrunsk;
            nsubsk_last = skhead_.nsubsk;
            if(skhead_.nrunsk!=nrunsk_last) darklf_(&skhead_.nrunsk);
		// Get the dark rate
		float darkmc;
		int sk_geometry = skheadg_.sk_geometry;
		if (sk_geometry==4) darkmc = mc->darkds*0.7880; // (1/1.269) lfallfit_sk4_mc.F::112
		else if (sk_geometry==5) darkmc = mc->darkds*0.7880; // (1/1.269) (tentative) lfallfit_sk4_mc.F::114
		else std::cout << "SK_GEOMETRY=" << sk_geometry << ". Not set up for earlier than SKIV" << std::endl;
        }
    }

    // clear variables for low & mu fitters
    // sets bs variables to 99999
    lfclear_all_();

	

    // apply lowfit
    int NHITCUT = (MC) ? 800 : 1000;  //  number of hit limit for clusfit (changed for relic analysis)

	// Iterate over the prompt and delayed events

	// TODO check how many AFT events there are
	// Also check that there is an AFT (or do this in TreeReader)
	
	int max_trig = 1280;

//	int has_she = skhead_.idtgsk&(1<<28);
//	int has_aft = skhead_.idtgsk&(1<0);
//  trigger_bits = skhead_.idtgsk;

//	if (trigger_bits.test(trig_id_she))
//	{
//		Log(toolName+" Has SHE",v_message,verbosity);
//	}
//	else if (!MC && trigger_bits.test(trig_id_aft))
//	{
//		Log(toolName+" Has AFT",v_message,verbosity);
//	}

	// TODO Now to extract the AFT info (for MC)
	if (MC)
	{
		int N10_thr = 7; // threshold for an AFT (hits)
		int twin = 200; // window size (ns)
		int aft_gate = 1300; // After trigger gate (ns)
		
		// CombinedFitter::FindAFTs();
		{
		// Search for the maximum number of hits in a 200 ns window
		// Return the maximum n200 and timing of the peak
		// CombinedFitter::GetMaxN200();

		// Do a neutron search - consider changing the method?????
		// Find N10 peaks and calculate n10, neff, etc for each peak
		// Treat peaks separated by more than 20 ns as different peaks
		// CombinedFitter::FindNeutronCandidates(); // option for MC
		}
	}

	else
	{
		// Get the number of AFTs that are neutron candidates
		// CombinedFitter::FindNeutronCandidates(); //option for data
	}

	// Now we have the number of AFTs for data and MC, do the combined fit

	// Use the goodness from the combined fit to select the best AFT
	
	for (int AFT = 0; AFT<1; AFT++)//m_data->HasAFT()+1; AFT++)
	{
		/*
		// TODO get the time sorted out.
		// This may be working for the prompt but not for the trigger. 
		// What are all of the different time data stored?
		// Trigger time? 
		int t_hour = skhead_.ntimsk[0];
		int t_min = skhead_.ntimsk[1];
		float t_sec = skhead_.ntimsk[2]+skhead_.ntimsk[3]/100.;

		int sub_48bitclock[3];
		// time difference for sub-events
		sub_48bitclock[0] = skhead_.nt48sk[0]+t_hour;
		sub_48bitclock[1] = skhead_.nt48sk[1]+t_min;
		sub_48bitclock[2] = skhead_.nt48sk[2]+t_sec;
		
		float t_nsec = (sub_48bitclock[2]+sub_48bitclock[1]*6.5536+sub_48bitclock[0]*4.294967296)*2.0;
//    Log(toolName+" Start of run "+toString(t_hour)+" h "+toString(t_min)+" m "+toString(t_sec)+" s",v_message,verbosity);
		float subtdiff = (t_nsec - prev_t_nsec)/1000.;
		printf("subevent time difference = %6.2f\n",subtdiff);
		*/
	
		//--------------Get the MC information we want----------------------//
		// NB in SKG4 the mc info is stored in reverse order:
		// neutron, positron, antineutrino for IBD events

		if (MC)
		{
			// mc vertex, direction, momentum, energy variables
			float *posmc = skroot_lowe_.posmc;
			float (*dirmc)[3] = skroot_lowe_.dirmc;
			float *pabsmc = skroot_lowe_.pabsmc;
			float *energymc = skroot_lowe_.energymc; 
			float p_x = 9999;
			float p_y = 9999;
			float p_z = 9999;
			float p_abs = 9999;
			double restmass;
			int pdg;
			
			if (!AFT)
			// Get the MC info for the prompt event
			{

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
				if (p_abs>0)
				{
					dirmc[0][0] = p_x/p_abs;
					dirmc[0][1] = p_y/p_abs;
					dirmc[0][2] = p_z/p_abs;
				}

				// Accessing the energy is different in skdetsim and SKG4
				if (mc->energy[0]!=0) 
				// In SKG4, the energy is already calculated
				{
					energymc[0] = mc->energy[1]; // positron energy
				}
				else 
				// In SKDetsim, mc energy is 0 so we have to calculate it
				{	
					//get the rest mass of the particle and calculate energy
					pdg = mc->ipvc[1]; // positron
					restmass = PdgToMass(pdg);
					energymc[0] = sqrt(pow(restmass,2)+pow(p_abs,2));
				}
			} // endif !aft
	
			else 
			// if AFT get the neutron MC info
			{

				// Here, the whole process is different for SKG4 and SKDetim
				if (mc->energy[0]!=0) 
				// SKG4 (energy is zero in SKDetsim)
				{
					// In SKG4, vertex and momentum are in secondaries branch
					SecondaryInfo* secondaries = mgr->GetSECONDARY();
					
					// Get the mc vertex of the neutron event
					posmc[0] = secondaries->vtxscnd[0][0];
					posmc[1] = secondaries->vtxscnd[0][1];
					posmc[2] = secondaries->vtxscnd[0][2];

					// Calculate the absolute momentum of the delayed
					p_x = mc->pvc[0][0];//x momentum
					p_y = mc->pvc[0][1];//y momentum
					p_z = mc->pvc[0][2];//z momentum
					p_abs = sqrt(pow(p_x,2)+pow(p_y,2)+pow(p_z,2));
					pabsmc[1] = p_abs;
					
					//calculate the mc direction
					if (p_abs>0)
					{
						dirmc[1][0] = p_x/p_abs;
						dirmc[1][1] = p_y/p_abs;
						dirmc[1][2] = p_z/p_abs;
					}
				
					//then get the energy of the neutron from the MC branch
					energymc[1] = mc->energy[0];
				}
				else
				// SKDetsim
				{
					// Get the mc x, y and z co-ordinates of the delayed
					posmc[0] = mc->pvtxvc[2][0];
					posmc[1] = mc->pvtxvc[2][1];
					posmc[2] = mc->pvtxvc[2][2];

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
				} // SKDetsim
			} // AFT	
		} // endif MC
		
		if (!AFT)
		{
			m_data->LoadSHE(readerName);
			evid = 0;
		}
		else if (AFT) 
		{ 
			m_data->LoadAFT(readerName);
			evid = 1;
		}
		else
		{
			Log(toolName+" not SHE or AFT",v_debug,verbosity);
		}
		//-------------------------------------------
		// Get the hit information in one of three ways
		//-------------------------------------------
		// TODO need to do something here to access the prompt and delayed hits separately.
		// In SKG4, all of the hits in a 535 us window are saved.
		// The prompt trigger is at t0_prompt = geant_t0 + MCInfo.prim_pret0
		// The MC time of the neutron capture is stored in SECONDARY.tscnd
		// but we don't really want to access hits using truth info.
		// Instead, we need to apply a secondary 'trigger' to find the hits from
		// the neutron capture. We can do this with a simple hit sum and look
		// for all secondary 'triggers' occurring after the prompt SHE trigger.
		
		float charges[skq_.nqisk];  	//final array of hit charges
		float times[skq_.nqisk];		//final array of hit times
		int cableIDs[skq_.nqisk];
		lbfset0(cableIDs,&skq_.nqisk);
		lbfset0(charges,&skq_.nqisk);
		lbfset0(times,&skq_.nqisk);
		int nhit=0;

		if (dataSrc == 0)
		{
			// Get the hit information stored in the common blocks (skt_/skq_)
			// ---------------------------------------------------
			nhit = skq_.nqisk;
			for (int i=0;i<skq_.nqisk; i++)
			{
				//skt_, skq_, skchnl_ common blocks
				cableIDs[i] = skchnl_.ihcab[i];
				charges[i] = skq_.qisk[skchnl_.ihcab[i]-1];
				times[i]   = skt_.tisk[skchnl_.ihcab[i]-1];
				if (i==skq_.nqisk-1)
				{ //time between first and last hit
					//std::cout << times[skq_.nqisk-1]-mc->prim << ", last time - first time = " << times[skq_.nqisk-1] - times[0] << std::endl;
				}
	//			std::cout << "cable " << skchnl_.ihcab[i] << "        " << charges[i] << "      " << times[i] << std::endl;
			}

		}

		else if (dataSrc==1)
		{
			// Get raw hit time and charge info from the common blocks (skt_/skq)
			nhit = sktqz_.nqiskz;
			for (int i=0;i<nhit; i++)
			{
				cableIDs[i] = sktqz_.icabiz[i];
				charges[i] = sktqz_.qiskz[i];
				times[i]   = sktqz_.tiskz[i];
				//std::cout << "cable " << cableIDs[i] << ": " << charges[i] << ", " << times[i] << std::endl;
			}
		}

		else 
		{
			//Get hit time and charge info from the tqreal branch
			//---------------------------------------------------
		
			/*  NB We can access hit time and charge info from the TQREAL branch
			*   but these are raw values and not corrected for 
			*   e.g. out-of-gate hits and bad channels .
			*   TODO Here, I have tried to go through the process
			*   which is usually carried out by tqrealsk ($SKOFL_ROOT/src/skrd/tqrealsk.F)
			*	when skread is called.
			*/
	//		myTreeReader->Get("TQREAL", tqreal);
			// 'raw' in the following means pre-skread i.e. info for all hits
			int nhits_raw = tqreal->nhits;
			Int_t *cableIDs_raw = new Int_t[MAXPM];//from TQREAL
			Float_t *charges_raw = new Float_t[MAXPM];//from TQREAL
			Float_t *times_raw = new Float_t[MAXPM];//from TQREAL
			Float_t *times_relative = new Float_t[MAXPM];//times relative to initial trigger
	   
			// TODO get some info about the trigger timing. Not really working at the moment. 
			float sub_trigger_time = tqreal->it0xsk;// at the moment this is only looking at the first trigger time
			float initial_trigger_time = skheadqb_.it0sk;
			float count_per_nsec = COUNT_PER_NSEC;
			//int ntrigsk = (sub_trigger_time-initial_trigger_time)/count_per_nsec/10.; //number of triggers - yields 0
			for (int i = 0;i<nhits_raw;i++){
				
				// first get the raw values from the tqreal branch
				// NB these are raw in the sense that they are for 
				// all of the hits, although in data some pre-processing
				// has been done by this stage e.g. to convert charge to pe.(by skrawread)
				--------------------------------------------------------------------------
				// get the charge and cable (PMT no.) data
				charges_raw[i]=tqreal->Q[i];
				cableIDs_raw[i] = tqreal->cables[i] & ( (1<<16)-1);//cable ID (i.e. PMT number) is stored in 1st 16 bits
				// then get the times relative to the trigger
				times_relative[i]=tqreal->T[i]-(sub_trigger_time-initial_trigger_time)/count_per_nsec;//(tqrealsk.F::122), count_per_nsec set in skhead.h;

				// then apply some selection to ensure the hits are in-gate (1.3us)
				// and in the inner detector	
				int in_gate = ( (tqreal->cables[i] >> 16) &1 ); // Upper 16 bits are the hit flags. We want the first of these.
				if (cableIDs_raw[i]>0 && cableIDs_raw[i]<=MAXPM && in_gate){
	//		   		TODO if BTEST( ISKBIT, 31-25).AND.IBAD(ICAB).NE.0) GOTO 120 ! skip badch if MC
	//				TODO mask bad channels
					//BTEST returns logical true if bit a POS in I is set: BTEST(I,POS)
					//std::bitset < sizeof(int)*32 > iskbit;
					//iskbit = skopt_.iskbit;
					//std::cout << "mc? " << iskbit.test(31-25) << std::endl;
					//if ( iskbit.test(31-25) ){
						// save the selected hits into arrays
						cableIDs[nhit]=cableIDs_raw[i];
						times[nhit] = times_relative[i];
						charges[nhit] = charges_raw[i];
						//std::cout << "cable " << cableIDs[nhit] << ": " << charges[i] << ", " << times[i] << std::endl;
						nhit++;
					//}
				}
			}
		}

		// Do the BONSAI fit
		//------------------
		// do the fit
		Float_t *bsvertex = skroot_lowe_.bsvertex;
		Float_t *bsresult = skroot_lowe_.bsresult;
		Float_t *bsgood = skroot_lowe_.bsgood;
		int nsel;
	//	Int_t bsnsel[2];
	  //float  bsdirks;
	  //float  bseffhit[12];
	  //float  bsenergy;
	  int    bsn50 = 9999;
	  //float  clvertex[4];
	  //float  clresult[4];
	  //float  cldir[3];
	  //float  clgoodness;
	  //float  cldirks;
	  //float  cleffhit[12];
	  //float  clenergy;
	  //int    cln50;
	  //int    islekeep;
	  //float  bspatlik;
	  //float  clpatlik;
	  //float  lwatert;
	  //int    lninfo;
	  //int32_t linfo[255];


		if (bonsaiSrc==0)
		{
			int* nhits = &skq_.nqisk;
			if (nhit>9 && nhit<=NHITCUT)
			{
				int nbf = bonsaifit_(&bsvertex[0],&bsresult[0],&bsgood[0],&nsel,&skq_.nqisk,&cableIDs[0],&times[0],&charges[0]);
			}
			
		}
		else 
		{
			if (nhit>NHITCUT)
			{
				std::cout << "Event " << ev << ", " << nhit << " nhits exceeded max." << std::endl;
				return false;
			}
			goodness *bshits = new goodness(bslike->sets(),bslike->chargebins(),bsgeom,nhit,cableIDs,times,charges);
			int nsel = bshits->nselected();
			if (bshits->nselected()<4) { 
				std::cout << "Event " << ev << ", " << bshits->nselected() << " selected hits not enough." << std::endl;
				return false;
			}
				fourhitgrid* bsgrid = new fourhitgrid(bsgeom->cylinder_radius(),bsgeom->cylinder_height(),bshits);

				bslike->set_hits(bshits);
				bslike->maximize(bsfit,bsgrid);
				int nbf = bslike->nfit();
				if (nbf >0){
		
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
					*bsgood = bsfit->maxq();								//max charge???
					bsfit->fitresult();										
					bslike->get_dir(bsresult);								//direction
					bsresult[5] = bslike->get_ll0();						//maximal log likelihood
					bslike->set_hits(NULL);				
				
					//std::cout << bsvertex[0] << ", " << bsvertex[1] << ", " << bsvertex[2] << ", " << bsvertex[3] << std::endl;
				}
	/*		else{
				std::cout << "not enough hits to fit" << std::endl;
			}*/
		}

		// get the other values based on the reconstructed vertex
		// TODO can we avoid using the fortran routines?

		// Calculate n50 (hits in (-10,40) ns around the fitted time of emission
		// lfallfit:283 lfnhit(50.,bsvertex,bsn50,bshit) 
		// lfnhit(twindow,vertex,nhit,ihitcab) 
		// input: twindow (width of timing window), vertex[3] (x,y,z vertex)
		// output: nhit (maximal number of hits in timing window), ihitcab[nhit] (array with cable numbers of these hits)
		if (bsvertex[0]<9999)
		{
			  Log(toolName+" calculating NX",v_debug,verbosity);
			  CalculateNX(50,bsvertex,cableIDs,times,bsn50);
		}


		// pass reconstructed variables from skroot_lowe_ common block to skroot file
		// TODO or do this without the common block
		skroot_set_lowe_(&lun,                      &skroot_lowe_.bsvertex[0], &skroot_lowe_.bsresult[0],
						 &skroot_lowe_.bsdir[0],    &skroot_lowe_.bsgood[0],   &skroot_lowe_.bsdirks,
						 &skroot_lowe_.bseffhit[0], &skroot_lowe_.bsenergy,    &bsn50,
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
		mgr->fill_tree(); // crashing for some reason
		mgr->Clear(); // zero out structures for next entry ! don't forget !

		//prev_t_nsec = t_nsec; // TODO reinstate when time is sorted
	} //endif AFT
    return true;
}



bool CombinedFitter::Finalise(){
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

int CombinedFitter::CalculateNX(int timewindow,float* vertex, int cableIDs[], float times[],int &bsn50)
{
	// Based on Fortran routine lfnhit ($skoflroot/lowe/sklowe/lfnhit.F)
	float cns2cm =21.58333; // speed of light in medium
	// Find tof subtracted times for all hits at reconstructed vertex
	vector<float> tof;
	Log(toolName+" getting tof subtracted times for all hits",v_debug,verbosity);
	for (int hit=0; hit<skq_.nqisk; hit++)
	{
		   tof.push_back(times[hit]-sqrt(pow((vertex[0]-geopmt_.xyzpm[cableIDs[hit]-1][0]),2)+pow((vertex[1]-geopmt_.xyzpm[cableIDs[hit]-1][1]),2)+pow((vertex[2]-geopmt_.xyzpm[cableIDs[hit]-1][2]),2))/cns2cm);
	}

	Log(toolName+" sorting tof subtracted times",v_debug,verbosity);
	// Sort in ascending order
	sort(tof.begin(),tof.end());
		   
	// Find the centre of the distribution
	Log(toolName+" finding the centre of the tof subtracted times distribution",v_debug,verbosity);
	bsn50 = 0;

	if (skq_.nqisk>0)
	{
			bsn50 = 1;
			int i = 0;
			while (i < skq_.nqisk-bsn50-1)
			{
					int j = i + bsn50;
					while ((j < skq_.nqisk) && (tof[j]-tof[i] <= timewindow))
					{   
							bsn50++;
							j++;
					}
			i++;
			}

	}

	else bsn50=0;

//	for (int i=0; i<skq_.nqisk; i++)
//	{
//		float dt = tof[i]-tof[skroot_lowe_.bsvertex[3]];
//		if ( (dt>=-timewindow/2.) && (dt<=timewindow/2.))
//		{
//			 bsn50++;
//		}
//	}
	return 0;
}

