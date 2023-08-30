/* vim:set noexpandtab tabstop=4 wrap */
#include "vectgen.h"

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"
#include "TRandom3.h"
#include <cmath> // floor

vectgen::vectgen():Tool(){}

bool vectgen::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	m_variables.Get("verbosity", m_verbose);           // how verbose to be
	m_variables.Get("output_file", output_file);       // location of output file.
	m_variables.Get("output_format", output_format);   // "zbs" or "DatTable" (SKG4)
	m_variables.Get("num_events", total_events);       // we may directly specify number of events...
	m_variables.Get("event_rate", event_rate);         // ... or calculate it from an event rate... [evts/minute]
	m_variables.Get("run_number", run_number);         // ... and a given run's duration.
	m_variables.Get("livetime_file", livetime_file);   // list of runs and their livetimes
	m_variables.Get("random_seed", seed);              // for random generator used in generating events
	m_variables.Get("min_energy", min_E);              // minimum neutrino energy (MeV)
	m_variables.Get("max_energy", max_E);              // maximum neutrino energy (MeV)
	int sk_geometry=6;
	m_variables.Get("sk_geometry",sk_geometry);        // not sure if this is used. SK Phase.
	get_ok =  m_variables.Get("pos_x", pos_x);         // position of the vertex
	get_ok &= m_variables.Get("pos_y", pos_y);         // if not specified, will be random
	get_ok &= m_variables.Get("pos_z", pos_z);
	if(get_ok){
		position[0][0] = pos_x;
		position[0][1] = pos_y;
		position[0][2] = pos_z;
	} else {
		random_position = 1;
		position[0][0] = -1000; // not sure if this is needed? Why is this done? (vectgen.F line 80)
	}
	
	// validity check
	if(output_format!="zbs" && output_format!="DatTable"){
		Log(m_unique_name+" Uknown output format "+output_format+"! Valid options are "
		    "'zbs' or 'DatTable'",v_error,m_verbose);
		return false;
	}
	
	// initialize the random number generator, using the fortran routine rluxgo
	if(seed==0){
		// for reference, sonia's genrand.sh used:
		// python -c "from random import SystemRandom; s = SystemRandom(); print(s.randrange(1000000000))"
		TRandom3 rand;
		seed = rand.Uniform(0, 1000000000);
	}
	// note the following just sets the seed, so we get no return.
	int zero = 0;
	int lux = 3;   // 3/4-star luxury. better than average random numbers, but not the best.
	rluxgo_(&lux,&seed,&zero,&zero);
	
	// set the geometry. Not sure where this is used...
	m_data->GeoSet(sk_geometry);
	
	// if user didn't specify the number of events to generate,
	// scan livetime file to calculate number of events based on duration of given run
	float total_livetime=0;
	if(total_events==0){
		if(run_number==0 || event_rate==0){
			Log(m_unique_name+": Please provide either a number of events to generate, "
			   "or an event rate and run number for livetime.",v_error,m_verbose);
			return false;
		}
		if(livetime_file==""){
			Log(m_unique_name+" Warning: No livetime file given - using SK-IV livetime file",v_warning,m_verbose);
			livetime_file = "/home/sklowe/realtime_sk4_rep/solar_apr19/timevent/timevent.r061525.r077958";
		}
		// calculate livetime, noting subrun start date+times and durations (for zbs output)
		total_livetime = CalculateLivetime(run_number, livetime_file, &date, &time, &duration);
		
		// loop over subruns. We do this because zbs file records run/subrun info in LOWMC bank
		// so we need to periodically update it as we generate events
		for(int srun = 0; srun < duration.size(); ++srun){
			// total livetime is in seconds, event rate is in minutes, so:
			float exact_num_events = event_rate*(duration.at(srun) / 60.f) + 1.f;
			int num_events = floor(exact_num_events);
			// to be very precise, since the number of events predicted based on livetime and rate
			// may be fractional, but the number of events generated is an integer, we account
			// for the fractional part by throwing a random number to see whether we'll generate
			// an additional event in the remaining partial minute or not.
			float remainder = ( exact_num_events - num_events );
			float rnd = rlu_();  // $SKOFL_ROOT/src/monlib/rlu.F
			if(rnd > remainder) ++num_events; // oh go on then, just one more...
			eventcounts.push_back(num_events);
			total_events += num_events;
		}
	} else {
		// else user specified number of events to generate.
		eventcounts.push_back(total_events);
		date.push_back(std::array<int,3>{0,0,0});
		time.push_back(std::array<int,3>{0,0,0});
		duration.push_back(total_events / event_rate);
	}
	
	logmessage = m_unique_name+" Will generate "+toString(total_events)+" IBD events "
	    "with positron energy in the range "+toString(min_E)+" -- "+toString(max_E)+" MeV";
	if(random_position) logmessage += " with random positions within the fiducial volume";
	else logmessage += " at position ("+toString(pos_x)+", "+toString(pos_y)+", "+toString(pos_z)+")";
	logmessage += " with random directions";
	Log(logmessage,v_debug,m_verbose);
	
	// initialize zbs. This is needed whether we're writing to zbs or not.
	skheadf_.sk_file_format = 0;    // set common block variable for ZBS format
	m_data->KZInit();
	
	// open the output file
	if(output_format=="zbs"){
		if(output_file=="vectgen_out") output_file += ".zbs"; // add an extension if using default filename
		Log(m_unique_name+" creating output ZBS file "+output_file,v_debug,m_verbose);
		// get a unique handle for the file
		zbs_LUN = m_data->GetNextLUN();
		// open output file
		set_rflist_zbs( zbs_LUN, output_file.c_str(), true );
		int ipt = 1;
		int ihndl=1;
		skopenf_( &zbs_LUN, &ipt, "Z",&get_ok, &ihndl );
		if(get_ok!=0){   // based on skdump_new, this should be 0 for success...
			Log(m_unique_name+" Error opening output ZBS file!",v_error,m_verbose);
		}
	} else {
		// output fstream
		if(output_file=="vectgen_out") output_file += ".dat";
		Log(m_unique_name+" creating output dat file "+output_file,v_debug,m_verbose);
		datout = new std::ofstream(output_file.c_str(),std::ios::out);
		if(not datout->is_open()){
			Log(m_unique_name+" Error! Could not open output file "+output_file+" for writing!"
			    +" Does the path exist?",v_debug,m_verbose);
			return false;
		}
		// write out preamble. The DatTable is only read after a line containing "#DATASTART" is found.
		// so before that, we can record anything we like.
		*datout << "DatTable generated with vectgen ToolFramework tool from git commit "
		       << "XXX" << "and random seed "<< seed << "\n";  // TODO retreive git commit
		*datout << "Run number: "<<run_number<<", livetime "<<total_livetime<<" seconds, "
		       << "event rate "<<event_rate<<" evts/min, corresponding to "
		       <<total_events <<" events were requested\n";
		// TODO record further info such as cross-section model, spectrum weighting if applicable etc.
		// for reference, note the following commands will need to be called (or put in the macro file)
		// to read the generated DatTable file as a source of primary particles
		*datout << "Please add the following to your SKG4 macro file to use this DatTable file\n";
		*datout << "/SKG4/Primary/ParticleType DatTable\n";
		*datout << "/SKG4/Primary/DatTable/TablePath "<<output_file<<"\n";
		*datout << "/SKG4/Primary/DatTable/TableType 4\n"; // 4 = fully specified info for all primaries
		*datout << "/SKG4/Primary/DatTable/TableNumParticle 3\n"; // each event has 3 primaries
		*datout << "the following line signals the start of SKG4 input data\n\n";
		*datout << "#DATASTART" << std::endl;
	}
	
	// notes:
	// positron energy & neutrino directions are random, then neutrino+neutron kinematics are
	// calculated from them, based on probabilities given by the Strumia-Vissani cross-section model.
	// positron energies are uniformly distributed over the specified range, so events should be
	// reweighted according to a suitable positron energy spectrum.
	// positions are thrown within the ID (possibly with a hard-coded offset from the wall, currently 0)
	// generating events over a range of positron energies of 1-90MeV, over the full SK-IV
	// run range (61525 -- 77958) and with a rate of 1 event per minute takes only seconds to do,
	// and should generate ~4 million events.
	
	return true;
}

bool vectgen::Execute(){
	
	// Generate next ibd event
	seed = event_num + 3; // give each call a new random seed... offset by 3 to match i1 in vectgen.F
	// spectrum_ generates
	int flag = 9; // not sure this is actually used.
//	std::cout<<"calling spectrum with subrun "<<subrun<<", date "<<date.at(subrun).size()
//	         <<", time "<<time.at(subrun).size()<<", position "<<position[0][0]
//	         <<", seed "<<seed<<", E range "<<min_E<<"--"<<max_E
//	         <<", and random position setting "<<random_position<<std::endl;
	// outputs
	spectrum_(&flag,date.at(subrun).data(),time.at(subrun).data(),&position[0][0],
	          &seed,&min_E,&max_E,&random_position);
	if(output_format=="zbs"){
		Log(m_unique_name+": writing primaries to zbs output",v_debug,m_verbose);
		// slmcmklow - $SKOFL_ROOT/lowe/sollib/slmcmklow.F - populates the LOWMC zbs bank
		// which contains the passed variables on run num, start date+time and duration.
		slmcmklow_(&run_number,&subrun,date.at(subrun).data(),time.at(subrun).data(),&duration.at(subrun));
		// write data banks to zbs file
		kzwrit_(zbs_LUN);
		// reset common blocks for next event
		kzeclr_();
	} else if(output_format=="DatTable"){
		Log(m_unique_name+": writing output primaries to DatTable file",v_debug,m_verbose);
		// DatTable doesn't support run information variables,
		// but we do need to write the generated IBD particles
		// from the VCWORK and VCVRTX banks to the output file.
		// VCVRTX is the primary vertices - we should only have one, for an isolated IBD event.
		// VCWORK is the primary particles - we should have 3; neutrino, positron and neutron.
		// (the neutrino shouldn't really be recorded as an output particle from the primary vertex,
		// but it doesn't impact the detector simulation and this is somewhere to keep its info)
		
		// VCVRTX - Primary Particle Bank
		// ==============================
		// sanity check - we expect 3 output particles per IBD event
		if(vcwork_.nvc!=3){
			std::cerr<<"We have "<<vcwork_.nvc<<" primary particles!"<<std::endl;
			m_verbose=10;
		}
		// the DatTable file needs the pdg code, position, time, kinetic energy and direction of each primary
		// loop over output particles
		for(int i=0; i<vcwork_.nvc; ++i){
			// calculate the energy from the momentum and mass.
			float momentum_squared =  // TODO check index order
				pow(vcwork_.pvc[i][0],2.f) + pow(vcwork_.pvc[i][1],2.f) + pow(vcwork_.pvc[i][2],2.f);
			float rest_mass = vcwork_.amasvc[i];
			float kinetic_energy = sqrt(momentum_squared + pow(rest_mass,2.f)) - rest_mass;
			// calculate momentum magnitude to normalize it to form the unit direction
			float momentum_mag = sqrt(momentum_squared);
			
			// output the data to the DatTable
			*datout << i << " "                             // primary index in this event
			       << vcwork_.ipvc[i] << " "                // pdg code
			       << kinetic_energy  << " "                // KE [MeV]
			       << vcwork_.posivc[i][0]*10.f << " "      // pos_x [mm]  TODO check this is populated
			       << vcwork_.posivc[i][1]*10.f << " "      // pos_y [mm]  or vcwork_.posvc or vcvrtx_.pvtxvc?
			       << vcwork_.posivc[i][2]*10.f << " "      // pos_z [mm]  TODO check index order
			       << vcwork_.pvc[i][0]/momentum_mag << " " // mom_x [unit vector]
			       << vcwork_.pvc[i][1]/momentum_mag << " " // mom_x [unit vector] TODO check index order
			       << vcwork_.pvc[i][2]/momentum_mag << " " // mom_x [unit vector]
			       << vcwork_.timvc[i] << " \n";            // creation time [ns]
			
		}
	}  // end DatTable output
	
	++event_num;
	// check if this should be our last event
	if(event_num==total_events){
		Log(m_unique_name+" Finished generating all events, terminating loop",v_debug,m_verbose);
		m_data->vars.Set("StopLoop",1);
	}
	// check if we're moving to a new subrun
	if(event_num==eventcounts.at(subrun)) ++subrun;
	
	return true;
}


bool vectgen::Finalise(){
	
	if(datout){
		datout->close();
		delete datout;
		datout=nullptr;
	}
	if(zbs_LUN>0){
		// close the file
		skclosef_(&zbs_LUN);  // $SKOFL_ROOT/src/iolib/skclosef.f
		// inform DataModel the LUN is again free
		m_data->FreeLUN(zbs_LUN);
	}
	
	return true;
}


// TODO move to Algorithms or something more general
float vectgen::CalculateLivetime(int run_number, std::string livetime_filename, std::vector<std::array<int,3>>* start_date, std::vector<std::array<int,3>>* start_time, std::vector<float>* srun_duration){
	// open the file of livetimes - 3 steps as usual
	// 1. acquire a LUN
	int LUN = m_data->GetNextLUN();   // get a unique handle for the file
	// 2. set rflist
//	set_rflist_zbs(LUN,livetime_filename.c_str(),false);
	// TODO currently the set_rflist_zbs helper hard-codes a format (s6) of 'recl=5670' (zbs file...?)
	// whereas we need to use a format of 'form=unformatted', so just do it manually:
	set_rflist_(&LUN, livetime_filename.c_str(),"LOCAL","","RED","","","form=unformatted status=old","","",
	                 livetime_filename.length(),  5,    0,  3,   0, 0,           27,                 0,0);
	// 3. open the file (it's a fortran binary file, fyi)
	int fileIndex = 1;
	int ihndl=1;
	skopenf_(&LUN, &fileIndex, "f", &get_ok, &ihndl);
	
	// ok ready to read the file contents.
	// define variables to store line contents into
	int runnum=run_number;
	int subrun=0;
	int runstartdate[3];  // yy,mm,dd
	int runstarttime[3];  // hh,mm,ss
	float livetime=0;     // duration of run in seconds
	float solardir[3];    // cosine of direction to sun for this run
	int istat=0;          // status of next file line read (EOF = -1)
	int is_bad;
	// carry over runstartdate and time between loops to check for duplicates
	int prvrunstartdate[3];
	int prvrunstarttime[3];
	bool duplicate;
	// accumulate livetime over all subruns
	float total_livetime = 0;
	int line_num=0;
	// initialize slredtimev, may not be strictly necessary.
	int minus1=-1;
	slredtimev_(&minus1, &runnum, &subrun, &runstartdate[0], &runstarttime[0], &livetime, &solardir[0], &istat);
	// scan the livetime file for entries matching the requested run
	while(istat==0){
		if(m_verbose>2 && (line_num%1000)==0) printf(".");
		// $SKOFL_ROOT/lowe/sollib/slredtimev.f
		slredtimev_(&LUN, &runnum, &subrun, &runstartdate[0], &runstarttime[0], &livetime, &solardir[0], &istat);
		++line_num;
		if(istat < 0) break;                   // end of file
		if(runnum<60000) runnum += 65536;      // "for SK-IV analysis" it says..??
		if(runnum > run_number) break;         // past the desired run, further entries are of no interest to us.
		if(runnum!=run_number) continue;       // not our desired run; we don't care about it.
		// check if this is a bad run
		is_bad = lfbadrun_(&runnum, &subrun);  // $SKOFL_ROOT/lowe/sklowe/lfbadrun.F
		if(is_bad!=0) continue;                // will not contribute to analysis
		// vectgen.F also checks for multiple entries with the same start date+time; all duplicates are skipped
		duplicate = true;
		for(int i=0; i<3; ++i){
			duplicate &= ((runstartdate[i]==prvrunstartdate[i]) && (runstarttime[i]==prvrunstarttime[i]));
		}
		if(duplicate) continue;                // skip duplicates (as per vectgen.F, we do use the first entry)
		
		// All checks passed! Increment our livetime by the duration of this subrun.
		total_livetime += livetime;
		// copy out start date, time and duration of each subrun if requested
		if(start_date)    start_date->push_back({runstartdate[0],runstartdate[1],runstartdate[2]});
		if(start_time)    start_time->push_back({runstarttime[0],runstarttime[1],runstarttime[2]});
		if(srun_duration) srun_duration->push_back(livetime);
	}
	// close the file
	cclose_(&LUN);  // $SKOFL_ROOT/src/iolib/cclose.f
	// inform DataModel the LUN is again free
	m_data->FreeLUN(LUN);
	Log(m_unique_name+" livetime for run "+toString(run_number)+" was "
	    +toString(total_livetime)+" seconds",v_debug,m_verbose);
	
	return total_livetime;
}


