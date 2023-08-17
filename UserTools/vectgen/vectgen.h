/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef vectgen_H
#define vectgen_H

#include <string>
#include <iostream>
#include <fstream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.

/**
* \class vectgen
* Generate the output particles from inverse beta decay events. based on the old vectgen fortran file.
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class vectgen: public Tool {
	
	public:
	vectgen();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	float CalculateLivetime(int runnum, std::string livetime_filename, std::vector<std::array<int,3>>* start_date=nullptr, std::vector<std::array<int,3>>* start_time=nullptr, std::vector<float>* srun_duration=nullptr);
	
	// tool variables
	// ==============
	std::string output_file = "vectgen_out";
	std::string output_format = "DatTable"; // or can be 'zbs'
	int num_events = 0;         // num IBD events to generate. if !=0, overrides use of run number + rate
	float event_rate = 0;       // events per minute
	int run_number = 0;         // run to use to look up livetime
	std::string livetime_file = ""; // list of runs and their run durations. binary fortran file.
	int seed = 0;               // if 0, will be generated on the fly.
	float min_E = 1;            // energy range of *positrons* [MeV]
	float max_E = 90;           // 
	float pos_x = 0;            // we may specify a fixed position
	float pos_y = 0;
	float pos_z = 0;
	int random_position = 0;
	int zbs_LUN = 0;            // LUN for output zbs file, if writing to zbs
	std::ofstream* datout = nullptr;  // ofstream for output if writing to DatTable
	std::vector<std::array<int,3>> date;  // date and ...
	std::vector<std::array<int,3>> time;  // time and ...
	std::vector<float> duration;          // duration of matching subruns
	std::vector<int> eventcounts;         // number of events for this subrun, based on duration and event rate
	int total_events = 0;                 // total events in this run, over all subruns
	float position[3][10];      // passed to subroutines that use this variable to hold position of all primaries
	
	int subrun = 0;             // counter for Execute loop.
	int event_num = 0;          // counter for Execute loop.
	
};


#endif
