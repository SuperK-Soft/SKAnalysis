/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef SimplifyTree_H
#define SimplifyTree_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.
#include "MTreeReader.h"

/**
* \class SimplifyTree
* A simple Tool to read a raw input file and convert it to a simple TTree that requires no libraries to read.
*
* $Author: M.O'Flaherty $
* $Date: 2022/01/15 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/

class SimplifyTree: public Tool {
	
	public:
	SimplifyTree();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	private:
	// functions
	// =========
	bool GetBranchValues();
	
	// tool variables
	// ==============
	std::string toolName;
	ConnectionTable* myConnectionTable=nullptr;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int iexecute=0;
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// variables to read in
	// ====================
	MTreeReader* iTreeReader = nullptr;
	int LUN;
	MTreeReader oTreeReader;
	std::string infilename = "";
	const Header* myHeader=nullptr;
	const TQReal* myTQReal=nullptr;
	const TQReal* myTQAReal=nullptr;
	int cableNumber;
	float charge;
	float time;
	long it0sk;
	long it0xsk;
	
	// variables to write out
	// ======================
	TFile* fout = nullptr;
	TTree* outtree = nullptr;
	std::string outputdir;
	float tubePosition[3];
	// simple output TTree branches
	int sk_phase;
	int run_num;
	int subrun_num;
	int event_num;
	std::string timestring;
	int trigger_flags;
	int event_flags;
	int readout_t0;
	int trigger_t0;
	int readout_len;
	std::vector<int> hit_id;        // PMT IDs
	std::vector<double> hit_q;      // charges
	std::vector<double> hit_t;      // times
	std::vector<int> hit_flags;    // hit flags
	std::vector<double> hit_x;      // hit x
	std::vector<double> hit_y;      // hit y
	std::vector<double> hit_z;      // hit z
	std::vector<double> hit_theta;  // hit theta
	std::vector<double> hit_loc;    // hit loc (kIDTop, kIDWall, kIDBot... from ConnectionTable.cc)
	float tubeR;
	float tubetheta;
	
};


#endif
