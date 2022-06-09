/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef TreeReader_H
#define TreeReader_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "SkrootHeaders.h" // MCInfo, Header etc.
#include "Constants.h"

#include "fortran_routines.h"

/**
* \class TreeReader
*
* FIXME
*
* $Author: M.O'Flaherty $
* $Date: 2019/05/28 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class TreeReader: public Tool {
	
	public:
	TreeReader();         ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute();   ///< Execute function used to perform Tool purpose.
	bool Finalise();  ///< Finalise funciton used to clean up resources.
	
	// we need to provide access to these functions via the DataModel...
	bool HasAFT();
	bool LoadAFT();
	bool LoadSHE();
	
	private:
	// functions
	// =========
	int ReadEntry(long entry_number, bool use_buffered=false);
	int AFTRead(long entry_number);
	int CheckForAFTROOT(long entry_number);
	int CheckForAFTZebra(long entry_number);
	int LoadAFTROOT();
	int LoadAFTZebra();
	void PrintTriggerBits();
	int LoadConfig(std::string configfile);
	void CloseLUN();
	const Header* myHeader=nullptr;
	void PrintSubTriggers();
	
	// tool variables
	// ==============
	std::string toolName;
	std::string inputFile="";
	std::string FileListName="InputFileList";
	std::string selectionsFile="";
	std::string cutName="";
	std::string treeName;
	std::string readerName;
	int maxEntries=-1;
	int firstEntry=0;
	int entrynum=0;
	int readEntries=0;                // count how many entries we've actually returned
	SKROOTMODE skrootMode=SKROOTMODE::READ;  // default to read
	int skreadMode=0;                 // 0=skread only, 1=skrawread only, 2=both
	int skreadUser=0;                 // 0=auto, 1=skread only, 2=skrawread only, 3=both
	int LUN=10;                       // Assumed 10 by some SK routines, only change if you know what you're doing!
	std::string skroot_options="31";  // 31 = read HEADER (required).
	int skroot_badopt=23;             // 23 = LOWE default (mask bad chs, dead chs, noisy ID chs and OD chs)
	int skroot_badch_ref_run=0;       // reference run for bad channel list for e.g. MC.
	int sk_geometry=4;                // TODO increment the default to 6.
	std::string outputFile="";        // for when using SKROOT copy mode
	int skip_ped_evts = 1;            // automatically skip pedestal/status TTree entries
	bool loadSheAftPairs=false;       // should we load and buffer the AFT for an SHE event, if there is one?
	bool onlyPairs=false;             // should we only return pairs of SHE+AFT events
	int entriesPerExecute=1;          // alternatively, read and buffer N entries per Execute call
	std::vector<int> triggerMasks;    // trigger bits required to return an entry
	
	std::vector<std::string> list_of_files;
	
	bool has_aft=false;    // do we have an AFT event buffered that matches this SHE event
	bool aft_loaded=false; // is the AFT loaded into the common blocks at present
	long buffered_entry = -1;
	long file_cursor=0;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
	// variables to read in
	// ====================
	MTreeReader myTreeReader;                  // the TTree reader
	MTreeSelection* myTreeSelections=nullptr;  // a set of entries one or more cuts
	std::vector<std::string> ActiveInputBranches;
	std::vector<std::string> ActiveOutputBranches;
	
	// functions involved in buffering common blocks
	// to load SHE+AFT pairs together
	int PushCommons();
	int PopCommons();
	int FlushCommons();
	bool LoadCommons(int buffer_i);
	bool LoadNextZbsFile();
	
	// common blocks to buffer
	// =======================
	// TODO trim down this list, almost certainly many of these are either
	// not populated by skread/skrawread, or are not used by reconstruction algorithms
	// and therefore do not need to be buffered.
	// TODO right now these do not really need to be vectors, since we only ever fill
	// them with at most one entry. Generlizing for storing buffering many entries,
	// but ... simplify if this is not useful.
	
	// event header - run, event numbers, trigger info...
	std::vector<skhead_common> skhead_vec;
	std::vector<skheada_common> skheada_vec;
	std::vector<skheadg_common> skheadg_vec;
	std::vector<skheadf_common> skheadf_vec;
	std::vector<skheadc_common> skheadc_vec;
	std::vector<skheadqb_common> skheadqb_vec;
	
	// low-e event variables
	std::vector<skroot_lowe_common> skroot_lowe_vec;
	std::vector<skroot_mu_common> skroot_mu_vec;
	std::vector<skroot_sle_common> skroot_sle_vec;
	
	// commons containing arrays of T, Q, ICAB....
	std::vector<skq_common> skq_vec;
	std::vector<skqa_common> skqa_vec;
	std::vector<skt_common> skt_vec;
	std::vector<skta_common> skta_vec;
	std::vector<skchnl_common> skchnl_vec;
	std::vector<skthr_common> skthr_vec;
	std::vector<sktqz_common> sktqz_vec;
	std::vector<sktqaz_common> sktqaz_vec;
	std::vector<rawtqinfo_common> rawtqinfo_vec;
	
	std::vector<sktrighit_common> sktrighit_vec;
	std::vector<skqv_common> skqv_vec;
	std::vector<sktv_common> sktv_vec;
	std::vector<skchlv_common> skchlv_vec;
	std::vector<skthrv_common> skthrv_vec;
	std::vector<skhitv_common> skhitv_vec;
	std::vector<skpdstv_common> skpdstv_vec;
	std::vector<skatmv_common> skatmv_vec;
	
	// OD mask....? nhits, charge, flag...
	std::vector<odmaskflag_common> odmaskflag_vec;
	
	// hardware trigger variables; counters, trigger words, prevt0...
	// spacer and trigger info.
	std::vector<skdbstat_common> skdbstat_vec;
	std::vector<skqbstat_common> skqbstat_vec;
	std::vector<skspacer_common> skspacer_vec;
	
	// gps word and time.
	std::vector<skgps_common> skgps_vec;
	std::vector<t2kgps_common> t2kgps_vec;
	
	// hw counter difference to previous event.
	std::vector<prevt0_common> prevt0_vec;
	std::vector<tdiff_common> tdiff_vec;
	std::vector<mintdiff_common> mintdiff_vec;
	
	// trigger hardware counters, word, spacer length...
	std::vector<sktrg_common> sktrg_vec;
	
	// MC particles and vertices, event-wise.
	std::vector<vcvrtx_common> vcvrtx_vec;
	std::vector<vcwork_common> vcwork_vec;
	
};


#endif
