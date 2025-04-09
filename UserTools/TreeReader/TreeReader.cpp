/* vim:set noexpandtab tabstop=4 wrap */
#include "TreeReader.h"
#include "TTree.h"
#include <set>
#include <bitset>
#include <algorithm> // std::reverse
#include <wordexp.h>  // wordexp

#include "Algorithms.h"
#include "Constants.h"
#include "type_name_as_string.h"
#include "MTreeSelection.h"
#include "TreeManagerMod.h"
#include "SuperWrapper.h"
#include "fortran_routines.h"
#include "SK_helper_functions.h"
#include <cstdlib> // atoi

TreeReader::TreeReader():Tool(){}

const std::vector<std::string> default_branches{
	"HEADER",
	"TQREAL",
	"TQAREAL",
	"LOWE",
	"ATMPD",
	"UPMU",
	"MU",
	"SLE",
	"SWTRGLIST",
	"MC",
	"SECONDARY",
	"TQLIST",
	"ODTQLIST",
	"HWTRGLIST",
	"PEDESTALS",
	"EVENTHEADER",
	"EVENTTRAILER",
	"SOFTWARETRG",
	"QBEESTATUS",
	"DBSTATUS",
	"SPACERS",
	"PREVT0",
	"MISMATCHEDHITS",
	"GPSLIST",
	"T2KGPSLIST",
	"IDODXTLK",
	"SN"
};

bool TreeReader::Initialise(std::string configfile, DataModel &data){
	
	m_data= &data;
	
	Log(m_unique_name+": Initializing",v_debug,m_verbose);
	
	// Get the Tool configuration variables
	// ------------------------------------
	LoadConfig(configfile);
	m_unique_name = "TreeReader "+readerName;
	m_data->tool_configs[m_unique_name] = &m_variables;
	myTreeReader.SetName(readerName);
	
	// safety check that we were given an input file
	if(inputFile=="" && FileListName==""){
		// unless we are working in SKROOT write mode...
		if(skrootMode!=SKROOTMODE::WRITE){
			Log(m_unique_name+" error! no InputFile or FileListName given!",v_error,m_verbose);
			m_data->vars.Set("StopLoop",1);
			return false;
		}
	} else if(inputFile!=""){
		// single file takes precedence
		list_of_files.emplace_back(inputFile);
	} else {
		get_ok = m_data->CStore.Get(FileListName, list_of_files);
		if(!get_ok){
			Log(m_unique_name+" error! Could not find file list "+FileListName+" in CStore!"
				+" Ensure LoadFileList tool is run before this tool!",v_error,m_verbose);
			m_data->vars.Set("StopLoop",1);
			return false;
		}
	}
	// check if we're reading we have at least one input file
	if(skrootMode!=SKROOTMODE::WRITE && list_of_files.empty()){
		Log(m_unique_name+" error! no files given!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// detect zebra files based on extention of first file
	if(skrootMode!=SKROOTMODE::WRITE){
		std::string firstfile = list_of_files.front();
		if( (firstfile.substr(firstfile.length()-4,firstfile.length())==".zbs") ||
		    (firstfile.substr(firstfile.length()-4,firstfile.length())==".dat") ){
			Log(m_unique_name+" using zebra mode",v_debug,m_verbose);
			skrootMode=SKROOTMODE::ZEBRA;
		}
	}
	
	// safety check that the requested name to associate to this reader is free
	get_ok = m_data->Trees.count(readerName);
	if(get_ok){
		Log(m_unique_name+" error! TreeReader tool used to open file with name "+readerName
			+" but this name is already taken! Each name must be unique!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// safety check we have an output file too, if working in SKROOT write or copy mode
	if((skrootMode==SKROOTMODE::WRITE || skrootMode==SKROOTMODE::COPY) && outputFile==""){
		logmessage = m_unique_name+" error! SKROOT mode is ";
		logmessage += ((skrootMode==SKROOTMODE::WRITE) ? "write" : "copy");
		logmessage += " but no outputFile specified!";
		Log(logmessage,v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// safety check that if asked to read an SKROOT file, it has a TTree called 'data'
	// if not, the TreeManager will segfault!
	// while we have it, also make a note of the tree branches, as they're not standard
	std::vector<std::string> present_branches;
	if(skrootMode==SKROOTMODE::READ || skrootMode==SKROOTMODE::COPY){
		// i guess we can only pracitcally check the first file
		// i don't think it'll seg as long as at least one file has a 'data' tree
		// XXX although, perhaps it would be better to check all of them?
		bool got_good_file=false;
		for(int first_valid_file_i=0; first_valid_file_i<list_of_files.size(); ++first_valid_file_i){
			TFile* f_temp = TFile::Open(list_of_files.at(first_valid_file_i).c_str(),"READ");
			if(!f_temp || f_temp->IsZombie()){
				Log(m_unique_name+" ERROR opening input file '"+list_of_files.at(first_valid_file_i)
				    +"'",v_error,m_verbose);
				continue;
			}
			TTree* t_temp= (TTree*)f_temp->Get(treeName.c_str());
			if(t_temp==nullptr){
				Log(m_unique_name+" ERROR! input file "+list_of_files.at(first_valid_file_i)
				    +" has no '"+treeName+"' TTree!",v_error,m_verbose);
				f_temp->Close();
				// might as well remove it from the list of files
				list_of_files.erase(std::next(list_of_files.begin(),first_valid_file_i));
				continue;
			}
			got_good_file=true;
			present_branches.resize(t_temp->GetListOfBranches()->GetEntriesFast());
			for(int i=0; i<present_branches.size(); ++i){
				present_branches.at(i) = t_temp->GetListOfBranches()->At(i)->GetName();
			}
			f_temp->Close();
		}
		if(!got_good_file){
			Log(m_unique_name+" ERROR! Found no good files in file list!",v_error,m_verbose);
			m_data->vars.Set("StopLoop",1);
			return false;
		}
	}
	
	// warning check: see if we're given an input when we're in WRITE mode
	if(skrootMode==SKROOTMODE::WRITE && (inputFile!="" || FileListName!="")){
		Log(m_unique_name+" warning! InputFile or FileListName given, but mode is skroot::write! "
					+"Inputs will be ignored: use another reader instance to read input from another file",
					v_error,m_verbose);
	}
	// warning check: see if we've been given an output when we're not in SKROOT::COPY or WRITE mode
	if((skrootMode!=SKROOTMODE::WRITE && skrootMode!=SKROOTMODE::COPY) && outputFile!=""){
		logmessage  = m_unique_name+" warning! outputFile given, but SKROOT mode is ";
		logmessage += ((skrootMode==SKROOTMODE::NONE) ? "not enabled. " : "READ only. ");
		logmessage += "outputFile will be ignored!";
		Log(logmessage,v_warning,m_verbose);
	}
	
	// open the input TFile and TTree
	// ------------------------------
	// a lot of SK algorithms retrieve data via skroot_get_* calls behind the scenes,
	// which means if we're processing skroot files we probably need a corresponding TreeManager.
	if(skrootMode!=SKROOTMODE::NONE){
		
		// skroot_open_* invokes the singleton SuperManager class to create a TreeManager
		// to be associated with a given SKROOT file. This TreeManager is just the usual
		// TTree wrapper created via ROOT's MakeClass. 
		// The SuperManager keeps a map of LUNs (unique IDs) for each TreeManager.
		// Unfortunately there is no way to know if a LUN is already in use or not.
		// Calling GetManager will return a nullptr if it doesn't, but then skroot_open_*
		// will fail if you try to subsequently open a new file with this LUN. Doh!
		// TODO fix the SuperManager.
		// For now we'll just keep our own list of LUNs in the DataModel
		LUN = m_data->GetNextLUN(readerName, LUN);
		
		// "initialize data structures".
		// Allegedly a ZEBRA thing, but we seem to need this even for ROOT files??
		m_data->KZInit();
		
		// slight change in initialization depending on SK root vs zebra
		if(not (skrootMode==SKROOTMODE::ZEBRA)){
			Log(m_unique_name+": callling skroot_open_*",v_debug,m_verbose);
			
			// There are 3 modes to the TreeManager:
			// skroot_open_read_ calls the TreeManager constructor with mode = 2;
			// Some branches may be skipped to optimize reading by using ZeroInBranch.
			
			// skroot_open_write_ calls the TreeManager constructor with mode =1;
			// It creates an output TTree with the full set of SKROOT branches.
			// You should call skroot_set_* or TreeManager::Set* methods to populate branches, then
			// call skroot_fill_tree (or TreeManager::fill_tree) to write a new output tree entry.
			// Note: ZeroOutBranch doesn't work in this mode, all output branches will be created.
			// Note: This mode is commented as "zbs 2 root" (or sometimes "root 2 zbs"),
			// but neither skroot functions nor the TreeManager provide any means of handling zbs files.
			// To do that, see e.g. $SKOFL_ROOT/examples/lowe/zbs2skroot.F
			
			// skroot_open_ calls the TreeManager constructor with mode = 0;
			// This allows input file reading (as per mode 2), but also creates an output file
			// with an SKROOT tree set up for copying events from input to output.
			// (i.e. both trees' branch addresses point to the same object).
			// Input read can be optimized via ZeroInBranch - these will not be read in or copied across.
			// Further branches may be omitted from the copy to output via ZeroOutBranch.
			// Output entries are written on each TreeManager::fill_tree() call,
			// so a subset of entries may be copied across.
			
			// create the treemanager, and in write mode, the output file
			switch(skrootMode){
				case SKROOTMODE::READ:  skroot_open_read_(&LUN); break;
				case SKROOTMODE::WRITE: skroot_open_write_(&LUN, outputFile.c_str(), outputFile.size()); break;
				case SKROOTMODE::COPY:  skroot_open_(&LUN, outputFile.c_str(), outputFile.size()); break;
				default: break; /* no action, just to quiet compiler warnings */
			}
			
			// register any duplicate LUNs requested
			for(int& dup : duplicate_luns){
				AddDuplicateLun(LUN, dup);
			}
			
			// the following are not relevant for WRITE mode
			if(skrootMode!=SKROOTMODE::WRITE){
				Log(m_unique_name+": adding files to read...",v_debug,m_verbose);
				// set the input file(s) to read.
				// These need to be paths to files - glob patterns are not supported.
				// Only the first will be used to retrieve the RareList, but multiple calls can
				// be made to add subsequent files to the TChain.
				// This doesn't yet open the files, it just adds them to a list of strings.
				for(std::string& fname_in : list_of_files){
					/*
					printf("doing path expansion on path '%s'\n",fname_in.c_str());
					wordexp_t p;
					char** w;
					wordexp(fname_in.c_str(), &p, 0);
					w = p.we_wordv;
					for(size_t i=0; i<p.we_wordc; ++i) printf("***%s***\n",w[i]);
					wordfree(&p);
					printf("end of expansion\n");
					// seem to be that this will expand an array, delimited by shell IFS, if necessary
					if(p.we_wordc<1){
						Log(m_unique_name+" expansion of input file "+fname_in+" produced no output?!",
							v_error,m_verbose);
						continue;
					}
					// might be useful but we probably only need the first entry
					fname_in = w[0];
					printf("next SKROOT file:XXXX%sXXXXX\n",fname_in.c_str());
					*/
					// above doesn't seem to be necessary; ROOT does this automatically
					Log(m_unique_name+": next SKROOT file "+fname_in,v_debug,m_verbose);
					skroot_set_input_file_(&LUN, fname_in.c_str(), fname_in.size());
				}
				
				// disable unused input branches.
				int io_dir = 0; // 0 for input branch, 1 for output branch
				Log(m_unique_name+": disabling unused branches...",v_debug,m_verbose);
				Log(m_unique_name+": we have "+toString(SkippedInputBranches.size())+" user skipped branches"
				                  ", "+toString(ActiveInputBranches.size())+" user-selected active branches"
				                  ", and "+toString(present_branches.size())+"/"+toString(default_branches.size())
				                  +" standard branches present", v_debug,m_verbose);
				
				// ok first disable branches via a list of inactive branches specified by the user
				for(auto&& abranch : SkippedInputBranches){
					if(abranch=="HEADER") continue; // always required
					Log(m_unique_name+": disabling user-skipped input branch "+abranch,v_debug,m_verbose);
					skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
				}
				
				// alternatively allow users to specify just the branches they want to keep,
				// and skip everything else.
				if(ActiveInputBranches.size()){
					for(auto&& abranch : present_branches){
						if(std::find(ActiveInputBranches.begin(),ActiveInputBranches.end(),abranch) ==
							ActiveInputBranches.end()){
							if(abranch=="HEADER") continue; // always required
							Log(m_unique_name+": disabling user-disallowed input branch "+abranch,v_debug,m_verbose);
							skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
						}
					}
				}
				
				// lastly, the TreeManager attempts to read data from a standard set of branches,
				// and segfaults if they don't exist... SO let's explicitly disable everything standard
				// if it's not actually present
				for(auto&& abranch : default_branches){
					if(std::find(present_branches.begin(), present_branches.end(), abranch)==present_branches.end()){
						Log(m_unique_name+": disabling absent input branch "+abranch,v_debug,m_verbose);
						skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
					}
				}
				
				Log(m_unique_name+": done disabling unused branches",v_debug,m_verbose);
			}
			// disable unwanted output branches. Only applicable to COPY mode
			if(skrootMode==SKROOTMODE::COPY){
				int io_dir = 1;
				// first from list of disabled branches from user
				for(auto&& abranch : SkippedOutputBranches){
					if(abranch=="HEADER") continue; // always required
					Log(m_unique_name+": disabling user-skipped output branch "+abranch,v_debug,m_verbose);
					skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
				}
				// alternatively from a list of enabled branches from user
				if(ActiveOutputBranches.size()){
					for(auto&& abranch : default_branches){
						if(std::find(ActiveOutputBranches.begin(),ActiveOutputBranches.end(),abranch) ==
							ActiveOutputBranches.end()){
							if(abranch=="HEADER") continue; // always required
							Log(m_unique_name+": disabling user-disallowed output branch "+abranch,v_debug,m_verbose);
							skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
						}
					}
				}
				/*
				// and anything else which isn't present in the input file? XXX not sure this is needed
				for(auto&& abranch : default_branches){
					if(std::find(present_branches.begin(), present_branches.end(), abranch)==default_branches.end()){
						Log(m_unique_name+": disabling output branch "+abranch+" as it is absent from input Tree",v_debug,m_verbose);
						skroot_zero_branch_(&LUN, &io_dir, abranch.c_str(), abranch.size());
					}
				}
				*/
			}
			
			// ok now we perform the actual file opening, TTree cloning, and branch address setting.
			// except in the case of 'write', where all necessary steps are done on construction
			if(skrootMode!=SKROOTMODE::WRITE){
				if(treeName=="data"){
					// this mostly just invokes TreeManager::Initialise()
					skroot_init_(&LUN);
				} else {
					// TreeManager::Initialise hard-codes the name of the tree to 'data'.
					// but what if we don't want that? We can hack it just a little.
					TreeManager* mgr = skroot_get_mgr(&LUN);
					TreeManagerMod* mgrmod = (TreeManagerMod*)(mgr);
					mgrmod->Initialize(treeName);
					skheadf_.sk_file_format = 1;
					skheadf_.root_id = LUN;
				}
			} else {
				// if the user asked for a different tree name than 'data', apply that now...
				TreeManager* mgr = skroot_get_mgr(&LUN);
				TTree* otree = mgr->GetOTree();
				otree->SetName(treeName.c_str());
			}
			
		} else {
			// else input files are ZEBRA files
			Log(m_unique_name+" doing zebra initialization",v_debug,m_verbose);
			skheadf_.sk_file_format = 0;    // set common block variable for ZBS format
			
			//zbsinit_();  THIS IS THE SAME ROUTINE AS KZINIT!!
			
			// Set rflist and open file
			// '$SKOFL_ROOT/iolib/set_rflist.F' is used for setting the input file to open.
			// (See also $RELICWORKDIR/mc_sim/inject_dummy/src/set_rflistC.F)
			// This copies the input filename into the RFLIST common block.
			// Despite the name i can't see how 'rflist' actually supports a list,
			// (possibly as some environmental variable..????)
			// so we'll have to invoke skopenf for each file as we go until we have none left.
			// One way we can do this is reverse the list and use pop_back until we run out.
			std::reverse(list_of_files.begin(), list_of_files.end());
			
			// Load the first file
			get_ok = LoadNextZbsFile();
			if(get_ok==0){
				Log(m_unique_name+" Error loading first ZBS file!",v_error,m_verbose);
				return false;
			}
		}
//		std::cout<<"get_ok = "<<get_ok<<std::endl;
//		assert(false);
//		exit(-1);
		
		// we can access the manager via:
		//TreeManager* mgr = skroot_get_mgr(&LUN);
		
		// the primary reason to use the TreeManager is to populate the fortran common blocks
		// required by many old SK algorithms. This is not done by the TreeManager itself,
		// but by the fortran routines `skread` and `skrawread`. These internally read data
		// from the ROOT files via the skroot_get_* functions (e.g. via headsk.F),
		// and so depend on there being an underlying TreeManager.
		// There are a few extra steps we need to do to initialize up the fortran common blocks...
		
		// skoptn 25 indicates that skread should mask bad channels.
		// As a base there is a manually maintained list of bad channels, $SKOFL_ROOT/const/badch.dat
		// skoptn 26 indicates that the set of bad channels should be looked up using
		// the run and subrun numbers of the current event. In this case we can instead (in addition?)
		// use a list of bad channels generated during online processing, and saved to:
		// $SKOFL_ROOT/const/badch/badch.XXXX where XXXX is the run number.
		// For MC or calibration files, skoptn 26 should *not* be given, since the run number in MC data
		// is a dummy and the corresponding bad channel file will not be found.
		// Instead skbadch() may be called with a reference run to use for the bad channel list.
		// We could rely on the user passing the right options, but it would be good to verify this.
		// However, this is a bit tricky because we don't know if the input file is MC
		// until we've opened it and read some entries. On the other hand we normally set the skoptns
		// before reading from the file....
		// For now, we'll omit the run-wise bad channel masking for the MC check,
		// and then re-call skoptn_ with the additional 26 flag if it's required and suitable to do so.
		// This does rely on repeated calls to skoptn_ being handled correctly.
		bool use_runwise_bad_ch_masking = skroot_options.find("26")!=std::string::npos;
		if(use_runwise_bad_ch_masking) skroot_options.erase(skroot_options.find("26"),2);
		bool isMC=false;
		
		if(skrootMode!=SKROOTMODE::WRITE){
			// need to set skgeometry in skheadg common block
			m_data->GeoSet(sk_geometry);
			// note: it seems like we do need to call this before SKRAWREAD,
			// but if the geometry is incorrect, SKRAWREAD will reload the correct geometry. Maybe.
			
			// options for what to read etc.
			skoptn_(const_cast<char*>(skroot_options.c_str()), skroot_options.size());
			
			// options for masking OD / dead / noisy channels
			skbadopt_(&skroot_badopt);
		}
		
		// we can provide an MTreeReader except in SKROOT::WRITE mode or when reading zebra files
		if(skrootMode!=SKROOTMODE::ZEBRA && skrootMode!=SKROOTMODE::WRITE){
			// since an MTreeReader provides passive file access
			// it doesn't interfere with the use of a TreeManager
			Log(m_unique_name+" creating MTreeReader in parallel to TreeManager",v_debug,m_verbose);
			TTree* intree = skroot_get_tree(&LUN);  // n.b. no trailing underscore for this one
			get_ok = myTreeReader.Load(intree);
			if(not get_ok){
				Log(m_unique_name+" failed to make MTreeReader",v_error,m_verbose);
				return false;
			}
			if(mTreeReaderVerbosity) myTreeReader.SetVerbosity(mTreeReaderVerbosity);
			// closing and deleting the file will be done by the TreeManager
			myTreeReader.SetOwnsFile(false);
			
		}
		
		
		if(skrootMode!=SKROOTMODE::ZEBRA && skrootMode!=SKROOTMODE::WRITE){
			
			// skread/tqreal etc detect whether a file is MC or data by checking
			// the `mdrnsk` ("SK run mode" ...) member of the HEADER branch.
			// Unfortunately this is only populated in physics event entries
			// (not pedestal or status entries), which means we need to scan
			// the file for the first suitable entry before we can do the check.
			/*
			// The below works, but it spits out: "error reading the root tree"
			// until it finds a data event. Safe to ignore, but misleading...
			while(true){
				skcread_(&LUN, &get_ok); // get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
				std::cout<<"entry returned "<<get_ok<<std::endl;
				if(get_ok==0 || get_ok==2) break;
			}
			isMC = (skhead_.mdrnsk==0 || skhead_.mdrnsk==999999);
			*/
			// We could instead just copy the checks for PDST / RUNINFO entries
			// from the top of headsk.F and avoid using skcread.
			int tmp_entry=0;
			while(true){
				get_ok = myTreeReader.GetEntry(tmp_entry);
				if(get_ok<=0){
					Log(m_unique_name+" error! Hit end of tree while checking if MC!",v_error,m_verbose);
					return false;
				}
				const Header* header = nullptr;
				myTreeReader.Get("HEADER", header);
				get_ok  = ((header->nrunsk==0 && header->mdrnsk!=0 && header->mdrnsk!= 999999) ||
				           (std::bitset<8*sizeof(int)>(header->ifevsk).test(19)) ||
				           (header->nrunsk==0 && header->sk_geometry==0 && header->ifevsk!=0) );
				if(get_ok==0){
					// physics entry!
					isMC = (header->mdrnsk==0 || header->mdrnsk==999999);
					break;
				}
				++tmp_entry;
			}
			Log(m_unique_name+" read "+toString(tmp_entry+1)+" entries before determining that isMC = "
			    +toString(isMC),v_debug,m_verbose);
			
			/*
			// A much simpler option is to check for the 'MC' branch, but this could be:
			// 1) dropped during analysis, if the user decides it's not needed
			// 2) created but not filled, if the file is made with skroot_open_write
			// so that's not perfect either.
			isMC = (intree->FindBranch("MC")!=nullptr);
			*/
			
			// put this MC flag into the MTreeReader
			myTreeReader.SetMCFlag(isMC);
			m_data->vars.Set("inputIsMC",isMC);
			
			// XXX assume that if we're processing MC ROOT files, it's SKG4
			// set to 0 for skdetsim (according to lowfit_sk6_mc_badrun)
			wtpar_.flag_g4_wtpar = 1;
			
			// if we're processing MC, we should probably only call `skread`. If we're processing
			// data, we probably need to call `skrawread` as well. (see ReadEntry for more info).
			// SKRAWREAD also requires the presence of a TQLIST branch.
			// We'll use this as a default, but allow ourselves to be overridden by a user.
			if(skreadUser==0){
				bool no_tqlist = (std::find(present_branches.begin(), present_branches.end(), "TQLIST")==present_branches.end());
				 // skreadMode 0=skread only, 1=skrawread only, 2=both
				skreadMode = (isMC || no_tqlist) ? 0 : 2;
			} else {
				skreadMode = skreadUser-1;
			}
			
		} else if(skrootMode==SKROOTMODE::ZEBRA){
			
			// fall back to using skread to scan for MC if we're not reading ROOT files
			while(true){
				skcread_(&LUN, &get_ok); // get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
				Log(m_unique_name + " next isMC scan entry returned "+toString(get_ok),v_debug,m_verbose);
				// doesn't seem to be sufficient to identify non-physics events?
				// manually add in checks as per headsk.F for ROOT, these checks not done for ZBS though...?? XXX
				bool chk = ((skhead_.nrunsk==0 && skhead_.mdrnsk!=0 && skhead_.mdrnsk!= 999999) ||
				           (std::bitset<8*sizeof(int)>(skhead_.ifevsk).test(19)) ||
				           (skhead_.nrunsk==0 && skheadg_.sk_geometry==0 && skhead_.ifevsk!=0) );
				if(get_ok==0 && chk) get_ok = 3; // FIXME is this correct? necessary?
				if(get_ok==0 || get_ok==2) break;
			}
			if(get_ok==2){
				Log(m_unique_name+" hit end of file when doing MC detection!",v_error,m_verbose);
				return false;
			}
			// extract isMC
			isMC = (skhead_.mdrnsk==0 || skhead_.mdrnsk==999999);
			m_data->vars.Set("inputIsMC",isMC);
			
			// XXX assume that if we're processing MC ZBS files, it's skdetsim
			// set to 0 for skdetsim (according to lowfit_sk6_mc_badrun)
			wtpar_.flag_g4_wtpar = 0;
			
			if(skreadUser==0){
				skreadMode = (isMC) ? 0 : 2;  // 0=skread only, 1=skrawread only, 2=both
			} else {
				skreadMode = skreadUser-1;
			}
			
			// if user's intended mode is to invoke both skread and skrawread,
			// do the follow-up skrawread as well. Use negative LUN to prevent advance.
			if(skreadMode>0){
				Log(m_unique_name+" finishing first load with skrawread",v_debug,m_verbose);
				int LUN2 = -LUN;
				skcrawread_(&LUN2, &get_ok);
			}
		} // else SKROOT::Write mode, cannot check whether it's MC or not.
		
		// make a note of is_mc in m_variables so it gets saved in m_data->tool_configs
		m_variables.Set("is_mc",isMC);
		
		// now we know if our file is MC or not, set the appropriate runwise bad channel masking if required
		if(skrootMode!=SKROOTMODE::WRITE && skroot_options.find("25")!=std::string::npos && use_runwise_bad_ch_masking){
			if(not isMC){
				skroot_options = skroot_options + " 26";
				skoptn_(const_cast<char*>(skroot_options.c_str()), skroot_options.size());
			} else {
				if(skroot_ref_run==-1){
					Log(m_unique_name+" Error! skbadoptn contains 25 (mask bad channels) but not 26 "
						+"(i.e. look up bad channels based on run number). In this case one needs to provide "
						+"a reference run to use for the bad channel list! Please specify a run to use in "
						+"option 'mcReferenceRun' in "+m_unique_name+" config",v_error,m_verbose);
					return false;
				}
				Log(m_unique_name+" masking bad channels with reference run "
					+toString(skroot_ref_run),v_debug,m_verbose);
				combad_.log_level_skbadch = 0;  // n.b. loglevel of 4 suppresses printouts
				int refSubRunNo = 1;   // lowe school suggested "normally use subrun 1"
				int istat = 0;
				skbadch_(&skroot_ref_run, &refSubRunNo, &istat);
				/*
				*        istat  ; -1 : error
				*        istat  ;  0 : normal end  no modification
				*        istat  ;  1 : normal end  read only badch.XXXXXX
				*        istat  ;  2 : normal end  read only badch.XXXXXX.XXXXXX
				*        istat  ;  3 : normal end  read only /skam/const/badch.dat
				*        istat  ;  4 : normal end  read only badch2/badch.XXXXXX
				*        istat  ;  5 : normal end  read only badch2/badch.XXXX00.XXXX99
				*        istat  ;+10 : normal end  additional read /skam/const/badch.dat
				*/
				if(istat<0){
					Log(m_unique_name+" Error applying skbadch with reference run "+
					    toString(skroot_ref_run),v_error,m_verbose);
					return false;
				}
			}
		}
		
		// FIXME is there an order to the following actions?
		
		// for MC, use reference run to set parameters
		if(skrootMode!=SKROOTMODE::WRITE && isMC && skroot_ref_run>0
		      && (skhead_.nrunsk==0 || skhead_.nrunsk == 999999)){
			
			// we should override nrunsk according to harada-san
			skhead_.nrunsk = skroot_ref_run;
			
			// update water transparency
			Log(m_unique_name+" using reference run "+toString(skhead_.nrunsk)
				+" for water transparency",v_warning,m_verbose);
			float watert;
			int days_to_run_start = skday_data_.relapse[skhead_.nrunsk]; // defined in skdayC.h
			lfwater_(&days_to_run_start, &watert);
			Log(m_unique_name+" loaded new water transparency value "+toString(watert)
			    +" for run "+toString(skhead_.nrunsk),v_debug,m_verbose);
			// pass to downstream tools
			m_data->vars.Set("watert",watert);
			
			// initialise dark rates
			Log(m_unique_name+" Updating dark rates",v_debug,m_verbose);
			darklf_(&skroot_ref_run);
			
		}
		
		// initialize water transparency table
		skrunday_();
		
		// which of these should be called when??
		if(skheadg_.sk_geometry==4) skwt_gain_corr_();
		else skwt_();
		
	} else {
		// else not SK ROOT or zebra file; just a plain ROOT file.
		
		Log(m_unique_name+" creating MTreeReader to read tree "+treeName,v_debug,m_verbose);
		
		get_ok = myTreeReader.Load(list_of_files, treeName);
		if(not get_ok){
			Log(m_unique_name+" failed to open reader on tree "+treeName,v_error,m_verbose);
			return false;
		}
		if(mTreeReaderVerbosity) myTreeReader.SetVerbosity(mTreeReaderVerbosity);
		
		// for efficiency of reading, only enable used branches
		Log(m_unique_name+" activating branches",v_debug,m_verbose);
		if(SkippedInputBranches.size()){
			get_ok = myTreeReader.OnlyDisableBranches(SkippedInputBranches);
			if(!get_ok){
				Log(m_unique_name+" Did not recognise some branches in skipped branches list!",
				    v_error,m_verbose);
			}
		}
		if(ActiveInputBranches.size()){
			// only disable unlisted branches if we have a non-empty list of active branches
			get_ok = myTreeReader.OnlyEnableBranches(ActiveInputBranches);
			if(!get_ok){
				Log(m_unique_name+" Did not recognise some branches in active branches list!",
				    v_error,m_verbose);
			}
		}
	}
	
	// put the reader into the DataModel
	// Also register functions to load SHE / AFT commons, for access by other Tools if relevant.
	// We use std::mem_fn and std::bind to abstract away knowledge of the TreeReader class;
	// the DataModel can then provide access to these functions without dependance on this class.
	std::function<bool()> hasAFT  = std::bind(std::mem_fn(&TreeReader::HasAFT),  std::ref(*this));
	std::function<bool()> loadSHE = std::bind(std::mem_fn(&TreeReader::LoadSHE), std::ref(*this));
	std::function<bool()> loadAFT = std::bind(std::mem_fn(&TreeReader::LoadAFT), std::ref(*this));
	std::function<bool(int)> loadCommons = std::bind(std::mem_fn(&TreeReader::LoadCommons), std::ref(*this), std::placeholders::_1);
	std::function<int(long)> getTreeEntry = std::bind(std::mem_fn(&TreeReader::ReadEntry), std::ref(*this), std::placeholders::_1, false);
	// TODO we could remove the first argument now that the MTreeReader knows its name
	m_data->RegisterReader(readerName, &myTreeReader, hasAFT, loadSHE, loadAFT, loadCommons, getTreeEntry);
	
	// get first entry to process
	if(firstEntry<0) firstEntry=0;
	entrynum = firstEntry;
	
	// if we were given a selections file, only read entries that pass the specified cut
	if(selectionsFile!=""){
		// sanity check if given one in write mode
		if(skrootMode==SKROOTMODE::WRITE){
			Log(m_unique_name+" warning! selectionsFile given but we are in SKROOT WRITE mode!",
				v_warning,m_verbose);
		} else {
			
			// make the MTreeSelection to read the TEntryList file
			myTreeSelections = new MTreeSelection(selectionsFile);
			m_data->Selectors.emplace(readerName,myTreeSelections);
			
			if(cutName=="") cutName = myTreeSelections->GetTopCut();
			Log(m_unique_name+" reading only entries passing cut "+cutName
				+" in selections file "+selectionsFile,v_debug,m_verbose);
			
			// scan to the first entry passing our specified cut
			Log(m_unique_name+" scanning to first passing entry",v_debug,m_verbose);
			do {
				entrynum = myTreeSelections->GetNextEntry(cutName);
			} while(entrynum>0 && entrynum<firstEntry);
			if(entrynum<0){
				Log(m_unique_name+" was given both a selections file and a firstEntry,"
					+" but no passing entries were found after the specified starting entry!",v_error,m_verbose);
				return false;
			}
			Log(m_unique_name+" reading from entry "+toString(entrynum),v_debug,m_verbose);
		}
	}
	
	if(skreadMode!=0){
		// if calling skrawread, it seems like it doesn't work to just jump straight
		// to a given entry. So, very inefficient, but read all the entries in the chain
		// up to the requested entry
		// one of the things this does is set skpdst2_.iqbcalmode variable to 1 in entry 1 (at least in one instance);
		// note that one of the lomufit_gd reduction cuts is iqbcalmode!=0, so this is necessary
		// as part of the PreLoweReconstructionCuts Tool.
		Log(m_unique_name+" Scanning TChain to firstEntry "+toString(firstEntry),v_warning,m_verbose);
		// turn off verbosity for this scan as we know we're not interested in these events
		// and excessive printouts slow things down considerably
		int tmp_verb = m_verbose;
		m_verbose=0;
		for(int i=0; i<firstEntry; ++i){
			if(i%1000) Log(m_unique_name+" entry "+toString(i)+"...",v_error,m_verbose);
			ReadEntry(i,true);
		}
		m_verbose = tmp_verb;
	}
	
	return true;
}

bool TreeReader::Execute(){
	
	// nothing to do in write mode
	if(skrootMode==SKROOTMODE::WRITE) return true;
	
	// also offer the option to only read entries on request, not automatically on every entry
	// this can be useful if we aren't processing one entry per Execute loop,
	// but still want to use the TreeReader tool to populate SK common blocks.
	if(!autoRead) return true;
	
	Log(m_unique_name+" getting entry "+toString(entrynum),v_debug,m_verbose);
	
	// optionally buffer N entries per Execute call
	// clear the buffers before we start, unless we're buffering events between loops
	if(not (loadSheAftPairs && skrootMode==SKROOTMODE::ZEBRA)) FlushCommons();
	for(int buffer_entry_i=0; buffer_entry_i<entriesPerExecute; ++buffer_entry_i){
		
		// There are numerous reasons we may wish to skip entries from the input file:
		// + most entries in data SKROOT files are pedestal or status entries that don't actually
		//   contain detector data relating to a physics event. We'll usually want to skip these.
		// + i don't know how to determine the number of entries in a zbs file, so we may read off
		//   the end of our current zbs file, and need to re-try the read after loading the next file.
		// + if we're applying a trigger mask we'll skip entries whose trigger doesn't match.
		// + if we're only looking for SHE+AFT pairs, we'll skip all entries that aren't SHE
		//   with a follow-up AFT
		// + if we're reading a selection file, we'll skip all entries not specified by the selection
		// 
		// So, we'll need a loop that will keep reading entries until we find one we like.
		do {
			
			// load next entry
			Log(m_unique_name+" Reading entry "+toString(entrynum),v_debug,m_verbose);
			get_ok = ReadEntry(entrynum, true);
			Log(m_unique_name+" ReadEntry returned "+toString(get_ok),v_debug,m_verbose);
			
			// if we're processing ZBS files and ran off the end of this file,
			// load the next file if we have one and re-try the read.
			if(get_ok==0 && skrootMode==SKROOTMODE::ZEBRA && list_of_files.size()>0){
				Log(m_unique_name+" hit end of this ZBS file, loading next one",v_debug,m_verbose);
				skclosef_(&LUN);
				get_ok = LoadNextZbsFile();
				Log(m_unique_name+" loaded next ZBS file, return was "+toString(get_ok),v_debug,m_verbose);
				if(get_ok==0){
					Log(m_unique_name+" failure loading next ZBS file! Ending toolchain",v_error,m_verbose);
				} else {
					// loaded next file ok, repeat the read
					get_ok = -999;
				}
			}
			
			// trigger type checks - get the trigger bits
			if(get_ok>0 && skrootMode!=SKROOTMODE::NONE){
				
				std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
				
				// debug prints
				PrintTriggerBits();
				
				// apply our general check for required bits in the trigger mask
				// skip all events matching trigger types in skippedTriggers
				for(int bit_i=0; bit_i<skippedTriggers.size(); ++bit_i){
					int test_bit = skippedTriggers.at(bit_i);
					if(trigger_bits.test(test_bit)){
						get_ok=-999; // skip this event
					}
				}
				// alternatively to specifying every type we don't want,
				// we may choose to specify only the types we do want
				bool skipit=(allowedTriggers.size()>0); // only apply if given
				for(int bit_i=0; bit_i<allowedTriggers.size(); ++bit_i){
					int test_bit = allowedTriggers.at(bit_i);
					if(trigger_bits.test(test_bit)){
						skipit=false;
					}
				}
				if(skipit) get_ok=-999; // skip this event
				
				// if we're reading *only* SHE+AFT pairs, skip the entry if it's not SHE
				if(get_ok>0 && onlyPairs && !trigger_bits.test(28)){
					Log(m_unique_name+" Prompt entry is not SHE",v_debug,m_verbose);
					// its not SHE. If we only want SHE+AFT pairs, skip this entry.
					Log(m_unique_name+" Re-starting read process",v_debug,m_verbose);
					get_ok=-999;
				}
				
				// if we're reading SHE+AFT pairs and this is an SHE event,
				// search for and load the following AFT if there is one
				if(get_ok>0 && loadSheAftPairs && trigger_bits.test(28)){
					get_ok = AFTRead(entrynum);
					// returns: 1 if AFT loaded
					// returns: -100 if not AFT, but we're not only processing pairs
					// returns: -999 if not AFT (or error reading AFT) and we're only processing pairs
					// returns: <=0  if error during AFT read and we're not only processing pairs
				} else if(get_ok>0 && loadSheAftPairs){
					Log(m_unique_name+" PairLoading mode on but prompt event is not SHE, skipping follow-up read",
					    v_debug,m_verbose);
				}
				
			}  // else not SKROOT mode or bad prompt read. Skip trigger checks.
			
			// get the index of the next entry to read.
			if(myTreeSelections==nullptr){
				++entrynum;
				if(has_aft) ++entrynum; // we read two entries this round.
			} else {
				// get the entry number of the next entry we want to read
				entrynum = myTreeSelections->GetNextEntry(cutName);
				
				// if we were searching for SHE+AFT pairs, we may have the next entry
				// already loaded in the common block buffers.
				// If we don't want to read this entry, then discard it from the buffer.
				if(buffered_entry>0 && buffered_entry != entrynum){
					Log(m_unique_name+" Discarding buffered SHE from AFT search, since it is "
					   +"not in our selection entry list",v_debug,m_verbose);
					PopCommons();
				}
				buffered_entry = -1;
			}
			
			// TODO when reading SKROOT files we can't jump to specific entries,
			// (i think it doesn't populate the common blocks correctly?)
			// so we might need to read linearly and use the TreeSelections values
			// to keep reading until we get to an entry we want...
			
		} while(get_ok==-999);
		
		// if reading MC, we should override run number
		if( get_ok>0 && skrootMode!=SKROOTMODE::NONE &&
		    (skhead_.nrunsk==999999 || skhead_.nrunsk==0) ){
			skhead_.nrunsk = skroot_ref_run;
			skhead_.nsubsk = 1;
		}
		
		if(entriesPerExecute>1) PushCommons();
		if(get_ok==0) break;  // end of file
	} // read and buffer loop
	
	if(get_ok<0){
		// FIXME uhm, not sure waht this represents tbh.
		Log(m_unique_name+" get_ok returned "+toString(get_ok),v_error,m_verbose);
		assert(false);
		exit(-1);
		return false;
	}
	
	// when processing SKROOT files we can't rely on LoadTree(next_entry)
	// to indicate that there are more events to process - all remaining entries
	// may be pedestals that get skipped! If we skipped until we hit the end
	// of the TTree, bail out here as downstream tools will have no data to process.
	if(get_ok==0){
		m_data->vars.Set("Skip",true);
		m_data->vars.Set("StopLoop",true);
	} else {
		m_data->vars.Set("newRun",false);
		m_data->vars.Set("newSubrun",false);
		// check if we've changed run or subrun, and if so, update things like
		// water transparency and bad channel masking
		if(skrootMode!=SKROOTMODE::NONE && skhead_.nrunsk != last_nrunsk){
			RunChange();
			SubrunChange();
			last_nrunsk=skhead_.nrunsk;
			last_nsubsk=skhead_.nsubsk;
		} else if(skrootMode!=SKROOTMODE::NONE && skhead_.nsubsk != last_nsubsk){
			SubrunChange();
			last_nsubsk=skhead_.nsubsk;
		}
	}
	
	++readEntries;      // keep track of the number of entries we've actually returned
	if(readEntries%1000==0){
		Log(m_unique_name+" Read "+toString(readEntries),v_warning,m_verbose);
		// the TreeManager only calls TFile::Write on destructor,
		// so it seems like if we crash, we can end up losing everything.
		// intermittently invoke write
		if(skrootMode==SKROOTMODE::WRITE || skrootMode==SKROOTMODE::COPY){
			Write();
		}
	}
	aft_loaded = false; // whether common blocks currently hold SHE or AFT; reset this on every execution
	
	// check if we've hit the user-requested limit on number of entries to read
	if( ((maxEntries>0)&&(readEntries>=maxEntries)) ||
	    ((maxEntry>0)&&(entrynum>maxEntry)) ){
		Log(m_unique_name+" hit max events, setting StopLoop",v_message,m_verbose);
		m_data->vars.Set("StopLoop",1);
	}
	// use LoadTree to check if the next entry is valid without loading it
	// (this checks whether we've hit the end of the TTree/TChain)
	/*
	else if(myTreeReader.GetTree()){  // only possible if we have a TreeReader
		if(myTreeReader.GetTree()->LoadTree(entrynum)<0){
			Log(m_unique_name+" reached end of TTree, setting StopLoop",v_message,m_verbose);
			m_data->vars.Set("StopLoop",1);
		}
	}
	disabled as not needed now that we have can Skip, and this could potentially result
	in the TreeReader returning incorrect results if they go via the held TTree;
	e.g. myTreeReader->GetTree()->GetCurrentFile() may report the next file, when
	in fact the current in-memory event data relates to the last entry from the previous file.
	*/
	
	if(skrootMode!=::SKROOTMODE::NONE){
		Log(m_unique_name+" Returning entry skhead_.nevsk " + toString(skhead_.nevsk),v_debug,m_verbose);
		//std::cout<<"entry "<<entrynum<<", nevsk "<<skhead_.nevsk<<std::endl;
	} else {
		Log(m_unique_name+" Returning entry "+toString(entrynum),v_debug,m_verbose);
	}
	
	return true;
}

// FIXME no idea how much of this runchange / subrun change stuff is required
// either for MC or for data
bool TreeReader::RunChange(){
	
	// skip for invalid MC run numbers.
	if(skhead_.nrunsk==999999 || skhead_.nrunsk<=0){
		return true;
	}
	
	Log(m_unique_name+" Run Change",v_debug,m_verbose);
	m_data->vars.Set("newRun",true);
	get_ok = true;
	
	// check if this run is bad
	Log(m_unique_name+": Checking bad run flag for run "+toString(skhead_.nrunsk),v_debug,m_verbose);
	int isbad = lfbadrun_(&skhead_.nrunsk, &skhead_.nsubsk);
	if(isbad){
		Log(m_unique_name+" run "+toString(skhead_.nrunsk)+" flagged as a bad run by lfbadrun!",
		    v_warning,m_verbose);
		if(skipbadruns) SkipThisRun();
	}
	
	// update water transparency
	float watert;
	Log(m_unique_name+": Checking days since SK start",v_debug,m_verbose); // guessing what days_to_run_start is
	int days_to_run_start = skday_data_.relapse[skhead_.nrunsk]; // defined in skdayC.h
	if(days_to_run_start==0){
		Log(m_unique_name+" skday_data_.relapse returned 0 for run "+toString(skhead_.nrunsk)
		   /*+", skipping lfwater call!"*/, v_warning,m_verbose);
		// actually lfwater does return some nominal value even for invalid inputs,
		// which are probably better than nothing...?
	}
	// else {
		Log(m_unique_name+": Getting water transparency for day "+toString(days_to_run_start),v_debug,m_verbose);
		lfwater_(&days_to_run_start, &watert);
		Log(m_unique_name+" loaded new water transparency value "+toString(watert)
		    +" for run "+toString(skhead_.nrunsk),v_debug,m_verbose);
		// pass to downstream tools
		m_data->vars.Set("watert",watert);
	//}
	
	// update dark rate
	// What is this doing/required by?
	Log(m_unique_name+" Updating dark rates",v_debug,m_verbose);
	darklf_(&skhead_.nrunsk);
	
	return get_ok;
}

bool TreeReader::SubrunChange(){
	
	// skip for invalid MC run numbers
	if(skhead_.nrunsk==999999 || skhead_.nrunsk<=0){
		return true;
	}
	
	Log(m_unique_name+" SubRun Change",v_debug,m_verbose);
	m_data->vars.Set("newSubrun",true);
	get_ok = true;
	
	// update bad channels
	// read badch info & puts it into combad_ common block
	Log(m_unique_name+" Updating bad channel list for run "+toString(skhead_.nrunsk)
	    +", subrun "+toString(skhead_.nsubsk),v_debug,m_verbose);
	combad_.log_level_skbadch = 0;
	int ierr;
	skbadch_(&skroot_ref_run,&skhead_.nsubsk,&ierr);
	if(ierr<0){
		Log(m_unique_name+" Error calling skbadch_ in SubrunChange!",v_error,m_verbose);
		get_ok = false;
	} else {
		Log(m_unique_name+" bad channel list updated",v_debug,m_verbose);
	}
	
	return get_ok;
}

bool TreeReader::SkipThisRun(){
	Log(m_unique_name+" SkipThisRun called for run "+toString(skhead_.nrunsk)+", scanning for next run...",v_error,m_verbose);
	
	// probably the most efficient way would be to scan the set of filenames in the TChain,
	// (i.e. the list_of_files vector), parsing the run number from the filename
	// until we find the next file with a different run number.
	// how we then jump to that file though...? dunno. What's its correspoding TChain entry number?
	// It also requires on knowing the filename format and how to parse it.
	
	// the next best thing is just to keep jumping to the next file,
	// until we find one that has a new run number
	int next_run_num = skhead_.nrunsk;
	do {
		// get the number of entries in this file (TTree)
		int entry_in_current_file = myTreeReader.GetTree()->LoadTree(entrynum);
		int entries_in_current_file = myTreeReader.GetCurrentTree()->GetEntriesFast();
		Log(m_unique_name+" File "+myTreeReader.GetFile()->GetName()
		    +" has "+toString(entries_in_current_file)+" entries",v_debug,m_verbose);
		// jump forward to the first entry of the next file
		entrynum += (entries_in_current_file - entry_in_current_file);
		Log(m_unique_name+" Jumping forward to entry "+toString(entrynum),v_debug,m_verbose);
		
		// read that entry (should be first entry of next file), and get run number from the Header
		//  - ah, but this is only populated in some entries! so we need to scan until one is populated.
		while(true){
			get_ok = myTreeReader.GetEntry(entrynum);
			if(get_ok<=0){
				Log(m_unique_name+" Hit end of tree while looking for next good run!",v_warning,m_verbose);
				return false;
			}
			const Header* header = nullptr;
			get_ok = myTreeReader.Get("HEADER", header);
			if(!get_ok){
				Log(m_unique_name+" Error getting HEADER while scanning for good run",v_error,m_verbose);
				return false;
			}
			next_run_num = header->nrunsk;
			if(next_run_num!=0){
				get_ok=true;
				break;
			}
			++entrynum;
		}
		
		Log(m_unique_name+" Next file "+myTreeReader.GetFile()->GetName()
		    +" is from run "+toString(next_run_num),v_debug,m_verbose);
	} while (next_run_num==skhead_.nrunsk);
	
	
	/*
	// disable all branches and scan forward until run number changes
	TObjArray* branches = chain->GetListOfBranches();
	int nbranches = branches->GetEntries();
	std::vector<bool> bstatus(nbranches);
	for(int i=0; i<nbranches; ++i){
		TString bname = branches->At(i)->GetName();
		bstatus.at(i) = chain->GetBranchStatus(bname);
		chain->SetBranchStatus(bname, false);
	}
	
	// outer loop to keep going until run number changes
	int current_run_num = skhead_.nrunsk;
	int next_run_num = current_run_num;
	do {
		// scan forward until the file changes
		int current_tree_number = chain->GetTreeNumber();
		do {
			get_ok = chain->GetEntry(++entrynum);
			if(get_ok<0) break;
		} while (chain->GetTreeNumber()==current_tree_number);
		if(get_ok<0) break;
		chain->SetBranchStatus("HEADER",true);
		chain->GetEntry(entrynum);
		chain->SetBranchStatus("HEADER",false);
		// is this roundabout route necessary?
		TBranch* br = chain->GetBranch("HEADER");
		TLeaf* lf = (TLeaf*)br->GetListOfLeaves()->At(0);
		// or can we just cast the br to a TBranchElement here?
		TBranchElement* bev = (TBranchElement*)lf->GetBranch();
		const Header* hed = reinterpret_cast<const Header*>(bev->GetObject());
		next_run_num = hed->nrunsk;
	} while(next_run_num==current_run_num);
	
	// reenable the previously active branches
	for(int i=0; i<nbranches; ++i){
		TString bname = branches->At(i)->GetName();
		chain->SetBranchStatus(bname, bstatus.at(i));
	}
	*/
	
	return get_ok;
}

bool TreeReader::Write(){
	TreeManager* mgr = skroot_get_mgr(&LUN);
	if(!mgr) return false;
	TTree* otree = mgr->GetOTree();
	if(!otree) return false;
	TFile* ofile = otree->GetDirectory()->GetFile();
	if(!ofile || ofile->IsZombie()) return false;
	int nbyteswritten = ofile->Write(0,TObject::kOverwrite);
	Log(m_unique_name+": Wrote "+toString(nbyteswritten)+" to output file "
	    +ofile->GetName(),v_debug,m_verbose);
	return (nbyteswritten>=0);
}

bool TreeReader::Finalise(){
	
	if(myTreeSelections) delete myTreeSelections;
	
	if(skrootMode==SKROOTMODE::WRITE){
		TreeManager* mgr = skroot_get_mgr(&LUN);
		TTree* otree = mgr->GetOTree();
		//std::cout<<m_unique_name<<" tree has "<<otree->GetEntries()<<" entries in Finalise"<<std::endl;
	}
	
	
	// remove duplicate luns to prevent double deletion
	for(int& dup : duplicate_luns){
		EraseLun(dup);
	}
	
	// We could check our lunlist to see if there are any remaining TreeReaders,
	// and if this is the last one, clean up the SuperManager.
	// The issue with that is we can't check whether the SuperManager has any
	// outstanding TreeManagers, so if a user has done skroot_open_* themselves,
	// (and didn't create a corresponding entry in the lunlist)
	// then deleting the SuperManager would orphan the corresponding TreeManagers.
	// The best option is probably to just leave this to the OS - i.e. let it leak.
	//skroot_end_();
	
	return true;
}

TreeReader::~TreeReader(){
	
	// FIXME this shouldn't be needed? is AddTree messing it up?
	// do we also need to do this in copy mode now? if we're adding a tree?
	if(skrootMode==SKROOTMODE::WRITE){
		TreeManager* mgr = skroot_get_mgr(&LUN);
		TTree* otree = mgr->GetOTree();
		std::cout<<m_unique_name<<" tree has "<<otree->GetEntries()<<" entries in destructor"<<std::endl;
		std::cout<<m_unique_name<<" calling Write() on output file"<<std::endl;
		mgr->write();
	}
}

int TreeReader::ReadEntry(long entry_number, bool use_buffered){
	int bytesread=1;
	// load next entry data from TTree
	
	// skip the very first read in zebra mode as we already loaded it when checking if MC in Initialize
	if(skrootMode==SKROOTMODE::ZEBRA && entry_number==firstEntry){
		Log(m_unique_name+" skipping very first read as we got it from Initialize",v_debug,m_verbose);
	} else if(skrootMode!=SKROOTMODE::NONE){
		Log(m_unique_name+" ReadEntry using SK fortran routines to load data into common blocks",v_debug,m_verbose);
		// Populating fortran common blocks with SKROOT entry data requires using
		// SKRAWREAD and/or SKREAD.
		// These functions call various skroot_get_* functions to retrieve branch data.
		// They then use that data to populate the fortran common blocks.
		// * SKRAWREAD only works on data, not MC
		// * They both only call SOME branch getters
		
		// supposedly skrawread only reads the head, tq and pedestal branches < check this.
		// Supposedly skrawread loads some "constant tables" ...
		// These seem to be needed to call e.g. bonsai on data; if only skread
		// is called bonsai complains `neighbsk: error - ISEQ is <0 or > MAXPM! 0`
		// I have no idea what skread does.
		// TODO: Find out exactly what both of these functions do!
		// what do they read, calculate, and populate, along with any side-effects.
		// what does or doesn't get invoked in SKREAD when LUN is negative, and what
		// does this mean for the user. Is there an equivalent behaviour with skrawread?
		
		//   This will result in a segfault in SKROOT copy mode unless ALL remaining
		//   branches are either disabled or loaded by the user (e.g. via skroot_get_entry)
		
		//                         *** IMPORTANT ***
		// Both SKRAWREAD and SKREAD will invoke `skroot_next_entry` before data retrieval
		// *IF* the passed LUN is >0! This increments an entry number member of
		// the TreeManager. Subsequent `skroot_get_*` calls will load data from the
		// TTree entry corresponding to the internal entry number.
		// This allows one to call both SKRAWREAD and SKREAD, but ONLY in the correct
		// order and if SKREAD is given a negative LUN!
		// The sign of LUN also impacts whether various initialization stuff
		// is done, at least in SKREAD (not sure about SKRAWREAD).
		// So you MUST set LUN to be positive for ONE AND ONLY ONE CALL...probably.
		
		//                       *** VERY IMPORTANT ***
		// `skroot_jump_entry_(N) will actually set the internal index to (N-1)!
		// One needs to be very careful with calls to `skroot_jump_entry_`, `skrawread` and `skread`
		// to ensure the entry being read is indeed the entry you want!
		//                      ************************
		
		// =============================================================
		
		// first, set the internal entry number of the TreeManager.
		// remember this will actually set it to (entry_number - 1)
		// but we will then use either SKRAWREAD or SKREAD with a positive LUN
		// to increment it to the one we actually want, while also doing any
		// necessary internal reinitialization of these functions (presumably)
		if(skrootMode!=SKROOTMODE::ZEBRA){
			int entry_temp = static_cast<int>(entry_number);
			skroot_jump_entry_(&LUN, &entry_temp, &get_ok);
			//skroot_next_entry_(&LUN, &get_ok);
			if(get_ok==1){
				// ran off end of TTree!
				bytesread=0;
			}
		} else {
			// no idea how to specify the entry to read in zebra mode... TODO
			// for now we'll only support sequential reads
		}
		
		if(loadSheAftPairs && skrootMode==SKROOTMODE::ZEBRA && use_buffered && skhead_vec.size()>0){
			Log(m_unique_name+" buffered ZEBRA entry, using in place of read",v_debug,m_verbose);
			// if we have a buffered entry in hand, but it is not marked as an AFT trigger
			// for the current readout, then the buffered entry is an unprocessed event.
			// bypass the read and just load in the buffered data into the common blocks.
			LoadCommons(0);
			// then pop off the buffered data
			PopCommons();
		} else {
			Log(m_unique_name+" reading next entry from file",v_debug,m_verbose);
			// use skread / skrawread to get the next TTree entry and populate Fortran common blocks
			// skreadMode: 0=skread only, 1=skrawread only, 2=both
			if(bytesread>0 && skreadMode>0){
				Log(m_unique_name+" calling SKRAWREAD",v_debug,m_verbose);
				skcrawread_(&LUN, &get_ok); // N.B. positive LUN (see above)
				// for ZBS this doesn't seem to flag non-physics events as per for SKROOT...?
				// manually add in checks as per headsk.F for ROOT ... FIXME ? is this appropriate?
				if(get_ok==1){
					Log(m_unique_name+" read error "+toString(get_ok)+" calling skcrawread ",v_error,m_verbose);
					// lf_allfit actually continues the read loop if this is encountered,
					// so perhaps this is a recoverable error, or just an error relating to this entry?
					// FIXME if so it may be better to continue to next entry instead of bailing
					bytesread = -1;
				} else if(get_ok==2){
					// this just indicates we've reached the end of the file
					// if we've hit this, it's not good, because downstream tools
					// won't have any valid data!
					bytesread = 0;
				} else if(get_ok!=0) {
					// pedestal or status entry, no detector data, not an actual event
					// 3 = pedestal entry, 4 = runinfo entry.
					//Log(m_unique_name+" skrawread pedestal or status event, skipping",v_debug,m_verbose);
					// this happens a lot...
					if(skip_ped_evts) bytesread = -999;
					//std::cout<<"SKRAWREAD skipping event with status "<<get_ok<<std::endl;
				} else if(skrootMode==SKROOTMODE::ZEBRA){
					// get_ok == 0 -> normal physics event... but let's check manually.
					bool non_phys_evt = ((skhead_.nrunsk==0 && skhead_.mdrnsk!=0 && skhead_.mdrnsk!= 999999) ||
						                (std::bitset<8*sizeof(int)>(skhead_.ifevsk).test(19)) ||
						                (skhead_.nrunsk==0 && skheadg_.sk_geometry==0 && skhead_.ifevsk!=0) );
					/*
					std::cout<<"skhead_.nrunsk="<<skhead_.nrunsk<<"\n"
					         <<"skhead_.mdrnsk="<<skhead_.mdrnsk<<"\n"
					         <<"bit19="<<(std::bitset<8*sizeof(int)>(skhead_.ifevsk).test(19))<<"\n"
					         <<"skheadg_.sk_geometry="<<skheadg_.sk_geometry<<"\n"
					         <<"skhead_.ifevsk="<<skhead_.ifevsk<<std::endl;
					*/
					if(non_phys_evt) bytesread = -999; // doesn't look physical. FIXME is this correct???
				}
			}
			if(bytesread>0 && skreadMode!=1){  // skip skread if skrawread had an error
				Log(m_unique_name+" calling SKREAD",v_debug,m_verbose);
				int LUN2 = LUN;
				if(skreadMode==2) LUN2 = -LUN;  // if we already called skrawread, use a negative LUN
				skcread_(&LUN2, &get_ok);
				// for ZBS this doesn't seem to flag non-physics events as per for SKROOT...?
				// manually add in checks as per headsk.F for ROOT ... XXX ? is this appropriate?
				if(get_ok==1){
					// error reading entry
					Log(m_unique_name+" read error "+toString(get_ok)+" calling skcread ",v_error,m_verbose);
					bytesread = -1;
				} else if(get_ok==2){
					// end of file
					bytesread = 0;
				} else if(get_ok!=0) {
					// pedestal or status entry
					if(skip_ped_evts) bytesread = -999;
					//std::cout<<"SKREAD skipping event with status "<<get_ok<<std::endl;
				} else if(skrootMode==SKROOTMODE::ZEBRA){
					// get_ok == 0 -> normal read... but let's double check.
					bool non_phys_evt = ((skhead_.nrunsk==0 && skhead_.mdrnsk!=0 && skhead_.mdrnsk!= 999999) ||
						                (std::bitset<8*sizeof(int)>(skhead_.ifevsk).test(19)) ||
						                (skhead_.nrunsk==0 && skheadg_.sk_geometry==0 && skhead_.ifevsk!=0) );
					if(non_phys_evt) bytesread = -999; // manual non-physics override. FIXME is this correct?
				}
			}
			
			// As mentioned above, neither of these load all TTree branches.
			// To do that we need to call skroot_get_entry.
			if(bytesread>0 && skrootMode!=SKROOTMODE::ZEBRA){
				Log(m_unique_name+" calling skroot_get_entry",v_debug,m_verbose);
				skroot_get_entry_(&LUN);
				//skroot_get_tqskz_(&LUN); // might want to fix this - god I hate it all so much
			}
			if(bytesread > 0 && skrootMode == SKROOTMODE::ZEBRA){
			  Log(m_unique_name+" calling nerdnebk to retrieve NEUT bank", v_debug, m_verbose);
			  std::array<float, 3> interaction_pos = {};
			  nerdnebk_(interaction_pos.data());
			}
		}
		
	}
	if(bytesread >0 && skrootMode!=SKROOTMODE::ZEBRA) {
		Log(m_unique_name+" using MTreeReader to get next TTree entry",v_debug,m_verbose);
		// if in SKROOT mode we've already read from disk, just want to update
		// the internal MTreeReader variables, so skip the actual TTree::GetEntry call
		bytesread = myTreeReader.GetEntry(entry_number, (skrootMode!=SKROOTMODE::NONE));
	}
	Log(m_unique_name+" bytesread is "+toString(bytesread),v_debug,m_verbose);
	
	// stop loop if we ran off the end of the tree
	if(bytesread==0){
		Log(m_unique_name+" entry "+toString(entry_number)+" off end of input file!",v_error,m_verbose);
	} else if(bytesread==-999){
		Log(m_unique_name+" skrawread pedestal or status event",v_debug+10,m_verbose);
	}
	// stop loop if we had an error of some kind
	else if(bytesread<0){
		 if(bytesread==-10) Log(m_unique_name+" AutoClear error loading next input entry!",v_error,m_verbose);
		 else Log(m_unique_name+" IO error "+toString(bytesread)+" loading next input entry!",v_error,m_verbose);
		 m_data->vars.Set("StopLoop",1);
	}
	
	return bytesread;
}

int TreeReader::AFTRead(long entry_number){
	
	Log(m_unique_name+" Prompt entry is SHE, checking next entry for AFT", v_debug,m_verbose);
	has_aft=false; // default assumption
	
	// do a pre-check to see if we need to read the next entry.
	if(skrootMode!=SKROOTMODE::ZEBRA){
		// For ROOT files we can peek at the next TTree entry to determine whether it's
		// an AFT trigger without doing a full read.
		get_ok = CheckForAFTROOT(entry_number);
		// returns -1 on error trying to read AFT
		// -100 if no error but next entry is not AFT
		// -103 if AFT
	} else {
		// For zebra files we don't have this capability, so this always returns -103
		get_ok = CheckForAFTZebra(entry_number);
		// always returns -103
	}
	// get_ok = -103 requests a follow-up read (FIXME better return value)
	// If a follow-up read is requested, the current contents of the fortran commons
	// will have been put in the buffers, so they don't get overwritten by the next read.
	
	// if the pre-check didn't ask for a follow up read, we're done
	if(get_ok!=-103){
		// either due to an error reading, because we ran off the end of the file,
		// or just because this was not AFT, if we only want pairs we can skip this event
		if(onlyPairs) return -999;
		else return get_ok;
	}
	
	// if the pre-check indicated we need to do a follow up read, do that now.
	Log(m_unique_name+" Re-Invoking ReadEntry to check next entry",v_debug,m_verbose);
	
	// i assume that if there's an AFT, it'll always be the next entry,
	// i.e. there won't be things like status entries in between the SHE and AFT.
	get_ok = ReadEntry(entrynum+1, false);
	Log(m_unique_name+" Follow-up read returned "+toString(get_ok),v_debug,m_verbose);
	
	PrintTriggerBits();
	
	// check the follow-up read didn't encounter any errors
	if(get_ok>0){
		
		// We didn't. Check the follow-up entry to see if it's AFT.
		if(skrootMode!=SKROOTMODE::ZEBRA){
			get_ok = LoadAFTROOT();
			// always returns 1; we already knew this entry would be AFT
		} else {
			get_ok = LoadAFTZebra();
			// return -100: It wasn't an AFT, but onlyPairs is not set.
			//              The old SHE event has been put back into active commons.
			//              The next entry has been left in the buffer to use as a prompt candidate
			//              on the next Execute call.
			// return 1:    It was an AFT, we'll now have an SHE entry in active commons
			//              and the AFT entry in buffered commons.
			//              (this is always the case with ROOT, since we already peeked)
			// return -999: It wasn't an AFT, and onlyPairs is set.
			//              The old SHE will have been discarded, and the next entry
			//              will have been placed into the buffer to use as a new prompt candidate.
		}
		
		if(get_ok==1){
			Log(m_unique_name+" Successfully found SHE+AFT pair",v_debug,m_verbose);
			has_aft=true;
		} else if(get_ok == -100){
			// not AFT, but we're noy only reading pairs
			// the next entry has been left in buffers as a possible prompt candidate
			// for the next Execute call. Make a note of its entry number,
			// in case we're selecting only specific entries from the file and don't want it.
			buffered_entry = entry_number+1;
		}
		
	} else {
		// something went wrong with our follow-up read!
		// FIXME get_ok return values <0 all need to do this, change while loop?
		// or handle negative get_oks after the while loop?
	}
	
	
	return get_ok;
}

int TreeReader::CheckForAFTROOT(long entry_number){
	
	// this event is SHE, and we're looking for SHE+AFT pairs.
	// Peek at the next TTree entry to see if it's an associated AFT.
	Log(m_unique_name+" prompt event is SHE, peeking at next entry for AFT check",v_debug,m_verbose);
	
	int retval=-1;
	
	std::bitset<sizeof(int)*8> next_trigger_bits = 0;
	// make sure we don't read off the end of the tree
	if(myTreeReader.GetTree()->LoadTree(entry_number+1)>=0){
		// try to get the next HEAD entry
		get_ok = myTreeReader.GetTree()->GetBranch("HEADER")->GetEntry(entry_number+1);
		if(get_ok==0){
			Log(m_unique_name+" Error peeking next HEAD branch to look for AFT!",v_error,m_verbose);
			return -1;  // error reading AFT
		} else {
			const Header* header=nullptr;
			myTreeReader.Get("HEADER", header);
			next_trigger_bits = header->idtgsk;
		}
	} else {
		Log(m_unique_name+" can't check for AFT, no further entries in HEADER branch",
			v_debug,m_verbose);
			return -100;  // no error reading but no AFT
	}
	
	if(next_trigger_bits.test(29)){
		Log(m_unique_name+" next entry is AFT, requesting follow-up read",v_debug,m_verbose);
		// The next entry is indeed an AFT. We need to read it in properly now,
		// so buffer the current SHE data...
		PushCommons();
		// ... and indicate that we want to re-run ReadEntry to get the AFT entry.
		retval=-103;
	} else {
		Log(m_unique_name+" next entry is not AFT, no AFT this time.",v_debug,m_verbose);
		if(!onlyPairs){
			// if we're not explicitly requesting pairs we'll still process this SHE event
			// rewind Header branch so that anyone using the MTreeReader gets the right data
			myTreeReader.GetTree()->GetBranch("HEADER")->GetEntry(entry_number);
		}
		return -100;
	}
	
	return retval;
}

int TreeReader::CheckForAFTZebra(long entry_number){
	
	// At this point we should have an SHE entry loaded in the common blocks
	// "check if the next entry in the file is an AFT event."
	// To do this we need to read the next file entry into the fortran common blocks.
	// So that we don't lose the SHE data, buffer the current entry...
	PushCommons();
	
	// then ask ReadEntry to be called again, where we'll check if the next entry is an AFT.
	return -103;
}

int TreeReader::LoadAFTROOT(){
	// we peeked, so we already know this is an AFT trigger.
	Log(m_unique_name+" Successfully found SHE+AFT pair",v_debug,m_verbose);
	has_aft=true;
	
	// We now we have an SHE in the buffer and an AFT currently loaded.
	// Swap the two, because it's more intuitive that the SHE comes first.
	LoadCommons(0);
	
	return 1;
}

int TreeReader::LoadAFTZebra(){
	// At this point we currently have an unprocessed SHE event in the common block buffers,
	// and we've just read the next zebra file entry into the fortran common blocks.
	// Let's now check if the next zebra file entry is an AFT associated to our buffered SHE.
	Log(m_unique_name+" we have the next entry in active commons "
		+"and a prompt entry buffered. Checking trigger word",v_debug,m_verbose);
	
	int bytesread=-1;
	
	std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
	if(trigger_bits.test(29)){
		Log(m_unique_name+" Successfully found an SHE+AFT pair",v_debug,m_verbose);
		// this means we have an SHE in buffer and an AFT in the common blocks right now.
		// swap the SHE event back into the common blocks and AFT into the buffer.
		LoadCommons(0);
		// and that's it. This is the desired pair loading.
		bytesread = 1;
		
	} else {
		Log(m_unique_name+" Follow-up entry is not AFT",v_debug,m_verbose);
		
		// we have two options for proceeding here.
		// If the user ONLY wants SHE+AFT pairs...
		if(onlyPairs){
			Log(m_unique_name+" Dropping old SHE since we only want pairs",v_debug,m_verbose);
			// The currently buffered SHE did not have an associated AFT, so we have no use for it.
			PopCommons();  // drop it from the buffer.
			// we'll start this read all over, so put the new entry into the buffer
			// where it'll be used as a new prompt candidate on the next attempt.
			PushCommons();
			bytesread = -999;
		} else {
			Log(m_unique_name+" Swapping back previous entry, keeping next entry for next Execute call",
				v_debug,m_verbose);
			// else the user wants an AFT if there is one, but will still accept SHE events
			// without one. In that case, we still want to process our buffered entry,
			// so load it back into the common blocks, and retain our next entry in the buffer.
			LoadCommons(0);
			// ... we'll look again at the buffered entry on the next Execute call,
			// instead of reading from file.
			bytesread = -100;
		}
	}
	
	return bytesread;
}

void TreeReader::PrintTriggerBits(){
	std::bitset<sizeof(int)*8> trigger_bits = skhead_.idtgsk;
	
	Log(m_unique_name+" Trigger word for the active entry is: "
		+trigger_bits.to_string(),v_debug,m_verbose);
	if(m_verbose>(v_debug+1)){
		for(int i=0; i<(sizeof(int)*8); ++i){
			if(trigger_bits.test(i)) std::cout<<"bit "<<i<<" set"<<std::endl;
		}
	}
	/*
	std::bitset<sizeof(int)*8> event_bits = skhead_.ifevsk;
	Log(m_unique_name+" Event word for the active entry is: "+event_bits.to_string(),
		v_debug,m_verbose);
	if(m_verbose>(v_debug+1)){
		for(int i=0; i<(sizeof(int)*8); ++i){
			if(event_bits.test(i)) std::cout<<"bit "<<i<<" set"<<std::endl;
		}
	}
	*/
}

void TreeReader::PrintSubTriggers(){
	printf("calling runinf\n");
	runinfsk_();  // gets run info from SKROOT or ZBS file, populates /RUNINF/
	// is this already called as part of skread?
	// get trigger configuration for this run
	printf("This was event %d and had %d hardware triggers that were...{ ",skhead_.nevsk, skheadqb_.numhwsk);
	for(int i=0; i<skheadqb_.numhwsk; ++i){
		printf("%d, ", skheadqb_.hwsk[i]);
	}
	printf("}.\n");
	printf("The software ID for this event was %d\n",skheadqb_.nevswsk);
	printf( "subtrigger X, ntrigsk = %d, original T0 of the event (it0sk) = %d, "
			"\nT0 for PMT hits (it0xsk) = %d, time of subtrigger (t0_sub) = X, "
			"\ntime of first PMT hit (tiskz[0]) = %fns,"
			"raw time of first PMT hit(?) (tisk[0]) = %f,"
			"\nnqisk = %d, nqiskz = %d\n",
			/*i,*/ skheadqb_.ntrigsk, skheadqb_.it0sk, skheadqb_.it0xsk,
			/*t0_sub[0],*/ sktqz_.tiskz[0], skt_.tisk[0], skq_.nqisk, sktqz_.nqiskz);
	return;
	
	/*
	printf("calling softtrg_set_cond_\n");
	softtrg_set_cond_(...);
	// in contrast to slides, in examples it seems this accepts arguments
	// that must be retrieved from zebra blocks...don't think we should be calling it
	*/
	printf("calling get_sub_triggers_\n");
	// apply software trigger, obtain list of subtriggers of the requested type....
	std::bitset<sizeof(int)*8> trigbits;
	//trigbits.set(29);  // AFT
	//trigbits.set(28);  // SHE
	trigbits.set(2);
	int trig_id = (int)(trigbits.to_ulong());
	int ntrig;
	int MAX_TRIG=500;
	int* t0_sub = new int[MAX_TRIG];
	get_sub_triggers_(&trig_id, &ntrig, t0_sub, &MAX_TRIG);
	//  "get_sub_triggers: trgid=   536870912  is not supported yet."
	printf("The event had %d software triggers of type %d\n",ntrig, trig_id);
	if(ntrig==0) return;
	for(int i=0; i<ntrig; ++i){
		printf( "\tsubtrigger %d, ntrigsk = %d, original T0 of the event (it0sk) = %d, "
				"\n\tT0 for PMT hits (it0xsk) = %d, time of subtrigger (t0_sub) = %d, "
				"\n\ttime of first PMT hit (tiskz[0]) = %fns,"
				"raw time of first PMT hit(?) (tisk[0]) = %f,"
				"\n\tnqisk = %d, nqiskz = %d\n",
				i, skheadqb_.ntrigsk, skheadqb_.it0sk, skheadqb_.it0xsk, 
				t0_sub[i], sktqz_.tiskz[0], skt_.tisk[0], skq_.nqisk, sktqz_.nqiskz);
		int newt0 = skheadqb_.it0sk + t0_sub[i];
		set_timing_gate_(&newt0);   // populate tq variables etc
		// reset nqisk etc based on this subtrigger window
		int NEGLUN=-LUN;
		printf("doing skcread\n");
		skcread_(&NEGLUN, &get_ok);
		// now "modify PREVT0 if needed" it says - what does that mean???
		printf( "\t\tsubtrigger %d, ntrigsk = %d, original T0 of the event (it0sk) = %d, "
				"\n\t\tT0 for PMT hits (it0xsk) = %d, time of subtrigger (t0_sub) = %d, "
				"\n\t\ttime of first PMT hit (tiskz[0]) = %fns,"
				"raw time of first PMT hit(?) (tisk[0]) = %f,"
				"\n\t\tnqisk = %d, nqiskz = %d\n",
				i, skheadqb_.ntrigsk, skheadqb_.it0sk, skheadqb_.it0xsk, 
				t0_sub[i], sktqz_.tiskz[0], skt_.tisk[0], skq_.nqisk, sktqz_.nqiskz);
	}
	//set_timing_gate(skheadqb_.it0sk);
	//skroot_clear_tq_(&LUN);  // when do we need to call this?
	//delete_outside_hits_();
}

int TreeReader::LoadConfig(std::string configfile){
	Log(m_unique_name+" reading configuration file "+configfile,v_debug,m_verbose);
	// read the config file
	std::ifstream fin (configfile.c_str());
	std::string Line;
	
	if(not fin.is_open()){
		Log(m_unique_name+" failed to read configuration file "+configfile,v_error,m_verbose);
		return -1;
	}
	
	bool settingActiveInputBranches=false;
	bool settingActiveOutputBranches=false;
	bool settingSkippedInputBranches=false;
	bool settingSkippedOutputBranches=false;
	bool skFile=false;
	std::string allowedTriggersString="";
	std::string skippedTriggersString="";
	
	// scan over lines in the config file
	while (getline(fin, Line)){
		Log(m_unique_name+" parsing config line \""+Line+"\"",v_debug,m_verbose);
		// skip empty lines
		if (Line.empty()) continue;
		std::string LineCopy = Line; // make a copy so we can print it in case of parsing error
		// trim preceding whitespace
		Line.erase(0,Line.find_first_not_of(" \t\015"));
		// skip comment lines
		if(Line[0] == '#') continue;
		// trim line end comments (everything after and including a '#' character)
		if(Line.find('#')!=std::string::npos) Line.erase(Line.find_first_of('#'),std::string::npos);
		// trim trailing whitespace
		if(Line.find_last_not_of(" \t\n\015\014\013")!=std::string::npos)
			Line.erase(Line.find_last_not_of(" \t\n\015\014\013")+1,std::string::npos);
		
		// split apart the key and value
		std::string thekey   = Line.substr(0,Line.find_first_of(" \t\n\015\014\013"));
		std::string thevalue = Line.substr(Line.find_first_of(" \t\n\015\014\013")+1,std::string::npos);
		bool push_variable=true;  // whether to record this in m_variables (skip branch names and flags)
		
		// first check if we're entering or leaving the list of input or output branches to activate
		if (thekey=="StartActiveInputBranches"){
		  Log(m_unique_name+" starting active input branches list", v_debug, m_verbose);
			settingActiveInputBranches = true;
			push_variable=false;
		}
		else if(thekey=="EndActiveInputBranches"){
		  Log(m_unique_name+" ending active input branches list", v_debug, m_verbose);
			settingActiveInputBranches = false;
			push_variable=false;
		}
		else if(settingActiveInputBranches){
			ActiveInputBranches.push_back(Line);
			push_variable=false;
		}
		else if (thekey=="StartActiveOutputBranches"){
			settingActiveOutputBranches = true;
			push_variable=false;
		}
		else if(thekey=="EndActiveOutputBranches"){
			settingActiveOutputBranches = false;
			push_variable=false;
		}
		else if(settingActiveOutputBranches){
			ActiveOutputBranches.push_back(Line);
			push_variable=false;
		}
		// or the list of input or output branches to skip
		else if (thekey=="StartSkippedInputBranches"){
			settingSkippedInputBranches = true;
			push_variable=false;
		}
		else if(thekey=="EndSkippedInputBranches"){
			settingSkippedInputBranches = false;
			push_variable=false;
		}
		else if(settingSkippedInputBranches){
			SkippedInputBranches.push_back(Line);
			push_variable=false;
		}
		else if (thekey=="StartSkippedOutputBranches"){
			settingSkippedOutputBranches = true;
			push_variable=false;
		}
		else if(thekey=="EndSkippedOutputBranches"){
			settingSkippedOutputBranches = false;
			push_variable=false;
		}
		else if(settingSkippedOutputBranches){
			SkippedOutputBranches.push_back(Line);
			push_variable=false;
		}
		
		// other variables
		else if(thekey=="verbosity") m_verbose = stoi(thevalue);
		else if(thekey=="mTreeReaderVerbosity") mTreeReaderVerbosity = stoi(thevalue);
		else if(thekey=="inputFile") inputFile = thevalue;
		else if(thekey=="outputFile") outputFile = thevalue; // when using SKROOT copy mode
		else if(thekey=="FileListName") FileListName = thevalue;
		else if(thekey=="treeName") treeName = thevalue;
		else if(thekey=="readerName") readerName = thevalue;
		else if(thekey=="firstEntry") firstEntry = stoi(thevalue);
		else if(thekey=="maxEntries") maxEntries = stoi(thevalue);
		else if(thekey=="maxEntry") maxEntry = stoi(thevalue);
		else if(thekey=="selectionsFile") selectionsFile = thevalue;
		else if(thekey=="cutName") cutName = thevalue;
		else if(thekey=="skFile") skFile = stoi(thevalue);
		else if(thekey=="skrootMode") skrootMode = SKROOTMODE(stoi(thevalue));
		else if(thekey=="skreadMode") skreadUser = stoi(thevalue);
		else if(thekey=="LUN") LUN = stoi(thevalue);
		else if(thekey=="skoptn") skroot_options = thevalue;
		else if(thekey=="skbadopt") skroot_badopt = stoi(thevalue);
		else if(thekey=="mcReferenceRun") skroot_ref_run = stoi(thevalue);
		//else if(thekey=="skbadchrun") skroot_badch_ref_run = stoi(thevalue);
		// this reference run is used bad channel masking, dark rates and water transparency
		// not sure it make sense to use multiple different ones, but moreover it seems like
		// skhead_.nevsk is used internally by some algos so must be overridden for every MC readout
		// (or at least at all places before such algos are called - a global override is more manageable).
		//else if(thekey=="referenceWaterRun") reference_watert_run = stoi(thevalue);
		else if(thekey=="SK_GEOMETRY") sk_geometry = stoi(thevalue);
		else if(thekey=="skipPedestals") skip_ped_evts = stoi(thevalue);
		else if(thekey=="readSheAftTogether") loadSheAftPairs = stoi(thevalue);
		else if(thekey=="onlySheAftPairs") onlyPairs = stoi(thevalue);
		else if(thekey=="entriesPerExecute") entriesPerExecute = stoi(thevalue);
		else if(thekey=="allowedTriggers") allowedTriggersString = thevalue;
		else if(thekey=="skippedTriggers") skippedTriggersString = thevalue;
		else if(thekey=="skipBadRuns") skipbadruns = stoi(thevalue);
		else if(thekey=="autoEntryRead") autoRead = stoi(thevalue);
		// support for adding duplicate LUN numbers. This is rather silly because some SKOFL / ATMPD routines
		// hard-code the LUN number they read from, and if it's not matched to the one we're using, they either
		// read the wrong file, or dereference a pointer to a non-existent file and seg. Trouble is, LOWE group
		// has hard-coded use of LUN 10, while ATMPD has hard-coded use of LUN 1, so we can't even just make sure
		// we're using the expected number! So, we will allow referencing the same file under multiple LUNs
		else if(thekey.substr(0,3)=="LUN" && atoi(&thekey.data()[3])!=0){
			duplicate_luns.push_back(stoi(thevalue));
		}
		else {
			Log(m_unique_name+" error parsing config file line: \""+LineCopy
				+"\" - unrecognised variable \""+thekey+"\"",v_error,m_verbose);
		}
		if(push_variable){ m_variables.Set(thekey,thevalue); }
	}
	
	if(skFile==false){ skrootMode = SKROOTMODE::NONE; }
	if(loadSheAftPairs){
		if(entriesPerExecute>1){
			// SHE+AFT "pairs" represent two common blocks, so it doesn't make much sense to ask
			// to buffer 10 entries, and at the same time load SHE+AFT pairs together - pairs already
			// are two array entries, which is no different than with loadSheAftPairs=false
			// Only way to do this would be to implement some prescription for combining SHE+AFT,
			// e.g. by simply merging T and Q arrays..? but...that needs more thought, at least.
			Log(m_unique_name+" Error! readSheAftTogether and entriesPerExecute>1 cannot "
				+"both be used at the same time. Setting entriesPerExecute=1.",v_warning,m_verbose);
			entriesPerExecute=1;
		} else if(onlyPairs){
			skip_ped_evts = true;
		}
	}
	
	// parse the allowedTriggersString for allowed triggers
	while(true){
		try{
			size_t nextchar=0;
			int next_bit = stoi(allowedTriggersString,&nextchar);
			allowedTriggers.push_back(next_bit);
			if(nextchar==allowedTriggersString.length()) break;
			allowedTriggersString = allowedTriggersString.substr(nextchar+1,std::string::npos);
		} catch (...) { break; }
	}
	// same for skippedTriggerString
	while(true){
		try{
			size_t nextchar=0;
			int next_bit = stoi(skippedTriggersString,&nextchar);
			skippedTriggers.push_back(next_bit);
			if(nextchar==skippedTriggersString.length()) break;
			skippedTriggersString = skippedTriggersString.substr(nextchar+1,std::string::npos);
		} catch (...) { break; }
	}
	
	// done parsing, close config file
	fin.close();
	
	return 1;
}

void TreeReader::CloseLUN(){
	
	// close input skroot files, delete the TTreeManager
	// the SuperManager has a nullptr check so shouldn't seg even if this LUN is invalid
	if(skrootMode!=SKROOTMODE::ZEBRA) skroot_close_(&LUN);
	else skclosef_(&LUN);
	// note that if the LUN wasn't valid, that LUN will now be created in the SuperManager's map
	// with a nullptr entry, which completely locks the LUN - it doesn't point to a valid TreeManager,
	// and it's not possible to remove it from the list. TODO fix the SuperManager.
	
	// free the LUN number
	m_data->FreeLUN(LUN, readerName);
	
}

bool TreeReader::LoadNextZbsFile(){
	// get the next file
	std::string next_file = list_of_files.back();
	list_of_files.pop_back();
	// resolve any environmental variables and symlinks
	std::string cmd = std::string("readlink -f ")+next_file;
	Log(m_unique_name+" getting return from command '"+cmd+"'",v_debug+1,m_verbose);
	//next_file = getOutputFromFunctionCall(system, cmd.c_str());  // was crashing???
	next_file = getOutputFromFunctionCall(safeSystemCall, cmd);
	Log(m_unique_name+": next ZBS file "+next_file,v_debug,m_verbose);
	
	// ok now actually open the ZBS file.
	/*
	set_rflist_(&LUN, next_file.c_str(), "LOCAL", "", "RED", "", "", "recl=5670 status=old", "", "",
		        next_file.length(), 5, 0, 3, 0, 0, 20, 0, 0); // lengths of all string args
	// need to change "LOCAL" to "DISK" and 5 to 4 if file is on /disk... instead of elsewhere....?
	
	int fileIndex = 1; // 1st file in rflist (??? does it support multiple ???)
	int ihndl=1;
	char ftype='z';
	skopenf_(&LUN, &fileIndex, &ftype, &get_ok, &ihndl);
	*/
	
	set_rflist_zbs( LUN, next_file.c_str(), false );
	int ipt = 1;
	int ihndl=1;
	skopenf_( &LUN, &ipt, "Z", &get_ok, &ihndl );
	
	if(get_ok!=0){
		Log(m_unique_name+" Error loading next ZBS file '"+next_file,v_error,m_verbose);
		return false;
	} else {
		Log(m_unique_name+" next ZBS file '"+next_file+"' has been loaded",v_debug,m_verbose);
	}
	return true;
	
	/* for reference:
	// we also have the wrapper fort_fopen, which seems to be a slightly simpler interface
	// See '$RELICWORKDIR/data_reduc/third/tools_nov19/root2zbs/fort_fopen.F' (only supports writing)
	// or  '$ATMPD_ROOT/src/analysis/tutorials/fort_fopen.F'                  (supports read and write)
	//  see '$RELICWORKDIR/lomufit/mufit/src.sonia/check_mc/zbs_double_events.F'
	//  for writing ZBS files with set_rflist.)
	// See '$SKOFL_ROOT/src/iolib/skopenf.f' for some options e.g. on the meaning of the ftype argument
	int rw = 0;                     // = 0 reading, = 1 for writing
	char ftype = 'z';               // 'z' = zebra is probably the only one we need
	fort_fopen_(&LUN , next_file.c_str(), &ftype, &rw, next_file.length());
	// also note fort_fopen doesn't have any error returns - it has no return type, nor can i find
	// a c++ wrapper that provides one (though one would be easy to make...)
	// moreover ALL of these eventually call skopenf, which will invoke STOP on error....
	//$SKOFL_ROOT/src/recon/ska/zform_open.c
	*/
	
}

int TreeReader::PushCommons(){
	if(loadSheAftPairs && skhead_vec.size()){
		std::cerr<<"PUSH COMMONS WITH ALREADY EXISTING ENTRY!"<<std::endl;
		exit(-1);
	}
	// make a buffered copy of the current state of event-wise fortran common blocks
	// so that the user may access both SHE and AFT (or potentially arbitrary) events
	// TODO remove any that are not set by skread and used by reco algorithms
	// XXX we could consider using or looking at `skroot_set_tree_(&lun);`
	// which populates the SKROOT branches based on common blocks.
	// Perhaps we could call this and then buffer the generated e.g. TQREAL objects?
	
	// event header - run, event numbers, trigger info...
	skhead_vec.emplace_back(skhead_);            // skhead_common
	skheada_vec.emplace_back(skheada_);          // skheada_common
	skheadg_vec.emplace_back(skheadg_);          // skheadg_common
	skheadf_vec.emplace_back(skheadf_);          // skheadf_common
	skheadc_vec.emplace_back(skheadc_);          // skheadc_common
	skheadqb_vec.emplace_back(skheadqb_);        // skheadqb_common
	
	// low-e event variables
	skroot_lowe_vec.emplace_back(skroot_lowe_);  // skroot_lowe_common
	skroot_mu_vec.emplace_back(skroot_mu_);      // skroot_mu_common
	skroot_sle_vec.emplace_back(skroot_sle_);    // skroot_sle_common
	
	// commons containing arrays of T, Q, ICAB....
	// not sure which of these may be populated by skread
	skq_vec.emplace_back(skq_);                  // skq_common
	skqa_vec.emplace_back(skqa_);                // skqa_common
	skt_vec.emplace_back(skt_);                  // skt_common
	skta_vec.emplace_back(skta_);                // skta_common
	skchnl_vec.emplace_back(skchnl_);            // skchnl_common
	skthr_vec.emplace_back(skthr_);              // skthr_common
	sktqz_vec.emplace_back(sktqz_);              // sktqz_common
	sktqaz_vec.emplace_back(sktqaz_);            // sktqaz_common
	rawtqinfo_vec.emplace_back(rawtqinfo_);      // rawtqinfo_common
	
	sktrighit_vec.emplace_back(sktrighit_);      // sktrighit_common
	skqv_vec.emplace_back(skqv_);                // skqv_common
	sktv_vec.emplace_back(sktv_);                // sktv_common
	skchlv_vec.emplace_back(skchlv_);            // skchlv_common
	skthrv_vec.emplace_back(skthrv_);            // skthrv_common
	skhitv_vec.emplace_back(skhitv_);            // skhitv_common
	skpdstv_vec.emplace_back(skpdstv_);          // skpdstv_common
	skatmv_vec.emplace_back(skatmv_);            // skatmv_common
	
	// OD mask....? nhits, charge, flag...?
	odmaskflag_vec.emplace_back(odmaskflag_);    // odmaskflag_common
	
	// hardware trigger variables; counters, trigger words, prevt0...
	// no idea how many, if any, of these are populated by skread,
	// or moreover how many are needed by reconstruction algorithms.
	
	// spacer and trigger info.
	skdbstat_vec.emplace_back(skdbstat_);        // skdbstat_common
	skqbstat_vec.emplace_back(skqbstat_);        // skqbstat_common
	skspacer_vec.emplace_back(skspacer_);        // skspacer_common
	
	// gps word and time.
	skgps_vec.emplace_back(skgps_);              // skgps_common
	t2kgps_vec.emplace_back(t2kgps_);            // t2kgps_common
	
	// hw counter difference to previous event.
	prevt0_vec.emplace_back(prevt0_);            // prevt0_common
//	tdiff_vec.emplace_back(tdiff_);              // tdiff_common
	mintdiff_vec.emplace_back(mintdiff_);        // mintdiff_common
	
	// trigger hardware counters, word, spacer length...
	sktrg_vec.emplace_back(sktrg_);              // sktrg_common
	
	// MC particles and vertices, event-wise.
	vcvrtx_vec.emplace_back(vcvrtx_);            // vcvrtx_common
	vcwork_vec.emplace_back(vcwork_);            // vcwork_common
	
	return skhead_vec.size();
}

int TreeReader::PopCommons(){
	// drop an entry from the buffered common blocks
	// TODO remove any that are not set by skread and used by reco algorithms
	
	skhead_vec.pop_back();
	skheada_vec.pop_back();
	skheadg_vec.pop_back();
	skheadf_vec.pop_back();
	skheadc_vec.pop_back();
	skheadqb_vec.pop_back();
	skroot_lowe_vec.pop_back();
	skroot_mu_vec.pop_back();
	skroot_sle_vec.pop_back();
	skq_vec.pop_back();
	skqa_vec.pop_back();
	skt_vec.pop_back();
	skta_vec.pop_back();
	skchnl_vec.pop_back();
	skthr_vec.pop_back();
	sktqz_vec.pop_back();
	sktqaz_vec.pop_back();
	rawtqinfo_vec.pop_back();
	sktrighit_vec.pop_back();
	skqv_vec.pop_back();
	sktv_vec.pop_back();
	skchlv_vec.pop_back();
	skthrv_vec.pop_back();
	skhitv_vec.pop_back();
	skpdstv_vec.pop_back();
	skatmv_vec.pop_back();
	odmaskflag_vec.pop_back();
	skdbstat_vec.pop_back();
	skqbstat_vec.pop_back();
	skspacer_vec.pop_back();
	skgps_vec.pop_back();
	t2kgps_vec.pop_back();
	prevt0_vec.pop_back();
//	tdiff_vec.pop_back();
	mintdiff_vec.pop_back();
	sktrg_vec.pop_back();
	vcvrtx_vec.pop_back();
	vcwork_vec.pop_back();
	
	return skhead_vec.size();
}

int TreeReader::FlushCommons(){
	// drop all entries from the buffered common blocks
	// TODO remove any that are not set by skread and used by reco algorithms
	
	skhead_vec.clear();
	skheada_vec.clear();
	skheadg_vec.clear();
	skheadf_vec.clear();
	skheadc_vec.clear();
	skheadqb_vec.clear();
	skroot_lowe_vec.clear();
	skroot_mu_vec.clear();
	skroot_sle_vec.clear();
	skq_vec.clear();
	skqa_vec.clear();
	skt_vec.clear();
	skta_vec.clear();
	skchnl_vec.clear();
	skthr_vec.clear();
	sktqz_vec.clear();
	sktqaz_vec.clear();
	rawtqinfo_vec.clear();
	sktrighit_vec.clear();
	skqv_vec.clear();
	sktv_vec.clear();
	skchlv_vec.clear();
	skthrv_vec.clear();
	skhitv_vec.clear();
	skpdstv_vec.clear();
	skatmv_vec.clear();
	odmaskflag_vec.clear();
	skdbstat_vec.clear();
	skqbstat_vec.clear();
	skspacer_vec.clear();
	skgps_vec.clear();
	t2kgps_vec.clear();
	prevt0_vec.clear();
//	tdiff_vec.clear();
	mintdiff_vec.clear();
	sktrg_vec.clear();
	vcvrtx_vec.clear();
	vcwork_vec.clear();
	
	return 1;
}

bool TreeReader::LoadCommons(int buffer_i){
	// check we have such a buffered entry
	if(buffer_i>=skhead_vec.size()){
		Log(m_unique_name+" Error! Asked to load common block buffer entry "+toString(buffer_i)
			+" out of range 0->"+skhead_vec.size()+"!",v_error,m_verbose);
		return false;
	}
	
	std::swap(skhead_, skhead_vec.at(buffer_i));
	std::swap(skheada_, skheada_vec.at(buffer_i));
	std::swap(skheadg_, skheadg_vec.at(buffer_i));
	std::swap(skheadf_, skheadf_vec.at(buffer_i));
	std::swap(skheadc_, skheadc_vec.at(buffer_i));
	std::swap(skheadqb_, skheadqb_vec.at(buffer_i));
	std::swap(skroot_lowe_,skroot_lowe_vec.at(buffer_i));
	std::swap(skroot_mu_,skroot_mu_vec.at(buffer_i));
	std::swap(skroot_sle_,skroot_sle_vec.at(buffer_i));
	std::swap(skq_, skq_vec.at(buffer_i));
	std::swap(skqa_, skqa_vec.at(buffer_i));
	std::swap(skt_, skt_vec.at(buffer_i));
	std::swap(skta_, skta_vec.at(buffer_i));
	std::swap(skchnl_, skchnl_vec.at(buffer_i));
	std::swap(skthr_, skthr_vec.at(buffer_i));
	std::swap(sktqz_, sktqz_vec.at(buffer_i));
	std::swap(sktqaz_, sktqaz_vec.at(buffer_i));
	std::swap(rawtqinfo_, rawtqinfo_vec.at(buffer_i));
	std::swap(sktrighit_, sktrighit_vec.at(buffer_i));
	std::swap(skqv_, skqv_vec.at(buffer_i));
	std::swap(sktv_, sktv_vec.at(buffer_i));
	std::swap(skchlv_, skchlv_vec.at(buffer_i));
	std::swap(skthrv_, skthrv_vec.at(buffer_i));
	std::swap(skhitv_, skhitv_vec.at(buffer_i));
	std::swap(skpdstv_, skpdstv_vec.at(buffer_i));
	std::swap(skatmv_, skatmv_vec.at(buffer_i));
	std::swap(odmaskflag_, odmaskflag_vec.at(buffer_i));
	std::swap(skdbstat_, skdbstat_vec.at(buffer_i));
	std::swap(skqbstat_, skqbstat_vec.at(buffer_i));
	std::swap(skspacer_, skspacer_vec.at(buffer_i));
	std::swap(skgps_, skgps_vec.at(buffer_i));
	std::swap(t2kgps_, t2kgps_vec.at(buffer_i));
	std::swap(prevt0_, prevt0_vec.at(buffer_i));
//	std::swap(tdiff_, tdiff_vec.at(buffer_i));   // segfauts?????
	std::swap(mintdiff_, mintdiff_vec.at(buffer_i));
	std::swap(sktrg_, sktrg_vec.at(buffer_i));
	std::swap(vcvrtx_, vcvrtx_vec.at(buffer_i));
	std::swap(vcwork_, vcwork_vec.at(buffer_i));
	
	return true;
}

bool TreeReader::HasAFT(){
	// check if this entry has an AFT event buffered
	return has_aft;
}

bool TreeReader::LoadAFT(){
	Log(m_unique_name+" LoadAFT called: has_aft="+toString(has_aft)+", aft_loaded="+toString(aft_loaded),
	    v_debug,m_verbose);
	if(has_aft && !aft_loaded){
		aft_loaded = LoadCommons(0);
		return aft_loaded;
	} // else either already loaded, or no AFT to load
	return aft_loaded;
}

bool TreeReader::LoadSHE(){
	Log(m_unique_name+" LoadSHE called: has_aft="+toString(has_aft)+", aft_loaded="+toString(aft_loaded),
	    v_debug,m_verbose);
	if(has_aft && aft_loaded){
		aft_loaded = !LoadCommons(0);
		return aft_loaded;
	} // else SHE already loaded
	return true;
}
