/* vim:set noexpandtab tabstop=4 wrap */
#ifndef RootReadTest_H
#define RootReadTest_H

#include <string>
#include <iostream>
#include <map>

#include "Tool.h"
#include "MTreeReader.h"

#include "SkrootHeaders.h" // MCInfo, Header etc.


/**
* \class RootReadTest
*
* A tool to test the MTreeReader
*
* $Author: M.O'Flaherty $
* $Date: 2020/08/18 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class RootReadTest: public Tool {
	
	public:
	
	RootReadTest(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	int ReadEntry(long entry_number);
	int GetBranchesNtuple();
	int GetBranchesSKROOT();
	int CheckEntryNtuple();
	int CheckEntrySKROOT();
	
	MTreeReader myTreeReader;
	std::string inputFile;
	std::string treeName;
	std::string testFileType;
	int entrynum=0;
	int maxEvents=-1;
	int get_ok;
	
	// what's annoying is that we can deduce these types at runtime from the branch title,
	// but construction of the wrapper class needs to know them at compile time.
	// note also: you cannot use basic_array<float**> instead of basic_array<float(*)[3]>
	// because they are fundamentally different types. Sorry, it's the best solution i can come up with.
	int n_secondaries_2;
	basic_array<int*> secondary_PDG_code_2;            // 1D array
	basic_array<float(*)[3]> secondary_start_vertex_2; // 2D array
	
	const MCInfo* mc_info=nullptr;
	const Header* file_header=nullptr;
	
};



#endif
