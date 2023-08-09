#ifndef ODCut_H
#define ODCut_H

#include <string>
#include <iostream>
#include <bitset>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "ConnectionTable.h"


/**
 * \class ODCut
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class ODCut: public Tool {
	
	
	public:
	
	ODCut(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	
	private:
	MTreeReader* myTreeReader = nullptr;
	
	const Header* myHeader = nullptr;
	
	std::string treeReaderName;
	int verbosity;
	int v_error = 0;
	
	std::bitset<sizeof(int)*8>  triggerID;
	
	int Nskipped = 0;
	
	
};


#endif
