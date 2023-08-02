#ifndef WallCut_H
#define WallCut_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "skroot.h"
#include "ConnectionTable.h"


/**
 * \class WallCut
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class WallCut: public Tool {
	
	
	public:
	
	WallCut(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	
	private:
	
	std::string treeReaderName;
	
	MTreeReader* myTreeReader = nullptr;
	
	const LoweInfo* myLowe=nullptr;
	
	int verbosity = 1;
	int m_verbose;
	int v_error = 0;
	
	basic_array<float> reconVertex;
	
	int Nskipped = 0;
};


#endif
