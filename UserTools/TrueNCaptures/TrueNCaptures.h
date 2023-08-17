#ifndef TrueNCaptures_H
#define TrueNCaptures_H

#include <string>
#include <iostream>

#include "Tool.h"

class TFile;
class TTree;


/**
 * \class TrueNCaptures
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class TrueNCaptures: public Tool {
	
	public:
	
	 TrueNCaptures(); ///< Simple constructor
	 bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	 bool Execute(); ///< Executre function used to perform Tool perpose. 
	 bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	bool PrintCaptures();
	bool MakePlots(int step);
	
	std::string plotsfile="";
	TFile* fplots=nullptr;
	TTree* tplots=nullptr;
	std::map<std::string, double> dbranchvars;
	std::map<std::string, int> ibranchvars;
	std::map<std::string, std::vector<double>> vbranchvars;
	
};


#endif
