#ifndef PlotHitTimes_H
#define PlotHitTimes_H

#include <string>
#include <iostream>
#include <vector>
#include <bitset>

#include "Tool.h"
#include "ColourWheel.h"


/**
 * \class PlotHitTimes
 *
 * This Tool plots the times of hits, highlighting in-gate and out-of-gate hits, along with hits in the main and subtriggers
*
* $Author: M.O'Flaherty $
* $Date: 2023/08/01 $
*/

class TFile;
class TCanvas;
class MTreeReader;

class PlotHitTimes: public Tool {
	
	public:
	
	PlotHitTimes(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	int GetSubtriggerFlags(int subtrigtype, std::vector<std::bitset<32>>& in_subtrigger_flags, int n_triggers);
	
	private:
	std::string m_unique_name="";
	
	MTreeReader* myTreeReader=nullptr;
	TFile* fout = nullptr;
	TCanvas* c_subtriggers = nullptr;
	ColourWheel colourwheel;
	
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
};


#endif
