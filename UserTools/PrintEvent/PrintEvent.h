#ifndef PrintEvent_H
#define PrintEvent_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
 * \class PrintEvent
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class PrintEvent: public Tool {
	
	public:
	
	PrintEvent(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	bool PrintTQRealInfo(bool ID, bool verbose=false);
	bool PrintSubTriggers(bool verbose=false);
	bool PrintSubTriggersMC(bool verbose=false);
	bool PrintLowEInfo();
	bool PrintATMPDInfo();
	bool PrintMuInfo();
	bool PrintUpMuInfo();
	bool PrintSLEInfo();
	bool PrintPrevT0();
	bool PrintHWTriggerInfo();
	bool PrintSWTriggerInfo();
	bool PrintHeaderInfo();
	bool PrintMCInfo();
	bool PrintRunInfo();
	bool PrintBadChannels();
	bool PrintDarkInfo();
	bool PrintTriggerInfo();
	bool Print_sktrg();
	bool PrintHits();
	bool PrintTQRealHits(bool ID, int nhits);
	bool PrintTQCommons(bool ID, int nhits);
	bool PrintTQZCommons(bool ID, int nhits);
	bool PrintTQRawCommons(bool ID, int nhits);
	
	MTreeReader* myTreeReader=nullptr;
	
};


#endif
