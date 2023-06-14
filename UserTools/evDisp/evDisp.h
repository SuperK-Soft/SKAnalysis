#ifndef evDisp_H
#define evDisp_H

#include <string>
#include <iostream>
#include <bitset>

#include "Tool.h"
#include "MTreeReader.h"


/**
 * \class evDisp
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class TH2D;
class TGraph2D;
class TCanvas;
class TVirtualPad;

#include "skroot.h"
#include "ConnectionTable.h"

class evDisp: public Tool {
	
	public:
	
	evDisp(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	bool GetData();
	int FindHits();
	bool BranchCheck();
	
	// tool variables
	// ==============
	std::string toolName;
	std::string treeReaderName;
	MTreeReader* myTreeReader=nullptr;
	
	int plotVar=0;    // 0=charge, 1=time.
	int dataSrc=0;    // 0=sktqz_ common block, 1=TQReal branch
	int evtSrc=0;     // 0=SHE, 1=AFT
	int plotStyle=0;  // 0=points, 1=histogram
	int histSet=0;     // 0=time and charge, 1=T-TOF and charge (temporary)
	int inGateOnly=0; // 0 = only plot in-gate hits
	
	float cap_rad=0;
	float max_barrel_x=0;
	float min_barrel_x=3.2;
	
	const Header* myHeader=nullptr;
	const TQReal* myTQReal=nullptr;
	const SecondaryInfo* mySecondary=nullptr;
	const LoweInfo* myLowe=nullptr;
	ConnectionTable* myConnectionTable=nullptr;
	int totalPMTsActivated=0;
	int cableNumber = 0;
	float charge = 0;
	float time = 0;
	bool in_gate = 0;
	float var = 0;
	float tubePosition[3];
	float tubeRadialCoordinate;
	float tubeAngularCoordinate;
	int zMax = 1810;
	int zMin = -1810;
	int runNum = 0;
	int subrunNum = 0;
	int eventNum = 0;
	int gateWidth = 0;
	int loweFlag = 0;
	unsigned int readoutT0;
	
	float recEnergy;
	basic_array<float> recVertex;
	basic_array<float> recDirection;
	
	std::bitset<sizeof(int)*8>  triggerID;
	std::bitset<sizeof(int)*8> it0xskbitset;
	
	TH2D* topCapHeatMap=nullptr;
	TH2D* bottomCapHeatMap=nullptr;
	TH2D* barrelHeatMap=nullptr;
	TGraph2D* topCapHitMap=nullptr;
	TGraph2D* bottomCapHitMap=nullptr;
	TGraph2D* barrelHitMap=nullptr;
	TH1F* hitsVTimeHist=nullptr;
	TH1F* hitsVChargeHist=nullptr;
	
	TCanvas* displayCanvas=nullptr;
	TVirtualPad* displayPad=nullptr;
	TVirtualPad* botDisplayPad=nullptr;
	
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage="";
	int get_ok=0;
	
};


#endif
