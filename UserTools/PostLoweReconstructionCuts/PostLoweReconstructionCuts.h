#ifndef PostLoweReconstructionCuts_H
#define PostLoweReconstructionCuts_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class PostLoweReconstructionCuts
*
* Initial relic reductions based on quality of lowe reconstruction. Based on `make_precut.F`
*
* $Author: M.O'Flaherty $
* $Date: 2023/08/30 $
* $Contact: moflaher@km.icrr.u-tokyo.ac.jp
*/

class PostLoweReconstructionCuts: public Tool {
	
	public:
	
	PostLoweReconstructionCuts();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	std::string selectorName;
	bool getLoweVarsFromFile=false;
};


#endif
