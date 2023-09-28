#ifndef PreLoweReconstructionCuts_H
#define PreLoweReconstructionCuts_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class PreLoweReconstructionCuts
*
* relic candidate cuts we can do before lowe reconstruction
*
* $Author: M.O'Flaherty $
* $Date: 2023/08/30 $
* $Contact: moflaher@km.icrr.u-tokyo.ac.jp
*/

class PreLoweReconstructionCuts: public Tool {
	
	public:
	
	PreLoweReconstructionCuts();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	std::string selectorName;
	
};


#endif
