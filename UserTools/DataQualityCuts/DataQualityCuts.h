#ifndef DataQualityCuts_H
#define DataQualityCuts_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class DataQualityCuts
*
* Cuts before any kind of reconstruction. Based on `make_precut` and incorporating `lowfit_gd.F`
*
* $Author: M.O'Flaherty $
* $Date: 2023/08/30 $
* $Contact: moflaher@km.icrr.u-tokyo.ac.jp
*/

class DataQualityCuts: public Tool {
	
	public:
	
	DataQualityCuts();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	std::string selectorName;
	
};


#endif
