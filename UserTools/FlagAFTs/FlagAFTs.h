#ifndef FlagAFTs_H
#define FlagAFTs_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class FlagAFTs
*
* Please fill out the descripton and author information.
*
* $Author: M.O'Flaherty $
* $Date: 07/08/2023 $
* $Contact: moflaher@km.icrr.u-tokyo.ac.jp
*/

class FlagAFTs: public Tool {
	
	public:
	
	FlagAFTs();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	
};


#endif
