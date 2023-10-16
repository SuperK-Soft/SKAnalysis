#ifndef WriteSkEvent_H
#define WriteSkEvent_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class WriteSkEvent
*
* make a new TTree entry in an SKROOT file by invoking skroot_fill_tree, and possibly associated functions.
*
* $Author: M.O'Flaherty $
* $Date: 2023/09/06 $
* $Contact: moflaher@km.icrr.u-tokyo.ac.jp
*/

class WriteSkEvent: public Tool {
	
	public:
	
	WriteSkEvent();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	std::string treeReaderName;
	int LUN=0;
	bool delete_outside_hits=true;
	bool require_save_flag=false;
};


#endif
