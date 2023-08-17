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
* $Author: ?.????? $
* $Date: ????/??/?? $
* $Contact: ???@km.icrr.u-tokyo.ac.jp
*/

class WriteSkEvent: public Tool {
	
	public:
	
	WriteSkEvent();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	
	int LUN=0;
	bool delete_outside_hits=true;
};


#endif
