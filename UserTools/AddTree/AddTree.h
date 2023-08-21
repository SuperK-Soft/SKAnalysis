#ifndef AddTree_H
#define AddTree_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
* \class AddTree
*
* Create a TTree with a TreeManager using skroot_open_write, but associate it with an existing file.
*
* $Author: M.O'F $
* $Date: 2023/08/15 $
*/

class TTree;
class MTreeReader;

class AddTree: public Tool {
	
	public:
	AddTree();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	MTreeReader* myTreeReader=nullptr;
	TTree* thistree=nullptr;
	TFile* ofile=nullptr;
};


#endif
