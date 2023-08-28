/* vim:set noexpandtab tabstop=4 wrap */
#ifndef LoadFileList_H
#define LoadFileList_H

#include <string>
#include <iostream>

#include "Tool.h"

class TApplication;

/**
* \class LoadFileList
*
* A simple tool to load a list of files for reading.
*
* $Author: M.O'Flahery $
* $Date: 2020/08/12 $
* Contact: marcus.o-flaherty@warwick.ac.uk
*/
class LoadFileList: public Tool {
	
	public:
	
	LoadFileList(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool purpose.
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	private:
	int GetFileList();
	
	std::string inputFile="";
	std::string inputDirectory=".";
	std::string fileList="";
	std::vector<std::string> list_of_files;
	std::string filePattern="";
	bool useRegex=false;
	std::string FileListName="InputFileList";
	
};


#endif
