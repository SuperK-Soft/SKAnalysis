/* vim:set noexpandtab tabstop=4 wrap */
#include "LoadFileList.h"
#include "FindFilesInDirectory.h"

#include "type_name_as_string.h"

LoadFileList::LoadFileList():Tool(){
	toolName=type_name<decltype(this)>(); toolName.pop_back();
}


bool LoadFileList::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	m_data= &data;
	
	// Get the Tool configuration variables
	// ====================================
	m_variables.Get("verbosity",verbosity);           // how verbose to be
	m_variables.Get("inputFile",inputFile);           // a single specific input file
	m_variables.Get("inputDirectory",inputDirectory); // input files dir
	m_variables.Get("fileList",fileList);             // a file containing a list of input filenames
	m_variables.Get("filePattern",filePattern);       // a pattern to match files in input directory
	m_variables.Get("useRegex",useRegex);             // is the pattern a glob or a regex
	m_variables.Get("FileListName",FileListName);     // what key to use to store the list in the CStore
	
	Log(toolName+": Initializing",v_debug,verbosity);
	
	// get list of files to process.
	int num_files = GetFileList();
	Log(toolName+" loaded "+toString(num_files)+" files to read",v_debug,verbosity);
	
	// set the files into the CStore
	m_data->CStore.Set(FileListName, list_of_files);
	
	return true;
}

int LoadFileList::GetFileList(){
	if(inputFile!=""){
		Log(toolName+" will process file "+inputFile,v_message,verbosity);
		list_of_files.emplace_back(inputFile);
		return 1;
	} else if(fileList!=""){
		Log(toolName+" loading list of input files from "+fileList,v_debug,verbosity);
		// user has provided a list of files. Load the list.
		int num_files = ReadListFromFile(fileList, list_of_files);
		// Check the first file to see if we have just filenames, or full paths
		std::string first_file = list_of_files.front();
		bool abs_path = (first_file.find('/')!=std::string::npos);
		
		// if we got the path separately, prepend it
		if(!abs_path){
			// check that we have a directory
			if(inputDirectory==""){
				Log(toolName+" Error! Input list of files does not contain absolute paths,"
					+" but input directory is not given!",v_error,verbosity);
				return -1;
			}
			// if ok, build absolute paths
			for(auto&& afile : list_of_files) afile = inputDirectory + "/" + afile;
		}
		return list_of_files.size();
	}
	
	// otherwise search a directory for files matching a given pattern.
	// quick sanity check that we have a directory
	if(inputDirectory==""){
		Log(toolName+" Error! Using input file pattern but input directory is not given!",v_error,verbosity);
		return -1;
	}
	Log(toolName+" searching for files in "+inputDirectory
		+" that match pattern '"+filePattern+"'",v_debug,verbosity);
	
	// This function is in DataModel/FindFilesInDirectory.cpp
	// int FindFilesInDirectory(std::string inputdir, std::string pattern, std::vector<std::string> &matches, bool case_sensitive=false, int max_subdir_depth=0, bool use_regex=false, std::vector<std::string>* filenames=nullptr, std::vector<std::vector<std::string>>* output_submatches=nullptr, bool verbose=false);
	bool case_sensitive = false;  // case insensitive (default false)
	int max_subdir_depth=1;       // this dir only, 0: all the way down (default 0)
	int num_files = FindFilesInDirectory(inputDirectory, filePattern, list_of_files, case_sensitive, max_subdir_depth, useRegex, nullptr, nullptr, (verbosity>=v_debug));
	
	return num_files;
}


bool LoadFileList::Execute(){
	return true;
}


bool LoadFileList::Finalise(){
	return true;
}
