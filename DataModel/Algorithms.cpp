/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#include "Algorithms.h"
#include "Constants.h"
//#include <libgen.h>  // dirname and basename
#include <sys/stat.h>  // dirname and basename
#include <sys/types.h> // for stat() test to see if file or folder
#include <unistd.h>    // for system()
#include <sys/wait.h>  // for wait()
#include <errno.h>     // for WEXITSTATUS etc
#include <stdlib.h>    // for fork()?
//#include <memory>
//#include <exception>
//#include <cstring>  // strncpy
#include <fstream>
#include <sstream>
#include <cctype> // ::tolower
#include <cxxabi.h>  // demangle

#include "TStyle.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObjectTable.h"

int ReadListFromFile(std::string filename, std::vector<std::string> &lines, char commentchar, bool trim_whitespace){
	// read each new line into a std::vector<string> and return
	std::ifstream fin (filename.c_str());
	// return if not found or can't be opened
	if(not (fin.is_open() && fin.good())){
		return -1;
	}
	std::string Line;
	// loop over file lines
	while (getline(fin, Line)){
		if (Line[0] == commentchar){
			continue;
		} else{
			// check for trailing comments
			if (Line.find(commentchar) != std::string::npos){
				Line.erase(Line.begin()+Line.find(commentchar), Line.end());
			}
			// trim trailing whitespace
			if(trim_whitespace){
				Line.erase(Line.find_last_not_of(" \t\n\015\014\013")+1);
			}
			// add to vector
			lines.push_back(Line);
		}
	}
	fin.close();
	// return number of lines added
	return lines.size();
}

std::string GetStdoutFromCommand(std::string cmd, int bufsize){
	/*
	  credit: Jeremy Morgan, source:
	  https://www.jeremymorgan.com/tutorials/c-programming/how-to-capture-the-output-of-a-linux-command-in-c/
	*/
	std::string data;
	FILE * stream;
	char* buffer = new char[bufsize];
	cmd.append(" 2>&1");
	
	stream = popen(cmd.c_str(), "r");
	if(stream){
		while(!feof(stream)){
			if (fgets(buffer, bufsize, stream) != NULL) data.append(buffer);
		}
		pclose(stream);
	}
	delete[] buffer;
	return data;
}

void SetRootColourPlotStyle(){
	const int NRGBs = 5;
	const int NCont = 255;
	
	double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	double red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

void PrintVector(TVector3& avec, bool newline){
	std::string nl = (newline) ? "\n" : "";
	std::cout<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<")"<<nl;
}

void PrintVector(TLorentzVector& avec, bool newline){
	std::string nl = (newline) ? "\n" : "";
	std::cout<<"("<<avec.X()<<", "<<avec.Y()<<", "<<avec.Z()<<avec.T()<<")"<<nl;
}

double MomentumToEnergy(basic_array<float>& mom, int pdg){
	double mass = PdgToMass(pdg);
	if(mass<0) return -1;
	double momsq = pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2);
	return sqrt(momsq+pow(mass,2));
}

double MomentumToEnergy(TVector3& mom, int pdg){
	double mass = PdgToMass(pdg);
	if(mass<0) return -1;
	double momsq = pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2);
	return sqrt(momsq+pow(mass,2));
}

double Mag2(basic_array<float>& mom){
	return pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2);
}

double Mag(basic_array<float>& mom){
	return sqrt(Mag(mom));
}

bool CheckPath(std::string path, std::string& type){
	struct stat s;
	if(stat(path.c_str(),&s)==0){
		if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory?? how does this work?
			type="d";  //it's a directory
			return true;
		} else if(s.st_mode & S_IFREG){ // mask to check if it's a file??
			type="f"; //it's a file
			return true;
		} else {
			// exists, but neither file nor directory?
			type="???";
			return false;
			//assert(false&&"Check input path: stat says it's neither file nor directory..?");
		}
	} else {
		// does not exist - could be a pattern, e.g. "/path/to/rootfiles_*.root"
		type="none";
		return false;
	}
	return false;
}

std::string ToLower(std::string astring){
	std::transform(astring.begin(), astring.end(), astring.begin(), ::tolower);  // why is the :: needed?
	return astring;
}

void PrintObjectTable(){
	// should work? will gObjectTable be recognised in the invoked context?
//	std::string objecttable = getOutputFromFunctionCall([](){gObjectTable->Print();});
	// safer?
	auto&& printobjtable = [](TObjectTable* gobjs){gobjs->Print();};
	std::string objecttable = getOutputFromFunctionCall(printobjtable, gObjectTable);
	// or could we even obtain a pointer to gObjectTable::Print(), and just pass that?
	
	// anyway. Scan through all the output and extract just the line relating to number of
	// TObjArrays, so we don't swamp the output.
	std::stringstream objectstream(objecttable);
	std::string aline;
	while(getline(objectstream,aline, '\n')){
		if(aline.find("TObjArray")!=std::string::npos) break;
	}
	std::cout<<aline<<std::endl;
}

int safeSystemCall(std::string cmd){
	int ret=-1;
	if (fork() == 0){
		system(cmd.c_str());
		exit(1);
	} else {
		int stat;
		wait(&stat);
		//ret = stat;
		ret = WIFEXITED(stat) ? WEXITSTATUS(stat) : WTERMSIG(stat);
	}
	return ret;
}

int safeSystemCallVerbose(std::string cmd){
	std::cout<<"getting return from '"+cmd+"'"<<std::endl;
	int ret=-1;
	// afaict this shouldn't be necessary when using 'system', rather than 'exec*' directly.
	// system is supposed to already handling forking and waiting. However, for some reason
	// this was not working correctly in its usage
	// (captured stdout was empty for a "readlink" call... maybe a synchronization issue?)
	// using fork and wait seemed to make it work, but we should have been using popen!
	pid_t pid = fork();
	if ( pid == 0){
		// child process.
		system(cmd.c_str());  // run the command.
		_exit(1);             // exit the child process - but do not close parent file descriptors
		                      // (i.e. call `_exit()`, not `exit()`)
	} else {
		// parent process. wait for the child to terminate.
		int stat;
		wait(&stat);
		// check the termination status
		ret = WIFEXITED(stat) ? WEXITSTATUS(stat) : WTERMSIG(stat);
		if(WIFEXITED(stat)){
			// returned "normally"
			int retval = WEXITSTATUS(stat);
			printf("Exit status: %d\n", retval);  // print return val
			// if not succssful, maybe it set a relevant error in errno?
			// print out the corresponding error based on the current value of errno
			if(retval!=0){
				perror(cmd.c_str());
				// prints its argument, then the descrption associated with current errno value
			}
		} else if(WIFSIGNALED(stat)){
			// returned due to receipt of a signal - print the received signal
			psignal(WTERMSIG(stat), "Exit signal");
		} // or many others, see http://man.yolinux.com/cgi-bin/man2html?cgi_command=waitpid(2)
	}
	return ret;
}

// FIXME could replace this with TClassEdit::STLKind(const char* type_as_string) which does the same.
// But ROOT 5 does not know about c++11 container types, so better to use this.
// may also want to use TClassEdit::ShortType(const char* type_as_string) to discard qualifiers
bool IsStlContainer(std::string type_as_string){
	size_t base_pos = type_as_string.find('<');
	if(base_pos==std::string::npos) return false;
	return (container_types.count(type_as_string.substr(0,base_pos)));
}
