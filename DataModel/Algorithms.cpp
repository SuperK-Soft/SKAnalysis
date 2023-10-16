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
#include <locale>      // std::isspace
#include <memory>
#include <string.h>
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
#include "TLegend.h"
#include "TH1.h"
#include "TPie.h"
#include "TROOT.h"

std::string toString(const std::string s){ return s; }

std::string toString(const EventType& ev){
	std::stringstream ss;
	ss << ev;
	return ss.str();
}

std::string toString(const TVector3& vec){
	std::string s = "("+toString(vec.X())+", "+toString(vec.Y())+", "+toString(vec.Z())+")";
	return s;
}

std::string toString(const TLorentzVector& vec){
	std::string s = "("+toString(vec.T())+", "+toString(vec.X())+", "+toString(vec.Y())+", "+toString(vec.Z())+")";
	return s;
}

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
		if(Line.empty()) continue;
		if (Line[0] == commentchar){
			continue;
		} else{
			// check for trailing comments
			if (Line.find(commentchar) != std::string::npos){
				Line.erase(Line.find(commentchar));
				if(Line.empty()) continue;
			}
			// trim trailing whitespace
			if(trim_whitespace){
				Line.erase(Line.find_last_not_of(" \t\n\015\014\013")+1);
				if(Line.empty()) continue;
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
	// pop off trailing newline
	if(std::isspace(data.back())) data.pop_back();
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
	// remove any potential std:: prefix
	if(type_as_string.substr(0,5)=="std::") type_as_string.erase(0,5);
	size_t base_pos = type_as_string.find('<');
	if(base_pos==std::string::npos) return false;
	return (container_types.count(type_as_string.substr(0,base_pos)));
}


std::unique_ptr<TPie> GeneratePieFromHisto(std::string histoname, int verbose){
	TH1F* histo = (TH1F*)gROOT->FindObject(histoname.c_str());
	if(histo==nullptr){
		std::cerr<<"GeneratePieFromHisto could not find histo "<<histoname<<std::endl;
		return nullptr;
	}
	return GeneratePieFromHisto(histo, verbose);
}

std::unique_ptr<TPie> GeneratePieFromHisto(TH1F* histo, int verbose){
	// I may not have realised TPie::TPie(const TH1* h) exists ....?
	std::string histoname = std::string(histo->GetName());
	if(verbose) std::cout<<"creating pie chart from histo "<<histoname<<", which has "
						 <<histo->GetNbinsX()<<" bins with contents: "<<std::endl;
	std::vector< std::pair<std::string,float> > histbins;
	for(int bini=0; bini<histo->GetNbinsX(); bini++){
		TString binlabel = histo->GetXaxis()->GetBinLabel(bini+1);
		float binconts = histo->GetBinContent(bini+1);
		if(binconts<0.01) binconts = 0.0f;  // round floats. useful if the histo has been scaled.
		if(verbose && binconts!=0.0f) std::cout<<binlabel.Data()<<" : "<<binconts<<std::endl;
		if(binconts<0) std::cerr<<"error converting "<<histoname<<" to pie chart: bin "<<binlabel.Data()
			<<" has "<<binconts<<" entries!"<<std::endl;
		if(binconts!=0) histbins.emplace_back(binlabel.Data(),binconts);
	}
	
	auto thepie = std::unique_ptr<TPie>(new TPie(TString::Format("%sPie",histoname.c_str()), TString::Format("%s",histoname.c_str()), histbins.size()));
	
	for(size_t bini=0; bini<histbins.size(); bini++){
		std::pair<std::string,float> abin = histbins.at(bini);
		std::string thebinlabel = abin.first;
		float thebincontents = abin.second;
		if(thebincontents<0) std::cerr<<"error converting "<<histoname<<" to pie chart: bin "<<thebinlabel
			<<" has "<<thebincontents<<" entries!"<<std::endl;
		thepie->SetEntryVal(bini,thebincontents);  // NO +1 - TPie's have no underflow bin!
		thepie->SetEntryLabel(bini,thebinlabel.c_str());
	}
	return thepie;
}

int SystemCall(std::string cmd, std::string& retstring){
	// execute command, capture return status and any output
	
	int retstatus;  // exit code of command
	retstring="";
	
	// gotcha: our particular command redirects stdout,
	// so if we append 2>&1 we get stderr into our output file too!
	// fortunately we can put redirects whereever we like,
	// so insert it before anything else
	cmd.insert(0,"2>&1; ");
	cmd.append(";");
	
	// fork, execute cmd with '/bin/sh -c' and open a read pipe connected to outout (r)
	FILE * stream = popen(cmd.c_str(), "r");
	// check the fork succeeded and we got a pipe to read from
	if(stream){
		// read any return
		unsigned int bufsize=255;
		char buffer[bufsize]; // temporary buffer for reading return
		while(!feof(stream)){
			// if we're very paranoid we can check for read error with (ferror(stream)!=0)
			if (fgets(buffer, bufsize, stream) != NULL) retstring.append(buffer);
		}
		// pop off trailing newline
		if(std::isspace(retstring.back())) retstring.pop_back();
		// close the pipe, and capture command return value
		int stat = pclose(stream);
		retstatus = (WIFEXITED(stat)) ? WEXITSTATUS(stat) : WTERMSIG(stat);
	} else {
		retstring = "SystemCall popen returned nullptr for command '"+cmd
		          +"'', error is: "+strerror(errno);
		return -1;
	}
	if(retstatus!=0){
		retstring="SystemCall with command "+cmd+" failed with error code "
		         +std::to_string(retstatus)+", stderr returned '"+retstring+"'";
	}
	
	return retstatus;
	
}
