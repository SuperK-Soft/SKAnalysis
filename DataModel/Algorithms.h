/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include "OutputRedirector.h"   // for CStdoutRedirector

#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>   // for ofstream

#include "basic_array.h"
#include "OutputRedirector.h"   // for CStdoutRedirector

// TODO what's the best way of pulling these both in without introducing circular dependencies?
#include "SK_helper_functions.h"

class TVector3;
class TLorentzVector;
class TPie;
class TH1F;

int ReadListFromFile(std::string filename, std::vector<std::string> &lines, char commentchar='#', bool trim_whitespace=true);
std::string GetStdoutFromCommand(std::string cmd, int bufsize=500);
int SystemCall(std::string cmd, std::string& errmsg);  // or use this one, maybe better?
void SetRootColourPlotStyle();
double MomentumToEnergy(basic_array<float[3]>& mom, int pdg);
double MomentumToEnergy(TVector3& mom, int pdg);
double Mag2(basic_array<float>& mom);
double Mag(basic_array<float>& mom);
bool CheckPath(std::string path, std::string& type);
std::string ToLower(std::string astring);
void PrintObjectTable();
int safeSystemCall(std::string cmd);
int safeSystemCallVerbose(std::string cmd);
void PrintVector(TVector3& avec, bool newline=false);
void PrintVector(TLorentzVector& avec, bool newline=false);
bool IsStlContainer(std::string type_as_string);
std::string toString(const TVector3& vec);
std::string toString(const TLorentzVector& vec);
std::unique_ptr<TPie> GeneratePieFromHisto(TH1F* histo, int verbose=0);
std::unique_ptr<TPie> GeneratePieFromHisto(std::string histoname, int verbose=0);

namespace algorithms{
	
} // end namespace algorithms

// helper function: to_string with a precision
// particularly useful for printing doubles and floats in the Log function
template <typename T>
std::string toString(const T a_value, const int n = 2){
	std::ostringstream out;
	out.precision(n);
	out << std::fixed << a_value;
	return out.str();
}

template <typename T>
std::string toString(T* a_ptr, bool deref=false){
	if(!deref){
		// print the address
		std::stringstream out;
		out<<a_ptr;
		return out.str();
	}
	// try to print the object
	if(a_ptr==nullptr) return "?";
	std::stringstream out;
	out<<*a_ptr;
	return out.str();
}

template <typename T>
bool solveQuadratic(const T &a, const T &b, const T &c, T &x0, T &x1){
	T discr = b*b - 4*a*c;
	if(discr < 0){ return false; }
	else if(discr==0){ x0=-0.5*(b/a); x1=x0; }
	else {
		T q = (b>0) ? -0.5*(b+sqrt(discr)) : -0.5*(b-sqrt(discr));
		x0 = q/a;
		x1 = c/q;
	}
	if (x0>x1) std::swap(x0, x1);
	
	return true;
}

template <typename T>
void printVals(T &container, int messagelevel, int verbosity, std::string preamble="", std::string postamble=""){
	// if messagelevel==0 (error), print to cerr, otherwise print to cout
	std::ofstream outputbuf;
	if(messagelevel==0){
		outputbuf.copyfmt(std::cerr);
		outputbuf.clear(std::cerr.rdstate());
		outputbuf.basic_ios<char>::rdbuf(std::cerr.rdbuf());
	} else {
		outputbuf.copyfmt(std::cout);
		outputbuf.clear(std::cout.rdstate());
		outputbuf.basic_ios<char>::rdbuf(std::cout.rdbuf());
	}
	outputbuf<<preamble<<" {";
	// the below is ok for vectors and arrays, but won't work for maps
	//for(auto it=container.begin(); it!=container.end(); ++it){ outputbuf<<(*it); }
	// instead use a templated version
	if(messagelevel) print_container_cout(container);
	else             print_container_cerr(container);
	outputbuf<<"} "<<postamble<<std::endl;
}

// generic container printer, from https://eli.thegreenplace.net/2014/variadic-templates-in-c/
template <template <typename, typename...> class ContainerType,
          typename ValueType, typename... Args>
void print_container(const ContainerType<ValueType, Args...>& c) {
	for (const auto& v : c) {
		std::cout << v << ' ';
	}
	std::cout << '\n';
}

// Implement << for pairs: this is needed to print out mappings where range
// iteration goes over (key, value) pairs.
template <typename T, typename U>
std::ostream& operator<<(std::ostream& out, const std::pair<T, U>& p) {
  out << "[" << p.first << ", " << p.second << "]";
  return out;
}

// using https://github.com/zpasternack/Redirector instead
// version that takes no arguments
template<typename FunctionType>
std::string getOutputFromFunctionCall(FunctionType&& FunctionToCall){
	// set up capturing of stdout, capturing their current nominal outputs
	CStdoutRedirector theRedirector;
	theRedirector.StartRedirecting();
	// invoke the requested function
	std::forward<FunctionType>(FunctionToCall)();
	// restore the streams, capturing the intermediate buffer
	std::string captured_out = theRedirector.GetOutput();
	theRedirector.ClearOutput();
	theRedirector.StopRedirecting();
	
	// FIXME using std::cout and std::endl in prints, or e.g. 'system(...)' calls
	// still seem to have the ToolAnalysis '[1]: ' prefix in printouts.
	// Figure out how to prevent this, the following is not a good solution.
	// find and remove these
	for(int i=0; i<5; ++i){
		std::string substr="["+toString(i)+"]: ";
		size_t pos = captured_out.find(substr);
		while(pos!=std::string::npos){
			captured_out.erase(pos,substr.length());
			pos = captured_out.find(substr);
		}
	}
	// trim any trailing newlines
	while(captured_out.back()=='\n'||captured_out.back()=='\r') captured_out.pop_back();
	
	return captured_out;
}

// version that accepts arbitrary arguments?
template<typename FunctionType, typename... Rest>
std::string getOutputFromFunctionCall(FunctionType&& FunctionToCall, Rest... rest){
	// set up capturing of stdout, capturing their current nominal outputs
	CStdoutRedirector theRedirector;
	theRedirector.StartRedirecting();
	// invoke the requested function
	std::forward<FunctionType>(FunctionToCall)(rest...);
	// restore the streams, capturing the intermediate buffer
	std::string captured_out = theRedirector.GetOutput();
	theRedirector.ClearOutput();
	theRedirector.StopRedirecting();
	
	// FIXME prevent, rather than remove (see above version of this function)
	for(int i=0; i<5; ++i){
		std::string substr="["+toString(i)+"]: ";
		size_t pos = captured_out.find(substr);
		while(pos!=std::string::npos){
			captured_out.erase(pos,substr.length());
			pos = captured_out.find(substr);
		}
	}
	// trim any trailing newlines
	while(captured_out.back()=='\n'||captured_out.back()=='\r') captured_out.pop_back();
	
	return captured_out;
}

template <typename T>
bool IsStlContainer(const std::type_info &type){
	std::string demangled = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);
	if(demangled.length()<5 || demangled.substr(0,5)!="std::") return false;
	return IsStlContainer(demangled.substr(5,std::string::npos));
}

#endif
