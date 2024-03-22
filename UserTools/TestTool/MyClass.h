#ifndef MYCLASS_H
#define MYCLASS_H

#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>

//#include "BStore.h"

class MyClass {
	public:
	MyClass(){ internal.str(""); /*bs=0;*/ };
	//~MyClass(){};   // in general don't define one if you don't need to (rule of 3)
	// https://en.cppreference.com/w/cpp/language/rule_of_three
	// because we have a stringstream object we must manually define the copy constructor
	// as copy construction of stringstream is not defined
	MyClass(const MyClass& other): internal(other.internal.str()) {}
	
	// XXX even though it isn't used, both arguments MUST be given a name,
	// otherwise this fails when called from gInterpreter with
	// "symbol 'in' not defined" - so first arg is now 'name'.
	template<typename T> bool Set(std::string name, T in){
		std::cout<<"MyClass::Set with type "<<typeid(T).name()<<std::endl;
		std::cout<<"setting internal to "<<in<<std::endl;
		internal.str(""); internal.clear();
		internal << in;
		std::cout<<"internal is now "<<internal.str()<<std::endl;
		//if(bs) bs->Set(name,in);
		return true;
	}
	
	template<typename T> bool Get(std::string name, T& out){
		std::cout<<"MyClass::Get with type "<<typeid(T).name()<<std::endl;
		 std::cout<<"passing out internal "<<internal.str()<<std::endl;
		 internal >> out;
		 std::cout<<"out is now "<<out<<std::endl;
		//if(bs) bs->Get(name,out);
		return true;
	}
	
//	void SetBsPtr(void* ptr){
//		std::cout<<"setting bs to "<<ptr<<std::endl;
//		bs=(BStore*)ptr;
//	}
	
	bool Print(){
		std::cout<<"internal is currently "<<internal.str()<<std::endl;
		//std::cout<<"bs is currently "<<bs<<std::endl;
		return true;
	}
	
	std::stringstream internal;
	//BStore* bs;
};

//template<>
//bool MyClass::Get<char>(char &out){
//	std::cout << "Get specialization for void type" << std::endl;
//	return false;
//}

//template<>
//bool MyClass::Set<char>(char in){
//	std::cout << "Set specialization for void type" << std::endl;
//	return false;
//}

#endif

//#ifdef __CINT__
//  #pragma link C++ function MyClass::Set(std::string, int);
//  #pragma link C++ function MyClass::Set(std::string, double);
//  #pragma link C++ function MyClass::Get(std::string, int&);
//  #pragma link C++ function MyClass::Get(std::string, double&);
//#else
//  template bool MyClass::Set(std::string, int);
//  template bool MyClass::Set(std::string, double);
//  template bool MyClass::Get(std::string, int&);
//  template bool MyClass::Get(std::string, double&);
//#endif


