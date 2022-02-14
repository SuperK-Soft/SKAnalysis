#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
//#pragma link C++ nestedtypedefs;
//#pragma link C++ class BStore::*
//#pragma link C++ class BinaryStream::*
//#pragma link C++ class PointerWrapper::*
//#pragma link C++ class SerialisableObject::*
#include "MyClass.h"
#pragma extra_include "MyClass.h"
#include "TVector3.h"
#pragma extra_include "TVector3.h"
#include "TLorentzVector.h"
#pragma extra_include "TLorentzVector.h"
#pragma link C++ defined_in "BStore.h";
#pragma link C++ function BStore::Set(std::string, int);
#pragma link C++ function BStore::Get(std::string, int&);
#pragma link C++ function BStore::Set(std::string, double);
#pragma link C++ function BStore::Get(std::string, double&);
#pragma link C++ function BStore::Set(std::string, std::string);
#pragma link C++ function BStore::Set(std::string, std::string*, bool);
#pragma link C++ function BStore::Get(std::string, std::string&);
#pragma link C++ function BStore::Get(std::string, std::string*&);
#pragma link C++ function BStore::Set(std::string, TString);
#pragma link C++ function BStore::Set(std::string, TString*, bool);
#pragma link C++ function BStore::Get(std::string, TString&);
#pragma link C++ function BStore::Get(std::string, TString*&);
#pragma link C++ function BStore::Set(std::string, std::vector<std::string, std::allocator<std::string> >&);
#pragma link C++ function BStore::Set(std::string, std::vector<std::string, std::allocator<std::string> >*&, bool);
#pragma link C++ function BStore::Get(std::string, std::vector<std::string, std::allocator<std::string> >&);
#pragma link C++ function BStore::Get(std::string, std::vector<std::string, std::allocator<std::string> >*&);
#pragma link C++ function BStore::Set(std::string, TLorentzVector);
#pragma link C++ function BStore::Set(std::string, TLorentzVector*, bool);
#pragma link C++ function BStore::Get(std::string, TLorentzVector&);
#pragma link C++ function BStore::Get(std::string, TLorentzVector*&);
#pragma link C++ function BStore::Set(std::string, TVector3);
#pragma link C++ function BStore::Set(std::string, TVector3*, bool);
#pragma link C++ function BStore::Get(std::string, TVector3&);
#pragma link C++ function BStore::Get(std::string, TVector3*&);
#pragma link C++ function BStore::Set(std::string, MyClass);
#pragma link C++ function BStore::Set(std::string, MyClass*, bool);
#pragma link C++ function BStore::Get(std::string, MyClass&);
#pragma link C++ function BStore::Get(std::string, MyClass*&);
