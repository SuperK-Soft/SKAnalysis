#ifdef __CINT__
  #pragma link off all globals;
  #pragma link off all classes;
  #pragma link off all functions;
  #pragma link C++ nestedclasses;
  #pragma link C++ defined_in "MyClass.h";
//  #pragma link C++ function MyClass::Set(std::string, int);
//  #pragma link C++ function MyClass::Set(std::string, double);
//  #pragma link C++ function MyClass::Get(std::string, int&);
//  #pragma link C++ function MyClass::Get(std::string, double&);
  #pragma link C++ function MyClass::Set(std::string, std::string);
  #pragma link C++ function MyClass::Get(std::string, std::string&);
  #pragma link C++ function MyClass::Set(std::string, int);
  #pragma link C++ function MyClass::Get(std::string, int&);
  #pragma link C++ function MyClass::Set(std::string, int);
  #pragma link C++ function MyClass::Get(std::string, int&);
#endif
