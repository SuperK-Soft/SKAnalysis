#include "TestTool.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "MyClass.h"

TestTool::TestTool():Tool(){}

// probably the start of trying to call a global function from TCint via Reflex... unused
template<typename T>
void* BSGet(void* storep, std::string name, T &obj){
   BStore* store = (BStore*)(storep);
   store->Get(name, obj);
   void* retp = (void*)(&obj);
   return retp;
}

bool TestTool::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  
  std::cout<<"making root stuff"<<std::endl;
  f = new TFile("trashme.root","recreate");
  t = new TTree("testtree","testtree");
  t2 = new TTree("t2","t2");    // CHECK, WE MAY ALSO ONLY SUPPORT TREES ON THE HEAP, NOT THE STACK
  // (it segfaults if the TInterpreter tries to modify stack memory because it basically is
  // one application trying to modify memory space of another application...)
  
  std::cout<<"making bstore"<<std::endl;
  BStore* mystore = new BStore(true); // IMPORTANT: AGAIN THIS ONLY WORKS ON STORES ON THE HEAP, NOT THE STACK
  //system("rm trash");
  //mystore->Initnew("trash", uncompressed, 1);
  int myint = 555;
  double mydouble = 44.4;
//  TString mystring = "potato";  // TTree has issues with TString??? a ROOT issue not us.
  std::string mystring = "potato";
  TVector3 myvec{1,2,3};
  std::vector<std::string> mystrings{"potato","salad"};
  TH1D myhist("myhist","title",100,0,10);
  myhist.FillRandom("gaus",100);
  TLorentzVector mylorentz(1,2,3,4);
  MyClass myclass;
//  MCInfo* mci = new MCInfo();  // would need a ROOT dictionary for MCInfo: i think there is one? add to Makefile
  mystore->Set("myint",myint);
  mystore->Set("mydouble",mydouble);
  mystore->Set("mystring",mystring);
  mystore->Set("myvec",myvec);
  mystore->Set("mystrings",mystrings);
  mystore->Set("mylorentz",mylorentz);
  mystore->Set("myclass",myclass);
//  mystore->Set("mci",*mci);
  //mystore->Set("myhist",myhist);
  // ok TH1D actually crashes on TTree->Show(0)... which i'm going to assume is a ROOT problem.
  // tbh it's unlikely someone will want to put a histogram in a ROOT branch anyway

//  std::cout<<"making store a"<<std::endl;
//  Store mystore2;
//  int mint = 22;
//  double mdub = 33.3;
//  std::string ms = "carrot";
//  mystore2.Set("mint",mint);
//  mystore2.Set("mdub",mdub);
//  mystore2.Set("ms",ms);
//  
//  std::cout<<"making store b"<<std::endl;
//  Store mystore3;
//  int ii = 99;
//  float ff=12.3f;
//  std::string ss = "store3";
//  mystore3.Set("ii",ii);
//  mystore3.Set("ff",ff);
//  mystore3.Set("ss",ss);
//  
//  std::cout<<"setting stores a and b into datamodel"<<std::endl;
//  m_data->tool_configs["store_a"] = &mystore2;
//  m_data->tool_configs["store_b"] = &mystore3;
  
  /*
  std::cout<<"pointed getter of tstring"<<std::endl;
  TString* ptr;
  mystore->Get("mystring",ptr);
  std::cerr<<"ptr points to "<<ptr<<std::endl;
  std::cerr<<"this holds "<<ptr->Data()<<std::endl;

  gSystem->Load("libBStore_RootDict.so");
  //gSystem->Load("lib/libBStore_RootDict.so");
  TString cmd = TString::Format("BStore* bs = (BStore*)0x%x; ", mystore);
  gInterpreter->ProcessLine(cmd.Data());
  std::cout<<"getting mystring"<<std::endl;
//  cmd = "TString* tsp; bs->Get(\"mystring\",tsp);";
//  gInterpreter->ProcessLine(cmd.Data());
//  gInterpreter->ProcessLine("printf(\"tsp points to %s\\n\",tsp->Data());");
//  cmd = "TString ts; bs->Get(\"mystring\",ts);";
//  gInterpreter->ProcessLine(cmd.Data());
//  gInterpreter->ProcessLine("printf(\"tsp reads %s\\n\",tsp.Data());");
  cmd = "int iss; bs->Get(\"myint\",iss);";
  gInterpreter->ProcessLine(cmd.Data());
  gInterpreter->ProcessLine("printf(\"iss reads %d\\n\",iss);");
  
  return true;
*/
  
  std::cout<<"calling MakeBranches on Bstore"<<std::endl;
  //MakeBranches(t,mystore);
  m_data->StoreConverter.MakeBranches(t, mystore);
//  std::cout<<"calling MakeBranches on store map"<<std::endl;
//  MakeBranches(t2,m_data->tool_configs);
  
  std::cout<<"tree 1 structure is:"<<std::endl;
  t->Print();
//  std::cout<<"tree 2 structure is:"<<std::endl;
//  t2->Print();
  
  //std::cout<<"calling SetBranches on Bstore"<<std::endl;
  //SetBranches(t,mystore);
  //std::cout<<"calling SetBranches on store map"<<std::endl;
  //SetBranches(t2,m_data->tool_configs);
  
  std::cout<<"calling FillBranches on Bstore"<<std::endl;
  //FillBranches(t,mystore);
  m_data->StoreConverter.FillBranches(t, mystore);
  
  std::cout<<"entry count: "<<t->GetEntries()<<std::endl;
  std::cout<<"entry 0 is: "<<std::endl;
  t->Show(0);
  
  TVector3* avp = new TVector3();
  std::cout<<"resetting addresses"<<std::endl;
  t->ResetBranchAddresses();
  std::cout<<"setting address"<<std::endl;
  t->SetBranchAddress("myvec",&avp);
  std::vector<std::string> mystrings2;
  std::vector<std::string>* mystrings2p = &mystrings2;
  t->SetBranchAddress("mystrings",&mystrings2p);
  std::cout<<"getting entry 0"<<std::endl;
  t->GetEntry(0);
  std::cout<<"myvec was {"<<avp->X()<<","<<avp->Y()<<","<<avp->Z()<<"}"<<std::endl;
  std::cout<<"mystrings were: {";
  for(auto&& as : mystrings2) std::cout<<as<<", ";
  std::cout<<"\b\b}"<<std::endl;
  
//  std::cout<<"calling FillBranches on store map"<<std::endl;
//  FillBranches(t2,m_data->tool_configs);

//  std::cout<<"entry count: "<<t2->GetEntries();
//  std::cout<<"entry 0 is: "<<std::endl; t2->Show(0);

  delete avp;
  delete mystore;
  
  return true;
}


bool TestTool::Execute(){

  return true;
}


bool TestTool::Finalise(){
  
  std::cout<<"writing file"<<std::endl;
  f->Write();
  std::cout<<"closing file"<<std::endl;
  f->Close();
  std::cout<<"deleting file"<<std::endl;
  delete f;
  
  // this makes the application segfault on exit. While harmless, it looks bad.
  // the culprit line is 'TClass::SetUnloaded' called by ROOT's libCore or something.
  // I think the easiest solution is just to not manually delete the interpreter. ¯\_(ツ)_/¯
  //if(meInterpreter) delete meInterpreter; 
  

  return true;
}

// functions to transfer Stores or BStores into a TTree
void TestTool::MakeBranches(TTree* tree, Store* store) {
    for (auto const& pair: GetMap(store)) {
        auto key = pair.first.c_str();
        if(tree->FindBranch(key)) continue;  // continue if it already exists
        float tmpNum;
        if (store->Get(key, tmpNum)) {
            tree->Branch(key, &tmpNum);
        } else {
            std::string tmpStr;
            tmpStr = pair.second;
            tree->Branch(key, &tmpStr);
        }
    }
}

void TestTool::MakeBranches(TTree* tree, BStore* store) {
    for (auto const& pair: GetMap(store)) {
        auto key = pair.first.c_str();
        std::cout<<"next branch: "<<key<<std::endl;
        if(tree->FindBranch(key)) continue;  // continue if it already exists
        std::cout<<"doesn't exist"<<std::endl;
        
        std::cout<<"store is "<<store<<std::endl;
        std::cout<<"store->Has("<<key<<") = "<<store->Has(key)<<std::endl;
        std::string thetype = store->Type(key);
        std::cout<<"mangled type: "<<thetype<<std::endl;
        thetype = abi::__cxa_demangle(thetype.c_str(), nullptr, nullptr, nullptr);
        std::cout<<"unmangled type: "<<thetype<<std::endl;
        // in order to not split objects, we have to know if it's a basic type
        // and then modify the branching command accordingly
        bool basic_type = (fundamental_types.count(thetype)>0);
        // ROOT 5 doesn't like nested containers - we need a dictionary for them.
        // dictionaries for some basic nested containers of fundamental types
        // are located in libRootStl.so, which we link into libBStore_RootDict.so
        // Users will need to generate their own dictionaries for containers of custom classes,
        // and link those to libBStore_RootDict as well.
        // Even with this, CINT doesn't like the use of 'std::' prefixes in nested containers:
        // root [0] std::vector<std::vector<double>> mydubs;
        // will fail with:
        // Error: Invalid type 'vector<std::vector<double>' in declaration of '>mydubs'
        // while doing:
        // root [0] vector<vector<double>> mydubs;
        // will succeed (provided we have a suitable dictionary entry).
        std::string subs="std::";
        size_t len = subs.length();
        while(thetype.find(subs)!=std::string::npos){ thetype.erase(thetype.find(subs),len); }
        std::cout<<"after std::prefix removals: "<<thetype<<std::endl;
        // also remove qualifiers like 'const' and pointers...but only for the outermost type, if nested.
        // XXX is this correct functionality for pointers?
        thetype = RemoveQualifiers(thetype);
        std::cout<<"after qualifier removals: "<<thetype<<std::endl;
        // if a user tries to store a container of pointers such as a vector<Hit*>
        // we will error. I suppose what the user may *like* would be to have us
        // dereference and store the elements, but that's going a step too far.
        if(thetype.find('*')!=std::string::npos){
            Log(m_unique_name+" Error! Can't make a TBranch for BStore entry "+key
                +" of type "+thetype+" which contains pointers.",v_error,m_verbose);
            continue;
        }
        std::cout<<"no pointers left, basic_type is: "<<basic_type<<std::endl;
        
        // BStores can store many different types.
        // Rather than having a switch statement based on the type string,
        // we can hack it using a TInterpreter, provided ROOT is aware of the type.
        // TODO we could probably load libraries into the interpreter to support
        // e.g. SKROOT class objects as well, if need be.
        TString command = TString::Format("TTree* t = (TTree*)0x%x; ", tree);
        if(basic_type){
            command += thetype + " " + key + "; TBranch* b = t->Branch(\"" + key + "\", &" + key + ");";
        } else {
            // splitlevel->0, do not fragment objects.
            command += thetype + "* "+ key + "=0; TBranch* b = t->Branch(\"" + key + "\", &" + key + ",32000,0);";
        }
        std::cout<<"Making Tree branch with command: "<<std::endl<<command.Data()<<std::endl;
        gROOT->ProcessLine(command.Data());
        // NOTE! the branch address is not yet correct, but the branch is now
        // made with the correct type. We will need to set the address before calling Fill,
        // but this can be set with a void*, so we still don't need to instantiate
        // a suitable object of arbitrary type.
        // I'm not sure if the address of the object in the BStore can be expected
        // to remain constant, so it will be safer to re-set the address before
        // each call to TBranch::Fill. If you're concerned about efficiency,
        // it's probably better to use the TTree in a normal way ^_^;
    }
}

// Stores are ascii based, so basically can only hold
// numerical values or strings. Try to get it into a number,
// otherwise fall back to string.
// TODO improve this by trying to identify e.g. int from float?
void TestTool::SetBranches(TTree* tree, Store* store) {
    for (auto const& pair: GetMap(store)) {
        auto key = pair.first.c_str();
        auto value = pair.second;
        float tmpNum;                       // XXX this function would require a persistent
        if (store->Get(key, tmpNum)) {      // member tmpNum, but it is unused.
            tree->SetBranchAddress(key, &tmpNum);
        } else {
            std::string tmpStr;             // XXX this function would require persistent members
            std::string* tmpStrPtr=&tmpStr; // tmpStr and tmpStrPtr but it is unused
            tmpStr = value;
            tree->SetBranchAddress(key, &tmpStrPtr);
        }
    }
}

// BStores can hold all sorts of things. But unlike Stores,
// which will internally hold a float as a string,
// BStores hold data internally in proper binary.
// It is sufficient therefore just to take the address.
void TestTool::SetBranches(TTree* tree, BStore* store) {
    std::cout<<"BStore SetBranches"<<std::endl;
    for (auto const& pair: GetMap(store)) {
        std::string key = pair.first;
        auto val = pair.second;
        std::string thetype = store->Type(key);
        thetype = abi::__cxa_demangle(thetype.c_str(), nullptr, nullptr, nullptr);
        std::string* val_buffer = &((*store)[key]->buffer);
        void* data_pointer = (thetype=="std::string") ? (void*)val_buffer : (void*)&val_buffer->operator[](0);
        tree->SetBranchAddress(key.c_str(), &data_pointer);
        std::cout<<"data_pointer is "<<data_pointer<<std::endl;
        void** addp = (void**)tree->GetBranch(key.c_str())->GetAddress();
        std::cout<<"branch set up to fill from object at "<<(*addp)<<std::endl;
        std::cout<<"key "<<key<<" holds value ";
        if(key=="myint"){ int* mi = (int*)(data_pointer); std::cout<<(*mi)<<std::endl; }
        if(key=="mydouble"){ double* mi = (double*)(data_pointer); std::cout<<(*mi)<<std::endl; }
        if(key=="mystring"){ TString* mi = (TString*)(data_pointer); std::cout<<(*mi)<<std::endl; }
        //if(key=="mystring"){ std::string* mi = (std::string*)(data_pointer); std::cout<<(*mi)<<std::endl; }
        if(key=="myvec"){ TVector3* mi = (TVector3*)(data_pointer);
        std::cout<<"{"<<mi->X()<<","<<mi->Y()<<","<<mi->Z()<<"}"<<std::endl; }
    }
}

void TestTool::MakeBranches(TTree* tree, std::map<std::string, Store*> &stores)
{
    for(auto&& astore : stores) MakeBranches(tree, astore.second);
}

void TestTool::SetBranches(TTree* tree, std::map<std::string, Store*> &stores)
{
     for(auto&& astore : stores) SetBranches(tree, astore.second);
}

bool TestTool::FillBranches(TTree* tree, std::map<std::string, Store*> &stores)
{
    bool ret=true;
    for(auto&& astore : stores) ret &= FillBranches(tree, astore.second);
    return ret;
}

std::string TestTool::RemoveQualifiers(std::string thetype, bool allpointers){
    std::string underlying_type = thetype;
    // remove pointerage
    if(allpointers){
        // e.g. reduce std::vector<Hit*>* to std::vector<Hit>*
        while(underlying_type.find('*')!=std::string::npos){
            underlying_type.erase(underlying_type.find('*'),1);
        }
    } else {
        // e.g. reduce std::vector<Hit*>* to std::vector<Hit*>
        if(underlying_type.back()=='*') underlying_type.pop_back();
    }
    // remove any const qualifiers
    while(underlying_type.find("const")!=std::string::npos){
        underlying_type.erase(underlying_type.find("const"),5);
    }
    // collapse double spaces
    while(underlying_type.find("  ")!=std::string::npos){
        underlying_type.replace(underlying_type.find(" "),2," ");
    }
    // trim
    underlying_type = underlying_type.substr(underlying_type.find_first_not_of(" "),
                      underlying_type.find_last_not_of(" ")-
                      underlying_type.find_first_not_of(" ")+1);
    return underlying_type;
}

bool TestTool::RegisterTypes(std::string thetype){
    bool no_error=true; // should we return early on error...?
    // check for templated types
    int start = thetype.find('<');
    if(start!=std::string::npos){
        std::string outertype=thetype.substr(0,start); // e.g. "std::map"
        // check if this is an stl container or a custom class we need to register (e.g. Cluster<Hit>)
        if(outertype.length()>5){
            // there might be a 'std::' prefix
            if(outertype.substr(0,5)!="std::"){
                // maybe it's still an stl container, without the prefix...?
                if(container_types.count(outertype)==0){
                    // not a container! register this type
                    no_error &= RegisterType(outertype);
                } // else it's an stl containter just without the 'std::' prefix
            } // else it lives in std namespace - it's an stl container
        } else {
            // it's not long enough to have a 'std::' prefix, but...
            if(container_types.count(outertype)==0){
                // not a container! register this type too
                no_error &= RegisterType(outertype);
            } // else it's an stl containter, just without the 'std::' prefix
        }
        //std::string innertype=thetype.substr(start+1,thetype.length()-start-1);  // e.g "int,std::vector<int>"
        // handle possible trailing spaces
        std::string innertype=thetype.substr(start+1,thetype.find_last_of('>')-start-1);
        // to handle containers of pairs and maps we need to search
        // for a ',' outside a matched pair of '<' and '>'
        if(innertype.find(',')){
            size_t splitpos=std::string::npos;
            int nestlevel=0;
            for(int i=0; i<innertype.length(); ++i){
                if(innertype[i]=='<') ++nestlevel;
                if(innertype[i]=='>') --nestlevel;
                if(innertype[i]==','&&nestlevel==0){
                    splitpos=i; break;
                }
            }
            if(splitpos!=std::string::npos){
                // innertype is a pair of two types! Register both
                no_error &= RegisterTypes(innertype.substr(0,splitpos));
                no_error &= RegisterTypes(innertype.substr(splitpos+1,std::string::npos));
            } else {
                // there are pairs somewhere here, but not at this level
                no_error &= RegisterTypes(innertype);
            }
        } else {
            // no commas, we can be sure innertype is one type
            no_error &= RegisterTypes(innertype);
        }
    } else {
        // no templates here. Check if its a fundamental type
        if(fundamental_types.count(thetype)==0){
            no_error &= RegisterType(thetype);
        } // else nothing to do for fundamental types
    }
    return no_error;
}

bool TestTool::RegisterType(std::string thetype){
    // we need to add two new lines to #include the header defining this class
    // (with the exception of (most?) stl containers).
    // For classes of which ROOT is aware we can determine this from the TClass:
    TClass* tc = TClass::GetClass(thetype.c_str());
    if(tc==nullptr){
        Log(m_unique_name+" Error! Could not build TClass from "+thetype
            +"! Have you built a ROOT dictionary for it? (See Makefile)",
            v_error,m_verbose);
        return false;
    }
    const char* declared_in = tc->GetDeclFileName();
    if(declared_in==nullptr || strcmp(declared_in,"")==0){
        // this could be e.g. because either it's unknown to ROOT
        // or because it's a combination of classes,
        // e.g. std::vector<TString>
        Log(m_unique_name+" Error! Could not determine header "+thetype
            +" resides in! Have you built a ROOT dictionary for it?"
            +" (See Makefile)",
            v_error,m_verbose);
        return false;
    }
    std::string decl_file = declared_in;
    if(decl_file.find("prec_stl")!=std::string::npos){
        std::cout<<"skipping includes for class "<<thetype
                 <<" as we think this is a stl container of fundamentals"<<std::endl;
    } else {
        // ok we have a header
        // strip any preceding path from the header file
        int startpos = decl_file.find_last_of('/');
        startpos = (startpos==std::string::npos) ? 0 : startpos+1;
        std::string headername = decl_file.substr(startpos, std::string::npos);
        std::cout<<thetype<<" is declared in "<<headername<<std::endl;
        // add the lines
        std::string cmd = "sed -i '10i #pragma extra_include \"" + headername
                + "\"' BStore_Linkdef.hh";
        get_ok |= system(cmd.c_str());
        cmd = "sed -i '10i #include \"" + headername + "\"' BStore_Linkdef.hh";
        get_ok |= system(cmd.c_str());
        
        // If it's a user class (not a ROOT class) we ALSO need to link in
        // the library with the class definition when building BStore_RootDict.
        // This is handled in the Makefile by linking the same libraries
        // defining user classes to both main and libBStore_RootDict.so
    } // else it's from stl, we can skip the #includes
    return true;
}

bool TestTool::FillBranches(TTree* tree, Store* store){
    //std::cout<<"calling SetBranches from FillBranches"<<std::endl;
    //SetBranches(tree, store);
    std::cout<<"calling TBranch::Fill"<<std::endl;
    int num_entries=-1;
    for (auto const& pair: GetMap(store)) {
        std::string key = pair.first;
        std::string value = pair.second;
        float tmpNum;
        std::string tmpStr;
        std::string* tmpStrPtr = &tmpStr;
        if(store->Get(key, tmpNum)) {
            std::cout<<"Store saving number "<<tmpNum<<", setting address to "<<&tmpNum<<std::endl;
            tree->SetBranchAddress(key.c_str(), &tmpNum);
            TBranch* p = tree->GetBranch(key.c_str());
            std::cout<<"branch at "<<p<<std::endl;
            void** add = (void**)p->GetAddress();
            std::cout<<" will fill with object at "<<(add)<<std::endl;
        } else {
            tmpStr = value;
            std::cout<<"Store saving string "<<tmpStr<<", setting address to "<<&tmpStrPtr<<std::endl;
            tree->SetBranchAddress(key.c_str(), &tmpStrPtr);
            TBranch* p = tree->GetBranch(key.c_str());
            std::cout<<"branch at "<<p<<std::endl;
            void** add = (void**)p->GetAddress();
            std::cout<<" will fill with object at "<<(*add)<<std::endl;
        }
        std::cout<<"calling fill on branch "<<key<<std::endl;
        tree->GetBranch(key.c_str())->Fill();
        std::cout<<"variable is set"<<std::endl;
        int n_entries = tree->GetBranch(key.c_str())->GetEntries();
        if(num_entries<n_entries) num_entries = n_entries;  // keep track of the max in any branch
        std::cout<<"done"<<std::endl;
    }
    tree->SetEntries(num_entries);
    tree->ResetBranchAddresses();
    return true;
}

bool TestTool::FillBranches(TTree* tree, BStore* store){
    int num_entries=-1;
    std::cout<<"calling TBranch::Fill"<<std::endl;
    for (auto const& pair: GetMap(store)) {
        std::string key = pair.first;
        auto val = pair.second;
        std::string thetype = store->Type(key);
        thetype = abi::__cxa_demangle(thetype.c_str(), nullptr, nullptr, nullptr);
        // have we encountered this type before: if so we can skip
        // the checks about whether the BStore is aware of its type
        bool is_known = (known_types.count(thetype)>0);
        bool basic_type;
        if(is_known){
            basic_type = known_types.at(thetype);
        } else {
            // we've not seen this type before.
            // to be able to check if this is a fundamental type we need to strip qualifiers
            std::string underlying_type = RemoveQualifiers(thetype, true);
            // check if it's a fundamental type
            basic_type = (fundamental_types.count(underlying_type)!=0);
            // record our results for the future
            known_types.emplace(thetype, basic_type);
            underlying_types.emplace(thetype, underlying_type);
        }
        if(basic_type){
            // remove qualifiers. XXX does this handle pointers correctly?
            thetype = underlying_types.at(thetype);
            std::cout<<"basic type "<<thetype<<std::endl;
            // for basic datatypes we can get the address of the underlying object
            // straight from the BinaryStream object's buffer member
            std::string* val_buffer = &((*store)[key]->buffer);
            // std::string is a "special" case in that, despite not being a fundamental type,
            // the BinaryStream buffer holds the address of a contiguous chunk of memory
            // that we can perform a direct copy from and give to TBranch::SetAddress
            // The only catch is we need to use &buffer[0] instead of &buffer.
            void* data_pointer = (thetype.find("std::string")!=std::string::npos)
                     ? (void*)val_buffer : (void*)&val_buffer->operator[](0);
            void* branch_pointer = (basic_type) ? data_pointer : &data_pointer;
            tree->SetBranchAddress(key.c_str(), branch_pointer);
            // the following is purely debug code - let's get the address back from the TBranch,
            // and try to cast it to the appropriate datatype, in order to check that
            // the data at that address is indeed interpetable as the right datatype and value.
            std::cout<<"data_pointer is "<<data_pointer<<std::endl;
            void** addp = (void**)tree->GetBranch(key.c_str())->GetAddress();
            std::cout<<"branch set up to fill from object at "<<(*addp)<<std::endl;
            std::cout<<"key "<<key<<" holds value ";
            // of course, to cast it back to the right type at runtime requires
            // an if statement for every type. But this is only for debug!
            if(key=="myint"){ int* mi = (int*)(data_pointer); std::cout<<(*mi)<<std::endl; }
            if(key=="mydouble"){ double* mi = (double*)(data_pointer); std::cout<<(*mi)<<std::endl; }
            // funnily enough although TBranch::Fill works with the address of the underlying
            // char array used by a std::string, we cannot cast this back to std::string,
            // or even char*...? So we can't check it like this!
            //if(key=="mystring"){ std::string* mi = (std::string*)(data_pointer); std::cout<<(*mi)<<std::endl; }
            // It actually isn't just fundamental types for which we can do this.
            // The determining factor is whether the memory layout of the object is contiguous.
            // Classes such as TVector3 and TString, which only contain primitive members
            // (no pointers/classes, can be copied with a shallow copy)
            // satisfy this condition, while containers such as std::vector will store a pointer
            // to memory space allocated on the heap, so cannot be read simply by TBranch::Fill
            // simply by setting a suitable address.
            // (This is why the BinaryStream class has to provide suitable container overloads
            //  of the streamer operator)
            if(key=="mystring"){ TString* mi = (TString*)(data_pointer); std::cout<<(*mi)<<std::endl; }
            if(key=="myvec"){ TVector3* mi = (TVector3*)(data_pointer);
            std::cout<<"{"<<mi->X()<<","<<mi->Y()<<","<<mi->Z()<<"}"<<std::endl; }
            // debug print, one last check our branch address is correct
            TBranch* p = tree->GetBranch(key.c_str());
            std::cout<<"branch at "<<p<<std::endl;
            void** add = (void**)p->GetAddress();
            std::cout<<" will fill with object at "<<(*add)<<std::endl;
            // Alright! Copy that data.
            p->Fill();
            // make a note of the number of entries the tree now has
            int n_entries = tree->GetBranch(key.c_str())->GetEntries();
            if(num_entries<n_entries) num_entries = n_entries;  // keep track of the max in any branch
            std::cout<<"done"<<std::endl;
        } else {
            // remove qualifiers, but leave internal pointers..?
            // XXX does this handle pointers correctly?
            thetype = RemoveQualifiers(thetype);
            std::cout<<"compound type "<<thetype<<std::endl;
            // else this is not a fundamental data type.
            // if this data type stores its data e.g. via a pointer to heap memory
            // (as in the case of e.g. stl containers), then realistically we need
            // to pass an appropriate object to BStore::Get() in order to have it
            // call the appropriate functions to read the data out of the BinaryStream.
            // This means instantiating an object of the appropriate type -
            // a type we only have as a runtime string! This is tricky.
            // We can use ROOT's TInterpreter to instantiate objects at runtime. So we
            // pass the TInterpreter a pointer to our BStore, invoke
            // mybstore->Get("name", temp_object) to extract the object from the store
            // into a temporary object. We then call branch->SetAddress(&temp_object)
            // followed by branch->Fill() to write it to the TTree. Simple!
            // ....and in ROOT 6, this works! (ROOT 6 uses Cling, not CINT).
            // In ROOT 5 there is a wrinkle - to invoke BStore::Get<MyType>, CINT must know
            // about the Get method for *that specific type* - and CINT dictionaries
            // do not support generic templated functions! Instead a generated dictionary
            // will only contain information to tell CINT about templated functions for
            // EXPLICITLY REQUESTED types - i.e. we need to have a line in the dictionary
            // `#pragma link C++ function BStore::Get<MyType>` to support every possible type
            // we want to retrieve. This could mean a LOT of #pragma lines!
            // Fortunately there is a 'partial' out - the CINT dictionary can be built
            // and loaded completely independently of the main application - 
            // we can even add new `#pragma` lines, rebuild the dictionary,
            // and reload it into CINT while the application is running!
            
            // we can skip all this song and dance if we've seen this type before.
            if(not is_known){
                // check if we have a a suitable `#pragma` line. We expect one of the following form:
                TString cmd = "#pragma link C++ function BStore::Get(std::string, " + thetype + "&);";
                cmd = "grep -sqF '" + cmd + "' BStore_Linkdef.hh";
                std::cout<<"grepping linkdef for line "<<cmd.Data()<<std::endl;
                get_ok = system(cmd.Data());
                if(get_ok!=0){
                    std::cout<<"Not found, checking for linkdef"<<std::endl;
                    // failed to find it. Double check the file exists.
                    get_ok = system("ls -1 BStore_Linkdef.hh > /dev/null");
                    if(get_ok!=0){
                        Log(m_unique_name+" Error! Failed to find BStore_Linkdef.hh!",v_error,m_verbose);
                        return false;
                    } else {
                        std::cout<<"linkdef exists, adding new line(s)"<<std::endl;
                        // Linkdef file exists, but we do not have an entry for this Type. Add one.
                        // make a backup before modification
                        system("cp -f BStore_Linkdef.hh  BStore_Linkdef_1.hh");
                        std::cout<<"backup made"<<std::endl;
                        // insert all necessary #include and #pragma extra_include lines
                        RegisterTypes(thetype);
                        
                        // OK, now link the specific templated function instance
                        cmd  = "echo '#pragma link C++ function BStore::Set(std::string, " 
                             + thetype + ");' >> BStore_Linkdef.hh"; // set by reference
                        std::cout<<"adding new line: "<<cmd<<std::endl;
                        get_ok |= system(cmd.Data());
                        std::cout<<"returned "<<get_ok<<std::endl;
                        cmd  = "echo '#pragma link C++ function BStore::Set(std::string, " 
                             + thetype + "*, bool);' >> BStore_Linkdef.hh"; // set by pointer
                        std::cout<<"adding new line: "<<cmd<<std::endl;
                        get_ok |= system(cmd.Data());
                        std::cout<<"returned "<<get_ok<<std::endl;
                        cmd = "echo '#pragma link C++ function BStore::Get(std::string, " 
                             + thetype + "&);' >> BStore_Linkdef.hh";  // get by reference
                        std::cout<<"adding line "<<cmd.Data()<<std::endl;
                        get_ok |= system(cmd.Data());
                        std::cout<<"returned "<<get_ok<<std::endl;
                        cmd = "echo '#pragma link C++ function BStore::Get(std::string, " 
                             + thetype + "*&);' >> BStore_Linkdef.hh";  // get by pointer
                        std::cout<<"adding line "<<cmd.Data()<<std::endl;
                        get_ok |= system(cmd.Data());
                        std::cout<<"returned "<<get_ok<<std::endl;
                        // rebuild the dictionary
                        std::cout<<"Please check Linkdef file"<<std::endl;
                        
                        std::cout<<"building lib"<<std::endl;
                        get_ok |= system("make lib/libBStore_RootDict.so");
                        std::cout<<"after building lib"<<std::endl;
                        if(get_ok!=0){
                            Log(m_unique_name+ " failed to update BStore dictionary to support new type "
                                + thetype + "! Please rebuild manually",v_error,m_verbose);
                            return false;
                        }
                        // if we had a TInterpreter, we need to delete it and make a new one
                        // otherwise the new method will not register.
                        if(meInterpreter) delete meInterpreter;
                        std::cout<<"making new interpreter"<<std::endl;
                        meInterpreter = new TCint("meInterpreter","title");
                        std::cout<<"loading library"<<std::endl;
                        meInterpreter->Load("lib/libBStore_RootDict.so");
                    }
                } // else it already existed in the Linkdef file. In either case, it should do now.
            } // else we've already encountered this type, so no need to modify the dictionary
            
            // ok so we should now have a dictionary supporting our type.
            // make an TInterpreter instance and load the dictionary
            std::cout<<"making interpreter"<<std::endl;
            if(!meInterpreter){
                meInterpreter = new TCint("meInterpreter","title");
                std::cout<<"loading library"<<std::endl;
                meInterpreter->Load("lib/libBStore_RootDict.so");
            }
            
            // construct a pointer to our compiled-space BStore in Interpreted space
            std::cout<<"constructing pointer to compile-time BStore"<<std::endl;
            TString cmd = TString::Format("BStore* bs = (BStore*)0x%x; ", store);
            meInterpreter->ProcessLine(cmd.Data(),&error);
            if(error!=0){
                Log(m_unique_name+" Error! Failed to instantiate BStore pointer in FillBranches!"
                    + " TInterpreter error was "+TInterpreterErrors.at(error),v_error,m_verbose);
                return false;
            }
            // instantiate a pointer to an instance of the desired object type.
            // (we could also instantiate an instance, but a pointer may be more efficient)
            // This should only fail if ROOT is unaware of the class.
            cmd = thetype +"* " + key + ";";
            std::cout<<"declaring temporary: "<<cmd.Data()<<std::endl;
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                Log(m_unique_name+" Error! Failed to instantiate pointer to type " + thetype
                    +" in TInterpreter! Have you made a ROOT dictionary for your class"
                    +" and placed a rootmap file on $LD_LIBRARY_PATH?",v_error,m_verbose);
                return false;
            }
            std::cout<<"Invoking BStore::Get"<<std::endl;
            cmd = "bs->Get(\""+key+"\", "+key+");";
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                Log(m_unique_name+" Error! Failed to invoke BStore::Get for type "+thetype,v_error,m_verbose);
            }
            std::cout<<"Get returned"<<std::endl;
            cmd = "printf(\"temporary now points to %p\\n\","+key+");";
            meInterpreter->ProcessLine(cmd.Data(), &error);
            // OK we should now have pulled the object from the BStore.
            // Well, that was a lot of work! Pass the TInterpreter our TTree
            cmd = TString::Format("TTree* t = (TTree*)0x%x;", tree);
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                Log(m_unique_name+" FillTree failed to set a tree pointer?",v_error,m_verbose);
                return false;
            }
            // set the branch address and fill the TBranch
            cmd = "TBranch* b = t->GetBranch(\""+key+"\");";
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                Log(m_unique_name+" FillTree failed to get branch "+key
                    +", did you call MakeBranches first?",v_error,m_verbose);
                return false;
            }  // FIXME that's a point, make sure new Store objects are not introduced??
            meInterpreter->ProcessLine("printf(\"b is at %p\\n\",b);");
            cmd = "b->SetAddress(&"+key+");";
            std::cout<<"executing: '"<<cmd.Data()<<"'"<<std::endl;
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                Log(m_unique_name+" Error! Failed to set address of branch "+key
                    + " to pointer to type " + thetype,v_error,m_verbose);
                return false;
            }
            tree->GetBranch(key.c_str())->Fill();
            int n_entries = tree->GetBranch(key.c_str())->GetEntries();
            if(num_entries<n_entries) num_entries = n_entries;  // keep track of the max in any branch
            // and that should be it!
        } // not a fundamental data type
    }
    tree->SetEntries(num_entries);
    tree->ResetBranchAddresses();
    return true;
}
