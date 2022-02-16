#include "StoreToTTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TClassEdit.h"

StoreToTTree::StoreToTTree(){
    if(verbosity>0) std::cout<<"Making interpreter"<<std::endl;
    meInterpreter = new TCint("meInterpreter","title");
    if(verbosity>0) std::cout<<"loading library"<<std::endl;
    meInterpreter->Load("lib/libBStore_RootDict.so");
}

StoreToTTree::~StoreToTTree(){
    // This makes the application segfault on exit. While harmless, it looks bad.
    // The culprit line is 'TClass::SetUnloaded' called by ROOT's libCore or something.
    // I think the easiest solution is just to not manually delete the interpreter. ¯\_(ツ)_/¯
    //if(meInterpreter) delete meInterpreter; meInterpreter=nullptr;
}

// functions to transfer Stores or BStores into a TTree
void StoreToTTree::MakeBranches(TTree* tree, Store* store) {
    for (auto const& pair: *(store->GetMap())) {
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

void StoreToTTree::MakeBranches(TTree* tree, BStore* store) {
    if(verbosity>1) std::cout<<"Scanning entries in BStore"<<std::endl;
    for (auto const& pair: store->m_variables) {
        auto key = pair.first.c_str();
        if(verbosity>1) std::cout<<"entry "<<key<<std::endl;
        if(tree->FindBranch(key)) continue;  // continue if it already exists
        
        if(verbosity>1) std::cout<<"no existing branch, getting type"<<std::endl;
        std::string thetype = store->Type(key);
        if(verbosity>1) std::cout<<"raw type '"<<thetype<<"'"<<std::endl;
        thetype = abi::__cxa_demangle(thetype.c_str(), nullptr, nullptr, nullptr);
        if(verbosity>1) std::cout<<"demangled type '"<<thetype<<"'"<<std::endl;
        // Resolve typedefs: e.g. if someone has done
        // typedef std::array<double,3> pos;
        // then TClassEdit::ResolveTypedef("pos") would return "array<double,3>"
        thetype = TClassEdit::ResolveTypedef(thetype.c_str());
        if(verbosity>1) std::cout<<"typedef resolved type '"<<thetype<<"'"<<std::endl;
        // check if it's an enum and decay it to an integer if so -
        // otherwise if this enum isn't known about in the gInterpreter, it'll fall over.
        // XXX what if the gInterpreter does know about the enum type and we wish to preserve it?
        bool is_enum = gInterpreter->ClassInfo_IsEnum(thetype.c_str());
        if(is_enum){
            if(verbosity>1) std::cout<<" decaying enum to int"<<std::endl;
            thetype = "int";  // TODO this functionality is UNTESTED!
        }
        
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
        if(verbosity>1) std::cout<<"removing 'std::' namespace prefixes"<<std::endl;
        while(thetype.find(subs)!=std::string::npos){ thetype.erase(thetype.find(subs),len); }
        // also remove qualifiers like 'const' and pointers...but only for the outermost type, if nested.
        // XXX is this correct functionality for pointers?
        if(verbosity>1) std::cout<<"removing qualifiers"<<std::endl;
        thetype = RemoveQualifiers(thetype);
        // if a user tries to store a container of pointers such as a vector<Hit*>
        // we will error. I suppose what the user may *like* would be to have us
        // dereference and store the elements, but that's going a step too far.
        if(verbosity>1) std::cout<<"checking if the type is a container of pointers"<<std::endl;
        if(thetype.find('*')!=std::string::npos){
            std::cerr<<toolName<<" Error! Can't make a TBranch for BStore entry "<<key
                     <<" of type "<<thetype<<" which contains pointers!"<<std::endl;
            continue;
        }
        // in order to not split objects, we have to know if it's a basic type
        // and then modify the branching command accordingly
        bool basic_type = (fundamental_types.count(thetype)>0);
        if(verbosity>1) std::cout<<"basic type = "<<basic_type<<std::endl;
        
        // BStores can store many different types.
        // Rather than having a switch statement based on the type string,
        // we can hack it using a TInterpreter, provided ROOT is aware of the type.
        // TODO we could probably load libraries into the interpreter to support
        // e.g. SKROOT class objects as well, if need be.
        TString command = TString::Format("TTree* t = (TTree*)%p; ", (void*)tree);
        if(basic_type){
            command += thetype + " " + key + "; TBranch* b = t->Branch(\"" + key + "\", &" + key + ");";
        } else {
            // splitlevel->0, do not fragment objects.
            command += thetype + "* "+ key + "=0; TBranch* b = t->Branch(\"" + key + "\", &" + key + ",32000,0);";
        }
        if(verbosity>0) std::cout<<"Making Tree branch with command: "<<std::endl<<command.Data()<<std::endl;
        if(!meInterpreter){
            if(verbosity>0) std::cout<<"Making interpreter"<<std::endl;
            meInterpreter = new TCint("meInterpreter","title");
            if(verbosity>0) std::cout<<"loading library"<<std::endl;
            meInterpreter->Load("lib/libBStore_RootDict.so");
        }
        meInterpreter->ProcessLine(command.Data());
        if(verbosity>1) std::cout<<"done"<<std::endl;
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

void StoreToTTree::MakeBranches(TTree* tree, std::map<std::string, Store*> &stores){
    for(auto&& astore : stores) MakeBranches(tree, astore.second);
}

bool StoreToTTree::FillBranches(TTree* tree, std::map<std::string, Store*> &stores){
    bool ret=true;
    for(auto&& astore : stores) ret &= FillBranches(tree, astore.second);
    return ret;
}

// FIXME we could probably replace this with TClassEdit::ShortType(const char* fulltypename, int mode)
std::string StoreToTTree::RemoveQualifiers(std::string thetype, bool allpointers){
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

bool StoreToTTree::RegisterTypes(std::string thetype){
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

bool StoreToTTree::RegisterType(std::string thetype){
    // we need to add two new lines to #include the header defining this class
    // (with the exception of (most?) stl containers).
    // For classes of which ROOT is aware we can determine this from the TClass:
    TClass* tc = TClass::GetClass(thetype.c_str());
    if(tc==nullptr){
        std::cerr<<toolName<<" Error! Could not build TClass from "<<thetype
                 <<"! Have you built a ROOT dictionary for it? (See Makefile)"<<std::endl;
        return false;
    }
    const char* declared_in = tc->GetDeclFileName();
    if(declared_in==nullptr || strcmp(declared_in,"")==0){
        // this could be e.g. because either it's unknown to ROOT
        // or because it's a combination of classes,
        // e.g. std::vector<TString>
        std::cerr<<toolName<<" Error! Could not determine header "<<thetype
                 <<" resides in! Have you built a ROOT dictionary for it? (See Makefile)"<<std::endl;
        return false;
    }
    std::string decl_file = declared_in;
    if(decl_file.find("prec_stl")!=std::string::npos){
        if(verbosity>0) std::cout<<"skipping includes for class "<<thetype
                 <<" as we think this is a stl container of fundamentals"<<std::endl;
    } else {
        // ok we have a header
        // strip any preceding path from the header file
        int startpos = decl_file.find_last_of('/');
        startpos = (startpos==std::string::npos) ? 0 : startpos+1;
        std::string headername = decl_file.substr(startpos, std::string::npos);
        if(verbosity>0) std::cout<<thetype<<" is declared in "<<headername<<std::endl;
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

bool StoreToTTree::FillBranches(TTree* tree, Store* store){
    // first make sure the branches are present
    if(allowRaggedBranches){
        if(verbosity>0) std::cout<<"StoreToTTree making branches"<<std::endl;
        MakeBranches(tree, store);
    }
    // copy data across
    if(verbosity>0) std::cout<<"StoreToTTree copying data"<<std::endl;
    int num_entries=-1;
    for (auto const& pair: *(store->GetMap())) {
        std::string key = pair.first;
        std::string value = pair.second;
        float tmpNum;
        std::string tmpStr;
        std::string* tmpStrPtr = &tmpStr;
        if(store->Get(key, tmpNum)) {
            if(verbosity>0) std::cout<<"Store saving number "<<tmpNum<<", setting address to "<<&tmpNum<<std::endl;
            tree->SetBranchAddress(key.c_str(), &tmpNum);
            /*
            TBranch* p = tree->GetBranch(key.c_str());
            std::cout<<"branch at "<<p<<std::endl;
            void** add = (void**)p->GetAddress();
            std::cout<<" will fill with object at "<<(add)<<std::endl;
            */
        } else {
            tmpStr = value;
            if(verbosity>0) std::cout<<"Store saving string "<<tmpStr
                                     <<", setting address to "<<&tmpStrPtr<<std::endl;
            tree->SetBranchAddress(key.c_str(), &tmpStrPtr);
            /*
            TBranch* p = tree->GetBranch(key.c_str());
            std::cout<<"branch at "<<p<<std::endl;
            void** add = (void**)p->GetAddress();
            std::cout<<" will fill with object at "<<(*add)<<std::endl;
            */
        }
        if(verbosity>0) std::cout<<"calling fill on branch "<<key<<std::endl;
        tree->GetBranch(key.c_str())->Fill();
        if(verbosity>0) std::cout<<"variable is set"<<std::endl;
        int n_entries = tree->GetBranch(key.c_str())->GetEntries();
        if(num_entries<n_entries) num_entries = n_entries;  // keep track of the max in any branch
        if(verbosity>0) std::cout<<"done"<<std::endl;
    }
    tree->SetEntries(num_entries);
    tree->ResetBranchAddresses();
    return true;
}

bool StoreToTTree::FillBranches(TTree* tree, BStore* store){
    // first make sure the branches are present
    if(allowRaggedBranches){
        if(verbosity>0) std::cout<<"StoreToTTree making branches"<<std::endl;
        MakeBranches(tree, store);
    }
    int num_entries=-1;
    if(verbosity>0) std::cout<<"StoreToTTree copying data"<<std::endl;
    for (auto const& pair: store->m_variables) {
        std::string key = pair.first;
        std::string thetype = store->Type(key);
        if(verbosity>0) std::cout<<"raw type "<<thetype<<std::endl;
        thetype = abi::__cxa_demangle(thetype.c_str(), nullptr, nullptr, nullptr);
        if(verbosity>0) std::cout<<"demangled type "<<thetype<<std::endl;
        thetype = TClassEdit::ResolveTypedef(thetype.c_str());
        if(verbosity>0) std::cout<<"typedef resolved type "<<thetype<<std::endl;
        
        // have we encountered this type before: if so we can skip
        // the checks about whether the BStore is aware of its type
        bool is_known = (known_types.count(thetype)>0);
        bool basic_type;
        if(is_known){
            basic_type = known_types.at(thetype);
        } else {
            // we've not seen this type before.
            // to be able to check if this is a fundamental type we need to strip qualifiers
            std::string underlying_type = RemoveQualifiers(thetype);
            // check if it's a fundamental type
            basic_type = (fundamental_types.count(underlying_type)!=0);
            // record our results for the future
            known_types.emplace(thetype, basic_type);
            underlying_types.emplace(thetype, underlying_type);
        }
        if(basic_type){
            // remove qualifiers. XXX does this handle pointers correctly?
            thetype = underlying_types.at(thetype);
            if(verbosity>0) std::cout<<"basic type "<<thetype<<std::endl;
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
            /*
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
            */
            // Alright! Copy that data.
            tree->GetBranch(key.c_str())->Fill();
            // make a note of the number of entries the tree now has
            int n_entries = tree->GetBranch(key.c_str())->GetEntries();
            if(num_entries<n_entries) num_entries = n_entries;  // keep track of the max in any branch
            if(verbosity>0) std::cout<<"done"<<std::endl;
        } else {
            // remove qualifiers, but leave internal pointers..?
            // XXX does this handle pointers correctly?
            thetype = RemoveQualifiers(thetype);
            if(verbosity>0) std::cout<<"compound type "<<thetype<<std::endl;
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
                if(verbosity>0) std::cout<<"grepping linkdef for line "<<cmd.Data()<<std::endl;
                get_ok = system(cmd.Data());
                if(get_ok!=0){
                    if(verbosity>0) std::cout<<"Not found, checking for linkdef"<<std::endl;
                    // failed to find it. Double check the file exists.
                    get_ok = system("ls -1 BStore_Linkdef.hh > /dev/null");
                    if(get_ok!=0){
                        std::cerr<<toolName<<" Error! Failed to find BStore_Linkdef.hh!"<<std::endl;
                        return false;
                    } else {
                        if(verbosity>0) std::cout<<"linkdef exists, adding new line(s)"<<std::endl;
                        // Linkdef file exists, but we do not have an entry for this Type. Add one.
                        // make a backup before modification
                        system("cp -f BStore_Linkdef.hh  BStore_Linkdef_1.hh");
                        if(verbosity>0) std::cout<<"backup made"<<std::endl;
                        // insert all necessary #include and #pragma extra_include lines
                        RegisterTypes(thetype);
                        
                        // OK, now link the specific templated function instance
                        cmd  = "echo '#pragma link C++ function BStore::Set(std::string, " 
                             + thetype + ");' >> BStore_Linkdef.hh"; // set by reference
                        if(verbosity>0) std::cout<<"adding new line: "<<cmd<<std::endl;
                        get_ok |= system(cmd.Data());
                        if(verbosity>0) std::cout<<"returned "<<get_ok<<std::endl;
                        cmd  = "echo '#pragma link C++ function BStore::Set(std::string, " 
                             + thetype + "*, bool);' >> BStore_Linkdef.hh"; // set by pointer
                        if(verbosity>0) std::cout<<"adding new line: "<<cmd<<std::endl;
                        get_ok |= system(cmd.Data());
                        if(verbosity>0) std::cout<<"returned "<<get_ok<<std::endl;
                        cmd = "echo '#pragma link C++ function BStore::Get(std::string, " 
                             + thetype + "&);' >> BStore_Linkdef.hh";  // get by reference
                        if(verbosity>0) std::cout<<"adding line "<<cmd.Data()<<std::endl;
                        get_ok |= system(cmd.Data());
                        if(verbosity>0) std::cout<<"returned "<<get_ok<<std::endl;
                        cmd = "echo '#pragma link C++ function BStore::Get(std::string, " 
                             + thetype + "*&);' >> BStore_Linkdef.hh";  // get by pointer
                        if(verbosity>0) std::cout<<"adding line "<<cmd.Data()<<std::endl;
                        get_ok |= system(cmd.Data());
                        if(verbosity>0) std::cout<<"returned "<<get_ok<<std::endl;
                        // rebuild the dictionary
                        if(verbosity>0) std::cout<<"Please check Linkdef file"<<std::endl;
                        
                        if(verbosity>0) std::cout<<"building lib"<<std::endl;
                        get_ok |= system("make lib/libBStore_RootDict.so");
                        if(verbosity>0) std::cout<<"after building lib"<<std::endl;
                        if(get_ok!=0){
                            std::cerr<<toolName<< " failed to update BStore dictionary to support new type "
                                     <<thetype<<"! Please rebuild manually"<<std::endl;
                            return false;
                        }
                        // if we had a TInterpreter, we need to delete it and make a new one
                        // otherwise the new method will not register.
                        if(meInterpreter) delete meInterpreter;
                        if(verbosity>0) std::cout<<"making new interpreter"<<std::endl;
                        meInterpreter = new TCint("meInterpreter","title");
                        //meInterpreter = (TCint*)gInterpreter;
                        if(verbosity>0) std::cout<<"loading library"<<std::endl;
                        meInterpreter->Load("lib/libBStore_RootDict.so");
                    }
                } // else it already existed in the Linkdef file. In either case, it should do now.
            } // else we've already encountered this type, so no need to modify the dictionary
            
            // ok so we should now have a dictionary supporting our type.
            // make an TInterpreter instance and load the dictionary
            if(verbosity>0) std::cout<<"making interpreter"<<std::endl;
            if(!meInterpreter){
                meInterpreter = new TCint("meInterpreter","title");
                //meInterpreter = (TCint*)gInterpreter;
                if(verbosity>0) std::cout<<"loading library"<<std::endl;
                meInterpreter->Load("lib/libBStore_RootDict.so");
            }
            
            // construct a pointer to our compiled-space BStore in Interpreted space
            if(verbosity>0) std::cout<<"constructing pointer to compile-time BStore"<<std::endl;
            TString cmd = TString::Format("BStore* bs = (BStore*)%p; ", (void*)store);
            meInterpreter->ProcessLine(cmd.Data(),&error);
            if(error!=0){
                std::cerr<<toolName<<" Error! Failed to instantiate BStore pointer in FillBranches!"<<std::endl;
                //       << " TInterpreter error was "<<TInterpreterErrors.at(error)<<std::endl;
                // TInterpreterErrors are not useful, they provide no info on the error that occurred.
                // TODO But really we should check if it's fatal and spin up a new TCint if required.
                return false;
            }
            // instantiate a pointer to an instance of the desired object type.
            // (we could also instantiate an instance, but a pointer may be more efficient)
            // This should only fail if ROOT is unaware of the class.
            cmd = thetype +"* " + key + ";";
            if(verbosity>0) std::cout<<"declaring temporary: "<<cmd.Data()<<std::endl;
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                std::cerr<<toolName<<" Error! Failed to instantiate pointer to type "<<thetype
                         <<" in TInterpreter! Have you provided a ROOT dictionary for your class?"
                         <<" (see Makefile)"<<std::endl;
                return false;
            }
            if(verbosity>0) std::cout<<"Invoking BStore::Get"<<std::endl;
            cmd = "bs->Get(\""+key+"\", "+key+");";
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                std::cerr<<toolName<<" Error! Failed to invoke BStore::Get for type "<<thetype<<std::endl;
            }
            if(verbosity>0) std::cout<<"Get returned"<<std::endl;
            cmd = "printf(\"temporary now points to %p\\n\","+key+");";
            if(verbosity>0) meInterpreter->ProcessLine(cmd.Data(), &error);
            // OK we should now have pulled the object from the BStore.
            // Well, that was a lot of work! Pass the TInterpreter our TTree
            cmd = TString::Format("TTree* t = (TTree*)%p;", (void*)tree);
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                std::cerr<<toolName<<" FillBranches failed to set a tree pointer?"<<std::endl;
                return false;
            }
            // set the branch address and fill the TBranch
            cmd = "TBranch* b = t->GetBranch(\""+key+"\");";
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                std::cerr<<toolName<<" FillBranches failed to get branch "<<key
                         <<", did you call MakeBranches first?"<<std::endl;
                return false;
            }
            if(verbosity>0) meInterpreter->ProcessLine("printf(\"b is at %p\\n\",b);");
            cmd = "b->SetAddress(&"+key+");";
            if(verbosity>0) std::cout<<"executing: '"<<cmd.Data()<<"'"<<std::endl;
            meInterpreter->ProcessLine(cmd.Data(), &error);
            if(error!=0){
                std::cerr<<toolName<<" Error! Failed to set address of branch "<<key
                         << " to pointer to type "<<thetype<<std::endl;
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
