#include "TreeManagerMod.h"
#include <string>

// this is literally just a copy-paste of TreeManager::Initialise,
// but takes a tree name instead of hard-coding it to 'data'
void TreeManagerMod::Initialize(std::string treename)
{
  theTree  = 0;
  NEntries = 0; 

  if ( mode ==1 )
    {
      std::cerr << " You are calling Initialize but mode=1 (ZBS2ROOT) !\n"
		<< " You should not be here.... check your code or ask an expert\n";
      exit(-1);
    }

  std::vector<std::string>::const_iterator IT;
  IT = ListOfFiles.begin();
  file = new TFile((*IT).c_str());
  TTree* tree = (TTree*)file->Get(treename.c_str());

  // Check for TChain or TTree usage
  if( ListOfFiles.size() > 1 )
  {
      // make the chain
     TChain* chain = new TChain(treename.c_str(),"");
     for ( IT = ListOfFiles.begin() ; IT != ListOfFiles.end() ; ++IT)
       NEntries += chain->Add((*IT).c_str());

     theTree  = chain;
   }
   else 
   { 
     // make the tree, if there is only one file
     theTree  = tree; 
   }    
   
  // initial branch association done only once, before tree cloning 
  // if theTree is a TTree and not a TChain -rvw
  if( !theTree->InheritsFrom(TChain::Class()) ) SetupBranches(theTree); 

  theTree->SetCacheSize(40*1024*1024);
  pos = -1;  // need to call Next() at the first
  current_tree=-1;

  // Write mode : an output file is created "root2root mode"
  if ( mode == 0 ) {
      outputFile = new TFile(outputName.c_str(),"RECREATE"); 

      // Turn off branches we don't need to read in
      SetBranches(theTree,InBranchesToSkip,0);

      // Branches we want to read in but don't want to write out:
      // turn them off in the input tree before cloning
      SetBranches(theTree,OutBranchesToSkip,0);

      // Create the Clone, 0 specifies that events are not to be copied now.
      theOTree = (TTree*) theTree->CloneTree(0);  

      // Now turn them on again in the input tree so they will be read in
      // This will create warnings from the cloneTree, which does not have
      // these branches, and can be ignored. -rvw
      SetBranches(theTree ,OutBranchesToSkip,1);

      theOTree->SetCacheSize(40*1042*1024);

      //------------------------------------------------------------------------
      TFileCacheWrite *cachew = new TFileCacheWrite(outputFile, 40*1024*1024);
      //------------------------------------------------------------------------
      theTree->GetEntry(); // added to avoid an error in fill_tree() 20080717 is this ok?? (y.t.)
      printf("Initialized: total = %lld\n", theTree->GetEntries() );

  }
  // no output written : "read root mode"
  else if ( mode == 2 ) {
      outputFile = 0;
      theOTree = 0;
      pos = -1;   // need to call Next() immediately
      printf("InitReadRoot: total = %lld\n", theTree->GetEntries());
  }

  // initial branch association done only once, after any tree cloning 
  // if theTree is a TChain and not a TTree -rvw
  if( theTree->InheritsFrom(TChain::Class()) ) SetupBranches(theTree); 

  // For multiple files, this check should be done after at least one 
  // file has been connected to the TChain  -rvw
  NEntries = theTree->GetEntries();
  // replaced <= with < to allow 0 event in a file. 20090526 y.t.
  // if (NEntries <= 0) {
  if (NEntries < 0) {
      std::cerr << "TreeManager::Initialize() Error! NEntries = " << NEntries << ".\n";
      exit(-2);
  }

  // moved tree search above
  TList* RareList = (TList*) tree->GetUserInfo()->FindObject("Rare");
  int numrare=0;
  if (RareList) {
      numrare = RareList->GetSize();
  }

  TIter next(RareList);
  TObject* rare;

  if (!RareList || numrare==0) {
      std::cout << " No Rare List!" << std::endl;
  }
  else {
      std::string rarename="RUNINFO";
      
      if ((rare = (TObject *)next())) {
	  std::string objname = rare->GetName();
	  if (rarename.compare(objname) == 0) {
	      std::cout << "RUNINFO found in the rare list" <<std::endl;
	      RINFO = (RunInfo*)rare;
	  }
      }
  }

  // clear (may not be needed... (y.t.)) 
  Clear();
}

