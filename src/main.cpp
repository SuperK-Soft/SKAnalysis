#include <string>
#include "ToolChain.h"
//#include "DummyTool.h"

int main(int argc, char* argv[]){

  std::string config_file;
  if (argc==1)config_file="configfiles/Dummy/ToolChainConfig";
  else config_file=argv[1];

  ToolChain tools(config_file, argc, argv);


  //DummyTool dummytool;    

  //tools.Add("DummyTool",&dummytool,"configfiles/DummyToolConfig");

  //int portnum=24000;
  //  tools.Remote(portnum);
  //tools.Interactive();
  
  //  tools.Initialise();
  // tools.Execute();
  //tools.Finalise();
  
  
  return 0;
  
}
