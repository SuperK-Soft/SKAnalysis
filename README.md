![toolapplication logo](https://user-images.githubusercontent.com/14093889/147496518-f3751cd6-0c57-4dd1-8517-3a02b61e59f5.png)

# ***Please clone/fork this repository to build your own ToolApplication in either ToolFramework or ToolDAQ***

ToolApplication is a template repository for building your own applications in either: 

  - ToolFramework:- an open source general modular C++ Framework.
  - ToolDAQ:- DAQ focused frameowrk built on ToolFramework


****************************
# Installation
****************************

There are a few choices for installation, as mentioned this is an application template repository so the idea would be to fork or clone this repo into your own repo first, or build with it as a base. Some containers are provided for installation or trying out the code as well.

1. Using Docker / Singularity

   - For testing a complete install of Toolframework and ToolDAQ can be used by either:
     - ``` docker run --name=ToolFramework -it toolframework/toolframeworkapp```
     - ``` docker run --name=ToolDAQ -it tooldaq/tooldaqapp```
   - If you want a container to use as a base for your own application container you can use either:
     - ```toolframework/centos7``` which is a lightweight centos build with the prerequisites to install any ToolApplication from source
     - ```toolframework/core``` which is the same as above but with ToolFrameworkCore already installed in /opt/. This is a useful base for building ToolFramework containers
     - ```tooldaq/core``` which is the same as above but with ToolFrameworkCore, and ToolDAQframework plus its dependencies of boost and zmq already installed in /opt/. This is a useful base for building ToolDAQ containers


2. Install from source

   - Install Prerequisites: 
     - RHEL/Centos... ``` yum install git make gcc-c++ zlib-devel dialog ```
     - Debian/Ubuntu.. ``` apt-get install git make g++ libz-dev dialog ```

   - Then clone the repo with ```git clone https://github.com/ToolFramework/ToolApplication.git``` or more likely your own fork

   - Once clonned please run either:

     - ```./GetToolFramework.sh``` to install dependances and files for creating a ToolFramework app
     - ```./GetToolDAQ.sh``` to install dependances and files for creating a ToolDAQ app

(Note: if your doing this from inside one of the pre prepared core containers the core components will already be in the containers. In which case you should instead do:

 - ```ln -s /opt ./Dependencies```

followed by either:

 - ```./GetToolFramework.sh  --Final```
 - ```./GetToolDAQ.sh --Final```

To set up your application )



****************************
# Concept
****************************

The main executable creates a ToolChain which is an object that holds Tools. Tools are added to the ToolChain and then the ToolChain can be told to Initialise Execute and Finalise each Tool in the chain.

The ToolChain also holds a uesr defined DataModel which each Tool has access too and can read ,update and modify. This is the method by which data is passed between Tools.

User Tools can be generated for use in the ToolChain by included scripts.

For more information consult the ToolFramework manual:

https://docs.google.com/document/d/18rgYMOAGt3XiW9i0qN9kfUK1ssbQgLV1gQgG3hyVUoA

or the Doxygen docs

- https://toolframework.github.io/ToolApplication/

- https://toolframework.github.io/ToolFrameworkCore/

- http://tooldaq.github.io/ToolDAQFramework/

Copyright (c) 2016 Benjamin Richards (benjamin.richards@warwick.ac.uk)
