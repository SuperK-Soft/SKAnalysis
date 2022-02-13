#!/bin/bash

init=1
tooldaq=1
boostflag=1
zmq=1
final=1
rootflag=0
setup=1
threads=`nproc --all`

while [ ! $# -eq 0 ]
do
    case "$1" in
	--help | -h)
	    echo "This script should be run once after initially cloning the ToolApplication repository. It retrieves the ToolFrameworkCore TooDAQFramework ZMQ and BOOST repositories that provides the core framework and dependancies on which your application will be built."
	    exit
	    ;;

	--with_root | -r)
	    echo "Installing ToolDAQ with root"
	    rootflag=1 
	    ;;
	
	--no_boost | -b)
            echo "Installing ToolDAQ without boost"
            boostflag=0
	    ;;
	
	--no_init )
	     echo "Installing ToolDAQ without creating ToolDAQ Folder"
	    init=0;
	    ;;

	--no_zmq )
            echo "Installing ToolDAQ without zmq"
            zmq=0
            ;;

	--no_tooldaq )
	    echo "Installing dependancies without ToolDAQ"
	    tooldaq=0
	    ;;

	--no_final )
            echo "Installing ToolDAQ without compiling ToolAnalysis"
            final=0
            ;;

	--ToolDAQ_ZMQ )
            echo "Installing ToolDAQ & ZMQ"
	    boostflag=0
	    rootflag=0
	    final=0
            ;;

	--Boost )
            echo "Installing Boost"
	    init=0
	    tooldaq=0
	    zmq=0
	    final=0
	    rootflag=0
            ;;

	--Root )
            echo "Installing Root"
	    init=0
	    tooldaq=0
	    boostflag=0
	    zmq=0
	    final=0
	    rootflag=1
            ;;
	
	
	--Final )
            echo "Compiling ToolDAQ"
	    init=0
	    tooldaq=0
	    boostflag=0
	    rootflag=0
	    zmq=0
            ;;

    esac
    shift
done

if [ $init -eq 1 ]
then
    
    mkdir Dependencies
fi

cd Dependencies

if [ $tooldaq -eq 1 ]
then
git clone https://github.com/ToolFramework/ToolFrameworkCore.git

cd ToolFrameworkCore
make clean
make -j $threads

export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
cd ../

fi

if [ $zmq -eq 1 ]
then
    git clone https://github.com/ToolDAQ/zeromq-4.0.7.git
    
    cd zeromq-4.0.7
    
    ./configure --prefix=`pwd`
    make -j $threads
    make install
    
    export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
    
    cd ../
fi

if [ $boostflag -eq 1 ]
then
    
    git clone https://github.com/ToolDAQ/boost_1_66_0.git
     
    cd boost_1_66_0

    rm -rf INSTALL    
    mkdir install 
    
    ./bootstrap.sh --prefix=`pwd`/install/  > /dev/null 2>/dev/null
    ./b2 install iostreams -j $threads
    
    export LD_LIBRARY_PATH=`pwd`/install/lib:$LD_LIBRARY_PATH
    cd ../
fi


if [ $rootflag -eq 1 ]
then
    
    wget https://root.cern.ch/download/root_v5.34.34.source.tar.gz
    tar zxvf root_v5.34.34.source.tar.gz
    rm -rf root_v5.34.34.source.tar.gz
    cd root
    
    ./configure --enable-rpath
    make -j $threads
    make install
    
    source ./bin/thisroot.sh
    
    cd ../
    
fi

if [ $tooldaq -eq 1 ]
then
    git clone https://github.com/ToolDAQ/ToolDAQFramework.git
    
    cd ToolDAQFramework
    make clean
    make -j $threads
    export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
    cd ../
    
fi


cd ../

if [ $final -eq 1 ]
then
    
    echo "current directory"
    echo `pwd`
if [ $setup -eq 1 ]
then   
    cp -r ./Dependencies/ToolFrameworkCore/DataModel/* ./DataModel
    cp -r ./Dependencies/ToolDAQFramework/DataModel/* ./DataModel
    cp -r ./Dependencies/ToolFrameworkCore/UserTools/* ./UserTools
    cp -r ./Dependencies/ToolDAQFramework/UserTools/template/* ./UserTools/template
    cp -r ./Dependencies/ToolDAQFramework/configfiles/* ./configfiles
    mkdir src
    cp -r ./Dependencies/ToolDAQFramework/src/main.cpp ./src/
    cp ./Dependencies/ToolDAQFramework/Application/* ./
    git add DataModel/*
    git add UserTools/*
    git add configfiles/*
    git add ./Makefile
    git add ./CMakeLists.txt
    git add ./Setup.sh
    git add ./src/main.cpp
    rm -f ./GetToolFramework.sh
    sed -i 's/setup=1/setup=0/' ./GetToolDAQ.sh
fi   
    make clean
    make -j $threads
    
    export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
fi
