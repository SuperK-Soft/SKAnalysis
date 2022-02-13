#!/bin/bash

if [ "$1" == "" ]
then
    if  command -v dialog &> /dev/null
    then

	dialog --inputbox "Enter your ToolChain config name:" 0 0 2>answer
	if [ $? -eq 0 ]
        then
            dir=`cat answer`
	    echo dir $dir
	    rm -f answer
	    if [ "$dir" == "" ]
	    then
		clear
		rm -f answer
		echo -e "\e[38;5;196mError!!! no name given \e[0m"
		exit
	    else
		clear
                rm -f answer
	    fi
	else
	    clear
	    rm -f answer
	    exit
	fi

    else
	echo -e "\e[38;5;196mError no name given \e[0m: usage = \"./Create_run_config.sh \e[38;5;226m <ToolChainNAME> \e[0m "
    fi
else

dir=$1

fi

if [ -d $dir ]
then
    echo -e "\e[38;5;196mERROR!!! ToolChain config directory already exists \e[0m"
else
    
     cp -R template $dir
     more $dir/ToolChainConfig | sed s:"configfiles/":"configfiles/"$dir"/": > $dir/tmp
     mv $dir/tmp $dir/ToolChainConfig
     more $dir/ToolsConfig | sed s:"configfiles/":"configfiles/"$dir"/": > $dir/tmp
     mv $dir/tmp $dir/ToolsConfig
     ln -s  configfiles/$dir/ToolChainConfig ../$dir
     
fi
