#!/bin/bash

tool=$1

if  command -v dialog &> /dev/null
then
    
    
    if [ "$tool" == "" ]
    then

	num=1
	
	
	dialog --radiolist "Choose Tool with <spacebar>, Enter for OK and ESC for Cancel:" 0 0 0 \
	    `for Tool in \`ls */*.cpp | grep -v template/ | sed s:/:' ': | awk '{print $2}' | sed s:.cpp:: | grep -v Factory\`
do
    fin=0
    for current in \`cat Factory/Factory.cpp | grep -v '/' |grep if| awk '{print $4}' | sed s:';':: \`
    do
        if [ $Tool == $current ]
        then
            fin=1
            echo "$Tool $num off "
        fi
    done

    if [ $fin -eq 0 ]
    then
        echo "$Tool $num off "
    fi
    num=$(expr 1 + $num)
done

#for Tool in \`ls InactiveTools/*/*.cpp | sed s:/:' ':g | awk '{print $3}' | sed s:.cpp:: \`
#do
#    echo "$Tool $num off "
#    num=$(expr 1 + $num)
#done` 2> tmptool
	

	if [ $? -eq 0 ]
	then
	    tool=`cat tmptool`
	    rm -f tmptool
	    
	    if [ "$tool" == "" ]
	    then
		rm -f tmptool 
		clear
		echo -e "\e[38;5;196mERROR!!! no Tool name given \e[0m"
		exit
	    fi
	    
	else
	    clear
	    exit
	fi
	
    fi
    
    if [ -d $tool ]
    then
	
	dialog --title "Message"  --yesno "Are you sure you want to perminenatly delete Tool $tool? \n \n Note: Tools can be temporarily disabled using ToolSelect.sh instead" 0 0
	
	if [ $? -eq 0 ]
	then
	    clear
	    echo -e "\e[38;5;226mDeleting Tool $tool  \e[0m"
	    rm -r $tool
	    cat Unity.h | grep -v $tool.h > Unity.new
	    mv Unity.new Unity.h
	    cat Factory/Factory.cpp | grep -v "$tool;" > Factory/Factory.new
	    mv Factory/Factory.new Factory/Factory.cpp
	else
	    clear
	    exit
	fi
    else
	clear 
	echo -e "\e[38;5;196mError Tool $tool doest exist \e[0m"	
	exit
    fi
else	
    
    
    if [ "$tool" == "" ]
    then
	echo -e "\e[38;5;196mError no name given\e[0m: usage = \" ./deleteTool.sh \e[38;5;226m <ToolNAME> \e[0m\" "	
    else
	if [ -d $tool ]
	then
	    echo -e "\e[38;5;196mCAUTION!!!! About to perminantly delete Tool $tool are you sure? (\e[38;5;226mYes\e[38;5;196m,\e[38;5;226m No\e[38;5;196m) \n\e[38;5;46mNote: Tools can be temporarily disabled using ToolSelect.sh instead \e[0m"
	    read answer
	    while [[ "$answer" != "Yes" && "$answer" != "No" ]]
	    do
		echo -e "\e[38;5;196mMust be \"Yes\" or \"No\"  \e[0m"
		read answer
	    done
	    if [ "$answer" == "Yes" ]
	    then
		echo -e "\e[38;5;226mDeleting Tool $tool  \e[0m"
		rm -r $tool
		cat Unity.h | grep -v $tool.h > Unity.new
		mv Unity.new Unity.h
		cat Factory/Factory.cpp | grep -v "$tool;" > Factory/Factory.new
		mv Factory/Factory.new Factory/Factory.cpp
	    else
		exit
	    fi
	else
	    echo -e "\e[38;5;196mError Tool $tool doest exist \e[0m"
	    exit
	fi
	
    fi

fi

    
