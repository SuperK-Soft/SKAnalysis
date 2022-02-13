#!/bin/bash

tool=$1
newname=$2

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

    if [ "$newname" == "" ]
    then
	dialog --inputbox "What would you like the new name to be?:" 0 0 2>answer
	
	if [ $? -eq 0 ]
	then
	    
	    newname=`cat answer`
	    rm -f answer

	else
	    rm -f answer
	    clear
	    exit
	fi
    fi
    


else


    if [ "$tool" == "" ]
    then
	echo -e "\e[38;5;196mError no name given\e[0m: usage = \" ./renameTool.sh \e[38;5;226m <ToolNAME>  <newNAME>\e[0m\" "	
	exit
    fi

    if [ "$newname" == "" ]
    then
	echo -e "\e[38;5;196mError no name given\e[0m: usage = \" ./renameTool.sh \e[38;5;226m <ToolNAME>  <newNAME>\e[0m\" "	
	exit
    fi
    
fi


if  [ -d $newname ]
then
    clear
    echo -e "\e[38;5;196mError Tool \e[38;5;226m$newname \e[38;5;196malready exists \e[0m"
    exit
fi



if [ -d $tool ]
then
    clear
    mv $tool $newname
    mv $newname/$tool.h $newname/$newname.h
    sed -i s:$tool:$newname:g $newname/$newname.h
    mv $newname/$tool.cpp $newname/$newname.cpp
    sed -i s:$tool:$newname:g $newname/$newname.cpp
    sed -i s:$tool.h:$newname.h:g Unity.h
    sed -i s:\"$tool\":\"$newname\":g Factory/Factory.cpp
    sed -i s:"$tool;":"$newname;":g Factory/Factory.cpp
    echo -e "\e[38;5;226mTool $tool renamted to $newname\e[0m"
else	
    clear
    echo -e "\e[38;5;196mError Tool \e[38;5;226m$tool \e[38;5;196m doesnt exists \e[0m"
    exit
fi

