#!/bin/bash

if  command -v dialog &> /dev/null
then
    
    if [ "$1" == "" ]
    then

	dialog --inputbox "Enter new Tool name:" 0 0 2>answer
	
	if [ $? -eq 0 ]
	then
	    dir=`cat answer`

	    if [ "$dir" == "" ]
		then
		 rm -f answer
                clear
		echo -e "\e[38;5;196mERROR!!! no Tool name given \e[0m"
                exit
	    fi
	    
	    dialog --radiolist "Select Tool Template:" 0 0 0 \
		"<Blank>" . on \
		`for template in \`ls template/* | grep .cpp | sed s:template/MyTool:: | sed s:.cpp::\`
	   do
	   echo $template . off
	   done
          ` 2>answer
	    
	    if [ $? -eq 0 ]
	    then
		template=`cat answer`
		
	    else
		rm -f answer
		clear
		exit
	    fi
	else
            rm -f answer
            clear
            exit
	fi
	rm -f answer
	clear
	
	if [ "$template" == "<Blank>" ]
	then
	    template=""
	fi
	
    else
	
	dir=$1
	template=$2
    fi

else
    
    dir=$1
    template=$2
fi


    if [ "$dir" != "" ]
    then
	
	if [ -d $dir ]
	then
	    echo -e "\e[38;5;196mERROR!!! Tool directory already exists \e[0m"
	else
	    
	    mkdir $dir
	    more template/MyTool$template.h | sed s:MyTool$template:$dir:g | sed s:MYTOOL${template}_H:${dir}_H:g > ./$dir/$dir.h
	    more template/MyTool$template.cpp | sed s:MyTool$template:$dir:g | sed s:MyTool$template\(\):$dir\(\):g > ./$dir/$dir.cpp
	    more template/README.md | sed s:MyTool:$dir:g | sed s:MyTool\(\):$dir\(\):g > ./$dir/README.md
	    echo "#include \"$dir.h\"" >>Unity.h
	    
	    while read line
	    do
		if [ "$line" == "return ret;" ]
		then
		    echo "  if (tool==\""$dir"\") ret=new "$dir";" >>Factory/Factory.cpp.tmp
		fi
		echo "$line" >>Factory/Factory.cpp.tmp
	    done < Factory/Factory.cpp
	    mv Factory/Factory.cpp.tmp Factory/Factory.cpp
	fi
    else
	
	echo -e "\e[38;5;196mError no name given \e[0m: usage = \"./newTool.sh \e[38;5;226m <ToolNAME> \e[38;5;46m <TemplateNAME> \e[0m\"  if <TemplateName> is blank then blank template is used"
	echo -e "Valid tools template names are:"
	echo -e "\e[38;5;46m<BLANK>"
	for name in `ls template/ |grep '\.h' |grep -v "h~" |sed s:"\.h"::|sed s:"MyTool"::`
	do
	    echo $name
	done    
	echo -e "\e[0m"
    fi
    
