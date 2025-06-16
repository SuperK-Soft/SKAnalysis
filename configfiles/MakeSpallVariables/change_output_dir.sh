#!/bin/bash

#OLDDIRNAME="run870"
if [ $# -lt 2 ]; then
	echo "usage: $0 <old_dir_name> <new_dir_name>"
	exit 1
fi

#for file in *; do
#        sed -i "s%/home/moflaher/%/home/moflaher/%g" ${file}
#done

#for file in *; do
#        sed -i "s%/disk02/usr6/moflaher/%/disk02/usr6/moflaher/%g" ${file}
#done
 
OLDDIRNAME=$1
NEWDIRNAME=$2
for file in *; do
        sed -i "s%/disk02/usr6/moflaher/${OLDDIRNAME}%/disk02/usr6/moflaher/${NEWDIRNAME}%g" ${file}
done
