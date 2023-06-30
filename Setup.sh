#!/bin/bash

#Application path location of applicaiton

Dependencies=`pwd`/Dependencies

export LD_LIBRARY_PATH=`pwd`/lib:${Dependencies}/TMVA/lib:$LD_LIBRARY_PATH

for folder in `ls -d ${PWD}/UserTools/*/ `
do
    export PYTHONPATH=$folder:${PYTHONPATH}
done

export SEGFAULT_SIGNALS="all"

# new location of secret sauce on new sukap
export RFA_ROOT='/opt/FJSVrdass/lib'
