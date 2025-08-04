#!/bin/bash

#Application path location of applicaiton

Dependencies=`pwd`/Dependencies

export LD_LIBRARY_PATH=`pwd`/lib:${Dependencies}/TMVA/lib:${Dependencies}/Kirk:${Dependencies}/relic_sk4:$LD_LIBRARY_PATH

for folder in `ls -d ${PWD}/UserTools/*/ `
do
    export PYTHONPATH=$folder:${PYTHONPATH}
done

export SEGFAULT_SIGNALS="all"

# new location of secret sauce on new sukap
export RFA_ROOT='/opt/FJSVrdass/lib'

# we had to add this after switching to the warwick ATMPD...
export SKPATH=${SKOFL_ROOT}/const:${SKOFL_ROOT}/const/lowe:${ATMPD_ROOT}/const
