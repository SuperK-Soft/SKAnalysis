#!/bin/bash

#Application path location of applicaiton

Dependencies=`pwd`/Dependencies

export LD_LIBRARY_PATH=`pwd`/lib:${Dependencies}/TMVA/lib:$LD_LIBRARY_PATH

for folder in `ls -d ${PWD}/UserTools/*/ `
do
    export PYTHONPATH=$folder:${PYTHONPATH}
done

export SEGFAULT_SIGNALS="all"
export LD_LIBRARY_PATH=$PWD/Dependencies/relic_sk4_ana/data_reduc/third/lib:$LD_LIBRARY_PATH
source /usr/local/sklib_gcc8/bashenv_gcc8_skofl_21a+atmpd_21a
source /usr/local/sklib_gcc8/skofl_21a/env.sh

export NEUTRON=1

export PYTHONPATH=$PYTHONPATH:/home/giampaol/miniconda3/envs/spectral/lib/python3.8/site-packages/
export PYTHONPATH=$PYTHONPATH:/home/giampaol/miniconda3/lib/python3.8/site-packages/
export PYTHONPATH=/home/mattnich/.local/lib/python3.8/site-packages:$PYTHONPATH
