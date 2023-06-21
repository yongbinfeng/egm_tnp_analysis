#!/bin/bash 

# if [[ `uname -n` == *"lxplus"* ]] || [[ `uname -n` == *"cmswmass2"* ]]; then 
# ## install gcc 4.9
#     echo sourcing all the things
#     optVer=x86_64-slc6-gcc49-opt
#     optVer=x86_64-slc6-gcc48-opt

#     . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh 

# ## add python lib
#     . /opt/rh/python27/enable
#     export LD_LIBRARY_PATH=/opt/rh/python27/root/usr/lib64/:$LD_LIBRARY_PATH

# fi

#export PYTHONPATH=.:$PYTHONPATH
# next line allows the setup to be run from anywhere, but will still use the location of the setup.sh file to set the python path
export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ):$PYTHONPATH"

#/data/shared/singularity/pythonrootdevf32.sif

## for python 2
#/data/shared/singularity/python2rootdevf32.sif

## setup for cmswmass2
##=============================================
## cmssw-cc7
## source /cvmfs/cms.cern.ch/cmsset_default.sh
## cmsrel CMSSW_11_2_0
## cd CMSSW_11_2_0
## cmsenv
#alias python=python3
#source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos8-gcc10-opt/setup.sh
