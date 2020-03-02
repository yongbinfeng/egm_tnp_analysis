#!bin/bash 

if [[ `uname -n` == *"lxplus"* ]]; then 
## install gcc 4.9
    optVer=x86_64-slc6-gcc49-opt
    optVer=x86_64-slc6-gcc48-opt

    . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh 

## add python lib
    . /opt/rh/python27/enable
    export LD_LIBRARY_PATH=/opt/rh/python27/root/usr/lib64/:$LD_LIBRARY_PATH
fi

export PYTHONPATH=.:$PYTHONPATH
