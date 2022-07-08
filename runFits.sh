#!/bin/bash 
## 
## for i in tracking idip trigger iso; do for j in plus minus both; do for k in BtoF GtoH BtoH; do echo bash runFits.sh -t $i -c $j -e $k; done ; done ; done

## tnptypes: idip, trigger, tracking, iso, isonotrig, reco
## charges : plus, minus
## eras    : BtoF, GtoH, BtoH
while getopts c:t:e: flag
do
    case "${flag}" in
        c) charge=${OPTARG};;
        t) tnptype=${OPTARG};;
        e) era=${OPTARG};;
    esac
done

#mctruth=''
mcTruth='--mcTruth'

echo "charge: $charge";
echo "tnptype: $tnptype";
echo "era: $era";

python3 tnpEGM_fitter.py --era ${era} --flag mu_${tnptype}_${charge} ${mcTruth} --createBins  
python3 tnpEGM_fitter.py --era ${era} --flag mu_${tnptype}_${charge} ${mcTruth} --createHists 
python3 tnpEGM_fitter.py --era ${era} --flag mu_${tnptype}_${charge} ${mcTruth} --doFit
python3 tnpEGM_fitter.py --era ${era} --flag mu_${tnptype}_${charge} ${mcTruth} --doFit --mcSig --altSig
python3 tnpEGM_fitter.py --era ${era} --flag mu_${tnptype}_${charge} ${mcTruth} --doFit         --altSig

# sum up everything
python3 tnpEGM_fitter.py --era ${era} --flag mu_${tnptype}_${charge} ${mcTruth} --sumUp
# 

