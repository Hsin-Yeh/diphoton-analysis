#!/bin/bash

echo "The script starts now."

clusterid=${1}
procid=${2}
THEYEAR=${3}
THEINPATH=${4}
THEMODEL=${5}
THETOYSTART=${6}
THETOYEND=$(($THETOYSTART + 1))

echo "System: "
uname -a

cd /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis
eval `scramv1 runtime -sh`
cd -

export PWD=`pwd`

models="pow expow invpow invpowlin moddijet " 

for wind in 500_550 550_600 600_650 650_700 700_750 750_800 800_900 900_1000 1000_1200 1200_1800 1800_2500 2500_3500 3500_4500 4500_5500;do  mkdir -p ${PWD}/biasfiles/${THEYEAR}/${THEMODEL}/${wind}; done;

mkdir -p output/bkg/biasfiles
cp /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/*.rs .

echo "bias_study.exe ${THEYEAR} \"${THEINPATH}\" \"${PWD}\" ${THEMODEL} ${THETOYSTART} ${THETOYEND}"
bias_study.exe ${THEYEAR} "${THEINPATH}" "${PWD}" ${THEMODEL} ${THETOYSTART} ${THETOYEND}

for wind in 500_550 550_600 600_650 650_700 700_750 750_800 800_900 900_1000 1000_1200 1200_1800 1800_2500 2500_3500 3500_4500 4500_5500;do  cp ${PWD}/biasfiles/${year}/${model}/${wind}/*.root /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/output/bkg/biasfiles/${year}/${model}/${wind}/.; done;

#cp output/bkg/biasfiles/*.root /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/Parallel/${year}/${model}/output/.



 
