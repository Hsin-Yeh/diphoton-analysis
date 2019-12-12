#!/bin/bash

echo "The script starts now."

clusterid=${1}
procid=${2}
THEYEAR=${3}
THEINOUTPATH=${4}
THEMODEL=${5}
THETOYSTART=${6}
THETOYEND=$(($THETOYSTART + 1))

echo "System: "
uname -a

cd /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis
eval `scramv1 runtime -sh`
cd -

export PWD=`pwd`

mkdir -p output/bkg/biasfiles
cp /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/*.rs .

echo "bias_study.exe ${THEYEAR} \"${THEINOUTPATH}\" ${THEMODEL} ${THETOYSTART} ${THETOYEND}"
bias_study.exe ${THEYEAR} "${THEINOUTPATH}" ${THEMODEL} ${THETOYSTART} ${THETOYEND}

#cp output/bkg/biasfiles/*.root /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/Parallel/${year}/${model}/output/.



 
