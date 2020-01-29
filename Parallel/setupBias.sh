#!/bin/bash

echo "The script starts now."

clusterid=${1}
procid=${2}
THEYEAR=${3}
THEINPATH=${4}
THEMODEL=${5}
THETOYSTART=${6}
THETOYEND=$(($THETOYSTART + 99))

echo "System: "
uname -a

cd /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis
eval `scramv1 runtime -sh`
cd -

export PWD=`pwd`

#models="pow expow invpow invpowlin moddijet" 
models="expow" 

mkdir -p ${PWD}/biasfiles/${THEYEAR}/${THEMODEL}
mkdir -p output/bkg/biasfiles
mkdir -p Tools/json

cp /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/Tools/json/windowsnorm_${THEYEAR}.json Tools/json
cp /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/HighMass-hgg_models_Bkg.bias.rs .
cp ${THEINPATH}/bkg_${THEMODEL}_${THEYEAR}.root ${PWD}/output/bkg/.

echo "bias_study.exe ${THEYEAR} \"${PWD}/output/bkg\" \"${PWD}/biasfiles/${THEYEAR}\" ${THEMODEL} ${THETOYSTART} ${THETOYEND} readJson"
bias_study.exe ${THEYEAR} "${PWD}/output/bkg" "${PWD}/biasfiles/${THEYEAR}" ${THEMODEL} ${THETOYSTART} ${THETOYEND} readJson

cp ${PWD}/biasfiles/${THEYEAR}/${THEMODEL}/*.root /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/output/bkg/biasfiles/${THEYEAR}/${THEMODEL}/.

#cp output/bkg/biasfiles/*.root /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/Parallel/${year}/${model}/output/.



 
