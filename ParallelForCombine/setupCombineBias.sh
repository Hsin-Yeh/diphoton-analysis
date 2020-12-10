#!/bin/bash

echo "The script starts now."

clusterid=${1}
procid=${2}
THEINOUTPATH=${3}
THEMODEL=${4}
THENTOYS=${5}
THECOUP=${6}
THEINSIGNAME=${7}
THESEED=${8}
THEYEAR=${9}
THEMUIN=${10}
THEMASS=${11}
THECOMBCHOICE=${12}
THENOMINALMODEL=${13}
CURRENTTOY=${14}

export mainpath="/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis"

cd ${mainpath}
eval `scramv1 runtime -sh`
cd -

export PWD=`pwd`

cp ${mainpath}/datacards/${THEINSIGNAME}_${THECOUP}_${THEMODEL}_${THEYEAR}.dat . 
cp ${mainpath}/datacards/bkg_${THEMODEL}_*_${THEYEAR}.root .
cp ${mainpath}/datacards/${THEINSIGNAME}_${THECOUP}_${THEYEAR}.root .
cp ${mainpath}/datacards/${THEINSIGNAME}_${THECOUP}_${THENOMINALMODEL}_${THEYEAR}.dat . 
cp ${mainpath}/datacards/bkg_${THENOMINALMODEL}_*_${THEYEAR}.root .


if [ ${THECOMBCHOICE} == "samefunfit" ]; then
#combine -M FitDiagnostics ${THEINSIGNAME}_${THECOUP}_${THEMODEL}_${THEYEAR}.dat -t ${THENTOYS} --toysFrequentist --saveToys --expectSignal ${THEMUIN} -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --seed THESEED
echo "combine -M FitDiagnostics ${THEINSIGNAME}_${THECOUP}_${THEMODEL}_${THEYEAR}.dat -t ${THENTOYS} --toysFrequentist --saveToys --expectSignal ${THEMUIN} -m ${THEMASS} -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --seed -1"

combine -M FitDiagnostics ${THEINSIGNAME}_${THECOUP}_${THEMODEL}_${THEYEAR}.dat -t ${THENTOYS} --toysFrequentist --saveToys --expectSignal ${THEMUIN} -m ${THEMASS} -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --seed -1 -v -1

#The final file will be of the form: 
mv fitDiagnostics_${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS}.root fitDiagnostics_${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS}_${CURRENTTOY}.root

cp fitDiagnostics_${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS}_${CURRENTTOY}.root ${THEINOUTPATH}/${THECOUP}/mu${THEMUIN}/${THEMODEL}/mass${THEMASS}/.

rm -rf ${THEINSIGNAME}_${THECOUP}_${THEMODEL}_${THEYEAR}.dat bkg_${THEMODEL}_${THEYEAR}.root ${THEINSIGNAME}_${THECOUP}_${THEYEAR}.root higgsCombine_*.root fitDiagnostics_*.root

fi 

if [ ${THECOMBCHOICE} == "diffunfit" ]; then

combine -M GenerateOnly ${THEINSIGNAME}_${THECOUP}_${THEMODEL}_${THEYEAR}.dat -m ${THEMASS} -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --saveToys --expectSignal ${THEMUIN} -t ${THENTOYS} --seed -1 -v -1 --bypassFrequentistFit
#combine -M GenerateOnly ${THEINSIGNAME}_${THECOUP}_${THEMODEL}_${THEYEAR}.dat -m ${THEMASS} -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --saveToys --expectSignal ${THEMUIN} -t ${THENTOYS} --seed -1 -v -1 --toysFrequentist

toysfile=`ls |grep GenerateOnly`

#combine -M FitDiagnostics ${THEINSIGNAME}_${THECOUP}_${THENOMINALMODEL}_${THEYEAR}.dat --setRobustFitTolerance=1 -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --toysFile ${toysfile} -t ${THENTOYS} --saveWorkspace --robustFit=1 -v 5 --seed -1 -m ${THEMASS} --setParameters "PhotonsMass_bkg_dijet_linc_cat0=5.0,PhotonsMass_bkg_dijet_linc_cat1=8.0,PhotonsMass_bkg_dijet_logc_cat0=-1.0,PhotonsMass_bkg_dijet_logc_cat1=-1.0"

combine -M FitDiagnostics ${THEINSIGNAME}_${THECOUP}_${THENOMINALMODEL}_${THEYEAR}.dat -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --toysFile ${toysfile} -t ${THENTOYS} --saveWorkspace -v -1 --seed -1 -m ${THEMASS} --cminDefaultMinimizerStrategy=0
#combine -M FitDiagnostics ${THEINSIGNAME}_${THECOUP}_${THENOMINALMODEL}_${THEYEAR}.dat -n _${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS} --toysFile ${toysfile} -t ${THENTOYS} --saveWorkspace -v -1 --seed -1 -m ${THEMASS} --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance 0.01 --minos poi
mv fitDiagnostics_${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS}.root fitDiagnostics_${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS}_${CURRENTTOY}.root

cp fitDiagnostics_${THEINSIGNAME}_mu${THEMUIN}_${THECOUP}_${THEMODEL}_${THEMASS}_${CURRENTTOY}.root ${THEINOUTPATH}/${THECOUP}/mu${THEMUIN}/${THEMODEL}/mass${THEMASS}/.

rm -rf *.root *.dat

# combine -M FitDiagnostics grav_kMpl001_Laurent_2017.dat -n _grav_mu1_kMpl001_Laurent_4000 --toysFile higgsCombine_grav_mu1_kMpl001_Laurent_4000.GenerateOnly.mH4000.-785808385.root -t 100
# combine -M GenerateOnly grav_kMpl001_Laurent_2017.dat -m 4000 -n _grav_mu1_kMpl001_Laurent_4000 --saveToys --expectSignal 1 -t 100 --seed -1 -v -1 --bypassFrequentistFit

fi



