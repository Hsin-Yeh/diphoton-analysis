#!/bin/bash

echo "The script starts now."

clusterid=${1}
procid=${2}
THEYEAR=${3}
THEINPATH=${4}
THESAMPLE=${5}

echo "System: "
uname -a

cd ${THEINPATH}
eval `scramv1 runtime -sh`
cd -

export PWD=`pwd`

mkdir -p ${PWD}/input/${THEYEAR}/workspaces
mkdir -p ${PWD}/input/${THEYEAR}/responsefcnfit
mkdir -p ${PWD}/input/${THEYEAR}/mgen

cp ${THEINPATH}/HighMassGG.rs .

echo "prepareWorkspaces.exe ${THEYEAR} \"${THEINPATH}/input/trees\" \"${PWD}/input\" ${THESAMPLE}"
prepareWorkspaces.exe ${THEYEAR} "${THEINPATH}/input/trees" "${PWD}/input" ${THESAMPLE}

cp ${PWD}/input/${THEYEAR}/workspaces/*_${THEYEAR}.root ${THEINPATH}/input/${THEYEAR}/workspaces/.
cp ${PWD}/input/${THEYEAR}/responsefcnfit/*.png ${THEINPATH}/input/${THEYEAR}/responsefcnfit/.
cp ${PWD}/input/${THEYEAR}/mgen/*.png ${THEINPATH}/input/${THEYEAR}/mgen/.
cp ${PWD}/input/${THEYEAR}/*.txt ${THEINPATH}/input/${THEYEAR}/.




