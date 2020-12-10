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

mkdir -p ${PWD}/output/${THEYEAR}/workspaces
mkdir -p ${PWD}/input/${THEYEAR}/FinalConvShapes

cp ${THEINPATH}/HighMassGG.rs .

echo "prepareShapes.exe ${THEYEAR} \"${THEINPATH}/input/trees\" \"${THEINPATH}/input/${THEYEAR}/workspaces\" \"${PWD}/output/${THEYEAR}/workspaces\" ${THESAMPLE}"
prepareShapes.exe ${THEYEAR} "${THEINPATH}/input/trees" "${PWD}/input/${THEYEAR}/workspaces" "${PWD}/output/${THEYEAR}/workspaces" ${THESAMPLE}

cp ${PWD}/output/${THEYEAR}/workspaces/*_${THEYEAR}.root ${THEINPATH}/output/${THEYEAR}/workspaces/.
cp ${PWD}/output/${THEYEAR}/FinalConvShapes/*.png ${THEINPATH}/output/${THEYEAR}/FinalConvShapes/.
