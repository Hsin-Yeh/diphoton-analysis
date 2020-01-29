#!/bin/tcsh

cd /afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/Parallel
eval `scramv1 runtime -csh`
cd -

#The final output files will be here
setenv finalout "/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/output/bkg/biasfiles/2017/finaloutput"

rm -rf ${finalout}
mkdir ${finalout}

#This is where is the input from the jobs
setenv input "/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/output/bkg/biasfiles/2017/expow"

# Models we are studying
#setenv models "pow expow invpow invpowlin moddijet " 
setenv models "expow" 

# Years 
setenv years "2017"

# Year first
foreach year ($years)
echo "------------------------"
echo "Year ${year}"

# Starting the loop through all models
foreach model ($models)
echo "===================================================================================="
echo "Model $model" 

hadd -f ${finalout}/tree_bias_${model}_${year}.root ${input}/tree_bias_*.root

end

end 





