#!/bin/tcsh

cd /afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis/Parallel
eval `scramv1 runtime -csh`
cd -

# Models we are studying
#setenv models "pow expow invpow invpowlin moddijet " 
#setenv models "expow" 
#setenv models "Laurent invpowlin Atlas Exponential Expow"
#setenv models "Laurent PowerLaw invpowlin Atlas Exponential Expow"
# setenv models "Laurent PowerLaw Atlas Exponential Expow invpow invpowlin"
#setenv models "Atlas Exponential Expow invpow invpowlin"
# setenv models "Laurent Exponential Expow"
#setenv models "Exponential invpow"
# setenv models "pow"
setenv models "Laurent pow PowerLaw Exponential VVdijet expow dijet"
# Year 
setenv years "2017"
#setenv years "2018"

# Year first
foreach year ($years)
echo "------------------------"
echo "Year ${year}"

#The final output files will be here
setenv finalout "/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis/output/${year}/bkg/biasfiles/finaloutput"

rm -rf ${finalout}
mkdir ${finalout}

#This is where is the input from the jobs
setenv input "/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis/output/${year}/bkg/biasfiles"

# Starting the loop through all models
foreach model ($models)
echo "===================================================================================="
echo "Model $model" 

hadd -f ${finalout}/tree_bias_${model}_${year}.root ${input}/${model}/tree_bias_*.root

end

end 

