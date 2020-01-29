#!/bin/tcsh

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

setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/Parallel/${year}/${model}/jobs"

setenv runumberslist ` ls ${workpath} | grep .sub `

foreach run  ($runumberslist)

#echo ${run}
chmod 755 ${workpath}/${run}

echo "Sending ${run}"
#bsub -q 8nh -o /tmp/junk ${workpath}/${run}
condor_submit ${workpath}/${run}
echo "condor_submit ${workpath}/${run}"
#bsub -q 8nh -o ${workpath}/../logs/${run}.txt ${workpath}/${run}
#echo "bsub -q 8nh -o /tmp/junk ${workpath}/${run}"

end

end 

end




