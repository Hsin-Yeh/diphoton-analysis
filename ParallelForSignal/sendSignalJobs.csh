#!/bin/tcsh

#Years to process
#No signal for 2016 for now
setenv years "2017 2018"

# Year first
foreach year ($years)
echo "------------------------"
echo "Year ${year}"

setenv workpath "/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis/ParallelForSignal/${year}/jobs"

setenv runumberslist ` ls ${workpath} | grep .sub `

foreach run  ($runumberslist)

#echo ${run}
chmod 755 ${workpath}/${run}

echo "Sending ${run}"
condor_submit ${workpath}/${run}
echo "condor_submit ${workpath}/${run}"

end

end 
