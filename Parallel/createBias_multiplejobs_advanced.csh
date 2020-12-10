#!/bin/tcsh

#It is advantageous to submit multiple jobs as a single cluster because:
# 1. Only one copy of the checkpoint file is needed to represent all jobs 
#    in a cluster until they begin execution. 
# 2. There is much less overhead involved for Condor to start the next job 
#    in a cluster than for Condor to start a new cluster. This can make a big 
#    difference when submitting lots of short jobs.

# Models we are studying
#setenv models "pow expow invpow invpowlin moddijet" 
#setenv models "expow" 
# setenv models "Laurent PowerLaw Atlas Exponential Expow invpow invpowlin"
#setenv models "Atlas Exponential Expow invpow invpowlin"
# setenv models "Laurent Exponential Expow"
#setenv models "Exponential invpow"
# setenv models "pow"
setenv models "Laurent pow PowerLaw Exponential VVdijet expow dijet"
# Years 
setenv years "2017"
#setenv years "2018"

# When I wanted to reduce the total jobs from 1000 to 100, that is 
# from 100 batches times 10 jobs per batch to 10 batches times
# 10 jobs, I just changed the foreach line below with the batches, 
# nothing else. 
setenv alljobs "10" #was 100 
#This is for the number of jobs per clusterid
set jobsperclusterchoice=10
#Since all jobs are 1000 and we want 10 jobs per cluster we split 
#1000 jobs in 100 batches. 
# SET BATCHES ON THE LOOP BY HAND

setenv PWD `pwd`

# Year first
foreach year ($years)
echo "------------------------"
echo "Year ${year}"

# in/out
setenv inoutpath "/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis/output/${year}/bkg"

#will need the range directory for so many files
# rm -rf ${inoutpath}/biasfiles

# Starting the loop through all models
foreach model ($models)
echo "===================================================================================="
echo "Model $model" 

#Create local structure for the output
rm -rf ${year}/${model}/output ${year}/${model}/jobs ${year}/${model}/logs
mkdir -p ${year}/${model}/output ${year}/${model}/jobs ${year}/${model}/logs
chmod 755 -R ${year}/${model}/output ${year}/${model}/jobs ${year}/${model}/logs

rm -rf ${inoutpath}/biasfiles/${model}
mkdir -p ${inoutpath}/biasfiles/${model}

#foreach file (`ls ${filelistpath}`)
#foreach batch (`seq 0 99`)
foreach batch (`seq 0 9`)

echo '+JobFlavour = "tomorrow" ' > bias_$batch.sub
#echo '+JobFlavour = "microcentury" ' > bias_$batch.sub
echo ' ' >> bias_$batch.sub
echo "executable  = ${PWD}/setupBias.sh" >> bias_$batch.sub
#echo "arguments   = "'$(ClusterID) $(ProcId)'" ${ncut} ${thick} ${file} ${thicknum} " >> bias_${file}.sub
echo "arguments   = "'$(ClusterID) $(ProcId)'" ${year} ${inoutpath} ${model} "'$(infile)'" " >> bias_$batch.sub
echo "output      = ${PWD}/${year}/${model}/logs/bias_"'$(infile)'".out " >> bias_$batch.sub
echo "error       = ${PWD}/${year}/${model}/logs/bias_"'$(infile)'".err " >> bias_$batch.sub
echo "log         = ${PWD}/${year}/${model}/logs/bias_"'$(infile)'"_htc.log " >> bias_$batch.sub
#echo 'requirements = (OpSysAndVer =?= "CentOS7") ' >> bias_${batch}.sub
# echo '+MaxJobExecutionTime = 1800' >> bias_$batch.sub
echo 'max_retries = 1' >> bias_$batch.sub

rm voodoo
touch voodoo
#foreach jobspercluster (`seq 1 ${jobsperclusterchoice}`)
foreach jobspercluster (`seq 1 100 1000`)
set num=`expr ${batch} \* 1000  + ${jobspercluster} ` 
echo -n "${num} " >> voodoo
end 

setenv batchfilelist `cat voodoo`
echo "queue infile in (${batchfilelist}) " >> bias_$batch.sub

mv bias_$batch.sub ${year}/${model}/jobs/bias_$batch.sub
chmod 755 ${year}/${model}/jobs/bias_$batch.sub

echo bias_$batch.sub

#rm voodoo

end

end

end
