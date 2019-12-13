#!/bin/tcsh

#It is advantageous to submit multiple jobs as a single cluster because:
# 1. Only one copy of the checkpoint file is needed to represent all jobs 
#    in a cluster until they begin execution. 
# 2. There is much less overhead involved for Condor to start the next job 
#    in a cluster than for Condor to start a new cluster. This can make a big 
#    difference when submitting lots of short jobs.

# Models we are studying
setenv models "pow expow invpow invpowlin moddijet " 

# Years 
setenv years "2017"

# in/out
setenv inoutpath "/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_10_2_13/src/diphoton-analysis/output/bkg"

# windows for the output structure
setenv windows "500_550 550_600 600_650 650_700 700_750 750_800 800_900 900_1000 1000_1200 1200_1800 1800_2500 2500_3500 3500_4500 4500_5500"

setenv alljobs "1000"
#This is for the number of jobs per clusterid
set jobsperclusterchoice=10
#Since all jobs are 1000 and we want 10 jobs per cluster we split 
#1000 jobs in 100 batches. 
# SET BATCHES ON THE LOOP BY HAND

setenv PWD `pwd`

#will need the range directory for so many files
rm -rf ${inoutpath}/biasfiles

# Year first
foreach year ($years)
echo "------------------------"
echo "Year ${year}"

# Starting the loop through all models
foreach model ($models)
echo "===================================================================================="
echo "Model $model" 

#Create local structure for the output
rm -rf ${year}/${model}/output ${year}/${model}/jobs ${year}/${model}/logs
mkdir -p ${year}/${model}/output ${year}/${model}/jobs ${year}/${model}/logs
chmod 755 -R ${year}/${model}/output ${year}/${model}/jobs ${year}/${model}/logs

foreach wind ($windows)
mkdir -p ${inoutpath}/biasfiles/${year}/${model}/${wind}
end

#foreach file (`ls ${filelistpath}`)
foreach batch (`seq 0 99`)

echo '+JobFlavour = "tomorrow" ' > bias_$batch.sub
echo ' ' >> bias_$batch.sub
echo "executable  = ${PWD}/setupBias.sh" >> bias_$batch.sub
#echo "arguments   = "'$(ClusterID) $(ProcId)'" ${ncut} ${thick} ${file} ${thicknum} " >> bias_${file}.sub
echo "arguments   = "'$(ClusterID) $(ProcId)'" ${year} ${inoutpath} ${model} "'$(infile)'" " >> bias_$batch.sub
echo "output      = ${PWD}${year}/${model}/logs/bias_"'$(infile)'".out " >> bias_$batch.sub
echo "error       = ${PWD}${year}/${model}/logs/bias_"'$(infile)'".err " >> bias_$batch.sub
echo "log         = ${PWD}${year}/${model}/logs/bias_"'$(infile)'"_htc.log " >> bias_$batch.sub
echo 'requirements = (OpSysAndVer =?= "CentOS7") ' >> bias_${batch}.sub
echo 'max_retries = 1' >> bias_$batch.sub

rm voodoo
touch voodoo
#foreach jobspercluster (`seq 1 ${jobsperclusterchoice}`)
foreach jobspercluster (`seq 1 10`)
set num=`expr ${batch} \* 10  + ${jobspercluster} ` 
echo -n "${num} " >> voodoo
end 

setenv batchfilelist `cat voodoo`
echo "queue infile in (${batchfilelist}) " >> bias_$batch.sub

mv bias_$batch.sub ${year}/${model}/jobs/bias_$batch.sub
chmod 755 ${year}/${model}/jobs/bias_$batch.sub

echo bias_$batch.sub

rm voodoo

end

end

end
~

