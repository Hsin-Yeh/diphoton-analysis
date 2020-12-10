#!/bin/tcsh

#setenv combmode "samefunfit" 
setenv combmode "diffunfit" 

setenv nominalmodel "dijet"

#setenv models "pow Laurent Exponential VVdijet expow invpow moddijet"
#2017 Models
setenv models "Laurent Exponential"
#2018 Models 
#setenv models "Laurent PowerLaw Exponential VVdijet"

#setenv masses `seq 600 100 5000`
setenv masses "600 700 800 900 1000 1100 1200 1500 1800 2100 2400 2700 3000 3500 4000 4500 5000 5500 6000" 
#setenv masses "600"

setenv couplings "kMpl001 kMpl01 kMpl02"
#setenv couplings "kMpl001"

#I will change it one at a time
setenv insigname "grav"

# Years 
setenv years "2017"
#setenv years "2018"

setenv musinjected `seq 1 5`
#setenv musinjected "5 10 15"
#setenv musinjected `seq 1 3`
#setenv musinjected `seq 1 1`
setenv ntoys 100
setenv theseed 397

#This is for the number of jobs per clusterid
set jobsperclusterchoice=10

setenv PWD `pwd`

# Year first
foreach year ($years)
echo "------------------------"
echo "Year ${year}"

# in/out
setenv inoutpath "/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis/output/${year}/combine_bias/${insigname}/${combmode}"

# Coupling now
foreach coup ($couplings)
echo "------------------------"
echo "Coupling ${coup}"

#mus injected
foreach muin ($musinjected)
echo "------------------------"
echo "Muin ${muin}"

# Starting the loop through all models
foreach model ($models)
echo "===================================================================================="
echo "Model $model"

#masses
foreach mass ($masses)
echo "------------------------"
echo "Mass ${mass}"

#Create local structure for the output
rm -rf ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/output ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/jobs ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/logs
mkdir -p ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/output ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/jobs ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/logs
chmod 755 -R ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/output ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/jobs ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/logs

#Output of the job will be here
rm -rf ${inoutpath}/${coup}/mu${muin}/${model}/mass${mass}
mkdir -p ${inoutpath}/${coup}/mu${muin}/${model}/mass${mass}

#foreach batch (`seq 0 100`)
foreach batch (`seq 0 0`)

echo '+JobFlavour = "tomorrow" ' > bias_$batch.sub
#echo '+JobFlavour = "microcentury" ' > bias_$batch.sub
echo ' ' >> bias_$batch.sub
echo "executable  = ${PWD}/setupCombineBias.sh" >> bias_$batch.sub
#echo "arguments   = "'$(ClusterID) $(ProcId)'" ${ncut} ${thick} ${file} ${thicknum} " >> bias_${file}.sub
echo "arguments   = "'$(ClusterID) $(ProcId)'" ${inoutpath} ${model} ${ntoys} ${coup} ${insigname} ${theseed} ${year} ${muin} ${mass} ${combmode} ${nominalmodel} "'$(ProcID)'" " >> bias_$batch.sub
echo "output      = ${PWD}/${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/logs/bias_"'$(infile)'".out " >> bias_$batch.sub
echo "error       = ${PWD}/${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/logs/bias_"'$(infile)'".err " >> bias_$batch.sub
echo "log         = ${PWD}/${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/logs/bias_"'$(infile)'"_htc.log " >> bias_$batch.sub
#echo "output      = ${PWD}/${year}/${insigname}/${combmode}/${coup}/bias_"'$(infile)'".out " >> bias_$batch.sub
#echo "error       = ${PWD}/${year}/${insigname}/${combmode}/${coup}/bias_"'$(infile)'".err " >> bias_$batch.sub
#echo "log         = ${PWD}/${year}/${insigname}/${combmode}/${coup}/bias_"'$(infile)'"_htc.log " >> bias_$batch.sub

#echo 'requirements = (OpSysAndVer =?= "CentOS7") ' >> bias_${batch}.sub
echo 'max_retries = 1' >> bias_$batch.sub

rm voodoo
touch voodoo
foreach jobspercluster (`seq 1 10`)
#foreach jobspercluster (`seq 1 20`)
set num=`expr ${batch} \* 10  + ${jobspercluster} `
echo -n "${num} " >> voodoo
end

setenv batchfilelist `cat voodoo`
echo "queue infile in (${batchfilelist}) " >> bias_$batch.sub

mv bias_$batch.sub ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/jobs/bias_$batch.sub
chmod 755 ${year}/${insigname}/${combmode}/${coup}/mu${muin}/${model}/mass${mass}/jobs/bias_$batch.sub

echo bias_$batch.sub

end 

end

end

end

end 

end












