#!/bin/tcsh

#It is advantageous to submit multiple jobs as a single cluster because:
# 1. Only one copy of the checkpoint file is needed to represent all jobs 
#    in a cluster until they begin execution. 
# 2. There is much less overhead involved for Condor to start the next job 
#    in a cluster than for Condor to start a new cluster. This can make a big 
#    difference when submitting lots of short jobs.

#Years to process
#No signal for 2016 for now
setenv years "2017 2018"

# Samples we are studying
setenv samples "RSGravitonToGammaGamma_kMpl001_M_750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_1000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_1250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_1500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_1750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_2000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_2250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_2500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_2750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_3000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_3250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_3500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_4000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl001_M_5000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_1000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_1250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_1500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_1750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_2000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_2250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_2500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_3000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_3500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_4000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_4250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_4500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_4750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_5000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_5250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_5500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_5750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_6000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_6500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_7000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl01_M_8000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_1000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_1250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_1500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_1750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_2000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_2250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_2500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_3000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_3500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_4000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_4500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_4750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_5000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_5250_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_5500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_5750_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_6000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_6500_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_7000_TuneCP2_13TeV_pythia8 RSGravitonToGammaGamma_kMpl02_M_8000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_1000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_1250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_1500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_1750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_2000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_2250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_2500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_2750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_3000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_3250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_3500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_4000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_4500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_0p014_M_5000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_1000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_1250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_1500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_1750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_2000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_2250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_2500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_3000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_3500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_4000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_4250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_4500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_4750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_1p4_M_5000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_1000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_1250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_1500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_1750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_2000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_2250_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_2500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_3000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_3500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_4000_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_4500_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_4750_TuneCP2_13TeV_pythia8 GluGluSpin0ToGammaGamma_W_5p6_M_5000_TuneCP2_13TeV_pythia8"

setenv inoutpath "/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis"

setenv PWD `pwd`

# Year first
foreach year ($years)
echo "------------------------"
echo "Year ${year}"

#Create local structure for jobs
rm -rf ${year}/output ${year}/jobs ${year}/logs
mkdir -p ${year}/output ${year}/jobs ${year}/logs
chmod 755 -R ${year}/output ${year}/jobs ${year}/logs

# Starting the loop through all models
foreach sam ($samples)
echo "===================================================================================="
echo "Sample $sam" 


echo '+JobFlavour = "workday" ' > $sam.sub
echo ' ' >> $sam.sub
echo "executable  = ${PWD}/setuSignalJobs.sh" >> $sam.sub
echo "arguments   = "'$(ClusterID) $(ProcId)'" ${year} ${inoutpath} ${sam}" >> $sam.sub
echo "output      = ${PWD}/${year}/logs/${sam}.out " >> $sam.sub
echo "error       = ${PWD}/${year}/logs/${sam}.err " >> $sam.sub
echo "log         = ${PWD}/${year}/logs/${sam}_htc.log " >> $sam.sub
echo 'requirements = (OpSysAndVer =?= "CentOS7") ' >> $sam.sub
echo 'max_retries = 1' >> $sam.sub
echo "queue 1 " >> $sam.sub

mv $sam.sub ${year}/jobs/$sam.sub
chmod 755 ${year}/jobs/$sam.sub

echo $sam.sub

end

end
