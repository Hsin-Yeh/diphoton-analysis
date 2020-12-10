#!/bin/tcsh

setenv year $1
setenv signal $2
setenv coupling $3

setenv datacardfile "/afs/cern.ch/work/h/hsinyeh/public/diphoton-analysis/CMSSW_10_2_13/src/diphoton-analysis/datacards/${signal}_${coupling}_dijet_${year}.dat"

#setenv masslist `seq 550 25 5000`
setenv masslist "550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000 1025 1050 1075 1100 1125 1150 1175 1200 1225 1250 1275 1300 1325 1350 1375 1400 1425 1450 1475 1500 1525 1550 1575 1600 1625 1650 1675 1700 1725 1750 1775 1800 1825 1850 1875 1900 1925 1950 1975 2000 2025 2050 2075 2100 2125 2150 2175 2200 2225 2250 2275 2300 2325 2350 2375 2400 2425 2450 2475 2500 2525 2550 2575 2600 2625 2650 2675 2700 2725 2750 2775 2800 2825 2850 2875 2900 2925 2950 2975 3000 3025 3050 3075 3100 3125 3150 3175 3200 3225 3250 3275 3300 3325 3350 3375 3400 3425 3450 3475 3500 3525 3550 3575 3600 3625 3650 3675 3700 3725 3750 3775 3800 3825 3850 3875 3900 3925 3950 3975 4000 4025 4050 4075 4100 4125 4150 4175 4200 4225 4250 4275 4300 4325 4350 4375 4400 4425 4450 4475 4500 4525 4550 4575 4600 4625 4650 4675 4700 4725 4750 4775 4800 4825 4850 4875 4900 4925 4950 4975 5000 5500 6000 6500 7000 7500 8000"

#setenv masslist "555 580 605 630"

# Method in this script are: AsymptoticLimits, ExpSignificance, ExpSignificanceWithPval, ObsSignificance, ObsSignificanceWithPval 
setenv methods "ObsSignificance"
# setenv methods "AsymptoticLimits"

foreach method ($methods)

rm -rf combineJobs13TeV/${year}/${signal}/${coupling}/${method}/All
mkdir -p combineJobs13TeV/${year}/${signal}/${coupling}/${method}/All

foreach mass ($masslist)

echo "====================================================================="
echo $mass

##
##    
##    Expected Significance 
##
##

if ($method == "ExpSignificance") then
  combine -d $datacardfile -M Significance -m $mass --signif --pval --cminDefaultMinimizerType=Minuit2 -t -1 --expectSignal=1.0 -n Expected
  mv higgsCombineExpected.Significance.mH$mass.root combineJobs13TeV/${year}/${signal}/${coupling}/${method}/All/.
endif

if ($method == "ObsSignificance") then
  combine -d $datacardfile -M Significance -m $mass --signif --pval --cminDefaultMinimizerType=Minuit2 -n Observed
  mv higgsCombineObserved.Significance.mH$mass.root combineJobs13TeV/${year}/${signal}/${coupling}/${method}/All/.
endif

#end of loop over masses
end

##
##    
##    Asymptotic Limits
##
##

if ($method == "AsymptoticLimits") then

rm finalResults
touch finalResults

foreach mass ($masslist)

echo "====================================================================="
echo $mass

if ($method == "AsymptoticLimits") then

#a priori limits. I see in the post-fit or a-posteriori expected limit weird 
#one and two sigma region above 1.2 TeV. So, I will go to a priori limits at the moment. 
combine -M AsymptoticLimits -m $mass -s -1 --bypassFrequentistFit $datacardfile > ${datacardfile}_results
mv higgsCombineTest.AsymptoticLimits.mH$mass.*.root combineJobs13TeV/${year}/${signal}/${coupling}/${method}/All/.

setenv obs    `cat ${datacardfile}_results  | grep  "Observed Limit:" | awk '{print $5}'`
setenv expM2s `cat ${datacardfile}_results  | grep  "Expected  2.5%:" | awk '{print $5}'`
setenv expM1s `cat ${datacardfile}_results  | grep  "Expected 16.0%:" | awk '{print $5}'`
setenv exp    `cat ${datacardfile}_results  | grep  "Expected 50.0%:" | awk '{print $5}'`
setenv expP1s `cat ${datacardfile}_results  | grep  "Expected 84.0%:" | awk '{print $5}'`
setenv expP2s `cat ${datacardfile}_results  | grep  "Expected 97.5%:" | awk '{print $5}'`

echo $mass $obs $expM2s $expM1s $exp $expP1s $expP2s
echo $mass $obs $expM2s $expM1s $exp $expP1s $expP2s >> finalResults

#end of loop over masses
end

cat finalResults | sort -n > finalResults2
mv finalResults2 combineJobs13TeV/${year}/${signal}/${coupling}/${method}/All/finalResults

endif

#end of loop over methods
end 




