# EXO DiPhoton Analysis Code

## To build

```
# set architecture
## bash
export SCRAM_ARCH=slc6_amd64_gcc493  
## tcsh
setenv SCRAM_ARCH slc6_amd64_gcc493

# setup cmssw release
cmsrel CMSSW_7_6_4  
cd CMSSW_7_6_4/src  
cmsenv  

# clone repository
git clone git@github.com:cms-exotica-diphotons/diphoton-analysis  

# compile code
cd diphoton-analysis  
scram b -j 16
```

## Fake Rate Analysis

Fake rate twiki  
https://twiki.cern.ch/twiki/bin/viewauth/CMS/ExoPhotonFakeRate13TeV

### EDAnalyzers  
ExoDiPhotonFakeRateAnalyzer  
ExoDiPhotonMCFakeRateRealTemplateAnalyzer  
ExoDiPhotonMCFakeRateClosureTestAnalyzer

### Analysis and plotting code  
FakeRateAnalysis

## Signal Analysis

### RSG signal

#### EDAnalyzer  
ExoDiPhotonRSGSignalAnalyzer

## Background Analysis
ExoDiPhotonBackgroundAnalyzer

## Misc.

### Tools

Scripts for general analysis
