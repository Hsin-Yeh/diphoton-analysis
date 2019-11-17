// creates combine datacard for the resonant analysis
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "TH1.h"
#include "TFile.h"

#include "diphoton-analysis/Tools/interface/sampleList.hh" // integrated luminosity only defined once
#include "diphoton-analysis/Tools/interface/utilities.hh"

class nuisance {

public:
  nuisance(std::string syst, std::string distribution, std::vector<std::string> contribution) :
    m_syst(syst),
    m_distribution(distribution),
    m_contribution(contribution) {}
  std::string m_syst;
  std::string m_distribution;
  std::vector<std::string> m_contribution;
};

double getYield(const std::string& region, const std::string& sample, const std::string& datacardYear, double& yieldError, const TF1 * scaleFactor = 0);
std::string getDiphotonYieldVariations(const std::string& region, const std::string& variation);
void makeDatacard(const std::string& signalPoint, const std::string& region, const std::string& datacardYear);

std::string datacardYear;

int main(int argc, char *argv[])
{

  if(argc!=2) {
    std::cout << "Syntax: makeDatacardRes.exe [2016/2017/2018] \n";
    return -1;
  }
  datacardYear = argv[1];

  std::vector<std::string> samples = getSampleListForFit();

  for(auto isample : samples) {
    
    //We will process one year each time
    if ( isample.find(datacardYear) == std::string::npos ) continue; 
    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;

    makeDatacard( isample, "EBEB", datacardYear);
    makeDatacard( isample, "EBEE", datacardYear);

  }

}

void makeDatacard(const std::string& signalPoint, const std::string& region, const std::string& datacardYear)
{

  //Systematics
  std::string lumiError = std::to_string(1 + luminosityErrorFrac[datacardYear]);
  //We have two processes, sig and bkg. bkg is estimated from data so 
  //a dash is there. 
  nuisance lumi("lumi", "lnN", {lumiError, "-"});
  nuisance pileup("pileup", "lnN", {"1.05", "-"});
  //nuisance eff("eff", "lnN", {"1.1", "1.05"});
  nuisance thetaSmear("thetaSmear$CHANNEL", "param", {"0.0", "1.0"});
  

  std::ofstream output;
  std::string filename("datacards/");
  filename+=signalPoint;
  filename+="_";
  filename+=datacardYear;
  filename+="_" + region;
  filename+=".dat";
  output.open(filename);
  if (output.is_open()) {

    // header
    output << "imax * " << nchannels << " number of channels" << std::endl;
    output << "jmax * " << nbackgrounds << " number of backgrounds" << std::endl;
    output << "kmax * " <<  << " number of nuisance parameters" << std::endl;
    output << "-------------------" << std::endl;

    output << "shapes data_obs * bkg_" << datacardYear << ".root " << " HLFactory_ws:Data_$CHANNEL" << endl;
    output << "shapes sig * SignalParametricShapes_ws.root " << " model_signal:SignalShape_kMpl001_$CHANNEL" << endl;
    output << "shapes bkg * bkg_"<< datacardYear << ".root " << " HLFactory_ws:PhotonsMassBkg_DiJet_$CHANNEL" << endl;

    output << "-------------------" << std::endl;

    output << "---------------" << std::endl;
    output << "bin          " << region << std::endl;
    output <<  "observation   -1"  <<  std::endl;
    output << "------------------------------" << std::endl;
    output << "bin          " << region << << "     " << region << std::endl;
    output << "process                 sig      bkg" << std::endl;
    output << "process                   0        1" << std::endl;
    output << "rate                     -1       -1 " << std::endl;

    output << "\n-------------------" << std::endl;
    
    // output systematics
    for(auto inuisance : allNuisances) {
      output << inuisance.m_syst  << "      " << inuisance.m_distribution << "  ";
      for( auto icontrib : inuisance.m_contribution) {
	output << icontrib << " ";
      }
      output << "\n";
    }

  } // closes is_open()
  else {
    std::cout << "Could not open output file " << filename << std::endl;
  }
}



