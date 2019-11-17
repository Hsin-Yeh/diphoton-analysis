#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "FWCore/Utilities/interface/Exception.h"
#include "diphoton-analysis/CommonClasses/interface/CrossSections.h"

//ROOT
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPaveText.h"

#include "diphoton-analysis/RooUtils/interface/RooSpline1D.h"
#include "RooRealVar.h"
#include "RooVoigtian.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"

//-----------------------------------------------------------------------------------
struct eff_reco {

  std::string kMpl;
  std::string M_bins;
  //Values for sel are: no_selection, isGood
  std::string sel; 
  //Values for cat are both, BB, BE. 
  std::string cat;
  //Efficiency of acceptance of efficiency times acceptance
  //depending on the the sel value
  double efforacc; 
  
};

//-----------------------------------------------------------------------------------
struct theGraphs{
  
  std::string coup;
  std::string cat;
  TGraphErrors* thegraph;
  TF1 * fun;
  
};

//-----------------------------------------------------------------------------------
struct xsec{
  
  std::string name;
  std::string M_bins;
  double val;
  double error ;
  
};

using namespace RooFit;

//-----------------------------------------------------------------------------------
//Declarations here definition after main
std::map<std::string, std::vector<eff_reco> > computefficiency(const std::string &year);
void plotefficiency(TCanvas* cc, std::map<std::string, std::vector<eff_reco> > coup_eff_reco, bool doave);
TGraphErrors* graph(std::string coup, std::map<std::string, std::vector<eff_reco> > quant);
std::vector<theGraphs> plot(TCanvas* cc, std::string coup, std::map<std::string, std::vector<eff_reco> > quantBB, std::map<std::string, std::vector<eff_reco> > quantBE , std::map<std::string, std::vector<eff_reco> > quantTotal, bool doave, std::string label);
std::map<std::string, std::vector<eff_reco> > reduce(std::map<std::string, std::vector<eff_reco> > coup_eff_reco, std::string cut);
TH1F * createHisto(TChain * newtree1, const std::string &histo_name, int nBins, double xMin, double xMax, std::string cut);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);
void writeToJson(std::map< std::string , std::vector<double> > valuestowrite, std::string outputfile);
std::map<std::string, std::vector<xsec> > loadXsections(const std::string &year);
TGraphErrors* getXsecGraph( std::vector<xsec> xsections);
RooSpline1D* graphToSpline(std::string name, TGraphErrors *graph, RooRealVar* MH, double xmin, double upperxmax);
void plotsplines(TCanvas* cc,  RooSpline1D * xsSplines, RooRealVar* MH, double xmin, double upperxmax, std::string label);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

  std::string region, year, inputdir, outputdir;

  if(argc!=2) {
    std::cout << "Syntax: signalNorm.exe [2016/2017/2018] " << std::endl;
    return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
  }

  //========================================================================
  //read the reduced input root trees
  // initForFit(inputdir);

  // include signal samples but not unskimmed data samples
  init(false, true);

  std::map<std::string, std::vector<xsec> > xsections = loadXsections(year);
  std::map<std::string, TGraphErrors*> grxs;
  std::map<std::string, RooSpline1D * > xsSplines;
  
  //Couplings for the moment
  std::vector<std::string> coups; coups.clear();
  coups.push_back("001");
  coups.push_back("01");
  coups.push_back("02");
  
  RooRealVar* MH = new RooRealVar("MH", "MH", 1000.);
  MH->setConstant();
  RooRealVar* kmpl = new RooRealVar("kmpl", "kmpl", 0.01);
  kmpl->setConstant();

  double upperxmax = 0.;
  double xmin = 230.;
  std::string plabel;
  
  std::map<std::string , TCanvas*> can;

  for(auto cp : coups) {

    can[cp] = new TCanvas(Form("can_coup_%s", cp.c_str()),Form("can_coup_%s", cp.c_str()));

    if ( cp == "001" ){upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4}";}
    else if ( cp == "01" ){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2}";}
    else if ( cp == "02" ){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2}";}
    else {
      std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
      exit(1);
    }
    
    grxs[cp] = getXsecGraph(xsections[cp]);
    RooSpline1D *xsSpline = graphToSpline(Form("fxs_%s",cp.c_str()), grxs[cp], MH , xmin, upperxmax);
    xsSplines[cp] = xsSpline;

    //For debugging to see what this spline looks like
    plotsplines(can[cp], xsSplines[cp],  MH ,xmin, upperxmax, plabel);
    can[cp]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/xsSplines_%s.png",cp.c_str()));

  }

  

  //========================================================================
  std::map<std::string, std::vector<eff_reco> > effreco = computefficiency(year);
  // std::map<std::string, std::vector<eff_reco> > accBB = reduce(effreco, "ABB");
  // std::map<std::string, std::vector<eff_reco> > accBE = reduce(effreco, "ABE");
  // std::map<std::string, std::vector<eff_reco> > accTotal = reduce(effreco, "ATotal");
  std::map<std::string, std::vector<eff_reco> > exABB = reduce(effreco, "exABB");
  std::map<std::string, std::vector<eff_reco> > exABE = reduce(effreco, "exABE");
  std::map<std::string, std::vector<eff_reco> > exATotal = reduce(effreco, "exATotal");

  TCanvas* cc001 = new TCanvas("cc001", "cc001");
  std::map<std::string, std::vector<theGraphs> > graphsofexA;

  graphsofexA["001"] = plot(cc001, "001", exABB, exABE, exATotal, false, "#varepsilon #otimes A");
  cc001->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exAsVsMass_001.png");

  TCanvas* cc01 = new TCanvas("cc01", "cc01");
  graphsofexA["01"] = plot(cc01, "01", exABB, exABE, exATotal, false, "#varepsilon #otimes A");
  cc01->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exAsVsMass_01.png");

  TCanvas* cc02 = new TCanvas("cc02", "cc02");
  graphsofexA["02"] = plot(cc02, "02", exABB, exABE, exATotal, false, "#varepsilon #otimes A");
  cc02->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exAsVsMass_02.png");

  //Graphs of exA to splines
  std::map<std::string , TCanvas*> cans;
  std::map<std::string , RooSpline1D *> exASplines;

  for(auto grs : graphsofexA) {
    for(auto gr : grs.second) {
      cans[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())] = new TCanvas(Form("%s_%s", gr.coup.c_str(), gr.cat.c_str()), Form("%s_%s", gr.coup.c_str(), gr.cat.c_str()));

      
      if ( gr.coup == "001" ){upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4}";}
      else if ( gr.coup == "01" ){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2}";}
      else if ( gr.coup == "02" ){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2}";}
      else {
	std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
	exit(1);
      }

      RooSpline1D *exASpline = graphToSpline(Form("fea_%s_%s",gr.coup.c_str(), gr.cat.c_str() ), gr.thegraph , MH , xmin, upperxmax);
      exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())] = exASpline;

      //For debugging to see what this spline looks like
      plotsplines(cans[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())], exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())],  MH ,xmin, upperxmax, plabel);

      cans[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exASplines_%s_%s.png",gr.coup.c_str(), gr.cat.c_str()));

		
      }//end of loop through graphs
  }//end of loop through couplings

  



  // plotefficiency(cc, effrecoBB, false);
  // cc->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exAsVsMass_001.png");
  // TCanvas* ccBE = new TCanvas("ccBE", "ccBE");
  // plotefficiency(ccBE, effrecoBE, false);
  // ccBE->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effVsMass_BE.png");
  // TCanvas* ccave = new TCanvas("ccave", "ccave");
  // plotefficiency(ccave, effreco, true);
  // ccave->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effVsMass_average.png");

  //========================================================================
  //Acceptance values
  std::vector<double> acc_val;
  acc_val.push_back(-1.057554393); 
  acc_val.push_back(-2.047333072); 
  acc_val.push_back(3.7439397255); 

  std::map< std::string, std::vector<double> > mapnamestoacc;
  mapnamestoacc["acc_EBEB_p2"] = acc_val;
  mapnamestoacc["acc_EBEB_p0"] = acc_val;

  writeToJson(mapnamestoacc, "Tools/json/acceptance.json");

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
std::map<std::string, std::vector<eff_reco> > computefficiency(const std::string &year)
{
  
  std::map<std::string, std::string> cuts;
  //Acceptance cuts
  cuts["ABB"] = "(Diphoton.Minv > 230 && Photon1.pt>125 && Photon2.pt>125 && Photon1.isEB && Photon2.isEB)";
  cuts["ABE"] = "(Diphoton.Minv > 330 && Photon1.pt>125 && Photon2.pt>125 && ( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE )))";
  cuts["ATotal"] = "(Diphoton.Minv > 230 && Photon1.pt>125 && Photon2.pt>125) && ( (Photon1.isEB && Photon2.isEB) || (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ) )";
  // std::cout << cuts["ATotal"] << std::endl;
  //Selection efficiency cuts. 
  cuts["eBB"] = "(Photon1.isEB && Photon2.isEB)*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  cuts["eBE"] = "( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ))*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  cuts["eTotal"] = "( (Photon1.isEB && Photon2.isEB) || (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ) )*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  // std::cout << cuts["eTotal"] << std::endl;
  //Acceptance times efficiency
  cuts["exABB"] = cuts["ABB"] + "*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)*weightAll";
  cuts["exABE"] = cuts["ABE"] + "*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)*weightAll";
  cuts["exATotal"] = "(" + cuts["ATotal"] + ")" + "*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)*weightAll";
  // std::cout << cuts["exATotal"] << std::endl;

  // cuts["BB"] = "isGood*(Diphoton.Minv > 230 && Photon1.pt>125 && Photon2.pt>125 && Photon1.isEB && Photon2.isEB)*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  // cuts["BE"] = "isGood*(Diphoton.Minv > 330 && Photon1.pt>125 && Photon2.pt>125 && ( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE )))*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";

  std::vector<std::string> samples = getSampleList();

  std::vector<int> stringScales = {3000, 3500, 4000, 4500, 5000, 5500, 6000};

  std::string coupling = "";
  std::string M_bins = "";
  double Ngen = 0.;

  std::map<std::string, std::vector<eff_reco> > the_eff_reco; 
  std::map<std::string, TH1F *> histograms;
  std::string histo_name;

  for(auto isample : samples) {
    //Run a single year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    //Run only on RS samples for now. 
    if( isample.find("RSGravitonToGammaGamma") == std::string::npos ) continue;

    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    std::cout << "xsec " << ExoDiPhotons::crossSection(getSampleBase(isample,year)) << std::endl;
    coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");

    int nBins = 2385;
    double xMin = 0.0;
    double xMax = 0.0; 
    if ( coupling == "001" ){xMax = 6000. ;}
    else if ( coupling == "01" ){xMax = 9000. ;}
    else if ( coupling == "02" ){xMax = 9000. ;}
    else {
      std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
      exit(1);
    }

    eff_reco tmpeffreco;

    std::string baseName(getSampleBase(isample, year));

    //========================================================================
    //Initial Distribution
    //========================================================================
 
    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName;
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, "no_selection");
    Ngen = histograms[histo_name]->Integral();
    // std::cout << "Ngen " << Ngen << std::endl;
    // std::cout << "Name " << getBase(isample) << " chains[getBase(isample)]->GetName() " << chains[getBase(isample)]->GetName() << " chains[getBase(isample)]->GetEntries() " << chains[getBase(isample)]->GetEntries() << std::endl;

    //========================================================================
    //Acceptance
    //========================================================================
    //------------------------------------------------------------------------
    //Sel: ABB
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "ABB";
    tmpeffreco.cat = "ABB";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_ABB";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["ABB"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   

    std::cout << "A BB " << tmpeffreco.efforacc << std::endl; 

    //------------------------------------------------------------------------
    //Sel: ABE
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "ABE";
    tmpeffreco.cat = "ABE";

    histo_name = baseName + "_ABE";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["ABE"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   
    std::cout << "A BE " << tmpeffreco.efforacc << std::endl; 
    //------------------------------------------------------------------------
    //Sel: ATotal
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "ATotal";
    tmpeffreco.cat = "ATotal";

    histo_name = baseName + "_ATotal";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["ATotal"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   
    std::cout << "A Total " << tmpeffreco.efforacc << std::endl; 

    std::cout << "-----------------------------------------------------------------" << std::endl;
    //========================================================================
    //Selection Efficiency 
    //========================================================================
    //------------------------------------------------------------------------
    //Sel: eBB
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "eBB";
    tmpeffreco.cat = "eBB";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_eBB";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["eBB"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   

    std::cout << "e BB " << tmpeffreco.efforacc << std::endl; 

    //------------------------------------------------------------------------
    //Sel: eBE
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "eBE";
    tmpeffreco.cat = "eBE";

    histo_name = baseName + "_eBE";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["eBE"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   
    std::cout << "e BE " << tmpeffreco.efforacc << std::endl; 
    //------------------------------------------------------------------------
    //Sel: eTotal
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "eTotal";
    tmpeffreco.cat = "eTotal";

    histo_name = baseName + "_eTotal";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["eTotal"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   
    std::cout << "e Total " << tmpeffreco.efforacc << std::endl; 

    //========================================================================
    //Efficiency x Acceptance
    //========================================================================
    //------------------------------------------------------------------------
    //Sel: exABB
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "exABB";
    tmpeffreco.cat = "exABB";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_exABB";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exABB"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * ExoDiPhotons::crossSection(getSampleBase(isample,year)) ) ;

    the_eff_reco[coupling].push_back(tmpeffreco);   

    std::cout << "exA BB " << tmpeffreco.efforacc << std::endl; 

    //------------------------------------------------------------------------
    //Sel: exABE
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "exABE";
    tmpeffreco.cat = "exABE";

    histo_name = baseName + "_exABE";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exABE"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * ExoDiPhotons::crossSection(getSampleBase(isample,year)) ) ;

    the_eff_reco[coupling].push_back(tmpeffreco);   
    std::cout << "exA BE " << tmpeffreco.efforacc << std::endl; 
    //------------------------------------------------------------------------
    //Sel: exATotal
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "exATotal";
    tmpeffreco.cat = "exATotal";

    histo_name = baseName + "_exATotal";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exATotal"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * ExoDiPhotons::crossSection(getSampleBase(isample,year)) ) ;

    the_eff_reco[coupling].push_back(tmpeffreco);   
    std::cout << "exA Total " << tmpeffreco.efforacc << std::endl; 

   
  }//end of loop through samples

  return the_eff_reco;

}

//-----------------------------------------------------------------------------------
TGraphErrors* graph(std::string coup, std::map<std::string, std::vector<eff_reco> > quant){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraphErrors* graph;

  unsigned int nP = quant[coup].size(); 

  double masses[nP];
  double massesErr[nP];
  double qu[nP];
  double quErr[nP];

  for(unsigned int iE =0; iE < nP; iE++){
    masses[iE] = std::stod(quant[coup][iE].M_bins); 
    massesErr[iE] = 0.;
    qu[iE] = quant[coup][iE].efforacc; 
    quErr[iE] = 0.;
  }
  
  graph = new TGraphErrors(nP, masses, qu, massesErr, quErr);

  return graph;

}

//-----------------------------------------------------------------------------------
TGraphErrors* getXsecGraph( std::vector<xsec> xsections)
{  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraphErrors* graph;

  unsigned int nP = xsections.size(); 

  double masses[nP];
  double massesErr[nP];
  double xsecval[nP];
  double xsecErr[nP];
  
  for(unsigned int iE =0; iE < nP; iE++){
    masses[iE] = std::stod(xsections[iE].M_bins); 
    massesErr[iE] = 0.;
    xsecval[iE] = xsections[iE].val; 
    xsecErr[iE] = xsections[iE].error;
  }
  
  graph = new TGraphErrors(nP, masses, xsecval, massesErr, xsecErr);

  return graph;

}
//-----------------------------------------------------------------------------------
void plotsplines( TCanvas* cs, RooSpline1D * Spline, RooRealVar* MH, double xmin, double upperxmax, std::string label){

  cs->cd();
  TGraph * graph = new TGraph();
  int point=0;
  
  for (double m =xmin; m<upperxmax; m=m+5.){
    MH->setVal(m);
    graph->SetPoint(point,m,Spline->getVal());
    point++;
  }

  TPaveText *pt2 = new TPaveText(.6,.85,.9,1.0,"NDC");
  pt2->AddText(Form("%s XS Spline",label.c_str()));
  graph->Draw("ALP");
  pt2->Draw();

}
//-----------------------------------------------------------------------------------
RooSpline1D* graphToSpline(std::string name, TGraphErrors *graph, RooRealVar* MH, double xmin, double upperxmax){
  
  std::vector<double> xValues, yValues;
  for (double mh=xmin; mh<(upperxmax+0.25); mh+=5.){
    xValues.push_back(mh);
    yValues.push_back(graph->Eval(mh));
  }
  RooSpline1D *res = new RooSpline1D(name.c_str(),name.c_str(),*MH,xValues.size(),&(xValues[0]),&(yValues[0]));
  return res;
}

//-----------------------------------------------------------------------------------
std::vector<theGraphs> plot(TCanvas* cc, std::string coup, std::map<std::string, std::vector<eff_reco> > quantBB, std::map<std::string, std::vector<eff_reco> > quantBE , std::map<std::string, std::vector<eff_reco> > quantTotal, bool doave, std::string label)
{
  
  std::vector<theGraphs> graphs; //BB, BE, Total
  theGraphs tmpgr; 
  
  std::map<std::string, double> upperxmax;
  upperxmax["001"] = 6000.;
  upperxmax["01"]  = 9000.;
  upperxmax["02"]  = 9000.;

  cc->cd();

  if (!doave){
    TLegend* legmc = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);

    if (coup == "001") {legmc->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-4}","C");}
    else if (coup == "01") {legmc->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-2}","C");}
    else if (coup == "02") {legmc->SetHeader("#frac{#Gamma}{m} = 5.6 #times 10^{-2}","C");}

    //BB
    tmpgr.coup = coup;
    tmpgr.cat = "BB";
    tmpgr.thegraph = graph(coup, quantBB);
    tmpgr.thegraph->GetYaxis()->SetRangeUser(0., 1.0);
    tmpgr.thegraph->SetMarkerColor(2);
    tmpgr.thegraph->SetLineColor(2);
    legmc->AddEntry( tmpgr.thegraph , "EBEB J=2" ,"pl" );
    
    tmpgr.fun = new TF1(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "pol2", 500,upperxmax[coup]);
    tmpgr.fun->SetLineColor(2);
    tmpgr.thegraph->Fit(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "R");
    tmpgr.thegraph->GetYaxis()->SetTitle(Form("%s",label.c_str()));
    tmpgr.thegraph->GetXaxis()->SetTitle("m_{X} [GeV]");
    tmpgr.thegraph->Draw("APE");
    tmpgr.fun->Draw("same");

    graphs.push_back(tmpgr);

    tmpgr.coup = coup;
    tmpgr.cat = "BE";
    tmpgr.thegraph = graph(coup, quantBE);
    tmpgr.thegraph->GetYaxis()->SetRangeUser(0., 1.0);
    tmpgr.thegraph->SetMarkerColor(4);
    tmpgr.thegraph->SetLineColor(4);
    legmc->AddEntry( tmpgr.thegraph , "EBEE J=2" ,"pl" );
    
    tmpgr.fun = new TF1(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "pol2", 500,upperxmax[coup]);
    tmpgr.fun->SetLineColor(4);
    tmpgr.thegraph->Fit(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "R");
    tmpgr.thegraph->GetYaxis()->SetTitle(Form("%s",label.c_str()));
    tmpgr.thegraph->GetXaxis()->SetTitle("m_{X} [GeV]");
    tmpgr.thegraph->Draw("PSE");
    tmpgr.fun->Draw("same");

    graphs.push_back(tmpgr);

    tmpgr.coup = coup;
    tmpgr.cat = "Total";
    tmpgr.thegraph = graph(coup, quantTotal);
    tmpgr.thegraph->GetYaxis()->SetRangeUser(0., 1.0);
    tmpgr.thegraph->SetMarkerColor(8);
    tmpgr.thegraph->SetLineColor(8);
    legmc->AddEntry( tmpgr.thegraph , "Total J=2" ,"pl" );
    
    tmpgr.fun = new TF1(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "pol2", 500,upperxmax[coup]);
    tmpgr.fun->SetLineColor(8);
    tmpgr.thegraph->Fit(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "R");
    tmpgr.thegraph->GetYaxis()->SetTitle(Form("%s",label.c_str()));
    tmpgr.thegraph->GetXaxis()->SetTitle("m_{X} [GeV]");
    tmpgr.thegraph->Draw("PSE");
    tmpgr.fun->Draw("same");

    graphs.push_back(tmpgr);

    legmc->Draw("same");
  }

  return graphs;

}

// //-----------------------------------------------------------------------------------
// void plot(TCanvas* cc, std::map<std::string, std::vector<eff_reco> > effreco, bool doave){
  
//   gStyle->SetOptStat(0);
//   gStyle->SetOptFit(0);

//   double masses001[effreco["001"].size()];
//   double massesErr001[effreco["001"].size()];
//   double eff001[effreco["001"].size()];
//   double effErr001[effreco["001"].size()];

//   double masses01[effreco["01"].size()];
//   double massesErr01[effreco["01"].size()];
//   double eff01[effreco["01"].size()];
//   double effErr01[effreco["01"].size()];

//   double masses02[effreco["02"].size()];
//   double massesErr02[effreco["02"].size()];
//   double eff02[effreco["02"].size()];
//   double effErr02[effreco["02"].size()];

//   //Since the final total selection efficiency is the average of the 3 couplings 
//   //we will save that also. 
//   //BE CAREFUL: For the moment we will only use points where we have values for all the 
//   //three couplings
//   std::map<std::string, std::vector<eff_reco> > effrecoBB;
//   std::map<std::string, std::vector<eff_reco> > effrecoBE;
//   std::map<std::string, std::vector<eff_reco> > effrecoTotal;
//   if (doave){
//     for(unsigned int iE =0; iE < effreco["01"].size(); iE++){
//       std::cout << effreco["01"][iE].M_bins << std::endl;
//       double curmass = std::stod(effreco["01"][iE].M_bins); 
//       if ( curmass > 5000. || curmass == 2750. || curmass == 3250. || curmass == 4250. || curmass == 4500. || curmass == 4750.  ){continue;}
//       if ( effreco["01"][iE].sel == "ABB"){effrecoBB["01"].push_back( effreco["01"][iE] );}
//       if ( effreco["01"][iE].sel == "ABE"){effrecoBE["01"].push_back( effreco["01"][iE] );}
//       if ( effreco["01"][iE].sel == "ATotal"){effrecoTotal["01"].push_back( effreco["01"][iE] );}
//     }

//     for(unsigned int iE =0; iE < effreco["02"].size(); iE++){
//       double curmass = std::stod(effreco["02"][iE].M_bins); 
//       if ( curmass > 5000. || curmass == 2750. || curmass == 3250. || curmass == 4250. || curmass == 4500. || curmass == 4750.  ){continue;}
//       if ( effreco["02"][iE].sel == "ABB"){effrecoBB["02"].push_back( effreco["02"][iE] );}
//       if ( effreco["02"][iE].sel == "ABE"){effrecoBE["02"].push_back( effreco["02"][iE] );}
//       if ( effreco["02"][iE].sel == "ATotal"){effrecoTotal["02"].push_back( effreco["02"][iE] );}
//     }

//     for(unsigned int iE =0; iE < effreco["001"].size(); iE++){
//       if ( effreco["001"][iE].sel == "BB"){effrecoBB["001"].push_back( effreco["001"][iE] );}
//       if ( effreco["001"][iE].sel == "BE"){effrecoBE["001"].push_back( effreco["001"][iE] );}
//     }


//   }
  
//   double massesBB[effrecoBB["01"].size()];
//   double massesErrBB[effrecoBB["01"].size()];
//   double effBB[effrecoBB["01"].size()];
//   double effErrBB[effrecoBB["01"].size()];

//   double massesBE[effrecoBE["01"].size()];
//   double massesErrBE[effrecoBE["01"].size()];
//   double effBE[effrecoBE["01"].size()];
//   double effErrBE[effrecoBE["01"].size()];

//   std::map<std::string, double> upperxmax;
//   upperxmax["001"] = 6000.;
//   upperxmax["01"]  = 9000.;
//   upperxmax["02"]  = 9000.;
  
//   for(unsigned int iE =0; iE < effreco["001"].size(); iE++){
//     masses001[iE] = std::stod(effreco["001"][iE].M_bins); 
//     massesErr001[iE] = 0.;
//     eff001[iE] = effreco["001"][iE].eff; 
//     effErr001[iE] = 0.;
//   }

//   for(unsigned int iE =0; iE < effreco["01"].size(); iE++){
//     masses01[iE] = std::stod(effreco["01"][iE].M_bins); 
//     massesErr01[iE] = 0.;
//     eff01[iE] = effreco["01"][iE].eff; 
//     effErr01[iE] = 0.;
//   }
  
//   for(unsigned int iE =0; iE < effreco["02"].size(); iE++){
//     masses02[iE] = std::stod(effreco["02"][iE].M_bins); 
//     massesErr02[iE] = 0.;
//     eff02[iE] = effreco["02"][iE].eff; 
//     effErr02[iE] = 0.;
//   }

//   if (doave){
//     //BB average
//     for(unsigned int iE =0; iE < effrecoBB["01"].size(); iE++){
//       std::cout << effrecoBB["01"].size() << " " << effrecoBB["001"].size() <<  " " << effrecoBB["02"].size() << std::endl;
//       std::cout << effrecoBB["01"][iE].M_bins << " " << effrecoBB["02"][iE].M_bins << " " << effrecoBB["001"][iE].M_bins << std::endl;
//       massesBB[iE] = std::stod(effrecoBB["01"][iE].M_bins); 
//       massesErrBB[iE] = 0.;
//       //The tweaks to divide by 3 and use all couplings are above
//       effBB[iE] = (effrecoBB["001"][iE].eff+effrecoBB["01"][iE].eff+effrecoBB["02"][iE].eff)/3.; 
//       effErrBB[iE] = 0.;
//     }
//     //BE average
//     for(unsigned int iE =0; iE < effrecoBE["01"].size(); iE++){
//       massesBE[iE] = std::stod(effrecoBE["01"][iE].M_bins); 
//       massesErrBE[iE] = 0.;
//       effBE[iE] = (effrecoBE["001"][iE].eff+effrecoBE["01"][iE].eff+effrecoBE["02"][iE].eff)/3.; 
//       effErrBE[iE] = 0.;
//     }
//   }
  
  
//   cc->cd();

//   if (!doave){
//     TLegend* legmc = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
//     legmc->SetTextFont(42);
//     legmc->SetBorderSize(0);
//     legmc->SetFillStyle(0);

//     TGraphErrors* graph001 = new TGraphErrors(effreco["001"].size(), masses001, eff001, massesErr001, effErr001);
//     graph001->GetYaxis()->SetRangeUser(0., 1.0);
//     graph001->SetMarkerColor(2);
//     graph001->SetLineColor(2);
//     legmc->AddEntry( graph001 , "#frac{#Gamma}{m} = 1.4 #times 10^{-4}" ,"pl" );

//     TGraphErrors* graph01 = new TGraphErrors(effreco["01"].size(), masses01, eff01, massesErr01, effErr01);
//     graph01->GetYaxis()->SetRangeUser(0., 1.0);
//     graph01->SetMarkerColor(4);
//     graph01->SetLineColor(4);
//     legmc->AddEntry( graph01 , "#frac{#Gamma}{m} = 1.4 #times 10^{-2}" ,"pl" );

//     TGraphErrors* graph02 = new TGraphErrors(effreco["02"].size(), masses02, eff02, massesErr02, effErr02);
//     graph02->GetYaxis()->SetRangeUser(0., 1.0);
//     graph02->SetMarkerColor(8);
//     graph02->SetLineColor(8);
//     legmc->AddEntry( graph02 , "#frac{#Gamma}{m} = 5.6 #times 10^{-2}" ,"pl" );

 
//     TF1 *f001 = new TF1("f001", "pol2", 500,upperxmax["001"]);
//     f001->SetLineColor(2);
//     graph001->Fit("f001", "R");
//     graph001->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
//     graph001->GetXaxis()->SetTitle("m_{G} [GeV]");
//     graph001->Draw("APE");
//     f001->Draw("same");

//     TF1 *f01 = new TF1("f01", "pol2", 500,upperxmax["01"]);
//     f01->SetLineColor(4);
//     graph01->Fit("f01", "R");
//     graph01->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
//     graph01->GetXaxis()->SetTitle("m_{G} [GeV]");
//     graph01->Draw("PSE");
//     f01->Draw("same");

//     TF1 *f02 = new TF1("f02", "pol2", 500,upperxmax["02"]);
//     f02->SetLineColor(8);
//     graph02->Fit("f02", "R");
//     graph02->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
//     graph02->GetXaxis()->SetTitle("m_{G} [GeV]");
//     graph02->Draw("PSE");
//     f02->Draw("same");

//     legmc->Draw("same");
//   } else{
//     TLegend* legmc = new TLegend(0.58, 0.70, 0.75, 0.9, "", "bNDC");
//     legmc->SetTextFont(42);
//     legmc->SetBorderSize(0);
//     legmc->SetFillStyle(0);

//     TGraphErrors* graphBB = new TGraphErrors(effrecoBB["01"].size(), massesBB, effBB, massesErrBB, effErrBB);
//     graphBB->GetYaxis()->SetRangeUser(0., 1.0);
//     graphBB->SetMarkerColor(2);
//     graphBB->SetLineColor(2);
//     legmc->AddEntry( graphBB , "BB" ,"pl" );

//     TGraphErrors* graphBE = new TGraphErrors(effrecoBE["01"].size(), massesBE, effBE, massesErrBE, effErrBE);
//     graphBE->GetYaxis()->SetRangeUser(0., 1.0);
//     graphBE->SetMarkerColor(4);
//     graphBE->SetLineColor(4);
//     legmc->AddEntry( graphBE , "BE" ,"pl" );

//     TF1 *fBB = new TF1("fBB", "pol2", 500,upperxmax["01"]);
//     fBB->SetLineColor(2);
//     graphBB->Fit("fBB", "R");
//     graphBB->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
//     graphBB->GetXaxis()->SetTitle("m_{G} [GeV]");
//     graphBB->Draw("APE");
//     fBB->Draw("same");

//     TF1 *fBE = new TF1("fBE", "pol2", 500,upperxmax["01"]);
//     fBE->SetLineColor(4);
//     graphBE->Fit("fBE", "R");
//     graphBE->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
//     graphBE->GetXaxis()->SetTitle("m_{G} [GeV]");
//     graphBE->Draw("PSE");
//     fBE->Draw("same");

//     legmc->Draw("same");

//   }
  




// }
//-----------------------------------------------------------------------------------
std::map<std::string, std::vector<eff_reco> > reduce(std::map<std::string, std::vector<eff_reco> > effreco, std::string cut){

  std::map<std::string, std::vector<eff_reco> > out; 
  
  for(unsigned int iE =0; iE < effreco["001"].size(); iE++){
    if (effreco["001"][iE].sel == cut){ out["001"].push_back(effreco["001"][iE]); }
  }

  for(unsigned int iE =0; iE < effreco["01"].size(); iE++){
    if (effreco["01"][iE].sel == cut){ out["01"].push_back(effreco["01"][iE]); }
  }

  for(unsigned int iE =0; iE < effreco["02"].size(); iE++){
    if (effreco["02"][iE].sel == cut){ out["02"].push_back(effreco["02"][iE]); }
  }

  return out; 
}

//-----------------------------------------------------------------------------------
TH1F * createHisto(TChain * newtree1, const std::string &histo_name, int nBins, double xMin, double xMax, std::string cut){
   
  // std::cout << "Making histograms with cut\n" << cut << std::endl;
  TH1F * hist = new TH1F(histo_name.c_str(), histo_name.c_str(), nBins, xMin, xMax);

  if (cut == "no_selection"){
    newtree1->Project(histo_name.c_str(), "Diphoton.Minv");
  } else {
    newtree1->Project(histo_name.c_str(), "Diphoton.Minv", cut.c_str());
  }
    // move overflow to last bin
  int nbins = hist->GetNbinsX();
  float lastBin = hist->GetBinContent(nbins);
  float overflow = hist->GetBinContent(nbins+1);
  hist->SetBinContent(nbins, lastBin + overflow);
  hist->SetBinContent(nbins+1, 0.0);

  return hist;

}

//-----------------------------------------------------------------------------------
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim)
{
  unsigned first_delim_pos = s.find(start_delim);
  unsigned end_pos_of_first_delim = first_delim_pos + start_delim.length();
  unsigned last_delim_pos = s.find(stop_delim);

  return s.substr(end_pos_of_first_delim,
		  last_delim_pos - end_pos_of_first_delim);
}


//-----------------------------------------------------------------------------------
// remove year
std::string getSampleBase(const std::string & sampleName, const std::string & year)
{
  std::string newString(sampleName);
  if( sampleName.find("_201") != std::string::npos) {
    newString.replace(newString.find("_201"), 5, "");
  }
  if(sampleName.find("_R2F2") != std::string::npos) {
    newString.replace(newString.find("_R2F2"), 5, "_diphotonkfactorScalesUp");
  }
  if(sampleName.find("_R0p5F0p5") != std::string::npos) {
    newString.replace(newString.find("_R0p5F0p5"), 9, "_diphotonkfactorScalesDown");
  }
  // "data_obs" is always the name of the data observation histogram
  std::string data("data_" + year);
  if( sampleName.compare(data) == 0) newString = "data_obs";
  return newString;
}

//-----------------------------------------------------------------------------------
// ignore variations to get dataset name
std::string getBase(const std::string & sampleName)
{
  if(sampleName.compare("gg_R2F2_2016") == 0 ) return "gg_2016";
  if(sampleName.compare("gg_R0p5F0p5_2016") == 0 ) return "gg_2016";
  if(sampleName.compare("gg_R2F2_2017") == 0 ) return "gg_2017";
  if(sampleName.compare("gg_R0p5F0p5_2017") == 0 ) return "gg_2017";
  if(sampleName.compare("gg_R2F2_2018") == 0 ) return "gg_2018";
  if(sampleName.compare("gg_R0p5F0p5_2018") == 0 ) return "gg_2018";
  return sampleName;
}

void writeToJson(std::map< std::string , std::vector<double> > valuestowrite, std::string outputfile)
{
  // Short alias for this namespace
  namespace pt = boost::property_tree;
  
  // Create a root
  pt::ptree root;

  for (auto &valmap : valuestowrite){
    pt::ptree values_node;
    for (auto &val : valmap.second){
      // Create an unnamed node containing the value
      pt::ptree val_node;
      val_node.put("", val);

      // Add this node to the list.
      values_node.push_back(std::make_pair("", val_node));
    }
    root.add_child(valmap.first, values_node);
  }
  pt::write_json(std::cout, root);
  pt::write_json(outputfile, root);


}

std::map<std::string , std::vector<xsec> > loadXsections(const std::string & year)
{
  
  std::map<std::string, std::vector<xsec> >thexsections; 

  std::vector<std::string> samples = getSampleList();

  xsec tmpxsec; 
  std::string coup;

  for(auto isample : samples) {
    //Run a single year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    //Run only on RS samples for now. 
    if( isample.find("RSGravitonToGammaGamma") == std::string::npos ) continue;

    tmpxsec.name = getSampleBase(isample,year);
    tmpxsec.val = ExoDiPhotons::crossSection(getSampleBase(isample,year));
    tmpxsec.error = 0.;

    std::cout << "xsec " << ExoDiPhotons::crossSection(getSampleBase(isample,year)) << std::endl;
    coup = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    tmpxsec.M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");

    thexsections[coup].push_back(tmpxsec);

  }
  
  return thexsections;

}
