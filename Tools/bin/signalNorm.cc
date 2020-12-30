#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "diphoton-analysis/CommonClasses/interface/CrossSections.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"

//ROOT
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPaveText.h"

//RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooVoigtian.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooPolyVar.h"
#include "RooCustomizer.h"
#include "RooBinning.h"

int nBinsMass = 1100;//150//1250

//-----------------------------------------------------------------------------------
struct eff_reco {

  std::string samplename;
  std::string kMpl;
  std::string M_bins;
  //Values for sel are: no_selection, isGood
  std::string sel;
  //Values for cat are both, BB, BE.
  std::string cat;
  //Efficiency or acceptance or efficiency times acceptance
  //depending on the the sel value
  double efforacc;

};

//-----------------------------------------------------------------------------------
struct theGraphs{

  std::string coup;
  std::string cat;
  TGraphErrors* thegraph;
  TF1 * fun;
  double funp0;
  double funp1;
  double funp2;
  double funp3;

};

//-----------------------------------------------------------------------------------
struct xsec{

  std::string name;
  double coup;
  std::string M_bins;
  double val;
  double error ;

};

//-----------------------------------------------------------------------------------
struct pdfs{

  std::string name;
  std::string coup;
  std::string cat;
  RooAbsPdf * pdf;

};

using namespace RooFit;

//-----------------------------------------------------------------------------------
//Declarations here definition after main
std::map<std::string, std::vector<eff_reco> > computefficiency(const std::string &year);
void plotefficiency(TCanvas* cc, std::map<std::string, std::vector<eff_reco> > coup_eff_reco, bool doave);
TGraphErrors* graph(std::string coup, std::map<std::string, std::vector<eff_reco> > quant);
void plotexAs(std::map<std::string, RooProduct *> exAs, RooRealVar* MH, RooRealVar* kmpl, std::vector<std::string> cats, std::vector<std::string> coups, std::string insigname, std::string outplots);
void plotexAs(std::map<std::string, std::map<std::string, RooProduct *> > exAs, RooRealVar* MH, RooRealVar* kmpl, std::vector<std::string> cats, std::vector<std::string> coups, std::string insigname, std::string outplots);
std::vector<theGraphs> plot(TCanvas* cc, std::string coup, std::map<std::string, std::vector<eff_reco> > quantBB, std::map<std::string, std::vector<eff_reco> > quantBE , std::map<std::string, std::vector<eff_reco> > quantTotal, bool doave, std::string label, std::string insigname);
std::map<std::string, std::vector<eff_reco> > reduce(std::map<std::string, std::vector<eff_reco> > coup_eff_reco, std::vector<std::string>coups, std::string cut);
void plotallcoups(TCanvas* ccallcoup, std::vector<theGraphs> graphsofeff001, std::vector<theGraphs> graphsofeff01, std::vector<theGraphs> graphsofeff02, std::string insigname);
TH1F * createHisto(TChain * newtree1, const std::string &histo_name, int nBins, double xMin, double xMax, std::string cut);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);
// void writeToJson(std::map< std::string , std::vector<double> > valuestowrite, std::string outputfile);
void writeToJson( std::map<std::string, std::vector<eff_reco> > valuestowrite, std::string outputfile);
std::map<std::string, std::vector<eff_reco> > readJson(std::string inputfile);
std::map<std::string, std::vector<xsec> > loadXsections(const std::string &year, bool basedonCoup);
TGraphErrors* getXsecGraph( std::vector<xsec> xsections, bool basedonCoup);
RooSpline1D* graphToSpline(std::string name, TGraphErrors *graph, RooRealVar* MH, double xmin, double upperxmax);
void plotsplines(TCanvas* cc,  RooSpline1D * xsSplines, RooRealVar* MH, double xmin, double upperxmax, std::string label);
std::vector<eff_reco> ave_eff(std::map<std::string, std::vector<eff_reco> > effreco, std::string insigname);
TGraphErrors* avegraph(std::vector<eff_reco> quant);
TGraphErrors * getgraph(std::vector<double> xvalues, std::vector<double> yvalues, std::vector<double> xvaluesErr, std::vector<double> yvaluesErr);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

  std::string region, year, inputworkspaces, outputworkspaces, outplots, insigname, Json;

  if(argc!=7) {
    std::cout << "Syntax: signalNorm.exe [2016/2017/2018] [Json] [inputworkspaces] [outputworkspaces] [outplots] [insigname]" << std::endl;
    return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    Json = argv[2];
    if (Json !="createJson" and Json!="readJson"){
      std::cout << "Only createJson and readJson are allowed " << std::endl;
    }
    inputworkspaces = argv[3];
    outputworkspaces = argv[4];
    outplots = argv[5];
    insigname = argv[6];
  }

  //========================================================================
  // include signal samples but not unskimmed data samples
  init(false, true);

  //========================================================================
  //========================================================================
  //ZERO PART : We give whatever input/output variables needed here.
  //========================================================================
  //========================================================================
  //========================================================================

  //There are 3 files, one for each coupling, containing the signal model.
  std::map<std::string , TFile*> fin;
  std::map<std::string , RooWorkspace*> wsin;

  //This is where we will save the output workspace
  std::map<std::string , TFile*> fout;
  std::map<std::string , RooWorkspace*> ws_out;

  //Couplings for the moment
  std::vector<std::string> coups; coups.clear();
  if (insigname == "grav"){
    coups.push_back("kMpl001");
    coups.push_back("kMpl01");
    coups.push_back("kMpl02");
  } else if (insigname == "heavyhiggs"){
    coups.push_back("0p014");
    coups.push_back("1p4");
    coups.push_back("5p6");
  } else {
    std::cout << "The signal you are looking for doesn't exist. " << std::endl;
    exit(1);
  }

  //Categories
  std::vector<std::string> cats; cats.clear();
  cats.push_back("EBEB");
  cats.push_back("EBEE");
  cats.push_back("All");

  // for (auto cp : coups){
    //Input
    // fin[cp] = TFile::Open( Form("%s/SignalParametricShapes_ws_%s.root", inputworkspaces.c_str(), cp.c_str()) );
    // wsin[cp] = (RooWorkspace*) fin[cp]->Get("ws_inputs");
    // std::cout<< wsin[cp]->var("mgg")->getBinning().numBins()<< std::endl;
    //Output
  //   fout[cp] = new TFile(Form("%s/%s_%s_%s.root", outputworkspaces.c_str(), insigname.c_str(),cp.c_str(),year.c_str()), "RECREATE");
  //   ws_out[cp] = new RooWorkspace("wtemplates","wtemplates");
  // }

  //Here we take note of the pdfs we are going to read
  std::vector<pdfs> thepdfs;

  pdfs tmppdf;
  for (auto cp : coups){
    for (auto cat : cats){

      tmppdf.name = Form("SignalShape_%s_%s", cp.c_str(), cat.c_str() );
      tmppdf.cat = cat;
      tmppdf.coup = cp;

      thepdfs.push_back(tmppdf);

    }
  }

  //Parameters in the input workspace.
  // std::map<std::string , std::vector<std::string> > params; //[cat][params]
  // std::map<std::string , std::string> reparam_by_cat;
  // for (auto cat : cats){
  //   params[cat].push_back( Form("thetaSmear%s",cat.c_str() ) );
  //   params[cat].push_back( "deltaSmear" );
  //   reparam_by_cat[ Form("thetaSmear%s",cat.c_str() ) ] = Form("thetaSmear%s_13TeV",cat.c_str() );
  //   reparam_by_cat[ "deltaSmear" ] = Form("deltaSmear%s",cat.c_str() );
  // }

  //Observables
  // std::map<std::string , RooRealVar * > obsCat;
  // obsCat["EBEE"] =  new RooRealVar("mggEBEE","mggEBEE", 2335, 330., 5000.);
  // obsCat["EBEB"] =  new RooRealVar("mggEBEB","mggEBEB", 2385, 230., 5000.);
  // obsCat["EBEE"] =  new RooRealVar("mgg","mgg", nBinsMass, 500., 6000.);//2250
  // // RooBinning* mggRooBinning = new RooBinning(nBinsMass, 500., 6000.,"mgg");
  // obsCat["EBEE"]->setBinning(*mggRooBinning,"mgg");

  // obsCat["EBEB"] =  obsCat["EBEE"];
  // obsCat["All"]  =  new RooRealVar("mggAll","mggAll", nBinsMass, 230., 6000.);//2385
  // // RooAbsBinning& mggBinsEBEE =  obsCat["EBEE"]->getBinning();
  // // obsCat["EBEE"]->setBinning(mggBinsEBEE);
  // // RooAbsBinning& mggBinsEBEB =  obsCat["EBEB"]->getBinning();
  // // obsCat["EBEB"]->setBinning(mggBinsEBEB);
  // obsCat["EBEB"]->setBinning(*mggRooBinning,"mgg");

  // std::cout<< "++++++++++++++++++++++++++++++++"<< std::endl;
  // std::cout<< obsCat["EBEE"]->getBinning().numBins()<< std::endl;
  // std::cout<< obsCat["EBEB"]->getBinning().numBins()<< std::endl;
  // std::cout<< obsCat["EBEB"]->hasBinning("mgg")<< std::endl;
  // std::cout<< "++++++++++++++++++++++++++++++++"<< std::endl;
  RooRealVar * obsCat =  new RooRealVar("mgg","mgg", nBinsMass, 500., 6000.);;
  obsCat->setBins(nBinsMass);

  // RooBinning* mggRooBinning = new RooBinning(nBinsMass, 500., 6000.,"mgg");
  // obsCat->setBinning(*mggRooBinning);
  std::cout<< "++++++++++++++++++++++++++++++++"<< std::endl;
  // std::cout<< obsCat->getBinning().numBins()<< std::endl;
  std::cout<< obsCat->getBins()<< std::endl;
  std::cout<< "++++++++++++++++++++++++++++++++"<< std::endl;

  //This is the uncertainty used on the scale <<<==== CHECK THIS
  // double unc = 0.01;
  // std::map<std::string , RooRealVar * > thetaScale;
  // std::map<std::string , RooRealVar * > deltaScale;
  // std::map<std::string , RooArgList * > rooList;
  // std::map<std::string , RooFormulaVar * > scaledMean;
  //========================================================================
  //========================================================================
  //FIRST PART: Create or read json file with efficiency (e), acceptance (A)
  //            and efficiency times acceptance (exA).
  //========================================================================
  //========================================================================

  //========================================================================
  std::map<std::string, std::vector<eff_reco> > effreco;
  //========================================================================
  //Will save the values in a json file for faster running.
  if (Json == "createJson"){
    effreco = computefficiency(year);
    writeToJson(effreco, Form("Tools/json/acceptance_%s.json", year.c_str()));
    return 0;
  } else if (Json == "readJson"){
    effreco = readJson(Form("Tools/json/acceptance_%s.json", year.c_str()));
  }

  //========================================================================
  //========================================================================
  //SECOND PART: Plot e, A, average e and save the graphs for later.
  //========================================================================
  //========================================================================

  //acceptance per coupling and mass point
  std::map<std::string, std::vector<eff_reco> > accBB = reduce(effreco, coups, "ABB");
  std::map<std::string, std::vector<eff_reco> > accBE = reduce(effreco, coups, "ABE");
  std::map<std::string, std::vector<eff_reco> > accTotal = reduce(effreco, coups, "ATotal");
  //efficiency per coupling and mass point
  std::map<std::string,std::vector<eff_reco> > effBB = reduce(effreco, coups, "eBB");
  std::map<std::string, std::vector<eff_reco> > effBE = reduce(effreco, coups, "eBE");
  std::map<std::string, std::vector<eff_reco> > effTotal = reduce(effreco, coups, "eTotal");
  //acceptance x efficiency per coupling and mass point
  std::map<std::string,std::vector<eff_reco> > accxeffBB = reduce(effreco, coups, "exABB");
  std::map<std::string, std::vector<eff_reco> > accxeffBE = reduce(effreco, coups, "exABE");
  std::map<std::string, std::vector<eff_reco> > accxeffTotal = reduce(effreco, coups, "exATotal");

  //========================================================================
  //Eff plots and graphs
  // TCanvas* cckMpl001 = new TCanvas("cckMpl001", "cckMpl001");
  std::map<std::string, TCanvas* > cck;
  std::map<std::string, std::vector<theGraphs> > graphsofeff;

  for (auto cp : coups){

    cck[cp] = new TCanvas(Form("cc%s",cp.c_str()), Form("cc%s",cp.c_str()));
    graphsofeff[cp] = plot(cck[cp], cp, effBB, effBE, effTotal, false, "(#varepsilon #otimes A)/A", insigname);
    cck[cp]->SaveAs( Form("%s/effsVsMass_%s.png", outplots.c_str(),cp.c_str()) );

  }


  // graphsofeff["kMpl001"] = plot(cckMpl001, "kMpl001", effBB, effBE, effTotal, false, "(#varepsilon #otimes A)/A");
  // cckMpl001->SaveAs( Form("%s", outplots.c_str())    "/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effsVsMass_kMpl001.png");

  // TCanvas* cckMpl01 = new TCanvas("cckMpl01", "cckMpl01");
  // graphsofeff["kMpl01"] = plot(cckMpl01, "kMpl01", effBB, effBE, effTotal, false, "(#varepsilon #otimes A)/A");
  // cckMpl01->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effsVsMass_kMpl01.png");

  // TCanvas* cckMpl02 = new TCanvas("cckMpl02", "cckMpl02");
  // graphsofeff["kMpl02"] = plot(cckMpl02, "kMpl02", effBB, effBE, effTotal, false, "(#varepsilon #otimes A)/A");
  // cckMpl02->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effsVsMass_kMpl02.png");

  //Now, since we will average the eff for the 3 coupling, let's plot all coupling together
  //in the same plot to see if this is a logical think to do.
  TCanvas* ccallcoup = new TCanvas("ccallcoup", "ccallcoup");
  if (insigname == "grav"){
    plotallcoups(ccallcoup, graphsofeff["kMpl001"], graphsofeff["kMpl01"], graphsofeff["kMpl02"], insigname);
  } else if (insigname == "heavyhiggs"){
    plotallcoups(ccallcoup, graphsofeff["0p014"], graphsofeff["1p4"], graphsofeff["5p6"], insigname);
  }
  ccallcoup->SaveAs( Form("%s/effsVsMass_allcoup_%s.png", outplots.c_str(),insigname.c_str()) );

  //========================================================================
  //Acceptance plots and graphs
  // TCanvas* ccacckMpl001 = new TCanvas("ccacckMpl001", "ccacckMpl001");
  std::map<std::string, TCanvas* > ccacck;
  std::map<std::string, std::vector<theGraphs> > graphsofacc;

  for (auto cp : coups){

    ccacck[cp] = new TCanvas(Form("ccacc%s",cp.c_str()), Form("ccacc%s",cp.c_str()));
    graphsofacc[cp] = plot(ccacck[cp], cp, accBB, accBE, accTotal, false, "A", insigname);
    ccacck[cp]->SaveAs( Form("%s/accsVsMass_%s.png", outplots.c_str(),cp.c_str()) );

  }



  // graphsofacc["kMpl001"] = plot(ccacckMpl001, "kMpl001", accBB, accBE, accTotal, false, "A");
  // ccacckMpl001->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/accsVsMass_kMpl001.png");

  // TCanvas* ccacckMpl01 = new TCanvas("ccacckMpl01", "ccacckMpl01");
  // graphsofacc["kMpl01"] = plot(ccacckMpl01, "kMpl01", accBB, accBE, accTotal, false, "A");
  // ccacckMpl01->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/accsVsMass_kMpl01.png");

  // TCanvas* ccacckMpl02 = new TCanvas("ccacckMpl02", "ccacckMpl02");
  // graphsofacc["kMpl02"] = plot(ccacckMpl02, "kMpl02", accBB, accBE, accTotal, false, "A");
  // ccacckMpl02->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/accsVsMass_kMpl02.png");

  //========================================================================
  //exA plots and graphs
  std::map<std::string, TCanvas* > ccexAk;
  std::map<std::string, std::vector<theGraphs> > graphsofexAk;

  for (auto cp : coups){

    ccexAk[cp] = new TCanvas(Form("ccexAk%s",cp.c_str()), Form("ccexAk%s",cp.c_str()));
    graphsofexAk[cp] = plot(ccexAk[cp], cp, accxeffBB, accxeffBE, accxeffTotal, false, "#varepsilon #otimes A", insigname);
    ccexAk[cp]->SaveAs( Form("%s/accxeffsVsMass_%s.png", outplots.c_str(),cp.c_str()) );

  }


  //========================================================================
  //Average of the efficiency over all couplings
  std::vector<eff_reco> effBB_ave = ave_eff(effBB, insigname);
  std::vector<eff_reco> effBE_ave = ave_eff(effBE, insigname);
  std::vector<eff_reco> effTotal_ave = ave_eff(effTotal, insigname);

  //Get the average graph to fit
  std::map<std::string , TGraphErrors* > ave_graphs;
  ave_graphs["EBEB"] = avegraph(effBB_ave);
  ave_graphs["EBEE"] = avegraph(effBE_ave);
  ave_graphs["All"] = avegraph(effTotal_ave);

  //========================================================================
  //========================================================================
  //THIRD PART: Build the Eff(mX) RooPolyVar object
  //========================================================================
  //========================================================================

  //Ready to create the RooPolyVar for the different categories
  std::map<std::string, RooPolyVar *> effMX; //[cat][ RooPolyVar *]
  std::map<std::string , RooSpline1D *> effMXSplines;//[cat][RooSpline1D  *]
  std::map<std::string, TF1 *> ff;
  RooArgList *eff_coef;

  RooRealVar* MH = new RooRealVar("MH", "MH", 1000.);
  MH->setConstant();
  RooRealVar* kmpl = new RooRealVar("kmpl", "kmpl", 0.01);
  kmpl->setConstant();

  std::map<std::string , TCanvas*> can;
  double upperxmax = 9000.;
  double xmin = 230.;

  for (auto ct : cats){

    std::string mname = "ave_eff_" + ct;
    can[mname] = new TCanvas(Form("can_cat_%s", ct.c_str()),Form("can_cat_%s", ct.c_str()));
    can[mname]->cd();


    if (ave_graphs[ct]->GetN() > 1){
      //Will go with pol0
      ff[mname] = new TF1(TString::Format("ff_%s",mname.c_str()), "pol0", xmin, upperxmax);
    } else {
      ff[mname] = new TF1(TString::Format("ff_%s",mname.c_str()), "pol0", xmin, upperxmax);
    }

    ave_graphs[ct]->Fit(TString::Format("ff_%s",mname.c_str()), "R");
    ave_graphs[ct]->Draw("APE");
    ff[mname]->Draw("same");

    RooSpline1D *effMXSpline = graphToSpline(Form("eff_%s",ct.c_str()), ave_graphs[ct] , MH , xmin, upperxmax);
    effMXSplines[ct] = effMXSpline;

    TPaveText *pt2 = new TPaveText(.6,.65,.9,0.8,"NDC");
    pt2->AddText(Form("%s GeV Eff(mX) Polynomial",mname.c_str()));
    pt2->Draw();

    //make the polynomial for roofit
    if (ave_graphs[ct]->GetN() > 1){
      // eff_coef = new RooArgList( RooFit::RooConst(ff[mname]->GetParameter(0)) , RooFit::RooConst(ff[mname]->GetParameter(1)) , RooFit::RooConst(ff[mname]->GetParameter(2)) );
      //Since we go with pol0
      eff_coef = new RooArgList( RooFit::RooConst(ff[mname]->GetParameter(0)) );
    } else {
      eff_coef = new RooArgList( RooFit::RooConst(ff[mname]->GetParameter(0)) );
    }

    effMX[ct] = new RooPolyVar( Form("eff_%s",ct.c_str()), Form("eff_%s",ct.c_str()), *MH, *eff_coef, 0);

    can[mname]->SaveAs(Form("%s/EFF_MX_%s_%s.png",outplots.c_str(), insigname.c_str(), mname.c_str()));


    // can[mname]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/EFF_MX_%s.png",mname.c_str()));

    // mname += "spline";
    // can[mname] = new TCanvas(Form("can_cat_%s_spline", ct.c_str()),Form("can_cat_%s_spline", ct.c_str()));
    // can[mname]->cd();

    // //For debugging to see what this spline looks like
    // plotsplines(can[mname], exASplines[ct],  MH ,xmin, upperxmax, plabel);
    // plotsplines( TCanvas* cs, RooSpline1D * Spline, RooRealVar* MH, double xmin, double upperxmax, ){

    // can[mname]->SaveAs(Form("%s/EFF_MX_%s_%s.png",outplots.c_str(), insigname.c_str(), mname.c_str()));

  //     cans[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/out


  }

  //========================================================================
  //========================================================================
  //FOURTH PART: Build the accept(mX,kpl) RooPolyVar object
  //========================================================================
  //========================================================================
  //========================================================================

  //Do the fit per coupling and then do the fit per parameter of that.
  std::map<std::string , TCanvas*> cans;
  //In the following parameter we will save the p0,p1,p2 graphs for each
  //category, so it is [cat][0] is p0 and so on.
  std::map<std::string , std::vector<TGraphErrors*> > acc_p;
  std::map<std::string , std::vector<TF1*> > fit_acc_p;
  std::map<std::string , std::map<std::string , RooSpline1D *> >accSplines;//[coup][cat][RooSpline1D  *]

  std::string plabel;
  std::map<std::string , std::vector<double> > p0vals, p1vals, p2vals, p3vals, p0valsErr, p1valsErr, p2valsErr, p3valsErr, xvals, xvalsErr;

  //remap names to match the ones in the signal pdf
  std::map<std::string , std::string > remapnames;
  remapnames["kMpl001"] = "kMpl001";
  remapnames["kMpl01"] = "kMpl01";
  remapnames["kMpl02"] = "kMpl02";
  remapnames["0p014"] = "0p014";
  remapnames["1p4"] = "1p4";
  remapnames["5p6"] = "5p6";
  remapnames["BB"] = "EBEB";
  remapnames["BE"] = "EBEE";
  remapnames["Total"] = "All";

  for(auto grs : graphsofacc) {
    for(auto gr : grs.second) {

      double curcoup = 0.;
      //6000,9000,9000 for 17/18
      if ( gr.coup == "kMpl001" || gr.coup == "0p014"){upperxmax = 3500 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4}"; curcoup = 1.4 * pow(10.,-4);}
      else if ( gr.coup == "kMpl01" || gr.coup == "1p4"){upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2}"; curcoup = 1.4 * pow(10.,-2);}
      else if ( gr.coup == "kMpl02" || gr.coup == "5p6"){upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2}"; curcoup = 5.6 * pow(10.,-2);}
      else {
	std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
	exit(1);
      }

      std::cout << "Cat " << gr.cat << " Coup " << gr.coup << " curcoup " << curcoup << " funp0 " << gr.funp0 << " funp1 " << gr.funp1 << " funp2 " << gr.funp2 << std::endl;

      p0vals[remapnames[gr.cat]].push_back( gr.funp0 );
      p1vals[remapnames[gr.cat]].push_back( gr.funp1 );
      p2vals[remapnames[gr.cat]].push_back( gr.funp2 );
      p3vals[remapnames[gr.cat]].push_back( gr.funp3 );
      xvals[remapnames[gr.cat]].push_back( curcoup );
      p0valsErr[remapnames[gr.cat]].push_back( 0. );
      p1valsErr[remapnames[gr.cat]].push_back( 0. );
      p2valsErr[remapnames[gr.cat]].push_back( 0. );
      p3valsErr[remapnames[gr.cat]].push_back( 0. );
      xvalsErr[remapnames[gr.cat]].push_back( 0. );

      accSplines[gr.coup][remapnames[gr.cat]] = graphToSpline(Form("acc_%s_%s",gr.coup.c_str(),remapnames[gr.cat].c_str()), gr.thegraph , MH , xmin, upperxmax);

    }
  }

  for (auto cat : cats){

    // A std::vector is at its heart an array. To get the array just get the address of the first element.
    acc_p[cat].push_back( new TGraphErrors( p0vals[cat].size(), &xvals[cat][0], &p0vals[cat][0], &xvalsErr[cat][0], &p0valsErr[cat][0] ) );
    acc_p[cat].push_back( new TGraphErrors( p1vals[cat].size(), &xvals[cat][0], &p1vals[cat][0], &xvalsErr[cat][0], &p0valsErr[cat][0] ) );
    acc_p[cat].push_back( new TGraphErrors( p2vals[cat].size(), &xvals[cat][0], &p2vals[cat][0], &xvalsErr[cat][0], &p0valsErr[cat][0] ) );
    acc_p[cat].push_back( new TGraphErrors( p3vals[cat].size(), &xvals[cat][0], &p3vals[cat][0], &xvalsErr[cat][0], &p0valsErr[cat][0] ) );

  }

  //This is for the vs coupling description of p0,p1,p2.
  std::map<std::string, std::vector<RooPolyVar *> > acckMpl; //[cat][ RooPolyVar *]
  std::map<std::string, std::vector<RooArgList *> > coeff_coeffs;
  //And this is the list that will take the above kMpl dependence
  std::map<std::string, RooArgList *> acc_coeffs;

  for (auto cat : cats){

    //This is hardcoded
    xmin = 1.39 * pow(10.,-4);
    upperxmax = 5.61 * pow(10.,-2);

    //p0
    cans[Form("%s_p0",cat.c_str())] = new TCanvas(Form("%s_p0",cat.c_str()), Form("%s_p0",cat.c_str()) );
    cans[Form("%s_p0",cat.c_str())]->cd();
    fit_acc_p[cat].push_back( new TF1(Form("fit_acc_%s_p0",cat.c_str()), "pol2", xmin, upperxmax ) );
    acc_p[cat.c_str()][0]->Fit(Form("fit_acc_%s_p0",cat.c_str()), "R");
    acc_p[cat.c_str()][0]->GetYaxis()->SetTitle("p0");
    acc_p[cat.c_str()][0]->GetXaxis()->SetTitle("kMpl");
    if (cat == "EBEB") {acc_p[cat.c_str()][0]->GetYaxis()->SetRangeUser(0.20, 0.29);}
    else if (cat == "EBEE"){acc_p[cat.c_str()][0]->GetYaxis()->SetRangeUser(0.20, 0.29);}
    else if (cat == "All") {acc_p[cat.c_str()][0]->GetYaxis()->SetRangeUser(0.46, 0.52);}
    acc_p[cat.c_str()][0]->Draw("APE");
    fit_acc_p[cat][0]->Draw("same");

    cans[Form("%s_p0",cat.c_str())]->SaveAs(Form("%s/acc_%s_%s_p0.png",outplots.c_str(),insigname.c_str(),cat.c_str()));

    // cans[Form("%s_p0",cat.c_str())]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/acc_%s_p0.png",cat.c_str()));

    coeff_coeffs[cat].push_back( new RooArgList( RooFit::RooConst(fit_acc_p[cat][0]->GetParameter(0)) , RooFit::RooConst(fit_acc_p[cat][0]->GetParameter(1)) , RooFit::RooConst(fit_acc_p[cat][0]->GetParameter(2)) ) );

    acckMpl[cat].push_back( new RooPolyVar( Form("acc_%s_0",cat.c_str()), Form("acc_%s_0",cat.c_str()), *kmpl, *coeff_coeffs[cat][0], 0) );

    //p1
    cans[Form("%s_p1",cat.c_str())] = new TCanvas(Form("%s_p1",cat.c_str()), Form("%s_p1",cat.c_str()) );
    cans[Form("%s_p1",cat.c_str())]->cd();
    fit_acc_p[cat].push_back( new TF1(Form("fit_acc_%s_p1",cat.c_str()), "pol2", xmin, upperxmax ) );
    acc_p[cat.c_str()][1]->Fit(Form("fit_acc_%s_p1",cat.c_str()), "R");
    acc_p[cat.c_str()][1]->GetYaxis()->SetTitle("p1");
    acc_p[cat.c_str()][1]->GetXaxis()->SetTitle("kMpl");
    if (cat == "EBEB") {acc_p[cat.c_str()][1]->GetYaxis()->SetRangeUser(0.11 * pow(10., -3.), 0.17 * pow(10., -3.) );}
    else if (cat == "EBEE") {acc_p[cat.c_str()][1]->GetYaxis()->SetRangeUser(-70. * pow(10., -6.), -58. * pow(10., -6.) );}
    else if (cat == "All") {acc_p[cat.c_str()][1]->GetYaxis()->SetRangeUser(60. * pow(10., -6.), 99. * pow(10., -6.) );}
    acc_p[cat.c_str()][1]->Draw("APE");
    fit_acc_p[cat][1]->Draw("same");

    cans[Form("%s_p1",cat.c_str())]->SaveAs(Form("%s/acc_%s_%s_p1.png",outplots.c_str(),insigname.c_str(),cat.c_str()));

    // cans[Form("%s_p1",cat.c_str())]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/acc_%s_p1.png",cat.c_str()));

    coeff_coeffs[cat].push_back( new RooArgList( RooFit::RooConst(fit_acc_p[cat][1]->GetParameter(0)) , RooFit::RooConst(fit_acc_p[cat][1]->GetParameter(1)) , RooFit::RooConst(fit_acc_p[cat][1]->GetParameter(2)) ) );

    acckMpl[cat].push_back( new RooPolyVar( Form("acc_%s_1",cat.c_str()), Form("acc_%s_1",cat.c_str()), *kmpl, *coeff_coeffs[cat][1], 0) );

    //p2
    cans[Form("%s_p2",cat.c_str())] = new TCanvas(Form("%s_p2",cat.c_str()), Form("%s_p2",cat.c_str()) );
    cans[Form("%s_p2",cat.c_str())]->cd();
    fit_acc_p[cat].push_back( new TF1(Form("fit_acc_%s_p2",cat.c_str()), "pol2", xmin, upperxmax ) );
    acc_p[cat.c_str()][2]->Fit(Form("fit_acc_%s_p2",cat.c_str()), "R");
    acc_p[cat.c_str()][2]->GetYaxis()->SetTitle("p2");
    acc_p[cat.c_str()][2]->GetXaxis()->SetTitle("kMpl");
    if (cat == "EBEB") {acc_p[cat.c_str()][2]->GetYaxis()->SetRangeUser(-16. * pow(10., -9.), -7. * pow(10., -9.) );}
    else if (cat == "EBEE") {acc_p[cat.c_str()][2]->GetYaxis()->SetRangeUser(3.8 * pow(10., -9.), 5.5 * pow(10., -9.) );}
    else if (cat == "All") {acc_p[cat.c_str()][2]->GetYaxis()->SetRangeUser(1. * pow(10., -10.), -8. * pow(10., -9.) );}

    acc_p[cat.c_str()][2]->Draw("APE");
    fit_acc_p[cat][2]->Draw("same");

    cans[Form("%s_p2",cat.c_str())]->SaveAs(Form("%s/acc_%s_%s_p2.png",outplots.c_str(),insigname.c_str(),cat.c_str()));

    // cans[Form("%s_p2",cat.c_str())]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/acc_%s_p2.png",cat.c_str()));

    coeff_coeffs[cat].push_back( new RooArgList( RooFit::RooConst(fit_acc_p[cat][2]->GetParameter(0)) , RooFit::RooConst(fit_acc_p[cat][2]->GetParameter(1)) , RooFit::RooConst(fit_acc_p[cat][2]->GetParameter(2)) ) );

    acckMpl[cat].push_back( new RooPolyVar( Form("acc_%s_2",cat.c_str()), Form("acc_%s_2",cat.c_str()), *kmpl, *coeff_coeffs[cat][2], 0) );

    //p3
    // cans[Form("%s_p3",cat.c_str())] = new TCanvas(Form("%s_p3",cat.c_str()), Form("%s_p3",cat.c_str()) );
    // cans[Form("%s_p3",cat.c_str())]->cd();
    // fit_acc_p[cat].push_back( new TF1(Form("fit_acc_%s_p3",cat.c_str()), "pol2", xmin, upperxmax ) );
    // acc_p[cat.c_str()][3]->Fit(Form("fit_acc_%s_p3",cat.c_str()), "R");
    // acc_p[cat.c_str()][3]->GetYaxis()->SetTitle("p3");
    // acc_p[cat.c_str()][3]->GetXaxis()->SetTitle("kMpl");
    // if (cat == "EBEB") {acc_p[cat.c_str()][3]->GetYaxis()->SetRangeUser(-16. * pow(10., -9.), -7. * pow(10., -9.) );}
    // else if (cat == "EBEE") {acc_p[cat.c_str()][3]->GetYaxis()->SetRangeUser(3.8 * pow(10., -9.), 5.5 * pow(10., -9.) );}
    // else if (cat == "All") {acc_p[cat.c_str()][3]->GetYaxis()->SetRangeUser(1. * pow(10., -10.), -8. * pow(10., -9.) );}

    // acc_p[cat.c_str()][3]->Draw("APE");
    // fit_acc_p[cat][3]->Draw("same");

    // cans[Form("%s_p3",cat.c_str())]->SaveAs(Form("%s/acc_%s_%s_p3.png",outplots.c_str(),insigname.c_str(),cat.c_str()));

    // coeff_coeffs[cat].push_back( new RooArgList( RooFit::RooConst(fit_acc_p[cat][3]->GetParameter(0)) , RooFit::RooConst(fit_acc_p[cat][3]->GetParameter(1)) , RooFit::RooConst(fit_acc_p[cat][3]->GetParameter(2)) ) );

    // acckMpl[cat].push_back( new RooPolyVar( Form("acc_%s_3",cat.c_str()), Form("acc_%s_3",cat.c_str()), *kmpl, *coeff_coeffs[cat][3], 0) );

    acc_coeffs[cat] =  new RooArgList( *acckMpl[cat][0], *acckMpl[cat][1], *acckMpl[cat][2]);


  }

  //========================================================================
  //========================================================================
  //FIFTH PART : HERE WE WILL BUILD Eff(mX) x accept(mX,kpl)
  //========================================================================
  //========================================================================

  //Check also slides 12 and 13 here:
  //https://indico.cern.ch/event/458780/contributions/1971886/attachments/1179583/1707075/Chiara_diphotOct30.pdf

  //Ready to create the exA RooPolyVar
  std::map<std::string, RooPolyVar *> acceptanceMHandkMpl; //[cat][ RooPolyVar *]
  std::map<std::string, RooProduct *> exAs;
  std::map<std::string, std::map<std::string, RooProduct *> > exAsSplines;

  for (auto cat : cats){

    acceptanceMHandkMpl[cat] = new RooPolyVar( Form("acc_%s",cat.c_str()), Form("acc_%s", cat.c_str()), *MH, *acc_coeffs[cat],0);

    exAs[cat] = new RooProduct( Form("eff_acc_%s", cat.c_str()), Form("eff_acc_%s", cat.c_str()) , RooArgList(*effMX[cat],*acceptanceMHandkMpl[cat]) );

    for (auto coup : coups){
      effMXSplines[cat]->Print();
      accSplines[coup][cat]->Print();
      exAsSplines[coup][cat] = new RooProduct( Form("eff_acc_%s_%s_splines", coup.c_str(), cat.c_str()), Form("eff_acc_%s_%s_splines", coup.c_str(), cat.c_str()) , RooArgList(*effMXSplines[cat],*accSplines[coup][cat]) );
    }

  }

  //Let's see if what we build makes sense
  plotexAs(exAs, MH, kmpl, cats, coups, insigname, outplots);
  plotexAs(exAsSplines, MH, kmpl, cats, coups, insigname, outplots);


  //========================================================================
  //========================================================================
  //SIXTH PART : Read the parametric signal model with its variables and
  //             write the output workspace with norm, pdf, nuisances and
  //             renamed variables.
  //========================================================================
  //========================================================================
  std::map<std::string, RooRealVar* > reparamVars;
  // RooRealVar* obsIn;
  // RooCustomizer* custom;
  std::map<std::string, RooFormulaVar* > fwhm_rooformula; //[cat][formula]

  //Save the normalization for later
  std::ofstream myfileNorm;
  //Will save on top directory everything and not per year
  myfileNorm.open( "SignalNorm.txt", std::ios::out | std::ios::app );

  for (auto curpdf : thepdfs){

    std::cout << " Cat and Coup " << curpdf.cat << " " << curpdf.coup << std::endl;

    // fin[curpdf.coup]->cd();

    std::string signame = insigname + "_" + curpdf.coup;

    double curcoup = 0.;

    if ( curpdf.coup == "kMpl001" || curpdf.coup == "0p014"){curcoup = 1.4 * pow(10.,-4);}
    else if ( curpdf.coup == "kMpl01" || curpdf.coup == "1p4"){curcoup = 1.4 * pow(10.,-2);}
    else if ( curpdf.coup == "kMpl02" || curpdf.coup == "5p6"){curcoup = 5.6 * pow(10.,-2);}
    else {
      std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
      exit(1);
    }


    kmpl->setVal(curcoup);

    //----------------------------------------------------------------------
    //Read variables and rename them
    // for (auto par : params[curpdf.cat]){
    //   RooRealVar* dst = new RooRealVar(reparam_by_cat[par].c_str(), reparam_by_cat[par].c_str(), 0.);
    //   dst->setConstant(true);
    //   reparamVars[par] = dst;
    //   ws_out[curpdf.coup]->import(*dst, RooFit::RecycleConflictNodes());
    // }

    // curpdf.pdf = wsin[curpdf.coup]->pdf(curpdf.name.c_str());
    // obsIn =  wsin[curpdf.coup]->var("mgg");
    // obsIn->setBins(nBinsMass);
    // RooAbsBinning& mggBins =  obsIn->getBinning();
    // obsIn->setBinning(mggBins);


    //Clone the current pdf
    // custom = new RooCustomizer(*curpdf.pdf,"");
    // // custom->replaceArg(*obsIn,*obsCat[curpdf.cat]);
    // custom->replaceArg(*obsIn,*obsCat);

    // for (auto par: reparamVars){
    //   RooRealVar* srcVar = wsin[curpdf.coup]->var(par.first.c_str());
    //   par.second->setVal( srcVar->getVal() );
    //   //par.first is the existing in workspace, while the reparam name is the par.second
    //   custom->replaceArg( *srcVar  , *par.second );
    // }

    //----------------------------------------------------------------------
    //Energy scale
    // RooAbsArg * mean_formula = wsin[curpdf.coup]->arg( Form("mean_%s", curpdf.cat.c_str() ) ) ;
    // rooList[curpdf.cat] = new RooArgList( *mean_formula );

    // thetaScale[curpdf.cat] = new RooRealVar(Form("thetaScale%s_13TeV", curpdf.cat.c_str() ), Form("thetaScale%s_13TeV", curpdf.cat.c_str() ), 0., -4, 4);
    // thetaScale[curpdf.cat]->setConstant(true);
    // deltaScale[curpdf.cat] = new RooRealVar(Form("deltaScale%s_13TeV", curpdf.cat.c_str() ), Form("deltaScale%s_13TeV", curpdf.cat.c_str() ), 1.*unc);
    // deltaScale[curpdf.cat]->setConstant(true);

    // rooList[curpdf.cat]->add(*thetaScale[curpdf.cat]);
    // rooList[curpdf.cat]->add(*deltaScale[curpdf.cat]);

    // scaledMean[curpdf.cat] = new RooFormulaVar(Form("scaled_mean_%s",curpdf.cat.c_str() ), Form("scaled_mean_%s",curpdf.cat.c_str() ), "@0*(1.+@0*@1)", *rooList[curpdf.cat] );

    // custom->replaceArg( *mean_formula ,*scaledMean[curpdf.cat] );

    // curpdf.pdf = (RooAbsPdf*) custom->build(true);

    // //----------------------------------------------------------------------
    // //Saving fhwm parametrization
    // int cc = 0;
    // if (curpdf.cat == "EBEB"){cc=0;}
    // else if (curpdf.cat == "EBEE"){cc=1;}
    // else if (curpdf.cat == "All"){cc=2;}
    // RooAbsArg * fhwm_formula = wsin[curpdf.coup]->arg( Form("FHWM_%s",curpdf.cat.c_str()) );
    // fhwm_formula->Print();
    // RooRealVar* pp0 = dynamic_cast<RooRealVar*>( wsin[curpdf.coup]->arg( Form("p0_cat%d", cc) )  );
    // RooRealVar* pp1 = dynamic_cast<RooRealVar*>( wsin[curpdf.coup]->arg( Form("p1_cat%d", cc) )  );
    // double p0 = pp0->getVal();
    // double p1 = pp1->getVal();

    // // fhwm parametrized vs MH
    // fwhm_rooformula[curpdf.cat] = new RooFormulaVar(Form("fwhm_rooformula_%s",curpdf.cat.c_str()), Form("fwhm_rooformula_%s",curpdf.cat.c_str()), Form("%f + %f*@0", p0, p1), RooArgList(*MH) );

    // fwhm_rooformula[curpdf.cat]->Print();
    // ws_out[curpdf.coup]->import(*fwhm_rooformula[curpdf.cat], RooFit::RecycleConflictNodes());

    //----------------------------------------------------------------------
    //Signal normalization and naming is here

    //Remember that in the creation of the json file with e and A values we have
    //used the weightAll cut to normalize the signal samples to 1000/pb luminosity.
    //So, here we will multiply with luminosity to get to the correct norm pdf.
    //No need for the spline aka Hgg creation below.

    // curpdf.pdf->SetName( Form("model_signal_%s_%s",signame.c_str(), curpdf.cat.c_str()) );
    // RooProduct *norm = new RooProduct(Form("model_signal_%s_%s_norm", signame.c_str(), curpdf.cat.c_str()), Form("model_signal_%s_%s_norm", signame.c_str(), curpdf.cat.c_str()), RooArgList( *exAs[curpdf.cat], RooFit::RooConst(luminosity[year]) ) );

    // ws_out[curpdf.coup]->import(*norm, RooFit::RecycleConflictNodes());
    // ws_out[curpdf.coup]->import(*curpdf.pdf, RooFit::RecycleConflictNodes());

    // make some debug checks
    for (int m =500; m<5000; m=m+500) {
      MH->setVal(m);
      std::cout << "---------------------------------------------------------" << std::endl;
      std::cout << "[INFO] MH " << m <<  " kmpl " << curcoup << " - exA "  <<  (exAs[curpdf.cat]->getVal())
		<< " intL = " << luminosity[year]
		<< " predicted events " <<  exAs[curpdf.cat]->getVal() * luminosity[year]
		<< " predicted events " <<  exAsSplines[curpdf.coup][curpdf.cat]->getVal() * luminosity[year] << std::endl;
    }

    for (double mh=500.; mh<(6000.+0.25); mh+=50.){
      MH->setVal(mh);
      //Year Coupling MassPoint Category Norm
      myfileNorm << Form("%s %s %i %s %f\n", year.c_str(), curpdf.coup.c_str(), int(mh), curpdf.cat.c_str(), exAs[curpdf.cat]->getVal() * luminosity[year] ) ;
    }


  }

  myfileNorm.close();

  // for (auto cp : coups){

  //   fout[cp]->cd();

  //   ws_out[cp]->Write();
  //   // ws_out[cp]->Print();

  //   fout[cp]->Write();
  //   fout[cp]->Close();

  // }


  //========================================================================
  //========================================================================
  //BACKUP PART : This part contains commented code with some effort to
  //              build splines for cross sections and a more Hgg style norm
  //              pdf but it doesn't include the kMpl so I comment it out.
  //========================================================================
  //========================================================================

  // std::map<std::string, std::vector<xsec> > xsections = loadXsections(year, true);
  // std::map<std::string, std::vector<xsec> > xsectionsKMpl = loadXsections(year, false);
  // std::map<std::string, TGraphErrors*> grxs;
  // std::map<std::string, RooSpline1D *> xsSplines;
  // // std::map<std::string, RooSpline1D *> xsSplineskMpl;
  // std::map<std::string, RooPolyVar *> xskMpl;
  // std::map<std::string, TGraphErrors*> kMplgrs; //[MH]
  // std::map<std::string, TF1 *> fm;

  // for(auto cp : coups) {

  //   can[cp] = new TCanvas(Form("can_coup_%s", cp.c_str()),Form("can_coup_%s", cp.c_str()));

  //   if ( cp == "kMpl001" || cp == "0p014"){upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4}"; }
  //   else if ( cp == "kMpl01" || cp == "1p4"){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2}"; }
  //   else if ( cp == "kMpl02" || cp == "5p6"){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2}"; }
  //   else {
  //     std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
  //     exit(1);
  //   }

  //   grxs[cp] = getXsecGraph(xsections[cp], true);

  //   std::ofstream myfile;
  //   myfile.open( Form("%s/XS_%s.txt", outplots.c_str(), cp.c_str()) , std::ios::out );

  //   for (double mh=500.; mh<(upperxmax+0.25); mh+=50.){
  //     std::cout << std::scientific << std::setprecision(5) << grxs[cp]->Eval(mh) << std::endl;
  //     myfile << Form("%f %.5e \n", mh, grxs[cp]->Eval(mh) ) ;

  //   }

  //   RooSpline1D *xsSpline = graphToSpline(Form("fxs_%s",cp.c_str()), grxs[cp], MH , xmin, upperxmax);

  //   // for (double mh=500.; mh<(upperxmax+0.25); mh+=100.){
  //   //   MH->setVal(mh);
  //   //   std::cout << xsSpline->getVal() << std::endl;
  //   //   myfile << std::setprecision(12) << Form("%f %f \n", mh, xsSpline->getVal() ) ;
  //   // }

  //   myfile.close();

  //   xsSplines[cp] = xsSpline;

  //   std::cout << outplots.c_str() << std::endl;
  //   //For debugging to see what this spline looks like
  //   plotsplines(can[cp], xsSplines[cp],  MH ,xmin, upperxmax, plabel);
  //   can[cp]->SaveAs(Form("%s/xsSplines_%s.png", outplots.c_str(), cp.c_str()));

  // }

  // RooArgList *xs_coef;
  // for(auto mb : xsectionsKMpl){

  //   std::string mname = mb.first;
  //   //This is hardcoded
  //   xmin = 1.39 * pow(10.,-4);
  //   upperxmax = 5.61 * pow(10.,-2);

  //   can[mname] = new TCanvas( Form("can%s",mname.c_str()), Form("can%s",mname.c_str()) );
  //   can[mname]->cd();

  //   kMplgrs[mname] = getXsecGraph(xsectionsKMpl[mname], false);
  //   // RooSpline1D *xsSplineKMpl = graphToSpline(Form("fxs_%s",mname.c_str()), kMplgrs[mname], kmpl , xmin, upperxmax);
  //   // xsSplineskMpl[mname] = xsSplineKMpl;
  //   // std::cout << kMplgrs[mname]->GetN() << std::endl;

  //   if (kMplgrs[mname]->GetN() > 1){
  //     fm[mname] = new TF1(TString::Format("fm_%s",mname.c_str()), "pol2", xmin, upperxmax);
  //   } else {
  //     fm[mname] = new TF1(TString::Format("fm_%s",mname.c_str()), "pol0", xmin, upperxmax);
  //   }

  //   kMplgrs[mname]->Fit(TString::Format("fm_%s",mname.c_str()), "R");
  //   kMplgrs[mname]->Draw("APE");
  //   fm[mname]->Draw("same");

  //   TPaveText *pt2 = new TPaveText(.6,.65,.9,0.8,"NDC");
  //   pt2->AddText(Form("%s GeV XS Polynomial",mname.c_str()));
  //   pt2->Draw();

  //   //make the polynomial for roofit
  //   if (kMplgrs[mname]->GetN() > 1){
  //     xs_coef = new RooArgList( RooFit::RooConst(fm[mname]->GetParameter(0)) , RooFit::RooConst(fm[mname]->GetParameter(1)) , RooFit::RooConst(fm[mname]->GetParameter(2)) );
  //   } else {
  //     xs_coef = new RooArgList( RooFit::RooConst(fm[mname]->GetParameter(0)) );
  //   }

  //   xskMpl[mname] = new RooPolyVar( Form("xs_%s",mname.c_str()), Form("xs_%s",mname.c_str()), *kmpl, *xs_coef, 0);


  //   //For debugging to see what this spline looks like
  //   // plotsplines(can[mname], xsSplineskMpl[mname],  kmpl ,xmin, upperxmax, Form("%s XS Spline",mname.c_str()) );

  //   can[mname]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/xs_%s.png",mname.c_str()));

  // }

  // //Graphs of exA to splines
  // std::map<std::string , RooSpline1D *> exASplines;
  // std::map<std::string , RooAbsReal *> pdf_norm;


  // for(auto grs : graphsofexA) {
  //   for(auto gr : grs.second) {
  //     cans[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())] = new TCanvas(Form("%s_%s", gr.coup.c_str(), gr.cat.c_str()), Form("%s_%s", gr.coup.c_str(), gr.cat.c_str()));


  //     if ( gr.coup == "001" ){upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4}";}
  //     else if ( gr.coup == "01" ){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2}";}
  //     else if ( gr.coup == "02" ){upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2}";}
  //     else {
  // 	std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
  // 	exit(1);
  //     }

  //     RooSpline1D *exASpline = graphToSpline(Form("fea_%s_%s",gr.coup.c_str(), gr.cat.c_str() ), gr.thegraph , MH , xmin, upperxmax);
  //     exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())] = exASpline;

  //     //For debugging to see what this spline looks like
  //     plotsplines(cans[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())], exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())],  MH ,xmin, upperxmax, plabel);

  //     cans[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exASplines_%s_%s.png",gr.coup.c_str(), gr.cat.c_str()));

  //     //Now we have the tools to create the norm pdf. We should give it the correct name also.
  //     std::string normpdfname = Form("SignalShape_%s_%s_norm",remapnames[gr.coup].c_str(), remapnames[gr.cat].c_str());
  //     RooAbsReal *finalNorm = new RooFormulaVar( normpdfname.c_str(), normpdfname.c_str(), "@0*@1*@2",RooArgList(*xsSplines[gr.coup.c_str()], *exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())], RooFit::RooConst(luminosity[year]) ) );

  //     pdf_norm[normpdfname] = finalNorm;

  //     // make some debug checks
  //     for (int m =500; m<9000; m=m+500){
  // 	MH->setVal(m);
  // 	std::cout << "[INFO] MH " << m <<  " - ea "  <<  (exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())]->getVal())
  // 		  << " intL= " << luminosity[year]
  // 		  << " xs " << xsSplines[gr.coup.c_str()]->getVal()
  // 		  << "norm " << (exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())]->getVal())*xsSplines[gr.coup.c_str()]->getVal()
  // 		  << "predicted events " <<  (exASplines[Form("%s_%s", gr.coup.c_str(), gr.cat.c_str())]->getVal())*xsSplines[gr.coup.c_str()]->getVal()*luminosity[year] <<  std::endl;
  //     }

  //     }//end of loop through graphs
  // }//end of loop through couplings

  // for (auto pdf: pdf_norm){
  //   ws_out->import(*pdf.second);
  // }

  // ws_out->Print();
  // ws_out->Write();

  // fout->Write();
  // fout->Close();

  // plotefficiency(cc, effrecoBB, false);
  // cc->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exAsVsMass_001.png");
  // TCanvas* ccBE = new TCanvas("ccBE", "ccBE");
  // plotefficiency(ccBE, effrecoBE, false);
  // ccBE->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effVsMass_BE.png");
  // TCanvas* ccave = new TCanvas("ccave", "ccave");
  // plotefficiency(ccave, effreco, true);
  // ccave->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effVsMass_average.png");

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
std::map<std::string, std::vector<eff_reco> > computefficiency(const std::string &year)
{

  std::map<std::string, std::string> cuts;
  //Acceptance cuts
  cuts["ABB"] = "(Diphoton.Minv > 230 && Photon1.pt>125 && Photon2.pt>125 && Photon1.isEB && Photon2.isEB)*weightAll";
  cuts["ABE"] = "(Diphoton.Minv > 330 && Photon1.pt>125 && Photon2.pt>125 && ( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE )))*weightAll";
  cuts["ATotal"] = "((Diphoton.Minv > 230 && Photon1.pt>125 && Photon2.pt>125) && ( (Photon1.isEB && Photon2.isEB) || (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ) ) )*weightAll";
  // std::cout << cuts["ATotal"] << std::endl;
  //Selection efficiency cuts.
  if ( year == "2016"){

    cuts["eBB"] = "(Photon1.isEB && Photon2.isEB)*isGood*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)*weightAll";
    cuts["eBE"] = "( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ))*isGood*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)*weightAll";
    cuts["eTotal"] = "( (Photon1.isEB && Photon2.isEB) || (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ) )*isGood*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)*weightAll";
    // std::cout << cuts["eTotal"] << std::endl;
    //Acceptance times efficiency
    cuts["exABB"] = cuts["ABB"] + "*isGood*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)";
    cuts["exABE"] = cuts["ABE"] + "*isGood*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)";
    cuts["exATotal"] = "(" + cuts["ATotal"] + ")" + "*isGood*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)";
    // std::cout << cuts["exATotal"] << std::endl;

    // cuts["BB"] = "isGood*(Diphoton.Minv > 230 && Photon1.pt>125 && Photon2.pt>125 && Photon1.isEB && Photon2.isEB)*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)";
    // cuts["BE"] = "isGood*(Diphoton.Minv > 330 && Photon1.pt>125 && Photon2.pt>125 && ( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE )))*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)";

  } else{
  cuts["eBB"] = "(Photon1.isEB && Photon2.isEB)*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)*weightAll";
  cuts["eBE"] = "( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ))*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)*weightAll";
  cuts["eTotal"] = "( (Photon1.isEB && Photon2.isEB) || (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE ) )*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)*weightAll";
  // std::cout << cuts["eTotal"] << std::endl;
  //Acceptance times efficiency
  cuts["exABB"] = cuts["ABB"] + "*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  cuts["exABE"] = cuts["ABE"] + "*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  cuts["exATotal"] = "(" + cuts["ATotal"] + ")" + "*isGood*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  // std::cout << cuts["exATotal"] << std::endl;

  // cuts["BB"] = "isGood*(Diphoton.Minv > 230 && Photon1.pt>125 && Photon2.pt>125 && Photon1.isEB && Photon2.isEB)*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  // cuts["BE"] = "isGood*(Diphoton.Minv > 330 && Photon1.pt>125 && Photon2.pt>125 && ( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE )))*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
  }

  std::vector<std::string> samples = getSampleList();

  std::vector<int> stringScales = {3000, 3500, 4000, 4500, 5000, 5500, 6000};

  std::string coupling = "";
  std::string M_bins = "";
  double NacBB = 0.;double NacBE = 0.;double NacTotal = 0.;

  std::map<std::string, std::vector<eff_reco> > the_eff_reco;
  std::map<std::string, TH1F *> histograms;
  std::string histo_name;

  for(auto isample : samples) {
    //Run a single year each time
    if ( isample.find(year) == std::string::npos ) continue;
    //Run only on RS samples for now.
    // if( isample.find("RSGravitonToGammaGamma") == std::string::npos ) continue;

    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    double xsecvalue = 0.;
    if (year == "2016"){
      xsecvalue = ExoDiPhotons::crossSection(getSampleBase(isample,year)+"_"+year);
      std::cout << "!!!!! 2016 xsec value!!!!!"<< std::endl;
    } else {
      xsecvalue = ExoDiPhotons::crossSection(getSampleBase(isample,year));
    }
    std::cout << "xsec " << xsecvalue << std::endl;

    if (isample.find("RSGraviton") != std::string::npos) {
      coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    } else if ( isample.find("GluGluSpin0") != std::string::npos ){
      coupling = get_str_between_two_str(getBase(isample), "GluGluSpin0ToGammaGamma_W_", "_M_");
    }

    M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");

    int nBins = 2385;
    double xMin = 0.0;
    double xMax = 0.0;
    if ( coupling == "001" ){xMax = 6000. ; coupling = "kMpl001";}
    else if ( coupling == "01" ){xMax = 9000. ; coupling = "kMpl01";}
    else if ( coupling == "02" ){xMax = 9000. ; coupling = "kMpl02";}
    else if ( coupling == "0p014" ){ xMax = 6000. ; }
    else if ( coupling == "1p4"){ xMax = 9000. ; }
    else if ( coupling == "5p6"){ xMax = 9000.; }
    else {
      std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
      exit(1);
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //2016 upper cut
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (coupling == "kMpl001" && year == "2016" && std::stof(M_bins) > 3900. ){continue;}

    eff_reco tmpeffreco;

    std::string baseName(getSampleBase(isample, year));

    tmpeffreco.samplename = baseName;

    //========================================================================
    //Initial Distribution
    //========================================================================

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName;
    // chains[getBase(isample)]->Print();
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, "no_selection");
    TCanvas *c1 = new TCanvas();
    histograms[histo_name]->Draw("HIST");
    c1->Update();
    c1->SaveAs("test/test_%s.png",histo_name);
    // Ngen = histograms[histo_name]->Integral();
    // std::cout << "Ngen " << Ngen << std::endl;
    // std::cout << "Name " << getBase(isample) << " chains[getBase(isample)]->GetName() " << chains[getBase(isample)]->GetName() << " chains[getBase(isample)]->GetEntries() " << chains[getBase(isample)]->GetEntries() << std::endl;

    //========================================================================
    //Acceptance (gen-level selection)
    //========================================================================
    //------------------------------------------------------------------------
    //Sel: ABB
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "ABB";
    tmpeffreco.cat = "ABB";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_ABB";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["ABB"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * xsecvalue );
    //Apart from the acceptance, we also want the denominator for the efficiency below
    NacBB = histograms[histo_name]->Integral();

    the_eff_reco[coupling].push_back(tmpeffreco);

    std::cout << "A BB " << tmpeffreco.efforacc << std::endl;

    //------------------------------------------------------------------------
    //Sel: ABE
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "ABE";
    tmpeffreco.cat = "ABE";

    histo_name = baseName + "_ABE";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["ABE"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * xsecvalue );
    //Apart from the acceptance, we also want the denominator for the efficiency below
    NacBE = histograms[histo_name]->Integral();

    the_eff_reco[coupling].push_back(tmpeffreco);
    std::cout << "A BE " << tmpeffreco.efforacc << std::endl;
    //------------------------------------------------------------------------
    //Sel: ATotal
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "ATotal";
    tmpeffreco.cat = "ATotal";

    histo_name = baseName + "_ATotal";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["ATotal"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * xsecvalue );
    //Apart from the acceptance, we also want the denominator for the efficiency below
    NacTotal = histograms[histo_name]->Integral();

    the_eff_reco[coupling].push_back(tmpeffreco);
    std::cout << "A Total " << tmpeffreco.efforacc << std::endl;

    std::cout << "-----------------------------------------------------------------" << std::endl;
    //========================================================================
    //Selection Efficiency (defined as exA/A and then averaging)
    //========================================================================
    //------------------------------------------------------------------------
    //Sel: eBB
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "eBB";
    tmpeffreco.cat = "eBB";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_eBB";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exABB"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / NacBB;

    the_eff_reco[coupling].push_back(tmpeffreco);

    std::cout << "e BB " << tmpeffreco.efforacc << std::endl;

    //------------------------------------------------------------------------
    //Sel: eBE
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "eBE";
    tmpeffreco.cat = "eBE";

    histo_name = baseName + "_eBE";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exABE"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / NacBE;

    the_eff_reco[coupling].push_back(tmpeffreco);
    std::cout << "e BE " << tmpeffreco.efforacc << std::endl;
    //------------------------------------------------------------------------
    //Sel: eTotal
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "eTotal";
    tmpeffreco.cat = "eTotal";

    histo_name = baseName + "_eTotal";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exATotal"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / NacTotal;

    the_eff_reco[coupling].push_back(tmpeffreco);
    std::cout << "e Total " << tmpeffreco.efforacc << std::endl;

    //========================================================================
    //Efficiency x Acceptance (reco selection): Just for crosscheck
    //========================================================================
    //------------------------------------------------------------------------
    //Sel: exABB
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "exABB";
    tmpeffreco.cat = "exABB";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_exABB";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exABB"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * xsecvalue ) ;

    the_eff_reco[coupling].push_back(tmpeffreco);

    std::cout << "exA BB " << tmpeffreco.efforacc << std::endl;

    //------------------------------------------------------------------------
    //Sel: exABE
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "exABE";
    tmpeffreco.cat = "exABE";

    histo_name = baseName + "_exABE";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exABE"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * xsecvalue ) ;

    the_eff_reco[coupling].push_back(tmpeffreco);
    std::cout << "exA BE " << tmpeffreco.efforacc << std::endl;
    //------------------------------------------------------------------------
    //Sel: exATotal
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "exATotal";
    tmpeffreco.cat = "exATotal";

    histo_name = baseName + "_exATotal";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["exATotal"]);
    tmpeffreco.efforacc = histograms[histo_name]->Integral() / (1000. * xsecvalue ) ;

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
TGraphErrors* avegraph(std::vector<eff_reco> quant){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraphErrors* graph;

  unsigned int nP = quant.size();

  double masses[nP];
  double massesErr[nP];
  double qu[nP];
  double quErr[nP];

  for(unsigned int iE =0; iE < nP; iE++){
    masses[iE] = std::stod(quant[iE].M_bins);
    massesErr[iE] = 0.;
    qu[iE] = quant[iE].efforacc;
    quErr[iE] = 0.;
  }

  graph = new TGraphErrors(nP, masses, qu, massesErr, quErr);

  return graph;

}

//-----------------------------------------------------------------------------------
TGraphErrors* getXsecGraph( std::vector<xsec> xsections, bool basedonCoup)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraphErrors* graph;

  unsigned int nP = xsections.size();

  double xval[nP];
  double xvalErr[nP];
  double xsecval[nP];
  double xsecErr[nP];

  for(unsigned int iE =0; iE < nP; iE++){
    if (basedonCoup) { xval[iE] = std::stod(xsections[iE].M_bins); }
    else {
      xval[iE] = xsections[iE].coup ;
      std::cout << "xsections[iE].coup  " << xsections[iE].coup << std::endl;
    }
    xvalErr[iE] = 0.;
    xsecval[iE] = xsections[iE].val;
    xsecErr[iE] = xsections[iE].error;
  }

  graph = new TGraphErrors(nP, xval, xsecval, xvalErr, xsecErr);

  return graph;

}
//-----------------------------------------------------------------------------------
void plotexAs( std::map<std::string, RooProduct *> exAs, RooRealVar* MH, RooRealVar* kmpl, std::vector<std::string> cats, std::vector<std::string> coups, std::string insigname, std::string outplots){

  //We need one graph per coupling per category
  std::map<std::string, TGraph* > graph;
  std::map<std::string, TCanvas* > cs;

  double upperxmax = 0.;
  double xmin = 0.;
  std::string plabel;

  std::map<std::string, TLegend* > legmc;

  for(auto cp : coups) {
    //A counter for the canvas draw
    int count = 0;
    cs[cp] = new TCanvas(Form("can_%s",cp.c_str()), Form("can_%s",cp.c_str()));
    legmc[cp] = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
    legmc[cp]->SetTextFont(42);
    legmc[cp]->SetBorderSize(0);
    legmc[cp]->SetFillStyle(0);

    double curcoup = 0.;

    if ( cp == "kMpl001" || cp == "0p014"){
      upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4}"; curcoup = 1.4 * pow(10.,-4);
      legmc[cp]->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-4}","C");
    } else if ( cp == "kMpl01" || cp == "1p4"){
      upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2}"; curcoup = 1.4 * pow(10.,-2);
      legmc[cp]->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-2}","C");
    } else if ( cp == "kMpl02" || cp == "5p6"){
      upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2}"; curcoup = 5.6 * pow(10.,-2);
      legmc[cp]->SetHeader("#frac{#Gamma}{m} = 5.6 #times 10^{-2}","C");
    } else {
      std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
      exit(1);
    }

    kmpl->setVal(curcoup);

    for (auto cat : cats){
      cs[cp]->cd();

      ++count;

      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] = new TGraph();
      int point=0;

      for (double m =xmin; m<upperxmax; m=m+50.){
	MH->setVal(m);
	// std::cout << exAs[cat]->getVal() << std::endl;
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetPoint(point,m,exAs[cat]->getVal());
	point++;
      }

      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->GetYaxis()->SetRangeUser(0., 1.0);
      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->GetYaxis()->SetTitle("#varepsilon #otimes A");
      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->GetXaxis()->SetTitle("m_{X} [GeV]");
      if (cp == "kMpl001" || cp == "0p014") {graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerStyle(kFullDiamond);}
      else if (cp == "kMpl01" || cp == "1p4") {graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerStyle(kFullTriangleUp);}
      else if (cp == "kMpl02" || cp == "5p6") {graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerStyle(kFullCircle);}

      int thespin = 10000;
      if (insigname == "grav") {thespin=2;}
      else if (insigname == "heavyhiggs") {thespin=0;}

      if ( cat == "EBEB"){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerColor(2);
	legmc[cp]->AddEntry( graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] , Form("EBEB J=%d",thespin) ,"p" );
      } else if ( cat == "EBEE"){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerColor(4);
	legmc[cp]->AddEntry( graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] , Form("EBEE J=%d",thespin) ,"p" );
      } else if ( cat == "All"){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerColor(8);
	legmc[cp]->AddEntry( graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] , Form("Total J=%d",thespin) ,"p" );
      }

      if (count==1){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->Draw("AP");
      } else {
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->Draw("PS");
      }
      if (count==3){
	legmc[cp]->Draw("same");
      }

    }// end of loop over cats

    cs[cp]->SaveAs(Form("%s/exAs_RooProduct_%s.png", outplots.c_str(), cp.c_str()));

    // cs[cp]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exAs_RooProduct_kMpl_%s.png", cp.c_str()));

  }//end of loops over couplings

}
//-----------------------------------------------------------------------------------
void plotexAs( std::map<std::string,std::map<std::string, RooProduct *> >exAs, RooRealVar* MH, RooRealVar* kmpl, std::vector<std::string> cats, std::vector<std::string> coups, std::string insigname, std::string outplots){

  //We need one graph per coupling per category
  std::map<std::string, TGraph* > graph;
  std::map<std::string, TCanvas* > cs;

  double upperxmax = 0.;
  double xmin = 0.;
  std::string plabel;

  std::map<std::string, TLegend* > legmc;

  for(auto cp : coups) {
    //A counter for the canvas draw
    int count = 0;
    cs[cp] = new TCanvas(Form("can_%s",cp.c_str()), Form("can_%s",cp.c_str()));
    legmc[cp] = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
    legmc[cp]->SetTextFont(42);
    legmc[cp]->SetBorderSize(0);
    legmc[cp]->SetFillStyle(0);

    double curcoup = 0.;

    if ( cp == "kMpl001" || cp == "0p014"){
      upperxmax = 6000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4}"; curcoup = 1.4 * pow(10.,-4);
      legmc[cp]->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-4}","C");
    } else if ( cp == "kMpl01" || cp == "1p4"){
      upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2}"; curcoup = 1.4 * pow(10.,-2);
      legmc[cp]->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-2}","C");
    } else if ( cp == "kMpl02" || cp == "5p6"){
      upperxmax = 9000 ; plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2}"; curcoup = 5.6 * pow(10.,-2);
      legmc[cp]->SetHeader("#frac{#Gamma}{m} = 5.6 #times 10^{-2}","C");
    } else {
      std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
      exit(1);
    }

    kmpl->setVal(curcoup);

    for (auto cat : cats){
      cs[cp]->cd();

      ++count;

      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] = new TGraph();
      int point=0;

      for (double m =xmin; m<upperxmax; m=m+5.){
	MH->setVal(m);
	// std::cout << exAs[cat]->getVal() << std::endl;
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetPoint(point,m,exAs[cp][cat]->getVal());
	point++;
      }

      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->GetYaxis()->SetRangeUser(0., 1.0);
      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->GetYaxis()->SetTitle("#varepsilon #otimes A");
      graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->GetXaxis()->SetTitle("m_{X} [GeV]");
      if (cp == "kMpl001" || cp == "0p014") {graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerStyle(kFullDiamond);}
      else if (cp == "kMpl01" || cp == "1p4") {graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerStyle(kFullTriangleUp);}
      else if (cp == "kMpl02" || cp == "5p6") {graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerStyle(kFullCircle);}

      int thespin = 10000;
      if (insigname == "grav") {thespin=2;}
      else if (insigname == "heavyhiggs") {thespin=0;}

      if ( cat == "EBEB"){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerColor(2);
	legmc[cp]->AddEntry( graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] , Form("EBEB J=%d",thespin) ,"p" );
      } else if ( cat == "EBEE"){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerColor(4);
	legmc[cp]->AddEntry( graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] , Form("EBEE J=%d",thespin) ,"p" );
      } else if ( cat == "All"){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->SetMarkerColor(8);
	legmc[cp]->AddEntry( graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())] , Form("Total J=%d",thespin) ,"p" );
      }

      if (count==1){
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->Draw("AP");
      } else {
	graph[Form("exA_%s_%s", cat.c_str(), cp.c_str())]->Draw("PS");
      }
      if (count==3){
	legmc[cp]->Draw("same");
      }

    }// end of loop over cats

    cs[cp]->SaveAs(Form("%s/exAs_RooProduct_%s_Splines.png", outplots.c_str(), cp.c_str()));

    // cs[cp]->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/exAs_RooProduct_kMpl_%s.png", cp.c_str()));

  }//end of loops over couplings

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
  pt2->AddText(Form("%s Spline",label.c_str()));
  graph->Draw("ALP");
  pt2->Draw();

}
//-----------------------------------------------------------------------------------
RooSpline1D* graphToSpline(std::string name, TGraphErrors *graph, RooRealVar* MH, double xmin, double upperxmax){

  std::vector<double> xValues, yValues;
  for (double mh=xmin; mh<(upperxmax+0.25); mh+=50.){
    xValues.push_back(mh);
    yValues.push_back(graph->Eval(mh));
  }
  RooSpline1D *res = new RooSpline1D(name.c_str(),name.c_str(),*MH,xValues.size(),&(xValues[0]),&(yValues[0]));
  return res;
}

//-----------------------------------------------------------------------------------
std::vector<theGraphs> plot(TCanvas* cc, std::string coup, std::map<std::string, std::vector<eff_reco> > quantBB, std::map<std::string, std::vector<eff_reco> > quantBE , std::map<std::string, std::vector<eff_reco> > quantTotal, bool doave, std::string label, std::string insigname)
{

  std::vector<theGraphs> graphs; //BB, BE, Total
  theGraphs tmpgr;

  std::map<std::string, double> upperxmax;
  upperxmax["kMpl001"] = 3500.;//6000 for 2017/18
  upperxmax["kMpl01"]  = 6000.;//9000 for 2017/18
  upperxmax["kMpl02"]  = 6000.;//9000 for 2017/18
  upperxmax["0p014"] = 6000.;//6000 for 2017/18
  upperxmax["1p4"]  = 9000.;
  upperxmax["5p6"]  = 9000.;

  int thespin = 10000;
  if (insigname == "grav") {thespin=2;}
  else if (insigname == "heavyhiggs") {thespin=0;}

  cc->cd();

  if (!doave){
    TLegend* legmc = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);

    if (coup == "kMpl001" || coup == "0p014") {legmc->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-4}","C");}
    else if (coup == "kMpl01" || coup == "1p4") {legmc->SetHeader("#frac{#Gamma}{m} = 1.4 #times 10^{-2}","C");}
    else if (coup == "kMpl02" || coup == "5p6") {legmc->SetHeader("#frac{#Gamma}{m} = 5.6 #times 10^{-2}","C");}

    //BB
    tmpgr.coup = coup;
    tmpgr.cat = "BB";
    tmpgr.thegraph = graph(coup, quantBB);
    tmpgr.thegraph->GetYaxis()->SetRangeUser(0., 1.0);
    tmpgr.thegraph->SetMarkerColor(2);
    if (coup == "kMpl001" || coup == "0p014") {tmpgr.thegraph->SetMarkerStyle(kFullDiamond);}
    else if (coup == "kMpl01" || coup == "1p4") {tmpgr.thegraph->SetMarkerStyle(kFullTriangleUp);}
    else if (coup == "kMpl02" || coup == "5p6") {tmpgr.thegraph->SetMarkerStyle(kFullCircle);}
    tmpgr.thegraph->SetLineColor(2);
    legmc->AddEntry( tmpgr.thegraph , Form("EBEB J=%d",thespin) ,"pl" );

    tmpgr.fun = new TF1(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "pol2", 500,upperxmax[coup]);
    tmpgr.fun->SetLineColor(2);
    tmpgr.fun->SetLineStyle(10);
    tmpgr.thegraph->Fit(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "R");
    tmpgr.thegraph->GetYaxis()->SetTitle(Form("%s",label.c_str()));
    tmpgr.thegraph->GetXaxis()->SetTitle("m_{X} [GeV]");
    tmpgr.thegraph->Draw("APE");
    tmpgr.fun->Draw("same");

    tmpgr.funp0 = tmpgr.fun->GetParameter(0);
    tmpgr.funp1 = tmpgr.fun->GetParameter(1);
    tmpgr.funp2 = tmpgr.fun->GetParameter(2);
    // tmpgr.funp3 = tmpgr.fun->GetParameter(3);

    graphs.push_back(tmpgr);

    tmpgr.coup = coup;
    tmpgr.cat = "BE";
    tmpgr.thegraph = graph(coup, quantBE);
    tmpgr.thegraph->GetYaxis()->SetRangeUser(0., 1.0);
    tmpgr.thegraph->SetMarkerColor(4);
    if (coup == "kMpl001" || coup == "0p014") {tmpgr.thegraph->SetMarkerStyle(kFullDiamond);}
    else if (coup == "kMpl01" || coup == "1p4") {tmpgr.thegraph->SetMarkerStyle(kFullTriangleUp);}
    else if (coup == "kMpl02" || coup == "5p6") {tmpgr.thegraph->SetMarkerStyle(kFullCircle);}
    tmpgr.thegraph->SetLineColor(4);
    legmc->AddEntry( tmpgr.thegraph , Form("EBEE J=%d",thespin) ,"pl" );

    tmpgr.fun = new TF1(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "pol2", 500,upperxmax[coup]);
    tmpgr.fun->SetLineColor(4);
    tmpgr.fun->SetLineStyle(5);
    tmpgr.thegraph->Fit(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "R");
    tmpgr.thegraph->GetYaxis()->SetTitle(Form("%s",label.c_str()));
    tmpgr.thegraph->GetXaxis()->SetTitle("m_{X} [GeV]");
    tmpgr.thegraph->Draw("PSE");
    tmpgr.fun->Draw("same");

    tmpgr.funp0 = tmpgr.fun->GetParameter(0);
    tmpgr.funp1 = tmpgr.fun->GetParameter(1);
    tmpgr.funp2 = tmpgr.fun->GetParameter(2);
    // tmpgr.funp3 = tmpgr.fun->GetParameter(3);

    graphs.push_back(tmpgr);

    tmpgr.coup = coup;
    tmpgr.cat = "Total";
    tmpgr.thegraph = graph(coup, quantTotal);
    tmpgr.thegraph->GetYaxis()->SetRangeUser(0., 1.0);
    tmpgr.thegraph->SetMarkerColor(8);
    if (coup == "kMpl001" || coup == "0p014") {tmpgr.thegraph->SetMarkerStyle(kFullDiamond);}
    else if (coup == "kMpl01" || coup == "1p4") {tmpgr.thegraph->SetMarkerStyle(kFullTriangleUp);}
    else if (coup == "kMpl02" || coup == "5p6") {tmpgr.thegraph->SetMarkerStyle(kFullCircle);}
    tmpgr.thegraph->SetLineColor(8);
    legmc->AddEntry( tmpgr.thegraph , Form("Total J=%d",thespin) ,"pl" );

    tmpgr.fun = new TF1(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "pol2", 500,upperxmax[coup]);
    tmpgr.fun->SetLineColor(8);
    tmpgr.fun->SetLineStyle(1);
    tmpgr.thegraph->Fit(Form("f%s_%s",coup.c_str(),tmpgr.cat.c_str()), "R");
    tmpgr.thegraph->GetYaxis()->SetTitle(Form("%s",label.c_str()));
    tmpgr.thegraph->GetXaxis()->SetTitle("m_{X} [GeV]");
    tmpgr.thegraph->Draw("PSE");
    tmpgr.fun->Draw("same");

    tmpgr.funp0 = tmpgr.fun->GetParameter(0);
    tmpgr.funp1 = tmpgr.fun->GetParameter(1);
    tmpgr.funp2 = tmpgr.fun->GetParameter(2);
    // tmpgr.funp2 = tmpgr.fun->GetParameter(3);

    graphs.push_back(tmpgr);

    legmc->Draw("same");
  }

  return graphs;

}


void plotallcoups(TCanvas* ccallcoup, std::vector<theGraphs> graphsofeffkMpl001, std::vector<theGraphs> graphsofeffkMpl01, std::vector<theGraphs> graphsofeffkMpl02, std::string insigname){

  ccallcoup->cd();

  int thespin = 10000;
  if (insigname == "grav") {thespin=2;}
  else if (insigname == "heavyhiggs") {thespin=0;}

  TLegend* legmc = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);

  std::map< std::string, std::string > remapcats;
  remapcats["BB"] = "EBEB"; remapcats["BE"] = "EBEE"; remapcats["Total"] = "Total";

  int i=0;

  for (auto gr : graphsofeffkMpl001){
    if (i==0){ gr.thegraph->Draw(); ++i;}
    else { gr.thegraph->Draw("same");}
    gr.fun->Draw("same");
    legmc->AddEntry( gr.thegraph , Form("#frac{#Gamma}{m} = 1.4 #times 10^{-4} , %s J=%d", remapcats[gr.cat].c_str(), thespin) ,"pl" );
  }
  for (auto gr : graphsofeffkMpl01){
    gr.thegraph->Draw("same");
    gr.fun->Draw("same");
    legmc->AddEntry( gr.thegraph , Form("#frac{#Gamma}{m} = 1.4 #times 10^{-2} , %s J=%d", remapcats[gr.cat].c_str(), thespin ) ,"pl" );
  }
  for (auto gr : graphsofeffkMpl02){
    gr.thegraph->Draw("same");
    gr.fun->Draw("same");
    legmc->AddEntry( gr.thegraph , Form("#frac{#Gamma}{m} = 5.6 #times 10^{-2} , %s J=%d", remapcats[gr.cat].c_str(), thespin ) ,"pl" );
  }

  legmc->Draw("same");

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
std::map<std::string, std::vector<eff_reco> > reduce(std::map<std::string, std::vector<eff_reco> > effreco, std::vector<std::string> coups, std::string cut){

  std::map<std::string, std::vector<eff_reco> > out;

  for (auto cp : coups){
    for(unsigned int iE =0; iE < effreco[cp].size(); iE++){
      if (effreco[cp][iE].sel == cut){ out[cp].push_back(effreco[cp][iE]); }
    }
  }

  // for(unsigned int iE =0; iE < effreco["kMpl001"].size(); iE++){
  //   if (effreco["kMpl001"][iE].sel == cut){ out["kMpl001"].push_back(effreco["kMpl001"][iE]); }
  // }

  // for(unsigned int iE =0; iE < effreco["kMpl01"].size(); iE++){
  //   if (effreco["kMpl01"][iE].sel == cut){ out["kMpl01"].push_back(effreco["kMpl01"][iE]); }
  // }

  // for(unsigned int iE =0; iE < effreco["kMpl02"].size(); iE++){
  //   if (effreco["kMpl02"][iE].sel == cut){ out["kMpl02"].push_back(effreco["kMpl02"][iE]); }
  // }

  return out;
}

//-----------------------------------------------------------------------------------
std::vector<eff_reco> ave_eff(std::map<std::string, std::vector<eff_reco> > effreco, std::string insigname){

  std::vector<eff_reco> out;
  out.clear();

  std::map<std::string , int> denompermasspoint;
  std::map<std::string , double> numeratorpermasspoint;
  //Initialization
  for (auto iE : effreco){
    for (auto tm : iE.second){
    denompermasspoint[tm.M_bins] = 0;
    numeratorpermasspoint[tm.M_bins] = 0.;
    }
  }

  //Couplings for the moment
  std::vector<std::string> coups; coups.clear();
  if (insigname == "grav"){
    coups.push_back("kMpl001");
    coups.push_back("kMpl01");
    coups.push_back("kMpl02");
  } else if (insigname == "heavyhiggs"){
    coups.push_back("0p014");
    coups.push_back("1p4");
    coups.push_back("5p6");
  } else {
    std::cout << "The signal you are looking for doesn't exist. " << std::endl;
    exit(1);
  }


  for (auto iE : effreco){
    for (auto tm : iE.second){
      ++denompermasspoint[tm.M_bins];
      numeratorpermasspoint[tm.M_bins] += tm.efforacc;
    }
  }

  //We will choose the one with the biggest number of samples
  eff_reco tmp;
  if (insigname == "grav"){
    for (auto iE : effreco["kMpl01"]){
      tmp.M_bins = iE.M_bins;
      tmp.sel    = iE.sel;
      tmp.cat    = iE.cat;
      tmp.efforacc = numeratorpermasspoint[iE.M_bins]/ ( (double) denompermasspoint[iE.M_bins]);

      out.push_back(tmp);

    }
  } else if (insigname== "heavyhiggs"){
    for (auto iE : effreco["1p4"]){
      tmp.M_bins = iE.M_bins;
      tmp.sel    = iE.sel;
      tmp.cat    = iE.cat;
      tmp.efforacc = numeratorpermasspoint[iE.M_bins]/ ( (double) denompermasspoint[iE.M_bins]);

      out.push_back(tmp);
    }
  }//end of if

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

void writeToJson(std::map<std::string, std::vector<eff_reco> > valuestowrite, std::string outputfile)
{
  // Short alias for this namespace
  namespace pt = boost::property_tree;

  // Create a root
  pt::ptree root;

  std::string sel_node;

  for (auto &valmap : valuestowrite){
    // pt::ptree values_node;
    std::string samname;
    for (auto &val : valmap.second){
      samname = val.samplename;
      sel_node = val.sel;

      root.put(Form("%s.%s.kMpl",samname.c_str(), sel_node.c_str() ), valmap.first);
      root.put(Form("%s.%s.M_bins",samname.c_str(), sel_node.c_str() ), val.M_bins);
      root.put(Form("%s.%s.sel",samname.c_str(), sel_node.c_str() ), val.sel);
      root.put(Form("%s.%s.cat",samname.c_str(), sel_node.c_str() ), val.cat);
      root.put(Form("%s.%s.efforacc",samname.c_str(), sel_node.c_str() ), std::to_string(val.efforacc));

    }
  }
  pt::write_json(std::cout, root);
  pt::write_json(outputfile, root);


}

std::map<std::string, std::vector<eff_reco> > readJson(std::string inputfile){

  // Short alias for this namespace
  namespace pt = boost::property_tree;
  using boost::property_tree::ptree;

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(inputfile, root);

  std::map<std::string, std::vector<eff_reco> > theeffreco;

  std::string samname;
  std::string sel_node;

  for (ptree::const_iterator it = root.begin(); it != root.end(); ++it) {
    eff_reco tmpeffreco;

    samname = it->first;

    for (ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {

      sel_node = jt->first;

      tmpeffreco.kMpl   = root.get<std::string>( Form("%s.%s.kMpl",samname.c_str(), sel_node.c_str() ) ) ;
      tmpeffreco.M_bins = root.get<std::string>( Form("%s.%s.M_bins",samname.c_str(), sel_node.c_str() ) ) ;
      tmpeffreco.sel    = root.get<std::string>( Form("%s.%s.sel",samname.c_str(), sel_node.c_str() ) ) ;
      tmpeffreco.cat    = root.get<std::string>( Form("%s.%s.cat",samname.c_str(), sel_node.c_str() ) ) ;
      tmpeffreco.efforacc = std::stod( root.get<std::string>( Form("%s.%s.efforacc",samname.c_str(), sel_node.c_str()) ) );

      theeffreco[ tmpeffreco.kMpl ].push_back( tmpeffreco );

      std::cout << tmpeffreco.kMpl << " " << tmpeffreco.M_bins << " " << tmpeffreco.sel << " "  << tmpeffreco.cat <<  " "  << tmpeffreco.efforacc << std::endl;
    }
  }

  return theeffreco;

}

std::map<std::string , std::vector<xsec> > loadXsections(const std::string & year, bool basedonCoup)
{

  std::map<std::string, std::vector<xsec> >thexsections;

  std::vector<std::string> samples = getSampleList();

  xsec tmpxsec;
  std::string coup;

  for(auto isample : samples) {
    //Run a single year each time
    if ( isample.find(year) == std::string::npos ) continue;
    //Run only on RS samples for now.
    // if( isample.find("RSGravitonToGammaGamma") == std::string::npos ) continue;

    double xsecvalue = 0.;
    if (year == "2016"){
      xsecvalue = ExoDiPhotons::crossSection(getSampleBase(isample,year)+"_"+year);
      std::cout << "!!!!! 2016 xsec value!!!!!"<< std::endl;
    } else{
      xsecvalue = ExoDiPhotons::crossSection(getSampleBase(isample,year));
    }
    std::cout << "xsec " << xsecvalue << std::endl;

    tmpxsec.name = getSampleBase(isample,year);
    tmpxsec.val = xsecvalue;
    tmpxsec.error = 0.;

    // std::cout << "xsec " << xsecvalue << std::endl;
    coup = get_str_between_two_str(getBase(isample), "kMpl", "_M_");

    if (isample.find("RSGraviton") != std::string::npos) {
      coup = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    } else if ( isample.find("GluGluSpin0") != std::string::npos ){
      coup = get_str_between_two_str(getBase(isample), "GluGluSpin0ToGammaGamma_W_", "_M_");
    }

    if ( coup == "001" ){ coup = "kMpl001";}
    else if ( coup == "01" ){ coup = "kMpl01";}
    else if ( coup == "02" ){ coup = "kMpl02";}
    //GluGlu coups are ok
    else if ( coup == "0p014" || coup == "1p4" || coup == "5p6"){ std::cout << "Coups for GluGlu ok "<<std::endl; }
    else {
      std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
      exit(1);
    }

    tmpxsec.M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");

    if (basedonCoup) { thexsections[coup].push_back(tmpxsec); }
    else {

      if ( coup == "kMpl001" || coup == "0p014" ){ tmpxsec.coup =  1.4 * pow(10.,-4); }
      else if ( coup == "kMpl01" || coup == "1p4" ) { tmpxsec.coup =  1.4 * pow(10.,-2); }
      else if ( coup == "kMpl02" || coup == "5p6" ) { tmpxsec.coup =  5.6 * pow(10.,-2); }
      else {
	std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
	exit(1);
      }
      thexsections[tmpxsec.M_bins].push_back(tmpxsec);

    }

  }

    return thexsections;

}

TGraphErrors * getgraph(std::vector<double> xvalues, std::vector<double> yvalues, std::vector<double> xvaluesErr, std::vector<double> yvaluesErr){

  TGraphErrors * thegraph = new TGraphErrors( yvalues.size(), &xvalues[0], &yvalues[0], &xvaluesErr[0], &yvaluesErr[0] );
  return thegraph;
}
