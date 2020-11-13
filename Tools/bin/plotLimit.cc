#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>

#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
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
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TF1.h"

#include "RooRealVar.h"

//-----------------------------------------------------------------------------------
struct xsec{
  
  std::string name;
  double coup;
  std::string M_bins;
  double val;
  double error ;
  
};

//-----------------------------------------------------------------------------------
//Declarations here definition after main
std::map<std::string, std::vector<xsec> > loadXsections(const std::string &year, bool basedonCoup, std::string signame, bool use_fb);
TGraphErrors* getXsecGraph( std::vector<xsec> xsections, bool basedonCoup);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);
RooSpline1D* graphToSpline(std::string name, TGraphErrors *graph, RooRealVar* MH, double xmin, double upperxmax);
TGraphErrors* SplineTograph(std::string name, RooSpline1D* thespline, RooRealVar* MH, double xmin, double upperxmax);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string input, coupling, year, outputdir, signame;
  if(argc!=6) {
    std::cout << "Syntax: plotLimit.exe [inputfile] [coupling] [year] [outputdir] [signalname]" << std::endl;
    return -1;
  }
  else {
    input = argv[1];
    coupling = argv[2];
    year = argv[3];
    outputdir = argv[4];
    signame = argv[5];
  }

  //========================================================================
  // include signal samples but not unskimmed data samples
  init(false, true);

  // gROOT->Reset();
  //  gROOT->ProcessLine(".L CMSstyle.C");
  //  CMSstyle();

  //  gROOT->ProcessLine(".L CMS1bACStyle.C");
  //  CMS1bACStyle();

  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(0);
  gSystem->Exec(Form("wc -l %s | awk '{print $1}' > pptt", input.c_str() ));
  std::ifstream numPoin("pptt");
  int NumPunkte = 0;//221;
  numPoin >> NumPunkte;
  std::cout<< NumPunkte << std::endl;
  std::vector<double> mass, obs, exp, expP1s, expP2s, expM1s, expM2s; 
  mass.resize(NumPunkte);
  obs.resize(NumPunkte);
  exp.resize(NumPunkte);
  expP1s.resize(NumPunkte);
  expP2s.resize(NumPunkte);
  expM1s.resize(NumPunkte);
  expM2s.resize(NumPunkte);
  
  // double mass[3]=  {0,0,0};
  // double obs[3]=   {0,0,0};
  // double exp[3]=   {0,0,0};
  // double expP1s[3]={0,0,0};
  // double expP2s[3]={0,0,0};
  // double expM1s[3]={0,0,0};
  // double expM2s[3]={0,0,0};
  
  std::ifstream file(input);

  TGraphErrors* expGraph_init = new TGraphErrors();
  TGraphErrors* exp1SGraph_init = new TGraphErrors();
  TGraphErrors* exp2SGraph_init = new TGraphErrors();
  TGraphErrors* obsGraph_init = new TGraphErrors();
  // TGraph* theoryGraph = new TGraph();
  // TGraph* obsCMSGraph = new TGraph();

  TGraphAsymmErrors * graphmede = new TGraphAsymmErrors();
  TGraphAsymmErrors * graph68up = new TGraphAsymmErrors();
  TGraphAsymmErrors * graph68dn = new TGraphAsymmErrors();
  TGraphAsymmErrors * graph95up = new TGraphAsymmErrors();
  TGraphAsymmErrors * graph95dn = new TGraphAsymmErrors();

  bool use_fb = true;
  std::map<std::string, std::vector<xsec> > xsections = loadXsections(year, true, signame, use_fb);
  std::map<std::string, TGraphErrors*> grxs;
  std::map<std::string, TGraphErrors*> grxs_spline;
  grxs[coupling] = getXsecGraph(xsections[coupling], true);
  // std::map<std::string, RooSpline1D *> xsSplines;
  RooRealVar* MH = new RooRealVar("MH", "MH", 1000.);
  MH->setConstant();
  double upperxmax = 9000.;
  double xmin = 230.;
  RooSpline1D *xsSpline = graphToSpline(Form("fxs_%s",coupling.c_str()), grxs[coupling], MH , xmin, upperxmax);
  std::map<std::string, TF1 *> fm;
  fm[coupling] = new TF1(TString::Format("fm_%s",coupling.c_str()), "pol6", xmin, upperxmax);
  grxs[coupling]->Fit(TString::Format("fm_%s",coupling.c_str()), "R");
  grxs_spline[coupling] = SplineTograph(Form("fxs_spline_%s",coupling.c_str()), xsSpline, MH , xmin, upperxmax);
  
  //Search for heavy, top-like quark pair production in the dilepton final state in pp collisions at sqrt(s) = 7 TeV CMS Collaboration
  // Submitted on 24 Mar 2012
  // CMS-EXO-11-050, CERN-PH-EP-2012-081
  // arXiv:1203.5410v1 [hep-ex]
  
  // theoryGraph->SetPoint(0,350,5.297);
  // theoryGraph->SetPoint(1,400,2.386);
  // theoryGraph->SetPoint(2,450,1.153);
  // theoryGraph->SetPoint(3,500,0.590);
  // theoryGraph->SetPoint(4,550,0.315);  
  // theoryGraph->SetPoint(5,600,0.174);  
  // theoryGraph->SetPoint(6,650,0.0999); 
  // theoryGraph->SetPoint(7,700,0.0585); 
  // theoryGraph->SetPoint(8,750,0.0350);   
  // theoryGraph->SetPoint(9,800,0.0213);   
  // theoryGraph->SetPoint(10,850,0.0132);   
  // theoryGraph->SetPoint(11,900,0.00828);
  // theoryGraph->SetPoint(12,950,0.00525);
  // theoryGraph->SetPoint(13,1000,0.00336);
  // theoryGraph->SetPoint(14,1400,0.000114);
  // theoryGraph->SetPoint(15,1500,0.0000499);

 
  // obsCMSGraph->SetPoint(0,350,0.47);
  // obsCMSGraph->SetPoint(1,400,0.26);
  // obsCMSGraph->SetPoint(2,450,0.22);
  // obsCMSGraph->SetPoint(3,500,0.18);
  // obsCMSGraph->SetPoint(4,550,0.16);  
  // obsCMSGraph->SetPoint(5,600,0.14);  


  for( int i = 0 ; i < NumPunkte ; ++i ) {
    file >> mass[i] >> obs[i] >> expM2s[i] >> expM1s[i] >> exp[i] >> expP1s[i] >> expP2s[i];
    expGraph_init->SetPoint(i,mass[i],exp[i]);
    obsGraph_init->SetPoint(i,mass[i],obs[i]);
    std::cout << mass[i]<< " " << obs[i]<< " " << expM2s[i]<< " " << expM1s[i]<< " " << exp[i]<< " " << expP1s[i]<< " " << expP2s[i]<<std::endl;
    
    exp1SGraph_init->SetPoint(i,mass[i],expM1s[i]);
    exp2SGraph_init->SetPoint(i,mass[i],expM2s[i]);

    graphmede->SetPoint(i,mass[i],exp[i]);
    graph68up->SetPoint(i,mass[i],expP1s[i]);
    graph68dn->SetPoint(i,mass[i],expM1s[i]);
    graph95up->SetPoint(i,mass[i],expP2s[i]);
    graph95dn->SetPoint(i,mass[i],expM2s[i]);
    
  }    
  // TGraph* expGraph = new TGraph(mass.size(), &mass[0], &exp);
  // TGraph* obsGraph = new TGraph(mass.size(), &mass[0],&obs);
  // TGraph* exp1SGraph = new TGraph(mass.size(), &mass[0],&expM1s);
  // TGraph* exp2SGraph = new TGraph(mass.size(), &mass[0],&expM2s);

   for( int i = 0 ; i < NumPunkte ; ++i ) {	
     exp1SGraph_init->SetPoint(NumPunkte+i,mass[NumPunkte-1-i],expP1s[NumPunkte-1-i]);
     exp2SGraph_init->SetPoint(NumPunkte+i,mass[NumPunkte-1-i],expP2s[NumPunkte-1-i]);
   }

   //To splines and back: This is to get a better interpolation which we do not need for the moment. 
   // RooSpline1D *expGraphSpline = graphToSpline(Form("fexpgr_%s",coupling.c_str()), expGraph_init, MH , xmin, upperxmax);
   // TGraphErrors* expGraph = SplineTograph(Form("fexpgr_spline_%s",coupling.c_str()), expGraphSpline, MH , xmin, upperxmax);
   // RooSpline1D *exp1SGraphSpline = graphToSpline(Form("fexpgr_%s",coupling.c_str()), exp1SGraph_init, MH , xmin, upperxmax);
   // TGraphErrors* exp1SGraph = SplineTograph(Form("fexpgr_spline_%s",coupling.c_str()), exp1SGraphSpline, MH , xmin, upperxmax);
   // RooSpline1D *exp2SGraphSpline = graphToSpline(Form("fexpgr_%s",coupling.c_str()), exp2SGraph_init, MH , xmin, upperxmax);
   // TGraphErrors* exp2SGraph = SplineTograph(Form("fexpgr_spline_%s",coupling.c_str()), exp2SGraphSpline, MH , xmin, upperxmax);
   // RooSpline1D *obsGraphSpline = graphToSpline(Form("fexpgr_%s",coupling.c_str()), obsGraph_init, MH , xmin, upperxmax);
   // TGraphErrors* obsGraph = SplineTograph(Form("fexpgr_spline_%s",coupling.c_str()), obsGraphSpline, MH , xmin, upperxmax);

   TGraphErrors* expGraph = (TGraphErrors*) expGraph_init->Clone();
   TGraphErrors* exp1SGraph = (TGraphErrors*) exp1SGraph_init->Clone();
   TGraphErrors* exp2SGraph = (TGraphErrors*) exp2SGraph_init->Clone();
   TGraphErrors* obsGraph = (TGraphErrors*) obsGraph_init->Clone();


   
   //smooth
   TString fitstring = "[0] + [1]*x*x + [2]*x*x*x +[3]*x*x*x*x + [4]*x";
   TF1 *medfunc  = new TF1("medfunc" , fitstring, 500., 8000.);
   TF1 *up68func = new TF1("up68func", fitstring, 500., 8000.);
   TF1 *dn68func = new TF1("dn68func", fitstring, 500., 8000.);
   TF1 *up95func = new TF1("up95func", fitstring, 500., 8000.);
   TF1 *dn95func = new TF1("dn95func", fitstring, 500., 8000.);

   // expGraph->Fit(medfunc,"R,M,EX0","Q");

   graphmede->Fit(medfunc,"R,M,EX0","Q");
   graph68up->Fit(up68func,"R,M,EX0","Q");
   graph68dn->Fit(dn68func,"R,M,EX0","Q");
   graph95up->Fit(up95func,"R,M,EX0","Q");
   graph95dn->Fit(dn95func,"R,M,EX0","Q");
   
   TCanvas *canv = new TCanvas("canv","Title",800,600); 
   canv->SetLogy();
   canv->SetLogx();
   canv->SetRightMargin(0.08);
   canv->SetLeftMargin(0.15);
  
   exp1SGraph->SetFillColor(kGreen);
   exp2SGraph->SetFillColor(kYellow);
   exp2SGraph->GetXaxis()->SetTitleSize(0.045);
   exp2SGraph->GetYaxis()->SetTitleSize(0.045);

   exp2SGraph->GetYaxis()->SetRangeUser(0.0001,30);
   exp2SGraph->GetYaxis()->SetTitle(signame == "grav" ?  "95% CL limit #sigma(pp#rightarrowG#rightarrow#gamma#gamma) (fb)" : "95% CL limit #sigma(pp#rightarrowS#rightarrow#gamma#gamma) (fb)" );
   exp2SGraph->GetXaxis()->SetTitle(signame == "grav" ? "m_{G} (GeV)" : "m_{S} (GeV)");
   exp2SGraph->GetXaxis()->SetMoreLogLabels();
   exp2SGraph->GetXaxis()->SetRangeUser(500,8000);


   exp2SGraph->Draw("AF");
   exp1SGraph->Draw("F");
   expGraph->SetLineColor(kBlue);
   expGraph->SetLineStyle(7);
   expGraph->SetLineWidth(4);
   expGraph->Draw("L");

   obsGraph->SetMarkerColor(1);
   obsGraph->SetMarkerStyle(20);
   obsGraph->SetMarkerSize(1);
   obsGraph->SetLineWidth(2);
   obsGraph->SetLineColor(kBlack);
   obsGraph->SetLineStyle(1);
   // obsGraph->Smooth();
// obsGraph->Draw("PL");
   // obsGraph->Draw("L");
   //   obsGraph->Draw("LC");

   // graphmede->Draw("same");
  
   // theoryGraph->SetLineWidth(3);
   // theoryGraph->SetLineColor(kRed);
   // theoryGraph->SetLineStyle(4);
   // theoryGraph->Draw("PL");

   // fm[coupling]->SetLineWidth(3);
   // fm[coupling]->SetLineColor(kRed);
   // fm[coupling]->SetLineStyle(4);
   // fm[coupling]->Draw("same");

   grxs_spline[coupling]->SetLineWidth(3);
   grxs_spline[coupling]->SetLineColor(kRed);
   grxs_spline[coupling]->SetLineStyle(4);
   grxs_spline[coupling]->Draw("same");
   
   // obsCMSGraph->SetMarkerColor(kViolet);
   // obsCMSGraph->SetMarkerStyle(20);
   // obsCMSGraph->SetMarkerSize(1);
   // obsCMSGraph->SetLineWidth(2);
   // obsCMSGraph->SetLineColor(kViolet);
   // obsCMSGraph->SetLineStyle(1);
   //obsCMSGraph->Draw("PL");
 
   TLegend *leg = new TLegend(0.2,0.20,0.5,0.50,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->SetTextSize(0.033);
   // leg->AddEntry(theoryGraph,"Theory NNLO","L");

   std::string plabel;
   if ( coupling == "kMpl001" ){ plabel = "#tilde{k}=0.01,J=2"; }
   else if ( coupling == "kMpl01" ){ plabel = "#tilde{k}=0.1,J=2";}
   else if ( coupling == "kMpl02" ){ plabel = "#tilde{k}=0.2,J=2";}
   else if ( coupling == "0p014"){ plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-4},J=0"; }
   else if ( coupling == "1p4"){ plabel = "#frac{#Gamma}{m} = 1.4 #times 10^{-2},J=0"; }
   else if ( coupling == "5p6"){ plabel = "#frac{#Gamma}{m} = 5.6 #times 10^{-2},J=0"; }
   else {
     std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
     exit(1);
   }

   leg->SetHeader(plabel.c_str(),"C");
   leg->AddEntry(expGraph,"expected Limit","L"); //L_{int}=36.4/pb
   leg->AddEntry(exp1SGraph,"#pm1#sigma","F");
   leg->AddEntry(exp2SGraph,"#pm2#sigma","F");
   //   leg->AddEntry(obsGraph,"observed Limit (Asymptotic)","L");
   leg->AddEntry(grxs_spline[coupling],"G_{RS}#rightarrow#gamma#gamma (LO)","l");
     
   leg->Draw();
     
   TLatex* cmsText=new TLatex(0.17,0.90, "CMS");
   cmsText->SetNDC(kTRUE);
   cmsText->SetTextFont(61);
   cmsText->SetLineColor(0);
   cmsText->SetLineStyle(1);
   cmsText->SetLineWidth(1);
   cmsText->SetTextSize(0.035);
   cmsText->Draw();

   TLatex* extraText=new TLatex(0.23,0.90, "Preliminary");
   extraText->SetNDC(kTRUE);
   extraText->SetTextFont(52);
   extraText->SetLineColor(0);
   extraText->SetLineStyle(1);
   extraText->SetLineWidth(1);
   extraText->SetTextSize(0.035);
   extraText->Draw();
   
   std::string thelumi = "";
   if (year == "fullRun2"){
     double fullRun2lumi = luminosity["2016"] + luminosity["2017"] + luminosity["2018"];
     
     std::stringstream stream;
     stream << std::fixed << std::setprecision(3) << fullRun2lumi;
     thelumi = stream.str();
   } else {
     thelumi = std::to_string(luminosity[year]);
   }
   
   TLatex* lumiText=new TLatex(0.70,0.90, Form("%s fb^{-1} (13 TeV)", thelumi.c_str() ) );
   lumiText->SetNDC(kTRUE);
   lumiText->SetTextFont(42);
   lumiText->SetLineColor(0);
   lumiText->SetLineStyle(1);
   lumiText->SetLineWidth(1);
   lumiText->SetTextSize(0.035);
   lumiText->Draw();

   canv->SaveAs( Form("%s/limitplot_%s_%s_%s.png", outputdir.c_str(), signame.c_str(), coupling.c_str(), year.c_str()) );
   

   // TGraph* exp1SGraph = new TGraph();
   // TGraph* exp2SGraph = new TGraph();
   // TGraph* obsGraph = new TGraph();

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
std::map<std::string , std::vector<xsec> > loadXsections(const std::string & year, bool basedonCoup, std::string signame, bool use_fb)
{
  std::map<std::string, std::vector<xsec> >thexsections; 

  std::vector<std::string> samples = getSampleList();

  xsec tmpxsec; 
  std::string coup;

  for(auto isample : samples) {
    //Run a single year each time
    //For the full Run 2, we will set the year to 2017
    //if ( isample.find(2017) == std::string::npos ) continue; 
    if ( isample.find("2017") == std::string::npos ) continue;  

    if (isample.find("RSGraviton") != std::string::npos && signame == "grav") {
      coup = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    } else if ( isample.find("GluGluSpin0") != std::string::npos && signame == "heavyhiggs"){
      coup = get_str_between_two_str(getBase(isample), "GluGluSpin0ToGammaGamma_W_", "_M_");
    } else{
      continue;
    }

    std::cout << coup << std::endl;
    if ( coup == "001" ){ coup = "kMpl001";}
    else if ( coup == "01" ){ coup = "kMpl01";}
    else if ( coup == "02" ){coup = "kMpl02";}
    // else if ( coup == "0p014" ){  }
    // else if ( coup == "1p4"){ }
    // else if ( coup == "5p6"){ }
    else {
      std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
      exit(1);
    }

    std::cout << coup << " " << isample << std::endl;
    double scl = use_fb ? 1000. : 1.;
    tmpxsec.name = getSampleBase(isample,year);
    tmpxsec.val = ExoDiPhotons::crossSection(getSampleBase(isample,year)) * scl;
    tmpxsec.error = 0.;
    if (year == "fullRun2"){ tmpxsec.val = 3. * tmpxsec.val;  }

    std::cout << "xsec " << ExoDiPhotons::crossSection(getSampleBase(isample,year)) << std::endl;
    std::cout << "xsec in fb " << tmpxsec.val << std::endl;
    tmpxsec.M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");
    
    if (basedonCoup) { thexsections[coup].push_back(tmpxsec); }
    else { 
      
      if ( coup == "kMpl001" ){ tmpxsec.coup =  1.4 * pow(10.,-4); }
      else if ( coup == "kMpl01" ) { tmpxsec.coup =  1.4 * pow(10.,-2); }
      else if ( coup == "kMpl02" ) { tmpxsec.coup =  5.6 * pow(10.,-2); }
      else {
	std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
	exit(1);
      }
      thexsections[tmpxsec.M_bins].push_back(tmpxsec); 
      
    }

  }
  
    return thexsections;
    
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
RooSpline1D* graphToSpline(std::string name, TGraphErrors *graph, RooRealVar* MH, double xmin, double upperxmax){

  std::cout << "Graph to Spline: " << name << std::endl; 
  std::vector<double> xValues, yValues;
  for (double mh=xmin; mh<(upperxmax+0.25); mh+=10.){
    xValues.push_back(mh);
    yValues.push_back(graph->Eval(mh));
  }
  RooSpline1D *res = new RooSpline1D(name.c_str(),name.c_str(),*MH,xValues.size(),&(xValues[0]),&(yValues[0]));
  return res;
}

//-----------------------------------------------------------------------------------
TGraphErrors* SplineTograph(std::string name, RooSpline1D* thespline, RooRealVar* MH, double xmin, double upperxmax){

  std::cout << "Spline to Graph: " << name << std::endl; 
  TGraphErrors* graph;

  std::vector<double> xValues, yValues, xValuesErr, yValuesErr;
  
  for (double mh=xmin; mh<(upperxmax+0.25); mh+=10.){
    MH->setVal(mh);
    xValues.push_back(mh);
    yValues.push_back(thespline->getVal());
    xValuesErr.push_back(0.);
    yValuesErr.push_back(0.);
  }

  graph = new TGraphErrors(yValues.size(), &xValues[0], &yValues[0], &xValuesErr[0], &yValuesErr[0]);
  graph->SetName(name.c_str());

  return graph;
}

//-----------------------------------------------------------------------------------
// TGraph * smoothen(TGraph *gr){
//   //hardBound=self.options.smoothen_boundary
//   double relwindow = 5.e-2;
//   Double_t ix, iy, jx, jy;
  
//   for(unsigned int ip =0; ip < gr->GetN(); ip++){

//     gr->GetPoint(ip,ix,iy);
//     double window = relwindow * ix;
//     std::vector<double> ipoints;
//     for(unsigned int jp =0; jp < gr->GetN(); jp++){
//       gr->GetPoint(jp,jx,jy);

//       if (fabs(ix-jx)<window){ ipoints.push_back(jp); }
//     }
//     if (ipoints.size() < 3) { continue; }
//     double minw = gr->GetX()[std::min(ipoints)];
//     double maxw = gr->GetX()[std::max(ipoints)];
//     double left = ix-window/2.;
//     double right = ix+window/2.;
//     if (minw>left) { left = minw; }
//     if (maxw<right) { right = maxw; }
//     TGraph *gr2 = new TGraph();
    
//         map(lambda y: gr2.SetPoint(y[0],gr.GetX()[y[1]],gr.GetY()[y[1]]), enumerate(filter(lambda x:  gr.GetX()[x]>=minw and gr.GetX()[x]<=maxw, ipoints )))
//         ## print "Smoothing ", ix, left, right
//         func = ROOT.TF1("f","[0]*pow(x,[1])")
//         ## gr2.Print()
//         slope = ( log(gr2.GetY()[gr2.GetN()-1])-log(gr2.GetY()[0]) ) / ( log(gr2.GetX()[gr2.GetN()-1])-log(gr2.GetX()[0]) )
//         intercept = pow(gr2.GetY()[0],-slope)
//         gr2.Fit(func,"Q")
//         points.append( (ip, ix, func.Eval(ix)) )

//       }
      
//   for (auto point : points){
//     ip, ix, iy = point
//       scl = iy / gr.GetY()[ip] 
//       gr.SetPoint(ip,ix,iy)
//       gr.SetPointEYlow( ip, gr.GetEYlow()[ip]*scl )
//       gr.SetPointEYhigh( ip, gr.GetEYhigh()[ip]*scl )

//       }
        
// }
