#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
#include "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h"

//RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooStats/HLFactory.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "RooStats/RooStatsUtils.h"

//ROOT
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

using namespace RooFit;
using namespace RooStats;

//-----------------------------------------------------------------------------------
static const Int_t NCAT = 2; //BB and BE for the moment 
Int_t MINmass= 310;
Int_t MINmassBE = 340;
Int_t MAXmass= 1000;
Int_t nBinsMass=120;
//-----------------------------------------------------------------------------------
//Declarations here definition after main
RooRealVar* buildRooVar(std::string name, std::string title, int nBins, double xMin, double xMax, std::string unit);
void AddBkgData(RooWorkspace* w, const std::string &isample, const std::string &year, const std::string &ws_dir);
void PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir);
std::vector<RooFitResult*> BkgModelFitDiJetFunc(RooWorkspace* w, bool blind, bool dobands, const std::string &year, const std::string &ws_dir);
std::vector<RooFitResult*> BkgModelFitPowFunc(RooWorkspace* w, bool blind, bool dobands, const std::string &year, const std::string &ws_dir);
std::vector<RooFitResult*> BkgModelFitExpPARFunc(RooWorkspace* w);
void runAllFits(const std::string &year, const std::string &ws_dir);
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample);
std::string getBase(const std::string & sampleName);
void SetConstantParams(const RooArgSet* params) ;
TPaveText* get_labelsqrt( int legendquadrant );
TPaveText* get_labelcms( int legendquadrant, std::string year, bool sim);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir;

  if(argc!=3) {
    std::cout << "Syntax: bkgFit.exe [2016/2017/2018] [input]" << std::endl;
      return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    inputdir = argv[2];
 }

  //========================================================================
  //read the reduced input root trees
  initForFit(inputdir);
  //========================================================================
  //Run the fits
  runAllFits(year,"output/bkg");

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void runAllFits(const std::string &year, const std::string &ws_dir){

  std::vector<std::string> samples = getSampleListForFit();

  for(auto isample : samples) {
    //We will process one year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    //We will only process data since bkg shape is data-driven
    if ( isample.find("data") == std::string::npos ) continue; 
    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    runfits(year, ws_dir, isample);
  }

}
//-----------------------------------------------------------------------------------
void AddBkgData(RooWorkspace* w, const std::string &isample, const std::string &year, const std::string &ws_dir) {

  Int_t ncat = NCAT;//BB and BE for the moment

  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;

  RooRealVar* mgg = buildRooVar("mgg", "M(gg)", 250, MINmass, MAXmass,"GeV");
  // RooRealVar* weight = buildRooVar("weight", "event weight", 0, 0 , 1000 , "nounits");
  RooRealVar* eventClass = buildRooVar("eventClass", "eventClass", 0, -10, 10, "nounits");
  
  RooArgSet* rooVars = new RooArgSet(*mgg, *eventClass); 

  //Cut on the region
  TString mainCut = TString::Format("mgg>=(%.1f) && mgg<=(%.1f)", minMassFit, maxMassFit);

  // Create dataset
  RooDataSet Data("Data","dataset",treesforfit[getBase(isample)], *rooVars, mainCut, "weight"); 
  std::cout << "Data reading to RooDataSet" << std::endl;
  Data.Print("v");
  std::cout << "---- nX:  " << Data.sumEntries() << std::endl;

  // split into NCAT categories
  RooDataSet* dataToFit[NCAT];  
  for (int c=0; c<ncat; ++c) {
    // int theCat = c+1;

    if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==0"));
    if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==1"));

    std::cout << "Data for category = " << c << std::endl;
    dataToFit[c]->Print("v");
    std::cout << "---- nX:  " << dataToFit[c]->sumEntries() << std::endl;

    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }

  std::cout << "data, no split" << std::endl;
  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut);
  w->import(*data, Rename("Data"));
    
  data->Print("v");
  std::cout << "---- nX:  " << data->sumEntries() << std::endl; 

  RooPlot* plotPhotonsMassBkgOnlyData[NCAT];
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Only Data",0,0,500,500);
  TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661,"brNDC");

  for (int c=0; c<ncat; ++c) {
    plotPhotonsMassBkgOnlyData[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);
    dataToFit[c]->plotOn(plotPhotonsMassBkgOnlyData[c], MarkerColor(c+2), LineColor(c+2) );  
    if (c==0){plotPhotonsMassBkgOnlyData[c]->Draw();}
    else if (c==1){plotPhotonsMassBkgOnlyData[c]->Draw("same");}
    legdata->AddEntry(plotPhotonsMassBkgOnlyData[c]->getObject(0),Form("Data_cat%d",c),"LPE");

  }  
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  // legdata->SetTextAlign(31);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");

  ctmp->SetLogy();
  ctmp->SaveAs( Form("%s/BkgOnlyData_%s_log.png", ws_dir.c_str(),year.c_str()) );
  ctmp->SetLogy(0);
  ctmp->SaveAs( Form("%s/BkgOnlyData_%s.png", ws_dir.c_str(),year.c_str()) );

  
  
}

//-----------------------------------------------------------------------------------
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample){

  std::cout << ws_dir << std::endl;

  TFile *fout = new TFile(Form("%s/bkg_%s.root", ws_dir.c_str(),year.c_str()), "recreate");

  TString fileBaseName("HighMassGG");    
  TString fileBkgName("HighMassGG.inputbkg");
  TString card_name("HighMass-hgg_models_Bkg.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  // range for the variables
  w->var("mgg")->setMin(MINmass);
  w->var("mgg")->setMax(MAXmass);
  w->Print("V");

  //========================================================================
  std::cout << "Adding bkg data" <<std::endl;
  AddBkgData(w, isample, year, ws_dir);
  //========================================================================
  
  //blind or not
  bool blind = true;
  //bands
  bool dobands = true;
  std::vector<RooFitResult*> fitresults_dijet = BkgModelFitDiJetFunc(w, blind, dobands, year, ws_dir); 
  //std::vector<RooFitResult*> fitresults_pow = BkgModelFitPowFunc(w, blind, year, ws_dir); 
  //std::vector<RooFitResult*> fitresults_ExpPAR = BkgModelFitExpPARFunc(w); 


  // TCanvas* canv1 = new TCanvas(Form("Canvas1_%s",isample.c_str()),Form("Canvas_%s",isample.c_str()));
  // canv1->cd();
  // RooPlot* p = w->var("mgg")->frame(RooFit::Range(230., 10000.), RooFit::Bins( 2385));

  // //========================================================================
  // std::cout << "Fit the response function" <<std::endl;
  // sigModelResponseFcnFit(w, stof(M_bins), coupling, year);
  // //========================================================================
  // std::cout << "Fit the gen function" <<std::endl;
  // sigModelGenFcnFit(w, stof(M_bins), coupling, year);
  //========================================================================
  //Save the resulting workspace
  // w->Print("V");
  fout->cd();
  w->Write();
  fout->Write();
  fout->Close();


}

std::vector<RooFitResult*> BkgModelFitDiJetFunc(RooWorkspace* w, bool blind, bool dobands, const std::string &year, const std::string &ws_dir) {

  Int_t ncat = NCAT;
  RooDataSet* data[NCAT];
  std::vector<RooFitResult*> fitresult;

  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
  Float_t minMassFit, maxMassFit;

  for (int c = 0; c < ncat; ++c) {
    if (c==0){ 
      minMassFit = MINmass;
      maxMassFit = MAXmass;
    } else if (c==1){
      minMassFit = MINmassBE;
      maxMassFit = MAXmass;    
    }

    
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_DiJet_cat%d",c), "TMath::Max(1e-50,pow(@0,@1+@2*log(@0)))" , RooArgList(*mgg, *w->var(TString::Format("PhotonsMass_bkg_dijet_linc_cat%d",c)), *w->var(TString::Format("PhotonsMass_bkg_dijet_logc_cat%d",c))) );
  
    //RooPowLogPdf *PhotonsMassBkgTmp0 = new RooPowLogPdf(TString::Format("PhotonsMassBkg_DiJet_cat%d",c), TString::Format("PhotonsMassBkg_DiJet_cat%d",c),  *mgg, *w->var(TString::Format("PhotonsMass_bkg_dijet_linc_cat%d",c)), *w->var(TString::Format("PhotonsMass_bkg_dijet_logc_cat%d",c))) ;

    fitresult.push_back( (RooFitResult* ) PhotonsMassBkgTmp0->fitTo(*data[c], RooFit::Minimizer("Minuit2"), RooFit::PrintLevel(2),SumW2Error(kTRUE), Range(minMassFit,maxMassFit), RooFit::Save(kTRUE)) );
    w->import(*PhotonsMassBkgTmp0);

    std::cout << TString::Format("******************************** Background Fit results DIJET cat %d***********************************", c) << std::endl;
    fitresult[c]->Print("V");
    
    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    PlotFitResult(w, ctmp, c, mgg, data[c], PhotonsMassBkgTmp0, minMassFit, maxMassFit, blind, dobands, 3, year, ws_dir);
    ctmp->SaveAs( Form("%s/Bkg_cat%d_DIJET_%s.png", ws_dir.c_str(),c,year.c_str()) );

  }
  
  return fitresult;

}


void PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir){

  RooPlot* plotPhotonsMassBkg[NCAT];

  plotPhotonsMassBkg[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);
  
  data->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
  PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit));
  //RooPlot::chiSquare() should give you the chi2/ndof.
  double chi2 = plotPhotonsMassBkg[c]->chiSquare(numoffittedparams);
  Int_t ndof = nBinsMass-numoffittedparams;
  std::cout<<"------> ndof"<< ndof<<std::endl;
  //https://root.cern.ch/root/html524/TMath.html#TMath:Prob   
  double prob = TMath::Prob(chi2*ndof, ndof);
  std::cout << prob << std::endl;

  if( blind ) {
    RooDataSet* data_blind = (RooDataSet*) data->reduce(*w->var("mgg"),"mgg < 1000");
    data_blind->plotOn(plotPhotonsMassBkg[c]); 
  } else {
    data->plotOn(plotPhotonsMassBkg[c]);    
  } 
  
  plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
  plotPhotonsMassBkg[c]->Draw();  

  TLegend *legdata = new TLegend(0.37,0.67,0.62,0.82, TString::Format("Category %d",c), "brNDC");
  legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
  legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: DiJet","L");
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  // legdata->SetTextAlign(31);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");

  TPaveText* label_cms = get_labelcms(0, year, false);
  TPaveText* label_sqrt = get_labelsqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  
  //write down the chi2 of the fit on the
 
  TPaveText* label_chi2 = new TPaveText(0.55,0.33,0.79,0.44, "brNDC");
  label_chi2->SetFillColor(kWhite);
  label_chi2->SetTextSize(0.035);
  label_chi2->SetTextFont(42);
  label_chi2->SetTextAlign(31); // align right
  label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
  label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
  label_chi2->Draw("same");

  //********************************************************************************//
  if (dobands) {

    RooAbsPdf *cpdf = PhotonsMassBkgTmp0;
    TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
    RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
    nlim->removeRange();
      
    RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
    for (int i=1; i<(plotPhotonsMassBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
      double lowedge = plotPhotonsMassBkg[c]->GetXaxis()->GetBinLowEdge(i);
      double upedge  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinUpEdge(i);
      double center  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinCenter(i);
	
      double nombkg = nomcurve->interpolate(center);
      nlim->setVal(nombkg);
      mgg->setRange("errRange",lowedge,upedge);
      RooAbsPdf *epdf = 0;
      epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
      RooAbsReal *nll = epdf->createNLL(*(data),Extended());
      RooMinimizer minim(*nll);
      minim.setStrategy(0);
      // double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
      double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
      minim.migrad();
      minim.minos(*nlim);
      printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
      onesigma->SetPoint(i-1,center,nombkg);
      onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
      minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
      // eventually if cl = 0.95 this is the usual 1.92!      
	
      minim.migrad();
      minim.minos(*nlim);
	
      twosigma->SetPoint(i-1,center,nombkg);
      twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
      delete nll;
      delete epdf;
    }

    mgg->setRange("errRange",minMassFit,maxMassFit);
      
    twosigma->SetLineColor(kGreen);
    twosigma->SetFillColor(kGreen);
    twosigma->SetMarkerColor(kGreen);
    twosigma->Draw("L3 SAME");
      
    onesigma->SetLineColor(kYellow);
    onesigma->SetFillColor(kYellow);
    onesigma->SetMarkerColor(kYellow);
    onesigma->Draw("L3 SAME");
      
    legdata->AddEntry(onesigma,"#pm1#sigma","F");
    legdata->AddEntry(twosigma,"#pm2#sigma","F");
    legdata->Draw("same");
    
    plotPhotonsMassBkg[c]->Draw("SAME"); 
  }
  
}

std::vector<RooFitResult*> BkgModelFitPowFunc(RooWorkspace* w, bool blind, bool dobands, const std::string &year, const std::string &ws_dir) {

  Int_t ncat = NCAT;
  RooDataSet* data[NCAT];
  std::vector<RooFitResult*> fitresult;

  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;


  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Pow_cat%d",c), "TMath::Max(1e-50,1.+ @0*@1 + pow(@0,2.)*@2 + pow(@0,3.)*@3)", RooArgList(*mgg, *w->var(TString::Format("PhotonsMass_bkg_pow_a0_cat%d",c)), *w->var(TString::Format("PhotonsMass_bkg_pow_a1_cat%d",c)), *w->var(TString::Format("PhotonsMass_bkg_pow_a2_cat%d",c)),*w->var(TString::Format("PhotonsMass_bkg_pow_a3_cat%d",c)) ) );
  
    fitresult.push_back( (RooFitResult* ) PhotonsMassBkgTmp0->fitTo(*data[c], RooFit::Minimizer("Minuit2"), RooFit::PrintLevel(2),SumW2Error(kTRUE), Range(minMassFit,maxMassFit), RooFit::Save(kTRUE)) );
    w->import(*PhotonsMassBkgTmp0);

    std::cout << TString::Format("******************************** Background Fit results POW cat %d***********************************", c) << std::endl;
    fitresult[c]->Print("V");

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    PlotFitResult(w, ctmp, c, mgg, data[c], PhotonsMassBkgTmp0, minMassFit, maxMassFit, blind, dobands, 1, year, ws_dir);
    ctmp->SaveAs( Form("%s/Bkg_cat%d_POW_%s.png", ws_dir.c_str(),c,year.c_str()) );
    
  }
  
  return fitresult;

}

std::vector<RooFitResult*> BkgModelFitExpPARFunc(RooWorkspace* w) {

  Int_t ncat = NCAT;
  RooDataSet* data[NCAT];
  std::vector<RooFitResult*> fitresult;

  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;
  
  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");
 
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    mgg->setRange("bkg range", MINmass, MAXmass);
  
    // fit con expol 
    RooFormulaVar *p1mod= new RooFormulaVar(TString::Format("par1ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_ExpPAR1_cat%d",c)));
    RooFormulaVar *p2mod= new RooFormulaVar(TString::Format("par2ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_ExpPAR2_cat%d",c)));

    RooAbsPdf* PhotonsMassBkg = new RooGenericPdf(TString::Format("PhotonsMassBkg_Lau_o2_cat%d",c), "exp(-@1*@0)*pow(@0, @2)", RooArgList(*mgg, *p1mod, *p2mod));
    
    fitresult.push_back( PhotonsMassBkg->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE)) );   
    std::cout<<TString::Format("******************************** Background Fit results cat %d ***********************************", c)<<std::endl;
    fitresult[c]->Print("V");

  }

  return fitresult;

}

//-----------------------------------------------------------------------------------
RooRealVar* buildRooVar(std::string name, std::string title, int nBins, double xMin, double xMax, std::string unit){

  RooRealVar* rooVar;
  
  if ( unit == "nounits"){
    rooVar =  new RooRealVar( name.c_str() , title.c_str(), xMin, xMax,  "");
  } else{
    rooVar =  new RooRealVar( name.c_str() , title.c_str(), xMin, xMax, unit.c_str() );

    rooVar->setMin(xMin); 
    rooVar->setMax(xMax);
    rooVar->setBins(nBins);
  }
  return rooVar;

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
void SetConstantParams(const RooArgSet* params) {
  
  std::cout << std::endl; std::cout << "Entering SetConstantParams" << std::endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}
//-----------------------------------------------------------------------------------
TPaveText* get_labelcms( int legendquadrant = 0 , std::string year="2016", bool sim=false) {

  if( legendquadrant!=0 && legendquadrant!=1 && legendquadrant!=2 && legendquadrant!=3 ) {
    std::cout << "warning! legend quadrant '" << legendquadrant << "' not yet implemented for cms label. using 2." << std::endl;
    legendquadrant = 2;
  }

  float x1 = 0.; float y1 = 0.; float x2 = 0.; float y2 = 0.;
  if( legendquadrant==1 ) {
    x1 = 0.63;
    y1 = 0.83;
    x2 = 0.8;
    y2 = 0.87;
  } else if( legendquadrant==0 ) {
    x1 = 0.175;
    y1 = 0.953;
    x2 = 0.6;
    y2 = 0.975;

  } else if( legendquadrant==3 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
  }
 
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brndc" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if( legendquadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
  
  std::string lefttext;
   
     
  if (sim)  lefttext = "cms simulation"; 
  else {
    lefttext = "cms preliminary, 19.5 fb^{-1}";
  }
  cmslabel->AddText(lefttext.c_str());
  return cmslabel;

}
//-----------------------------------------------------------------------------------
TPaveText* get_labelsqrt( int legendquadrant ) {

  if( legendquadrant!=0 && legendquadrant!=1 && legendquadrant!=2 && legendquadrant!=3 ) {
    std::cout << "warning! legend quadrant '" << legendquadrant << "' not yet implemented for sqrt label. using 2." << std::endl;
    legendquadrant = 2;
  }


  float x1 = 0.; float y1 = 0.; float x2 = 0.; float y2 = 0.;
  if( legendquadrant==1 ) {
    x1 = 0.63;
    y1 = 0.78;
    x2 = 0.8;
    y2 = 0.82;
  } else if( legendquadrant==2 ) {
    x1 = 0.25;
    y1 = 0.78;
    x2 = 0.42;
    y2 = 0.82;
  } else if( legendquadrant==3 ) {
    x1 = 0.25;
    y1 = 0.16;
    x2 = 0.42;
    y2 = 0.2;
  } else if( legendquadrant==0 ) {
    x1 = 0.65;
    y1 = 0.953;
    x2 = 0.87;
    y2 = 0.975;
  }


  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "brndc");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); 
  label_sqrt->AddText("#sqrt{s} = 8 tev");
  return label_sqrt;

}

