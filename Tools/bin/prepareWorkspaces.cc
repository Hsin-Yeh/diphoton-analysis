#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
// #include "HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"

//RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooStats/HLFactory.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "RooVoigtian.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooExponential.h"

//ROOT
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;
using namespace RooStats;

//-----------------------------------------------------------------------------------
static const Int_t NCAT = 2; //BB and BE for the moment 
Int_t MINmass= 300;
Int_t MAXmass= 10000;
bool gofwithtoys=false;

struct gof {
  double prob;
  double chi2overndof;
};

struct theFitResult {
  RooFitResult* fitres;
  double minNll;
  gof gofresults;
  int minimizestatus;
  int hessestatus;
  int minosstatus;

};

//-----------------------------------------------------------------------------------
//Declarations here definition after main
RooRealVar* buildRooVar(std::string name, std::string title, int nBins, double xMin, double xMax, std::string unit);
void runAllFits(const std::string &year, const std::string &ws_dir, const std::string &checksample);
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample, FILE *resFilegen, FILE *resFileresp);
void AddSigData(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year, const std::string &isample);
theFitResult theFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, double minMassFit, double maxMassFit, std::string thefitrange);
void sigModelResponseFcnFit(RooWorkspace* w,Float_t mass, std::string coupling, const std::string &year, const std::string &ws_dir, const std::string &samplename, FILE *resFileresp);
void sigModelGenFcnFit(RooWorkspace* w,Float_t mass,std::string coupling, const std::string &year, const std::string &ws_dir, const std::string &samplename, FILE *resFilegen);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);
void SetConstantParams(const RooArgSet* params) ;
TPaveText* get_labelsqrt( int legendquadrant );
TPaveText* get_labelcms( int legendquadrant, std::string year, bool sim);
float widthtonum(std::string coupling);
double computeHistFHWM(TH1F*  hist);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir, outputdir, checksample;

  if(argc!=5) {
    std::cout << "Syntax: prepareWorkspaces.exe [2016/2017/2018] [input] [output] [checksample]" << std::endl;
      return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    inputdir = argv[2];
    outputdir = argv[3];
    checksample = argv[4];
 }

  //========================================================================
  //read the reduced input root trees
  initForFit(inputdir);
  //========================================================================
  //Run the fits
  runAllFits(year,outputdir,checksample);

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void runAllFits(const std::string &year, const std::string &ws_dir, const std::string &checksample){

  std::vector<std::string> samples = getSampleListForFit();
  FILE *resFilegen;
  resFilegen = fopen(Form("%s/%s/gen_fits_%s_%s.txt", ws_dir.c_str(), year.c_str(), year.c_str(), checksample.c_str() ),"w");
  // fprintf(resFilegen,"\\hline\n");
  // fprintf(resFilegen,"\\hline\n");
  // fprintf(resFilegen,"Sample & Mass & Width & Category & gof$(\\chi^{2}/ndof)$ & gof(prob) & (minim,hesse,minos) \\\\\n");
  // fprintf(resFilegen,"\\hline\n");
  FILE *resFileresp;
  resFileresp = fopen(Form("%s/%s/resp_fits_%s_%s.txt", ws_dir.c_str(), year.c_str(), year.c_str(), checksample.c_str() ),"w");
  // fprintf(resFileresp,"\\hline\n");
  // fprintf(resFileresp,"\\hline\n");
  // fprintf(resFileresp,"Sample & Mass & Width & Category & gof$(\\chi^{2}/ndof)$ & gof(prob) & (minim,hesse,minos) \\\\\n");
  // fprintf(resFileresp,"\\hline\n");

  for(auto isample : samples) {
    //We will process one year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    if ( isample.find(checksample) == std::string::npos && checksample!="all") continue;    
    // if ( checksample!="all" ) continue;    
    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    runfits(year, ws_dir, isample, resFilegen, resFileresp);
  }

  fclose(resFilegen);
  fclose(resFileresp);


}
//-----------------------------------------------------------------------------------
void AddSigData(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &isample) {

  Int_t ncat = NCAT;//BB and BE for the moment

  RooRealVar* mgg = buildRooVar("mgg", "m_{#gamma#gamma}", 2385, 230., 10000., "GeV");
  RooRealVar* mggGen = buildRooVar("mggGen", "m_{#gamma#gamma} gen", 2385, 230., 10000, "GeV");
  // RooRealVar* weight = buildRooVar("weight", "event weight", 0, 0 , 1000 , "nounits");
  RooRealVar* eventClass = buildRooVar("eventClass", "eventClass", 0, -10, 10, "nounits");
  
  RooArgSet* rooVars = new RooArgSet(*mgg, *mggGen, *eventClass); 

  //Cut on the signal region
  TString CutSignalRegion = TString::Format("mgg>=300&&mgg<=10000&&mggGen>=300&&mggGen<=10000");  
  //Cut on the reduced mass
  TString CutReducedMass = TString::Format("(mgg-mggGen)>-600.&&(mgg-mggGen)<600."); 

  //================================================================
  //Only mgg
  // // treesforfit[getBase(isample)]->Print();
  // RooDataSet* sigWeightedK = new RooDataSet("sigWeightedK","datasetK", treesforfit[getBase(isample)], *rooVars, CutSignalRegion, "weight" );
  // // RooDataSet* sigWeightedK = new RooDataSet("sigWeightedK","datasetK", treesforfit[getBase(isample)], *rooVars, CutSignalRegion, "" );
  // sigWeightedK->Print();

  // RooDataSet* signalK[ncat];
  // RooDataSet* signalAllK;
  
  // for (int c=0; c<ncat; ++c) {
  //   if (c==0) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),CutReducedMass+TString::Format("&& eventClass==0"));
  //   if (c==1) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),CutReducedMass+TString::Format("&& eventClass==1"));

  //   TString myCut;
  //   w->import(*signalK[c],RooFit::Rename(TString::Format("SigWeightK_cat%d",c)));
  //   std::cout << "cat " << c << ", signalK[c]: " << std::endl;
  //   signalK[c]->Print("V");
  //   std::cout << "---- for category " << c << ", nX for signal[c]:  " << signalK[c]->sumEntries() << std::endl; 
  // }
    
  // // Create full weighted signal data set without categorization
  // signalAllK = (RooDataSet*) sigWeightedK->reduce(RooArgList(*w->var("mgg")),CutReducedMass);
  // w->import(*signalAllK, RooFit::Rename("SigWeightK"));
  // std::cout << "now signalAllK" << std::endl;
  // signalAllK->Print("V");
  // std::cout << "---- nX for signalAll:  " << signalAllK->sumEntries() << std::endl; 
  // std::cout << "==================================================================" << std::endl;
  //================================================================

  RooDataSet sigWeighted("sigWeighted","dataset",treesforfit[getBase(isample)], *rooVars, CutSignalRegion, "weight"); 
  // -------------------------
  // reduced mass
  RooFormulaVar *massReduced_formula = new RooFormulaVar("massReduced_formula","","@0-@1",RooArgList(*w->var("mgg"),*w->var("mggGen")));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(-600., 600.);
    
  // split in categories, wrt mgg - this is the dataset to be used for the convolution
  std::cout << std::endl;
  std::cout << "preparing dataset with observable mgg " << std::endl;
  RooDataSet* signal[NCAT];
    
    for (int c=0; c<ncat; ++c) {
      if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),CutReducedMass+TString::Format("&& eventClass==0"));
      if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),CutReducedMass+TString::Format("&& eventClass==1"));
      
      TString myCut;
      w->import(*signal[c],RooFit::Rename(TString::Format("SigWeight_cat%d",c)));
      std::cout << "cat " << c << ", signal[c]: " << std::endl;
      signal[c]->Print("V");
      std::cout << "---- for category " << c << ", nX for signal[c]:  " << signal[c]->sumEntries() << std::endl; 
      std::cout << std::endl;
    }
    
    // Create full weighted signal data set without categorization
    RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(*w->var("mgg"),CutReducedMass);
    w->import(*signalAll, RooFit::Rename("SigWeight"));
    std::cout << "now signalAll" << std::endl;
    signalAll->Print("V");
    std::cout << "---- nX for signalAll:  " << signalAll->sumEntries() << std::endl; 
    std::cout << std::endl;
    



  bool wantResponse = true;
  // -------------------------
  // split in categories, wrt massReduced - to study the detector response
  if (wantResponse) {
    std::cout << std::endl;
    std::cout << "preparing dataset with observable massReduced" << std::endl;
    RooDataSet* signalR[NCAT];
    for (int c=0; c<ncat; ++c) {
      if (c==0) signalR[c] = (RooDataSet*) sigWeighted.reduce(RooArgSet(*w->var("massReduced")),CutReducedMass+TString::Format("&& eventClass==0"));
      if (c==1) signalR[c] = (RooDataSet*) sigWeighted.reduce(RooArgSet(*w->var("massReduced")),CutReducedMass+TString::Format("&& eventClass==1"));
	
      TString myCut;
      w->import(*signalR[c],RooFit::Rename(TString::Format("SigWeightReduced_cat%d",c)));
    }
    std::cout << std::endl;
    // Create full weighted signal data set without categorization
    RooDataSet* signalRAll = (RooDataSet*) sigWeighted.reduce(RooArgSet(*w->var("massReduced")),CutReducedMass);
    w->import(*signalRAll, RooFit::Rename("SigWeightReduced"));
    std::cout << "now signalRAll" << std::endl;
    signalRAll->Print("V");
    std::cout << "---- nX for signalRAll:  " << signalRAll->sumEntries() << std::endl; 
    std::cout << std::endl;
    std::cout << "==================================================================" << std::endl;
  
  }

  bool wantGenLevel=true;
  RooDataSet sigWeightedGen("sigWeightedGen","datasetGen",treesforfit[getBase(isample)],*rooVars, CutSignalRegion,"weight");   
  // -------------------------
  // split in categories, wrt genMass - to study the theory width
  if (wantGenLevel) { 
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "preparing dataset with observable mggGen, no split in categories since they're all the same" << std::endl;
    RooDataSet* signalG[NCAT];
    for (int c=0; c<ncat; ++c) {
      if (c==0) signalG[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),TString::Format("eventClass==0"));
      if (c==1) signalG[c] = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"),TString::Format("eventClass==1"));
      
      TString myCut;
      w->import(*signalG[c],RooFit::Rename(TString::Format("SigWeightGen_cat%d",c)));
    }
    RooDataSet* signalGAll = (RooDataSet*) sigWeightedGen.reduce(*w->var("mggGen"));
    w->import(*signalGAll, RooFit::Rename("SigWeightGen"));
    std::cout << "now signalGAll" << std::endl;
    signalGAll->Print("V");
    std::cout << "---- nX for signalGAll:  " << signalGAll->sumEntries() << std::endl; 
    std::cout << std::endl;
    std::cout << std::endl;
  }
  //w->writeToFile("data_001_5000.root");
  std::cout << "workspace summary" << std::endl;
  w->Print();


}
//-----------------------------------------------------------------------------------
theFitResult theFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, double minMassFit, double maxMassFit,std::string thefitrange){

  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  // params_test->printLatex();
  std::cout << "--------------------- BEFORE ITERATIONS-------------------------------" << std::endl;
  int stat=1;
  double minnll=10e8;
  RooFitResult *fitTest = new RooFitResult("","");
  theFitResult tmpfit; 

  double nllval = 0.;
  int minosstatus = 1;
  int hessestatus = 1;

  while (stat!=0){
    if (ntries>=MaxTries) break;
    // params_test->printLatex();
    //std::cout << "current try " << ntries << " stat=" << stat << " minnll=" << minnll << std::endl;
    std::cout << "--------------------- FITTING-------------------------------" << std::endl;
    // fitTest = pdf->fitTo(*data, RooFit::Minimizer("Minuit2","minimize"), RooFit::Offset(kTRUE), RooFit::Strategy(2), RooFit::PrintLevel(3), RooFit::Warnings(false), RooFit::SumW2Error(kTRUE), RooFit::Range(minMassFit,maxMassFit), RooFit::Save(kTRUE));
    // fitTest = pdf->fitTo(*data, RooFit::Minimizer("Minuit2","minimize"), RooFit::Offset(kTRUE), RooFit::Strategy(2), RooFit::PrintLevel(3), RooFit::Warnings(false), RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));
    // RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Offset(kTRUE),RooFit::Strategy(2));   

    //-------------------
    RooNLLVar *nll;
    if (thefitrange == "nofitrange"){  
      nll=new RooNLLVar("nll","nll",*pdf,*data);
    } else{
      nll=new RooNLLVar("nll","nll",*pdf,*data, RooFit::Range(thefitrange.c_str()));
    }
    RooMinimizer *minuit_fitTest = new RooMinimizer(*nll);
    // minuit_fitTest->setOffsetting(kTRUE);
    minuit_fitTest->setStrategy(1);
    minuit_fitTest->setPrintLevel(-1000);
    minuit_fitTest->minimize("Minuit2","minimize");
    std::cout << "Now running hesse" << std::endl;
    hessestatus = minuit_fitTest->hesse();
    std::cout << "Running minos" << std::endl;
    minosstatus = minuit_fitTest->minos();
    fitTest = minuit_fitTest->save("fitTest","fitTest");
    // offset= nll->offset();
    // std::cout << nll->isOffsetting() << std::endl;
    // std::cout << nll->offsetCarry() << std::endl;

    // When Offset(True) is used, the value returned by RooFitResult.minNll() is the
    // result of the minimization, which is not the value of the negative log
    // likelihood at the minimum. Storing the RooNLLVar and evaluating it at the
    // minimum actually evaluates it at the currently set parameter values, which will
    // give the minimum even if an offset is used in the minimization for numerical
    // stability.

    // I do not know why they used the two lines below in the past. 
    // They give minus the result we get with nll->getVal() and while 
    // with getVal we finally get some decent p-value, with the two lines
    // below we always get zero for p-value. 
    // minnll_woffset=fitTest->minNll();
    // minnll=-offset-minnll_woffset;
    nllval = nll->getVal();
    stat=fitTest->status();
    //-------------------

    if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
    ntries++; 
 
    tmpfit.fitres = fitTest;
    //this is with offset taking into account
    tmpfit.minNll = nllval;//minnll;
    
    tmpfit.minimizestatus = stat;
    tmpfit.hessestatus = hessestatus;
    tmpfit.minosstatus = minosstatus;
    
    
  }
  // std::cout << "------------------------OFFSET-----------------------------" << std::endl;
  // std::cout << "end of runFit stat=" << stat << " offset=" << offset << " minnll with offset=" << minnll_woffset << " diff= " << minnll<< " nll->getVal() " << nllval << std::endl;
  std::cout << "end of runFit stat=" << stat << " minnll with nll->getVal() " << nllval << std::endl;
  *stat_t = stat;
  *NLL = minnll;

  return tmpfit;
}

//-----------------------------------------------------------------------------------
void sigModelResponseFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year, const std::string &ws_dir, const std::string &samplename, FILE *resFileresp) {
  Int_t ncat = NCAT;
  RooDataSet* signal;
  RooDCBShape* responseadd[NCAT+1];
  int iMass = (int)  abs(mass);
  TCanvas* canv2 = new TCanvas( Form("Canvas2_M%f_k%s", mass , coupling.c_str()) , Form("Canvas2_M%f_k%s", mass , coupling.c_str()) );
  canv2->cd();
  // TPaveText* label_cms = get_labelcms(0, year , true);
  // TPaveText* label_sqrt = get_labelsqrt(0);

  double massMin = -600.; 
  double massMax = 600.;
  std::string thefitrange = "nofitrange";
  if (coupling == "1p4" || coupling == "5p6"){
    massMin = -100.;
    massMax = 100.;
    thefitrange = "fitrange";
  }
  // if(coupling=="001" || coupling == "0p014"){
  //   massMin=0.8*mass;
  //   massMax=1.2*mass;
  // }else if(coupling=="01" || coupling == "1p4"){
  //   massMin=0.6*mass;
  //   massMax=1.4*mass;
  // }else if(coupling=="02" || coupling == "5p6"){
  //   massMin=0.4*mass;
  //   massMax=1.6*mass;
  //   // massMin=0.8*mass;
  //   // massMax=1.2*mass;
  // }
  // if(massMin<300)massMin=300;
  // if(massMax>9000)massMax=9000;
  // if (iMass==750){massMin=600;massMax=850;}
  
  for(int c = 0; c<ncat+1; c++){
    std::cout << "RESPONSE FUNCTION FIT FOR CATEGORY " << c<< std::endl; 
    
    // if(c==1)continue;
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

    TString myCut;
    if(c==0||c==1)signal = (RooDataSet*)w->data(TString::Format("SigWeightReduced_cat%d",c));
    if(c==2)signal = (RooDataSet*)w->data("SigWeightReduced");

    //We will fit each distribution of the larger couplings up to 2*FWHM.
    // TH1F* hist = (TH1F*) signal->createHistogram( Form("reducedmass_sigHist_cat%d",c), *w->var("massReduced"), RooFit::Binning(120, massMin,massMax) );
    // // int lastnotoverflow = hist->GetXaxis()->GetNbins();
    // int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
    // int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
    // double fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
    // std::cout << "FWHM for cat " << c << " is "
    // 	      << fwhm << " = " << hist->GetBinCenter(bin2) << " - " << hist->GetBinCenter(bin1)
    // 	      << " with bin2 " << bin2 << " with bin1 " << bin1
    // 	      << " and maximum at bin " << hist->GetMaximum()
    // 	      << " with maximum value " << hist->GetBinCenter(hist->GetMaximum()) << std::endl;

    // fwhm = computeHistFHWM(hist);
    // std::cout << "FWHM for cat " << c << " is " << fwhm << std::endl;
    // //In case of large couplings for spin 0 samples don't trust events above 2*FWHM. 
    // if(coupling == "5p6" || coupling == "1p4"){
    //   massMin = hist->GetBinCenter(hist->GetMaximum()) - 2.*fwhm;
    //   massMax = hist->GetBinCenter(hist->GetMaximum()) + 2.*fwhm;
    // }

    //cb pos                                                                                                                     
    RooFormulaVar cbpos_mean(TString::Format("reducedmass_cbpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("reducedmass_sig_mean_cat%d",c)));
    RooFormulaVar cbpos_sigma(TString::Format("reducedmass_cbpos_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("reducedmass_sig_sigma_cat%d",c)));
    RooFormulaVar cbpos_alphacb(TString::Format("reducedmass_cbpos_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbpos_cat%d",c)));
    RooFormulaVar cbpos_n(TString::Format("reducedmass_cbpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_npos_cat%d",c)));     
    //cb neg
    RooFormulaVar cbneg_n(TString::Format("reducedmass_cbneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_nneg_cat%d",c)));
    RooFormulaVar cbneg_alphacb(TString::Format("reducedmass_cbneg_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbneg_cat%d",c)));   
   
    responseadd[c]= new  RooDCBShape(TString::Format("responseaddpdf_cat%d",c),TString::Format("responseaddpdf_cat%d",c) , *w->var("massReduced"), cbpos_mean, cbpos_sigma,  cbneg_alphacb, cbpos_alphacb,  cbneg_n,cbpos_n) ;

    w->import(*responseadd[c], RooFit::RecycleConflictNodes() );
    bool fixPar = true;
    if(mass==4000 && c==0 && samplename == "GluGluSpin0" && fixPar){
      w->var(TString::Format("reducedmass_sig_mean_cat%d",c))->setVal(-18.);
      w->var(TString::Format("reducedmass_sig_sigma_cat%d",c))->setVal(38.);
      w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c))->setVal(-1.);
      w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c))->setVal(1.);
      w->var(TString::Format("reducedmass_sig_npos_cat%d",c))->setVal(5);
      w->var(TString::Format("reducedmass_sig_nneg_cat%d",c))->setVal(5);
    }// else if(mass==1500 && fixPar){
    //   w->var(TString::Format("reducedmass_sig_mean_cat%d",c))->setVal(-12.);
    //   w->var(TString::Format("reducedmass_sig_sigma_cat%d",c))->setVal(25);
    //   w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c))->setVal(-2);
    //   w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c))->setVal(1);
    //   w->var(TString::Format("reducedmass_sig_npos_cat%d",c))->setVal(3);
    //   w->var(TString::Format("reducedmass_sig_nneg_cat%d",c))->setVal(7);
    // }
    
    int fitStatus = 0;
    double thisNll = 0.;
    if (coupling == "1p4" || coupling == "5p6"){
      w->var("massReduced")->setRange("fitrange",massMin, massMax);
    }
    theFitResult fitresults = theFit(responseadd[c], signal, &thisNll, &fitStatus, /*max iterations*/ 3, massMin, massMax, thefitrange) ;//"nofitrange" for normal range

    // RooFitResult* fitresults = (RooFitResult* ) responseadd[c]->fitTo(*signal,RooFit::SumW2Error(kTRUE), RooFit::Range(-600, 600), RooFit::Save(kTRUE), RooFit::Minos(kTRUE), RooFit::Strategy(2) );
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults.fitres->Print("V");
   
    double prob;
    RooRealVar norm("norm","norm",signal->sumEntries(),0,10E6);
    //norm.removeRange();
    RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*responseadd[c],norm);
    int nBinsForMassreduced = 120;
    int np = pdf->getParameters(*signal)->getSize();
 
    RooRealVar* massreduced = (RooRealVar*) w->var("massReduced");
    RooPlot* plotg = massreduced->frame(Range(-600., 600.),Title("mass reduced"),Bins(nBinsForMassreduced));
    signal->plotOn(plotg,Name("signal"));

    responseadd[c]->plotOn(plotg, LineColor(kBlue), Name("pdf"));
    // responseadd[c]->plotOn(plotg,Components(TString::Format("responsecbneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
    //responseadd[c]->plotOn(plotg,Components(TString::Format("responsecbpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    double chi2 = plotg->chiSquare("pdf","signal",np);
    std::cout << "[INFO] Calculating GOF for pdf " << responseadd[c]->GetName() << ", using " <<np << " fitted parameters" <<std::endl;
    prob = TMath::Prob(chi2*(nBinsForMassreduced-np),nBinsForMassreduced-np);
    //chi2/ndof is chi2
    std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMassreduced-np) << std::endl;
    std::cout << "[INFO] p-value  =  " << prob << std::endl;
    
    plotg->GetYaxis()->SetRangeUser(0.01,plotg->GetMaximum()*10 );
    plotg->GetXaxis()->SetTitle("#Delta m [GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    // plotg->GetXaxis()->SetTitlesize(0.07);
    // plotg->GetXaxis()->SetTitleoffset(1.40);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
    // legmc->SetTextSize(0.0286044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    // legmc->SetHeader( "#splitline{m_{x}=150 gev}{#splitline{}{class 0}}");
    legmc->AddEntry(plotg->getObject(0),"Simulation","lp");    
   
    legmc->AddEntry(plotg->getObject(1),"Double-Sided CB ","l");
    // legmc->AddEntry(plotg->getObject(2),"cb 1","l");   
    //legmc->AddEntry(plotg->getObject(3),"cb 2","l");   

    
   
    plotg->Draw();
    
    // lat->draw("same");
    // latex->draw("same");
    legmc->Draw("same");
    //    int ipos=11 ;
    //cms_lumi( c1,true,ipos );
    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    canv2->SetLogy(0);


    canv2->SaveAs(TString::Format("%s/%s/responsefcnfit/%s_responsefcnfitcb_cat%d_M%d_k_%s.png", ws_dir.c_str(),year.c_str(),samplename.c_str(),c,iMass,coupling.c_str())); 
    canv2->SetLogy();
    
    canv2->SaveAs(TString::Format("%s/%s/responsefcnfit/%s_responsefcnfitcb_cat%d_M%d_k_%s_log.png", ws_dir.c_str(),year.c_str(),samplename.c_str(),c,iMass,coupling.c_str())); 

    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    canv2->SetLogy(0);  
	
    w->defineSet(TString::Format("responseaddpdfparam_cat%d",c),RooArgSet(*w->var(TString::Format("reducedmass_sig_sigma_cat%d",c)), 
    									  *w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c)),
    									  *w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c)),
    									  *w->var(TString::Format("reducedmass_sig_npos_cat%d",c)),
    									  *w->var(TString::Format("reducedmass_sig_nneg_cat%d",c)),	   
    									  *w->var(TString::Format("reducedmass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("responseaddpdfparam_cat%d",c)));

    //Let's save the results to file
    float coup = widthtonum(coupling);
    fprintf(resFileresp,"%s & %10.0f & %10.1e & %d & %10.2f & %10.2f & (%d,%d,%d) \\\\\n", samplename.c_str() , mass, coup, c, chi2, prob, fitresults.minimizestatus ,fitresults.hessestatus ,fitresults.minosstatus );

  }//end of loop over categories
  
}

//-----------------------------------------------------------------------------------
void sigModelGenFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year, const std::string &ws_dir, const std::string &samplename, FILE *resFilegen) {
  int ncat = NCAT;
  RooDataSet* signal;
  //  RooBreitWigner* bw[NCAT+1];
  // RooGenericPdf* bw[NCAT+1];
  // RooVoigtian* mgenadd[NCAT+1];
  // RooDCBShape* mgenadd_dcb[NCAT+1];
  // RooGaussian* mgenadd_gaus[NCAT+1];
  // RooExponential* mgenadd_expo[NCAT+1];
  // RooAddPdf* mgenadd[NCAT+1];
  RooDCBShape* mgenadd[NCAT+1];
  // RooFFTConvPdf* mgenadd[NCAT+1];

  int iMass =  (int) abs(mass);
  TCanvas* canv3 = new TCanvas( Form("Canvas3_M%f_k%s", mass , coupling.c_str()) , Form("Canvas3_M%f_k%s", mass , coupling.c_str()) );
  canv3->cd();
  // TPaveText* label_cms = get_labelcms(0, "2016", true);
  // TPaveText* label_sqrt = get_labelsqrt(0);

  double massMin = 0.; 
  double massMax = 0.;
  // double curcoup = 0.;
  if(coupling=="001" || coupling == "0p014"){
    massMin=0.8*mass;
    massMax=1.2*mass;
    // curcoup = 1.4 * pow(10.,-4);
  }else if(coupling=="01" || coupling == "1p4"){
    massMin=0.8*mass;
    massMax=1.2*mass;
    // curcoup = 1.4 * pow(10.,-2);
  }else if(coupling=="02" || coupling == "5p6"){
    massMin=0.8*mass;
    massMax=1.2*mass;
    // curcoup = 5.6 * pow(10.,-2);
    // massMin=0.8*mass;
    // massMax=1.2*mass;
  }
  if(massMin<300)massMin=300;
  if(massMax>9000)massMax=9000;
  // if (iMass==750){massMin=600;massMax=850;}

  
  for(int c = 2; c<ncat+1; c++){
  // for(int c = 0; c<ncat; c++){

    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

    TString myCut;
    signal = (RooDataSet*)w->data("SigWeightGen");
    //We will fit each distribution of the larger couplings up to 2*FWHM.
    //Be careful to open a small window and not the whole range because for
    //spin 0 the distribution rises at the beginning. 
    TH1F* hist = (TH1F*) signal->createHistogram( Form("mgen_sigHist_cat%d",c), *w->var("mggGen"), RooFit::Binning(1000, massMin,massMax) );
    int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
    int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
    double fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
    std::cout << "FWHM for cat " << c << " is " << fwhm << " = " << hist->GetBinCenter(bin2) << " - " << hist->GetBinCenter(bin1) <<std::endl;
    // bool weirdshapelow = false;
    // bool weirdshapehigh = false;
    
    // if (mass - hist->GetBinCenter(bin1) > 500.){weirdshapelow=true;}
    // while (weirdshapelow){
    //   std::cout << "YYYYYYYYYYYYYYYYYYYYYY" << std::endl; 
    //   weirdshapelow = (mass - hist->GetBinCenter(bin1) > 500.);
    //   ++bin1;
    // }
    // if (hist->GetBinCenter(bin2) - mass > 500.){weirdshapehigh=true;}
    // while (weirdshapehigh){
    //   weirdshapehigh = (hist->GetBinCenter(bin2) - mass > 500.);
    //   --bin2;
    // }
    // if (weirdshapelow || weirdshapehigh) {fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);}
    // if (weirdshapelow){
    //   fwhm = 2*(mass - hist->GetBinCenter(bin2));
    // } else if (weirdshapehigh){
    //   fwhm = 2*hist->GetBinCenter(bin1);
    // }

    w->var(TString::Format("mgen_sig_mean_cat%d",c))->setVal(mass);
    // w->var(TString::Format("mgen_sig_mean_cat%d",c))->setRange(mass*0.8, mass*1.2);
    std::cout << coupling << " "  << year << " " << mass << std::endl;
    //------------------------------------------------------------
    //Some tweaking first for RS samples
    if (coupling == "02" && year == "2017"){
      if (mass == 2000. || mass == 2250.){
	w->var(TString::Format("mgen_sig_sigma_cat%d",c))->setVal(15.);
      } else if (mass == 5250. || mass == 6500.){
	w->var(TString::Format("mgen_sig_sigma_cat%d",c))->setVal(30.);
      } else if (mass == 8000.){
	w->var(TString::Format("mgen_sig_sigma_cat%d",c))->setVal(150.);
      }
    }//end of if over RS 2017 with 02 coup
    //------------------------------------------------------------
    //Now to GluGlu
    if (coupling == "1p4" ){
      if (mass >= 1000.){
	w->var(TString::Format("mgen_sig_sigma_cat%d",c))->setVal(50.);
      }
      if(mass == 4000. && year == "2017"){
	w->var(TString::Format("mgen_sig_nneg_cat%d",c))->setVal(1e-05);
      }
      if( (mass == 4750. || mass == 4250.) && year == "2018" ) {
	w->var(TString::Format("mgen_sig_nneg_cat%d",c))->setVal(1e-03);
	w->var(TString::Format("mgen_sig_alphacbneg_cat%d",c))->setVal(0.2);
      }
    }//end of if over GluGlu 1p4
    if (coupling == "5p6"){
      // if (mass >= 2000. && mass <= 2500. ){
      if (mass == 2000.){
	w->var(TString::Format("mgen_sig_nneg_cat%d",c))->setVal(1e-03);
      } else if (mass >= 1000. && mass < 3500. ){
	// w->var(TString::Format("mgen_sig_nneg_cat%d",c))->setVal(1e-03);
	w->var(TString::Format("mgen_sig_sigma_cat%d",c))->setVal(50.);
	// w->var(TString::Format("mgen_sig_alphacbneg_cat%d",c))->setVal(10.);
      } else if (mass == 3500. ){
	w->var(TString::Format("mgen_sig_sigma_cat%d",c))->setVal(100.);
	// w->var(TString::Format("mgen_sig_alphacbneg_cat%d",c))->setVal(10.);
      }
    }//end of if over GluGlu 5p6 
    if (coupling == "0p014"){
      if (mass >= 2000.){
	// w->var(TString::Format("mgen_sig_nneg_cat%d",c))->setVal(1e-03);
	w->var(TString::Format("mgen_sig_sigma_cat%d",c))->setVal(50.);
      }
    }//end of if over GluGlu 0p014

    //In case of large couplings for spin 0 samples don't trust events above 2*FWHM. 
    if(coupling == "5p6" || coupling == "1p4"){
      massMin = mass - 2.*fwhm;
      massMax = mass + 2.*fwhm;
    }
    
    //Below are the DCB parameters, but will use the same with voigtian where possible
    //cb pos                                                                                                                   
    RooFormulaVar cbpos_mean(TString::Format("mgen_cbpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("mgen_sig_mean_cat%d",c)));
    RooFormulaVar cbpos_sigma(TString::Format("mgen_cbpos_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("mgen_sig_sigma_cat%d",c)));
    RooFormulaVar cbpos_alphacb(TString::Format("mgen_cbpos_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbpos_cat%d",c)));
    RooFormulaVar cbpos_n(TString::Format("mgen_cbpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_npos_cat%d",c)));    
    //cb neg
    RooFormulaVar cbneg_n(TString::Format("mgen_cbneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_nneg_cat%d",c)));
    RooFormulaVar cbneg_alphacb(TString::Format("mgen_cbneg_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbneg_cat%d",c)));    
    mgenadd[c]=  new RooDCBShape(TString::Format("mgenaddpdf_cat%d",c),TString::Format("mgenaddpdf_cat%d",c) , *w->var("mggGen"), cbpos_mean, cbpos_sigma,  cbneg_alphacb, cbpos_alphacb,  cbneg_n,cbpos_n) ;
    // mgenadd_dcb[c] =  new  RooDCBShape(TString::Format("mgenadd_dcbpdf_cat%d",c),TString::Format("mgenadd_dcbpdf_cat%d",c), *w->var("mggGen"), cbpos_mean, cbpos_sigma,  cbneg_alphacb, cbpos_alphacb,  cbneg_n,cbpos_n) ;
    // RooRealVar mg("mg","mg", mass - 500., 0. ,  mass + 1000.) ;
    // RooRealVar sg("sg","sg",100,0.,200.) ;
    // RooFormulaVar expo_lambda(TString::Format("mgen_expolambda_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_expolambda_cat%d",c)));
    // RooRealVar lambda("lambda", "slope", -2., -10., 0.);
    // mgenadd_gaus[c] = new RooGaussian(TString::Format("mgenadd_gauspdf_cat%d",c),TString::Format("mgenadd_gauspdf_cat%d",c),*w->var("mggGen"),mg,sg);
    // mgenadd_expo[c] = new RooExponential(TString::Format("mgenadd_expopdf_cat%d",c),TString::Format("mgenadd_expopdf_cat%d",c),*w->var("mggGen"), expo_lambda);
    // RooRealVar* var = new RooRealVar("var", "var", massMin, massMax);
    // w->var("mggGen")->setBins(4000, "cache");
    //Convolution
    // mgenadd[c] = new RooFFTConvPdf(TString::Format("mgenaddpdf_cat%d",c),TString::Format("mgenaddpdf_cat%d",c) , *w->var("mggGen"), *mgenadd_dcb[c], *mgenadd_gaus[c]);

    //add the DCB with an exponential
    // RooArgList *pdfs_holder = new RooArgList();
    // pdfs_holder->add(*mgenadd_dcb[c]);
    // pdfs_holder->add(*mgenadd_expo[c]);
    // RooArgList *coeffs_holder= new RooArgList();
    // coeffs_holder->add( *w->var(TString::Format("mgen_sig_frac_cat%d",c)) );
    // mgenadd[c] = new RooAddPdf(TString::Format("mgenaddpdf_cat%d",c),TString::Format("mgenaddpdf_cat%d",c),*pdfs_holder,*coeffs_holder, true);
    mgenadd[c]->Print("V");
    
    // // add the DCB and the Gaussian
    // RooArgList *pdfs_holder = new RooArgList();
    // pdfs_holder->add(*mgenadd_dcb[c]);
    // pdfs_holder->add(*mgenadd_gaus[c]);
    // RooArgList *coeffs_holder= new RooArgList();
    // coeffs_holder->add( *w->var(TString::Format("mgen_sig_frac_cat%d",c)) );
    // mgenadd[c] = new RooAddPdf(TString::Format("mgenaddpdf_cat%d",c),TString::Format("mgenaddpdf_cat%d",c),*pdfs_holder,*coeffs_holder, true);
    // w->var(TString::Format("mgen_width_sig_cat%d",c))->setVal(curcoup);
    // RooFormulaVar width(TString::Format("mgen_voiwidth_sig_cat%d",c),"","@0", *w->var( TString::Format("mgen_width_sig_cat%d",c)) );
    
    // mgenadd[c]=  new  RooVoigtian(TString::Format("mgenaddpdf_cat%d",c),TString::Format("mgenaddpdf_cat%d",c) , *w->var("mggGen"), cbpos_mean, cbpos_sigma, width );
    
    
    w->import(*mgenadd[c]);
    std::cout<<mass<<" "<<signal->sumEntries()<<std::endl;
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;

    int fitStatus = 0;
    double thisNll = 0.;
    w->var("mggGen")->setRange("fitrange",massMin, massMax) ;
    // RooFitResult* fitresult = (RooFitResult* ) mgenadd[c]->fitTo(*signal,SumW2Error(kTRUE), RooFit::Range("fitrange"), RooFit::Minos(kTRUE), RooFit::Save(kTRUE));  
    theFitResult fitresults = theFit(mgenadd[c], signal, &thisNll, &fitStatus, /*max iterations*/ 3, massMin, massMax, "fitrange") ;

    std::cout << "massMin " << massMin << " massMax "<< massMax << std::endl;
    fitresults.fitres->Print("V");
    // fitresult->Print("V");
    
    // theFitResult fitresults;
    // fitresults.fitres = fitresult;

    int nBinsForMassreduced = 120;
    w->var("mggGen")->setRange("plotrange",mass*0.8,mass*1.2) ;
    RooPlot* plotg = w->var("mggGen")->frame(RooFit::Range("plotrange"),Title("mass generated"),Bins(nBinsForMassreduced));
    signal->plotOn(plotg,Name("signal"));
    mgenadd[c]->plotOn(plotg, LineColor(kBlue), Name("pdf"), RooFit::Range("plotrange"));

    double prob;
    RooRealVar norm("norm","norm",signal->sumEntries(),0,10E6);
    //norm.removeRange();
    RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mgenadd[c],norm);
    int np = pdf->getParameters(*signal)->getSize();
 
    double chi2 = plotg->chiSquare("pdf","signal",np);
    std::cout << "[INFO] Calculating GOF for pdf " << mgenadd[c]->GetName() << ", using " <<np << " fitted parameters" <<std::endl;
    prob = TMath::Prob(chi2*(nBinsForMassreduced-np),nBinsForMassreduced-np);
    //chi2/ndof is chi2
    std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMassreduced-np) << std::endl;
    std::cout << "[INFO] p-value  =  " << prob << std::endl;

   
    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma}^{gen} [GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotg->getObject(0),"Gen Simulation","lp");    
    legmc->AddEntry(plotg->getObject(1),"Double-Sided CB ","l");
    //  legmc->AddEntry(plotg->getObject(2),"cb 1","l");   
    // legmc->AddEntry(plotg->getObject(3),"cb 2","l");   
 
    // bw[c]->plotOn(plotg, LineColor(kBlue));
    plotg->GetYaxis()->SetRangeUser(0.01,plotg->GetMaximum()*1.2 );
    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    // plotg->GetXaxis()->SetTitlesize(0.07);
    // plotg->GetXaxis()->SetTitleoffset(1.40);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    plotg->Draw();
    legmc->Draw("same");
    canv3->SetLogy(0);


    canv3->SaveAs(TString::Format("%s/%s/mgen/%s_mgen_cat%d_M%d_k_%s.png", ws_dir.c_str(),year.c_str(),samplename.c_str(),c,iMass,coupling.c_str())); 
    canv3->SetLogy();
    
    canv3->SaveAs(TString::Format("%s/%s/mgen/%s_mgen_cat%d_M%d_k_%s_log.png", ws_dir.c_str(),year.c_str(),samplename.c_str(),c,iMass,coupling.c_str())); 
    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    canv3->SetLogy(0);  
	
    w->defineSet(TString::Format("mgenpdfparam_cat%d",c),RooArgSet(*w->var(TString::Format("mgen_sig_sigma_cat%d",c)), 
    								   *w->var(TString::Format("mgen_sig_alphacbpos_cat%d",c)),
    								   *w->var(TString::Format("mgen_sig_alphacbneg_cat%d",c)),
    								   *w->var(TString::Format("mgen_sig_npos_cat%d",c)),
    								   *w->var(TString::Format("mgen_sig_nneg_cat%d",c)),	   
    								   *w->var(TString::Format("mgen_sig_frac_cat%d",c)),  
    								   *w->var(TString::Format("mgen_sig_mean_cat%d",c))));
    // w->defineSet(TString::Format("mgenpdfparam_cat%d",c),RooArgSet(*w->var(TString::Format("mgen_sig_sigma_cat%d",c)), 
    // 								   *w->var(TString::Format("mgen_width_sig_cat%d",c)),  
    // 								   *w->var(TString::Format("mgen_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("mgenpdfparam_cat%d",c)));

    //Let's save the results to file
    float coup = widthtonum(coupling);
    fprintf(resFilegen,"%s & %10.0f & %10.1e & %d & %10.2f & %10.2f & (%d,%d,%d) \\\\\n", samplename.c_str() , mass, coup, c, chi2, prob, fitresults.minimizestatus ,fitresults.hessestatus ,fitresults.minosstatus );
    // fprintf(resFilegen,"%s & %10.0f & %10.1e & %d & %10.2f & %10.2f & %d \\\\\n", samplename.c_str() , mass, coup, c, chi2, prob, fitresult->status() );


    
  }

}
//-----------------------------------------------------------------------------------
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample, FILE *resFilegen, FILE *resFileresp){

  std::string coupling= ""; 
  if (isample.find("RSGraviton") != std::string::npos) {
    coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
  } else if ( isample.find("GluGluSpin0") != std::string::npos ){
    coupling = get_str_between_two_str(getBase(isample), "GluGluSpin0ToGammaGamma_W_", "_M_");
  }
  std::string M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");
  std::string samplename = isample.substr(0,getBase(isample).find("ToGammaGamma")); 
  std::cout << ws_dir << std::endl;

  TFile *fout = new TFile(Form("%s/%s/workspaces/ws_%s_ResponseAndGen_M%s_k%s_%s.root", ws_dir.c_str(), year.c_str(), samplename.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()), "recreate");

  HLFactory hlf("HLFactory", "HighMassGG.rs", false);
  RooWorkspace* w = hlf.GetWs();
  // range for the variables
  w->var("mgg")->setMin(MINmass);
  w->var("mgg")->setMax(MAXmass);
  w->var("mggGen")->setMin(MINmass);
  w->var("mggGen")->setMax(MAXmass);
  // w->Print("V");

  //========================================================================
  std::cout << "Adding signal data for mass " << M_bins << " and coupling " << coupling <<std::endl;
  AddSigData(w, stof(M_bins), coupling, isample);
  //========================================================================

  TCanvas* canv1 = new TCanvas(Form("Canvas1_%s",isample.c_str()),Form("Canvas_%s",isample.c_str()));
  canv1->cd();
  // RooPlot* p = w->var("mgg")->frame(RooFit::Range(230., 10000.), RooFit::Bins( 2385));

  //========================================================================
  std::cout << "Fit the response function" <<std::endl;
  sigModelResponseFcnFit(w, stof(M_bins), coupling, year, ws_dir, samplename, resFileresp);
  //========================================================================
  std::cout << "Fit the gen function" <<std::endl;
  sigModelGenFcnFit(w, stof(M_bins), coupling, year, ws_dir, samplename, resFilegen);
  //========================================================================
  //Save the resulting workspace
  w->Print("V");
  fout->cd();
  w->Write();
  fout->Write();
  fout->Close();


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

float widthtonum(std::string coupling){
  float coup = 0.; 

  if ( coupling == "001" || coupling == "0p014"){coup = 1.4 * pow(10.,-4);}
  else if ( coupling == "01" || coupling == "1p4" ){coup = 1.4 * pow(10.,-2);}
  else if ( coupling == "02" || coupling == "5p6" ){coup = 5.6 * pow(10.,-2);}

  return coup;


}


double computeHistFHWM(TH1F*  hist){

  double  halfMaxVal = 0.5*hist->GetMaximum();
  int  maxBin = hist->GetMaximumBin();
        
  int  binLeft = 0; 
  int  binRight = 0;
  double xLeft = 0.;
  double xRight = 0.;
  double xWidth;
      
  int  last = 1;
  for(int ibin = 1; ibin<=maxBin;ibin++){
    double   binVal = hist->GetBinContent(ibin);
    if (binVal >= halfMaxVal){
      binLeft = last;
      break;
    }
    if (binVal>0) last=ibin;
  }
  last = hist->GetXaxis()->GetNbins()+1;
  for(int ibin = hist->GetXaxis()->GetNbins()+1; ibin>maxBin;ibin--){ 
    double   binVal = hist->GetBinContent(ibin);
    if (binVal >= halfMaxVal){
      binRight = last;
      break;
    }
    if (binVal>0) last=ibin;
  }
  for(int ibin = 1; ibin<=maxBin;ibin++){ 
    double  binVal = hist->GetBinContent(ibin);
    if (binVal >= halfMaxVal){
      binLeft = ibin;
      break;
    }
  }
  for(int ibin = maxBin+1; ibin<=hist->GetXaxis()->GetNbins()+1;ibin++){  
    double binVal = hist->GetBinContent(ibin);
    if (binVal < halfMaxVal){
      binRight = ibin-1;
      break;
    }
  } 
  xWidth = 0.;
  if (binLeft > 0 && binRight > 0 ){
    xLeft = hist->GetXaxis()->GetBinCenter(binLeft);
    xRight = hist->GetXaxis()->GetBinCenter(binRight);
    xWidth = xRight-xLeft;
    std::cout<<TString::Format("FWHM = %f",xWidth)<<std::endl;
  }  else{
    std::cout<<"Did not succeed to compute the FWHM"<<std::endl;
  }

  return xWidth;
}
