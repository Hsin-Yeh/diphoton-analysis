#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"

//RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooStats/HLFactory.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooBinning.h"
#include "RooVoigtian.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "RooExtendPdf.h"

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
void runAllFits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &checksample);
void runfits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &isample, FILE *resFileasimovfit, FILE *resFilerecofit);
theFitResult theFit(RooAbsPdf *pdf, RooDataHist *data, double *NLL, int *stat_t, int MaxTries, double minMassFit, double maxMassFit);
void sigModelShapeFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year, const std::string &ws_outdir, const std::string &samplename, bool doPlot, FILE *resFilerecofit);
void asimovDatasetFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &ws_outdir, const std::string &samplename, FILE *resFileasimovfit);
RooDataHist* throwAsimov( double nexp, RooAbsPdf *pdf, RooRealVar *x, RooDataHist *asimov );
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);
void SetConstantParams(const RooArgSet* params) ;
TPaveText* get_labelsqrt( int legendquadrant );
TPaveText* get_labelcms( int legendquadrant, std::string year, bool sim);
float widthtonum(std::string coupling);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputsamples, inputworkspaces, outputworkspaces, checksample;

  if(argc!=6) {
    std::cout << "Syntax: prepareShapes.exe [2016/2017/2018] [inputsamples] [inputworkspaces] [outputworkspaces]" << std::endl;
      return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    inputsamples = argv[2];
    inputworkspaces = argv[3];
    outputworkspaces = argv[4];
    checksample = argv[5];

 }

  //========================================================================
  //read the reduced input root trees
  initForFit(inputsamples);
  //========================================================================
  //Run the fits
  runAllFits(year, inputworkspaces, outputworkspaces,checksample);

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void sigModelShapeFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year, const std::string &ws_outdir, const std::string &samplename, bool doPlot, FILE *resFilerecofit) {
  int ncat = NCAT;
  RooDataSet* signalK;
  RooDCBShape* reducedmass[NCAT+1];
  RooDCBShape* mgen[NCAT+1];
  RooFFTConvPdf* conv[NCAT+1];
  int iMass =  (int) abs(mass);
  TCanvas* canv1 = new TCanvas( Form("Canvas1_M%f_k%s", mass , coupling.c_str()) , Form("Canvas1_M%f_k%s", mass , coupling.c_str()) );
  canv1->cd();
  // TPaveText* label_cms = get_labelcms(0, "2016", true);
  // TPaveText* label_sqrt = get_labelsqrt(0);
  
  double massMin = 0.; 
  double massMax = 0.;
  if(coupling=="001" || coupling == "0p014"){
    massMin=0.8*mass;
    massMax=1.2*mass;
  }else if(coupling=="01" || coupling == "1p4"){
    massMin=0.6*mass;
    massMax=1.4*mass;
  }else if(coupling=="02" || coupling == "5p6"){
    massMin=0.4*mass;
    massMax=1.6*mass;
  }

  for(int c = 0; c<ncat+1; c++){
    std::cout << "FINAL SHAPE FIT FOR CATEGORY " << c<< std::endl; 
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();
    TString myCut;
      
    //pickup the functions from the ws
    RooRealVar* var = new RooRealVar("var", "var", massMin, massMax);
    RooRealVar* MH = new RooRealVar("MH", "MH",mass);
    w->import(*MH);
    MH->setConstant();
    //response
    // double mh =100.;
    //gen shape
    //cb pos                                                                                                                     
    RooFormulaVar cbpos_mean_mgen("cbpos_sig_mean_cat2_mgen","", "@0", *w->var("mgen_sig_mean_cat2"));
    RooFormulaVar cbpos_sigma_mgen("cbpos_sig_sigma_cat2_mgen", "", "sqrt(@0*@0)", *w->var("mgen_sig_sigma_cat2"));
    RooFormulaVar cbpos_alphacb_mgen("cbpos_sig_alphacb_cat2_mgen","", "@0", *w->var( "mgen_sig_alphacbpos_cat2"));
    RooFormulaVar cbpos_n_mgen("cbpos_sig_n_cat2_mgen","", "@0", *w->var( "mgen_sig_npos_cat2"));
    //cb neg
    RooFormulaVar cbneg_n_mgen("cbneg_sig_n_cat2_mgen","", "@0", *w->var( "mgen_sig_nneg_cat2"));
    RooFormulaVar cbneg_alphacb_mgen("cbneg_sig_alphacb_cat2_mgen","", "@0", *w->var( "mgen_sig_alphacbneg_cat2"));

    // RooFormulaVar cbpos_mean_mgen(Form("cbpos_sig_mean_cat%d_mgen",c),"", "@0", *w->var(Form("mgen_sig_mean_cat%d",c)));
    // RooFormulaVar cbpos_sigma_mgen(Form("cbpos_sig_sigma_cat%d_mgen",c), "", "sqrt(@0*@0)", *w->var(Form("mgen_sig_sigma_cat%d",c)));
    // RooFormulaVar cbpos_alphacb_mgen(Form("cbpos_sig_alphacb_cat%d_mgen",c),"", "@0", *w->var( Form("mgen_sig_alphacbpos_cat%d",c)));
    // RooFormulaVar cbpos_n_mgen(Form("cbpos_sig_n_cat%d_mgen",c),"", "@0", *w->var( Form("mgen_sig_npos_cat%d",c)));
    // //cb neg
    // RooFormulaVar cbneg_n_mgen(Form("cbneg_sig_n_cat%d_mgen",c),"", "@0", *w->var( Form("mgen_sig_nneg_cat%d",c)));
    // RooFormulaVar cbneg_alphacb_mgen(Form("cbneg_sig_alphacb_cat%d_mgen",c),"", "@0", *w->var( Form("mgen_sig_alphacbneg_cat%d",c)));


    mgen[c]= new  RooDCBShape(TString::Format("mgen_cat%d",c),TString::Format("mgen_cat%d",c) , *var, cbpos_mean_mgen, cbpos_sigma_mgen,  cbneg_alphacb_mgen, cbpos_alphacb_mgen,  cbneg_n_mgen,cbpos_n_mgen) ;
    //response shape
    //cb pos                                  
    // double scale=1;
    RooFormulaVar cbpos_mean_reducedmass(TString::Format("cbpos_sig_mean_cat%d_reducedmass",c),"","@0",RooArgList(*w->var(TString::Format("reducedmass_sig_mean_cat%d",c)),*w->var("MH"))  );
    RooFormulaVar cbpos_sigma_reducedmass(TString::Format("cbpos_sig_sigma_cat%d_reducedmass",c), "", "sqrt(@0*@0)", RooArgList(*w->var(TString::Format("reducedmass_sig_sigma_cat%d",c)),*w->var(TString::Format("mgen_sig_mean_cat%d",c))));
    RooFormulaVar cbpos_alphacb_reducedmass(TString::Format("cbpos_sig_alphacb_cat%d_reducedmass",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbpos_cat%d",c)));
    RooFormulaVar cbpos_n_reducedmass(TString::Format("cbpos_sig_n_cat%d_reducedmass",c),"", "@0", *w->var( TString::Format("reducedmass_sig_npos_cat%d",c)));
    //cb neg
    RooFormulaVar cbneg_n_reducedmass(TString::Format("cbneg_sig_n_cat%d_reducedmass",c),"", "@0", *w->var( TString::Format("reducedmass_sig_nneg_cat%d",c)));
    RooFormulaVar cbneg_alphacb_reducedmass(TString::Format("cbneg_sig_alphacb_cat%d_reducedmass",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbneg_cat%d",c)));
    reducedmass[c]= new  RooDCBShape(TString::Format("addpdf_cat%d_reducedmass",c), TString::Format("addpdf_cat%d_reducedmass",c), *var, cbpos_mean_reducedmass, cbpos_sigma_reducedmass,  cbneg_alphacb_reducedmass, cbpos_alphacb_reducedmass,  cbneg_n_reducedmass,cbpos_n_reducedmass) ;
    
    //convolution
    var->setBins(4000, "cache");
    conv[c] = new RooFFTConvPdf(TString::Format("convpdf_cat%d",c),TString::Format("convpdf_cat%d",c), *var,*mgen[c],*reducedmass[c]);
    conv[c]->Print("V");
    std::cout<<"******************************************************************************************"<<std::endl;
    RooArgSet* params = conv[c]->getParameters(*w->var("mgg")) ;
    params->Print("v") ;
   
    //asimov dataset
    // double min;
    // double max;
    // min=mass*0.8;
    // max=mass*1.2;
    if(massMin<300)massMin=300;
    if(massMax>9000)massMax=9000;
    
    RooBinning binning(160, massMin, massMax, "binning");
    TH1F* dummy = new TH1F("dummy", "dummy", 160, massMin, massMax);
    RooDataHist* datahist_asimov = new RooDataHist("datahist_asimov","datahist_asimov",*var, dummy, 1);
    // RooDataHist* datahist_asimov_clone = new RooDataHist("datahist_asimov","datahist_asimov",*w->var("mgg"), dummy, 1);
      
    //  RooDataHist* datahist_asimov; //= mgen[c]->generateBinned(RooArgSet(*var),RooFit::NumEvents(signal->numEntries()), RooFit::Asimov());
    //datahist_asimov = signal->binnedClone();
    throwAsimov(20000, conv[c], var,datahist_asimov); 

    if(c<2)    datahist_asimov->SetName(TString::Format("signal_asimov_cat%d",c));
    if(c==2)    datahist_asimov->SetName("signal_asimov");
    std::cout<<"--------------------------------------------------- "<<datahist_asimov->sumEntries()<<std::endl;
    w->import(*datahist_asimov);

    if(doPlot){
      //plotting...
      if(c==0||c==1)signalK = (RooDataSet*)w->data(TString::Format("SigWeight_cat%d",c));
      if(c==2)signalK = (RooDataSet*)w->data("SigWeight");
 
      w->Print("v");
      
      int nBinsForMassreco = 160;
      RooPlot* plotgg = (w->var("mgg"))->frame(Range(massMin,massMax),Title("mass generated"),Bins(nBinsForMassreco));
      signalK->Print("v");
      signalK->plotOn(plotgg,Name("signalreco"));
      std::cout<<"flag3"<<std::endl;
      RooPlot* plotg = var->frame(Range(massMin,massMax),Title("mass generated"),Bins(nBinsForMassreco));
      RooDataSet* fakedata= conv[c]->generate(*var, signalK ->sumEntries());
      fakedata->plotOn(plotg, RooFit::Invisible());
      // reducedmass[c]->plotOn(plotg, LineColor(kRed));
      mgen[c]->plotOn(plotg, LineColor(kGreen));
      conv[c]->plotOn(plotg, LineColor(kBlue), Name("convpdf"));

      double prob;
      RooRealVar norm("norm","norm",signalK->sumEntries(),0,10E6);
      //norm.removeRange();
      RooExtendPdf *convpdf = new RooExtendPdf("ext","ext",*conv[c],norm);
      int np = convpdf->getParameters(*signalK)->selectByAttrib("Constant",kFALSE)->getSize();
 
      double chi2 = plotg->chiSquare("convpdf","signalreco",np);
       std::cout << "[INFO] Calculating GOF for pdf " << conv[c]->GetName() << ", using " <<np << " fitted parameters" <<std::endl;
      prob = TMath::Prob(chi2*(nBinsForMassreco-np),nBinsForMassreco-np);
      //chi2/ndof is chi2
      std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMassreco-np) << std::endl;
      std::cout << "[INFO] p-value  =  " << prob << std::endl;
 
      // datahist_asimov->plotOn(plotg, MarkerColor(kRed), MarkerSize(0.5));
      //check parameters 
      RooArgSet* model_params = conv[c]->getParameters(*w->var("mgg")) ;
      model_params->Print("v") ;
      std::cout<<"flag4"<<std::endl;
    
      plotg->GetXaxis()->SetTitle("m_{#gamma#gamma}[GeV]");
      plotg->GetXaxis()->SetTitleFont(42);
      plotg->GetYaxis()->SetTitle("arbitrary scale");
      plotg->GetYaxis()->SetTitleFont(42);
      plotg->GetYaxis()->SetTitleSize(0.04);
      TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->AddEntry(plotgg->getObject(0),"Reco Simulation","lp");    
      //    legmc->AddEntry(plotg->getObject(0),"Response ","l");
      legmc->AddEntry(plotg->getObject(1),"Gen shape","l");   
      legmc->AddEntry(plotg->getObject(2),"Convolution","l");   
    
      // bw[c]->plotOn(plotg, LineColor(kBlue));
      plotg->GetYaxis()->SetRangeUser(0.01,30000 );
      plotg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
      plotg->GetXaxis()->SetTitleFont(42);
      // plotg->GetXaxis()->SetTitlesize(0.07);
      // plotg->GetXaxis()->SetTitleoffset(1.40);
      plotg->GetYaxis()->SetTitle("arbitrary scale");
      plotg->GetYaxis()->SetTitleFont(42);
      plotg->GetYaxis()->SetTitleSize(0.04);
      plotg->Draw();
      plotgg->Draw("same");
      legmc->Draw("same");
      canv1->SetLogy(0);
    

      canv1->SaveAs(TString::Format("%s/../FinalConvShapes/%s_mreco_cat%d_M%d_k_%s.png", ws_outdir.c_str(), samplename.c_str(), c, iMass, coupling.c_str())); 

      canv1->SetLogy();			        
    					        
      canv1->SaveAs(TString::Format("%s/../FinalConvShapes/%s_mreco_cat%d_M%d_k_%s_log.png", ws_outdir.c_str(), samplename.c_str(), c, iMass,coupling.c_str())); 

      canv1->SetLogy(0);

      //Let's save the results to file
      float coup = widthtonum(coupling);
      fprintf(resFilerecofit,"%s & %10.0f & %10.1e & %d & %10.2f & %10.2f  \\\\\n", samplename.c_str() , mass, coup, c, chi2, prob);


      
    }//end of doPlot
    
  }//end of loop over categories
 
}

//-----------------------------------------------------------------------------------
void runAllFits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &checksample){

  std::vector<std::string> samples = getSampleListForFit();

  FILE *resFileasimovfit;

  resFileasimovfit = fopen(Form("%s/../FinalConvShapes/asimov_fits_%s.txt", ws_outdir.c_str(), year.c_str()),"w");
  fprintf(resFileasimovfit,"\\hline\n");
  fprintf(resFileasimovfit,"\\hline\n");
  fprintf(resFileasimovfit,"Sample & Mass & Width & Category & gof$(\\chi^{2}/ndof)$ & gof(prob) & (minim,hesse,minos) \\\\\n");
  fprintf(resFileasimovfit,"\\hline\n");

  FILE *resFilerecofit;

  resFilerecofit = fopen(Form("%s/../FinalConvShapes/reco_fits_%s.txt", ws_outdir.c_str(), year.c_str()),"w");
  fprintf(resFilerecofit,"\\hline\n");
  fprintf(resFilerecofit,"\\hline\n");
  fprintf(resFilerecofit,"Sample & Mass & Width & Category & gof$(\\chi^{2}/ndof)$ & gof(prob) \\\\\n");
  fprintf(resFilerecofit,"\\hline\n");

  
  for(auto isample : samples) {
    //We will process one year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    if ( isample.find(checksample) == std::string::npos && checksample!="all") continue;
    // if (isample.find("RSGraviton") == std::string::npos) { continue; }

    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    runfits(year, ws_indir, ws_outdir, isample, resFileasimovfit, resFilerecofit);
  }

  fclose(resFileasimovfit);
  fclose(resFilerecofit);

}

//-----------------------------------------------------------------------------------
void runfits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &isample, FILE *resFileasimovfit, FILE *resFilerecofit){

  std::string coupling= ""; 
  if (isample.find("RSGraviton") != std::string::npos) {
    coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
  } else if ( isample.find("GluGluSpin0") != std::string::npos ){
    coupling = get_str_between_two_str(getBase(isample), "GluGluSpin0ToGammaGamma_W_", "_M_");
  }
  std::string M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");
  std::string samplename = isample.substr(0,getBase(isample).find("ToGammaGamma")); 
  
  std::cout << "Input dir " << ws_indir << std::endl;
  std::cout << "Output dir " << ws_outdir << std::endl;

  TFile *fresandgen = TFile::Open(Form("%s/ws_%s_ResponseAndGen_M%s_k%s_%s.root", ws_indir.c_str(), samplename.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()));

  RooWorkspace* ws = (RooWorkspace*) fresandgen->Get("HLFactory_ws"); 
  //========================================================================
  std::cout << "Make convolution of the response and the gen level shape" <<std::endl;
  bool doPlot = true; //For the final plot
  sigModelShapeFcnFit(ws, stof(M_bins), coupling, year, ws_outdir, samplename, doPlot, resFilerecofit);
  //========================================================================
  std::cout << "Throw an asimov dataset from the DCB and fit it" <<std::endl;
  asimovDatasetFcnFit(ws, stof(M_bins), coupling, ws_outdir, samplename, resFileasimovfit);

  //========================================================================
  //Save the resulting workspace
  TFile *fout = new TFile(Form("%s/ws_%s_ResponseAndGen_M%s_k%s_%s_final.root", ws_outdir.c_str(), samplename.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()), "recreate");
  fout->cd();
  // ws->Print("V");
  ws->Write();
  fout->Write();
  fout->Close();

}

//-----------------------------------------------------------------------------------
theFitResult theFit(RooAbsPdf *pdf, RooDataHist *data, double *NLL, int *stat_t, int MaxTries, double minMassFit, double maxMassFit){

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
    // fitTest = pdf->fitTo(*data, RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(2), RooFit::PrintLevel(3), RooFit::Warnings(false), RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE), RooFit::Range(minMassFit,maxMassFit));
    // RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Offset(kTRUE),RooFit::Strategy(2));   

    //-------------------
    RooNLLVar *nll=new RooNLLVar("nll","nll",*pdf,*data, RooFit::Range(minMassFit,maxMassFit));
    RooMinimizer *minuit_fitTest = new RooMinimizer(*nll);
    // minuit_fitTest->setOffsetting(kTRUE);
    minuit_fitTest->setStrategy(1);
    minuit_fitTest->setPrintLevel(-1000);
    minuit_fitTest->minimize("Minuit2","minimize");
    std::cout << "Now running hesse" << std::endl;
    hessestatus = minuit_fitTest->hesse();
    std::cout << "Running minos" << std::endl;
    minosstatus = minuit_fitTest->minos();
    //-------------------
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
void asimovDatasetFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &ws_outdir, const std::string &samplename, FILE *resFileasimovfit) {
  int ncat = NCAT;
  RooDataHist* signal;
  // RooVoigtian* voigt[NCAT+1];
  RooDCBShape* cb[NCAT+1];
  TCanvas* canv2 = new TCanvas( Form("Canvas2_M%f_k%s", mass , coupling.c_str()) , Form("Canvas2_M%f_k%s", mass , coupling.c_str()) );
  canv2->cd();
  // TPaveText* label_cms = get_labelcms(0, "2012", true);
  // TPaveText* label_sqrt = get_labelsqrt(0);
  int iMass = (int)abs(mass);
  for(int c = 0; c<ncat+1; c++){
    
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();
    TString myCut;
    
    //pickup the functions from the ws

    // RooRealVar* gH = new RooRealVar("gH", "gH", 0, 1000);
    
    if(c<2)    signal = (RooDataHist*)w->data(TString::Format("signal_asimov_cat%d",c));
    if(c==2)   signal = (RooDataHist*)w->data ("signal_asimov");
   
 
    //crystallball
    //cb pos                                  

    RooRealVar* cb_deltaMean = new RooRealVar(TString::Format("cb_deltaMean_cat%d",c),TString::Format("cb_deltaMean_cat%d",c), 0., -300., 300.);
    RooFormulaVar* cb_mean = new RooFormulaVar(TString::Format("cb_mean_cat%d",c),"@0+@1", RooArgList(*cb_deltaMean, *w->var("MH")));
    RooRealVar* cb_sigma = new RooRealVar(TString::Format("cb_sigma_cat%d",c),TString::Format("cb_sigma_cat%d",c),20, 0., 200);//90
    RooRealVar* cb_alphaL = new RooRealVar(TString::Format("cb_alphaL_cat%d",c),TString::Format("cb_alphaL_cat%d",c),1.5, 0.7, 1.7);
    RooRealVar* cb_alphaR = new RooRealVar(TString::Format("cb_alphaR_cat%d",c),TString::Format("cb_alphaR_cat%d",c),2, 1.5, 4);
    //cb_alphaR->setConstant();
    //cb_alphaL->setConstant();
    RooRealVar* cb_nL = new RooRealVar(TString::Format("cb_nL_cat%d",c),TString::Format("cb_nL_cat%d",c), 3.3, 3, 4.5);
    RooRealVar* cb_nR = new RooRealVar(TString::Format("cb_nR_cat%d",c),TString::Format("cb_nR_cat%d",c), 4., 3, 6);
    // if(coupling!="001") cb_nR->setRange(0.,15);
    // if(coupling!="001") cb_nR->setVal(2);
    // if(coupling!="001")cb_nL->setRange(0., 15);
    // if(coupling!="001")cb_nL->setVal(2);
    // if(coupling!="001") cb_alphaR->setRange(0.,5);
    // if(coupling!="001") cb_alphaR->setVal(2);
    // if(coupling!="001")cb_alphaL->setRange(0., 5);
    // if(coupling!="001")cb_alphaL->setVal(2);
    
    // if(mass>3200&&c==0&&coupling=="001") cb_nR->setRange(4.,5);
    // if(mass>3200&&c==0&&coupling=="001") cb_nR->setVal(4.3);
    // if(mass>3200&&c==0&&coupling=="001")cb_nL->setRange(4.5,6.);
    // if(mass>3200&&c==0&&coupling=="001")cb_nL->setVal(5);
    //if(mass>3200&&c==1)cb_nL->setRange(3,4.5);
    //cb_nR->setConstant();
    //cb_nL->setConstant();

    // if (iMass==750){cb_alphaL->setVal(0.5);cb_alphaR->setVal(0.5);cb_nL->setVal(1);cb_nR->setVal(1.);}
    
    cb[c] =  new RooDCBShape(TString::Format("cbshape_cat%d",c),TString::Format("cbshape_cat%d",c) , *w->var("var"), *cb_mean, *cb_sigma,  *cb_alphaL, *cb_alphaR,  *cb_nL,*cb_nR) ;
    w->import(*cb[c]);
   
    double massMin = 0.; 
    double massMax = 0.;
    if(coupling=="001" || coupling == "0p014"){
      massMin=0.8*mass;
      massMax=1.2*mass;
    }else if(coupling=="01" || coupling == "1p4"){
      massMin=0.6*mass;
      massMax=1.4*mass;
    }else if(coupling=="02" || coupling == "5p6"){
      massMin=0.4*mass;
      massMax=1.6*mass;
      // massMin=0.8*mass;
      // massMax=1.2*mass;
    }
    if(massMin<300)massMin=300;
    if(massMax>9000)massMax=9000;
    // if (iMass==750){massMin=600;massMax=850;}

    //fitting
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    // RooFitResult* fitresults = (RooFitResult* ) cb[c]->fitTo(*signal,SumW2Error(kTRUE), Range(massMin,massMax), RooFit::Save(kTRUE));
    // fitresults->Print("V");

    int fitStatus = 0;
    double thisNll = 0.;
    theFitResult fitresults = theFit(cb[c], signal, &thisNll, &fitStatus, /*max iterations*/ 3, massMin, massMax) ;

    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;

    // fitresults.fitres->Print("V");
  
    //plotting...
    double prob = 0.;
    double chi2 = 0.;
    int nBinsForMass = 160;
    RooDCBShape*  dcb;
    RooPlot* plotg;
    // while(prob<0.01){
    // RooPlot* plotg = (w->var("var"))->frame(Range(massMin,massMax),Title("mass"),Bins(nBinsForMass));
    plotg = (w->var("var"))->frame(Range(massMin,massMax),Title("mass"),Bins(nBinsForMass));
    signal->plotOn(plotg,Name("signal"));
    
    // voigt[c]->plotOn(plotg, LineColor(kBlue));
    // cb[c]->plotOn(plotg, LineColor(kBlue));
    //check parameters 
    // RooArgSet* model_params = voigt[c]->getParameters(*w->var("mgg")) ;
    //model_params->Print("v") ;
    RooArgList model_params = fitresults.fitres->floatParsFinal();
    // RooArgList model_params = fitresults->floatParsFinal();
    model_params.Print("v") ;
    //save a new function in the ws as vittorio said

    RooRealVar* dcb_dm = new RooRealVar(TString::Format("dcb_dm_cat%d",c),TString::Format("dcb_dm_cat%d",c), cb_deltaMean->getVal());
    RooRealVar* dcb_m = new RooRealVar(TString::Format("dcb_m_cat%d",c),TString::Format("dcb_m_cat%d",c), cb_mean->getVal());
    RooRealVar* dcb_s = new RooRealVar(TString::Format("dcb_s_cat%d",c),TString::Format("dcb_s_cat%d",c), cb_sigma->getVal());
    RooRealVar* dcb_aL = new RooRealVar(TString::Format("dcb_aL_cat%d",c),TString::Format("dcb_aL_cat%d",c),cb_alphaL->getVal());
    RooRealVar* dcb_aR = new RooRealVar(TString::Format("dcb_aR_cat%d",c),TString::Format("dcb_aR_cat%d",c),cb_alphaR->getVal());
    RooRealVar* dcb_nL = new RooRealVar(TString::Format("dcb_nL_cat%d",c),TString::Format("dcb_nL_cat%d",c),cb_nL->getVal() );
    RooRealVar* dcb_nR = new RooRealVar(TString::Format("dcb_nR_cat%d",c),TString::Format("dcb_nR_cat%d",c),cb_nR->getVal() );
    dcb_dm->setError(cb_deltaMean->getError());
    dcb_m->setError(cb_deltaMean->getError());
    dcb_s->setError(cb_sigma->getError());
    dcb_aL->setError(cb_alphaL->getError());
    dcb_aR->setError(cb_alphaR->getError());
    dcb_nL->setError(cb_nL->getError());
    dcb_nR->setError(cb_nR->getError());
    // RooDCBShape*  dcb =  new RooDCBShape(TString::Format("dcbshape_cat%d",c),TString::Format("dcbshape_cat%d",c) , *w->var("var"), *dcb_m, *dcb_s,  *dcb_aL, *dcb_aR,  *dcb_nL,*dcb_nR) ;
    dcb =  new RooDCBShape(TString::Format("dcbshape_cat%d",c),TString::Format("dcbshape_cat%d",c) , *w->var("var"), *dcb_m, *dcb_s,  *dcb_aL, *dcb_aR,  *dcb_nL,*dcb_nR) ;
    w->import(*dcb);
    w->import(*dcb_dm);
    dcb->plotOn(plotg, LineColor(kBlue), Name("pdf"));

    RooRealVar norm("norm","norm",signal->sumEntries(),0,10E6);
    RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*dcb,norm);
    int np = pdf->getParameters(*signal)->getSize();
    chi2 = plotg->chiSquare("pdf","signal",np);

    std::cout << "[INFO] Calculating GOF for pdf " << dcb->GetName() << ", using " <<np << " fitted parameters" <<std::endl;
    prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np);
    //chi2/ndof is chi2
    std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
    std::cout << "[INFO] p-value  =  " << prob << std::endl;
    // }

    RooArgSet* model_params2 = dcb->getParameters(*w->var("var")) ;
    model_params2->Print("v") ;

    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma}[GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    TLegend* legmc = new TLegend(0.58, 0.54, 0.91, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->AddEntry(plotg->getObject(0),"Asimov Data","lp");    
    //  legmc->AddEntry(plotg->getObject(1),"Voigtian Fit ","l");
    legmc->AddEntry(plotg->getObject(1),"DCB Fit ","l");
    
    // bw[c]->plotOn(plotg, LineColor(kBlue));
    plotg->GetYaxis()->SetRangeUser(0.1,30000 );
    plotg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    plotg->GetXaxis()->SetTitleFont(42);
    // plotg->GetXaxis()->SetTitlesize(0.07);
    // plotg->GetXaxis()->SetTitleoffset(1.40);
    plotg->GetYaxis()->SetTitle("arbitrary scale");
    plotg->GetYaxis()->SetTitleFont(42);
    plotg->GetYaxis()->SetTitleSize(0.04);
    plotg->Draw();
    legmc->Draw("same");
    canv2->SetLogy(0);
    
    canv2->SaveAs(TString::Format("%s/../FinalConvShapes/%s_fitasimov_cat%d_M%d_k_%s.png", ws_outdir.c_str(), samplename.c_str(), c,iMass,coupling.c_str())); 

    canv2->SetLogy();
   
    canv2->SaveAs(TString::Format("%s/../FinalConvShapes/%s_fitasimov_cat%d_M%d_k_%s_log.png", ws_outdir.c_str(), samplename.c_str(), c,iMass,coupling.c_str())); 

    canv2->SetLogy(0);  

    //Let's save the results to file
    float coup = widthtonum(coupling);
    fprintf(resFileasimovfit,"%s & %d & %10.1e & %d & %10.2f & %10.2f & (%d,%d,%d) \\\\\n", samplename.c_str() , iMass, coup, c, chi2, prob, fitresults.minimizestatus ,fitresults.hessestatus ,fitresults.minosstatus);

    // fprintf(resFileasimovfit,"%s & %d & %10.1e & %d & %10.2f & %10.2f & (%d,%d,%d) \\\\\n", samplename.c_str() , iMass, coup, c, chi2, prob, -1 , -1 , -1);

  }
		 
}
//-----------------------------------------------------------------------------------
RooDataHist* throwAsimov( double nexp, RooAbsPdf *pdf, RooRealVar *x, RooDataHist *asimov)
{
  RooBinning binning(360, 300, 9000, "binning");
    RooArgSet mset( *x );
    if( asimov != 0 ) {
        asimov->reset();
    } else {
      asimov = new RooDataHist(Form("asimov_dataset_%s",pdf->GetName()),Form("asimov_dataset_%s",pdf->GetName()),mset, "binning");
    }
    pdf->fillDataHist( asimov, &mset, 1, false, true );

    for( int i = 0 ; i < asimov->numEntries() ; i++ ) {
        asimov->get( i ) ;

        // Expected data, multiply p.d.f by nEvents
        Double_t w = asimov->weight() * nexp;
        asimov->set( w, sqrt( w ) );
	std::cout<<i<<" "<<w<<std::endl;
    }
    std::cout<<nexp<<" "<<asimov->sumEntries()<<std::endl;
    Double_t corr = nexp / asimov->sumEntries();
    for( int i = 0 ; i < asimov->numEntries() ; i++ ) {
        RooArgSet theSet = *( asimov->get( i ) );
        asimov->set( asimov->weight()*corr, sqrt( asimov->weight()*corr ) );
    }
    
    return asimov;
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

