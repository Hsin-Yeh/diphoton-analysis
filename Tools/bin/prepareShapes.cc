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

//-----------------------------------------------------------------------------------
//Declarations here definition after main
void runAllFits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir);
void runfits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &isample);
void sigModelShapeFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year, bool doPlot);
void asimovDatasetFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &yeart);
RooDataHist* throwAsimov( double nexp, RooAbsPdf *pdf, RooRealVar *x, RooDataHist *asimov );
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);
void SetConstantParams(const RooArgSet* params) ;
TPaveText* get_labelsqrt( int legendquadrant );
TPaveText* get_labelcms( int legendquadrant, std::string year, bool sim);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir;

  if(argc!=3) {
    std::cout << "Syntax: prepareShapes.exe [2016/2017/2018] [input]" << std::endl;
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
  runAllFits(year,"input/workspaces", "output/workspaces");

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void sigModelShapeFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year, bool doPlot) {
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
  
  for(int c = 0; c<ncat+1; c++){
    std::cout << "FINAL SHAPE FIT FOR CATEGORY " << c<< std::endl; 
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();
    TString myCut;
      
    //pickup the functions from the ws
    double massMin = 0.; 
    double massMax = 0.;
    if(coupling=="001"){
      massMin=0.8*mass;
      massMax=1.2*mass;
    }else if(coupling=="01"){
      massMin=0.6*mass;
      massMax=1.4*mass;
    }else if(coupling=="02"){
      massMin=0.4*mass;
      massMax=1.6*mass;
    }
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
    if(massMax>6000)massMax=6000;
    
    RooBinning binning(160, massMin, massMax, "binning");
    TH1F* dummy = new TH1F("dummy", "dummy", 160, massMin, massMax);
    RooDataHist* datahist_asimov = new RooDataHist("datahist_asimov","datahist_asimov",*var, dummy, 1);
    // RooDataHist* datahist_asimov_clone = new RooDataHist("datahist_asimov","datahist_asimov",*w->var("mgg"), dummy, 1);
      
    //  RooDataHist* datahist_asimov; //= mgen[c]->generateBinned(RooArgSet(*var),RooFit::NumEvents(signal->numEntries()), RooFit::Asimov());
    //datahist_asimov = signal->binnedClone();
    throwAsimov(10000, conv[c], var,datahist_asimov); 
   
    if(c<2)    datahist_asimov->SetName(TString::Format("signal_asimov_cat%d",c));
    if(c==2)    datahist_asimov->SetName("signal_asimov");
    std::cout<<"--------------------------------------------------- "<<datahist_asimov->sumEntries()<<std::endl;
    w->import(*datahist_asimov);

    bool doPlot=false;
     
    if(doPlot){
      //plotting...
      if(c==0||c==1)signalK = (RooDataSet*)w->data(TString::Format("SigWeightK_cat%d",c));
      if(c==2)signalK = (RooDataSet*)w->data("SigWeightK_cat0");
      std::cout<<"flag1"<<std::endl;
      w->Print("v");
      std::cout<<"flag2"<<std::endl;
      RooPlot* plotgg = (w->var("mgg"))->frame(Range(massMin,massMax),Title("mass generated"),Bins(160));
      signalK->Print("v");
      signalK->plotOn(plotgg);
      std::cout<<"flag3"<<std::endl;
      RooPlot* plotg = var->frame(Range(massMin,massMax),Title("mass generated"),Bins(160));
      RooDataSet* fakedata= conv[c]->generate(*var, signalK ->sumEntries());
      fakedata->plotOn(plotg, RooFit::Invisible());
      // reducedmass[c]->plotOn(plotg, LineColor(kRed));
      mgen[c]->plotOn(plotg, LineColor(kGreen));
      conv[c]->plotOn(plotg, LineColor(kBlue));
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
    

      canv1->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalConvShapes/mreco_cat%d_M%d_k_%s.png",c,iMass,coupling.c_str())); 

      canv1->SetLogy();			        
    					        
      canv1->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalConvShapes/mreco_cat%d_M%d_k_%s_log.png",c,iMass,coupling.c_str())); 

      canv1->SetLogy(0);  
    }
    
  }
 
}

//-----------------------------------------------------------------------------------
void runAllFits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir){

  std::vector<std::string> samples = getSampleListForFit();

  for(auto isample : samples) {
    //We will process one year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    runfits(year, ws_indir, ws_outdir, isample);
  }

}

//-----------------------------------------------------------------------------------
void runfits(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &isample){

  std::string coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
  std::string M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");
  std::cout << "Input dir " << ws_indir << std::endl;
  std::cout << "Output dir " << ws_outdir << std::endl;

  TFile *fresandgen = TFile::Open(Form("%s/ws_ResponseAndGen_M%s_k%s_%s.root", ws_indir.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()));

  RooWorkspace* ws = (RooWorkspace*) fresandgen->Get("HLFactory_ws"); 
  //========================================================================
  std::cout << "Make and fit the convolution of the response and the gen level shape" <<std::endl;
  bool doPlot = true; //For the final plot
  sigModelShapeFcnFit(ws, stof(M_bins), coupling, year, doPlot);
  //========================================================================
  std::cout << "Throw an asimov dataset from the DCB" <<std::endl;
  asimovDatasetFcnFit(ws, stof(M_bins), coupling, year);

  //========================================================================
  //Save the resulting workspace
  TFile *fout = new TFile(Form("%s/ws_ResponseAndGen_M%s_k%s_%s_final.root", ws_outdir.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()), "recreate");
  fout->cd();
  // ws->Print("V");
  ws->Write();
  fout->Write();
  fout->Close();

}

//-----------------------------------------------------------------------------------
void asimovDatasetFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &yeart) {
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
    if(coupling!="001") cb_nR->setRange(0.,15);
    if(coupling!="001") cb_nR->setVal(2);
    if(coupling!="001")cb_nL->setRange(0., 15);
    if(coupling!="001")cb_nL->setVal(2);
    if(coupling!="001") cb_alphaR->setRange(0.,5);
    if(coupling!="001") cb_alphaR->setVal(2);
    if(coupling!="001")cb_alphaL->setRange(0., 5);
    if(coupling!="001")cb_alphaL->setVal(2);
    
    if(mass>3200&&c==0&&coupling=="001") cb_nR->setRange(4.,5);
    if(mass>3200&&c==0&&coupling=="001") cb_nR->setVal(4.3);
    if(mass>3200&&c==0&&coupling=="001")cb_nL->setRange(4.5,6.);
    if(mass>3200&&c==0&&coupling=="001")cb_nL->setVal(5);
    //if(mass>3200&&c==1)cb_nL->setRange(3,4.5);
    //cb_nR->setConstant();
    //cb_nL->setConstant();
    cb[c] =  new RooDCBShape(TString::Format("cbshape_cat%d",c),TString::Format("cbshape_cat%d",c) , *w->var("var"), *cb_mean, *cb_sigma,  *cb_alphaL, *cb_alphaR,  *cb_nL,*cb_nR) ;
    w->import(*cb[c]);
   
    //fitting
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    RooFitResult* fitresults = (RooFitResult* ) cb[c]->fitTo(*signal,SumW2Error(kTRUE), Range(mass*0.7,mass*1.3), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
  
    //plotting...
    double massMin = 0.; 
    double massMax = 0.;
    if(coupling=="001"){
      massMin=0.8*mass;
      massMax=1.2*mass;
    }else if(coupling=="01"){
      massMin=0.6*mass;
      massMax=1.4*mass;
    }else if(coupling=="02"){
      massMin=0.4*mass;
      massMax=1.6*mass;
    }
    if(massMin<300)massMin=300;
    if(massMax>6000)massMax=6000;
    RooPlot* plotg = (w->var("var"))->frame(Range(massMin,massMax),Title("mass"),Bins(60));
    signal->plotOn(plotg);
    
    // voigt[c]->plotOn(plotg, LineColor(kBlue));
    // cb[c]->plotOn(plotg, LineColor(kBlue));
    //check parameters 
    // RooArgSet* model_params = voigt[c]->getParameters(*w->var("mgg")) ;
    //model_params->Print("v") ;
    RooArgList model_params = fitresults->floatParsFinal();
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
    RooDCBShape*  dcb =  new RooDCBShape(TString::Format("dcbshape_cat%d",c),TString::Format("dcbshape_cat%d",c) , *w->var("var"), *dcb_m, *dcb_s,  *dcb_aL, *dcb_aR,  *dcb_nL,*dcb_nR) ;
    w->import(*dcb);
    w->import(*dcb_dm);
    dcb->plotOn(plotg, LineColor(kBlue));
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
    

    canv2->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalConvShapes/fitasimov_cat%d_M%d_k_%s.png",c,iMass,coupling.c_str())); 

    canv2->SetLogy();
   
    canv2->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalConvShapes/fitasimov_cat%d_M%d_k_%s_log.png",c,iMass,coupling.c_str())); 

    canv2->SetLogy(0);  
  }
		 
}
//-----------------------------------------------------------------------------------
RooDataHist* throwAsimov( double nexp, RooAbsPdf *pdf, RooRealVar *x, RooDataHist *asimov )
{
  RooBinning binning(360, 300, 6000, "binning");
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

