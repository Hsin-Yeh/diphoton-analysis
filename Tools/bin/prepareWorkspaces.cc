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
RooRealVar* buildRooVar(std::string name, std::string title, int nBins, double xMin, double xMax, std::string unit);
void runAllFits(const std::string &year, const std::string &ws_dir);
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample);
void AddSigData(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &isample);
void sigModelResponseFcnFit(RooWorkspace* w,Float_t mass, std::string coupling, const std::string &year);
void sigModelGenFcnFit(RooWorkspace* w,Float_t mass,std::string coupling, const std::string &year);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
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
    std::cout << "Syntax: makeParametricFits.exe [2016/2017/2018] [input]" << std::endl;
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
  runAllFits(year,"input/workspaces");

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void runAllFits(const std::string &year, const std::string &ws_dir){

  std::vector<std::string> samples = getSampleListForFit();

  for(auto isample : samples) {
    //We will process one year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    runfits(year, ws_dir, isample);
  }

}
//-----------------------------------------------------------------------------------
void AddSigData(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &isample) {

  Int_t ncat = NCAT;//BB and BE for the moment

  RooRealVar* mgg = buildRooVar("mgg", "m_{#gamma#gamma}", 2385, 230., 10000., "GeV");
  RooRealVar* mggGen = buildRooVar("mggGen", "m_{#gamma#gamma} gen", 2385, 230., 10000, "GeV");
  RooRealVar* weight = buildRooVar("weight", "event weight", 0, 0 , 1000 , "nounits");
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
void sigModelResponseFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year) {
  Int_t ncat = NCAT;
  RooDataSet* signal;
  RooDCBShape* responseadd[NCAT+1];
  int iMass = (int)  abs(mass);
  TCanvas* canv2 = new TCanvas( Form("Canvas2_M%f_k%s", mass , coupling.c_str()) , Form("Canvas2_M%f_k%s", mass , coupling.c_str()) );
  canv2->cd();
  TPaveText* label_cms = get_labelcms(0, year , true);
  TPaveText* label_sqrt = get_labelsqrt(0);

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
    //cb pos                                                                                                                     
    RooFormulaVar cbpos_mean(TString::Format("reducedmass_cbpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("reducedmass_sig_mean_cat%d",c)));
    RooFormulaVar cbpos_sigma(TString::Format("reducedmass_cbpos_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("reducedmass_sig_sigma_cat%d",c)));
    RooFormulaVar cbpos_alphacb(TString::Format("reducedmass_cbpos_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbpos_cat%d",c)));
    RooFormulaVar cbpos_n(TString::Format("reducedmass_cbpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_npos_cat%d",c)));     
    //cb neg
    RooFormulaVar cbneg_n(TString::Format("reducedmass_cbneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_nneg_cat%d",c)));
    RooFormulaVar cbneg_alphacb(TString::Format("reducedmass_cbneg_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("reducedmass_sig_alphacbneg_cat%d",c)));   
   
    responseadd[c]= new  RooDCBShape(TString::Format("responseaddpdf_cat%d",c),TString::Format("responseaddpdf_cat%d",c) , *w->var("massReduced"), cbpos_mean, cbpos_sigma,  cbneg_alphacb, cbpos_alphacb,  cbneg_n,cbpos_n) ;

    w->import(*responseadd[c]);
    bool fixPar = false;
    if(mass==2000 && fixPar){
      w->var(TString::Format("reducedmass_sig_mean_cat%d",c))->setVal(-3.);
      w->var(TString::Format("reducedmass_sig_sigma_cat%d",c))->setVal(25);
      w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c))->setVal(1);
      w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c))->setVal(1);
      w->var(TString::Format("reducedmass_sig_npos_cat%d",c))->setVal(5);
      w->var(TString::Format("reducedmass_sig_nneg_cat%d",c))->setVal(5);
    }else if(mass==1500 && fixPar){
      w->var(TString::Format("reducedmass_sig_mean_cat%d",c))->setVal(-12.);
      w->var(TString::Format("reducedmass_sig_sigma_cat%d",c))->setVal(25);
      w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c))->setVal(-2);
      w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c))->setVal(1);
      w->var(TString::Format("reducedmass_sig_npos_cat%d",c))->setVal(3);
      w->var(TString::Format("reducedmass_sig_nneg_cat%d",c))->setVal(7);
    }
    RooFitResult* fitresults = (RooFitResult* ) responseadd[c]->fitTo(*signal,RooFit::SumW2Error(kTRUE), RooFit::Range(-600, 600), RooFit::Save(kTRUE));
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   
    RooRealVar* massreduced = (RooRealVar*) w->var("massReduced");
    RooPlot* plotg = massreduced->frame(Range(-600., 600.),Title("mass reduced"),Bins(120));
    signal->plotOn(plotg);
   

    responseadd[c]->plotOn(plotg, LineColor(kBlue));
    // responseadd[c]->plotOn(plotg,Components(TString::Format("responsecbneg_cat%d",c)), LineColor(kGreen), LineStyle(kDotted));
    //responseadd[c]->plotOn(plotg,Components(TString::Format("responsecbpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
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


    canv2->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/responsefcnfitcbcb_cat%d_M%d_k_%s.png",c,iMass,coupling.c_str())); 
    canv2->SetLogy();
    
    canv2->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/responsefcnfitcbcb_cat%d_M%d_k_%s_log.png",c,iMass,coupling.c_str())); 

    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    canv2->SetLogy(0);  
	
    w->defineSet(TString::Format("responseaddpdfparam_cat%d",c),RooArgSet(*w->var(TString::Format("reducedmass_sig_sigma_cat%d",c)), 
    									  *w->var(TString::Format("reducedmass_sig_alphacbpos_cat%d",c)),
    									  *w->var(TString::Format("reducedmass_sig_alphacbneg_cat%d",c)),
    									  *w->var(TString::Format("reducedmass_sig_npos_cat%d",c)),
    									  *w->var(TString::Format("reducedmass_sig_nneg_cat%d",c)),	   
    									  *w->var(TString::Format("reducedmass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("responseaddpdfparam_cat%d",c)));
  }//end of loop over categories
  
}

//-----------------------------------------------------------------------------------
void sigModelGenFcnFit(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &year) {
  int ncat = NCAT;
  RooDataSet* signal;
  //  RooBreitWigner* bw[NCAT+1];
  RooGenericPdf* bw[NCAT+1];
  RooDCBShape* mgenadd[NCAT+1];
  int iMass =  (int) abs(mass);
  TCanvas* canv3 = new TCanvas( Form("Canvas3_M%f_k%s", mass , coupling.c_str()) , Form("Canvas3_M%f_k%s", mass , coupling.c_str()) );
  canv3->cd();
  TPaveText* label_cms = get_labelcms(0, "2016", true);
  TPaveText* label_sqrt = get_labelsqrt(0);
  
  for(int c = 0; c<ncat+1; c++){
    
    TLatex *lat  = new TLatex(0.6,0.9,TString::Format("cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

    TString myCut;
    signal = (RooDataSet*)w->data("SigWeightGen");
    w->var(TString::Format("mgen_sig_mean_cat%d",c))->setVal(mass);
    w->var(TString::Format("mgen_sig_mean_cat%d",c))->setRange(mass*0.8, mass*1.2);
    //cb pos                                                                                                                     
    RooFormulaVar cbpos_mean(TString::Format("mgen_cbpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("mgen_sig_mean_cat%d",c)));
    RooFormulaVar cbpos_sigma(TString::Format("mgen_cbpos_sig_sigma_cat%d",c), "", "sqrt(@0*@0)", *w->var(TString::Format("mgen_sig_sigma_cat%d",c)));
    RooFormulaVar cbpos_alphacb(TString::Format("mgen_cbpos_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbpos_cat%d",c)));
    RooFormulaVar cbpos_n(TString::Format("mgen_cbpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_npos_cat%d",c)));    
    //cb neg
    RooFormulaVar cbneg_n(TString::Format("mgen_cbneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_nneg_cat%d",c)));
    RooFormulaVar cbneg_alphacb(TString::Format("mgen_cbneg_sig_alphacb_cat%d",c),"", "@0", *w->var( TString::Format("mgen_sig_alphacbneg_cat%d",c)));    
    mgenadd[c]=  new  RooDCBShape(TString::Format("mgenaddpdf_cat%d",c),TString::Format("mgenaddpdf_cat%d",c) , *w->var("mggGen"), cbpos_mean, cbpos_sigma,  cbneg_alphacb, cbpos_alphacb,  cbneg_n,cbpos_n) ;
    w->import(*mgenadd[c]);
    std::cout<<mass<<" "<<signal->sumEntries()<<std::endl;
    std::cout<<TString::Format("******************************** signal fit results cb+cb  mass %f cat %d***********************************", mass, c)<<std::endl;
    RooFitResult* fitresults = (RooFitResult* ) mgenadd[c]->fitTo(*signal,SumW2Error(kTRUE), Range(mass*0.8, mass*1.2), RooFit::Save(kTRUE));   
    fitresults->Print("V");

  
    RooPlot* plotg = w->var("mggGen")->frame(Range(mass*0.8,mass*1.2),Title("mass generated"),Bins(120));
    signal->plotOn(plotg);
    mgenadd[c]->plotOn(plotg, LineColor(kBlue));
   
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


    canv3->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/mgen_cat%d_M%d_k_%s.png",c,iMass,coupling.c_str())); 
    canv3->SetLogy();
    
    canv3->SaveAs(TString::Format("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/mgen_cat%d_M%d_k_%s_log.png",c,iMass,coupling.c_str())); 
    plotg->GetYaxis()->SetRangeUser(0.0001,plotg->GetMaximum()*0.12 );
    canv3->SetLogy(0);  
	
    w->defineSet(TString::Format("mgenpdfparam_cat%d",c),RooArgSet(*w->var(TString::Format("mgen_sig_sigma_cat%d",c)), 
								   *w->var(TString::Format("mgen_sig_alphacbpos_cat%d",c)),
								   *w->var(TString::Format("mgen_sig_alphacbneg_cat%d",c)),
								   *w->var(TString::Format("mgen_sig_npos_cat%d",c)),
								   *w->var(TString::Format("mgen_sig_nneg_cat%d",c)),	   
								   *w->var(TString::Format("mgen_sig_frac_cat%d",c)),  
								   *w->var(TString::Format("mgen_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("mgenpdfparam_cat%d",c)));
  }

}
//-----------------------------------------------------------------------------------
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample){

  std::string coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
  std::string M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");
  std::cout << ws_dir << std::endl;

  TFile *fout = new TFile(Form("%s/ws_ResponseAndGen_M%s_k%s_%s.root", ws_dir.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()), "recreate");

  TString fileBaseName("HighMassGG");    
  TString fileBkgName("HighMassGG.inputbkg");
  HLFactory hlf("HLFactory", "HighMassGG.rs", false);
  RooWorkspace* w = hlf.GetWs();
  // range for the variables
  w->var("mgg")->setMin(MINmass);
  w->var("mgg")->setMax(MAXmass);
  w->var("mggGen")->setMin(MINmass);
  w->var("mggGen")->setMax(MAXmass);
  w->Print("V");

  //========================================================================
  std::cout << "Adding signal data for mass " << M_bins << " and coupling " << coupling <<std::endl;
  AddSigData(w, stof(M_bins), coupling, isample);
  //========================================================================

  TCanvas* canv1 = new TCanvas(Form("Canvas1_%s",isample.c_str()),Form("Canvas_%s",isample.c_str()));
  canv1->cd();
  RooPlot* p = w->var("mgg")->frame(RooFit::Range(230., 10000.), RooFit::Bins( 2385));

  //========================================================================
  std::cout << "Fit the response function" <<std::endl;
  sigModelResponseFcnFit(w, stof(M_bins), coupling, year);
  //========================================================================
  std::cout << "Fit the gen function" <<std::endl;
  sigModelGenFcnFit(w, stof(M_bins), coupling, year);
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

