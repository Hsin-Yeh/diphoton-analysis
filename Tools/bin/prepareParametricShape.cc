#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
#include <string>  

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
#include "TGraphErrors.h"

using namespace RooFit;
using namespace RooStats;

//-----------------------------------------------------------------------------------
struct theTree {

  std::string kMpl;
  std::string M_bins;
  std::string name;
  TFile * fresandgen_final;
  RooWorkspace * ws;

};

//-----------------------------------------------------------------------------------
static const Int_t NCAT = 2; //BB and BE for the moment 

//-----------------------------------------------------------------------------------
//Declarations here definition after main
double computePdfFHWM(RooDCBShape pdf,RooRealVar roobs, double MH, bool plot, std::string name);
void plotAllSignalsAsimov(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &couplingIn);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir, thecoupling;

  if(argc!=4) {
    std::cout << "Syntax: prepareParametricShape.exe [2016/2017/2018] [input] [kMplxxx]" << std::endl;
      return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    inputdir = argv[2];
    thecoupling = argv[3];
 }

  //========================================================================
  //read the reduced input root trees
  initForFit(inputdir);
  //========================================================================
  //Run the fits
  plotAllSignalsAsimov(year,"output/workspaces", "output/FinalParametricShape/workspaces", thecoupling);

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
double computePdfFHWM(RooDCBShape pdf,RooRealVar roobs, double MH, bool plot, std::string name){
  TCanvas* ccc = new TCanvas(Form("ccc_%s", name.c_str()), Form("ccc_%s", name.c_str()));
  ccc->cd();
  int nBins = 1000;
  double mean = MH;
  double  sigma  = mean*0.2;
  TH1F*  hist = (TH1F*) pdf.createHistogram(Form("sigHist_%s",name.c_str()),roobs, RooFit::Binning(nBins,mean-4.*sigma,mean+4.*sigma) );
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
  if(plot){
    hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetBinCenter(maxBin)-5*xWidth,hist->GetXaxis()->GetBinCenter(maxBin)+5*xWidth);
    hist->Draw("HIST");
       
    TLine* xL = new TLine(xLeft,0,xLeft,halfMaxVal);
    TLine* xR =  new TLine(xRight,0,xRight,halfMaxVal);
    xL->SetLineColor(kRed);
    xR->SetLineColor(kRed);
    xL->Draw("same");  
    xR->Draw("same");
    ccc->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/test_%s.png",name.c_str()));
  }

  return xWidth;
}


//-----------------------------------------------------------------------------------
void plotAllSignalsAsimov(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &couplingIn){

  std::vector<std::string> samples = getSampleListForFit();
  bool plot = true;

  //========================================================================
  std::cout << "Calculate the final parametric shape" <<std::endl;
  std::cout << "Input dir " << ws_indir << std::endl;
  std::cout << "Output dir " << ws_outdir << std::endl;
  std::string coupling = "";
  std::string M_bins = "";

  double upperxmax = 0.; 
  if ( couplingIn == "kMpl001" ){upperxmax = 6000 ;}
  else if ( couplingIn == "kMpl01" ){upperxmax = 9000 ;}
  else if ( couplingIn == "kMpl02" ){upperxmax = 9000 ;}
  else {
    std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
    exit(1);
  }

  RooDCBShape final_shape[5];

  std::vector<theTree> theInput; 

  for(auto isample : samples) {

    //BE CAREFUL
    //We will process one year each time and one coupling each time
    if ( isample.find(year) == std::string::npos ) continue; 
    if ( isample.find(couplingIn) == std::string::npos ) continue; 
    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;

    theTree tmpin;

    coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");

    tmpin.kMpl = coupling;
    tmpin.M_bins = M_bins;
    tmpin.name = getBase(isample);
    tmpin.fresandgen_final = TFile::Open(Form("%s/ws_ResponseAndGen_M%s_k%s_%s_final.root", ws_indir.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()));  
    tmpin.ws = (RooWorkspace*) tmpin.fresandgen_final->Get("HLFactory_ws");
 
    theInput.push_back(tmpin);

  }
  
  double m[theInput.size()];
  double mErr[theInput.size()];
  double s[theInput.size()];
  double sErr[theInput.size()];
  double aL[theInput.size()];
  double aLErr[theInput.size()];
  double aR[theInput.size()];
  double aRErr[theInput.size()];
  double nL[theInput.size()];
  double nLErr[theInput.size()];
  double nR[theInput.size()];
  double nRErr[theInput.size()];
  double masses[theInput.size()];
  double massesErr[theInput.size()];

  TString svar = "responseaddpdf";
  TString spdf = "responseaddpdf";
  
  RooRealVar* MH = new RooRealVar("MH", "MH", 300, upperxmax);
  MH->setConstant();
  //RooWorkspace* ws_out = new RooWorkspace( Form("%s/model_signal", ws_outdir.c_str()) );
  // RooWorkspace* ws_out = new RooWorkspace( "model_signal" );
  RooWorkspace* ws_out = new RooWorkspace( "ws_inputs" );
  RooDCBShape* sigshape[NCAT+1];
  TF1* fm[NCAT+1];
  TF1* fs[NCAT+1];
  TF1* faL[NCAT+1];
  TF1* faR[NCAT+1];
  TF1* fnL[NCAT+1];
  TF1* fnR[NCAT+1];
  TF1* ffhmw[NCAT+1];

  

  for(int c =0; c< NCAT+1; c++){

    std::map< std::string , RooDCBShape* > res;
    std::map< std::string , TH1F* > h;
    std::map< std::string , RooDataHist* > resdata; 
    std::map< std::string , RooPlot* > pres; 
    std::map< std::string , RooArgSet* > model_params;
    double fhwm[theInput.size()];

    TLegend* legmc = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);

    TCanvas* cc1 = new TCanvas(Form("cc1_cat%d", c), Form("cc1_cat%d", c));

    for(unsigned int iM =0; iM < theInput.size(); iM++){

      theTree tmpin = theInput[iM];
      masses[iM] = std::stod(tmpin.M_bins); 
      massesErr[iM] = 0.;

      res[tmpin.name] = (RooDCBShape*) tmpin.ws->pdf(TString::Format("dcbshape_cat%d",c));
      //compute FWHM
      fhwm[iM] = computePdfFHWM(*res[tmpin.name],*tmpin.ws->var("var"), std::stoi(tmpin.M_bins), plot, Form("%d_%d",iM,c));
      std::cout << "FWHM " << fhwm[iM] <<" for mass " << tmpin.M_bins << std::endl;

      cc1->cd();
      h[tmpin.name] = new TH1F(Form("h%d_%d",iM,c), Form("h%d_%d",iM,c),460, 300,upperxmax);

      h[tmpin.name]->SetLineColor(iM+1);
      h[tmpin.name]->SetMarkerColor(iM+1);

      legmc->AddEntry( h[tmpin.name] , Form("M_{X} = %d GeV", std::stoi(tmpin.M_bins) ),"pl" );

      if(c<2){  resdata[tmpin.name] = (RooDataHist*)tmpin.ws->data(TString::Format("signal_asimov_cat%d",c));}
      if(c==2){ resdata[tmpin.name] = (RooDataHist*)tmpin.ws->data("signal_asimov"); }

      pres[tmpin.name] = tmpin.ws->var("var")->frame(Range(300,upperxmax),Title("mass asimov"),Bins(420));
      resdata[tmpin.name]->plotOn(pres[tmpin.name],MarkerColor(iM+1),LineColor(iM+1));
      res[tmpin.name]->plotOn(pres[tmpin.name],LineColor(iM+1));

      pres[tmpin.name]->GetXaxis()->SetTitle("#Delta m [GeV]");
      pres[tmpin.name]->GetYaxis()->SetTitle("a.u.");
      // pres[tmpin.name]->GetYaxis()->SetRangeUser(0.001,upperxmax);

      if (iM == 0){
	pres[tmpin.name]->Draw();
      } else{
	pres[tmpin.name]->Draw("same");
      }

      legmc->Draw("same");
   
      //save parameters 
      model_params[tmpin.name] = res[tmpin.name]->getParameters(*tmpin.ws->var("var")) ;
      model_params[tmpin.name]->Print("v") ;

      m[iM] = ((RooRealVar*) tmpin.ws->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)tmpin.ws->var("MH"))->getVal();
      mErr[iM] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
      std::cout<<m[iM]<<" "<<mErr[iM]<<std::endl;
      
      s[iM] = ((RooRealVar*) tmpin.ws->var(TString::Format("dcb_s_cat%d",c)))->getVal();
      sErr[iM] =((RooRealVar*) tmpin.ws->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
      std::cout<<s[iM]<<" "<<sErr[iM]<<std::endl;

      aL[iM] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
      aLErr[iM] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
      std::cout<<aL[iM]<<" "<<aLErr[iM]<<std::endl;
      
      aR[iM] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
      aRErr[iM] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
      std::cout<<aR[iM]<<" "<<aRErr[iM]<<std::endl;
      
      nR[iM] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
      nRErr[iM] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
      std::cout<<nR[iM]<<" "<<nRErr[iM]<<std::endl;

      nL[iM] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
      nLErr[iM] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
      std::cout<<nL[iM]<<" "<<nLErr[iM]<<std::endl;

   } //end of loop over signal samples of the same coupling
   
   
   cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/asimovAllMasses_%s_cat%d.png",couplingIn.c_str(), c));
   cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/asimovAllMasses_%s_cat%d.root",couplingIn.c_str(), c));

   // cc1->SetLogy();
   // cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/asimovAllMasses_%s_cat%d_log.png",couplingIn.c_str(), c));

   // cc1->SetLogy(0);
 
   TGraph* gr = new TGraph(theInput.size(),masses,fhwm);
   TGraphErrors* gm = new TGraphErrors(theInput.size(), masses, m, massesErr, mErr);
   TGraphErrors* gs = new TGraphErrors(theInput.size(), masses, s, massesErr, sErr);
   TGraphErrors* gaL = new TGraphErrors(theInput.size(), masses, aL, massesErr, aLErr);
   TGraphErrors* gaR = new TGraphErrors(theInput.size(), masses, aR, massesErr, aRErr);
   TGraphErrors* gnR = new TGraphErrors(theInput.size(), masses, nR, massesErr, nRErr);
   TGraphErrors* gnL = new TGraphErrors(theInput.size(), masses, nL, massesErr, nLErr);

   ffhmw[c] = new TF1(TString::Format("ffhmw_cat%d",c), "pol1", 500,upperxmax);
   gr->Fit(TString::Format("ffhmw_cat%d",c), "R");
   gr->GetXaxis()->SetTitle("m_X[GeV]");
   gr->GetYaxis()->SetTitle("FHWM [GeV]");
   gr->Draw("AP");   
   cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/FHWM_%s_cat%d.png",couplingIn.c_str(), c));


   fm[c] = new TF1(TString::Format("fm_cat%d",c), "pol2", 500,upperxmax);
   gm->Fit(TString::Format("fm_cat%d",c), "R");
   gm->GetYaxis()->SetTitle("#Delta m= m - m_{H} [GeV]");
   gm->GetXaxis()->SetTitle("m_{X}[GeV]");
   gm->Draw("APE");
   fm[c]->Draw("same");
   cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/meanVsMass_%s_cat%d.png",couplingIn.c_str(), c));

   
    gs->GetYaxis()->SetTitle("#sigma [GeV]");
    gs->GetXaxis()->SetTitle("m_{X}[GeV]");
    fs[c] = new TF1(TString::Format("fs_cat%d",c), "pol2", 500,upperxmax);
    gs->Fit(TString::Format("fs_cat%d",c), "R");
    gs->Draw("APE");
    fs[c]->Draw("same");
    cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/sigmaVsMass_%s_cat%d.png",couplingIn.c_str(), c));
    

    gaL->GetYaxis()->SetTitle("#alpha_{L} [GeV]");
    gaL->GetXaxis()->SetTitle("m_{X}[GeV]");
    faL[c] = new TF1(TString::Format("faL_cat%d",c), "pol2", 500,upperxmax);
    gaL->Fit(TString::Format("faL_cat%d",c), "R");
    gaL->Draw("APE");
    faL[c]->Draw("same");
    cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/aLVsMass_%s_cat%d.png",couplingIn.c_str(), c));

  
    gaR->GetYaxis()->SetTitle("#alpha_{R} [GeV]");
    gaR->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    faR[c] = new TF1(TString::Format("faR_cat%d",c), "pol2", 500,upperxmax);
    gaR->Fit(TString::Format("faR_cat%d",c), "R");
    gaR->Draw("APE");
    faR[c]->Draw("same");
    cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/aRVsMass_%s_cat%d.png",couplingIn.c_str(), c));
    
 
    gnR->GetYaxis()->SetTitle("n_{R} [GeV]");
    gnR->GetXaxis()->SetTitle("m_{X}[GeV]");
    fnR[c] = new TF1(TString::Format("fnR_cat%d",c), "pol2", 500,upperxmax);
    // if (c==0){
    //   fnR[c]->SetParameter(0,3.4);
    // } else if (c==1){

    // } else if (c==2){
    // }
    gnR->Fit(TString::Format("fnR_cat%d",c), "R");
    gnR->Draw("APE");
    fnR[c]->Draw("same");   
    cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/nRVsMass_%s_cat%d.png",couplingIn.c_str(), c));

    gnL->GetYaxis()->SetTitle("n_{L} [GeV]");
    gnL->GetXaxis()->SetTitle("m_{X} [GeV]");
    fnL[c] = new TF1(TString::Format("fnL_cat%d",c), "pol2", 500,upperxmax);
    gnL->Fit(TString::Format("fnL_cat%d",c), "R");
    gnL->Draw("APE");
    fnL[c]->Draw("same");
    cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/nLVsMass_%s_cat%d.png",couplingIn.c_str(), c));
   
    //build parametric model
    //TF1: p0+p1*x+p2*x*x
    
    RooRealVar* MH = new RooRealVar("MH", "MH", 300, upperxmax);
    MH->setConstant();
    RooRealVar* p0m = new RooRealVar(TString::Format("p0m_cat%d",c), TString::Format("p0m_cat%d",c),fm[c]->GetParameter(0));
    RooRealVar* p1m = new RooRealVar(TString::Format("p1m_cat%d",c), TString::Format("p1m_cat%d",c),fm[c]->GetParameter(1));
    RooRealVar* p2m = new RooRealVar(TString::Format("p2m_cat%d",c), TString::Format("p2m_cat%d",c),fm[c]->GetParameter(2));
   		    		                     	                    
    RooRealVar* p0s = new RooRealVar(TString::Format("p0s_cat%d",c), TString::Format("p0s_cat%d",c),fs[c]->GetParameter(0));
    RooRealVar* p1s = new RooRealVar(TString::Format("p1s_cat%d",c), TString::Format("p1s_cat%d",c),fs[c]->GetParameter(1));
    RooRealVar* p2s = new RooRealVar(TString::Format("p2s_cat%d",c), TString::Format("p2s_cat%d",c),fs[c]->GetParameter(2));
     		    
    RooRealVar* p0aL = new RooRealVar(TString::Format("p0aL_cat%d",c), TString::Format("p0aL_cat%d",c),faL[c]->GetParameter(0));
    RooRealVar* p1aL = new RooRealVar(TString::Format("p1aL_cat%d",c), TString::Format("p1aL_cat%d",c),faL[c]->GetParameter(1));
    RooRealVar* p2aL = new RooRealVar(TString::Format("p2aL_cat%d",c), TString::Format("p2aL_cat%d",c),faL[c]->GetParameter(2));
     		    					                      
    RooRealVar* p0aR = new RooRealVar(TString::Format("p0aR_cat%d",c), TString::Format("p0aR_cat%d",c),faR[c]->GetParameter(0));
    RooRealVar* p1aR = new RooRealVar(TString::Format("p1aR_cat%d",c), TString::Format("p1aR_cat%d",c),faR[c]->GetParameter(1));
    RooRealVar* p2aR = new RooRealVar(TString::Format("p2aR_cat%d",c), TString::Format("p2aR_cat%d",c),faR[c]->GetParameter(2));  
    		    
    RooRealVar* p0nL = new RooRealVar(TString::Format("p0nL_cat%d",c), TString::Format("p0nL_cat%d",c),fnL[c]->GetParameter(0));
    RooRealVar* p1nL = new RooRealVar(TString::Format("p1nL_cat%d",c), TString::Format("p1nL_cat%d",c),fnL[c]->GetParameter(1));
    RooRealVar* p2nL = new RooRealVar(TString::Format("p2nL_cat%d",c), TString::Format("p2nL_cat%d",c),fnL[c]->GetParameter(2));
   		    					                      
    RooRealVar* p0nR = new RooRealVar(TString::Format("p0nR_cat%d",c), TString::Format("p0nR_cat%d",c),fnR[c]->GetParameter(0));
    RooRealVar* p1nR = new RooRealVar(TString::Format("p1nR_cat%d",c), TString::Format("p1nR_cat%d",c),fnR[c]->GetParameter(1));
    RooRealVar* p2nR = new RooRealVar(TString::Format("p2nR_cat%d",c), TString::Format("p2nR_cat%d",c),fnR[c]->GetParameter(2));
    p0m ->setConstant();
    p1m ->setConstant();
    p2m ->setConstant();
    
    p0s ->setConstant();
    p1s ->setConstant();
    p2s ->setConstant();
    
    p0aL->setConstant();
    p1aL->setConstant();
    p2aL->setConstant();
    
    p0aR->setConstant();
    p1aR->setConstant();
    p2aR->setConstant();
    
    p0nL->setConstant();
    p1nL->setConstant();
    p2nL->setConstant();
    
    p0nR->setConstant();
    p1nR->setConstant();
    p2nR->setConstant();
   
  
    std::cout<<"m: "<<fm[0]->GetParameter(0)<<" "<<fm[0]->GetParameter(1)<<" "<<fm[0]->GetParameter(2)<<" "<<fm[0]->GetParameter(0)+fm[0]->GetParameter(1)*1750+fm[0]->GetParameter(2)*1750*1750+1750<<std::endl;
    std::cout<<"s: "<<fs[0]->GetParameter(0)<<" "<<fs[0]->GetParameter(1)<<" "<<fs[0]->GetParameter(2)<<" "<<fs[0]->GetParameter(0)+fs[0]->GetParameter(1)*1750+fs[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"aL: "<<faL[0]->GetParameter(0)<<" "<<faL[0]->GetParameter(1)<<" "<<faL[0]->GetParameter(2)<<" "<<faL[0]->GetParameter(0)+faL[0]->GetParameter(1)*1750+faL[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"aR: "<<faR[0]->GetParameter(0)<<" "<<faR[0]->GetParameter(1)<<" "<<faR[0]->GetParameter(2)<<faR[0]->GetParameter(0)+faR[0]->GetParameter(1)*1750+faR[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"nL: "<<fnL[0]->GetParameter(0)<<" "<<fnL[0]->GetParameter(1)<<" "<<fnL[0]->GetParameter(2)<<fnL[0]->GetParameter(0)+fnL[0]->GetParameter(1)*1750+fnL[0]->GetParameter(2)*1750*1750<<std::endl;
    std::cout<<"nR: "<<fnR[0]->GetParameter(0)<<" "<<fnR[0]->GetParameter(1)<<" "<<fnR[0]->GetParameter(2)<<fnR[0]->GetParameter(0)+fnR[0]->GetParameter(1)*1750+fnR[0]->GetParameter(2)*1750*1750<<std::endl;
  
    RooRealVar* mgg = new RooRealVar("mgg", "mgg", 300, 6000);
    RooRealVar* thetaSmearEBEB = new RooRealVar("thetaSmearEBEB", "thetaSmearEBEB", 0., -1, 1);
    RooRealVar* thetaSmearEBEE = new RooRealVar("thetaSmearEBEE", "thetaSmearEBEE", 0., -1, 1);
    RooRealVar* thetaSmearAll = new RooRealVar("thetaSmearAll", "thetaSmearAll", 0., -1, 1);
    // RooRealVar* thetaScaleEBEB = new RooRealVar("thetaSmearEBEB", "thetaSmearEBEB", 0., -1, 1);
    // RooRealVar* thetaScaleEBEE = new RooRealVar("thetaSmearEBEE", "thetaSmearEBEE", 0., -1, 1);
    // RooRealVar* thetaScaleAll = new RooRealVar("thetaSmearAll", "thetaSmearAll", 0., -1, 1);
    RooRealVar* deltaSmear = new RooRealVar("deltaSmear", "deltaSmear",  0.001);
    // RooRealVar* deltaScale = new RooRealVar("deltaSmear", "deltaSmear",  0.001);
    deltaSmear->setConstant();
    RooRealVar* s0_EBEB = new RooRealVar("smear0EBEB", "smear0EBEB",  0.01);
    RooRealVar* s0_EBEE = new RooRealVar("smear0EBEE", "smear0EBEE",  0.015);
    RooRealVar* s0_All = new RooRealVar("smear0All", "smear0All",  0.01);
    s0_EBEB->setConstant();
    s0_EBEE->setConstant();

    RooFormulaVar* DeltaSmearEBEB = new RooFormulaVar("DeltaSmearEBEB", "2*@0*@1*@2*@2",  RooArgList(*s0_EBEB,*deltaSmear,*MH));
    RooFormulaVar* DeltaSmearEBEE = new RooFormulaVar("DeltaSmearEBEE", "2*@0*@1*@2*@2",  RooArgList(*s0_EBEE,*deltaSmear,*MH));
    RooFormulaVar* DeltaSmearAll = new RooFormulaVar("DeltaSmearAll", "2*@0*@1*@2*@2",  RooArgList(*s0_All,*deltaSmear,*MH));
    RooPlot* plot = (RooPlot*)mgg->frame(Range(300, upperxmax));
    RooDCBShape* fin_shape[theInput.size()];
    RooFormulaVar* mean; 
    RooFormulaVar* sigma0;
    RooFormulaVar* sigma; 
    RooFormulaVar* aL; 
    RooFormulaVar* aR; 
    RooFormulaVar* nL; 
    RooFormulaVar* nR; 

    if(c==0){
      mean= new RooFormulaVar("mean_EBEB", "(@0+@1*@3+@2*@3*@3)+@3", RooArgList(*p0m, *p1m, *p2m, *MH));
      sigma0= new RooFormulaVar("sigma0_EBEB", "@0+@1*@3+@2*@3*@3", RooArgList(*p0s, *p1s, *p2s, *MH));
      sigma= new RooFormulaVar("sigma_EBEB", "sqrt(@0*@0+@1*@2)", RooArgList(*sigma0,*thetaSmearEBEB,*DeltaSmearEBEB));
      aL= new RooFormulaVar("aL_EBEB","@0+@1*@3+@2*@3*@3", RooArgList(*p0aL, *p1aL, *p2aL, *MH));
      aR= new RooFormulaVar("aR_EBEB", "@0+@1*@3+@2*@3*@3", RooArgList(*p0aR, *p1aR, *p2aR, *MH));
      nL= new RooFormulaVar("nL_EBEB","@0+@1*@3+@2*@3*@3", RooArgList(*p0nL, *p1nL, *p2nL, *MH));
      nR= new RooFormulaVar("nR_EBEB","@0+@1*@3+@2*@3*@3", RooArgList(*p0nR, *p1nR, *p2nR, *MH));
    }else if(c==1){
      mean= new RooFormulaVar("mean_EBEE", "(@0+@1*@3+@2*@3*@3)+@3", RooArgList(*p0m, *p1m, *p2m, *MH));
      sigma0= new RooFormulaVar("sigma0_EBEE", "@0+@1*@3+@2*@3*@3", RooArgList(*p0s, *p1s, *p2s, *MH));
      sigma= new RooFormulaVar("sigma_EBEE", "sqrt(@0*@0+@1*@2)", RooArgList(*sigma0,*thetaSmearEBEE,*DeltaSmearEBEE));
      aL= new RooFormulaVar("aL_EBEE","@0+@1*@3+@2*@3*@3", RooArgList(*p0aL, *p1aL, *p2aL, *MH));
      aR= new RooFormulaVar("aR_EBEE", "@0+@1*@3+@2*@3*@3", RooArgList(*p0aR, *p1aR, *p2aR, *MH));
      nL= new RooFormulaVar("nL_EBEE","@0+@1*@3+@2*@3*@3", RooArgList(*p0nL, *p1nL, *p2nL, *MH));
      nR= new RooFormulaVar("nR_EBEE","@0+@1*@3+@2*@3*@3", RooArgList(*p0nR, *p1nR, *p2nR, *MH));

    }else if(c==2){
      mean= new RooFormulaVar("mean_All", "(@0+@1*@3+@2*@3*@3)+@3", RooArgList(*p0m, *p1m, *p2m, *MH));
      sigma0= new RooFormulaVar("sigma0_All", "@0+@1*@3+@2*@3*@3", RooArgList(*p0s, *p1s, *p2s, *MH));
      sigma= new RooFormulaVar("sigma_All", "sqrt(@0*@0+@1*@2)", RooArgList(*sigma0,*thetaSmearAll,*DeltaSmearAll));
      aL= new RooFormulaVar("aL_All","@0+@1*@3+@2*@3*@3", RooArgList(*p0aL, *p1aL, *p2aL, *MH));
      aR= new RooFormulaVar("aR_All", "@0+@1*@3+@2*@3*@3", RooArgList(*p0aR, *p1aR, *p2aR, *MH));
      nL= new RooFormulaVar("nL_All","@0+@1*@3+@2*@3*@3", RooArgList(*p0nL, *p1nL, *p2nL, *MH));
      nR= new RooFormulaVar("nR_All","@0+@1*@3+@2*@3*@3", RooArgList(*p0nR, *p1nR, *p2nR, *MH));

    }

    for(unsigned int M =0; M<theInput.size(); M++){
      
      MH->setVal(masses[M]);
      // MH->setConstant();
      std::cout<<MH->getVal()<<std::endl;
      
      std::cout<<mean->getVal(*MH)<<std::endl;
      std::cout<<sigma->getVal(*MH)<<std::endl;
      std::cout<<aL->getVal(*MH)<<std::endl;
      std::cout<<aR->getVal(*MH)<<std::endl;
      std::cout<<nL->getVal(*MH)<<std::endl;
      std::cout<<nR->getVal(*MH)<<std::endl;
   
      fin_shape[M] =  new RooDCBShape(Form("SignalShape_%s_cat%d_M%f",couplingIn.c_str(), c, masses[M]),TString::Format("final_shape_%s_cat%d_M%f",couplingIn.c_str(),c, masses[M] ) ,*mgg, *mean, *sigma,  *aL, *aR,  *nL, *nR) ;

      fin_shape[M]->plotOn(plot, LineColor(M+1));
     
    }

    plot->GetYaxis()->SetTitle("a.u.");
    plot->GetYaxis()->SetRangeUser(0.001, 5);
    plot->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    plot->Draw();
    legmc->Draw("same");
    cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/final_shapes_%s_cat%d.png",couplingIn.c_str(), c));
    cc1->SetLogy(); 
    cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/final_shapes_%s_cat%d_log.png",couplingIn.c_str(), c));

    if(c==0)    sigshape[c] =  new RooDCBShape(Form("SignalShape_%s_EBEB",couplingIn.c_str()), Form("SignalShape_%s_EBEB",couplingIn.c_str()), *mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    if(c==1)    sigshape[c] =  new RooDCBShape(Form("SignalShape_%s_EBEE",couplingIn.c_str()), Form("SignalShape_%s_EBEE",couplingIn.c_str()), *mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    if(c==2)    sigshape[c] =  new RooDCBShape(Form("SignalShape_%s_All",couplingIn.c_str()), Form("SignalShape_%s_All",couplingIn.c_str()) ,*mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    ws_out->import(*sigshape[c]);
    
  } //end of loop over categories
   
   RooRealVar* p0_cat0 = new RooRealVar("p0_cat0", "p0_cat0", ffhmw[0]->GetParameter(0) );
   RooRealVar* p1_cat0 = new RooRealVar("p1_cat0", "p1_cat0", ffhmw[0]->GetParameter(1) );
   RooRealVar* p0_cat1 = new RooRealVar("p0_cat1", "p0_cat1", ffhmw[1]->GetParameter(0) );
   RooRealVar* p1_cat1 = new RooRealVar("p1_cat1", "p1_cat1", ffhmw[1]->GetParameter(1) );

   RooFormulaVar* FHWM_EBEB= new RooFormulaVar("FHWM_EBEB", "(@0+@1*@2)", RooArgList(*p0_cat0,*p1_cat0, *MH));
   RooFormulaVar* FHWM_EBEE= new RooFormulaVar("FHWM_EBEE", "(@0+@1*@2)", RooArgList(*p0_cat1,*p1_cat1, *MH));
   ws_out->import(*FHWM_EBEB);
   ws_out->import(*FHWM_EBEE);
  
   TFile* fout= new TFile(Form("%s/SignalParametricShapes_ws_%s.root", ws_outdir.c_str(), couplingIn.c_str() ) , "RECREATE");
   fout->cd();
   fm[0]->Write();       fm[1]->Write(); 
   fs[0]->Write(); 	 fs[1]->Write(); 
   faL[0]->Write();	 faL[1]->Write();
   faR[0]->Write();	 faR[1]->Write();
   fnL[0]->Write();	 fnL[1]->Write();
   fnR[0]->Write();	 fnR[1]->Write();
   ffhmw[0]->Write();	 ffhmw[1]->Write();

   ws_out->Print();
   
   ws_out->Write();
   fout->Write();
   fout->Close();
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
