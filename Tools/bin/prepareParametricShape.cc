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
#include "TStyle.h"

using namespace RooFit;
using namespace RooStats;

int nBinsMass = 1100;//150//1250

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
void plotAllSignalsAsimov(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &couplingIn, std::vector<std::string> cats);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);
std::string widthtonum(std::string coupling);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year,  inputsamples, inputworkspaces, outputworkspaces, thecoupling;

  if(argc!=6) {
    std::cout << "Syntax: prepareParametricShape.exe [2016/2017/2018] [inputsamples] [inputworkspaces] [outputworkspaces] [kMplxxx]" << std::endl;
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
    thecoupling = argv[5];
 }

  //Categories
  std::vector<std::string> cats; 
  cats.clear(); 
  cats.push_back("EBEB");
  cats.push_back("EBEE");

  //========================================================================
  //read the reduced input root trees
  initForFit(inputsamples);
  //========================================================================
  //Run the fits
  plotAllSignalsAsimov(year,inputworkspaces, outputworkspaces, thecoupling, cats);

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
    ccc->SaveAs(Form("%s.png",name.c_str()));
  }

  return xWidth;
}


//-----------------------------------------------------------------------------------
void plotAllSignalsAsimov(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &couplingIn, std::vector<std::string> cats){

  std::vector<std::string> samples = getSampleListForFit();
  bool plot = true;
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptStat(1111);
  m_gStyle->SetOptFit(1);

  std::string coup = widthtonum(couplingIn);

  //========================================================================
  std::cout << "Calculate the final parametric shape" <<std::endl;
  std::cout << "Input dir " << ws_indir << std::endl;
  std::cout << "Output dir " << ws_outdir << std::endl;
  std::string coupling = "";
  std::string M_bins = "";

  double MassMin = 0.; 
  double MassMax = 0.; 
  if ( couplingIn == "kMpl001" || couplingIn == "0p014"){MassMin = 740.; MassMax = 5000. ;}
  else if ( couplingIn == "kMpl01" || couplingIn == "1p4"){MassMin = 740.; MassMax = 5000. ;}
  else if ( couplingIn == "kMpl02" || couplingIn == "5p6"){MassMin = 740.; MassMax = 5000. ;}
  else {
    std::cout << "Only 'kMpl001', 'kMpl01', 'kMpl02', '0p014', '1p4' and '5p6' are allowed. " << std::endl;
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

    if (isample.find("RSGraviton") != std::string::npos) {
      coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    } else if ( isample.find("GluGluSpin0") != std::string::npos ){
      coupling = get_str_between_two_str(getBase(isample), "GluGluSpin0ToGammaGamma_W_", "_M_");
    }

    M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");
    std::string samplename = isample.substr(0,getBase(isample).find("ToGammaGamma")); 

    tmpin.kMpl = coupling;
    tmpin.M_bins = M_bins;
    tmpin.name = getBase(isample);
    tmpin.fresandgen_final = TFile::Open(Form("%s/ws_%s_ResponseAndGen_M%s_k%s_%s_final.root", ws_indir.c_str(), samplename.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()));  
    tmpin.ws = (RooWorkspace*) tmpin.fresandgen_final->Get("HLFactory_ws");
 
    theInput.push_back(tmpin);

  }

  unsigned int inputsize = 1;
  for(unsigned int iM =0; iM < theInput.size(); iM++){
    theTree tmpin = theInput[iM];
    if ( std::stod(tmpin.M_bins) > 4001. ){continue;}//HARDCODED: Just keep it in mind. 
    ++inputsize;

  }

  
  double m[inputsize];
  double mErr[inputsize];
  double s[inputsize];
  double sErr[inputsize];
  double aL[inputsize];
  double aLErr[inputsize];
  double aR[inputsize];
  double aRErr[inputsize];
  double nL[inputsize];
  double nLErr[inputsize];
  double nR[inputsize];
  double nRErr[inputsize];
  double masses[inputsize];
  double massesErr[inputsize];

  TString svar = "responseaddpdf";
  TString spdf = "responseaddpdf";
  
  RooRealVar* MH = new RooRealVar("MH", "MH", 300, 9000);
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

    FILE *resFile;
    resFile = fopen(Form("%s/../plots/paramsvsMX_fits_%s_%s_cat%d.txt", ws_outdir.c_str(),year.c_str(),couplingIn.c_str(),c),"w");
    if (couplingIn=="kMpl001"){
      fprintf(resFile,"\\hline\n");
      fprintf(resFile,"\\hline\n");
      fprintf(resFile,"\\multicolumn{6}{c}{Year: %s } \\\\\n", year.c_str());
      fprintf(resFile,"\\hline\n");
      fprintf(resFile,"Parameter (GeV) & Width & $p_{0}$ & $p_{1}$ & $p_{2}$ & $\\chi^{2}/N  $ \\\\\n");
      fprintf(resFile,"\\hline\n");
      fprintf(resFile,"\\multicolumn{6}{c}{Category: %s } \\\\\n", cats[c].c_str());
      fprintf(resFile,"\\hline\n");
    }

    std::map< std::string , RooDCBShape* > res;
    std::map< std::string , TH1F* > h;
    std::map< std::string , RooDataHist* > resdata; 
    std::map< std::string , RooPlot* > pres; 
    std::map< std::string , RooArgSet* > model_params;
    double fhwm[inputsize];

    TLegend* legmc = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);

    TCanvas* cc1 = new TCanvas(Form("cc1_cat%d", c), Form("cc1_cat%d", c));
    unsigned int iMindex = 0;
    
    for(unsigned int iM =0; iM < theInput.size(); iM++){

      theTree tmpin = theInput[iM];

      if ( std::stod(tmpin.M_bins) > 4001. ){continue;}//HARDCODED: Just keep it in mind. 

      //Should be after the continue to count only the points we want for the
      //chi^2/ndof table later. 
      ++iMindex;
      masses[iMindex] = std::stod(tmpin.M_bins); 
      massesErr[iMindex] = 0.;

      res[tmpin.name] = (RooDCBShape*) tmpin.ws->pdf(TString::Format("dcbshape_cat%d",c));
      //compute FWHM
      fhwm[iMindex] = computePdfFHWM(*res[tmpin.name],*tmpin.ws->var("var"), std::stoi(tmpin.M_bins), plot, Form("%s/../plots/test_%d_%d", ws_outdir.c_str(),iMindex,c));
      std::cout << "FWHM " << fhwm[iMindex] <<" for mass " << tmpin.M_bins << std::endl;

      cc1->cd();
      h[tmpin.name] = new TH1F(Form("h%d_%d",iMindex,c), Form("h%d_%d",iMindex,c),460, 300,MassMax);

      h[tmpin.name]->SetLineColor(iMindex+1);
      h[tmpin.name]->SetMarkerColor(iMindex+1);

      legmc->AddEntry( h[tmpin.name] , Form("M_{X} = %d GeV", std::stoi(tmpin.M_bins) ),"pl" );

      if(c<2){  resdata[tmpin.name] = (RooDataHist*)tmpin.ws->data(TString::Format("signal_asimov_cat%d",c));}
      if(c==2){ resdata[tmpin.name] = (RooDataHist*)tmpin.ws->data("signal_asimov"); }

      pres[tmpin.name] = tmpin.ws->var("var")->frame(Range(300,MassMax+1000.),Title("mass asimov"),Bins(420));
      resdata[tmpin.name]->plotOn(pres[tmpin.name],MarkerColor(iMindex+1),LineColor(iMindex+1));
      res[tmpin.name]->plotOn(pres[tmpin.name],LineColor(iMindex+1));

      pres[tmpin.name]->GetXaxis()->SetTitle("#Delta m [GeV]");
      pres[tmpin.name]->GetYaxis()->SetTitle("a.u.");
      // pres[tmpin.name]->GetYaxis()->SetRangeUser(0.001,MassMax);

      if (iMindex == 0){
	pres[tmpin.name]->Draw();
      } else{
	pres[tmpin.name]->Draw("same");
      }

      legmc->Draw("same");
   
      //save parameters 
      model_params[tmpin.name] = res[tmpin.name]->getParameters(*tmpin.ws->var("var")) ;
      model_params[tmpin.name]->Print("v") ;

      m[iMindex] = ((RooRealVar*) tmpin.ws->var(TString::Format("dcb_m_cat%d",c)))->getVal()- ((RooRealVar*)tmpin.ws->var("MH"))->getVal();
      mErr[iMindex] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_m_cat%d",c)))->getError() ;
      std::cout<<m[iMindex]<<" "<<mErr[iMindex]<<std::endl;
      
      s[iMindex] = ((RooRealVar*) tmpin.ws->var(TString::Format("dcb_s_cat%d",c)))->getVal();
      sErr[iMindex] =((RooRealVar*) tmpin.ws->var(TString::Format("dcb_s_cat%d",c)))->getError() ;
      std::cout<<s[iMindex]<<" "<<sErr[iMindex]<<std::endl;

      aL[iMindex] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aL_cat%d",c)))->getVal();
      aLErr[iMindex] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aL_cat%d",c)))->getError() ;
      std::cout<<aL[iMindex]<<" "<<aLErr[iMindex]<<std::endl;
      
      aR[iMindex] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aR_cat%d",c)))->getVal();
      aRErr[iMindex] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_aR_cat%d",c)))->getError() ;
      std::cout<<aR[iMindex]<<" "<<aRErr[iMindex]<<std::endl;
      
      nR[iMindex] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nR_cat%d",c)))->getVal();
      nRErr[iMindex] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nR_cat%d",c)))->getError() ;
      std::cout<<nR[iMindex]<<" "<<nRErr[iMindex]<<std::endl;

      nL[iMindex] = ((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nL_cat%d",c)))->getVal();
      nLErr[iMindex] =((RooRealVar*)tmpin.ws->var(TString::Format("dcb_nL_cat%d",c)))->getError() ;
      std::cout<<nL[iMindex]<<" "<<nLErr[iMindex]<<std::endl;

   } //end of loop over signal samples of the same coupling
   
   
   cc1->SaveAs(Form("%s/../plots/asimovAllMasses_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));
   cc1->SaveAs(Form("%s/../plots/asimovAllMasses_%s_cat%d.root", ws_outdir.c_str(),couplingIn.c_str(), c));
   // cc1->SetLogy();
   // cc1->SaveAs(Form("%s/../plots/asimovAllMasses_%s_cat%d_log.png", ws_outdir.c_str(),couplingIn.c_str(), c));

   // cc1->SetLogy(0);
 
   TGraph* gr = new TGraph(inputsize,masses,fhwm);
   TGraphErrors* gm = new TGraphErrors(inputsize, masses, m, massesErr, mErr);
   TGraphErrors* gs = new TGraphErrors(inputsize, masses, s, massesErr, sErr);
   TGraphErrors* gaL = new TGraphErrors(inputsize, masses, aL, massesErr, aLErr);
   TGraphErrors* gaR = new TGraphErrors(inputsize, masses, aR, massesErr, aRErr);
   TGraphErrors* gnR = new TGraphErrors(inputsize, masses, nR, massesErr, nRErr);
   TGraphErrors* gnL = new TGraphErrors(inputsize, masses, nL, massesErr, nLErr);

   ffhmw[c] = new TF1(TString::Format("ffhmw_cat%d",c), "pol1", MassMin,MassMax);
   gr->Fit(TString::Format("ffhmw_cat%d",c), "R");
   gr->GetXaxis()->SetTitle("m_X [GeV]");
   gr->GetYaxis()->SetTitle("FHWM [GeV]");
   
   gr->SetMarkerStyle(20);
   gr->SetTitle(" ");
   gr->GetXaxis()->SetRangeUser(250., MassMax + 100.);
   gr->GetYaxis()->SetTitleSize(0.045);
   gr->GetXaxis()->SetTitleSize(0.045);
   gr->GetYaxis()->SetTitleOffset(0.90);
   gr->GetXaxis()->SetTitleOffset(0.90);

   gr->Draw("AP");
   ffhmw[c]->Draw("same");
   // cc1->SaveAs(Form("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/FinalParametricShape/FHWM_%s_cat%d.png",couplingIn.c_str(), c));
   cc1->SaveAs(Form("%s/../plots/FHWM_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));
   // cc1->SaveAs(Form("%s/../plots/FHWM_%s_cat%d.root", ws_outdir.c_str(),couplingIn.c_str(), c));

   fm[c] = new TF1(TString::Format("fm_cat%d",c), "pol2", MassMin,MassMax);
   gm->Fit(TString::Format("fm_cat%d",c), "R");
   gm->GetYaxis()->SetTitle("#Delta m = m - m_{H} [GeV]");
   gm->GetXaxis()->SetTitle("m_{X}[GeV]");
   gm->SetMarkerStyle(20);
   gm->SetTitle(" ");
   gm->GetXaxis()->SetRangeUser(250., MassMax + 100.);
   gm->GetYaxis()->SetTitleSize(0.045);
   gm->GetXaxis()->SetTitleSize(0.045);
   gm->GetYaxis()->SetTitleOffset(0.90);
   gm->GetXaxis()->SetTitleOffset(0.90);
   gm->Draw("APE");
   fm[c]->Draw("same");
   cc1->SaveAs(Form("%s/../plots/meanVsMass_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));

   
   gs->GetYaxis()->SetTitle("#sigma [GeV]");
   gs->GetXaxis()->SetTitle("m_{X}[GeV]");
   fs[c] = new TF1(TString::Format("fs_cat%d",c), "pol2", MassMin,MassMax);
   gs->Fit(TString::Format("fs_cat%d",c), "R");
   gs->SetMarkerStyle(20);
   gs->SetTitle(" ");
   gs->GetXaxis()->SetRangeUser(250., MassMax + 100.);
   gs->GetYaxis()->SetTitleSize(0.045);
   gs->GetXaxis()->SetTitleSize(0.045);
   gs->GetYaxis()->SetTitleOffset(0.90);
   gs->GetXaxis()->SetTitleOffset(0.90);
   gs->Draw("APE");
   fs[c]->Draw("same");
   cc1->SaveAs(Form("%s/../plots/sigmaVsMass_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));


   gaL->GetYaxis()->SetTitle("#alpha_{L} [GeV]");
   gaL->GetXaxis()->SetTitle("m_{X}[GeV]");
   faL[c] = new TF1(TString::Format("faL_cat%d",c), "pol2", MassMin,MassMax);
   gaL->Fit(TString::Format("faL_cat%d",c), "R");
   gaL->SetMarkerStyle(20);
   gaL->SetTitle(" ");
   gaL->GetXaxis()->SetRangeUser(250., MassMax + 100.);
   gaL->GetYaxis()->SetTitleSize(0.045);
   gaL->GetXaxis()->SetTitleSize(0.045);
   gaL->GetYaxis()->SetTitleOffset(0.90);
   gaL->GetXaxis()->SetTitleOffset(0.90);
   gaL->Draw("APE");
   faL[c]->Draw("same");
   cc1->SaveAs(Form("%s/../plots/aLVsMass_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));
  
   gaR->GetYaxis()->SetTitle("#alpha_{R} [GeV]");
   gaR->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
   faR[c] = new TF1(TString::Format("faR_cat%d",c), "pol2", MassMin,MassMax);
   gaR->Fit(TString::Format("faR_cat%d",c), "R");
   gaR->SetMarkerStyle(20);
   gaR->SetTitle(" ");
   gaR->GetXaxis()->SetRangeUser(250., MassMax + 100.);
   gaR->GetYaxis()->SetTitleSize(0.045);
   gaR->GetXaxis()->SetTitleSize(0.045);
   gaR->GetYaxis()->SetTitleOffset(0.90);
   gaR->GetXaxis()->SetTitleOffset(0.90);
   gaR->Draw("APE");
   faR[c]->Draw("same");
   cc1->SaveAs(Form("%s/../plots/aRVsMass_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));

 
   gnR->GetYaxis()->SetTitle("n_{R} [GeV]");
   gnR->GetXaxis()->SetTitle("m_{X}[GeV]");
   fnR[c] = new TF1(TString::Format("fnR_cat%d",c), "pol2", MassMin,MassMax);
   // if (c==0){
   //   fnR[c]->SetParameter(0,3.4);
   // } else if (c==1){

   // } else if (c==2){
   // }
   gnR->Fit(TString::Format("fnR_cat%d",c), "R");
   gnR->SetMarkerStyle(20);
   gnR->SetTitle(" ");
   gnR->GetXaxis()->SetRangeUser(250., MassMax + 100.);
   gnR->GetYaxis()->SetTitleSize(0.045);
   gnR->GetXaxis()->SetTitleSize(0.045);
   gnR->GetYaxis()->SetTitleOffset(0.90);
   gnR->GetXaxis()->SetTitleOffset(0.90);
   gnR->Draw("APE");
   fnR[c]->Draw("same");   
   cc1->SaveAs(Form("%s/../plots/nRVsMass_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));

   gnL->GetYaxis()->SetTitle("n_{L} [GeV]");
   gnL->GetXaxis()->SetTitle("m_{X} [GeV]");
   fnL[c] = new TF1(TString::Format("fnL_cat%d",c), "pol2", MassMin,MassMax);
   gnL->Fit(TString::Format("fnL_cat%d",c), "R");
   gnL->SetMarkerStyle(20);
   gnL->SetTitle(" ");
   gnL->GetXaxis()->SetRangeUser(250., MassMax + 100.);
   gnL->GetYaxis()->SetTitleSize(0.045);
   gnL->GetXaxis()->SetTitleSize(0.045);
   gnL->GetYaxis()->SetTitleOffset(0.90);
   gnL->GetXaxis()->SetTitleOffset(0.90);
   gnL->Draw("APE");
   fnL[c]->Draw("same");
   cc1->SaveAs(Form("%s/../plots/nLVsMass_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));

   //Let's make the table
   fprintf(resFile," $\\Delta m = m - m_{H} $ & $ %s $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f $ \\\\\n", coup.c_str(), fm[c]->GetParameter(0), fm[c]->GetParError(0), fm[c]->GetParameter(1), fm[c]->GetParError(1), fm[c]->GetParameter(2), fm[c]->GetParError(2), fm[c]->GetChisquare()/fm[c]->GetNDF() );
   fprintf(resFile," $\\sigma $ & $ %s $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f $ \\\\\n", coup.c_str(), fs[c]->GetParameter(0), fs[c]->GetParError(0), fs[c]->GetParameter(1), fs[c]->GetParError(1), fs[c]->GetParameter(2), fs[c]->GetParError(2), fs[c]->GetChisquare()/fs[c]->GetNDF() );
   fprintf(resFile," $\\alpha_{L} $ & $ %s $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f $ \\\\\n", coup.c_str(), faL[c]->GetParameter(0), faL[c]->GetParError(0), faL[c]->GetParameter(1), faL[c]->GetParError(1), faL[c]->GetParameter(2), faL[c]->GetParError(2), faL[c]->GetChisquare()/faL[c]->GetNDF() );
   fprintf(resFile," $\\alpha_{R} $ & $ %s $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f $ \\\\\n", coup.c_str(), faR[c]->GetParameter(0), faR[c]->GetParError(0), faR[c]->GetParameter(1), faR[c]->GetParError(1), faR[c]->GetParameter(2), faR[c]->GetParError(2), faR[c]->GetChisquare()/faR[c]->GetNDF() );
   fprintf(resFile," $ n_{L} $ & $ %s $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f $ \\\\\n", coup.c_str(), fnL[c]->GetParameter(0), fnL[c]->GetParError(0), fnL[c]->GetParameter(1), fnL[c]->GetParError(1), fnL[c]->GetParameter(2), fnL[c]->GetParError(2), fnL[c]->GetChisquare()/fnL[c]->GetNDF() );
   fprintf(resFile," $ n_{R} $ & $ %s $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f \\pm %10.3f $ & $ %10.3f $ \\\\\n", coup.c_str(), fnR[c]->GetParameter(0), fnR[c]->GetParError(0), fnR[c]->GetParameter(1), fnR[c]->GetParError(1), fnR[c]->GetParameter(2), fnR[c]->GetParError(2), fnR[c]->GetChisquare()/fnR[c]->GetNDF() );
    fprintf(resFile,"\\hline\n");

   fclose(resFile);

   //build parametric model
   //TF1: p0+p1*x+p2*x*x
    
   RooRealVar* MH = new RooRealVar("MH", "MH", 300, MassMax);
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
  
    RooRealVar* mgg = new RooRealVar("mgg", "mgg", nBinsMass, 500., 6000.);
    mgg->setBins(nBinsMass);

    // RooBinning* mggRooBinning = new RooBinning(nBinsMass, 500., 6000.,"mgg");
    // mgg->setBinning(*mggRooBinning);
    // std::cout<< "++++++++++++++++++++++++++++++++"<< std::endl; 
    // std::cout<< mgg->getBinning().numBins()<< std::endl;
    std::cout<< mgg->getBins()<< std::endl;
    // std::cout<< "++++++++++++++++++++++++++++++++"<< std::endl; 
    
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
    RooPlot* plot = (RooPlot*)mgg->frame(Range(300, MassMax+1000.));
    RooDCBShape* fin_shape[inputsize];
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

    for(unsigned int M =0; M<inputsize; M++){
      
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
    cc1->SaveAs(Form("%s/../plots/final_shapes_%s_cat%d.png", ws_outdir.c_str(),couplingIn.c_str(), c));
    cc1->SetLogy(); 
    cc1->SaveAs(Form("%s/../plots/final_shapes_%s_cat%d_log.png", ws_outdir.c_str(),couplingIn.c_str(), c));
    cc1->SetLogy(0); 


    if(c==0)    sigshape[c] =  new RooDCBShape(Form("SignalShape_%s_EBEB",couplingIn.c_str()), Form("SignalShape_%s_EBEB",couplingIn.c_str()), *mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    if(c==1)    sigshape[c] =  new RooDCBShape(Form("SignalShape_%s_EBEE",couplingIn.c_str()), Form("SignalShape_%s_EBEE",couplingIn.c_str()), *mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    if(c==2)    sigshape[c] =  new RooDCBShape(Form("SignalShape_%s_All",couplingIn.c_str()), Form("SignalShape_%s_All",couplingIn.c_str()) ,*mgg, *mean, *sigma,  *aL, *aR,  *nL,*nR) ;
    ws_out->import(*sigshape[c]);
    
  } //end of loop over categories
   
   RooRealVar* p0_cat0 = new RooRealVar("p0_cat0", "p0_cat0", ffhmw[0]->GetParameter(0) );
   RooRealVar* p1_cat0 = new RooRealVar("p1_cat0", "p1_cat0", ffhmw[0]->GetParameter(1) );
   RooRealVar* p0_cat1 = new RooRealVar("p0_cat1", "p0_cat1", ffhmw[1]->GetParameter(0) );
   RooRealVar* p1_cat1 = new RooRealVar("p1_cat1", "p1_cat1", ffhmw[1]->GetParameter(1) );
   RooRealVar* p0_cat2 = new RooRealVar("p0_cat2", "p0_cat2", ffhmw[2]->GetParameter(0) );
   RooRealVar* p1_cat2 = new RooRealVar("p1_cat2", "p1_cat2", ffhmw[2]->GetParameter(1) );

   RooFormulaVar* FHWM_EBEB= new RooFormulaVar("FHWM_EBEB", "(@0+@1*@2)", RooArgList(*p0_cat0,*p1_cat0, *MH));
   RooFormulaVar* FHWM_EBEE= new RooFormulaVar("FHWM_EBEE", "(@0+@1*@2)", RooArgList(*p0_cat1,*p1_cat1, *MH));
   RooFormulaVar* FHWM_All= new RooFormulaVar("FHWM_All", "(@0+@1*@2)", RooArgList(*p0_cat2,*p1_cat2, *MH));
   ws_out->import(*FHWM_EBEB);
   ws_out->import(*FHWM_EBEE);
   ws_out->import(*FHWM_All);

   TFile* fout= new TFile(Form("%s/SignalParametricShapes_ws_%s.root", ws_outdir.c_str(), couplingIn.c_str() ) , "RECREATE");
   fout->cd();
   fm[0]->Write();       fm[1]->Write(); 
   fs[0]->Write(); 	 fs[1]->Write(); 
   faL[0]->Write();	 faL[1]->Write();
   faR[0]->Write();	 faR[1]->Write();
   fnL[0]->Write();	 fnL[1]->Write();
   fnR[0]->Write();	 fnR[1]->Write();
   ffhmw[0]->Write();	 ffhmw[1]->Write(); ffhmw[2]->Write();

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

//-----------------------------------------------------------------------------------
std::string widthtonum(std::string coupling){
  std::string coup = ""; 
 
  if ( coupling == "0p014" ){coup = "1.4 \\times 10^{-4}";}
  else if ( coupling == "1p4" ){coup = "1.4 \\times 10^{-2}";}
  else if ( coupling == "5p6" ){coup = "5.6 \\times 10^{-2}";}
  else if ( coupling == "kMpl001" ){coup = "0.01";}
  else if ( coupling == "kMpl01" ){coup = "0.1";}
  else if ( coupling == "kMpl02" ){coup = "0.2";}

  return coup;


}
