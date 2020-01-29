#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
#include "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

//RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooStats/HLFactory.h"
#include "RooPlot.h"
#include "RooConstVar.h"
#include "RooPolyVar.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "RooStats/RooStatsUtils.h"
#include "RooNLLVar.h"
#include "TRandom3.h"

//ROOT
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TBox.h"

using namespace RooFit;
using namespace RooStats;

static const Int_t NCAT = 2; //BB and BE for the moment 
float MINmass = 500.;
float MINmassBE = 500.;
//Will go up to 6000 for the moment. 
float MAXmass = 6000.;
//In the full range (0,6000) will have 20 GeV bins same as in the past. 
int nBinsMass = 300;
//Here 5 GeV
// int nBinsMass = 1200;

bool gofwithtoys=false;

TRandom3 *RandomGen = new TRandom3();

struct window {
  std::string name; 
  float low;
  float high;
  double norm;
  std::string truemodel;
};

struct thepdf {

  std::string cat; 
  int catnum; 
  std::string model; 
  RooAbsPdf* pdf; 
  int order; 
  std::vector<RooDataSet*> toys;
  RooWorkspace* ws; 
  double nominalnorm;
  std::vector<window> truenorms;

};

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
struct biasfitres {
  std::string cat; 
  std::string model; 
  int toynum;
  std::string windname; 
  double windnorm;
  double nomnorm;
  int minos;
  double hesseerr; 
  double fiterrh; 
  double fiterrl; 
  double bias; 
  double minMassFit; 
  double maxMassFit;


};

//-----------------------------------------------------------------------------------
//Declarations here definition after main
std::map<std::string, std::vector<window> > computewindnorm(const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, std::vector<std::string> models);
void writeToJson( std::map<std::string, std::vector<window> > valuestowrite, std::string outputfile);
std::map<std::string, std::vector<window> > readJson(std::string inputfile, std::string model);
std::vector<thepdf> throwtoys(RooWorkspace* w,const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, std::vector<std::string> models, int ntoystart, int ntoyend, std::map<std::string, std::vector<window> > windows);
std::vector<biasfitres> fitToys(RooWorkspace* w, std::vector<thepdf> thetoys, std::string familymodel, bool blind, bool dobands, const std::string &year, const std::string &ws_dir, const std::string &out_dir, std::vector<std::string> cats, std::vector<std::string> models, bool approx_minos);
void makebiasrootfile(std::vector<biasfitres> bfr, const std::string &out_dir, std::string familymodel, std::string year, int ntoystart, int ntoyend);
void analyzeBias(const std::string &ws_dir, const std::string &out_dir, std::vector<std::string> models, std::vector<std::string> cats, std::string year, int ntoystart, int ntoyend);
theFitResult theFit(RooAbsPdf *pdf, RooExtendPdf *epdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, bool fittoys);
gof PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, std::string model, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir, int order);
TPaveText* get_labelsqrt( int legendquadrant );
TPaveText* get_labelcms( int legendquadrant, std::string year, bool sim);
gof getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name, bool gofToys);
std::map<TString , TString > buildBiasTerm(std::string year, std::vector<std::string> cats, std::string paper, bool plotto1invfb);
void plotBiasTerm(const std::string &ws_dir, std::string year, std::vector<std::string> cats);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir, outputdir, themodel, ntoystart, ntoyend, Json;
  

  if(argc!=8) {
    std::cout << "Syntax: bias_study.exe [2016/2017/2018] [input] [output] [model] [toystart] [toysend] [Json]" << std::endl;
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
    themodel = argv[4];
    ntoystart = argv[5];
    ntoyend = argv[6];
    Json = argv[7];
    if (Json !="createJson" and Json!="readJson"){
      std::cout << "Only createJson and readJson are allowed " << std::endl;
    }

 }

  //========================================================================
  //========================================================================
  //ZERO PART : We give whatever input/output variables needed here.  
  //========================================================================
  //========================================================================
  //========================================================================

  //Categories
  std::vector<std::string> cats; 
  cats.clear(); 
  cats.push_back("EBEB");
  cats.push_back("EBEE");

  //Models
  std::vector<std::string> models;
  // models.push_back("pow");
  // models.push_back("expow");
  // models.push_back("invpow");
  // models.push_back("invpowlin");
  // models.push_back("moddijet");
  //Due to speed issues will run each model separately and sent to condor. 
  models.push_back(themodel);

  bool blind = true;
  bool dobands = false;
  bool approxminos = true;

  //========================================================================
  //========================================================================
  //FIRST PART: Create or read json file containing the number of events 
  //            based on the true underlying distribution h(mγγ) for all mass 
  //            windows. I was forced to save that to a file due to a 
  //            segmentation fault coming from the createIntegral on the mgg, 
  //            which whould be loaded from the intput file workspace and not 
  //            from the new created workspace contained the dijet function. 
  //========================================================================
  //========================================================================

  //========================================================================
  std::map<std::string, std::vector<window> > windows;
  //========================================================================
  //Save the values in a json file 
  if (Json == "createJson"){
    windows = computewindnorm(year, inputdir, cats, models);
    writeToJson(windows, Form("Tools/json/windowsnorm_%s.json", year.c_str()));
    return 0; 
  } else if (Json == "readJson"){
    windows = readJson(Form("%s/../../Tools/json/windowsnorm_%s.json", inputdir.c_str(), year.c_str()), themodel);
  }

  //========================================================================
  //========================================================================
  //SECOND PART: Throw toys. The toys thrown from the alternative models will 
  //             be fit with the dijet ansantz.
  //========================================================================
  //========================================================================
  TString card_name = "HighMass-hgg_models_Bkg.bias.rs" ; 
  HLFactory *hlf = new HLFactory("HLFactory", card_name, false);
  RooWorkspace *w = hlf->GetWs();
  w->Print("V");
  std::vector<thepdf> thetoys = throwtoys(w, year, inputdir, cats, models, std::stoi(ntoystart), std::stoi(ntoyend), windows);

  //========================================================================
  //========================================================================
  //THIRD PART: Fit toys and save necessary info for the bias study. 
  //========================================================================
  //========================================================================
  std::vector<biasfitres> thebiasres = fitToys(w, thetoys, themodel, blind, dobands, year, inputdir, outputdir, cats, models, approxminos);
  makebiasrootfile(thebiasres, outputdir, themodel, year, std::stoi(ntoystart), std::stoi(ntoyend));

  //========================================================================
  //========================================================================
  //FOURTH PART: Analyze the bias results. 
  //========================================================================
  //========================================================================
  analyzeBias(inputdir, outputdir, models, cats, year, std::stoi(ntoystart), std::stoi(ntoyend));



}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
std::map<std::string, std::vector<window> > computewindnorm(const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, std::vector<std::string> models){

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  //This is what we will return, the true norms. 
  //It contains: 
  //name,low,high,norm,truemodel 
  std::map<std::string, std::vector<window> > truenormwinds;//[c][window]
  truenormwinds.clear();

  //We build a helpful variabel to calculate the norms. 
  std::vector<window> windows;
  
  window tmpwind;
  tmpwind.low = 500.; tmpwind.high = 550.; tmpwind.name = "500_550";
  windows.push_back(tmpwind);
  tmpwind.low = 550.; tmpwind.high = 600.; tmpwind.name = "550_600";
  windows.push_back(tmpwind);
  tmpwind.low = 600.; tmpwind.high = 650.; tmpwind.name = "600_650";
  windows.push_back(tmpwind);
  tmpwind.low = 650.; tmpwind.high = 700.; tmpwind.name = "650_700";
  windows.push_back(tmpwind);
  tmpwind.low = 700.; tmpwind.high = 750.; tmpwind.name = "700_750";
  windows.push_back(tmpwind);
  tmpwind.low = 750.; tmpwind.high = 800.; tmpwind.name = "750_800";
  windows.push_back(tmpwind);
  tmpwind.low = 800.; tmpwind.high = 900.; tmpwind.name = "800_900";
  windows.push_back(tmpwind);
  tmpwind.low = 900.; tmpwind.high = 1000.; tmpwind.name = "900_1000";
  windows.push_back(tmpwind);
  tmpwind.low = 1000.; tmpwind.high = 1200.; tmpwind.name = "1000_1200";
  windows.push_back(tmpwind);
  tmpwind.low = 1200.; tmpwind.high = 1800.; tmpwind.name = "1200_1800";
  windows.push_back(tmpwind);
  tmpwind.low = 1800.; tmpwind.high = 2500.; tmpwind.name = "1800_2500";
  windows.push_back(tmpwind);
  tmpwind.low = 2500.; tmpwind.high = 3500.; tmpwind.name = "2500_3500";
  windows.push_back(tmpwind);
  tmpwind.low = 3500.; tmpwind.high = 4500.; tmpwind.name = "3500_4500";
  windows.push_back(tmpwind);
  tmpwind.low = 4500.; tmpwind.high = 6000.; tmpwind.name = "4500_6000";
  windows.push_back(tmpwind);

  //Order
  std::map<std::string, std::vector<int> > modelorder; //[model][order] for cats
  //These are set from the bkgFit.cc study
  modelorder["pow"].push_back(4);
  modelorder["pow"].push_back(4);

  modelorder["expow"].push_back(1);
  modelorder["expow"].push_back(1);

  modelorder["invpow"].push_back(2);
  modelorder["invpow"].push_back(2);

  modelorder["invpowlin"].push_back(2);
  modelorder["invpowlin"].push_back(2);

  modelorder["moddijet"].push_back(1);
  modelorder["moddijet"].push_back(1);

  modelorder["dijet"].push_back(1);
  modelorder["dijet"].push_back(1);

  std::map<std::string , TFile *> fin; //[model][file]
  std::map<std::string , RooWorkspace* > win; //[model][ws]

  RooRealVar* nBackground[NCAT];

  for (auto model : models){

    std::cout << "MODEL " << model << std::endl; 

    fin[model] = TFile::Open(Form("%s/bkg_%s_%s.root", ws_dir.c_str(), model.c_str(), year.c_str()));
    fin[model]->cd();
    win[model] = (RooWorkspace*) fin[model]->Get("HLFactory_ws");

    RooRealVar* mgg = win[model]->var("mgg");

    for (unsigned int c = 0; c < cats.size(); ++c) {

      std::cout << "Cat " << cats[c] << std::endl; 

      Float_t minMassFit = 0.;
      Float_t maxMassFit = 0.;
      if (c==0){ 
	minMassFit = MINmass;
	maxMassFit = MAXmass;
      } else if (c==1){
	minMassFit = MINmassBE;
	maxMassFit = MAXmass;    
      }

      RooAbsPdf* thetmppdf = (RooGenericPdf*) win[model]->pdf(TString::Format("PhotonsMassBkg_%s%d_%s",model.c_str(), modelorder[model][c], cats[c].c_str())); 

      RooArgSet* theparams = new RooArgSet();
      theparams = thetmppdf->getParameters(*mgg);
      theparams->printLatex();
      win[model]->loadSnapshot(TString::Format("%s%d_fit_params_cat%s",model.c_str(), modelorder[model][c], cats[c].c_str())); 
      std::cout << "ALWAYS CHECK THAT THE PARAMETERS LOADED ARE THE ONE EXPECTED" << std::endl; 
      theparams->printLatex();

      nBackground[c] = (RooRealVar*) win[model]->var(TString::Format("PhotonsMassBkg_%s%d_%s_norm",model.c_str(), modelorder[model][c], cats[c].c_str()));

      RooArgSet *set_mgg;

      mgg->setRange("origRange", minMassFit, maxMassFit);
      set_mgg = new RooArgSet(*mgg);
      double nominalnorm = nBackground[c]->getVal();
      double renorm = thetmppdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,RooFit::Range("origRange"))->getVal() / nominalnorm;
      
      //Time for the number of events predicted by the alternative true underline distribution in windows
      for (auto wind : windows){
	window tmpwi;
	tmpwi.name = wind.name;
	tmpwi.low = wind.low;
	tmpwi.high = wind.high;

	mgg->setRange(wind.name.c_str(), wind.low, wind.high);
	set_mgg = new RooArgSet(*mgg);

       	tmpwi.norm = (thetmppdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,RooFit::Range(tmpwi.name.c_str()))->getVal() / renorm  )  ; 

	tmpwi.truemodel = model;
	
	std::cout<< " renorm " << renorm << " range " << tmpwi.name << " norm orig " << thetmppdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,RooFit::Range("origRange"))->getVal() << " norm region "<< thetmppdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,RooFit::Range(tmpwi.name.c_str()))->getVal() << " true norm in window " << tmpwi.norm << std::endl;

	truenormwinds[cats[c]].push_back( tmpwi );
      }

    } //end of loop over categories
    
  } //end of loop over models

  return truenormwinds;

}

//-----------------------------------------------------------------------------------
std::vector<thepdf> throwtoys(RooWorkspace* w, const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, std::vector<std::string> models, int ntoystart, int ntoyend, std::map<std::string, std::vector<window> > windows){

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  
  //Order
  std::map<std::string, std::vector<int> > modelorder; //[model][order] for cats
  //These are set from the bkgFit.cc study
  modelorder["pow"].push_back(4);
  modelorder["pow"].push_back(4);

  modelorder["expow"].push_back(1);
  modelorder["expow"].push_back(1);

  modelorder["invpow"].push_back(2);
  modelorder["invpow"].push_back(2);

  modelorder["invpowlin"].push_back(2);
  modelorder["invpowlin"].push_back(2);

  modelorder["moddijet"].push_back(1);
  modelorder["moddijet"].push_back(1);

  modelorder["dijet"].push_back(1);
  modelorder["dijet"].push_back(1);

  std::map<std::string , TFile *> fin; //[model][file]
  std::map<std::string , RooWorkspace* > win; //[model][ws]
  std::vector<thepdf> pdfs;

  std::vector<RooDataSet*> data_toys; 

  RooRealVar* nBackground[NCAT];

  for (auto model : models){

    std::cout << "MODEL " << model << std::endl; 

    fin[model] = TFile::Open(Form("%s/bkg_%s_%s.root", ws_dir.c_str(), model.c_str(), year.c_str()));
    fin[model]->cd();
    win[model] = (RooWorkspace*) fin[model]->Get("HLFactory_ws");

    RooPlot* plotPhotonsMassBkg_toy[NCAT];
    RooRealVar* mgg = w->var("mgg");
    mgg->setBins(nBinsMass);

    for (unsigned int c = 0; c < cats.size(); ++c) {

      std::cout << "Cat " << cats[c] << std::endl; 

      data_toys.clear();

      Float_t minMassFit = 0.;
      Float_t maxMassFit = 0.;
      if (c==0){ 
	minMassFit = MINmass;
	maxMassFit = MAXmass;
      } else if (c==1){
	minMassFit = MINmassBE;
	maxMassFit = MAXmass;    
      }

      thepdf tmppdf; 
      tmppdf.cat = cats[c];
      tmppdf.catnum = c;
      tmppdf.model = model; 
      tmppdf.pdf = (RooGenericPdf*) win[model]->pdf(TString::Format("PhotonsMassBkg_%s%d_%s",model.c_str(), modelorder[model][c], cats[c].c_str()));
      tmppdf.order = modelorder[model][c];
      // tmppdf.ws = win[model];

      // //Some debugging print out
      // std::cout << << << std::endl; 
      //Parameters loaded should be the one expected. 
      std::cout << tmppdf.pdf->GetName() << std::endl; 

      RooArgSet* theparams = new RooArgSet();
      theparams = tmppdf.pdf->getParameters(*mgg);
      theparams->printLatex();

      // RooArgSet * thesnap = (RooArgSet *) win[model]->getSnapshot(TString::Format("%s%d_fit_params_cat%s",model.c_str(), modelorder[model][c], cats[c].c_str()));

      // TIterator* iter = thesnap->createIterator() ;
      // RooAbsArg* arg ;
      // while((arg=(RooAbsArg*)iter->Next())) {
      // 	arg->Print();
      // 	// arg->getVal();
      // }

      win[model]->loadSnapshot(TString::Format("%s%d_fit_params_cat%s",model.c_str(), modelorder[model][c], cats[c].c_str())); 
      // theparams->assignValueOnly(*thesnap);
      // *theparams = *thesnap ;
      // thesnap->printLatex();
      std::cout << "ALWAYS CHECK THAT THE PARAMETERS LOADED ARE THE ONE EXPECTED" << std::endl; 
      theparams->printLatex();

      nBackground[c] = (RooRealVar*) win[model]->var(TString::Format("PhotonsMassBkg_%s%d_%s_norm",model.c_str(), modelorder[model][c], cats[c].c_str()));

      tmppdf.nominalnorm = nBackground[c]->getVal();
      
      tmppdf.truenorms = windows[cats[c]];
      
      // thetmppdf->Print();
      std::cout << tmppdf.nominalnorm << " " << gRandom->Poisson(  nBackground[c]->getVal() ) << std::endl; 
      // //Back to the full range
      // mgg->setRange("origRange", minMassFit, maxMassFit);
      // set_mgg = new RooArgSet(*mgg);
      
      // fout[model]->cd();
    
      for (int toy = ntoystart; toy < ntoyend; ++toy){
	// RooDataSet* data = tmppdf.pdf->generate( RooArgSet(*mgg) , TRandom::Poisson(  nBackground[c]->getVal() ) );
	gRandom->SetSeed(0);
	RooDataSet* data = tmppdf.pdf->generate( RooArgSet(*mgg) , gRandom->Poisson(  nBackground[c]->getVal() ) );
	// RooDataSet* data = tmppdf.pdf->generate( RooArgSet(*mgg) ,   nBackground[c]->getVal() );
	data->SetName(Form("toy_%s_%s_%d",model.c_str(), cats[c].c_str(), toy));
	data->SetTitle(Form("toy_%s_%s_%d",model.c_str(), cats[c].c_str(), toy));
	data_toys.push_back(data);
	
	// data->Write();

	w->import(*data);
	
      }

      //====================================================================
      //Debugging plot for the last toy
      TCanvas* ctmp = new TCanvas(Form("ctmp_%d", c),"PhotonsMass Background From Toy");
      TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661,"brNDC");
          
      nBinsMass=(int) round( (maxMassFit-minMassFit)/20. );

      plotPhotonsMassBkg_toy[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);
      data_toys.back()->plotOn(plotPhotonsMassBkg_toy[c], MarkerColor(c+2), LineColor(c+2) );  
      if (c==0){plotPhotonsMassBkg_toy[c]->Draw();}
      else if (c==1){plotPhotonsMassBkg_toy[c]->Draw("same");}
      legdata->AddEntry(plotPhotonsMassBkg_toy[c]->getObject(0),Form("Data_Toy_%s",cats[c].c_str()),"LPE");

      legdata->SetTextSize(0.035);
      legdata->SetTextFont(42);
      // legdata->SetTextAlign(31);
      legdata->SetBorderSize(0);
      legdata->SetFillStyle(0);
      legdata->Draw("same");

      ctmp->SetLogy();
      ctmp->SaveAs( Form("%s/Bkg_Toy_%s_%s_log.png", ws_dir.c_str(), cats[c].c_str() ,year.c_str() ) );
      ctmp->SetLogy(0);
      ctmp->SaveAs( Form("%s/Bkg_Toy_%s_%s.png", ws_dir.c_str(), cats[c].c_str(), year.c_str()) );
      //====================================================================
      
      tmppdf.toys = data_toys;
      pdfs.push_back(tmppdf);

    }//end of loop over categories

    // delete win[model];
    // fin[model]->Close();
     // fout[model]->Close();
 
  }//end of loop over models

  return pdfs;

}
//-----------------------------------------------------------------------------------
std::vector<biasfitres> fitToys(RooWorkspace* w,std::vector<thepdf> thepdfwithtoys, std::string familymodel, bool blind, bool dobands, const std::string &year, const std::string &ws_dir, const std::string &out_dir, std::vector<std::string> cats, std::vector<std::string> models, bool approx_minos){

  // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  // Int_t ncat = NCAT;

  std::vector<biasfitres> thebiasresult; 
  std::vector<theFitResult> fitresult; 
  
  std::vector<TCanvas*> ctmp;

  RooArgList* coefList = new RooArgList();
  std::string formula = ""; 

  Float_t minMassFit = 0.;
  Float_t maxMassFit = 0.;

  RooRealVar* mgg = w->var("mgg");   

  for (auto pdf : thepdfwithtoys){

    std::cout << "===============================================================" << std::endl; 
    std::cout << "TOYS FROM MODEL " << pdf.model << std::endl; 

    biasfitres tmpbfr;
    tmpbfr.cat = pdf.cat; 
    tmpbfr.model = pdf.model; 

    // pdf.ws->Print("V");
 
    if (pdf.catnum==0){ 
      minMassFit = MINmass;
      maxMassFit = MAXmass;
    } else if (pdf.catnum==1){
      minMassFit = MINmassBE;
      maxMassFit = MAXmass;    
    }

    mgg->setUnit("GeV");
    mgg->setRange(minMassFit,maxMassFit);
    //Will use roughly the same bins as in the past (20 GeV here) although for the non reso 
    //analysis they use 25 GeV.  
    nBinsMass = (int) round( (maxMassFit-minMassFit)/20. );
    mgg->setBins(nBinsMass);

    coefList = new RooArgList(*mgg, 
    			      *w->var(TString::Format("PhotonsMass_bkg_dijet_linc_cat%d", pdf.catnum)), 
    			      *w->var(TString::Format("PhotonsMass_bkg_dijet_logc_cat%d", pdf.catnum)));
    
    formula = "TMath::Max(1e-50,pow(@0,@1+@2*log(@0)))";
    
    //The toys thrown from the alternative models will be fit with the dijet ansantz
    //So here we build the dijet ansantz
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(Form("PhotonsMassBkg_dijet_%s", cats[pdf.catnum].c_str()), Form("PhotonsMassBkg_dijet_%s", cats[pdf.catnum].c_str()), formula.c_str(),  *coefList );

    RooDataSet* dset;
    RooDataSet* edset;
    RooArgSet *set_mgg;

    for (auto toy : pdf.toys){

      // toy = (RooDataSet*) w->data(toy->GetName());
      
      //Before studying each toy we should clear the fitresult parameter
      //Update: No need to clear. I will use the back operator. 
      // fitresult.clear();    
      // ctmp.clear();          

      std::string thetoyname = toy->GetName();
      std::cout << "TOY " << thetoyname << std::endl;
  
      std::size_t found = thetoyname.find_last_of("_");
      int toynum = std::stoi(thetoyname.substr(found+1));
            
      // ==========================================================
      //Here we obtain ghat(mgg) after a fit to the toy from h(mgg). 
      int fitStatus = 0;
      double thisNll = 0.;
      std::cout << "Fitting true underline model " <<  pdf.model.c_str() << " with " << PhotonsMassBkgTmp0->GetName() << std::endl; 
      std::cout << "True underline model (toy) sumEntries " << toy->sumEntries() << std::endl; 
      fitresult.push_back( theFit(PhotonsMassBkgTmp0, NULL, toy, &thisNll, &fitStatus, /*max iterations*/ 3, false ) );
      
      std::cout << TString::Format("Background Fit results: Toy from %s cat %s", pdf.model.c_str(), pdf.cat.c_str()) << std::endl;
      //I use back because I do not know if pdf.catnum follows the 0 to 1 order for categories. 
      fitresult.back().fitres->Print("V");

      //************************************************
      // PhotonsMass background fit results per categories 
      // TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
      ctmp.push_back( new TCanvas(Form("ctmp_%s_%s",pdf.model.c_str(), toy->GetName()),"PhotonsMass Background Categories")  );

      gof gofresult = PlotFitResult(w, ctmp.back(), pdf.catnum, mgg, toy, "dijet", PhotonsMassBkgTmp0, minMassFit, maxMassFit, blind, dobands, fitresult.back().fitres->floatParsFinal().getSize(), year, ws_dir, 0);
      fitresult.back().gofresults = gofresult;
 
      ctmp.back()->SaveAs( Form("%s/Bkg_%s_%s.png", ws_dir.c_str(),toy->GetName(), year.c_str() ) );
      // ==========================================================

      RooArgSet* theparams = PhotonsMassBkgTmp0->getParameters(*mgg);
      w->saveSnapshot(Form("dijet_fit_params_cat%s", cats[pdf.catnum].c_str() ), *theparams, kTRUE);

      //For the pull
      RooRealVar *roonorm;
      RooFitResult *fitTest;
      for (auto wind : pdf.truenorms){

	//load the succesful fit param values in the whole range
	w->loadSnapshot(Form("dijet_fit_params_cat%s", cats[pdf.catnum].c_str() ) );
	theparams->printLatex();

	mgg->setRange(wind.name.c_str(), wind.low, wind.high);
	set_mgg = new RooArgSet(*mgg);
	RooAbsReal * integral = PhotonsMassBkgTmp0->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg), wind.name.c_str());
	dset = (RooDataSet*) toy->reduce(Form("mgg > %f && mgg < %f", wind.low, wind.high ) );

	double nomnorm = integral->getVal()*dset->sumEntries();
	if (nomnorm == 0.) {continue;}
	std::cout << "nomnorm " << nomnorm << " Range " << wind.name << " true norm " << wind.norm << std::endl;

	double largeNum = std::max(0.1,nomnorm*50.); 
	roonorm =  new RooRealVar(Form("norm_%s_%s_%s", pdf.model.c_str(), pdf.cat.c_str(),  wind.name.c_str()), Form("norm_%s_%s_%s", pdf.model.c_str(), pdf.cat.c_str(),  wind.name.c_str()) , nomnorm,-largeNum,largeNum);
	roonorm->setConstant(false);

	RooExtendPdf *epdf = new RooExtendPdf(Form("%s_%s_%s", pdf.model.c_str(), pdf.cat.c_str(), wind.name.c_str()), Form("%s_%s_%s", pdf.model.c_str(), pdf.cat.c_str(), wind.name.c_str()), *PhotonsMassBkgTmp0, *roonorm, wind.name.c_str());

	edset = (RooDataSet*) toy;

	// ======================================================
	// fitStatus = 0;
	// thisNll = 0.;
	// std::cout << "Fitting true underline model " <<  pdf.model.c_str() << " with " << PhotonsMassBkgTmp0->GetName() << std::endl; 
	// std::cout << "True underline model (toy) sumEntries (full range)" << toy->sumEntries() << std::endl; 
	// fitresult.push_back( theFit(PhotonsMassBkgTmp0, epdf, edset, &thisNll, &fitStatus, /*max iterations*/ 3, true ) );
      
	// std::cout << TString::Format("Background Fit results: Toy from %s cat %s", pdf.model.c_str(), pdf.cat.c_str()) << std::endl;
	// //I use back because I do not know if pdf.catnum follows the 0 to 1 order for categories. 
	// fitresult.back().fitres->Print("V");
	// ======================================================


	fitTest = new RooFitResult("","");

	RooAbsReal *nll = epdf->createNLL( *edset, Extended());
	RooMinimizer *minim = new RooMinimizer(*nll);
	minim->setEps(1000);
	minim->setOffsetting(true);
	minim->setMinimizerType("Minuit2");
	minim->setMaxIterations(15);
	minim->setMaxFunctionCalls(100);
	minim->setPrintLevel(0);
	std::cout << "Running minimize" << std::endl;
	//Strategy { Speed =0, Balance =1, Robustness =2 }
	minim->setStrategy(1);
	minim->minimize("Minuit2","minimize");
	fitTest = minim->save("fitTest","fitTest");

	// std::cout << "Running migrad" << std::endl;
	// int migrad = minim->migrad();
	// std::cout << "int migrad" << migrad << std::endl;
	// if (migrad != 0){ continue; }

	int minim_status = fitTest->status();
	// if (minim_status!=0){
	//   minim->setStrategy(0);
	//   minim->minimize("Minuit2","minimize");
	//   fitTest = minim->save("fitTest","fitTest");
	// }

	std::cout << "---------------------------------------------------" << std::endl;
	std::cout << "int minim_status" << minim_status << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;
	if (minim_status != 0){ continue; }

	nomnorm = roonorm->getVal();
                        
	std::cout << "Now running hesse" << std::endl;
	int hes = minim->hesse();
	std::cout << "---------------------------------------------------" << std::endl;
	std::cout << "int hesse" << hes << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;
	double hesseerr = roonorm->getError();
	double fiterrh  = roonorm->getErrorHi();
	double fiterrl  = roonorm->getErrorLo();

	int minos = -999;
	if (!approx_minos){
	  std::cout << "Running minos" << std::endl;
	  
	  minos = minim->minos( RooArgSet(*roonorm) ) ;
	  if (minos == 0){ 
	    if (roonorm->getErrorHi() != 0.){ fiterrh = roonorm->getErrorHi();}
	    if (roonorm->getErrorLo() != 0.){ fiterrl = roonorm->getErrorLo();}
	  }
	} else {
	  std::cout << "Computing approximate minos errors" << std::endl;
	  double fitval  = roonorm->getVal();
	  fiterrh = abs(roonorm->getErrorHi()/2.);
	  fiterrl = abs(roonorm->getErrorLo()/2.);

	  std::cout << "Computing NLL at minimum" << std::endl;
	  double minll  = nll->getVal();
	  if (fiterrl < fitval){ roonorm->setVal(fitval-fiterrl);}
	  else { roonorm->setVal(0.1); fiterrl = fitval - 0.1;}

	  std::cout << "evaluating NLL at " << roonorm->getVal() << std::endl;
	  roonorm->setConstant(true); 
	  RooMinimizer *minimm = new RooMinimizer(*nll);
	  minimm->setPrintLevel( 0 );
	  minimm->setMaxIterations(15);
	  minimm->setMaxFunctionCalls(100);                            
	  minimm->setStrategy(1);
	  minimm->setEps(1000);
	  minimm->setOffsetting(true);
	  minimm->setMinimizerType("Minuit2");
	  minimm->migrad();
	  double nllm = nll->getVal();
    
	  roonorm->setVal(fitval+fiterrh);
	  roonorm->setConstant(true); 
	
	  std::cout << "evaluating NLL at " << roonorm->getVal() << std::endl;
	  RooMinimizer *minimp = new RooMinimizer(*nll);
	  minimp->setPrintLevel( 0 );
	  minimp->setMaxIterations(15);
	  minimp->setMaxFunctionCalls(100);                            
	  minimp->setStrategy(1);
	  minimp->setEps(1000);
	  minimp->setOffsetting(true);
	  minimp->setMinimizerType("Minuit2");
	  minimp->migrad();
	  double nllp =  nll->getVal();
                
	  if ( (nllm-minll > 0.) && (nllp-minll > 0.) ){
	    fiterrh = std::max(hesseerr,fiterrh / sqrt(2.*(nllp-minll))); 
	    fiterrl = std::max(hesseerr,fiterrl / sqrt(2.*(nllm-minll)));
	    minos = 0;
	  } else {
	    minos = 1;
	  }

	  std::cout << "approx minos status " << minos << std::endl;

	}//end of approx minos
	double errh = (fiterrh != 0.) ? fiterrh : hesseerr;
	double errl = (fiterrl != 0.) ? fiterrl : hesseerr;
	
	double bias = 0.;
	if (nomnorm > wind.norm){                            
	  bias = (nomnorm-wind.norm)/abs(errl);
	} else{
	  bias = (nomnorm-wind.norm)/abs(errh);
	}

	//Will need to save the info for later processing. 
	//If the code reached here it means that it didn't stack and 
	//we can have some output, avoiding corrupted root files. 

	tmpbfr.toynum = toynum;
	tmpbfr.windname = wind.name;
	tmpbfr.windnorm = wind.norm;
	tmpbfr.nomnorm = nomnorm;
	tmpbfr.minos = minos;
	tmpbfr.hesseerr = hesseerr;
	tmpbfr.fiterrh = fiterrh;
	tmpbfr.fiterrl = fiterrl;
	tmpbfr.bias = bias;
	tmpbfr.minMassFit = minMassFit;
	tmpbfr.maxMassFit = maxMassFit;

	thebiasresult.push_back(tmpbfr);

      } //end of loop over windows

    } //end of loop over toys
    
  }//end of loop over pdfs 
  
  return thebiasresult;

}

//-----------------------------------------------------------------------------------
void makebiasrootfile(std::vector<biasfitres> bfr, const std::string &out_dir, std::string familymodel, std::string year, int ntoystart, int ntoyend){
  
  TFile* biasesfiles = new TFile(Form("%s/%s/tree_bias_%s_%s_%d_%d.root", out_dir.c_str(), familymodel.c_str(), familymodel.c_str(), year.c_str(), ntoystart, ntoyend ), "recreate");
 
  //We will create the ntuple to hold the bias results. 
  std::map<std::string, TNtuple* > biases; //[cat+windows][ntuple]
  biases.clear();   

  std::vector<window> windows;
  window tmpwind;
  tmpwind.low = 500.; tmpwind.high = 550.; tmpwind.name = "500_550";
  windows.push_back(tmpwind);
  tmpwind.low = 550.; tmpwind.high = 600.; tmpwind.name = "550_600";
  windows.push_back(tmpwind);
  tmpwind.low = 600.; tmpwind.high = 650.; tmpwind.name = "600_650";
  windows.push_back(tmpwind);
  tmpwind.low = 650.; tmpwind.high = 700.; tmpwind.name = "650_700";
  windows.push_back(tmpwind);
  tmpwind.low = 700.; tmpwind.high = 750.; tmpwind.name = "700_750";
  windows.push_back(tmpwind);
  tmpwind.low = 750.; tmpwind.high = 800.; tmpwind.name = "750_800";
  windows.push_back(tmpwind);
  tmpwind.low = 800.; tmpwind.high = 900.; tmpwind.name = "800_900";
  windows.push_back(tmpwind);
  tmpwind.low = 900.; tmpwind.high = 1000.; tmpwind.name = "900_1000";
  windows.push_back(tmpwind);
  tmpwind.low = 1000.; tmpwind.high = 1200.; tmpwind.name = "1000_1200";
  windows.push_back(tmpwind);
  tmpwind.low = 1200.; tmpwind.high = 1800.; tmpwind.name = "1200_1800";
  windows.push_back(tmpwind);
  tmpwind.low = 1800.; tmpwind.high = 2500.; tmpwind.name = "1800_2500";
  windows.push_back(tmpwind);
  tmpwind.low = 2500.; tmpwind.high = 3500.; tmpwind.name = "2500_3500";
  windows.push_back(tmpwind);
  tmpwind.low = 3500.; tmpwind.high = 4500.; tmpwind.name = "3500_4500";
  windows.push_back(tmpwind);
  tmpwind.low = 4500.; tmpwind.high = 6000.; tmpwind.name = "4500_6000";
  windows.push_back(tmpwind);
  
  for (auto wind : windows){
    biases[ Form("EBEB_%s", wind.name.c_str())] = new TNtuple(Form("tree_bias_EBEB_%s_%s", familymodel.c_str(), wind.name.c_str()), Form("tree_bias_EBEB_%s_%s", familymodel.c_str(),wind.name.c_str()), "toy:truth:fit:minos:errhe:errp:errm:bias:fitmin:fitmax" );
    biases[ Form("EBEE_%s", wind.name.c_str())] = new TNtuple(Form("tree_bias_EBEE_%s_%s", familymodel.c_str(), wind.name.c_str()), Form("tree_bias_EBEE_%s_%s", familymodel.c_str(),wind.name.c_str()), "toy:truth:fit:minos:errhe:errp:errm:bias:fitmin:fitmax" );
  }

  for (auto res : bfr){
    biases[Form("%s_%s", res.cat.c_str(), res.windname.c_str())]->Fill( res.toynum, res.windnorm, res.nomnorm,  res.minos, res.hesseerr, res.fiterrh, res.fiterrl, res.bias, res.minMassFit, res.maxMassFit );
  }

  biasesfiles->Write();
  biasesfiles->Close();

}
//-----------------------------------------------------------------------------------
void analyzeBias(const std::string &ws_dir, const std::string &out_dir, std::vector<std::string> models, std::vector<std::string> cats, std::string year, int ntoystart, int ntoyend){
  
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptStat(1111);
  m_gStyle->SetOptFit(1);

  std::map<std::string , TFile *> fin; //[model][file]
  std::map<TString , TH1F *> hb;  
  std::map<TString , TH1F *> hc;  
  std::map<TString , TH1F *> hd;  
  std::map<TString , TH1F *> hcb;  
  std::map<TString , TF1 *> gaus;  
  std::map<TString , TF1 *> gausd;  
  std::map<TString , TF1 *> gauscb;  
  std::map<TString , TF1 *> bias_func; //[cat]
  std::map<TString , TString > bias_formula; //[cat+windname]
  std::map<TString , TCanvas *> bcan;  
  std::map<std::string , std::map<TString , TGraphErrors *> > profiles; //[model][cat] 
  std::map<std::string , std::map<TString , TGraphErrors *> > bprofiles; //[model][cat] 
  std::map<std::string , std::map<TString , TGraphErrors *> > cprofiles; //[model][cat] 
  std::map<std::string , std::map<TString , TGraphErrors *> > cbprofiles; //[model][cat] 

  std::map<std::string , std::map<std::string , std::vector<double> > > bpullvals, bpullvalsErr, pullvals, pullvalsErr, bcpullvals, bcpullvalsErr, windvals, windvalsErr;
  
  bias_formula = buildBiasTerm(year, cats, "EXO-17-017", false); //false is to scale to year luminosity and not to 1/fb-1
  plotBiasTerm(ws_dir, year, cats);

  float xfirst = 1e6;
  float xlast  = 0.;

  for (auto model : models){
    
    std::cout << "MODEL " << model << std::endl; 
    fin[model] = TFile::Open(Form("%s/%s/tree_bias_%s_%s_%d_%d.root", out_dir.c_str(), model.c_str(), model.c_str(), year.c_str(), ntoystart, ntoyend));
    fin[model]->cd();
    
    TDirectory *curdir = gDirectory;
    TIter next(curdir->GetListOfKeys());
    TKey *key;

    xfirst = 1e6;
    xlast  = 0.;
    const Int_t nq = 1;
    Double_t prb[nq];  // position where to compute the quantiles in [0,1], median here. 
    Double_t med[nq];  // array to contain the median 
    Double_t medd[nq];  // array to contain the median 
    Double_t medcb[nq];  // array to contain the median 
    Double_t qtl[nq];  // array to contain the quantile 

    while ((key = (TKey*)next())) {
      // TClass *cl = gROOT->GetClass(key->GetClassName());
      TString name = key->GetName();
      std::cout << name << std::endl;
      if (name.Contains("tree_bias")){

    	//The first 10 characters should be true_bias_
    	name = name(10,name.Length());
    	//From the start to the first _ we have the category. 
    	std::string cat = name(0,name.Index("_"));
    	//The window is (in the form e.g. 500_550)
    	window wind; 
    	wind.name = name(name.Index(model)+model.length()+1,name.Length());
    	wind.low =  std::stod(wind.name.substr(0, wind.name.find("_")));
    	wind.high = std::stod(wind.name.substr(wind.name.find("_")+1, wind.name.length()));

    	std::cout << "model " << model << " cat " << cat << " wind.name " << wind.name << " wind.low " << wind.low << " wind.high " << wind.high <<std::endl;

    	TNtuple *tree = (TNtuple*)key->ReadObj();
    	TString slabel = Form("%s_%s_%s", cat.c_str(), model.c_str(), wind.name.c_str());

    	bcan[slabel] = new TCanvas(slabel,slabel);
    	bcan[slabel]->cd();

    	hb[slabel] = new TH1F(Form("h_bias_%s",slabel.Data()), Form("h_bias_%s",slabel.Data()), 501, -5.005, 5.005);
    	tree->Draw(Form("bias>>%s", hb[slabel]->GetName()) );
    	hb[slabel]->Fit("gaus","L+Q");
    	Int_t nentries = hb[slabel]->GetEntries();

    	std::cout << "hb[slabel] " << nentries << std::endl; 

    	hb[slabel]->Draw();
    	bcan[slabel]->SaveAs( Form("%s/bias_%s_%s.png", ws_dir.c_str(), slabel.Data() ,year.c_str() ) );

    	// hb[slabel]->GetListOfFunctions()->Print();
    	gaus[slabel] = (TF1*) hb[slabel]->GetListOfFunctions()->At(0);
    	prb[0] = 0.5;
    	med[0] = 0.;
    	hb[slabel]->GetQuantiles(nq,med,prb);

    	hc[slabel] = new TH1F(Form("h_coverage_%s",slabel.Data()), Form("h_coverage_%s",slabel.Data()), 501, 0., 5.01);
    	tree->Draw(Form("abs(bias)>>%s", hc[slabel]->GetName()) );

    	prb[0] = 0.683;
    	qtl[0] = 0.;
    	hc[slabel]->GetQuantiles(nq,qtl,prb);

    	hd[slabel] = new TH1F(Form("h_deviation_%s",slabel.Data()), Form("h_deviation_%s",slabel.Data()), 501, -100.2, 100.2 );
    	tree->Draw(Form("fit-truth>>%s", hd[slabel]->GetName()) );
    	hd[slabel]->Fit("gaus","L+Q");

    	gausd[slabel] = (TF1*) hd[slabel]->GetListOfFunctions()->At(0);
    	medd[0] = 0.;
    	hd[slabel]->GetQuantiles(nq,medd,prb);

    	xfirst = std::min(wind.low,xfirst);
    	xlast = std::max(wind.high,xlast);

    	pullvals[model][cat].push_back( fabs(medd[0])/(wind.high-wind.low) );
    	pullvalsErr[model][cat].push_back( gaus[slabel]->GetParameter(2)/(wind.high-wind.low) );

    	bpullvals[model][cat].push_back( gaus[slabel]->GetParameter(1)/(gaus[slabel]->GetParameter(2)) );
    	windvals[model][cat].push_back( 0.5*(wind.high+wind.low) );
    	bpullvalsErr[model][cat].push_back( 0. );
    	windvalsErr[model][cat].push_back( 0.5*(wind.high-wind.low) );

    	tree->GetEntry(0);
    	//Corrected bias
    	TString formulakey = Form("%s_%s", cat.c_str(), wind.name.c_str());
    	// bias_func[slabel] = new TF1( Form("err_correction_%s", slabel.Data()), bias_formula[formulakey], 0., 2e+6);
    	std::cout << formulakey << " " << bias_formula[formulakey] << std::endl;
    	bias_func[slabel] = new TF1( Form("err_correction_%s", slabel.Data()), Form("%s", bias_formula[formulakey].Data()), 0., 2e+6);

    	tree->SetAlias("berr",Form("(fit-truth)/bias*%f" , std::max(1.,gaus[slabel]->GetParameter(2))));
    	tree->SetAlias("corr_bias",Form("(fit-truth)/sqrt(berr^2+%f^2)", (bias_func[slabel]->Integral(wind.low,wind.high)*1.0) ) );

    	hcb[slabel] = new TH1F(Form("h_corr_bias_%s",slabel.Data()), Form("h_corr_bias_%s",slabel.Data()), 501, -5.005, 5.005);
    	tree->Draw(Form("corr_bias>>%s", hcb[slabel]->GetName()) );
    	hcb[slabel]->Fit("gaus","L+Q");

    	gauscb[slabel] = (TF1*) hcb[slabel]->GetListOfFunctions()->At(0);
    	medcb[0] = 0.;
    	hcb[slabel]->GetQuantiles(nq,medcb,prb);

    	bcpullvals[model][cat].push_back( gauscb[slabel]->GetParameter(1) );
    	bcpullvalsErr[model][cat].push_back( 0. );
	  

      } //end of loop through bias ntuples

    }//end of loop through keys of file
  }//end of loop through models

  //Time for plot
  m_gStyle->SetOptFit(0);

  std::map<TString , int > colors; 
  colors["expow"] = 2; //kRed
  colors["invpow"] = 4; //kBlue
  colors["invpowlin"] = 417; //kGreen+1
  colors["moddijet"] = 433; //kCyan+1
  colors["pow"] = 616; //kMagenta
  // colors[""] = ; //kYellow
  // colors[""] = ; //kGray
  // colors[""] = ; //kOrange
  
  //For pull
  std::map<std::string , TCanvas *> bcanv; //[cat]
  std::map<std::string ,TLegend *> bleg; //[cat]
  //For bias
  std::map<std::string , TCanvas *> canv; //[cat]
  std::map<std::string ,TLegend *> leg; //[cat]
  //For corrected pull bias
  std::map<std::string , TCanvas *> cbcanv; //[cat]
  std::map<std::string ,TLegend *> cbleg; //[cat]

  std::map<std::string , TH2F *> frame;
  std::map<std::string , TBox *> box;

  //First loop through cats since we want bias plots per category. 
  for (auto cat : cats){
    
    bcanv[cat] = new TCanvas(Form("profile_pull_%s",cat.c_str()),Form("profile_pull_%s",cat.c_str()));
    bcanv[cat]->SetLogx();
    bcanv[cat]->SetGridy();
    bcanv[cat]->SetGridx();

    bleg[cat]= new TLegend(0.6,0.12,0.95,0.47);
    bleg[cat]->SetFillStyle(0);
      
    frame[cat] = new TH2F(Form("frame_%s",cat.c_str()),Form("frame_%s",cat.c_str()),100,xfirst,xlast,100,-4,2);
    frame[cat]->SetStats(false);
    bcanv[cat]->cd();
    frame[cat]->Draw();
    frame[cat]->GetXaxis()->SetTitle("mass [GeV]");
    frame[cat]->GetXaxis()->SetMoreLogLabels();
    frame[cat]->GetYaxis()->SetTitle("( n_{fit} - n_{true} )/ #sigma_{fit}");

    box[cat]= new TBox(xfirst,-0.5,xlast,0.5);
    box[cat]->SetFillColor(kGray);
    box[cat]->Draw("same");       

    //------------------------------------------------------------------------------------------------------
    //First with the pull bias plot
    for (auto model : models){
      // A std::vector is at its heart an array. To get the array just get the address of the first element.
      bprofiles[model][cat] = new TGraphErrors( bpullvals[model][cat].size(), &windvals[model][cat][0], &bpullvals[model][cat][0], &windvalsErr[model][cat][0], &bpullvalsErr[model][cat][0] ) ; 

      bprofiles[model][cat]->GetXaxis()->SetRangeUser(xfirst,xlast);
      bprofiles[model][cat]->GetXaxis()->SetTitleOffset( 0.9 );
      bprofiles[model][cat]->SetMarkerColor(colors[model]);
      bprofiles[model][cat]->SetMarkerStyle(kFullCircle);
      bleg[cat]->AddEntry( bprofiles[model][cat] , Form("%s",model.c_str()) ,"pe" );
      bleg[cat]->SetHeader(Form("dijet %s",cat.c_str()));
      bprofiles[model][cat]->Draw("PSE");
      bleg[cat]->Draw("same");
      bcanv[cat]->RedrawAxis();
      bcanv[cat]->Modified();
      bcanv[cat]->Update();

    }//end of loop through models
      
    bcanv[cat]->SaveAs((Form("%s/profile_pull_%s_%s.png", ws_dir.c_str(), cat.c_str(), year.c_str())) );

    //------------------------------------------------------------------------------------------------------
    //Now with the deviation, nfit-ntrue/GeV
    canv[cat] = new TCanvas(Form("profile_bias_%s",cat.c_str()),Form("profile_bias_%s",cat.c_str()));
    canv[cat]->SetLogx();
    canv[cat]->SetGridy();
    canv[cat]->SetGridx();

    leg[cat]= new TLegend(0.6,0.12,0.95,0.47);
    leg[cat]->SetFillStyle(0);
      
    frame[cat]->SetStats(false);
    canv[cat]->cd();
    frame[cat]->Draw();
    frame[cat]->GetXaxis()->SetTitle("mass [GeV]");
    frame[cat]->GetXaxis()->SetMoreLogLabels();
    frame[cat]->GetYaxis()->SetTitle("| n_{fit} - n_{true} | / GeV");

    box[cat]->SetFillColor(kGray);
    box[cat]->Draw("same");       

    for (auto model : models){

      // A std::vector is at its heart an array. To get the array just get the address of the first element.
      profiles[model][cat] = new TGraphErrors( pullvals[model][cat].size(), &windvals[model][cat][0], &pullvals[model][cat][0], &windvalsErr[model][cat][0], &pullvalsErr[model][cat][0] ) ; 

      profiles[model][cat]->GetXaxis()->SetRangeUser(0.001, 0.3);
      profiles[model][cat]->GetXaxis()->SetTitleOffset( 0.9 );
      profiles[model][cat]->SetMarkerColor(colors[model]);
      profiles[model][cat]->SetMarkerStyle(kFullCircle);
      leg[cat]->AddEntry( profiles[model][cat] , Form("%s",model.c_str()) ,"pe" );
      leg[cat]->SetHeader(Form("dijet %s",cat.c_str()));
      profiles[model][cat]->Draw("PSE");
      leg[cat]->Draw("same");
      canv[cat]->RedrawAxis();
      canv[cat]->Modified();
      canv[cat]->Update();

    }//end of loop through models

    canv[cat]->SaveAs((Form("%s/profile_bias_%s_%s.png", ws_dir.c_str(), cat.c_str(), year.c_str())) );

   //------------------------------------------------------------------------------------------------------
    //Now with the corrected pull, nfit-ntrue/GeV
    cbcanv[cat] = new TCanvas(Form("profile_corr_pull_%s",cat.c_str()),Form("profile_corr_pull_%s",cat.c_str()));
    cbcanv[cat]->SetLogx();
    cbcanv[cat]->SetGridy();
    cbcanv[cat]->SetGridx();

    cbleg[cat]= new TLegend(0.6,0.12,0.95,0.47);
    cbleg[cat]->SetFillStyle(0);
      
    frame[cat]->SetStats(false);
    cbcanv[cat]->cd();
    frame[cat]->Draw();
    frame[cat]->GetXaxis()->SetTitle("mass [GeV]");
    frame[cat]->GetXaxis()->SetMoreLogLabels();
    frame[cat]->GetYaxis()->SetTitle("( n_{fit} - n_{true} )/ ( #sigma_{fit} #oplus bias )");

    box[cat]->SetFillColor(kGray);
    box[cat]->Draw("same");       

    for (auto model : models){

      // A std::vector is at its heart an array. To get the array just get the address of the first element.
      cbprofiles[model][cat] = new TGraphErrors( bcpullvals[model][cat].size(), &windvals[model][cat][0], &bcpullvals[model][cat][0], &windvalsErr[model][cat][0], &bcpullvalsErr[model][cat][0] ) ; 

      cbprofiles[model][cat]->GetXaxis()->SetRangeUser(0.001, 0.3);
      cbprofiles[model][cat]->GetXaxis()->SetTitleOffset( 0.9 );
      cbprofiles[model][cat]->SetMarkerColor(colors[model]);
      cbprofiles[model][cat]->SetMarkerStyle(kFullCircle);
      cbleg[cat]->AddEntry( cbprofiles[model][cat] , Form("%s",model.c_str()) ,"pe" );
      cbleg[cat]->SetHeader(Form("dijet %s",cat.c_str()));
      cbprofiles[model][cat]->Draw("PSE");
      cbleg[cat]->Draw("same");
      cbcanv[cat]->RedrawAxis();
      cbcanv[cat]->Modified();
      cbcanv[cat]->Update();

    }//end of loop through models

    cbcanv[cat]->SaveAs((Form("%s/profile_corr_pull_%s_%s.png", ws_dir.c_str(), cat.c_str(), year.c_str())) );

  }//end of loop through cats
    


																			    }

//-----------------------------------------------------------------------------------
theFitResult theFit(RooAbsPdf *pdf, RooExtendPdf *epdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, bool fittoys){

  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  // params_test->printLatex();
  std::cout << "--------------------- BEFORE ITERATIONS-------------------------------" << std::endl;
  int stat=1;
  double minnll=10e8;
  // double offset=10e8;
  // double minnll_woffset=10e8;
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
    // stat = fitTest->status();
    // minnll = fitTest->minNll();

    //Including offset 
    //-------------------
    RooNLLVar *nll;
    if(!fittoys){
      nll = new RooNLLVar("nll","nll",*pdf,*data);
    } else{
      nll = new RooNLLVar("nll","nll",*epdf,*data, Extended());
    }

    RooMinimizer *minuit_fitTest = new RooMinimizer(*nll);
    minuit_fitTest->setOffsetting(kTRUE);
    minuit_fitTest->setStrategy(2);
    // minuit_fitTest->setEps(1000);
    minuit_fitTest->setPrintLevel(5);
    if (fittoys){
      minuit_fitTest->setMaxIterations(15);
      minuit_fitTest->setMaxFunctionCalls(100);                            
    }
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
  std::cout << "------------------------OFFSET-----------------------------" << std::endl;
  // std::cout << "end of runFit stat=" << stat << " offset=" << offset << " minnll with offset=" << minnll_woffset << " diff= " << minnll<< " nll->getVal() " << nllval << std::endl;
  std::cout << "end of runFit stat=" << stat << " minnll with nll->getVal() " << nllval << std::endl;
  *stat_t = stat;
  *NLL = minnll;

  return tmpfit;
}

//-----------------------------------------------------------------------------------
gof PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, std::string model, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir, int order){

  RooPlot* plotPhotonsMassBkg[NCAT];

  std::string name = Form("%s/Bkg_cat%d_%s_%s.png", ws_dir.c_str(),c, model.c_str(), year.c_str() );
  std::cout << "GOODNESS OF FIT WITH/OUT TOYS" << std::endl; 
  gof gofresult = getGoodnessOfFit(mgg, PhotonsMassBkgTmp0, data, name, gofwithtoys);
  ctmp->cd();

  //We should get the chisquare from the full range fit
  //2 GeV bins
  //nBinsMass = (int) round( (maxMassFit-minMassFit)/2. );

  //roughly 20 GeV bins
  // nBinsMass = (int) round( (maxMassFit-minMassFit)/20. );
  // plotPhotonsMassBkg[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);
  // data->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
  // PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit));
  // //RooPlot::chiSquare() should give you the chi2/ndof.
  // //We should take the chisquare from the full range fit. 
  // double chi2 = plotPhotonsMassBkg[c]->chiSquare(numoffittedparams);
  // Int_t ndof = nBinsMass-numoffittedparams;
  // std::cout<<"------> ndof"<< ndof<<std::endl;
  // //https://root.cern.ch/root/html524/TMath.html#TMath:Prob   
  // double prob = TMath::Prob(chi2*ndof, ndof);
  // std::cout << prob << std::endl;

  if( blind ) { 
    maxMassFit = 1000.; 
  }
  
  nBinsMass = (int) round( (maxMassFit-minMassFit)/20. );
  plotPhotonsMassBkg[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);

  data->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
  PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit));

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
  if (order == 0){
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),Form("Parametric Model: %s", model.c_str()),"L");
  } else{
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),Form("Parametric Model: %s Order %d", model.c_str(), order),"L");
  }
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  // legdata->SetTextAlign(31);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");

  TPaveText* label_cms = get_labelcms(1, year, false);
  TPaveText* label_sqrt = get_labelsqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  
  //write down the chi2 of the fit on the
 
  TPaveText* label_chi2 = new TPaveText(0.55,0.33,0.79,0.44, "brNDC");
  label_chi2->SetFillColor(kWhite);
  label_chi2->SetTextSize(0.035);
  label_chi2->SetTextFont(42);
  label_chi2->SetTextAlign(31); // align right
  label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", gofresult.chi2overndof));
  if (gofwithtoys){
    label_chi2->AddText(TString::Format("Chi square Prob (toys) = %.3f", gofresult.prob));
  } else{
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", gofresult.prob));
  }
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

  return gofresult;
  
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

//-----------------------------------------------------------------------------------
gof getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name, bool gofToys){
	
  gof gofresult;
  double prob;
  int ntoys = 50;
  int nBinsForMass=mass->getBinning().numBins();
  std::cout << "mass bins " << nBinsForMass << std::endl;

  // Routine to calculate the goodness of fit. 
  name+="_gofTest";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();
  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Name("data"));

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
    if ( ((double)data->sumEntries()/nBinsForMass < 5) || gofToys ){ 

    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
    //  std::cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      if (itoy%10==0){std::cout << "toy " << itoy << std::endl;}
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *toy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
      //RooDataSet *toy = pdf->generate(RooArgSet(*mass),nToyEvents,0,1);
   //   pdf->fitTo(*toy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(2),RooFit::Offset(kTRUE)); 
      pdf->fitTo(*toy,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Strategy(2),RooFit::Offset(kTRUE));

      RooPlot *plot_t = mass->frame();
      toy->plotOn(plot_t);
      pdf->plotOn(plot_t);//,RooFit::NormRange("fitdata_1,fitdata_2"));

      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForMass-np));
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    TCanvas *can = new TCanvas();
    can->cd();
	double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));}
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForMass-np),toyhist.GetMaximum(),chi2*(nBinsForMass-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(Form("%s.png",name.c_str()));
    can->SaveAs(Form("%s.pdf",name.c_str()));

    // back to best fit 	
    params->assignValueOnly(preParams);
	std::cout << "[INFO] Probability from toys " << prob << " Probability from TMath::Prob " << TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np)<< std::endl;
  	} else {  prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np); }
  std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
  std::cout << "[INFO] p-value  =  " << prob << std::endl;

  gofresult.prob = prob;
  gofresult.chi2overndof = chi2;

  delete pdf;
  return gofresult;

}

//-----------------------------------------------------------------------------------
void writeToJson(std::map<std::string, std::vector<window> > valuestowrite, std::string outputfile)
{
  // Short alias for this namespace
  namespace pt = boost::property_tree;
  
  // Create a root
  pt::ptree root;

  std::string sel_node; 
  std::string sel_subnode; 

  for (auto &valmap : valuestowrite){
    // pt::ptree values_node;
    std::string model; 
    for (auto &val : valmap.second){
      // model = val.truemodel; 
      // sel_node = val.name;

      // root.put(Form("%s.%s.windlow", model.c_str(), sel_node.c_str() ) , std::to_string(val.low));
      // root.put(Form("%s.%s.windhigh",model.c_str(), sel_node.c_str() ), std::to_string(val.high));
      // root.put(Form("%s.%s.windnorm",model.c_str(), sel_node.c_str() ), std::to_string(val.norm));
      // root.put(Form("%s.%s.windcat",model.c_str(), sel_node.c_str() ), valmap.first);

      model = val.truemodel; 
      sel_node = valmap.first;
      sel_subnode = val.name  ;

      root.put(Form("%s.%s.%s.windlow", model.c_str(), sel_node.c_str(),sel_subnode.c_str()  ) , std::to_string(val.low));
      root.put(Form("%s.%s.%s.windhigh",model.c_str(), sel_node.c_str(),sel_subnode.c_str() ), std::to_string(val.high));
      root.put(Form("%s.%s.%s.windnorm",model.c_str(), sel_node.c_str(),sel_subnode.c_str() ), std::to_string(val.norm));

    }
  }
  pt::write_json(std::cout, root);
  pt::write_json(outputfile, root);


}

//-----------------------------------------------------------------------------------
std::map<std::string, std::vector<window> > readJson(std::string inputfile, std::string model){

  // Short alias for this namespace
  namespace pt = boost::property_tree;
  using boost::property_tree::ptree;
 
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(inputfile, root);

  std::map<std::string, std::vector<window> > windows;

  std::string modelin; 
  std::string sel_node; //category
  std::string sel_subnode; //windows

  for (ptree::const_iterator it = root.begin(); it != root.end(); ++it) {

    modelin = it->first;
    if ( modelin != model){continue;}

    for (ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {

      //This is the category
      sel_node = jt->first;

      for (ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt) {

	sel_subnode = kt->first;
	window tmpwin;


	tmpwin.low = std::stod( root.get<std::string>( Form("%s.%s.%s.windlow",modelin.c_str(), sel_node.c_str(), sel_subnode.c_str() ) ) ); 
	tmpwin.high = std::stod( root.get<std::string>( Form("%s.%s.%s.windhigh",modelin.c_str(), sel_node.c_str(), sel_subnode.c_str() ) ) ); 
	tmpwin.norm    = std::stod( root.get<std::string>( Form("%s.%s.%s.windnorm",modelin.c_str(), sel_node.c_str(), sel_subnode.c_str() ) ) ) ; 
	tmpwin.name    = sel_subnode; 
	tmpwin.truemodel = modelin;
    
	windows[sel_node].push_back( tmpwin  );
	
	std::cout << tmpwin.low << " " << tmpwin.high << " " << tmpwin.norm << " "  << tmpwin.name <<  " "  << tmpwin.truemodel << " "  <<sel_node  << std::endl;
      } //end of loop over windows
    }// end of loop over categories
  }//end of loop over models
  
  return windows;

}
//-----------------------------------------------------------------------------------
//Options for paper are: EXO-17-017 or EXO-16-027
std::map<TString , TString > buildBiasTerm(std::string year, std::vector<std::string> cats, std::string paper, bool plotto1invfb){

  std::map<TString , TString > formula; 
 
  std::vector<window> windows;
  window tmpwind;
  tmpwind.low = 500.; tmpwind.high = 550.; tmpwind.name = "500_550";
  windows.push_back(tmpwind);
  tmpwind.low = 550.; tmpwind.high = 600.; tmpwind.name = "550_600";
  windows.push_back(tmpwind);
  tmpwind.low = 600.; tmpwind.high = 650.; tmpwind.name = "600_650";
  windows.push_back(tmpwind);
  tmpwind.low = 650.; tmpwind.high = 700.; tmpwind.name = "650_700";
  windows.push_back(tmpwind);
  tmpwind.low = 700.; tmpwind.high = 750.; tmpwind.name = "700_750";
  windows.push_back(tmpwind);
  tmpwind.low = 750.; tmpwind.high = 800.; tmpwind.name = "750_800";
  windows.push_back(tmpwind);
  tmpwind.low = 800.; tmpwind.high = 900.; tmpwind.name = "800_900";
  windows.push_back(tmpwind);
  tmpwind.low = 900.; tmpwind.high = 1000.; tmpwind.name = "900_1000";
  windows.push_back(tmpwind);
  tmpwind.low = 1000.; tmpwind.high = 1200.; tmpwind.name = "1000_1200";
  windows.push_back(tmpwind);
  tmpwind.low = 1200.; tmpwind.high = 1800.; tmpwind.name = "1200_1800";
  windows.push_back(tmpwind);
  tmpwind.low = 1800.; tmpwind.high = 2500.; tmpwind.name = "1800_2500";
  windows.push_back(tmpwind);
  tmpwind.low = 2500.; tmpwind.high = 3500.; tmpwind.name = "2500_3500";
  windows.push_back(tmpwind);
  tmpwind.low = 3500.; tmpwind.high = 4500.; tmpwind.name = "3500_4500";
  windows.push_back(tmpwind);
  tmpwind.low = 4500.; tmpwind.high = 6000.; tmpwind.name = "4500_6000";
  windows.push_back(tmpwind);

  double luminosityscale;
  
  if (plotto1invfb){ luminosityscale = 1.;}
  else {luminosityscale = luminosity[year];}

  for (auto cat : cats){

    for (auto wind : windows){


      //The key of the formula will be cat+windname
      TString thekey = Form("%s_%s", cat.c_str(), wind.name.c_str()); 

      if (paper == "EXO-17-017"){
      	if (cat=="EBEB"){
	  formula[thekey] = Form("x<=650. ? 0.125/%f : pow(x,2.0-0.36*log(x))/%f ", luminosityscale, luminosityscale);
      	  // //Up to 650 GeV
      	  // if (wind.high <= 650.){ formula[thekey] = Form("0.125/%f", luminosityscale); } 
      	  // //Above 650 GeV
      	  // else{ formula[thekey] = Form("(pow(x,2.0-0.36*log(x)))/%f", luminosityscale);}
      	} else if (cat=="EBEE"){
	  formula[thekey] = Form("x<=750. ? 0.2/%f : (pow(x/600.,-3.5-0.2*log(x/600.))+15e-5)/%f", luminosityscale, luminosityscale);
      	  // //Up to 750 GeV
      	  // if (wind.high <= 750.){ formula[thekey] = Form("0.2/%f", luminosityscale); } 
      	  // //Above 750 GeV
      	  // else{ formula[thekey] = Form("(pow(x/600.,-3.5-0.2*log(x))+15e-5)/%f", luminosityscale);}
      	}
      } //end of EXO-17-017 choice
      else if (paper == "EXO-16-027"){
      	if (cat=="EBEB"){ formula[thekey] = Form("pow(x,2.2-0.4*log(x))/%f", luminosityscale);}
      	else if (cat=="EBEE"){formula[thekey] = Form("(0.1*(x/600.)^(-5)+2e-5)/%f", luminosityscale);}
      } //end of EXO-16-027 choice

    }//end of loop over windows
  }//end of loop over categories

  return formula;
  
}
//-----------------------------------------------------------------------------------
//We will here plot the functions we are using. 
void plotBiasTerm(const std::string &ws_dir, std::string year, std::vector<std::string> cats){
  
  std::map<TString , TCanvas *> can;  
  std::map<std::string , TH1F *> frame;
  // std::map<std::string , TBox *> box;
  std::map<TString , TLegend *> funleg;  
  std::map<TString , TF1 * > funsEXO17017; 
  std::map<TString , TF1 * > funsEXO16027; 

  std::map<TString , TString > formEXO17017 = buildBiasTerm(year, cats, "EXO-17-017", true); //true is to scale to 1/fb-1
  std::map<TString , TString > formEXO16027 = buildBiasTerm(year, cats, "EXO-16-027", true); //true is to scale to 1/fb-1

  std::vector<window> windows;
  window tmpwind;
  tmpwind.low = 500.; tmpwind.high = 550.; tmpwind.name = "500_550";
  windows.push_back(tmpwind);
  tmpwind.low = 550.; tmpwind.high = 600.; tmpwind.name = "550_600";
  windows.push_back(tmpwind);
  tmpwind.low = 600.; tmpwind.high = 650.; tmpwind.name = "600_650";
  windows.push_back(tmpwind);
  tmpwind.low = 650.; tmpwind.high = 700.; tmpwind.name = "650_700";
  windows.push_back(tmpwind);
  tmpwind.low = 700.; tmpwind.high = 750.; tmpwind.name = "700_750";
  windows.push_back(tmpwind);
  tmpwind.low = 750.; tmpwind.high = 800.; tmpwind.name = "750_800";
  windows.push_back(tmpwind);
  tmpwind.low = 800.; tmpwind.high = 900.; tmpwind.name = "800_900";
  windows.push_back(tmpwind);
  tmpwind.low = 900.; tmpwind.high = 1000.; tmpwind.name = "900_1000";
  windows.push_back(tmpwind);
  tmpwind.low = 1000.; tmpwind.high = 1200.; tmpwind.name = "1000_1200";
  windows.push_back(tmpwind);
  tmpwind.low = 1200.; tmpwind.high = 1800.; tmpwind.name = "1200_1800";
  windows.push_back(tmpwind);
  tmpwind.low = 1800.; tmpwind.high = 2500.; tmpwind.name = "1800_2500";
  windows.push_back(tmpwind);
  tmpwind.low = 2500.; tmpwind.high = 3500.; tmpwind.name = "2500_3500";
  windows.push_back(tmpwind);
  tmpwind.low = 3500.; tmpwind.high = 4500.; tmpwind.name = "3500_4500";
  windows.push_back(tmpwind);
  tmpwind.low = 4500.; tmpwind.high = 6000.; tmpwind.name = "4500_6000";
  windows.push_back(tmpwind);

  for (auto cat : cats){

    float xfirst = 1e6;
    float xlast  = 0.;
   
    can[cat] = new TCanvas(Form("BiasTermFunction_%s",cat.c_str()),Form("BiasTermFunction_%s",cat.c_str()));
    can[cat]->cd();
    can[cat]->SetLogx();
    can[cat]->SetLogy();
    funleg[cat] = new TLegend(0.45, 0.75, 0.8, 0.9, "" , "brNDC");
    funleg[cat]->SetFillStyle(0);

    // frame[cat] = new TH2F(Form("frame_%s",cat.c_str()),Form("frame_%s",cat.c_str()),100,xfirst,xlast,100,1e-6,1e-3);
    // frame[cat] = new TH1F(Form("frame_%s",cat.c_str()),Form("frame_%s",cat.c_str()),100,1e-6,1e-3);
    // frame[cat]->SetStats(false);
    // frame[cat]->Draw();
    // frame[cat]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    // frame[cat]->GetXaxis()->SetMoreLogLabels();
    // frame[cat]->GetYaxis()->SetTitle("#beta_{#gamma#gamma} (Events/fb^{-1})");

    // box[cat]= new TBox(xfirst,-0.5,xlast,0.5);
    // box[cat]->SetFillColor(kGray);
    // box[cat]->Draw("same");       

    TString thekey = "";

    int counter = 0;
    // bprofiles[model][cat]->GetXaxis()->SetRangeUser(xfirst,xlast);
    //   bprofiles[model][cat]->GetXaxis()->SetTitleOffset( 0.9 );
    //   bprofiles[model][cat]->SetMarkerColor(colors[model]);
  
    for (auto wind : windows){
  
      can[cat]->cd();
      thekey = Form("%s_%s", cat.c_str(), wind.name.c_str()); 
      xfirst = std::min(wind.low,xfirst);
      xlast = std::max(wind.high,xlast);

      // funsEXO17017[thekey] = new TF1( Form("err_correction_%s", thekey.Data()), Form("%s", formEXO17017[thekey].Data()), wind.low, wind.high);
      // funsEXO16027[thekey] = new TF1( Form("err_correction_%s", thekey.Data()), Form("%s", formEXO16027[thekey].Data()), wind.low, wind.high);
      funsEXO17017[thekey] = new TF1( Form("err_correction_%s", thekey.Data()), Form("%s", formEXO17017[thekey].Data()), 500., 1e6);
      funsEXO16027[thekey] = new TF1( Form("err_correction_%s", thekey.Data()), Form("%s", formEXO16027[thekey].Data()), 500., 1e6);
       
      funsEXO16027[thekey]->SetLineColor(2);
      funsEXO17017[thekey]->SetLineColor(4);
       
      funsEXO17017[thekey]->GetHistogram()->GetYaxis()->SetTitle("#beta_{#gamma#gamma} (Events/fb^{-1})");
      funsEXO17017[thekey]->GetHistogram()->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");

      funsEXO16027[thekey]->GetHistogram()->GetYaxis()->SetTitle("#beta_{#gamma#gamma} (Events/fb^{-1})");
      funsEXO16027[thekey]->GetHistogram()->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");

      //Range in Y common
      funsEXO17017[thekey]->GetHistogram()->GetYaxis()->SetRangeUser(1e-6,1.);
      funsEXO16027[thekey]->GetHistogram()->GetYaxis()->SetRangeUser(1e-6,1.);
      
      funsEXO17017[thekey]->GetHistogram()->GetXaxis()->SetRangeUser(0.,5e3);
      funsEXO16027[thekey]->GetHistogram()->GetXaxis()->SetRangeUser(0.,1e4);

      


      // funsEXO17017[thekey]->GetHistogram()->GetYaxis()->SetRangeUser(1e-6,1e-3);
      // funsEXO16027[thekey]->GetHistogram()->GetXaxis()->SetRangeUser(xfirst,xlast);
      if (counter == 0) {funsEXO17017[thekey]->Draw();}
      else {funsEXO17017[thekey]->Draw("same");}
      funsEXO16027[thekey]->Draw("same");
      can[cat]->Update();
      can[cat]->RedrawAxis();
      // can[cat]->Modified();
      can[cat]->Update();
      can[cat]->RedrawAxis();

      ++counter;

    }//end of loop over windows
    
    funleg[cat]->SetHeader(Form("Bias function - %s for 1 fb-1 ",cat.c_str()));
    funleg[cat]->AddEntry( funsEXO16027[thekey] , "EXO-16-027, bias from MC" ,"l" );
    funleg[cat]->AddEntry( funsEXO17017[thekey] , "EXO-17-017, bias from data" ,"l" );
  
    funleg[cat]->Draw("same");
    // can[cat]->RedrawAxis();
    // can[cat]->Modified();
    // can[cat]->Update();

    can[cat]->SaveAs((Form("%s/BiasTermFunctions_%s_%s.png", ws_dir.c_str(), cat.c_str(), year.c_str())) );
    
  }//end of loop over categories
  
 

}
