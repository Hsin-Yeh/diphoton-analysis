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
#include "RooConstVar.h"
#include "RooPolyVar.h"
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
#include "TRandom.h"

using namespace RooFit;
using namespace RooStats;

static const Int_t NCAT = 2; //BB and BE for the moment 
float MINmass = 320.;
float MINmassBE = 360.;
float MAXmass = 1000.;
Int_t nBinsMass= 120;

struct window {
  std::string name; 
  float low;
  float high;
  double norm;
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


//-----------------------------------------------------------------------------------
//Declarations here definition after main
std::vector<thepdf> throwtoys(const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, int ntoys, std::vector<window> windows);
void fitToys(std::vector<thepdf> thetoys, bool blind, bool dobands, const std::string &year, const std::string &ws_dir, std::vector<std::string> cats);
RooAbsPdf* buildPdf(std::string model, std::string name, int catnum, RooRealVar* xvar, RooPolyVar* polymgg, RooWorkspace* w, std::vector<std::string> cats);
void PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, std::string model, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir, int order);
TPaveText* get_labelsqrt( int legendquadrant );
TPaveText* get_labelcms( int legendquadrant, std::string year, bool sim);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir, order, ntoys;
  

  if(argc!=4) {
    std::cout << "Syntax: bkgFit.exe [2016/2017/2018] [input] [ntoys]" << std::endl;
      return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    inputdir = argv[2];
    ntoys = argv[3];
 }

  //Categories
  std::vector<std::string> cats; 
  cats.clear(); 
  cats.push_back("EBEB");
  cats.push_back("EBEE");

  bool blind = true;
  bool dobands = false;

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
  tmpwind.low = 4500.; tmpwind.high = 5500.; tmpwind.name = "4500_5500";
  windows.push_back(tmpwind);

  //========================================================================
  //throw toys 
  std::vector<thepdf> thetoys = throwtoys(year, inputdir, cats, std::stoi(ntoys), windows);
  fitToys(thetoys, blind, dobands, year, inputdir, cats);

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
std::vector<thepdf> throwtoys(const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, int ntoys, std::vector<window> windows){

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    //Models
  std::vector<std::string> models;
  models.push_back("pow");
  models.push_back("expow");
  models.push_back("invpow");
  models.push_back("invpowlin");
  models.push_back("moddijet");
  
  //Order
  std::map<std::string, std::vector<int> > modelorder; //[model][order] for cats
  //These are set from the bkgFit.cc study
  modelorder["pow"].push_back(1);
  modelorder["pow"].push_back(1);

  modelorder["expow"].push_back(2);
  modelorder["expow"].push_back(1);

  modelorder["invpow"].push_back(3);
  modelorder["invpow"].push_back(1);

  modelorder["invpowlin"].push_back(1);
  modelorder["invpowlin"].push_back(1);

  modelorder["moddijet"].push_back(1);
  modelorder["moddijet"].push_back(1);

  std::map<std::string , TFile *> fin; //[model][file]
  std::map<std::string , RooWorkspace* > win; //[model][ws]
  // std::map<std::string , RooAbsPdf* > pdfs; //
  std::vector<thepdf> pdfs;

  std::vector<RooDataSet*> data_toys; 

  RooRealVar* nBackground[NCAT];

  for (auto model : models){

    std::cout << "MODEL " << model << std::endl; 

    fin[model] = TFile::Open(Form("%s/bkg_%s_%s.root", ws_dir.c_str(), model.c_str(), year.c_str()));
    fin[model]->cd();
    win[model] = (RooWorkspace*) fin[model]->Get("HLFactory_ws");

    RooRealVar* mgg = win[model]->var("mgg");  

    for (unsigned int c = 0; c < cats.size(); ++c) {

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
      tmppdf.pdf = (RooGenericPdf*) win[model]->pdf(TString::Format("PhotonsMassBkg_%s_%s",model.c_str(), cats[c].c_str()));
      tmppdf.order = modelorder[model][c];
      tmppdf.ws = win[model];

      nBackground[c] = (RooRealVar*) win[model]->var(TString::Format("PhotonsMassBkg_%s%d_%s_norm",model.c_str(), modelorder[model][c], cats[c].c_str()));

      tmppdf.nominalnorm = nBackground[c]->getVal();

      RooArgSet *set_mgg;

      mgg->setRange("origRange", minMassFit, maxMassFit);
      set_mgg = new RooArgSet(*mgg);
      double renorm = tmppdf.pdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,Range("origRange"))->getVal() / tmppdf.nominalnorm;

      //Time for the number of events predicted by the alternative true underline distribution in windows
      std::vector<window> truenormwinds;
      truenormwinds.clear();
      for (auto wind : windows){
	window tmpwi;
	tmpwi.name = wind.name;
	tmpwi.low = wind.low;
	tmpwi.high = wind.high;

	mgg->setRange(wind.name.c_str(), wind.low, wind.high);
	set_mgg = new RooArgSet(*mgg);

       	tmpwi.norm = (tmppdf.pdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,Range(tmpwi.name.c_str()))->getVal() / renorm  )  ; 
	std::cout<< " renorm " << renorm << " range " << tmpwi.name << " norm orig " << tmppdf.pdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,Range("origRange"))->getVal() << " norm region "<< tmppdf.pdf->createIntegral(RooArgSet(*mgg), RooFit::NormSet(*set_mgg) ,Range(tmpwi.name.c_str()))->getVal() << " true norm in window " << tmpwi.norm << std::endl;

	truenormwinds.push_back( tmpwi );
      }
      tmppdf.truenorms = truenormwinds;
      
      // tmppdf.pdf->Print();
      std::cout << tmppdf.nominalnorm << " " << gRandom->Poisson(  nBackground[c]->getVal() ) << std::endl; 

      for (int toy = 0; toy < ntoys; ++toy){
	// RooDataSet* data = tmppdf.pdf->generate( RooArgSet(*mgg) , TRandom::Poisson(  nBackground[c]->getVal() ) );
	RooDataSet* data = tmppdf.pdf->generate( RooArgSet(*mgg) , gRandom->Poisson(  nBackground[c]->getVal() ) );
	data->SetName(Form("toy_%s_%s_%d",model.c_str(), cats[c].c_str(), toy));
	data->SetTitle(Form("toy_%s_%s_%d",model.c_str(), cats[c].c_str(), toy));
	data_toys.push_back(data);
      }

      tmppdf.toys = data_toys;
      pdfs.push_back(tmppdf);

    }//end of loop over categories

  }//end of loop over models

  return pdfs;

}
//-----------------------------------------------------------------------------------
void fitToys(std::vector<thepdf> thepdfwithtoys, bool blind, bool dobands, const std::string &year, const std::string &ws_dir, std::vector<std::string> cats){

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooPolyVar* polymgg[NCAT];
  RooPolyVar * OneMinusMgg[NCAT];

  //For the building of the models
  //Models
  std::vector<std::string> models;
  models.push_back("pow");
  models.push_back("expow");
  models.push_back("invpow");
  models.push_back("invpowlin");
  models.push_back("moddijet");

  std::map<std::string, TString > card_name; 
  std::map<std::string, HLFactory *> hlf;
  std::map<std::string, RooWorkspace* > w; 

  std::map<std::string, std::vector<RooFitResult*> > fitresult; 
  
  std::map<std::string, std::vector<TCanvas*> > ctmp;

  for (auto model : models){
    card_name[model] = Form( "HighMass-hgg_models_Bkg.%s.rs" , model.c_str() );
    hlf[model] = new HLFactory(Form("HLFactory_%s", model.c_str()), card_name[model], false);
    w[model] = hlf[model]->GetWs();
  }

  for (auto pdf : thepdfwithtoys){

    std::cout << "===============================================================" << std::endl; 
    std::cout << "MODEL " << pdf.model << std::endl; 

    Float_t minMassFit = 0.;
    Float_t maxMassFit = 0.;
    if (pdf.catnum==0){ 
      minMassFit = MINmass;
      maxMassFit = MAXmass;
    } else if (pdf.catnum==1){
      minMassFit = MINmassBE;
      maxMassFit = MAXmass;    
    }

    RooRealVar* mgg = w[pdf.model]->var("mgg");   
    mgg->setUnit("GeV");
    mgg->setMin(minMassFit);
    mgg->setMax(maxMassFit);

    RooArgList *coeffs = new RooArgList();
    
    for (int i=0; i<=pdf.order; ++i){
      coeffs->add( *w[pdf.model]->var(TString::Format("PhotonsMass_bkg_polymgg_a%d_cat%d", i, pdf.catnum))  );
    }

    polymgg[pdf.catnum] = new RooPolyVar(Form("PhotonsMassBkg_polymgg%d_%s", pdf.order, pdf.cat.c_str()), Form("PhotonsMassBkg_polymgg%d_%s", pdf.order, pdf.cat.c_str()), *mgg, *coeffs, 0 );

    //In case of the dijet model we want the polynomial of the 1-mgg
    if (pdf.model == "moddijet"){
      OneMinusMgg[pdf.catnum] = new RooPolyVar(Form("OneMinusMgg_%s", pdf.cat.c_str()), Form("OneMinusMgg_%s", pdf.cat.c_str()), *mgg, RooArgList( RooFit::RooConst(1.), RooFit::RooConst(-1.) ), 0 );
      polymgg[pdf.catnum] = new RooPolyVar(Form("PhotonsMassBkg_polymgg%d_%s", pdf.order, pdf.cat.c_str()), Form("PhotonsMassBkg_polymgg%d_%s", pdf.order, pdf.cat.c_str()), *OneMinusMgg[pdf.catnum], *coeffs, 0 );
    }

    RooAbsPdf* PhotonsMassBkgTmp0 = buildPdf(pdf.model, "", pdf.catnum, mgg, polymgg[pdf.catnum], w[pdf.model], cats);

    for (auto toy : pdf.toys){

      std::cout << "TOY " << toy->GetName() << std::endl;

      fitresult[pdf.model].push_back ( (RooFitResult* ) PhotonsMassBkgTmp0->fitTo(*toy, RooFit::Minimizer("Minuit2"), RooFit::PrintLevel(-1000),RooFit::Warnings(false),SumW2Error(kTRUE), Range(minMassFit,maxMassFit), RooFit::Save(kTRUE)) ) ;

      fitresult[pdf.model].back()->Print("V");

      //************************************************
      // PhotonsMass background fit results per categories 
      // TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
      ctmp[pdf.model].push_back( new TCanvas(Form("ctmp_%s_%s",pdf.model.c_str(), toy->GetName()),"PhotonsMass Background Categories",0,0,500,500)  );
      PlotFitResult(w[pdf.model], ctmp[pdf.model].back(), pdf.catnum, mgg, toy, pdf.model, PhotonsMassBkgTmp0, minMassFit, maxMassFit, blind, dobands, fitresult[pdf.model].back()->floatParsFinal().getSize(), year, ws_dir, pdf.order);
      ctmp[pdf.model].back()->SaveAs( Form("%s/Bkg_%s_%d_%s.png", ws_dir.c_str(),toy->GetName(), pdf.order, year.c_str() ) );
 

    } //end of loop over toys
    
  }//end of loop over pdfs 

}
//-----------------------------------------------------------------------------------
RooAbsPdf* buildPdf(std::string model, std::string name, int catnum, RooRealVar* xvar, RooPolyVar* polymgg, RooWorkspace* w, std::vector<std::string> cats){

  RooAbsPdf* thepdf = new RooGenericPdf(); 
  RooArgList *coefList = new RooArgList();
  std::string formula;
  //------------------------------------------------------------------------
  //dijet model
  if ( model == "dijet" ){
    
    coefList = new RooArgList(*xvar, *w->var(TString::Format("PhotonsMass_bkg_dijet%s_linc_cat%d", name.c_str(), catnum)), *w->var(TString::Format("PhotonsMass_bkg_dijet%s_logc_cat%d", name.c_str(), catnum)));
    formula = "TMath::Max(1e-50,pow(@0,@1+@2*log(@0)))";
    
    //------------------------------------------------------------------------
    //Pow model
  } else if ( model == "pow") {
    // coefList = new RooArgList(*xvar, *w->var(TString::Format("PhotonsMass_bkg_pow_a_cat%d", catnum)));
    coefList = new RooArgList(*polymgg, *w->var(TString::Format("PhotonsMass_bkg_pow_a_cat%d", catnum)));
    formula = "TMath::Max(1e-50,pow(@0,@1))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);
    
    //------------------------------------------------------------------------
    //expow model
  } else if ( model == "expow") {
    //Below the normal coef and formula
    // coefList = new RooArgList(*xvar, *w->var(TString::Format("PhotonsMass_bkg_expow_lam_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)));
    // coefList = new RooArgList(*polymgg, *w->var(TString::Format("PhotonsMass_bkg_expow_lam_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)));
    coefList = new RooArgList(*polymgg, *xvar, *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)));
    // formula = "exp(@3*@0)*pow(@1,@2)";
    formula = "exp(@0)*pow(@1,@2)";

    // w->var(TString::Format("PhotonsMass_bkg_expow_lam_cat%d", catnum))->setConstant(1.);

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);

    //------------------------------------------------------------------------
    //invpow model
  } else if ( model == "invpow") {
    // coefList = new RooArgList(*xvar, *w->var(TString::Format("PhotonsMass_bkg_invpow_slo_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum)) );
    // coefList = new RooArgList(*polymgg, *w->var(TString::Format("PhotonsMass_bkg_invpow_slo_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum)) );
    coefList = new RooArgList(*polymgg, *w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum)) );

    // formula = "pow(1+@0*@1,@2)";
    formula = "pow(1-@0,@1)";
    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

    //------------------------------------------------------------------------
    //invpowlin model
  } else if ( model == "invpowlin") {
    // coefList = new RooArgList(*xvar, *w->var(TString::Format("PhotonsMass_bkg_invpowlin_slo_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_invpowlin_alp_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_invpowlin_bet_cat%d", catnum)) );
    coefList = new RooArgList(*xvar, *w->var(TString::Format("PhotonsMass_bkg_invpowlin_slo_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_invpowlin_alp_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_invpowlin_bet_cat%d", catnum)) );
    // coefList = new RooArgList(*xvar, *polymgg, *w->var(TString::Format("PhotonsMass_bkg_invpowlin_slo_cat%d", catnum)) );

    formula = "pow(1+@0*@1,@2+@3*@0)";
    // formula = "pow(1-@0*@2,@1)";
    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

    //------------------------------------------------------------------------
    //moddijet model
  } else if ( model == "moddijet" ) {
    // coefList = new RooArgList(*xvar, *w->var(TString::Format("PhotonsMass_bkg_moddijet_lina_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_loga_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_linb_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_sqrb_cat%d", catnum)) ) ;
    // coefList = new RooArgList(*polymgg, *w->var(TString::Format("PhotonsMass_bkg_moddijet_lina_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_loga_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_linb_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_sqrb_cat%d", catnum)) ) ;
    coefList = new RooArgList(*xvar, *polymgg, *w->var(TString::Format("PhotonsMass_bkg_moddijet_lina_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_loga_cat%d", catnum)), *w->var(TString::Format("PhotonsMass_bkg_moddijet_linb_cat%d", catnum)) ) ;

    // w->var(TString::Format("PhotonsMass_bkg_moddijet_sqrb_cat%d", catnum))->setConstant(1.); 
    
    // formula = "TMath::Max(1e-50,pow(@0,@1+@2*log(@0))*pow(1.-@0*@4,@3))";
    formula = "TMath::Max(1e-50,pow(@0,@2+@3*log(@0))*pow(@1.,@4))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

  }

  thepdf = new RooGenericPdf(Form("PhotonsMassBkg_%s%s_%s", model.c_str(), name.c_str(), cats[catnum].c_str()), formula.c_str(),  *coefList);
      
  return thepdf;
}

//-----------------------------------------------------------------------------------
void PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, std::string model, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir, int order){

  RooPlot* plotPhotonsMassBkg[NCAT];

  if( blind ) { maxMassFit = 1000.; }
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

