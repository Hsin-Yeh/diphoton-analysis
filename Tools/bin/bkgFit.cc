#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
#include "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h"
#include "diphoton-analysis/Utils/interface/PdfModelBuilder.h"
#include "diphoton-analysis/Utils/interface/CMS_lumi.h"

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
#include "RooDataHist.h"
#include "RooMinimizer.h"
#include "RooStats/RooStatsUtils.h"
#include "RooNumIntConfig.h"
#include "RooNLLVar.h"
#include "TRandom3.h"
#include "RooHist.h"

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
#include "TArrow.h"

using namespace RooFit;
using namespace RooStats;

//-----------------------------------------------------------------------------------
static const Int_t NCAT = 2; //BB and BE for the moment 
// float MINmass, MAXmass, MINmassBE;
// std::map<std::string, float> MINmass, MAXmass;
float MINmass = 500.; //320
float MINmassBE = 500.; //340
//Will go up to 3000 for the moment. 
float MAXmass = 3000.;
float MAXmassBE = 2500.;
//In the full range (0,3000) will have 20 GeV bins same as in the past. 
int nBinsMass = 150;
int nBinsMassBE = 125;

bool gofwithtoys=false;

TRandom3 *RandomGen = new TRandom3();

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
void AddBkgData(RooWorkspace* w, const std::string &isample, const std::string &year, const std::string &ws_dir, bool blind, std::vector<std::string> cats);
theFitResult theFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries);
RooAbsPdf* buildPdf(std::string model, std::string name, int catnum, RooRealVar* xvar, RooWorkspace* w, std::vector<std::string> cats);
gof PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, std::string model, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir, int order);
theFitResult BkgModelFitFunc(RooWorkspace* w, std::string model, bool blind, bool dobands, const std::string &year, const std::string &ws_dir, int order, int c, std::vector<std::string> cats, bool findorder, std::vector<int> orderforcategory, RooArgSet* theparamsout);
void runAllFits(const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, int order);
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample, std::vector<std::string> cats, int order);
std::string getBase(const std::string & sampleName);
void SetConstantParams(const RooArgSet* params) ;
TPaveText* get_labelsqrt( int legendquadrant );
TPaveText* get_labelcms( int legendquadrant, std::string year, bool sim);
gof getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name, bool gofToys);
RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, std::string type, int order, const char* ext);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir, order;
  

  if(argc!=4) {
    std::cout << "Syntax: bkgFit.exe [2016/2017/2018] [input] [order]" << std::endl;
      return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    inputdir = argv[2];
    order = argv[3];
 }

  //========================================================================
  //read the reduced input root trees
  initForFitBkg(inputdir);

  //Categories
  std::vector<std::string> cats; 
  cats.clear(); 
  cats.push_back("EBEB");
  cats.push_back("EBEE");


  //========================================================================
  //Run the fits
  runAllFits(year,"output/bkg",cats,std::stoi(order));

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void runAllFits(const std::string &year, const std::string &ws_dir, std::vector<std::string> cats, int order){

  std::vector<std::string> samples = getSampleListForFit();

  for(auto isample : samples) {
    //We will process one year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    //We will only process data since bkg shape is data-driven
    if ( isample.find("data") == std::string::npos ) continue; 
    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    runfits(year, ws_dir, isample, cats, order);
  }

}
//-----------------------------------------------------------------------------------
void AddBkgData(RooWorkspace* w, const std::string &isample, const std::string &year, const std::string &ws_dir, bool blind, std::vector<std::string> cats) {

  Int_t ncat = NCAT;//BB and BE for the moment

  Float_t minMassFit, maxMassFit;
  //we want to see the whole range. Later we will limit the fit to our ranges of interest. 
  minMassFit = 0.;
  //MAXMASS is the greatest for BB so we are ok. 
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

    if (c==0){ w->import(*dataToFit[c],Rename("Data_EBEB") ); }
    if (c==1){ w->import(*dataToFit[c],Rename("Data_EBEE") ); }

  }

  std::cout << "data, no split" << std::endl;
  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut);
  w->import(*data, Rename("Data"));
    
  data->Print("v");
  std::cout << "---- nX:  " << data->sumEntries() << std::endl; 

  RooPlot* plotPhotonsMassBkgOnlyData[NCAT];
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Only Data");
  TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661,"brNDC");

  //We change here the max value since we want all data but for the plotting we will 
  //be blinded. We should use all data although for the fits. 
  if (blind){
    maxMassFit=1000.;
    nBinsMass=(int) round( (maxMassFit-minMassFit)/20. );
  }

  for (int c=0; c<ncat; ++c) {
    plotPhotonsMassBkgOnlyData[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);
    dataToFit[c]->plotOn(plotPhotonsMassBkgOnlyData[c], MarkerColor(c+2), LineColor(c+2) );  
    if (c==0){plotPhotonsMassBkgOnlyData[c]->Draw();}
    else if (c==1){plotPhotonsMassBkgOnlyData[c]->Draw("same");}
    legdata->AddEntry(plotPhotonsMassBkgOnlyData[c]->getObject(0),Form("Data_%s",cats[c].c_str()),"LPE");

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
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample, std::vector<std::string> cats, int order){

  std::cout << ws_dir << std::endl;

  TString fileBaseName("HighMassGG");    
  TString fileBkgName("HighMassGG.inputbkg");
  TString card_name("HighMass-hgg_models_Bkg.rs");
  // HLFactory hlf("HLFactory", card_name, false);
  // RooWorkspace* w = hlf.GetWs();
  // // range for the variables
  // w->var("mgg")->setMin(MINmass);
  // w->var("mgg")->setMax(MAXmass);
  // w->Print("V");
  std::map<std::string , HLFactory * >  hlf;
  std::map<std::string , RooWorkspace* > w;

  //First we find the order of the functions (true) and then we save those pdfs 
  bool findorder = true; 
  //blind or not
  bool blind = true;

  //========================================================================
  // std::cout << "Adding bkg data" <<std::endl;
  // AddBkgData(w, isample, year, ws_dir, blind, cats);
  //========================================================================
  
  //bands
  bool dobands = false;
  //dijet
  // std::vector<RooFitResult*> fitresults_dijet = BkgModelFitDiJetFunc(w, blind, dobands, year, ws_dir, cats); 
  std::map<std::string, std::map< int, double > > curminnll; //[model][order]
  std::map<std::string, std::map< int, int > > curdof; //[model][order]
  //the prevminll we have to know
  std::map<std::string, double > prevminnll; //[model]
  std::map<std::string, int > prevdof; //[model]
  std::map<std::string, std::map< int, double > > dnll; //[model][order]
  std::map<std::string, std::map< int, int > > diffdof; //[model][order]
  //Fit results
  // std::map< int, std::vector<RooFitResult*> > fitresults;  //[model][order][cats]
  std::map<std::string, std::map< int, theFitResult > > fitresults;  //[model][order]

  //Models
  std::vector<std::string> models;
  //Chosen ones
  // models.push_back("Laurent");
  // models.push_back("PowerLaw");
  // models.push_back("Atlas");
  // models.push_back("Exponential");
  // models.push_back("Expow");
  // models.push_back("invpowlin");

  // models.push_back("pow");
  // models.push_back("Laurent");
  // models.push_back("PowerLaw");
  // models.push_back("Atlas");
  // models.push_back("Exponential");
  // models.push_back("Expow");
  // models.push_back("Chebychev");
  // models.push_back("DijetSimple");
  // models.push_back("Dijet");
  // models.push_back("VVdijet");

  // models.push_back("expow");
  // models.push_back("invpow");
  // models.push_back("invpowlin");
  // models.push_back("moddijet");
  models.push_back("dijet");

  //Latex models
  std::map< std::string, std::string > latexModels;
  latexModels["pow"] =  "$p(m_{\\gamma\\gamma})^a$";
  latexModels["Laurent"] = "$\\sum_{i=0}^{N} p_{i}/m_{\\gamma\\gamma}^{n_{i}}$";
  latexModels["PowerLaw"] = "$\\sum_{i=0}^{N} p_{2i} m_{\\gamma\\gamma}^{p_{2i+1}} $";
  latexModels["Atlas"] = "$(1-m_{\\gamma\\gamma}^{1/3})^{b} m_{\\gamma\\gamma}^{\\sum_{j=0}^{k} a_{j}(logm_{\\gamma\\gamma})^{j}}$";
  latexModels["Exponential"] = "$\\sum_{i=0}^{N} e^{-p_{i} m_{\\gamma\\gamma} }$";
  latexModels["Expow"] = "$e^{p(m_{\\gamma\\gamma})} \\times m_{\\gamma\\gamma}^a$";
  latexModels["Chebychev"] = "$\\sum_{i=0}^{N} p_{i} m_{\\gamma\\gamma}^{i}$";
  latexModels["DijetSimple"] = "$m_{\\gamma\\gamma}^{a+b\\ log(m_{\\gamma\\gamma})}$";
  latexModels["Dijet"] = "$e^{\\sum_{i=0}^{N} p_{i} (log(m_{\\gamma\\gamma}))^{i}}$";
  latexModels["VVdijet"] = "$\\frac{p_{0} ((1-m_{\\gamma\\gamma})/\\sqrt(s))^{p_1}}{(m_{\\gamma\\gamma}/\\sqrt(s))^{p_{2}}}$";
  latexModels["expow"] = "$e^{p(m_{\\gamma\\gamma})} \\times m_{\\gamma\\gamma}^a$";
  latexModels["invpow"] =  "$(1 - p(m_{\\gamma\\gamma}))^a$";
  latexModels["invpowlin"] =  "$(1-m_{\\gamma\\gamma})^{p(m_{\\gamma\\gamma})} $";
  latexModels["moddijet"] =  "$m_{\\gamma\\gamma}^{a+b \\cdot log(m_{\\gamma\\gamma})} \\times p(1-m_{\\gamma\\gamma})^c$";

  //We will have a workspace and file for each final pdf, 
  //final meaning after choosing the correct order
  std::map<std::string, TFile *> fout; //[model][file]
  // TFile * fout; //[model][file]
  std::map<std::string, std::vector<int> > modelorder; //[model][order] for cats
    
  //SET YOUR FINAL ORDER CHOICES HERE (as usual c=0 EBEB c=1 EBEE)
  modelorder["pow"].push_back(4);
  modelorder["pow"].push_back(4);

  modelorder["Laurent"].push_back(1);//2
  modelorder["Laurent"].push_back(1);//2

  modelorder["PowerLaw"].push_back(1);//3
  modelorder["PowerLaw"].push_back(1);//3

  modelorder["Atlas"].push_back(1);//2
  modelorder["Atlas"].push_back(1);//2

  modelorder["Exponential"].push_back(1);//3
  modelorder["Exponential"].push_back(1);//3

  // modelorder["Chebychev"].push_back(3);
  // modelorder["Chebychev"].push_back(3);

  modelorder["DijetSimple"].push_back(2);//2
  modelorder["DijetSimple"].push_back(2);//2

  modelorder["Dijet"].push_back(3);//2
  modelorder["Dijet"].push_back(3);//2

  modelorder["VVdijet"].push_back(1);//1
  modelorder["VVdijet"].push_back(1);//1

  // modelorder["expow"].push_back(2);//1
  // modelorder["expow"].push_back(2);//1

  modelorder["Expow"].push_back(1);//2
  modelorder["Expow"].push_back(1);//2

  modelorder["invpow"].push_back(2);//1
  modelorder["invpow"].push_back(2);//1

  modelorder["invpowlin"].push_back(1);//2
  modelorder["invpowlin"].push_back(1);//2

  modelorder["moddijet"].push_back(2);
  modelorder["moddijet"].push_back(2);

  modelorder["dijet"].push_back(1);
  modelorder["dijet"].push_back(1);

  // fitresults[order] = BkgModelFitFunc(w, models[4], blind, dobands, year, ws_dir, order, cats); 
  
  std::vector<int> dummy;
  FILE *resFile ;
  resFile = fopen(Form("%s/fTestResults_%s.txt",ws_dir.c_str(),year.c_str()),"w");
  // fprintf(resFile,"Family label & Functional form & order & $ -2 \\Delta NLL_{N+1}=-2(NLL_{N+1} - NLL_{N}) $ & gof(\\chi^{2}/ndof) & gof(prob) & p-value & (mimin,hesse,minos) \\\\\n");
  fprintf(resFile,"\\hline\n");
  fprintf(resFile,"\\hline\n");
  fprintf(resFile,"Family label & Functional form & order & \\# params & $ -2 \\Delta NLL_{N+1} $ & gof$(\\chi^{2}/ndof)$ & gof(prob) & p-value & (minim,hesse,minos) \\\\\n");
  fprintf(resFile,"\\hline\n");
  fprintf(resFile,"\\hline\n");

  //We will save all models and orders in a file initially. 
  // fout = new TFile(Form("%s/bkg_%s.root", ws_dir.c_str(), year.c_str()), "recreate");
  // fout->cd();

  if (findorder){
    //Let's loop through all the models
    for (auto themodel : models){

      std::cout << themodel << std::endl; 

      hlf[themodel] = new HLFactory(Form("HLFactory_%s", themodel.c_str() ), card_name, false);
      w[themodel] = hlf[themodel]->GetWs();
      w[themodel]->var("mgg")->setMin(MINmass);
      w[themodel]->var("mgg")->setMax(MAXmass);
      w[themodel]->Print("V");
 
      //========================================================================
      std::cout << "Adding bkg data" <<std::endl;
      AddBkgData(w[themodel], isample, year, ws_dir, blind, cats);
      //========================================================================


      //In the first time, when searching one by one, change below to the 
      //pdf that you are studying. 
      // if (themodel != "dijet"){continue;}

      // if (themodel != "invpow"){continue;}
      // if (themodel != "expow") {continue;}
      // if (themodel != "invpowlin"){continue;}
      // if (themodel != "moddijet"){continue;}

      // //We will save all models and orders in a file initially. 
      fout[themodel] = new TFile(Form("%s/bkg_%s_%s.root", ws_dir.c_str(),themodel.c_str(), year.c_str()), "recreate");
      fout[themodel]->cd();

      RooArgSet* theparamsout = new RooArgSet();

      Int_t ncat = NCAT;
      for (int c = 0; c < ncat; ++c) {
	std::cout << "========================================================"<< std::endl;
	std::cout << "Cat " << cats[c] << std::endl;

	fprintf(resFile,"\\multicolumn{9}{c}{Category %s} \\\\\n",cats[c].c_str());
	fprintf(resFile,"\\hline\n");

	for (int i=1; i<=order; ++i){

	  //We won't go above the order we found. 
	  if (i > modelorder[themodel][c]){continue;}

	  //In power law only odd number of parameters are allowed. 
	  if (themodel == "PowerLaw" && i%2==0 ){continue;}

	  //In exponential only odd number of params allowed
	  if (themodel == "Exponential" && i%2==0 ){continue;}
	  //In exponential needs to be at least of order 3 and smaller than order 5
	  if (themodel == "Exponential" && (i<1 || i > 3) ){continue;}
	  
	  //In VVDijet only one order is allowed
	  if (themodel == "VVDijet" && (i<1 || i > 2) ){continue;}

	  //In DijetSimple only one order is allowed
	  if (themodel == "DijetSimple" && (i!=2) ){continue;}

	  //In Dijet we need at least order 2
	  if (themodel == "Dijet" && (i<2) ){continue;}

	  if (i==1){
	    prevminnll[themodel] = 0.;
	    prevdof[themodel] = 0;
	  } else{
	    prevminnll[themodel] = curminnll[themodel][i-1];
	    prevdof[themodel] = curdof[themodel][i-1];
	  }
	  // int counter = 0;
	  // while (TMath::Prob(  -dnll[themodel][i][c], diffdof[themodel][i][c] ) < 0.05 && counter < 3){
	  // fitresults[themodel][i].push_back( BkgModelFitFunc(w, themodel, blind, dobands, year, ws_dir, i, c, findorder, modelorder[themodel]) ); 
	  fitresults[themodel][i] = BkgModelFitFunc(w[themodel], themodel, blind, dobands, year, ws_dir, i, c, cats, findorder, modelorder[themodel], theparamsout) ; 
	  //   ++counter;
	  // }
	  curminnll[themodel][i] =  fitresults[themodel][i].minNll;
	  curdof[themodel][i]    =  fitresults[themodel][i].fitres->floatParsFinal().getSize();
	  dnll[themodel][i]      = 2.*(curminnll[themodel][i] - prevminnll[themodel]);
	  diffdof[themodel][i]   = curdof[themodel][i] - prevdof[themodel];

	  if ( (-dnll[themodel][i] < 0) && (i>1) ){
	    dnll[themodel][i] = 0.;
	    std::cout << "[WARNING] difference (prevNll-thisNll) < 0 --> chi2=0 --> prob=1" << std::endl;
	  }
	  std::cout << "Order " << i << " -2*dnll " << -dnll[themodel][i] << " dof difference " << diffdof[themodel][i]  << std::endl;
	  std::cout << "TMath::Prob(-2dnll,order-prev_order) " <<  TMath::Prob(  -dnll[themodel][i], diffdof[themodel][i] ) << std::endl;
	  fitresults[themodel][i].fitres->Print("V");

	  //Print the results to file, dijet not in the alternative functions
	  if (themodel != "dijet"){
	    fprintf(resFile,"%s & %s & %d & %d & %10.2f & %10.2f & %10.2f & %10.2f & (%d,%d,%d) \\\\\n", themodel.c_str(), latexModels[themodel].c_str(),i,curdof[themodel][i],-dnll[themodel][i], fitresults[themodel][i].gofresults.chi2overndof, fitresults[themodel][i].gofresults.prob, TMath::Prob(  -dnll[themodel][i], diffdof[themodel][i] ) , fitresults[themodel][i].minimizestatus ,fitresults[themodel][i].hessestatus ,fitresults[themodel][i].minosstatus );
	  }


	  //In the zero order no sense to take the diff
	  // if (i==1){
	  //   diffdof[i][0] = -1000;
	  //   diffdof[i][1] = -1000;
	  // }

	} //end of loop through orders

	  //Get ready for the next category
	  // prevminnll and prevdof will be initialized in the i=1 above;
	for (int i=1; i<=order; ++i){
	  curminnll[themodel][i] = 0.;
	  curdof[themodel][i] = 0;
	  dnll[themodel][i] = 0.;
	  diffdof[themodel][i] = 0;
	}
 

      
      } //end of loop through categories

     
      //Save the resulting workspace with all orders
      //w->Print("V");
      w[themodel]->Write();
      fout[themodel]->Write();
      fout[themodel]->Close();
    } //end of loop through models

    // w->Write();
    // fout->Write();
    // fout->Close();
    
    //Now let's print the results to see what order to choose
    //We want the results per category and then per model. 
    // for (unsigned int c = 0; c < cats.size(); ++c) {
    //   std::cout << "========================================================"<< std::endl;
    //   std::cout << "Cat " << cats[c] << std::endl;

    //   fprintf(resFile,"\\multicolumn{8}{c}{Category %s} \\\\\n",cats[c].c_str());
    //   fprintf(resFile,"\\hline\n");

    //   for (auto themodel : models){
	
    // 	std::cout << themodel << std::endl; 
	
    // 	for (int i=1; i<=order; ++i){

    // 	  //Until we solve the EBEE problem with fits we will use the EBEB max order. 
    // 	  if (i > modelorder[themodel][0]){continue;}

    // 	  // fitresult[c]->floatParsFinal().getSize()
    // 	  std::cout << "Order " << i << " -2*dnll " << -dnll[themodel][i][c] << " dof difference " << diffdof[themodel][i][c]  << std::endl;
    // 	  std::cout << "TMath::Prob(-2dnll,order-prev_order) " <<  TMath::Prob(  -dnll[themodel][i][c], diffdof[themodel][i][c] ) << std::endl;
    // 	  fitresults[themodel][i][c].fitres->Print("V");

    // 	  //Print the results to file, dijet not in the alternative functions
    // 	  if (themodel != "dijet"){
    // 	    fprintf(resFile,"%s & %s & %d & %10.2f & %10.2f & %10.2f & %10.2f & (%d,%d,%d) \\\\\n", themodel.c_str(), latexModels[themodel].c_str(),i,-dnll[themodel][i][c], fitresults[themodel][i][c].gofresults.chi2overndof, fitresults[themodel][i][c].gofresults.prob, TMath::Prob(  -dnll[themodel][i][c], diffdof[themodel][i][c] ) , fitresults[themodel][i][c].minimizestatus ,fitresults[themodel][i][c].hessestatus ,fitresults[themodel][i][c].minosstatus );
    // 	  }
    // 	} //end of loop through orders
    //   } //end of loop through models
    // } //end of loop through categories

  }//end of findorder

  fclose(resFile);
  
}

//-----------------------------------------------------------------------------------
theFitResult theFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){

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
    RooNLLVar *nll=new RooNLLVar("nll","nll",*pdf,*data);
    RooMinimizer *minuit_fitTest = new RooMinimizer(*nll);
    minuit_fitTest->setOffsetting(kTRUE);
    minuit_fitTest->setStrategy(2);
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
  std::cout << "------------------------OFFSET-----------------------------" << std::endl;
  // std::cout << "end of runFit stat=" << stat << " offset=" << offset << " minnll with offset=" << minnll_woffset << " diff= " << minnll<< " nll->getVal() " << nllval << std::endl;
  std::cout << "end of runFit stat=" << stat << " minnll with nll->getVal() " << nllval << std::endl;
  *stat_t = stat;
  *NLL = minnll;

  return tmpfit;
}

//-----------------------------------------------------------------------------------
RooAbsPdf* buildPdf(std::string model, std::string name, int catnum, RooRealVar* xvar, RooWorkspace* w, std::vector<std::string> cats){

  RooAbsPdf* thepdf = new RooGenericPdf(); 
  RooArgList *coefList = new RooArgList();
  RooFormulaVar* asq = new RooFormulaVar(TString::Format("asq_cat%d", catnum),"","-1*@0*@0", *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)) );
  RooFormulaVar* powlinasq = new RooFormulaVar(TString::Format("powlinasq_cat%d", catnum),"","-1*@0*@0", *w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum)));
  std::string formula;
  //------------------------------------------------------------------------
  //dijet model
  if ( model == "dijet" ){
    
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_dijet_linc_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_dijet_logc_cat%d", catnum)));
    formula = "TMath::Max(1e-50,pow(@0,@1+@2*log(@0)))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);

    //------------------------------------------------------------------------
    //Pow model
  } else if ( model == "pow1") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_a_cat%d", catnum))
			      );

    formula = "TMath::Max(1e-50,pow(@0,-1*@1*@1))";
    
    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);
    
    //------------------------------------------------------------------------
    //Pow model
  } else if ( model == "pow2") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam0_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam1_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_a_cat%d", catnum))
			      );

    formula = "TMath::Max(1e-50,pow( (@0*@1) + (@0*@0)*@2, -1*@3*@3))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);
    

    //------------------------------------------------------------------------
    //Pow model
  } else if ( model == "pow3") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam0_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam1_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam2_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_a_cat%d", catnum))
			      );

    formula = "TMath::Max(1e-50,pow( (@0*@1) + (@0*@0)*@2 + (@0*@0*@0)*@3, -1*@4*@4))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);
    
    //------------------------------------------------------------------------
    //Pow model
  } else if ( model == "pow4") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam0_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam1_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam2_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam3_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_a_cat%d", catnum))
			      );

    formula = "TMath::Max(1e-50,pow( (@0*@1) + (@0*@0)*@2 + (@0*@0*@0)*@3 + (@0*@0*@0*@0)*@4, -1*@5*@5))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);
    
    //------------------------------------------------------------------------
    //Pow model
  } else if ( model == "pow5") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam0_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam1_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam2_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam3_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_lam4_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_pow_a_cat%d", catnum))
			      );

    formula = "TMath::Max(1e-50,pow( (@0*@1) + (@0*@0)*@2 + (@0*@0*@0)*@3 + (@0*@0*@0*@0)*@4 + (@0*@0*@0*@0*@0)*@5, -1*@6*@6))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);

    //------------------------------------------------------------------------
    //expow model
  } else if ( model == "expow1") {
    //Below the normal coef and formula
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam0_cat%d", catnum)), 
			      // *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)));
			      *asq);
    // formula = "exp(@1*@0)*pow(@0, -1*@2*@2)";
    formula = "exp(@1*@0)*pow(@0, @2)";

    w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum))->setMin(0.0);  
    if (catnum==0){
      w->var(TString::Format("PhotonsMass_bkg_expow_lam0_cat%d", catnum))->setVal(0.);
      w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum))->setVal(1.5);
    } else if (catnum==1){
      w->var(TString::Format("PhotonsMass_bkg_expow_lam0_cat%d", catnum))->setVal(0.);
      w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum))->setVal(3.5);
    }

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);

    //------------------------------------------------------------------------
    //expow2 model
  } else if ( model == "expow2") {
    // RooArgList *bla = new RooArgList(*w->var(TString::Format("PhotonsMass_bkg_expow_lam0_cat%d", catnum)), 
    // 				     *w->var(TString::Format("PhotonsMass_bkg_expow_lam1_cat%d", catnum)));
    
    // RooFormulaVar *hmax = new RooFormulaVar(Form("hmax_cat%d", catnum),"( @1 != 0. ? 0. : @0*3500.)", *bla );
    // RooFormulaVar *hmax = new RooFormulaVar(Form("hmax_cat%d", catnum),"@0*3500.", *bla );
    // RooFormulaVar *hmax = new RooFormulaVar(Form("hmax_cat%d", catnum),"( @1 != 0. ? (-@0/(4.*@1)>340. && -@0/(4.*@1)<3500. ? @0*@0/(4.*@1+@1) : TMath::Max(@0*3500+2*@1*3500.*3500,@0*3500+2*@1*300.*300)) : @0*3500.)", *bla );
    // RooFormulaVar *hmax = new RooFormulaVar(Form("hmax_cat%d", catnum),"(@1 != 0. ? TMath::Max(@0*3500+2*@1*3500.*3500,@0*3500+2*@1*300.*300) : @0*3500.)", *bla );
    // RooFormulaVar *hmax = new RooFormulaVar(Form("hmax_cat%d", catnum),"( @1 != 0. ? 0.*@1 : @0*3500.)", *bla );

    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam0_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam1_cat%d", catnum)), 
			      // *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)) );
			      *asq);//, 
			      //*hmax);
    // formula = "exp( @1*@0 + @2*@0*@0)*pow( @0, -1*@3*@3 )";
    // formula = "exp( @1*@0 + @2*@0*@0)*pow( @0, @3 + @4)";
    formula = "exp( @1*@0 + @2*@0*@0)*pow( @0, @3)";

    w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum))->setMin(0.0);    

    //------------------------------------------------------------------------
    //expow3 model
  } else if ( model == "expow3") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam0_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam1_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam2_cat%d", catnum)), 
			      // *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)) );
			      *asq);
    // formula = "exp( @1*@0 + @2*@0*@0 + @3*@0*@0*@0 )*pow( @0, -1*@4*@4 )";
    formula = "exp( @1*@0 + @2*@0*@0 + @3*@0*@0*@0 )*pow( @0, @4 )";

    w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum))->setMin(0.0);    

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);

    //------------------------------------------------------------------------
    //expow4 model
  } else if ( model == "expow4") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam0_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam1_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam2_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_expow_lam3_cat%d", catnum)), 
			      *asq);
			      // *w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum)) );
    // formula = "exp( @1*@0 + @2*@0*@0 + @3*@0*@0*@0 + @4*@0*@0*@0*@0 )*pow( @0, -1*@5*@5)";
    formula = "exp( @1*@0 + @2*@0*@0 + @3*@0*@0*@0 + @4*@0*@0*@0*@0 )*pow( @0, @5 )";

    w->var(TString::Format("PhotonsMass_bkg_expow_alp_cat%d", catnum))->setMin(0.0);    

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;
    // exit(-1);
    //------------------------------------------------------------------------

    //invpow model
  } else if ( model == "invpow") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_invpow_slo_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpow_bet_cat%d", catnum)) );
    // formula = "pow(1+@0*@1,@2)";
    // formula = "pow(1-@0,@1)";
    formula = "pow(1+@0*@1,@2+@3*@0)";
    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

    //------------------------------------------------------------------------

    //invpow1 model
  } else if ( model == "invpow1") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_invpow_slo_cat%d", catnum)), 
			      // *w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum)) 
			      *powlinasq
			      );
    formula = "pow(1+@0*@1,@2)";

    w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setMin(0.0);   
    if (catnum==0){
      w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setVal(2.5);
    } else if (catnum==1){
      w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setVal(0.5);
    }

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

    //------------------------------------------------------------------------

    //invpow2 model
  } else if ( model == "invpow2") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_invpow_slo_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpow_qua_cat%d", catnum)), 
			      // *w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum)) 
			      *powlinasq
			      );

    formula = "pow(1+@0*@1+@2*@0*@0,@3)";

    w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setMin(0.0);    
    if (catnum==0){
      w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setVal(2.5);
    } else if (catnum==1){
      w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setVal(0.5);
    }

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

    //------------------------------------------------------------------------
    //invpowlin model
  } else if ( model == "invpowlin1") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_invpowlin_slo_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpowlin_alp_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpowlin_bet_cat%d", catnum)) 
			      );

    formula = "pow(1+@0*@1,@2+@3*@0)";
    // formula = "pow(1-@0*@2,@1)";
    
    if (catnum==0){
      w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setVal(2.5);
    } else if (catnum==1){
      w->var(TString::Format("PhotonsMass_bkg_invpow_alp_cat%d", catnum))->setVal(0.5);
    }

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

    //------------------------------------------------------------------------
    //invpowlin model
  } else if ( model == "invpowlin2") {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_invpowlin_slo_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpowlin_alp_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpowlin_bet_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_invpowlin_qua_cat%d", catnum)) 
			      );

    formula = "pow(1+@0*@1,@2+@3*@0+@4*@0*@0)";
    // formula = "pow(1-@0*@2,@1)";
    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

    //------------------------------------------------------------------------
    //moddijet model
  } else if ( model == "moddijet1" ) {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_lina_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_loga_cat%d", catnum)), 
			      // *linbsq,
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_linb_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_sqrb_cat%d", catnum)) 
			      ) ;

    w->var(TString::Format("PhotonsMass_bkg_moddijet_sqrb_cat%d", catnum))->setConstant(1); 
    
    formula = "TMath::Max(1e-50,pow(@0,@1+@2*log(@0))*pow(1.-@0*@4,@3))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

  
    //------------------------------------------------------------------------
    //moddijet model
  } else if ( model == "moddijet2" ) {
    coefList = new RooArgList(*xvar, 
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_lina_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_loga_cat%d", catnum)), 
			      // *linbsq,
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_linb_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_sqrb_cat%d", catnum)), 
			      *w->var(TString::Format("PhotonsMass_bkg_moddijet_quab_cat%d", catnum)) 
			      ) ;

    w->var(TString::Format("PhotonsMass_bkg_moddijet_sqrb_cat%d", catnum))->setConstant(1); 
    
    formula = "TMath::Max(1e-50,pow(@0,@1+@2*log(@0))*pow((1.-@0*@4)+(@5*(1.-@0*@4)^2),@3))";

    std::cout << "====================================================================" << std::endl;
    std::cout << formula << std::endl;

  }


  thepdf = new RooGenericPdf(Form("PhotonsMassBkg_%s%s_%s", model.c_str(), name.c_str(), cats[catnum].c_str()), Form("PhotonsMassBkg_%s%s_%s", model.c_str(), name.c_str(), cats[catnum].c_str()), formula.c_str(),  *coefList);
     
  std::cout << thepdf->GetName() << std::endl;

 
  return thepdf;
}

//-----------------------------------------------------------------------------------
theFitResult BkgModelFitFunc(RooWorkspace* w, std::string model, bool blind, bool dobands, const std::string &year, const std::string &ws_dir, int order, int c, std::vector<std::string> cats, bool findorder, std::vector<int> orderforcategory, RooArgSet* theparamsout) {

  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  RooDataSet* data;
  RooRealVar* nBackground; 
  theFitResult fitresult;

  std::map< std::string, RooAbsPdf* > PhotonsMassBkgTmp0; //[model][cat]
  std::map< std::string, RooArgSet  > thepreParams;
  
  RooRealVar* mgg = w->var("mgg");  
  mgg->setUnit("GeV");

  // mgg->setConstant(false);
  Float_t minMassFit = 0.;
  Float_t maxMassFit = 0.;

  if (c==0){ 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  } else if (c==1){
    minMassFit = MINmassBE;
    maxMassFit = MAXmassBE;    
  }
  mgg->setRange(minMassFit,maxMassFit);
  //Will use roughly the same bins as in the past (roughly 20 GeV here) although for the non reso 
  //analysis they use 25 GeV.  
  nBinsMass = (int) round( (maxMassFit-minMassFit)/20. );
  mgg->setBins(nBinsMass);
  
  //In case we do not want to find the correct order then 
  //we set the order we want for the relevant category
  if ( !findorder ) { order = orderforcategory[c];}

  data = (RooDataSet*) w->data(Form("Data_%s", cats[c].c_str()));
  
  std::string prevmodel = " ";
  //Will use the ones from Pasquale for expow
  if ((model=="pow") && (order==1)){model="pow1";prevmodel="pow";}
  if ((model=="pow") && (order==2)){model="pow2";prevmodel="pow";}
  if ((model=="pow") && (order==3)){model="pow3";prevmodel="pow2";}
  if ((model=="pow") && (order==4)){model="pow4";prevmodel="pow3";}
  if ((model=="pow") && (order==5)){model="pow5";prevmodel="pow4";}
  if ((model=="expow") && (order==1)){model="expow1";}
  if ((model=="expow") && (order==2)){model="expow2";prevmodel="expow";}
  if ((model=="expow") && (order==3)){model="expow3";prevmodel="expow2";}
  if ((model=="expow") && (order==4)){model="expow4";prevmodel="expow3";}
  if ((model=="invpow") && (order==1)){model="invpow1";}
  if ((model=="invpow") && (order==2)){model="invpow2";}
  if ((model=="invpowlin") && (order==1)){model="invpowlin1";}
  if ((model=="invpowlin") && (order==2)){model="invpowlin2";}
  if ((model=="moddijet") && (order==1)){model="moddijet1";}
  if ((model=="moddijet") && (order==2)){model="moddijet2";}

  PdfModelBuilder pdfsModel;
  pdfsModel.setObsVar(mgg);
  pdfsModel.setWorkspace(w);

  if (model!="Laurent" && model!="PowerLaw" && model!="Atlas" && model!="Exponential" && 
      model!="Chebychev" && model!="DijetSimple" && model!="Dijet" && model!="VVdijet" &&
      model!="Expow"
      ){
    PhotonsMassBkgTmp0[model] =  buildPdf(model, "", c, mgg, w, cats)  ;
  } else {
    PhotonsMassBkgTmp0[model] = getPdf(pdfsModel,model,order,Form("ftest_pdf_%s",cats[c].c_str() ));
  } 
  
  
  // Make list of model parameters
  // RooArgSet* theparams = PhotonsMassBkgTmp0->getParameters(*mgg);

  // RooArgSet* theparams = PhotonsMassBkgTmp0[model][c]->getParameters((const RooArgSet*)(0));
  RooArgSet* theparams = PhotonsMassBkgTmp0[model]->getParameters(*mgg);
  if (c==0){ 
    theparamsout->add(*theparams);
  } else if (c==1){
    theparamsout->add(*theparams);
    w->defineSet("parameters", *theparamsout);
    w->defineSet("observables", *mgg);
  }
 
  
  if (order==1){ w->saveSnapshot(Form("initial_%s_fit_params_cat%s", model.c_str(), cats[c].c_str() ), *theparams, kTRUE); }
  // else {w->loadSnapshot(Form("initial_%s_fit_params_cat%s", model.c_str(), cats[c].c_str() ) );}

  // w->saveSnapshot(Form("fit_params_cat%s", cats[c].c_str() ), *params, kTRUE);
  // w->loadSnapshot(Form("fit_params_cat%s", cats[c].c_str() ) );
  // Save snapshot of current parameters because the gof below will change them. 
  // RooArgSet curparams = (RooArgSet) theparams->snapshot() ;
  // RooArgSet preParams;
  // theparams->snapshot(preParams) ;

  nBackground = new RooRealVar(Form("PhotonsMassBkg_%s_%s_norm", model.c_str(), cats[c].c_str()), "nbkg", data->sumEntries(),0,3*data->sumEntries());

  // RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Pow_cat%d",c), "TMath::Max(1e-50,1.+ @0*@1 + pow(@0,2.)*@2 + pow(@0,3.)*@3)", RooArgList(*mgg, *w->var(TString::Format("PhotonsMass_bkg_pow_a0_cat%d",c)), *w->var(TString::Format("PhotonsMass_bkg_pow_a1_cat%d",c)), *w->var(TString::Format("PhotonsMass_bkg_pow_a2_cat%d",c)),*w->var(TString::Format("PhotonsMass_bkg_pow_a3_cat%d",c)) ) );
  
  //From https://roostatsworkbook.readthedocs.io/en/latest/modelbuilding.html :
  //For reasons of backward compatibility, likelihood offsetting is not applied by 
  //default in RooFit, It can be activated by adding Offset(kTRUE) to RooAbsReal::createNLL() 
  //or RooAbsReal::fitTo(), and users are recommended to do this for all non-trivial fits.
    
  //RooFit::Minos(kTRUE), RooFit::Strategy(2), RooFit::Offset(true),
  // ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-4); 
  // ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-8); 
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000); 
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");   
  //https://root.cern.ch/root/html/tutorials/roofit/rf901_numintconfig.C.html
  // Print current global default configuration for numeric integration strategies
  // RooAbsReal::defaultIntegratorConfig()->Print("v") ;
  // RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8);
  // RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8);
 
  //This is the current command but will change that
  // fitresult.push_back( (RooFitResult* ) PhotonsMassBkgTmp0->fitTo(*data[c], RooFit::Minimizer("Minuit2"), RooFit::Offset(kTRUE), RooFit::Strategy(2), RooFit::PrintLevel(3), RooFit::Warnings(false), SumW2Error(kTRUE), Range(minMassFit,maxMassFit), RooFit::Save(kTRUE)) );
  //to this one
  int fitStatus = 0;
  double thisNll = 0.;
  fitresult = theFit(PhotonsMassBkgTmp0[model], data, &thisNll, &fitStatus, /*max iterations*/ 3) ;
    
  //After the fit for a specific model the results for the next order will be better with 
  //initial values the previous order ones but when changing models this is an issue. So, 
  //we set them independently for each model
 
  // theparams->snapshot(thepreParams[model]) ;

  w->saveSnapshot(Form("%s_fit_params_cat%s", model.c_str(), cats[c].c_str() ), *theparams, kTRUE);
 
  if (order == orderforcategory[c]){
    w->import(*PhotonsMassBkgTmp0[model], RecycleConflictNodes());
    w->import(*nBackground, RecycleConflictNodes());
    theparams->printLatex();
  }

  std::cout << TString::Format("******************************** Background Fit results %s cat %s***********************************", model.c_str(), cats[c].c_str()) << std::endl;
  fitresult.fitres->Print("V");

  //************************************************
  // Plot PhotonsMass background fit results per categories 
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",550,600);
  gof gofresult = PlotFitResult(w, ctmp, c, mgg, data, model, PhotonsMassBkgTmp0[model], minMassFit, maxMassFit, blind, dobands, fitresult.fitres->floatParsFinal().getSize(), year, ws_dir, order);
  fitresult.gofresults = gofresult;
  ctmp->SetLogy();
  // ctmp->SaveAs( Form("%s/Bkg_cat%d_%s_%s_log.png", ws_dir.c_str(),c, model.c_str(), year.c_str() ) );
  ctmp->SaveAs( Form("%s/Bkg_cat%d_%s%d_%s.png", ws_dir.c_str(),c, model.c_str(), order, year.c_str() ) );
  ctmp->SetLogy(0);
  ctmp->SaveAs( Form("%s/Bkg_cat%d_%s%d_%s.png", ws_dir.c_str(),c, model.c_str(), order, year.c_str() ) );
  // ctmp->SaveAs( Form("%s/Bkg_cat%d_%s_%s.png", ws_dir.c_str(),c, model.c_str(), year.c_str() ) );

  //Due to the toys fit it changes the fit params. So, I have to load the snapshot. 
  // if (order>=2){
  //   w->loadSnapshot(Form("%s_fit_params_cat%s", prevmodel.c_str(), cats[c].c_str() ) );
  // }
    
  // theparams->assignValueOnly(preParams[model]);
  // *theparams = *curparams;
    
  // delete PhotonsMassBkgTmp0;

  return fitresult;

}

//-----------------------------------------------------------------------------------
gof PlotFitResult(RooWorkspace* w, TCanvas* ctmp, int c, RooRealVar* mgg, RooDataSet* data, std::string model, RooAbsPdf* PhotonsMassBkgTmp0, float minMassFit, float maxMassFit, bool blind, bool dobands, int numoffittedparams, const std::string &year, const std::string &ws_dir, int order){
  
  RooPlot* plotPhotonsMassBkg[NCAT];
  RooPlot* resid[NCAT]; 
  
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

  resid[c] = mgg->frame(minMassFit, maxMassFit,nBinsMass);

  TLine *lineone = new TLine(minMassFit,0,maxMassFit,0);
  lineone->SetLineColor(kBlue);
  lineone->SetLineWidth(2);
  resid[c]->addObject(lineone);

  RooCurve *fitc   = (RooCurve *) plotPhotonsMassBkg[c]->getObject(int(plotPhotonsMassBkg[c]->numItems()-2));
  RooHist* hist   = (RooHist*) plotPhotonsMassBkg[c]->getObject(int(plotPhotonsMassBkg[c]->numItems()-1));
  std::cout << hist->GetName() << " " << fitc->GetName() << std::endl; 
  //According to https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars 
  //we should use asymmetric vertical bars with "correct coverage". Since this is 
  //RooFit and hist->SetBinErrorOption(TH1::kPoisson) doesn't work we manually loop. 
  double alpha = 1 - 0.6827;
  for (int i = 0; i < hist->GetN(); ++i) {
    int N = 0;
    // if( blind && 
    // 	(hist->GetX()[i]-hist->GetErrorXlow(i)) >= minMassFit && 
    // 	(hist->GetX()[i]+hist->GetErrorXhigh(i)) <= maxMassFit
    // 	){ 
    //   hist->SetPoint(i,hist->GetX()[i],0.);
    //   N = fitc->Eval(hist->GetX()[i]);
    // } else {
    //   N = hist->GetY()[i];
    // }
    N = hist->GetY()[i];
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
    hist->SetPointEYlow(i, N-L);
    hist->SetPointEYhigh(i, U-N);
  }

  RooHist* hresid = plotPhotonsMassBkg[c]->residHist(hist->GetName(),fitc->GetName(),true);
  if (!dobands) {resid[c]->addPlotable(hresid,"PE");}

  ctmp->Divide(1,2);

  ctmp->cd(1);
  gPad->SetPad(0.,0.38,1.,0.95);
  gPad->SetTopMargin(0.0015);
  gPad->SetBottomMargin(0.02);
  gPad->SetLogy(true);
  gPad->SetFillStyle(0);
  gPad->SetTickx();
        
  ctmp->cd(2);        
  gPad->SetPad(0.,0.,1.,0.38);
  gPad->SetFillStyle(0);
  gPad->SetTopMargin(0.0015);
  gPad->SetBottomMargin(0.32);
  gPad->SetFillStyle(0);
  gPad->SetTickx();

  ctmp->cd(1);

  plotPhotonsMassBkg[c]->SetTitle("");
  plotPhotonsMassBkg[c]->GetYaxis()->SetTitleSize(0.09);
  plotPhotonsMassBkg[c]->GetXaxis()->SetLabelSize( 1.2*plotPhotonsMassBkg[c]->GetXaxis()->GetLabelSize() );
  plotPhotonsMassBkg[c]->GetXaxis()->SetTitleSize( 1.2*plotPhotonsMassBkg[c]->GetXaxis()->GetTitleSize() );
  plotPhotonsMassBkg[c]->GetXaxis()->SetTitleOffset( 1.15 );
  // double ymax = fitc->interpolate(plotPhotonsMassBkg[c]->GetXaxis()->GetXmin())*2.;
  // double ymin = fitc->interpolate(plotPhotonsMassBkg[c]->GetXaxis()->GetXmax())*0.3;
  // ymin = std::max(1.1e-1,ymin);
  // plotPhotonsMassBkg[c]->GetYaxis()->SetRangeUser(ymin,ymax);
  plotPhotonsMassBkg[c]->GetXaxis()->SetMoreLogLabels();
  plotPhotonsMassBkg[c]->GetYaxis()->SetLabelSize( plotPhotonsMassBkg[c]->GetXaxis()->GetLabelSize() * ctmp->GetWh() / gPad->GetWh() * 1.3 );
  plotPhotonsMassBkg[c]->GetYaxis()->SetTitleSize( plotPhotonsMassBkg[c]->GetXaxis()->GetTitleSize() * ctmp->GetWh() / gPad->GetWh() * 1.3 );
  plotPhotonsMassBkg[c]->GetYaxis()->SetTitleOffset( 0.75 );

  plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  plotPhotonsMassBkg[c]->SetAxisRange(0.11,plotPhotonsMassBkg[c]->GetMaximum()*2.0,"Y");

  if (!dobands){ plotPhotonsMassBkg[c]->Draw(); } 

  TLegend *legdata = new TLegend(0.2,0.37,0.35,0.62, "", "brNDC");
  legdata->SetFillColor(kWhite);
  legdata->SetFillStyle(0);
  legdata->SetShadowColor(kWhite);
  legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
  if (order == 0 || model == "dijet"){
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),Form("Fit Model: %s", model.c_str()),"L");
  } else{
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),Form("Fit Model: %s Order %d", model.c_str(), order),"L");
  }
  // if (dobands){
  //   legdata->AddEntry(onesigma,"#pm 1 s.d.","f");
  //   legdata->AddEntry(twosigma,"#pm 2 s.d.","f");
  // }

  legdata->SetTextSize(0.05);
  legdata->SetTextFont(42);
  // legdata->SetTextAlign(31);
  legdata->SetBorderSize(0);
  if (!dobands) {legdata->Draw("same");}

  TPaveText* label_cms = get_labelcms(1, year, false);
  TPaveText* label_sqrt = get_labelsqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  //write down the chi2 of the fit on the
 
  TPaveText* label_chi2 = new TPaveText(0.55,0.33,0.79,0.44, "brNDC");
  label_chi2->SetFillColor(kWhite);
  label_chi2->SetTextSize(0.05);
  label_chi2->SetTextFont(42);
  label_chi2->SetTextAlign(31); // align right
  label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", gofresult.chi2overndof));
  if (gofwithtoys){
    label_chi2->AddText(TString::Format("Chi square Prob (toys) = %.3f", gofresult.prob));
  } else{
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", gofresult.prob));
  }
  if (!dobands) {label_chi2->Draw("same");}


  ctmp->cd(2);
  gPad->SetGridy();
  resid[c]->SetTitle("");
  resid[c]->GetXaxis()->SetMoreLogLabels();
  resid[c]->GetXaxis()->SetNdivisions(510);
  resid[c]->GetYaxis()->SetNdivisions(505);
  resid[c]->GetYaxis()->CenterTitle();
  resid[c]->GetYaxis()->SetTitleSize( plotPhotonsMassBkg[c]->GetYaxis()->GetTitleSize() * 1.4 );
  resid[c]->GetYaxis()->SetTitleOffset( plotPhotonsMassBkg[c]->GetYaxis()->GetTitleOffset() * 0.6 ); 
  resid[c]->GetYaxis()->SetLabelSize( plotPhotonsMassBkg[c]->GetYaxis()->GetLabelSize() * 1.3 );
  resid[c]->GetXaxis()->SetTitleSize( plotPhotonsMassBkg[c]->GetXaxis()->GetTitleSize() * 2. * 1.3/1.2 );
  resid[c]->GetXaxis()->SetTitleOffset( plotPhotonsMassBkg[c]->GetXaxis()->GetTitleOffset() );
  resid[c]->GetXaxis()->SetLabelSize( plotPhotonsMassBkg[c]->GetXaxis()->GetLabelSize() * 6.5/3.5 );
  resid[c]->GetYaxis()->SetRangeUser( -3, 3 );
  resid[c]->GetYaxis()->SetTitle("(data-fit)/#sigma_{stat}");
  resid[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  if (!dobands) {resid[c]->Draw();}

  ctmp->cd(1);
  float margin = gPad->GetBottomMargin()+gPad->GetTopMargin();
  gPad->SetTopMargin(0.1*margin);
  gPad->SetBottomMargin(0.1*margin);

  ctmp->cd(2);
  margin = gPad->GetBottomMargin()+gPad->GetTopMargin();
  gPad->SetBottomMargin(margin);
  gPad->SetTopMargin(0.1*margin);
                
  plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("");
  plotPhotonsMassBkg[c]->GetXaxis()->SetLabelSize(0.);
  
  ctmp->cd();
  if (c==0){
    TLatex *ptCat0= new TLatex(0.2,0.44175, "EBEB");
    ptCat0->SetTextFont(61);
    ptCat0->SetTextSize(0.051);
    ptCat0->Draw("same");
  } else if (c==1){
    TLatex *ptCat1= new TLatex(0.2,0.44175, "EBEE");
    ptCat1->SetTextFont(61);
    ptCat1->SetTextSize(0.051);
    ptCat1->Draw("same");
  }
   
  CMS_lumi( (TPad*) gPad, 4, 0 );

  //********************************************************************************//
  if (dobands) {

    ctmp->cd(1);
    RooAbsPdf *cpdf = PhotonsMassBkgTmp0;
    TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *bias = new TGraphAsymmErrors();

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
      minim.setMinimizerType("Minuit2");
      minim.setStrategy(0);
      minim.setPrintLevel(-1);
      // double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
      // double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
      minim.migrad();

      minim.hesse();
      RooFitResult* result = minim.lastMinuitFit();
      double errm = nlim->getPropagatedError(*result);
      onesigma->SetPoint(i-1,center,nombkg);
      onesigma->SetPointError(i-1,lowedge,upedge,-errm, errm);

      // minim.minos(*nlim);
      // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
      // onesigma->SetPoint(i-1,center,nombkg);
      // onesigma->SetPointError(i-1,lowedge,upedge,-nlim->getErrorLo(),nlim->getErrorHi());
	
      // eventually if cl = 0.95 this is the usual 1.92!      
      // minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
      //if on the other hand cl = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0)= 0.95449974 this is 2. 
      minim.setErrorLevel(2); 
      minim.migrad();

      minim.hesse();
      result = minim.lastMinuitFit();
      errm = nlim->getPropagatedError(*result);
      twosigma->SetPoint(i-1,center,nombkg);
      // twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
      twosigma->SetPointError(i-1,lowedge,upedge,-nlim->getErrorLo(),nlim->getErrorHi());

      // minim.minos(*nlim);
      // printf("lowedge, = %5f, upper edge = %5f, errlo = %5f, errhi = %5f\n",lowedge,upedge,nlim->getErrorLo(),nlim->getErrorHi());
      // twosigma->SetPoint(i-1,center,nombkg);
      // twosigma->SetPointError(i-1,lowedge,upedge,-nlim->getErrorLo(),nlim->getErrorHi());
	
      delete nll;
      delete epdf;
    }

    mgg->setRange("errRange",minMassFit,maxMassFit);

    // plotPhotonsMassBkg[c]->Draw(); 
    plotPhotonsMassBkg[c]->remove(hist->GetName(), false);
    plotPhotonsMassBkg[c]->remove(fitc->GetName(), false);

    twosigma->SetLineColor(kYellow);
    twosigma->SetFillColor(kYellow);
    twosigma->SetMarkerColor(kYellow);
    // twosigma->Draw("E2 SAME");
    // twosigma->Draw("E2");

    onesigma->SetLineColor(kGreen);
    onesigma->SetFillColor(kGreen);
    onesigma->SetMarkerColor(kGreen);
    // onesigma->Draw("E2 SAME");
      
    bias->SetLineColor(kOrange);
    bias->SetFillColor(kOrange);
    bias->SetMarkerColor(kOrange);
    // bias->Draw("E2 SAME");

    plotPhotonsMassBkg[c]->addObject(twosigma,"LE3");
    plotPhotonsMassBkg[c]->addObject(onesigma,"LE3");
    plotPhotonsMassBkg[c]->addObject(fitc);
    plotPhotonsMassBkg[c]->addObject(hist);
   
    plotPhotonsMassBkg[c]->Draw(); 


    legdata->AddEntry(onesigma,"#pm 1 s.d.","F");
    legdata->AddEntry(twosigma,"#pm 2 s.d.","F");
    // legdata->AddEntry(bias,"bias","F");
    legdata->Draw("same");

    label_chi2->Draw("same");

    ctmp->cd(2);
    //error in  pulls
    TGraphAsymmErrors *ronesigma = (TGraphAsymmErrors *) onesigma->Clone();
    TGraphAsymmErrors *rtwosigma = (TGraphAsymmErrors *) twosigma->Clone();

    for (int i = 0; i < ronesigma->GetN(); ++i) {
      double px = ronesigma->GetX()[i];
      double py = ronesigma->GetY()[i];
      ronesigma->SetPoint(i,px,0.);
      rtwosigma->SetPoint(i,px,0.);
      //double hx = hist->GetX()[i];
      double hy = hist->GetY()[i];

      double oerrp = ronesigma->GetErrorYhigh(i);
      double oerrm = ronesigma->GetErrorYlow(i);
      double terrp = rtwosigma->GetErrorYhigh(i);
      double terrm = rtwosigma->GetErrorYlow(i);
      double herrp = hist->GetErrorYhigh(i);
      double herrm = hist->GetErrorYlow(i);

      if (py < hy){
	if (herrm == 0.) { continue; }
	oerrp /= herrm;
	terrp /= herrm;
	oerrm /= herrm;
	terrm /= herrm;
      } else{
	if (herrp == 0.) { continue; }
	oerrp /= herrp;
	terrp /= herrp;
	oerrm /= herrp;
	terrm /= herrp;
      }

      ronesigma->SetPointEYhigh(i,oerrp);
      ronesigma->SetPointEYlow(i,oerrm);
      rtwosigma->SetPointEYhigh(i,terrp);
      rtwosigma->SetPointEYlow(i,terrm);
      
    } //end of loop over sigma points

    resid[c]->addObject(rtwosigma,"E2");
    resid[c]->addObject(ronesigma,"E2");
    resid[c]->addPlotable(hresid,"PE");
    // resid[c]->addPlotable(rtwosigma,"E2");
    // resid[c]->addPlotable(ronesigma,"E2");

    resid[c]->Draw();

    // onesigma->Draw("same");
    
    
  }

  return gofresult;
  
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
RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, std::string type, int order, const char* ext){
  if (type=="DijetSimple") return pdfsModel.getDijetSimple(Form("%s_dijetsimple%d",ext,order),order); 
  else if (type=="Dijet") return pdfsModel.getDijet(Form("%s_dijet%d",ext,order),order); 
  else if (type=="VVdijet") return pdfsModel.getVVdijet(Form("%s_vvdijet%d",ext,order),order); 
  else if (type=="Atlas") return pdfsModel.getAtlas(Form("%s_atlas%d",ext,order),order); 
  else if (type=="Expow") return pdfsModel.getExpow(Form("%s_expow%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else {
    std::cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << std::endl;
    return NULL;
  }
}
