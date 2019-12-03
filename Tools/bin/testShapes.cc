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
  TFile * freco;
  TFile * fparamshape;
  RooWorkspace * wsreco;
  RooWorkspace * wsparamshape;

};

//-----------------------------------------------------------------------------------
static const Int_t NCAT = 2; //BB and BE for the moment 

//-----------------------------------------------------------------------------------
//Declarations here definition after main
void testTheShapes(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &couplingIn);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  std::string year, inputdir, thecoupling;

  if(argc!=4) {
    std::cout << "Syntax: testShapes.exe [2016/2017/2018] [input] [kMplxxx]" << std::endl;
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
  testTheShapes(year,"input/workspaces", "output/FinalParametricShape/workspaces", thecoupling);

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void testTheShapes(const std::string &year, const std::string &ws_indir, const std::string &ws_outdir, const std::string &couplingIn){

  std::vector<std::string> samples = getSampleListForFit();

  //========================================================================
  std::cout << "Check the final parametric shape" <<std::endl;
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
    tmpin.freco = TFile::Open(Form("%s/ws_ResponseAndGen_M%s_k%s_%s.root", ws_indir.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()));  
    tmpin.wsreco = (RooWorkspace*) tmpin.freco->Get("HLFactory_ws");
    tmpin.fparamshape = TFile::Open(Form("%s/SignalParametricShapes_ws_kMpl%s.root", ws_outdir.c_str(), coupling.c_str()));
    tmpin.wsparamshape = (RooWorkspace*) tmpin.fparamshape->Get("ws_inputs");

    theInput.push_back(tmpin);

  }
  
  double masses[theInput.size()];
  
  for(int c =0; c< NCAT; c++){

    TCanvas* cc1 = new TCanvas(Form("cc1_cat%d", c), Form("cc1_cat%d", c));

    std::map< std::string , RooDataSet* > data;
    std::map< std::string , RooAbsPdf* > shape;
    std::map< std::string , RooRealVar* > mgg ;
    std::map< std::string , RooRealVar* > MH ;
    RooPlot* p[theInput.size()]; 

    for(unsigned int iM =0; iM < theInput.size(); iM++){

      theTree tmpin = theInput[iM];
      masses[iM] = std::stod(tmpin.M_bins); 
      std::cout << "--------------------------------- " << std::endl;
      std::cout << "Mass point: " << masses[iM] << std::endl;

      mgg[tmpin.name] = tmpin.wsparamshape->var("mgg");
      MH[tmpin.name] = tmpin.wsparamshape->var("MH");

      MH[tmpin.name]->setVal(masses[iM]); 

      data[tmpin.name] = (RooDataSet*) tmpin.wsreco->data(TString::Format("SigWeight_cat%d",c));
      data[tmpin.name]->Print("v");

      p[iM] = (RooPlot*) mgg[tmpin.name]->frame(RooFit::Range(300, upperxmax ) );
      data[tmpin.name]->plotOn(p[iM], MarkerColor(1+iM));

      if(c==0) {
	shape[tmpin.name] = tmpin.wsparamshape->pdf(Form("SignalShape_kMpl%s_EBEB", coupling.c_str()));
	shape[tmpin.name]->plotOn(p[iM], LineColor(iM+1));				    
      }
      if(c==1){ 
	shape[tmpin.name] = tmpin.wsparamshape->pdf(Form("SignalShape_kMpl%s_EBEE", coupling.c_str()));
	shape[tmpin.name]->plotOn(p[iM], LineColor(iM+1));
      }				    



      if (iM == 0){
	p[iM]->GetXaxis()->SetTitle("m_{#gamma#gamma}");
	p[iM]->GetYaxis()->SetTitle("a.u.");
	cc1->cd();
	p[iM]->Draw();
      } else{
	cc1->cd();
	p[iM]->Draw("same");
      }

    } //end of loop over signal samples of the same coupling

    cc1->SaveAs(Form("%s/MassesCheck_%s_cat%d.png", ws_outdir.c_str(), couplingIn.c_str(), c));
    cc1->SetLogy(); 
    cc1->SaveAs(Form("%s/MassesCheck_%s_cat%d_log.png", ws_outdir.c_str(), couplingIn.c_str(), c));
    cc1->SetLogy(0);

  } //end of loop over categories

   
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
