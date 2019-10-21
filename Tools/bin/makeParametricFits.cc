#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"

//RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

//ROOT
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLeaf.h"

//-----------------------------------------------------------------------------------
struct theTree {

  std::string kMpl;
  int M_bins;
  std::string name;
  TTree * selectedTree;

};

//-----------------------------------------------------------------------------------
//Declarations here definition after main
RooRealVar* buildRooVar(std::string name, std::string title, int nBins, double xMin, double xMax, std::string unit);
void runAllFits(const std::string &year, const std::string &ws_dir);
void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample);
void AddSigData(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &isample);
std::vector<theTree> readtrees(const std::string &year);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);

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
  //========================================================================
  //create datasets
  // createdatasets(year,"input/workspaces");
  //========================================================================



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
void AddSigData(RooWorkspace* w, Float_t mass, std::string coupling, const std::string &isample) {

  Int_t ncat = 2;//BB and BE for the moment

  RooRealVar* mgg = buildRooVar("mgg", "m_{#gamma#gamma}", 2385, 230., 10000., "GeV");
  RooRealVar* mggGen = buildRooVar("mggGen", "m_{#gamma#gamma} gen", 2385, 230., 10000, "GeV");
  RooRealVar* weight = buildRooVar("weight", "event weight", 0, 0 , 1000 , "nounits");
  RooRealVar* eventClass = buildRooVar("eventClass", "eventClass", 0, -10, 10, "nounits");
  
  RooArgSet* rooVars = new RooArgSet(*mgg, *mggGen, *eventClass, *weight); 

  //Cut on the signal region
  TString CutSignalRegion = TString::Format("mgg>=300&&mgg<=10000&&mggGen>=300&&mggGen<=10000");  
  //Cut on the reduced mass
  TString CutReducedMass = TString::Format("(mgg-mggGen)>-600.&&(mgg-mggGen)<600."); 
  
  // RooDataSet* sigWeightedK = new RooDataSet("sigWeightedK","datasetK", treesforfit[getBase(isample)], *rooVars, CutSignalRegion, "weight" );
  RooDataSet* sigWeightedK = new RooDataSet("sigWeightedK","datasetK", treesforfit[getBase(isample)], *rooVars, CutSignalRegion, "" );

  RooDataSet* signalK[ncat];
  RooDataSet* signalAllK;
  
  for (int c=0; c<ncat; ++c) {
    if (c==0) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),CutReducedMass+TString::Format("&& eventClass==0"));
    if (c==1) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),CutReducedMass+TString::Format("&& eventClass==1"));
       
    TString myCut;
    w->import(*signalK[c],RooFit::Rename(TString::Format("SigWeightK_cat%d",c)));
    std::cout << "cat " << c << ", signalK[c]: " << std::endl;
    signalK[c]->Print("V");
    std::cout << "---- for category " << c << ", nX for signal[c]:  " << signalK[c]->sumEntries() << std::endl; 
  }
    
    // Create full weighted signal data set without categorization
    signalAllK = (RooDataSet*) sigWeightedK->reduce(RooArgList(*w->var("mgg")),CutReducedMass);
    w->import(*signalAllK, RooFit::Rename("SigWeightK"));
    std::cout << "now signalAllK" << std::endl;
    signalAllK->Print("V");
    std::cout << "---- nX for signalAll:  " << signalAllK->sumEntries() << std::endl; 
    std::cout << "==================================================================" << std::endl;


}

void runfits(const std::string &year, const std::string &ws_dir, const std::string &isample){
  
  std::string coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
  std::string M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");
  std::cout << ws_dir << std::endl;
  TFile *fout = new TFile(Form("%s/ws_ResponseAndGen_M%s_k%s_%s.root", ws_dir.c_str(), M_bins.c_str(), coupling.c_str(),year.c_str()), "recreate");
  RooWorkspace* workspace = new RooWorkspace( "HighMassDiphoton", "HighMassDiphoton" );
  
  std::cout << "Adding signal data for mass " << M_bins << " and coupling " << coupling <<std::endl;

  AddSigData(workspace, stof(M_bins), coupling, isample);

  TCanvas* canv = new TCanvas("canv","c",1);
  canv->cd();
  RooPlot* p = workspace->var("mgg")->frame(RooFit::Range(230., 10000.), RooFit::Bins( 2385));

  // TString fileBaseName("HighMassGG");    
  // TString fileBkgName("HighMassGG.inputbkg");
  // HLFactory hlf("HLFactory", "HighMassGG.rs", false);
  // RooWorkspace* w = hlf.GetWs();
  // int iMass = (int) abs(mass);   
  // // range for the variables
  // w->var("mgg")->setMin(MINmass);
  // w->var("mgg")->setMax(MAXmass);
  // w->var("mggGen")->setMin(MINmass);
  // w->var("mggGen")->setMax(MAXmass);
  // w->Print("V");
  
  // cout << endl; 
  // cout << "Now add signal data" << endl;
  // std::cout<<"------> "<<coupling<<std::endl;



  // AddSigData(w, mass, coupling, whichRel);
  
  // TCanvas* canv = new TCanvas("canv","c",1);
  // canv->cd();
  // RooPlot* p = w->var("mgg")->frame(Range(300,3000), Bins(260));

  // bool makeWs=false; 
  // if(makeWs){
  //   if(coupling=="001")sigModelResponseFcnFit(w, mass, coupling, whichRel);
  //   sigModelGenFcnFit(w, mass, coupling, whichRel);
  //   w->Print("V");
  //   TFile* f = new TFile(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k"+coupling+".root", iMass), "RECREATE");
  //   f->cd();
  //   w->Write();
  //   f->Write();
  //   f->Close();
  // }else{
  //   TFile* f001 = TFile::Open(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k001.root", iMass));
  //   TFile* f = TFile::Open(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k"+coupling+".root", iMass));

  //   RooWorkspace* ws001 = (RooWorkspace*) f001->Get("HLFactory_ws"); 
  //   RooWorkspace* ws = (RooWorkspace*) f->Get("HLFactory_ws"); 
  //   sigModelShapeFcnFit(ws001,ws, mass, coupling, whichRel);
  //   asimovDatasetFcnFit(ws, mass, coupling, whichRel);
  //   TFile* fout = new TFile(TString::Format(whichRel+"/ws_ResponseAndGen_M%d_k"+coupling+"_final.root", iMass), "RECREATE");
  //   fout->cd();
  //   ws->Write();
  //   fout->Write();
  //   fout->Close();
  //   }
 
}



// void createdatasets(const std::string &year, const std::string &ws_dir)
// {

//   Int_t ncat = 2;//BB and BE for the moment

//   TFile *fout = new TFile(Form("%s/HighMassDiphoton%s.root", ws_dir.c_str(), year.c_str()), "recreate");
//   RooWorkspace* workspace = new RooWorkspace( "HighMassDiphoton", "HighMassDiphoton" );

//   std::map<std::string, RooDataSet * > sigWeightedK;
 
//   RooRealVar* mgg = buildRooVar("mgg", "m_{#gamma#gamma}", 2385, 230., 6000., "GeV");
//   RooRealVar* mggGen = buildRooVar("mggGen", "m_{#gamma#gamma} gen", 2385, 230., 6000, "GeV");
//   RooRealVar* weight = buildRooVar("weight", "event weight", 0, 0 , 1000 , "nounits");
//   RooRealVar* eventClass = buildRooVar("eventClass", "eventClass", 0, -10, 10, "nounits");
  
//   RooArgSet* rooVars = new RooArgSet(*mgg, *mggGen, *eventClass, *weight); 

//   //Cuts
//   TString CutSignalRegion = TString::Format("mgg>=300&&mgg<=6000&&mggGen>=300&&mggGen<=6000");  
//   //Cut on the reduced mass
//   TString CutReducedMass = TString::Format("(mgg-mggGen)>-600.&&(mgg-mggGen)<600."); 

//   std::vector<std::string> samples = getSampleListForFit();

//   for(auto isample : samples) {
//     //We will process one year each time
//     if ( isample.find(year) == std::string::npos ) continue; 
//     std::cout << "Creating datasets for sample " << isample << std::endl;
    

//     //========================================================================
//     //Build datasets
//     //========================================================================
//     datasets[getBase(isample)] = new RooDataSet( getBase(isample).c_str(), getBase(isample).c_str(), treesforfit[getBase(isample)], *rooVars, CutSignalRegion, "weight" );
//     // std::cout << datasets[getBase(isample)]->GetName() <<  " numEntries " << datasets[getBase(isample)]->numEntries() << std::endl;
    

//     std::cout << "preparing dataset with observable mgg correct k" << std::endl;
//     RooDataSet* signalK[NCAT];
//     std::cout<<mass<<std::endl;
//     RooDataSet* signalAllK;
//     std::cout<<coupling<<std::endl;
    
//      for (int c=0; c<ncat; ++c) {
//       if (c==0) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==0"));
//       if (c==1) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==1"));
//       if (c==2) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==2"));
//       if (c==3) signalK[c] = (RooDataSet*) sigWeightedK->reduce(*w->var("mgg"),mainCut+TString::Format("&& eventClass==3"));
      
//       TString myCut;
//       w->import(*signalK[c],Rename(TString::Format("SigWeightK_cat%d",c)));
//       cout << "cat " << c << ", signalK[c]: " << endl;
//       signalK[c]->Print("V");
//       cout << "---- for category " << c << ", nX for signal[c]:  " << signalK[c]->sumEntries() << endl; 
//       cout << endl;
//     }
    
//     // Create full weighted signal data set without categorization
//     signalAllK = (RooDataSet*) sigWeightedK->reduce(RooArgList(*w->var("mgg")),mainCut);
//     w->import(*signalAllK, Rename("SigWeightK"));
//     cout << "now signalAllK" << endl;
//     signalAllK->Print("V");
//     cout << "---- nX for signalAll:  " << signalAllK->sumEntries() << endl; 
//     cout << endl;





   
//   }

//   fout->Write();
//   fout->Close();


//   // for (const auto& tr : theSelectedTrees) {
//   //   datasets[tr.name] = new RooDataSet( tr.name.c_str(), tr.name.c_str(), tr.selectedTree, *rooVars, CutSignalRegion, "weight");
//   // }
  
//   // //========================================================================
//   // //Build datasets
//   // //========================================================================
//   // // datasets[getBase(isample)] = new RooDataSet( getBase(isample).c_str(), getBase(isample).c_str(), RooArgSet(mgg), RooFit::Import( *finaltree ) );
//   // // std::cout << datasets[getBase(isample)]->GetName() <<  " numEntries " << datasets[getBase(isample)]->numEntries() << std::endl;

//   // fout->Write();
//   // fout->Close();

// }



// std::vector<theTree> readtrees(const std::string &year)
// {
  
//   std::vector<theTree> thtr; 
//   thtr.clear();
  
//   std::vector<std::string> samples = getSampleListForFit();

//   for(auto isample : samples) {
//     std::cout << isample << std::endl;
//   }

//   theTree tmpfortree;
//   for(auto isample : samples) {

//     std::cout << "Creating datasets for sample " << isample << std::endl;


//     tmpfortree.kMpl = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
//     tmpfortree.M_bins = std::stoi( get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_") );
//     tmpfortree.name = getBase(isample);
//     tmpfortree.selectedTree = treesforfit[getBase(isample)];

//     thtr.push_back(tmpfortree);
  
//   }

//   return thtr;

// }

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

