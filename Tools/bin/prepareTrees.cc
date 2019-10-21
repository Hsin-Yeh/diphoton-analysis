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
//Declarations here definition after main
void prepare(const std::string &region, const std::string &year, const std::string &outputdir);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

  std::string region, year, outputdir;

  if(argc!=4) {
    std::cout << "Syntax: prepareTrees.exe [BB/BE] [2016/2017/2018] [OutputDir]" << std::endl;
    return -1;
  }
  else {
    region = argv[1];
    if(region!="BB" and region!="BE") {
      std::cout << "Only 'BB' and 'BE' are allowed regions. " << std::endl;
      return -1;
    }
    year = argv[2];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
    outputdir = argv[3];
  }

  // include signal samples but not unskimmed data samples
  init(false, true);

  //========================================================================
  prepare(region, year, outputdir);

}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
void prepare(const std::string &region, const std::string &year, const std::string &outputdir)
{
  
  std::map<std::string, std::string> cuts;
  cuts["BB"] = "isGood*(Diphoton.Minv > 230 && Diphoton.deltaR > 0.45 && Photon1.pt>125 && Photon2.pt>125 && Photon1.isEB && Photon2.isEB)";
  cuts["BE"] = "isGood*(Diphoton.Minv > 330 && Diphoton.deltaR > 0.45 && Photon1.pt>125 && Photon2.pt>125 && ( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE )))";

  std::vector<std::string> samples = getSampleList();

  std::vector<int> stringScales = {3000, 3500, 4000, 4500, 5000, 5500, 6000};

  for(auto isample : samples) {
    std::cout << isample << std::endl;
  }

  for(auto isample : samples) {
    std::string sampleCut = cuts[region];
    if ( isample.find(year) == std::string::npos ) continue; 
    //Run only on RS samples for now. 
    if( isample.find("RSGravitonToGammaGamma") == std::string::npos ) continue;
    // apply weights for all samples except data
    bool is2015or2016 = isample.find("2015") != std::string::npos || isample.find("2016") != std::string::npos;
    if( isample.find("ADD") != std::string::npos
    	|| isample.find("gg70") != std::string::npos
    	|| is2015or2016 ) {
      sampleCut += "*(HLT_DoublePhoton60>0 || HLT_ECALHT800>0)";
    }
    else {
      sampleCut += "*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)";
    }
    if( isample.find("data") == std::string::npos ) {
      sampleCut+="*weightAll*" + std::to_string(luminosity[year]);
      // need to increase selection for ADD cuts to avoid negative weights
      // from background subtraction
      if( isample.find("ADD") != std::string::npos
	  or isample.find("gg70") != std::string::npos ) {
	sampleCut += "*(Diphoton.Minv > 600)";
      }
    }
    else {
      sampleCut += "*(Diphoton.Minv < 1000)";
    }
    // apply k-factor to Sherpa GG sample
    if( isample.find("gg_R2F2_") != std::string::npos) {
      if( is2015or2016 ) sampleCut += "*" + kfactorString(region, "R2F2_125GeV_CT10");
      else sampleCut += "*" + kfactorString(region, "R2F2_125GeV_NNPDF");
    }
    else if( isample.find("gg_R0p5F0p5_") != std::string::npos) {
      if( is2015or2016 ) sampleCut += "*" + kfactorString(region, "R0p5F0p5_125GeV_CT10");
      else sampleCut += "*" + kfactorString(region, "R0p5F0p5_125GeV_NNPDF");
    }
    else if( isample.find("gg_") != std::string::npos) {
      if( is2015or2016 ) sampleCut += "*" + kfactorString(region, "R1F1_125GeV_CT10");
      else sampleCut += "*" + kfactorString(region, "R1F1_125GeV_NNPDF");
    }

    std::cout << "Creating tree for sample " << isample << " with cut\n" << sampleCut << std::endl;
    std::string baseName(getSampleBase(isample, year));
    // std::cout << "---> "  << getBase(isample) << " " << baseName.c_str() << " " << sampleCut.c_str() << std::endl;
    // std::cout << "Chain " << chains[getBase(isample)]->	GetNtrees() << std::endl;

    //========================================================================
    //Skim
    //========================================================================
    //faster to skim the initial chain first a little bit
    chains[getBase(isample)]->SetBranchStatus("*",0);
    chains[getBase(isample)]->SetBranchStatus("Event",1);
    chains[getBase(isample)]->SetBranchStatus("Diphoton",1);
    chains[getBase(isample)]->SetBranchStatus("GenDiphoton",1);
    chains[getBase(isample)]->SetBranchStatus("isGood",1);
    chains[getBase(isample)]->SetBranchStatus("Photon1",1);
    chains[getBase(isample)]->SetBranchStatus("Photon2",1);

    TTree *newtree1 = chains[getBase(isample)]->CloneTree(0);
    newtree1->CopyEntries(chains[getBase(isample)]);
    // newtree1->Print();
    
    //========================================================================
    //Selection
    //========================================================================
    TTree *newtree2 = newtree1->CopyTree(cuts[region].c_str());
    // newtree2->Print();

    newtree2->SetBranchStatus("*",0);
    newtree2->SetBranchStatus("Event",1);
    newtree2->SetBranchStatus("Diphoton",1);
    newtree2->SetBranchStatus("GenDiphoton",1);

    TTree *finaltree = newtree2->CloneTree(0);
    finaltree->CopyEntries(newtree2);
    // finaltree->Print();

    //========================================================================
    //RooFit preparation
    //========================================================================
    //For the dataset creation in RooFit we need to have a tree-branch structure
    //where the variable must be a branch. So, we will create the for each of 
    //initial branches that we want (Diphoton, GenDiphoton, Event etc) a tree
    //with its variables as branches. 
    TFile *fout = new TFile(Form("%s/%s_%s.root", outputdir.c_str(), isample.c_str(), region.c_str()), "recreate");
    TTree *HighMassDiphotonTree = new TTree("HighMassDiphoton" , "A tree with all the info needed for RooFit");

    //Event branch
    Long64_t run, LS, evnum, processid, bx, orbit;
    Float_t ptHat, alphaqcd, alphaqed, qscale, x1, x2, pdf1, pdf2, weight0, weight, weightPuUp, weightPu, weightPuDown, weightLumi, weightAll; 
    Int_t interactingParton1PdgId, interactingParton2PdgId, pdf_id1, pdf_id2, npv_true;
    Bool_t beamHaloIDLoose, beamHaloIDTight, beamHaloIDTight2015;
    //Diphoton branch
    Double_t mgg,qt,deltaPhi,deltaEta,deltaR,cosThetaStar,cosThetaStar_old,chiDiphoton;
    Bool_t isEBEB,isEBEE,isEEEB,isEEEE;
    //GenDiphoton branch
    Double_t mggGen,Genqt,GendeltaPhi,GendeltaEta,GendeltaR,GencosThetaStar,GencosThetaStar_old,GenchiDiphoton;
    Bool_t GenisEBEB,GenisEBEE,GenisEEEB,GenisEEEE;
    //One more branch for the category
    // BB: 0 , BE:1
    Int_t eventClass; 
    
    //setBranches
    //Event branches
    HighMassDiphotonTree->Branch("run",&run,"run/L");
    HighMassDiphotonTree->Branch("LS",&LS,"LS/L");
    HighMassDiphotonTree->Branch("evnum",&evnum,"evnum/L");
    HighMassDiphotonTree->Branch("processid",&processid,"processid/L");
    HighMassDiphotonTree->Branch("bx",&bx,"bx/L");
    HighMassDiphotonTree->Branch("orbit",&orbit,"orbit/L");
    HighMassDiphotonTree->Branch("ptHat",&ptHat,"ptHat/F");
    HighMassDiphotonTree->Branch("alphaqcd",&alphaqcd,"alphaqcd/F");
    HighMassDiphotonTree->Branch("alphaqed",&alphaqed,"alphaqed/F");
    HighMassDiphotonTree->Branch("qscale",&qscale,"qscale/F");
    HighMassDiphotonTree->Branch("x1",&x1,"x1/F");
    HighMassDiphotonTree->Branch("x2",&x2,"x2/F");
    HighMassDiphotonTree->Branch("pdf1",&pdf1,"pdf1/F");
    HighMassDiphotonTree->Branch("pdf2",&pdf2,"pdf2/F");
    HighMassDiphotonTree->Branch("weight0",&weight0,"weight0/F");
    HighMassDiphotonTree->Branch("weight",&weight,"weight/F");
    HighMassDiphotonTree->Branch("weightPuUp",&weightPuUp,"weightPuUp/F");
    HighMassDiphotonTree->Branch("weightPu",&weightPu,"weightPu/F");
    HighMassDiphotonTree->Branch("weightPuDown",&weightPuDown,"weightPuDown/F");
    HighMassDiphotonTree->Branch("weightLumi",&weightLumi,"weightLumi/F");
    HighMassDiphotonTree->Branch("weightAll",&weightAll,"weightAll/F");
    HighMassDiphotonTree->Branch("interactingParton1PdgId",&interactingParton1PdgId,"interactingParton1PdgId/I");
    HighMassDiphotonTree->Branch("interactingParton2PdgId",&interactingParton2PdgId,"interactingParton2PdgId/I");
    HighMassDiphotonTree->Branch("pdf_id1",&pdf_id1,"pdf_id1/I");
    HighMassDiphotonTree->Branch("pdf_id2",&pdf_id2,"pdf_id2/I");
    HighMassDiphotonTree->Branch("npv_true ",&npv_true,"npv_true/I");
    HighMassDiphotonTree->Branch("beamHaloIDLoose",&beamHaloIDLoose,"beamHaloIDLoose/I");
    HighMassDiphotonTree->Branch("beamHaloIDTight",&beamHaloIDTight,"beamHaloIDTight/I");
    HighMassDiphotonTree->Branch("beamHaloIDTight2015",&beamHaloIDTight2015,"beamHaloIDTight2015/I");
    //Diphoton branches
    HighMassDiphotonTree->Branch("mgg",&mgg,"mgg/D");
    HighMassDiphotonTree->Branch("qt",&qt,"qt/D");
    HighMassDiphotonTree->Branch("deltaPhi",&deltaPhi,"deltaPhi/D");
    HighMassDiphotonTree->Branch("deltaEta",&deltaEta,"deltaEta/D");
    HighMassDiphotonTree->Branch("deltaR",&deltaR,"deltaR/D");
    HighMassDiphotonTree->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/D");
    HighMassDiphotonTree->Branch("cosThetaStar_old",&cosThetaStar_old,"cosThetaStar_old/D");
    HighMassDiphotonTree->Branch("chiDiphoton",&chiDiphoton,"chiDiphoton/D");
    HighMassDiphotonTree->Branch("isEBEB",&isEBEB,"isEBEB/O");
    HighMassDiphotonTree->Branch("isEBEE",&isEBEE,"isEBEE/O");
    HighMassDiphotonTree->Branch("isEEEB",&isEEEB,"isEEEB/O");
    HighMassDiphotonTree->Branch("isEEEE",&isEEEE,"isEEEE/O");
    //GenDiphoton branches
    HighMassDiphotonTree->Branch("mggGen",&mggGen,"mggGen/D");
    HighMassDiphotonTree->Branch("Genqt",&Genqt,"Genqt/D");
    HighMassDiphotonTree->Branch("GendeltaPhi",&GendeltaPhi,"GendeltaPhi/D");
    HighMassDiphotonTree->Branch("GendeltaEta",&GendeltaEta,"GendeltaEta/D");
    HighMassDiphotonTree->Branch("GendeltaR",&GendeltaR,"GendeltaR/D");
    HighMassDiphotonTree->Branch("GencosThetaStar",&GencosThetaStar,"GencosThetaStar/D");
    HighMassDiphotonTree->Branch("GencosThetaStar_old",&GencosThetaStar_old,"GencosThetaStar_old/D");
    HighMassDiphotonTree->Branch("GenchiDiphoton",&GenchiDiphoton,"GenchiDiphoton/D");
    HighMassDiphotonTree->Branch("GenisEBEB",&GenisEBEB,"GenisEBEB/O");
    HighMassDiphotonTree->Branch("GenisEBEE",&GenisEBEE,"GenisEBEE/O");
    HighMassDiphotonTree->Branch("GenisEEEB",&GenisEEEB,"GenisEEEB/O");
    HighMassDiphotonTree->Branch("GenisEEEE",&GenisEEEE,"GenisEEEE/O");
    //Category branch 
    HighMassDiphotonTree->Branch("eventClass",&eventClass,"eventClass/I");

    for (Int_t i=0;i<finaltree->GetEntries();i++) {
      finaltree->GetEntry(i);
      //Event branch
      run			= finaltree->GetBranch("Event")->GetLeaf("run")->GetValue(0);
      LS			= finaltree->GetBranch("Event")->GetLeaf("LS")->GetValue(0);
      evnum			= finaltree->GetBranch("Event")->GetLeaf("evnum")->GetValue(0);
      processid	         	= finaltree->GetBranch("Event")->GetLeaf("processid")->GetValue(0);
      bx			= finaltree->GetBranch("Event")->GetLeaf("bx")->GetValue(0);
      orbit			= finaltree->GetBranch("Event")->GetLeaf("orbit")->GetValue(0);
      ptHat			= finaltree->GetBranch("Event")->GetLeaf("ptHat")->GetValue(0);
      alphaqcd		        = finaltree->GetBranch("Event")->GetLeaf("alphaqcd")->GetValue(0);
      alphaqed		        = finaltree->GetBranch("Event")->GetLeaf("alphaqed")->GetValue(0);
      qscale			= finaltree->GetBranch("Event")->GetLeaf("qscale")->GetValue(0);
      x1			= finaltree->GetBranch("Event")->GetLeaf("x1")->GetValue(0);
      x2			= finaltree->GetBranch("Event")->GetLeaf("x2")->GetValue(0);
      pdf1			= finaltree->GetBranch("Event")->GetLeaf("pdf1")->GetValue(0);
      pdf2			= finaltree->GetBranch("Event")->GetLeaf("pdf2")->GetValue(0);
      weight0			= finaltree->GetBranch("Event")->GetLeaf("weight0")->GetValue(0);
      weight			= finaltree->GetBranch("Event")->GetLeaf("weight")->GetValue(0);
      weightPuUp		= finaltree->GetBranch("Event")->GetLeaf("weightPuUp")->GetValue(0);
      weightPu		        = finaltree->GetBranch("Event")->GetLeaf("weightPu")->GetValue(0);
      weightPuDown		= finaltree->GetBranch("Event")->GetLeaf("weightPuDown")->GetValue(0);
      weightLumi		= finaltree->GetBranch("Event")->GetLeaf("weightLumi")->GetValue(0);
      weightAll 		= finaltree->GetBranch("Event")->GetLeaf("weightAll")->GetValue(0);
      interactingParton1PdgId	= finaltree->GetBranch("Event")->GetLeaf("interactingParton1PdgId")->GetValue(0);
      interactingParton2PdgId	= finaltree->GetBranch("Event")->GetLeaf("interactingParton2PdgId")->GetValue(0);
      pdf_id1			= finaltree->GetBranch("Event")->GetLeaf("pdf_id1")->GetValue(0);
      pdf_id2			= finaltree->GetBranch("Event")->GetLeaf("pdf_id2")->GetValue(0);
      npv_true		        = finaltree->GetBranch("Event")->GetLeaf("npv_true")->GetValue(0);
      beamHaloIDLoose		= finaltree->GetBranch("Event")->GetLeaf("beamHaloIDLoose")->GetValue(0);
      beamHaloIDTight		= finaltree->GetBranch("Event")->GetLeaf("beamHaloIDTight")->GetValue(0);
      beamHaloIDTight2015       = finaltree->GetBranch("Event")->GetLeaf("beamHaloIDTight2015")->GetValue(0);
      //Diphoton branch
      mgg	        	= finaltree->GetBranch("Diphoton")->GetLeaf("Minv")->GetValue(0);
      qt		        = finaltree->GetBranch("Diphoton")->GetLeaf("qt")->GetValue(0);
      deltaPhi	                = finaltree->GetBranch("Diphoton")->GetLeaf("deltaPhi")->GetValue(0);
      deltaEta	                = finaltree->GetBranch("Diphoton")->GetLeaf("deltaEta")->GetValue(0);
      deltaR		        = finaltree->GetBranch("Diphoton")->GetLeaf("deltaR")->GetValue(0);
      cosThetaStar	        = finaltree->GetBranch("Diphoton")->GetLeaf("cosThetaStar")->GetValue(0);
      cosThetaStar_old          = finaltree->GetBranch("Diphoton")->GetLeaf("cosThetaStar_old")->GetValue(0);
      chiDiphoton	        = finaltree->GetBranch("Diphoton")->GetLeaf("chiDiphoton")->GetValue(0);
      isEBEB		        = finaltree->GetBranch("Diphoton")->GetLeaf("isEBEB")->GetValue(0);
      isEBEE		        = finaltree->GetBranch("Diphoton")->GetLeaf("isEBEE")->GetValue(0);
      isEEEB		        = finaltree->GetBranch("Diphoton")->GetLeaf("isEEEB")->GetValue(0);
      isEEEE                    = finaltree->GetBranch("Diphoton")->GetLeaf("isEEEE")->GetValue(0);
      //GenDiphoton branch
      mggGen	        	= finaltree->GetBranch("GenDiphoton")->GetLeaf("Minv")->GetValue(0);
      Genqt		        = finaltree->GetBranch("GenDiphoton")->GetLeaf("qt")->GetValue(0);
      GendeltaPhi	        = finaltree->GetBranch("GenDiphoton")->GetLeaf("deltaPhi")->GetValue(0);
      GendeltaEta	        = finaltree->GetBranch("GenDiphoton")->GetLeaf("deltaEta")->GetValue(0);
      GendeltaR		        = finaltree->GetBranch("GenDiphoton")->GetLeaf("deltaR")->GetValue(0);
      GencosThetaStar	        = finaltree->GetBranch("GenDiphoton")->GetLeaf("cosThetaStar")->GetValue(0);
      GencosThetaStar_old       = finaltree->GetBranch("GenDiphoton")->GetLeaf("cosThetaStar_old")->GetValue(0);
      GenchiDiphoton	        = finaltree->GetBranch("GenDiphoton")->GetLeaf("chiDiphoton")->GetValue(0);
      GenisEBEB		        = finaltree->GetBranch("GenDiphoton")->GetLeaf("isEBEB")->GetValue(0);
      GenisEBEE		        = finaltree->GetBranch("GenDiphoton")->GetLeaf("isEBEE")->GetValue(0);
      GenisEEEB		        = finaltree->GetBranch("GenDiphoton")->GetLeaf("isEEEB")->GetValue(0);
      GenisEEEE                 = finaltree->GetBranch("GenDiphoton")->GetLeaf("isEEEE")->GetValue(0);
      //Event category
      if (region == "BB"){
	eventClass = 0;
      } else if (region == "BE"){
	eventClass = 1;
      }
      
      //Fill the final tree
      HighMassDiphotonTree->Fill();
    }

      
    fout->Write();
    fout->Close();

  }

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

