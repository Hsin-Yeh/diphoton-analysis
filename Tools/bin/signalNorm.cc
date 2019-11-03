#include "diphoton-analysis/Tools/interface/sampleList.hh"
#include "diphoton-analysis/Tools/interface/utilities.hh"

//ROOT
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraphErrors.h"

//-----------------------------------------------------------------------------------
struct eff_reco {

  std::string kMpl;
  std::string M_bins;
  //Values for sel are: no_selection, isGood
  std::string sel; 
  //Values for cat are both, BB, BE. 
  std::string cat; 
  double eff; 
  
};

//-----------------------------------------------------------------------------------
//Declarations here definition after main
std::map<std::string, std::vector<eff_reco> > computefficiency(const std::string &year);
void plotefficiency(TCanvas* cc, std::map<std::string, std::vector<eff_reco> > coup_eff_reco, bool doave);
std::map<std::string, std::vector<eff_reco> > reduce(std::map<std::string, std::vector<eff_reco> > coup_eff_reco, std::string cut);
TH1F * createHisto(TChain * newtree1, const std::string &histo_name, int nBins, double xMin, double xMax, std::string cut);
std::string getSampleBase(const std::string & sampleName, const std::string & year);
std::string getBase(const std::string & sampleName);
std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

  std::string region, year, outputdir;

  if(argc!=2) {
    std::cout << "Syntax: signalNorm.exe [2016/2017/2018] " << std::endl;
    return -1;
  }
  else {
    year = argv[1];
    if(year!="2016" and year!="2017" and year!="2018") {
      std::cout << "Only '2016', 2017' and '2018' are allowed years. " << std::endl;
      return -1;
    }
  }

  // include signal samples but not unskimmed data samples
  init(false, true);

  //========================================================================
  std::map<std::string, std::vector<eff_reco> > effreco = computefficiency(year);
  std::map<std::string, std::vector<eff_reco> > effrecoBB = reduce(effreco, "BB");
  std::map<std::string, std::vector<eff_reco> > effrecoBE = reduce(effreco, "BE");
  TCanvas* cc = new TCanvas("cc", "cc");
  plotefficiency(cc, effrecoBB, false);
  cc->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effVsMass_BB.png");
  TCanvas* ccBE = new TCanvas("ccBE", "ccBE");
  plotefficiency(ccBE, effrecoBE, false);
  ccBE->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effVsMass_BE.png");
  TCanvas* ccave = new TCanvas("ccave", "ccave");
  plotefficiency(ccave, effreco, true);
  ccave->SaveAs("/afs/cern.ch/work/a/apsallid/CMS/Hgg/exodiphotons/CMSSW_9_4_13/src/diphoton-analysis/output/signalNorm/effVsMass_average.png");


}

//-----------------------------------------------------------------------------------
//Definitions
//-----------------------------------------------------------------------------------
std::map<std::string, std::vector<eff_reco> > computefficiency(const std::string &year)
{
  
  //We will follow the efficiency throughout the series of cuts 
  std::map<std::string, std::string> cuts;
  cuts["isGood"] = "isGood"; 
  cuts["isGood_Diphoton.deltaR"] = "isGood*(Diphoton.deltaR > 0.45)"; 
  cuts["isGood_Diphoton.deltaR_Photon1.pt"] = "isGood*(Diphoton.deltaR > 0.45 && Photon1.pt>125)"; 
  cuts["isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt"] = "isGood*(Diphoton.deltaR > 0.45 && Photon1.pt>125 && Photon2.pt>125)"; 
  cuts["isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt_HLT"] = "isGood*(Diphoton.deltaR > 0.45 && Photon1.pt>125 && Photon2.pt>125)*(HLT_DoublePhoton70>0 || HLT_ECALHT800>0)"; 

  cuts["BB"] = "isGood*(Diphoton.Minv > 230 && Diphoton.deltaR > 0.45 && Photon1.pt>125 && Photon2.pt>125 && Photon1.isEB && Photon2.isEB)";
  cuts["BE"] = "isGood*(Diphoton.Minv > 330 && Diphoton.deltaR > 0.45 && Photon1.pt>125 && Photon2.pt>125 && ( (Photon1.isEB && Photon2.isEE) || (Photon2.isEB &&  Photon1.isEE )))";

  std::vector<std::string> samples = getSampleList();

  std::vector<int> stringScales = {3000, 3500, 4000, 4500, 5000, 5500, 6000};

  std::string coupling = "";
  std::string M_bins = "";
  double Ngen = 0.;

  std::map<std::string, std::vector<eff_reco> > the_eff_reco; 
  std::map<std::string, TH1F *> histograms;
  std::string histo_name;

  for(auto isample : samples) {
    //Run a single year each time
    if ( isample.find(year) == std::string::npos ) continue; 
    //Run only on RS samples for now. 
    if( isample.find("RSGravitonToGammaGamma") == std::string::npos ) continue;

    std::cout << "Prosessing sample " << isample << " for year " << year << std::endl;
    coupling = get_str_between_two_str(getBase(isample), "kMpl", "_M_");
    M_bins = get_str_between_two_str(getBase(isample), "_M_", "_TuneCP2_13TeV_");

    // if (std::stod(M_bins) < 4000.){continue;}

    int nBins = 2385;
    double xMin = 0.0;
    double xMax = 0.0; 
    if ( coupling == "001" ){xMax = 6000. ;}
    else if ( coupling == "01" ){xMax = 9000. ;}
    else if ( coupling == "02" ){xMax = 9000. ;}
    else {
      std::cout << "Only 'kMpl001', 'kMpl01' and 'kMpl02' are allowed. " << std::endl;
      exit(1);
    }

    eff_reco tmpeffreco;

    //========================================================================
    //Skim
    //========================================================================
    //faster to skim the initial chain first a little bit
    // chains[getBase(isample)]->SetBranchStatus("*",0);
    // chains[getBase(isample)]->SetBranchStatus("Event",1);
    // chains[getBase(isample)]->SetBranchStatus("TriggerBit",1);
    // chains[getBase(isample)]->SetBranchStatus("TriggerPrescale",1);
    // chains[getBase(isample)]->SetBranchStatus("Diphoton",1);
    // chains[getBase(isample)]->SetBranchStatus("GenDiphoton",1);
    // chains[getBase(isample)]->SetBranchStatus("isGood",1);
    // chains[getBase(isample)]->SetBranchStatus("Photon1",1);
    // chains[getBase(isample)]->SetBranchStatus("Photon2",1);

    // TTree *newtree1 = chains[getBase(isample)]->CloneTree(0);
    // newtree1->CopyEntries(chains[getBase(isample)]);
    // newtree1->Print();

    //========================================================================
    //Selection and cutflow through efficiency
    //========================================================================
    std::string baseName(getSampleBase(isample, year));

    //------------------------------------------------------------------------
    //No selection
    // tmpeffreco.M_bins = M_bins;
    // tmpeffreco.sel = "no_selection";
    // tmpeffreco.cat = "both";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName;
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, "no_selection");
    Ngen = histograms[histo_name]->Integral();

    // tmpeffreco.eff = 1.;

    // the_eff_reco[coupling].push_back(tmpeffreco);   
    // //------------------------------------------------------------------------
    // //Sel: isGood
    // tmpeffreco.M_bins = M_bins;
    // tmpeffreco.sel = "isGood";
    // tmpeffreco.cat = "both";

    // std::cout << "-----------------------------------------------------------------" << std::endl;
    // histo_name = baseName + "_isGood";
    // histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["isGood"]);
    // tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    // the_eff_reco[coupling].push_back(tmpeffreco);   
    // //------------------------------------------------------------------------
    // //Sel: isGood_Diphoton.deltaR
    // tmpeffreco.M_bins = M_bins;
    // tmpeffreco.sel = "isGood_Diphoton.deltaR";
    // tmpeffreco.cat = "both";

    // std::cout << "-----------------------------------------------------------------" << std::endl;
    // histo_name = baseName + "_isGood_Diphoton.deltaR";
    // histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["isGood_Diphoton.deltaR"]);
    // tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    // the_eff_reco[coupling].push_back(tmpeffreco);   
    // //------------------------------------------------------------------------
    // //Sel: isGood_Diphoton.deltaR_Photon1.pt
    // tmpeffreco.M_bins = M_bins;
    // tmpeffreco.sel = "isGood_Diphoton.deltaR_Photon1.pt";
    // tmpeffreco.cat = "both";

    // std::cout << "-----------------------------------------------------------------" << std::endl;
    // histo_name = baseName + "_isGood_Diphoton.deltaR_Photon1.pt";
    // histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["isGood_Diphoton.deltaR_Photon1.pt"]);
    // tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    // the_eff_reco[coupling].push_back(tmpeffreco);   
    // //------------------------------------------------------------------------
    // //Sel: isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt
    // tmpeffreco.M_bins = M_bins;
    // tmpeffreco.sel = "isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt";
    // tmpeffreco.cat = "both";

    // std::cout << "-----------------------------------------------------------------" << std::endl;
    // histo_name = baseName + "_isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt";
    // histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt"]);
    // tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    // the_eff_reco[coupling].push_back(tmpeffreco);   
    // //------------------------------------------------------------------------
    // //Sel: isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt_HLT
    // tmpeffreco.M_bins = M_bins;
    // tmpeffreco.sel = "isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt_HLT";
    // tmpeffreco.cat = "both";

    // std::cout << "-----------------------------------------------------------------" << std::endl;
    // histo_name = baseName + "_isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt_HLT";
    // histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["isGood_Diphoton.deltaR_Photon1.pt_Photon2.pt_HLT"]);
    // tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    // the_eff_reco[coupling].push_back(tmpeffreco);   

    //------------------------------------------------------------------------
    //Sel: BB
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "BB";
    tmpeffreco.cat = "BB";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_BB";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["BB"]);
    tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   

    std::cout << tmpeffreco.eff << std::endl; 

    //------------------------------------------------------------------------
    //Sel: BE
    tmpeffreco.M_bins = M_bins;
    tmpeffreco.sel = "BE";
    tmpeffreco.cat = "BE";

    std::cout << "-----------------------------------------------------------------" << std::endl;
    histo_name = baseName + "_BE";
    histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts["BE"]);
    tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    the_eff_reco[coupling].push_back(tmpeffreco);   

    // tmpeffreco.M_bins = M_bins;
    // tmpeffreco.sel = region;
    // tmpeffreco.cat = region;

    // std::cout << "-----------------------------------------------------------------" << std::endl;
    // histo_name = baseName + region;
    // histograms[histo_name] = createHisto(chains[getBase(isample)], histo_name, nBins, xMin, xMax, cuts[region]);
    // tmpeffreco.eff = histograms[histo_name]->Integral() / Ngen;

    // the_eff_reco[coupling].push_back(tmpeffreco);   

  }//end of loop through samples

  return the_eff_reco;

}

//-----------------------------------------------------------------------------------
void plotefficiency(TCanvas* cc, std::map<std::string, std::vector<eff_reco> > effreco, bool doave){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  double masses001[effreco["001"].size()];
  double massesErr001[effreco["001"].size()];
  double eff001[effreco["001"].size()];
  double effErr001[effreco["001"].size()];

  double masses01[effreco["01"].size()];
  double massesErr01[effreco["01"].size()];
  double eff01[effreco["01"].size()];
  double effErr01[effreco["01"].size()];

  double masses02[effreco["02"].size()];
  double massesErr02[effreco["02"].size()];
  double eff02[effreco["02"].size()];
  double effErr02[effreco["02"].size()];

  //Since the final total selection efficiency is the average of the 3 couplings 
  //we will save that also. 
  //BE CAREFUL: For the moment we will only use points where we have values for all the 
  //three couplings
  std::map<std::string, std::vector<eff_reco> > effrecoBB;
  std::map<std::string, std::vector<eff_reco> > effrecoBE;
  if (doave){
    for(unsigned int iE =0; iE < effreco["01"].size(); iE++){
      std::cout << effreco["01"][iE].M_bins << std::endl;
      double curmass = std::stod(effreco["01"][iE].M_bins); 
      if ( curmass > 5000. || curmass == 2750. || curmass == 3250. || curmass == 4250. || curmass == 4500. || curmass == 4750.  ){continue;}
      if ( effreco["01"][iE].sel == "BB"){effrecoBB["01"].push_back( effreco["01"][iE] );}
      if ( effreco["01"][iE].sel == "BE"){effrecoBE["01"].push_back( effreco["01"][iE] );}
    }

    for(unsigned int iE =0; iE < effreco["02"].size(); iE++){
      double curmass = std::stod(effreco["02"][iE].M_bins); 
      if ( curmass > 5000. || curmass == 2750. || curmass == 3250. || curmass == 4250. || curmass == 4500. || curmass == 4750.  ){continue;}
      if ( effreco["02"][iE].sel == "BB"){effrecoBB["02"].push_back( effreco["02"][iE] );}
      if ( effreco["02"][iE].sel == "BE"){effrecoBE["02"].push_back( effreco["02"][iE] );}
    }

    for(unsigned int iE =0; iE < effreco["001"].size(); iE++){
      if ( effreco["001"][iE].sel == "BB"){effrecoBB["001"].push_back( effreco["001"][iE] );}
      if ( effreco["001"][iE].sel == "BE"){effrecoBE["001"].push_back( effreco["001"][iE] );}
    }


  }
  
  double massesBB[effrecoBB["01"].size()];
  double massesErrBB[effrecoBB["01"].size()];
  double effBB[effrecoBB["01"].size()];
  double effErrBB[effrecoBB["01"].size()];

  double massesBE[effrecoBE["01"].size()];
  double massesErrBE[effrecoBE["01"].size()];
  double effBE[effrecoBE["01"].size()];
  double effErrBE[effrecoBE["01"].size()];

  std::map<std::string, double> upperxmax;
  upperxmax["001"] = 6000.;
  upperxmax["01"]  = 9000.;
  upperxmax["02"]  = 9000.;
  
  for(unsigned int iE =0; iE < effreco["001"].size(); iE++){
    masses001[iE] = std::stod(effreco["001"][iE].M_bins); 
    massesErr001[iE] = 0.;
    eff001[iE] = effreco["001"][iE].eff; 
    effErr001[iE] = 0.;
  }

  for(unsigned int iE =0; iE < effreco["01"].size(); iE++){
    masses01[iE] = std::stod(effreco["01"][iE].M_bins); 
    massesErr01[iE] = 0.;
    eff01[iE] = effreco["01"][iE].eff; 
    effErr01[iE] = 0.;
  }
  
  for(unsigned int iE =0; iE < effreco["02"].size(); iE++){
    masses02[iE] = std::stod(effreco["02"][iE].M_bins); 
    massesErr02[iE] = 0.;
    eff02[iE] = effreco["02"][iE].eff; 
    effErr02[iE] = 0.;
  }

  if (doave){
    //BB average
    for(unsigned int iE =0; iE < effrecoBB["01"].size(); iE++){
      std::cout << effrecoBB["01"].size() << " " << effrecoBB["001"].size() <<  " " << effrecoBB["02"].size() << std::endl;
      std::cout << effrecoBB["01"][iE].M_bins << " " << effrecoBB["02"][iE].M_bins << " " << effrecoBB["001"][iE].M_bins << std::endl;
      massesBB[iE] = std::stod(effrecoBB["01"][iE].M_bins); 
      massesErrBB[iE] = 0.;
      //The tweaks to divide by 3 and use all couplings are above
      effBB[iE] = (effrecoBB["001"][iE].eff+effrecoBB["01"][iE].eff+effrecoBB["02"][iE].eff)/3.; 
      effErrBB[iE] = 0.;
    }
    //BE average
    for(unsigned int iE =0; iE < effrecoBE["01"].size(); iE++){
      massesBE[iE] = std::stod(effrecoBE["01"][iE].M_bins); 
      massesErrBE[iE] = 0.;
      effBE[iE] = (effrecoBE["001"][iE].eff+effrecoBE["01"][iE].eff+effrecoBE["02"][iE].eff)/3.; 
      effErrBE[iE] = 0.;
    }
  }
  
  
  cc->cd();

  if (!doave){
    TLegend* legmc = new TLegend(0.58, 0.34, 0.85, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);

    TGraphErrors* graph001 = new TGraphErrors(effreco["001"].size(), masses001, eff001, massesErr001, effErr001);
    graph001->GetYaxis()->SetRangeUser(0., 1.0);
    graph001->SetMarkerColor(2);
    graph001->SetLineColor(2);
    legmc->AddEntry( graph001 , "#frac{#Gamma}{m} = 1.4 #times 10^{-4}" ,"pl" );

    TGraphErrors* graph01 = new TGraphErrors(effreco["01"].size(), masses01, eff01, massesErr01, effErr01);
    graph01->GetYaxis()->SetRangeUser(0., 1.0);
    graph01->SetMarkerColor(4);
    graph01->SetLineColor(4);
    legmc->AddEntry( graph01 , "#frac{#Gamma}{m} = 1.4 #times 10^{-2}" ,"pl" );

    TGraphErrors* graph02 = new TGraphErrors(effreco["02"].size(), masses02, eff02, massesErr02, effErr02);
    graph02->GetYaxis()->SetRangeUser(0., 1.0);
    graph02->SetMarkerColor(8);
    graph02->SetLineColor(8);
    legmc->AddEntry( graph02 , "#frac{#Gamma}{m} = 5.6 #times 10^{-2}" ,"pl" );

 
    TF1 *f001 = new TF1("f001", "pol2", 500,upperxmax["001"]);
    f001->SetLineColor(2);
    graph001->Fit("f001", "R");
    graph001->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
    graph001->GetXaxis()->SetTitle("m_{G} [GeV]");
    graph001->Draw("APE");
    f001->Draw("same");

    TF1 *f01 = new TF1("f01", "pol2", 500,upperxmax["01"]);
    f01->SetLineColor(4);
    graph01->Fit("f01", "R");
    graph01->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
    graph01->GetXaxis()->SetTitle("m_{G} [GeV]");
    graph01->Draw("PSE");
    f01->Draw("same");

    TF1 *f02 = new TF1("f02", "pol2", 500,upperxmax["02"]);
    f02->SetLineColor(8);
    graph02->Fit("f02", "R");
    graph02->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
    graph02->GetXaxis()->SetTitle("m_{G} [GeV]");
    graph02->Draw("PSE");
    f02->Draw("same");

    legmc->Draw("same");
  } else{
    TLegend* legmc = new TLegend(0.58, 0.70, 0.75, 0.9, "", "bNDC");
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);

    TGraphErrors* graphBB = new TGraphErrors(effrecoBB["01"].size(), massesBB, effBB, massesErrBB, effErrBB);
    graphBB->GetYaxis()->SetRangeUser(0., 1.0);
    graphBB->SetMarkerColor(2);
    graphBB->SetLineColor(2);
    legmc->AddEntry( graphBB , "BB" ,"pl" );

    TGraphErrors* graphBE = new TGraphErrors(effrecoBE["01"].size(), massesBE, effBE, massesErrBE, effErrBE);
    graphBE->GetYaxis()->SetRangeUser(0., 1.0);
    graphBE->SetMarkerColor(4);
    graphBE->SetLineColor(4);
    legmc->AddEntry( graphBE , "BE" ,"pl" );

    TF1 *fBB = new TF1("fBB", "pol2", 500,upperxmax["01"]);
    fBB->SetLineColor(2);
    graphBB->Fit("fBB", "R");
    graphBB->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
    graphBB->GetXaxis()->SetTitle("m_{G} [GeV]");
    graphBB->Draw("APE");
    fBB->Draw("same");

    TF1 *fBE = new TF1("fBE", "pol2", 500,upperxmax["01"]);
    fBE->SetLineColor(4);
    graphBE->Fit("fBE", "R");
    graphBE->GetYaxis()->SetTitle("(#varepsilon #otimes A) / A");
    graphBE->GetXaxis()->SetTitle("m_{G} [GeV]");
    graphBE->Draw("PSE");
    fBE->Draw("same");

    legmc->Draw("same");

  }
  




}
std::map<std::string, std::vector<eff_reco> > reduce(std::map<std::string, std::vector<eff_reco> > effreco, std::string cut){

  std::map<std::string, std::vector<eff_reco> > out; 
  
  for(unsigned int iE =0; iE < effreco["001"].size(); iE++){
    if (effreco["001"][iE].sel == cut){ out["001"].push_back(effreco["001"][iE]); }
  }

  for(unsigned int iE =0; iE < effreco["01"].size(); iE++){
    if (effreco["01"][iE].sel == cut){ out["01"].push_back(effreco["01"][iE]); }
  }

  for(unsigned int iE =0; iE < effreco["02"].size(); iE++){
    if (effreco["02"][iE].sel == cut){ out["02"].push_back(effreco["02"][iE]); }
  }

  return out; 
}

//-----------------------------------------------------------------------------------
TH1F * createHisto(TChain * newtree1, const std::string &histo_name, int nBins, double xMin, double xMax, std::string cut){
   
  // std::cout << "Making histograms with cut\n" << cut << std::endl;
  TH1F * hist = new TH1F(histo_name.c_str(), histo_name.c_str(), nBins, xMin, xMax);

  if (cut == "no_selection"){
    newtree1->Project(histo_name.c_str(), "Diphoton.Minv");
  } else {
    newtree1->Project(histo_name.c_str(), "Diphoton.Minv", cut.c_str());
  }
    // move overflow to last bin
  int nbins = hist->GetNbinsX();
  float lastBin = hist->GetBinContent(nbins);
  float overflow = hist->GetBinContent(nbins+1);
  hist->SetBinContent(nbins, lastBin + overflow);
  hist->SetBinContent(nbins+1, 0.0);

  return hist;

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

