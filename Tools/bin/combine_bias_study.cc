// Bias Stydy: More info here: 
// https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/625.html

#include<string> 
#include<iostream>
#include<fstream>
#include<sstream>
#include <sys/stat.h>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

struct theFitResult {
  double meanFitE;
  double rmsFitE;
  double meanFitEerr;
  double rmsFitEerr;
 
};

//-----------------------------------------------------------------------------------
//Declaration here definition after main
theFitResult fitpulls(TH1 * histo, std::string outputdir, std::string scenario, std::string scenarioforplot);
std::string widthtonum(std::string coupling);
bool exists (const std::string& name);
double computeHistFHWM(TH1*  hist);

//-----------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

  std::string year, inputdir, outputdir;

  if(argc!=4) {
    std::cout << "Syntax: combine_bias_study.exe [2016/2017/2018] [input] [output] " << std::endl;
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
  }

  
  //========================================================================
  //========================================================================
  //ZERO PART : We give whatever input/output variables needed here.  
  //========================================================================
  //========================================================================
  //========================================================================

  //The mode of the task. It is referring to either generating and fitting
  //with the same function or with a different one. 
  std::string combmode = "diffunfit"; //"diffunfit"samefunfit
  //The type of signal we are running
  std::string insigname = "grav";
  //Injected signal strength
  std::vector<std::string> muinjected;
  muinjected.push_back("1");
  muinjected.push_back("2");
  muinjected.push_back("3");
  //muinjected.push_back("4");
  //muinjected.push_back("5");
  
  //Couplings for the moment
  std::vector<std::string> coups;
  coups.clear();
  if (insigname == "grav"){
    coups.push_back("kMpl001");
    coups.push_back("kMpl01");
    coups.push_back("kMpl02");
  } else if (insigname == "heavyhiggs"){
    coups.push_back("0p014");
    coups.push_back("1p4");
    coups.push_back("5p6");
  } else {
    std::cout << "The signal you are looking for doesn't exist. " << std::endl;
    exit(1);
  }

  //Models
  std::vector<std::string> models;
  // 2017: Chosen ones
  if ( year == "2017"){

    // models.push_back("pow"); 
    // models.push_back("Laurent");
    // models.push_back("Exponential"); 
    // models.push_back("VVdijet"); 
    models.push_back("expow1");
    models.push_back("invpow1");
    models.push_back("invpowlin1");
    // models.push_back("moddijet");

    //Add the dijet only in case of samefunfit
    if ( combmode == "samefunfit" ) {models.push_back("dijet");} 
  }

  // 2018: Chosen ones
  if ( year == "2018"){
    models.push_back("Laurent"); 
    models.push_back("PowerLaw"); 
    //models.push_back("Exponential");  
    models.push_back("VVdijet"); 
    //Add the dijet only in case of samefunfit
    if ( combmode == "samefunfit" ) {models.push_back("dijet");} 
  }
  
  //Masses run in this bias study
  std::vector<std::string> masses;
  masses.push_back("600");
  masses.push_back("700");
  masses.push_back("800");
  masses.push_back("900");
  masses.push_back("1000");
  masses.push_back("1100");
  masses.push_back("1200");
  masses.push_back("1500");
  masses.push_back("1800");
  masses.push_back("2100");
  masses.push_back("2400");
  masses.push_back("2700");
  masses.push_back("3000");
  masses.push_back("3500");
  masses.push_back("4000");
  masses.push_back("4500");
  masses.push_back("5000");
  masses.push_back("5500");
  masses.push_back("6000");

  std::vector<std::string> cats;
  cats->clear();
  cats.push_back("EBEE");
  cats.push_back("EBEB");

  //We will make a string to hold the file under investigation
  std::string scenario = "";
  std::string scenarioforplot = "";

  std::map<std::string , std::map<std::string , std::map<TString , TGraphErrors *> > > bprofiles; //[model][coup][muin] 
  std::map<std::string , std::map<std::string, std::map<std::string , std::vector<double> > > > bpullvals, bpullvalsErr, massvals, massvalsErr;
  
  for (auto model : models){
    for (auto cp : coups){
      for (auto muin : muinjected){
        for (auto mass : masses){
          for (auto cat : cats){

            std::cout << "MODEL " << model << " COUP " << cp << " MU " << muin << " MASS " << mass << std::endl;
            scenario = "tree_" + combmode + "_" + insigname + "_mu" + muin + "_" + cp + "_" + cat + "-" model + "_mass" + mass;
            scenarioforplot =  cp + " " + cat + " " + model + " " + mass + " GeV";

            //-----------------------------------------------------------------------------------
            //Read the file and the tree
            TString filename =  inputdir + "/" + scenario  + ".root";
            std::string filename_s = inputdir + "/" + scenario  + ".root";
            TString trname = "tree_fit_sb";
            //TString trname = "tree_fit_b";

            //Variables in the trees we want to plot
            Int_t fit_status;
            double  mu;
            double  muErr;
            double  muLoErr;
            double  muHiErr;

            //Getting trees from files
            std::cout << "Getting tree from file " << filename << std::endl;
            if ( !exists(filename_s) ){continue;}

            TFile* infile = TFile::Open(filename);
            //The tree in the file that we want to get
            TTree *tr = (TTree*) infile->Get(trname);

            tr->SetBranchAddress("fit_status", &fit_status);
            tr->SetBranchAddress("r", &mu);
            tr->SetBranchAddress("rErr", &muErr);
            tr->SetBranchAddress("rLoErr", &muLoErr);
            tr->SetBranchAddress("rHiErr", &muHiErr);

            //-----------------------------------------------------------------------------------
            //Make the histos from the tree

            double pullxlow = combmode == "samefunfit" ? -3.5 : -10.;
            double pullylow = combmode == "samefunfit" ? 3.5: 10.;
            //Histos
            TH1F * h_mu = new TH1F( "h_mu", "Signal strength from fit", 100, 0., 3.);
            TH1F * h_muErr = new TH1F( "h_muErr", "Signal strength error from fit", 100, 0., 16.);
            TH1F * h_muLoErr = new TH1F( "h_muLoErr", "Signal strength low error from fit", 100, 0., 3.);
            TH1F * h_muHiErr = new TH1F( "h_muHiErr", "Signal strength high error from fit", 100, 0., 16.);
            TH1F * h_pulls_mu = new TH1F("h_pulls_mu", "mu_{true} - mu_{fit} / err", 50, pullxlow, pullylow);

            std::cout << "Running bias study " << std::endl;
 
            //For the pull value
            double pull_mu = 0.;
            double err_mu = 0.;
            float muinj = std::stof(muin);

            std::cout << " Starting to loop over total entries of tree " << tr->GetEntries() << std::endl;
            for (Int_t i=0; i<tr->GetEntries(); i++) {
              tr->GetEntry(i);

              // if (fit_status==0){
              if (fit_status==0 && ((muHiErr+mu) < 19) && ((mu-muLoErr)>-19) ){
                // if ((muHiErr+mu) < 19 && (mu-muLoErr)>-19 ){
                // if (true){
                // std::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT " << fit_status << std::endl;
                if (mu > muinj){
                  // err_mu = muLoErr;
                  err_mu = muHiErr;
                  // err_mu = muErr;
                } else{
                  err_mu = muHiErr;
                }

                //HACK HERE
                err_mu = 0.5 * (muLoErr+muHiErr);
                // if (err_mu ==0){
                // 	pull_mu = 0.;
                // } else{
                // 	pull_mu = (muinj - mu )/err_mu;
                // }
                pull_mu = (muinj - mu )/err_mu;
                // if(fabs(pull_mu) > 10. ){
                // 	std::cout <<  "muinj= " << muinj << std::endl;
                // 	std::cout <<  "mu_fit - mu_inj  = " << (mu - muinj)<< std::endl;
                // 	std::cout <<  "err_mu = " << err_mu << std::endl;
                // 	std::cout <<  "pull_mu = " << pull_mu << std::endl;
                // }

                // if(fabs(pull_mu) < 0.01){
                // 	std::cout <<  "muinj= " << muinj << std::endl;
                // 	std::cout <<  "mu_fit - mu_inj  = " << (mu - muinj) << std::endl;
                // 	std::cout <<  "err_mu = " << err_mu << std::endl;
                // 	std::cout <<  "pull_mu = " << pull_mu << std::endl;

                // }

                h_mu->Fill(mu);
                h_muErr->Fill(muErr);
                h_muLoErr->Fill(muLoErr);
                h_muHiErr->Fill(muHiErr);

                h_pulls_mu->Fill( pull_mu );

              }

            } // End of loop over entries

            //-----------------------------------------------------------------------------------
            //Plot and save the final pull
            theFitResult tmpfitpull = fitpulls(h_pulls_mu, outputdir + "/" + combmode, scenario, scenarioforplot);

            bpullvals[model][cp][muin].push_back(tmpfitpull.meanFitE);
            bpullvalsErr[model][cp][muin].push_back(tmpfitpull.meanFitEerr);
            massvals[model][cp][muin].push_back( std::stod(mass) );
            massvalsErr[model][cp][muin].push_back( 0. );

            infile->Close();
          }
        }//end of loop over masses
      }//end of loop over muin
    }//end of loop over couplings
  }//end of loop over models

  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptStat(1111);
  m_gStyle->SetOptFit(1);

  std::map<TString , int > colors;
  colors["pow"] = 2;
  colors["Laurent"] = 3;
  colors["PowerLaw"] = 4;
  colors["Atlas"] = 5;
  colors["VVdijet"] = 6;
  colors["Exponential"] = 7;
  colors["invpowlin"] = 8;
  colors["expow"] = 9;
  colors["invpow"] = 28;
  colors["invpowlin"] = 34;
  colors["moddijet"] = 46;
  colors["dijet"] = 1;

  // colors["Laurent"] = 616; //kMagenta
  // colors["PowerLaw"] = kGray; //kGray
  // colors["Atlas"] = 433; //kCyan+1
  // colors["Expow"] = 2; //kRed
  // colors["invpow"] = kOrange; //kGray
  // colors["invpowlin"] = 417; //kGreen+1
  // colors["Exponential"] = 4; //kBlue
  // colors["invpow"] = kYellow; //kYellow
  // colors["dijet"] = kBlack; //kYellow
  // colors["moddijet"] = ; //kGray
  // colors["pow"] = ; //kOrange
  
  //For bias
  std::map<std::string , std::map<std::string , TCanvas *> > bcanv; //[cat]
  std::map<std::string , std::map<std::string ,TLegend *> > bleg; //[cat]

  std::map<std::string , std::map<std::string , TH2F *> > frame;
  std::map<std::string , std::map<std::string , TBox *> > box;

  float xfirst = 500.;
  float xlast  = 6500.;
  float yfirst = -3.;
  float ylast = 2.;
  
  for (auto cp : coups){
    std::string couplabel = widthtonum(cp);      
    for (auto muin : muinjected){
      std::string muinlabel = "#mu_{true} = " + muin ;

      std::string blabel = cp + "_" + muin;

      bcanv[cp][muin] = new TCanvas(Form("profile_pull_%s_%s",cp.c_str(), muin.c_str()), Form("profile_pull_%s_%s",cp.c_str(), muin.c_str()) );
      bcanv[cp][muin]->SetLogx();
      bcanv[cp][muin]->SetGridy();
      bcanv[cp][muin]->SetGridx();

      bleg[cp][muin]= new TLegend(0.6,0.12,0.95,0.47);
      bleg[cp][muin]->SetFillStyle(0);
      bleg[cp][muin]->SetBorderSize(0);

      frame[cp][muin] = new TH2F( Form("frame_%s_%s", cp.c_str(), muin.c_str() ), Form("frame_%s_%s", cp.c_str(), muin.c_str() ),100,xfirst,xlast,100,yfirst,ylast);
      frame[cp][muin]->SetStats(false);
      bcanv[cp][muin]->cd();
      frame[cp][muin]->Draw();
      frame[cp][muin]->SetTitle("");
      frame[cp][muin]->GetXaxis()->SetTitle("mass [GeV]");
      frame[cp][muin]->GetXaxis()->SetMoreLogLabels();
      frame[cp][muin]->GetYaxis()->SetTitle("( #mu_{fit} - #mu_{true} )/ #sigma_{fit}");

      box[cp][muin]= new TBox(xfirst,-0.5,xlast,0.5);
      box[cp][muin]->SetFillColor(kGray);
      box[cp][muin]->Draw("same");       

      //------------------------------------------------------------------------------------------------------
      //First with the pull bias plot
      for (auto model : models){
        // A std::vector is at its heart an array. To get the array just get the address of the first element.
        bprofiles[model][cp][muin] = new TGraphErrors( bpullvals[model][cp][muin].size(), &massvals[model][cp][muin][0], &bpullvals[model][cp][muin][0], &massvalsErr[model][cp][muin][0], &bpullvalsErr[model][cp][muin][0] ) ;

        bprofiles[model][cp][muin]->GetXaxis()->SetRangeUser(xfirst,xlast);
        bprofiles[model][cp][muin]->GetXaxis()->SetTitleOffset( 0.9 );
        bprofiles[model][cp][muin]->SetMarkerColor(colors[model]);
        bprofiles[model][cp][muin]->SetMarkerStyle(kFullCircle);
        bleg[cp][muin]->AddEntry( bprofiles[model][cp][muin] , Form("%s",model.c_str()) ,"pe" );
        // bleg[cp][muin]->SetHeader(Form("dijet %s",cat.c_str()));
        bprofiles[model][cp][muin]->Draw("PSE");
        bleg[cp][muin]->Draw("same");
        bcanv[cp][muin]->RedrawAxis();
        bcanv[cp][muin]->Modified();
        bcanv[cp][muin]->Update();

      }//end of loop through models

      char buf[500];
      TLatex lat;
      double xmi = xfirst; double xma = xlast; 
      double latx = xmi+(xma-xmi)/20.;
      double laty = ylast;

      sprintf(buf,"%s        %s   ", couplabel.c_str(), muinlabel.c_str() );
      lat.DrawLatex(latx,laty*0.85,buf);
      // sprintf(buf,"%s", muinlabel.c_str() );
      // lat.DrawLatex(latx,laty*0.6,buf);

      
      bcanv[cp][muin]->SaveAs((Form("%s/%s/profile_pull_%s_%s.png", outputdir.c_str(), combmode.c_str(), cp.c_str(), muin.c_str() )) );
    }//end of loop over muin
  }//end of loop over couplings
 
  
}

//======================================================================================
theFitResult fitpulls(TH1 * histo, std::string outputdir, std::string scenario, std::string scenarioforplot){

  theFitResult thefitres;
  
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("pulls","pulls",600,600);
  // c1->Divide(2,3);
  char buf[500];
  TLatex lat;
  lat.SetTextColor(kRed);

  // histo->SetMaximum(70);
  histo->Draw("pe");

  double meanE = histo->GetMean();
  double rmsE = histo->GetRMS();
  double meanEerr = histo->GetMeanError();
  double rmsEerr = histo->GetRMSError();
  // double fwhm = computeHistFHWM(histo);
  int bin1 = histo->FindFirstBinAbove(histo->GetMaximum()/2);
  int bin2 = histo->FindLastBinAbove(histo->GetMaximum()/2);
  double fwhm = histo->GetBinCenter(bin2) - histo->GetBinCenter(bin1);
  
  std::cout << " --- " << histo->GetName() << " = entries " << histo->GetEntries()
            << " nbins " << histo->GetNbinsX()
            << " " << histo->GetBinLowEdge(1)
            << " " << histo->GetBinLowEdge(histo->GetNbinsX())
            << " mean " << meanE
            << " rms " << rmsE
            << " MeanError " << meanEerr<< " RMSError " <<  rmsEerr
            << " underflows " << histo->GetBinContent(0)
            << " overflows " << histo->GetBinContent(histo->GetNbinsX()+1)
            << " FWHM " << fwhm
            << std::endl;

  //fit
  // double fitregionlow = ;
  histo->Fit("gaus","LR0","", //gauss
             meanE-2*rmsE, //meanE-2*rmsE,2*fwhm
             meanE+2*rmsE); //meanE+2*rmsE2*fwhm
  // -2., //meanE-2*rmsE,2*fwhm
  // 2.); //meanE+2*rmsE2*fwhm

  TF1 *fitResult = histo->GetFunction("gaus");//gaus

  fitResult->SetLineColor(2);

  double meanFitE = fitResult->GetParameter(1);
  double rmsFitE = fitResult->GetParameter(2);                                                                                  
  double meanFitEerr = fitResult->GetParError(1);
  double rmsFitEerr = fitResult->GetParError(2);    

  thefitres.meanFitE = meanFitE;
  thefitres.rmsFitE  = rmsFitE;
  thefitres.meanFitEerr = meanFitEerr;
  thefitres.rmsFitEerr = rmsFitEerr;
  
  fitResult->Draw("same");
    
  std::cout << " Gauss fit: mean " << meanFitE << " RMS " << rmsFitE
            << " Gauss fit: mean error " <<  meanFitEerr << " RMS error" << rmsFitEerr
            << std::endl;

  double xmi = 0.; double xma = 0.; 
  xmi = histo->GetXaxis()->GetXmin();
  xma = histo->GetXaxis()->GetXmax();

  double latx = xmi+(xma-xmi)/20.;
  double laty = histo->GetMaximum();
  lat.DrawLatex(latx,laty*0.9,scenarioforplot.c_str());
  sprintf(buf,"<meanfit> = %3.3f +/- %3.3f",fitResult->GetParameter(1),fitResult->GetParError(1));
  lat.DrawLatex(latx,laty*0.8,buf);
  sprintf(buf,"RMSfit = %3.3f +/- %3.3f",fitResult->GetParameter(2),fitResult->GetParError(2));
  lat.DrawLatex(latx,laty*0.7,buf);
  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  lat.DrawLatex(latx,laty*0.6,buf);
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  lat.DrawLatex(latx,laty*0.5,buf);

  TString savename = outputdir + "/" + scenario;
  c1->Print(savename+"_pull.pdf");
  c1->Print(savename+"_pull.png");

  return thefitres; 

 
}

//-----------------------------------------------------------------------------------
std::string widthtonum(std::string coupling){
  std::string coup = ""; 
 
  if ( coupling == "0p014" ){coup = "1.4 \\times 10^{-4}";}
  else if ( coupling == "1p4" ){coup = "1.4 \\times 10^{-2}";}
  else if ( coupling == "5p6" ){coup = "5.6 \\times 10^{-2}";}
  else if ( coupling == "kMpl001" ){coup = " #tilde{k} = #frac{k}{M_{pl}} = 0.01 ";}
  else if ( coupling == "kMpl01" ){ coup = " #tilde{k} = #frac{k}{M_{pl}} = 0.1";}
  else if ( coupling == "kMpl02" ){ coup = " #tilde{k} = #frac{k}{M_{pl}} = 0.2 ";}

  return coup;


}
//-----------------------------------------------------------------------------------
bool exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
//-----------------------------------------------------------------------------------
double computeHistFHWM(TH1*  hist){

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

  return xWidth;
}
