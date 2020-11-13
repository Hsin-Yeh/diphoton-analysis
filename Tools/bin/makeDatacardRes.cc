// creates combine datacard for the resonant analysis
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "TH1.h"
#include "TFile.h"

#include "diphoton-analysis/Tools/interface/sampleList.hh" // integrated luminosity only defined once
#include "diphoton-analysis/Tools/interface/utilities.hh"

class nuisance {

public:
  nuisance(std::string syst, std::string distribution, std::vector<std::string> contribution) :
    m_syst(syst),
    m_distribution(distribution),
    m_contribution(contribution) {}
  std::string m_syst;
  std::string m_distribution;
  std::vector<std::string> m_contribution;
};

double getYield(const std::string& region, const std::string& sample, const std::string& datacardYear, double& yieldError, const TF1 * scaleFactor = 0);
std::string getDiphotonYieldVariations(const std::string& region, const std::string& variation);
// void makeDatacard( const std::string& region, const std::string& datacardYear);
void makeDatacard( const std::string& model, const std::string& coup, const std::string& insigname, std::vector<std::string> cats, const std::string& datacardYear);
std::map<std::string, std::vector<int> > buildmodelorder (const std::string& datacardYear);

std::string datacardYear, insigname;

int main(int argc, char *argv[])
{

  if(argc!=3) {
    std::cout << "Syntax: makeDatacardRes.exe [2016/2017/2018] [insigname] \n";
    return -1;
  }
  datacardYear = argv[1];
  insigname = argv[2];

  //Categories
  std::vector<std::string> cats; 
  cats.clear(); 
  cats.push_back("EBEB");
  cats.push_back("EBEE");

  //Couplings for the moment
  std::vector<std::string> coups; coups.clear();
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
  // 2016: Chosen ones
  //EBEB: {"pow","Laurent","PowerLaw","Exponential","dijet"}
  //EBEE: {"pow","Laurent","Exponential","VVdijet","dijet"}
  //Below put the greatest one of the two above
  std::vector<std::string> exceptions;
  exceptions.clear();

  if ( datacardYear == "2016"){
    // models.push_back("Laurent"); 
    // models.push_back("pow"); 
    // models.push_back("PowerLaw");
    // models.push_back("Exponential"); 
    // models.push_back("VVdijet"); 
    models.push_back("dijet");

    //Here we put the exceptions of the models that do not describe the data in certain
    //categories. 
    exceptions.push_back("PowerLaw_EBEE");
    exceptions.push_back("VVdijet_EBEB");
  }

  // 2017: Chosen ones
  //EBEB: {"pow","Laurent","PowerLaw","Exponential","dijet"}
  //EBEE: {"Laurent","Exponential","VVdijet","expow","dijet"}
  //Below put the greatest one of the two above
  if ( datacardYear == "2017"){
    models.push_back("Laurent"); 
    // models.push_back("pow"); 
    // models.push_back("PowerLaw"); 
    models.push_back("Exponential"); 
    // models.push_back("VVdijet"); 
    // models.push_back("expow"); 
    models.push_back("dijet");

    //Here we put the exceptions of the models that do not describe the data in certain
    //categories. 
    exceptions.push_back("pow_EBEE");
    exceptions.push_back("PowerLaw_EBEE");
    exceptions.push_back("VVdijet_EBEB");
    exceptions.push_back("expow_EBEB");
  }

  //EBEB: {"pow","Laurent","PowerLaw","Exponential","VVdijet","dijet"}
  //EBEE: {"Laurent","PowerLaw","Exponential","VVdijet","invpowlin","dijet"}
  //Below put the greatest one of the two above
  if ( datacardYear == "2018"){
    models.push_back("Laurent"); 
    //models.push_back("pow"); 
    models.push_back("PowerLaw");
    models.push_back("Exponential"); 
    models.push_back("VVdijet"); 
    //models.push_back("invpowlin"); 
    //models.push_back("dijet");

    //Here we put the exceptions of the models that do not describe the data in certain
    //categories. 
    exceptions.push_back("pow_EBEE");
    exceptions.push_back("invpowlin_EBEB");
  }

  //We will process one year each time
  // makeDatacard( "EBEB", datacardYear);
  // makeDatacard( "EBEE", datacardYear);

  for (auto model : models){
    for (auto cp : coups){
      
      std::cout << "MODEL " << model << std::endl;
      makeDatacard(model, cp, insigname, cats, datacardYear);

    }//end of loop over couplings
  }//end of loop over models
    

}

// void makeDatacard(const std::string& region, const std::string& datacardYear)
void makeDatacard(const std::string& model, const std::string& coup, const std::string& insigname, std::vector<std::string> cats, const std::string& datacardYear)
{

  //Systematics
  std::vector<nuisance> allNuisances;  
  std::string lumiError = std::to_string(1 + luminosityErrorFrac[datacardYear]);
  //We have two processes, sig and bkg. bkg is estimated from data so 
  //a dash is there. 
  nuisance lumi("lumi", "lnN", {lumiError, "-", lumiError, "-"});
  nuisance eff("eff", "lnN", {"1.06", "-", "1.06", "-"});
  nuisance PDFs("PDFs", "lnN", {"1.06", "-", "1.06", "-"});
  nuisance thetaSmear("thetaSmear$CHANNEL_13TeV", "param", {"0.0", "1.0"});
  nuisance thetaScale("thetaScale$CHANNEL_13TeV", "param", {"0.0", "1.0"});
  nuisance dijetflatparam1("PhotonsMass_bkg_dijet_linc_cat0", "flatParam", {" "," "}); 
  nuisance dijetflatparam2("PhotonsMass_bkg_dijet_logc_cat0", "flatParam", {" "," "}); 
  nuisance dijetflatparam3("PhotonsMass_bkg_dijet_linc_cat1", "flatParam", {" "," "}); 
  nuisance dijetflatparam4("PhotonsMass_bkg_dijet_logc_cat1", "flatParam", {" "," "}); 

  allNuisances.push_back(lumi);
  allNuisances.push_back(eff);
  allNuisances.push_back(PDFs);
  allNuisances.push_back(thetaSmear);
  allNuisances.push_back(thetaScale);
  allNuisances.push_back(dijetflatparam1);
  allNuisances.push_back(dijetflatparam2);
  allNuisances.push_back(dijetflatparam3);
  allNuisances.push_back(dijetflatparam4);

  std::ofstream output;
  std::string filename("datacards/"+ insigname + "_" + coup + "_" + model+ "_");
  filename+=datacardYear;
  // filename+="_" + region;
  filename+=".dat";
  output.open(filename);

  //Order
  std::map<std::string, std::vector<int> > modelorder; //[model][order] for cats
  modelorder.clear();
  //Order of the truth/alternative models from BkgFit.cc study.
  modelorder = buildmodelorder(datacardYear);
  std::cout << modelorder["Laurent"].size() << std::endl;

  std::map<std::string,std::string> modelnametopdf;
  modelnametopdf["pow"] = "pow";
  modelnametopdf["Laurent"] = "lau";
  modelnametopdf["PowerLaw"] = "pow";
  modelnametopdf["Atlas"] = "atlas";
  modelnametopdf["VVdijet"] = "vvdijet";
  modelnametopdf["Exponential"] = "exp";
  modelnametopdf["invpowlin"] = "invpowlin";
  modelnametopdf["expow"] = "expow";
  modelnametopdf["invpow"] = "invpow";
  modelnametopdf["invpowlin"] = "invpowlin";
  modelnametopdf["moddijet"] = "moddijet";
  // modelnametopdf[""] = ;
    
  std::vector<std::string> thepdfname;
  thepdfname.clear();

  if (model!="Laurent" && model!="PowerLaw" && model!="Atlas" && model!="Exponential" && 
      model!="Chebychev" && model!="DijetSimple" && model!="Dijet" && model!="VVdijet" &&
      model!="Expow" && model!="dijet"){
    thepdfname.push_back( Form("PhotonsMassBkg_%s%d_EBEB", model.c_str(), modelorder[model][0]) );
    thepdfname.push_back( Form("PhotonsMassBkg_%s%d_EBEE", model.c_str(), modelorder[model][1]) );
  } else if ( model == "dijet"  ){
    thepdfname.push_back( "PhotonsMassBkg_" + model + "_EBEB" );
    thepdfname.push_back( "PhotonsMassBkg_" + model + "_EBEE" );
  } else {
    std::cout << Form("ftest_pdf_EBEB_%s" , modelnametopdf[model].c_str()) << std::endl;
    std::cout << modelnametopdf[model] << std::endl;
    std::cout << modelorder[model].size() << std::endl;
    std::cout << modelorder[model][0]<< std::endl;
    std::cout << Form("ftest_pdf_EBEB_%s%d" , modelnametopdf[model].c_str(), modelorder[model][0]) << std::endl;
    thepdfname.push_back( Form("ftest_pdf_EBEB_%s%d" , modelnametopdf[model].c_str(), modelorder[model][0]) );
    thepdfname.push_back( Form("ftest_pdf_EBEE_%s%d" , modelnametopdf[model].c_str(), modelorder[model][1]) );
  }

  if (output.is_open()) {

    // header
    output << "imax * number of channels" << std::endl;
    output << "jmax * number of backgrounds" << std::endl;
    output << "kmax * number of nuisance parameters" << std::endl;
    output << "-------------------" << std::endl;

    //output << "shapes data_obs * data_" << datacardYear << ".root " << " data_ws:Data_$CHANNEL" << std::endl;
    //output << "shapes bkg * bkg_dijet_001_" << datacardYear << ".root " << " bkg_ws:model_dijet_001_$CHANNEL" << std::endl;
    //output << "shapes sig * SignalParametricShapes_ws.root " << " model_signal:SignalShape_kMpl001_$CHANNEL" << endl;

    //data shape
    output << "shapes data_obs * bkg_" << model << "_" << "$CHANNEL_" << datacardYear << ".root " << " HLFactory_" << model << "_ws:Data_$CHANNEL" << std::endl;
    //bkg shape
    // output << "shapes bkg * bkg_" << model << "_" << "$CHANNEL_" << datacardYear << ".root " << " HLFactory_" << model << "_ws:" << thepdfname << std::endl;
    output << "shapes bkg EBEB bkg_" << model << "_" << "EBEB_" << datacardYear << ".root " << " HLFactory_" << model << "_ws:" << thepdfname[0] << std::endl;
     output << "shapes bkg EBEE bkg_" << model << "_" << "EBEE_" << datacardYear << ".root " << " HLFactory_" << model << "_ws:" << thepdfname[1] << std::endl;
   //signal shape
    output << "shapes sig * " << insigname << "_" << coup << "_" << datacardYear << ".root " << " wtemplates:model_signal_" << insigname << "_" << coup  << "_$CHANNEL" << std::endl; 
    output << "------------------------------" << std::endl;
    
    output << "bin          " <<  cats[0] <<  " " << cats[1] << std::endl;
    output << "observation   -1    -1"  <<  std::endl;
    output << "------------------------------" << std::endl;
    output << "bin          " <<  cats[0] << " " << cats[0] <<  " " << cats[1] <<  " " << cats[1] << std::endl;
    output << "process                 sig      bkg    sig      bkg" << std::endl;
    output << "process                   0        1      0        1" << std::endl;
    output << "rate                     1       1    1      1 " << std::endl;

    output << "\n-------------------" << std::endl;
    
    // output systematics
    for(auto inuisance : allNuisances) {
      output << inuisance.m_syst  << "      " << inuisance.m_distribution << "  ";
      for( auto icontrib : inuisance.m_contribution) {
	output << icontrib << " ";
      }
      output << "\n";
    }

  } // closes is_open()
  else {
    std::cout << "Could not open output file " << filename << std::endl;
  }
}



std::map<std::string, std::vector<int> > buildmodelorder (const std::string& datacardYear)
{

  std::map<std::string, std::vector<int> > modelorder;
  std::cout << datacardYear << std::endl;
  //-------------------------------------------------------------------------------
  if ( datacardYear == "2016"){

    modelorder["pow"].push_back(1);//2
    modelorder["pow"].push_back(1);//4

    modelorder["Laurent"].push_back(1);//2
    modelorder["Laurent"].push_back(1);//3

    modelorder["PowerLaw"].push_back(1);//1 
    modelorder["PowerLaw"].push_back(1);//1 No good fit

    modelorder["Atlas"].push_back(1);//1 No good fit
    modelorder["Atlas"].push_back(1);//1 No good fit

    modelorder["Exponential"].push_back(1);//3
    modelorder["Exponential"].push_back(1);//3

    // modelorder["Chebychev"].push_back(3);
    // modelorder["Chebychev"].push_back(3);

    modelorder["DijetSimple"].push_back(1);//2
    modelorder["DijetSimple"].push_back(1);//2

    modelorder["Dijet"].push_back(1);//2
    modelorder["Dijet"].push_back(1);//2

    modelorder["VVdijet"].push_back(1);//2
    modelorder["VVdijet"].push_back(1);//2

    modelorder["expow"].push_back(1);//1 No good fit
    modelorder["expow"].push_back(1);//1

    // modelorder["Expow"].push_back(2);//1
    // modelorder["Expow"].push_back(3);//1 No good fit

    modelorder["invpow"].push_back(1);//2 
    modelorder["invpow"].push_back(1);//2

    modelorder["invpowlin"].push_back(1);//2 No good fit
    modelorder["invpowlin"].push_back(1);//2 

    modelorder["moddijet"].push_back(1); //2
    modelorder["moddijet"].push_back(1); //2

    modelorder["dijet"].push_back(1);
    modelorder["dijet"].push_back(1);
  }

  //-------------------------------------------------------------------------------
  if ( datacardYear == "2017"){

    modelorder["pow"].push_back(1);//2
    modelorder["pow"].push_back(1);//4

    modelorder["Laurent"].push_back(1);//2
    modelorder["Laurent"].push_back(2);//3

    modelorder["PowerLaw"].push_back(1);//1 
    modelorder["PowerLaw"].push_back(1);//1 No good fit

    modelorder["Atlas"].push_back(1);//1 No good fit
    modelorder["Atlas"].push_back(1);//1 No good fit

    modelorder["Exponential"].push_back(1);//3
    modelorder["Exponential"].push_back(1);//3

    // modelorder["Chebychev"].push_back(3);
    // modelorder["Chebychev"].push_back(3);

    modelorder["DijetSimple"].push_back(1);//2
    modelorder["DijetSimple"].push_back(1);//2

    modelorder["Dijet"].push_back(1);//2
    modelorder["Dijet"].push_back(1);//2

    modelorder["VVdijet"].push_back(1);//2
    modelorder["VVdijet"].push_back(1);//2

    modelorder["expow"].push_back(1);//1 No good fit
    modelorder["expow"].push_back(1);//1

    // modelorder["Expow"].push_back(2);//1
    // modelorder["Expow"].push_back(3);//1 No good fit

    modelorder["invpow"].push_back(1);//2 
    modelorder["invpow"].push_back(1);//2

    modelorder["invpowlin"].push_back(1);//2 No good fit
    modelorder["invpowlin"].push_back(1);//2 

    modelorder["moddijet"].push_back(1); //2
    modelorder["moddijet"].push_back(1); //2

    modelorder["dijet"].push_back(1);
    modelorder["dijet"].push_back(1);
  }

  //-------------------------------------------------------------------------------
  if ( datacardYear == "2018"){

    std::cout << " 222222 "<< std::endl; 

    modelorder["pow"].push_back(1);//2
    modelorder["pow"].push_back(1);//4

    modelorder["Laurent"].push_back(1);//2
    modelorder["Laurent"].push_back(2);//3

    modelorder["PowerLaw"].push_back(1);//1 
    modelorder["PowerLaw"].push_back(1);//1 No good fit

    modelorder["Atlas"].push_back(1);//1 No good fit
    modelorder["Atlas"].push_back(1);//1 No good fit

    modelorder["Exponential"].push_back(1);//3
    modelorder["Exponential"].push_back(1);//3

    // modelorder["Chebychev"].push_back(3);
    // modelorder["Chebychev"].push_back(3);

    modelorder["DijetSimple"].push_back(1);//2
    modelorder["DijetSimple"].push_back(1);//2

    modelorder["Dijet"].push_back(1);//2
    modelorder["Dijet"].push_back(1);//2

    modelorder["VVdijet"].push_back(1);//2
    modelorder["VVdijet"].push_back(1);//2

    modelorder["expow"].push_back(1);//1 No good fit
    modelorder["expow"].push_back(1);//1

    // modelorder["Expow"].push_back(2);//1
    // modelorder["Expow"].push_back(3);//1 No good fit

    modelorder["invpow"].push_back(1);//2 
    modelorder["invpow"].push_back(1);//2

    modelorder["invpowlin"].push_back(1);//2 No good fit
    modelorder["invpowlin"].push_back(1);//2 

    modelorder["moddijet"].push_back(1); //2
    modelorder["moddijet"].push_back(1); //2

    modelorder["dijet"].push_back(1);
    modelorder["dijet"].push_back(1);
  }

  std::cout << modelorder["Laurent"].size() << std::endl;
  return modelorder;
}
