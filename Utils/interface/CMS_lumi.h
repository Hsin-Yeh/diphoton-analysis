#include <iostream>
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
// /#include "TASImage.h"

//
// Global variables
//

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

TString magText     = "";

bool disableExtraText = false;
bool writeExtraText = true;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.58;//THIS IS WHAT DEFINES the CMS text and the luminosity
float cmsTextSize      = 0.50;
float cmsTextOffset    = -0.5;  // only used in outOfFrame version
float extraTextOffset   = 0.05;

float relPosX    = 0.045;
float relPosY    = 0.;//0.035
float relExtraDY = -0.03;//AND THIS IS WHAT DEFINES THE preliminary position. 
float relExtraDX = -1.5;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "41.527 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";

bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10 );

