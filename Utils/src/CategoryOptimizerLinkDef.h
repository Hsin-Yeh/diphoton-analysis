#include "diphoton-analysis/Utils/interface/CategoryOptimizer.h"
#include "diphoton-analysis/Utils/interface/FunctionHelpers.h"
#include "diphoton-analysis/Utils/interface/NaiveCategoryOptimization.h"
#include "diphoton-analysis/Utils/interface/SimpleShapeCategoryOptimization.h"
#include "diphoton-analysis/Utils/interface/DataSetFiller.h"
#include "diphoton-analysis/Utils/interface/DataSetMixer.h"

#include "RVersion.h"

#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
#include "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h"
#include "diphoton-analysis/RooUtils/interface/RooSlicePdf.h"
#include "diphoton-analysis/RooUtils/interface/RooStarMomentMorph.h"

#pragma link C++ defined_in "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h";
#pragma link C++ defined_in "diphoton-analysis/RooUtils/interface/RooSlicePdf.h";
#pragma link C++ defined_in "diphoton-analysis/RooUtils/interface/RooStarMomentMorph.h";
#endif

#pragma link C++ defined_in "diphoton-analysis/Utils/interface/DataSetFiller.h";
#pragma link C++ defined_in "diphoton-analysis/Utils/interface/DataSetMixer.h";
#pragma link C++ defined_in "diphoton-analysis/Utils/interface/CategoryOptimizer.h";
#pragma link C++ defined_in "diphoton-analysis/Utils/interface/FunctionHelpers.h";
#pragma link C++ defined_in "diphoton-analysis/Utils/interface/NaiveCategoryOptimization.h";
#pragma link C++ defined_in "diphoton-analysis/Utils/interface/SimpleShapeCategoryOptimization.h";

#pragma link C++ class std::list<RooAbsData*>::iterator;

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

