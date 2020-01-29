#include "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h"
#include "diphoton-analysis/RooUtils/interface/RooSlicePdf.h"
#include "diphoton-analysis/RooUtils/interface/RooStarMomentMorph.h"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
#include "diphoton-analysis/RooUtils/interface/RooSpline1D.h"
#include "diphoton-analysis/RooUtils/interface/RooPowerLaw.h"
#include "diphoton-analysis/RooUtils/interface/RooPowerLawSum.h"

namespace  {
    struct dictionary {
            RooPowLogPdf pl;
	    RooSlicePdf  sl;
	    RooStarMomentMorph  sm;
            RooDCBShape dcb;
            RooPowerLaw plo;
            RooPowerLawSum plosu;
    };
}

