#include "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h"
#include "diphoton-analysis/RooUtils/interface/RooSlicePdf.h"
#include "diphoton-analysis/RooUtils/interface/RooStarMomentMorph.h"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"

namespace  {
    struct dictionary {
            RooPowLogPdf pl;
	    RooSlicePdf  sl;
	    RooStarMomentMorph  sm;
            RooDCBShape dcb;
    };
}

