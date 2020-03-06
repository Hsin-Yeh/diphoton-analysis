#include "diphoton-analysis/RooUtils/interface/RooPowLogPdf.h"
#include "diphoton-analysis/RooUtils/interface/RooSlicePdf.h"
#include "diphoton-analysis/RooUtils/interface/RooDCBShape.h"
#include "diphoton-analysis/RooUtils/interface/RooPowerLawSum.h"

namespace  {
    struct dictiona {
            RooPowLogPdf pl;
	    RooSlicePdf  sl;
            RooDCBShape dcb;
            RooPowerLawSum plosu;
    };
}

