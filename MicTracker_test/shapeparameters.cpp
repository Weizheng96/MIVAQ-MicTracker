#include "shapeparameters.h"

ShapeParameters::ShapeParameters()
{

}
void ShapeParameters::initialize(MicTrackerRefineByTempCons refineResult)
{
    intDiff_mu=refineResult.intDiff_mu;
    intDiff_sigma=refineResult.intDiff_sigma;
    sizeChange_mu=refineResult.sizeChange_mu;
    sizeChange_sigma=refineResult.sizeChange_sigma;
    sizeChangeRatio_mu=refineResult.sizeChangeRatio_mu;
    sizeChangeRatio_sigma=refineResult.sizeChangeRatio_sigma;
}
