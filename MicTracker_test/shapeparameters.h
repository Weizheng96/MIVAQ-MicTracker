#ifndef SHAPEPARAMETERS_H
#define SHAPEPARAMETERS_H

#include "mictrackerrefinebytempcons.h"

class ShapeParameters
{
public:
    ShapeParameters();
    void initialize(MicTrackerRefineByTempCons refineResult);

    float intDiff_mu;
    float intDiff_sigma;
    float sizeChange_mu;
    float sizeChange_sigma;
    float sizeChangeRatio_mu;
    float sizeChangeRatio_sigma;
};

#endif // SHAPEPARAMETERS_H
