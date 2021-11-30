#ifndef AUXILIARYPARAFORLINK_H
#define AUXILIARYPARAFORLINK_H

#include "motionparameters.h"
#include "shapeparameters.h"

using namespace std;

class AuxiliaryParaForLink
{
public:
    AuxiliaryParaForLink();
    void initialize(vector<int> data_rows_cols_slices, long time_points, int cellNumInAllFrames, MotionParameters xRpMtPara, MotionParameters xStaticPara, ShapeParameters apprPara);
    // size
    size_t nnT;
    vector<int> nnYXZ;

    // parameters for jump
    size_t maxJump=5;
    // false postive rate
    float FPdtctRate=10.0e-3;
    // link threshold pdf
    double pLinkThr;
    // number of detections in all frames
    size_t nCellRep;

    // source & sink
    size_t source;
    size_t sink;


};

#endif // AUXILIARYPARAFORLINK_H
