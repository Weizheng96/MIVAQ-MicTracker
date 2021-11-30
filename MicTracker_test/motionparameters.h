#ifndef MOTIONPARAMETERS_H
#define MOTIONPARAMETERS_H
#include "mictrackerrefinebytempcons.h"
#include <random>
#include <boost/math/special_functions/gamma.hpp>

class MotionParameters
{
public:
    MotionParameters();

    void Initialize(MicTrackerRefineByTempCons refineResult);
    MotionParameters static2motile(MotionParameters staticPara);

    vector<float> DsplcmtMu3d;
    vector<float> DsplcmtSigma3d;
    vector<float> distPhat;
    float invOvrlpRtPhat;
    float nbIncldRtPhat;

};

#endif // MOTIONPARAMETERS_H
