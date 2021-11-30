#include "motionparameters.h"

MotionParameters::MotionParameters()
{

}

void MotionParameters::Initialize(MicTrackerRefineByTempCons refineResult)
{
    DsplcmtMu3d.assign(refineResult.DsplcmtMu3d,refineResult.DsplcmtMu3d+3);
    DsplcmtSigma3d.assign(refineResult.DsplcmtSigma3d,refineResult.DsplcmtSigma3d+3);
    distPhat=refineResult.distPhat;
    invOvrlpRtPhat=refineResult.invOvrlpRtPhat;
    nbIncldRtPhat=refineResult.nbIncldRtPhat;

}



MotionParameters MotionParameters::static2motile(MotionParameters staticPara)
{
    MotionParameters motilePara(staticPara);

    for(auto& v:motilePara.DsplcmtSigma3d){v=v*pow(15,2);}

    int N=10000;
    default_random_engine generator;
    normal_distribution<float> distribution(0.0,1.0);

    float dist=0;
    vector<float> temp;
    temp.resize(N);


    for (int i=0; i<N; ++i) {
        dist=0;
        for(int j=0;j<3;++j){
            float num = distribution(generator);
            num=num*motilePara.DsplcmtSigma3d[j]+motilePara.DsplcmtMu3d[j];
            dist+=pow(num,2);
        }
        temp[i]=sqrt(dist);
    }

    MicTrackerRefineByTempCons::gamfit(temp, motilePara.distPhat);
    return motilePara;

}


