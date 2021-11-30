#include "auxiliaryparaforlink.h"

AuxiliaryParaForLink::AuxiliaryParaForLink()
{

}

void AuxiliaryParaForLink::initialize(vector<int> data_rows_cols_slices,long time_points,int cellNumInAllFrames,MotionParameters xRpMtPara, MotionParameters xStaticPara, ShapeParameters apprPara)
{


    // size
    nnT=time_points;
    nnYXZ=data_rows_cols_slices;

    // center position change pdf
    double pCtrPos=gammapdf((double)xRpMtPara.distPhat[0]*(double)xRpMtPara.distPhat[1],
            (double)xRpMtPara.distPhat[0],(double)xRpMtPara.distPhat[1] );
    // intensity difference pdf
    double pIntst=normalPDF(sqrt((double)apprPara.intDiff_sigma)*2,(double)apprPara.intDiff_mu,(double)apprPara.intDiff_sigma);
    // size difference pdf
    double pSize=normalPDF(sqrt((double)apprPara.sizeChangeRatio_sigma)*2,(double)apprPara.sizeChangeRatio_mu,(double)apprPara.sizeChangeRatio_sigma);
    // overlap ratio pdf
    double pOvlpRt=exppdf((double)1, (double)xStaticPara.invOvrlpRtPhat);
    // neighbor include Ratio pdf
    double pNbIncldRt=exppdf((double) 1, (double)xStaticPara.nbIncldRtPhat);
    // jump time pdf
    double pT=exppdf((double)(1+maxJump),(double)(maxJump+1)/2);

    pLinkThr=min(pIntst*pSize*pCtrPos*pOvlpRt*pNbIncldRt*pT,(double)FPdtctRate/(1-FPdtctRate));

    // cell number in all frames
    nCellRep=cellNumInAllFrames;
}
