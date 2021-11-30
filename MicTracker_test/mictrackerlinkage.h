#ifndef MICTRACKERLINKAGE_H
#define MICTRACKERLINKAGE_H
#include "mictrackerrefinebytempcons.h"
#include "motionparameters.h"
#include "shapeparameters.h"
#include "auxiliaryparaforlink.h"
#include "NetworkFlowSolver.hpp"

class MicTrackerLinkage
{
public:
    ///////////////////////////////////////////////
    /// functions about initializtion
    ///////////////////////////////////////////////
    MicTrackerLinkage();
    MicTrackerLinkage(MicTrackerRefineByTempCons refineResult);
    void Initialize(MicTrackerRefineByTempCons refineResult);
    void setIndLUT();

public:
    ///////////////////////////////////////////////
    /// functions about MCF
    ///////////////////////////////////////////////
    void getIntegratedMotionModelMCFlinking();
    void getPosteriorProbability();
    void buildGraph();
    void getTraces();
    void getFullTrace(const auto &flowMap, const size_t &cur, const size_t &t, vector<QVector4D> &tempTrace);

public:
    ///////////////////////////////////////////////
    /// functions for test
    ///////////////////////////////////////////////
    void test();


public:
    ///////////////////////////////////////////////
    /// parameters from refine result
    ///////////////////////////////////////////////

    // label map
    std::vector<cv::Mat> cell_label_maps;

    // video size
    vector<int> data_rows_cols_slices;
    long time_points = 0;
    int cellNumInAllFrames=0;

    // current cell info
    vector<size_t> number_cells; // time -> object num
    vector<vector<vector<size_t>>> VoxelIdxList; // time -> cell -> pixel->idx
    vector<vector<float>> avgItstyVec; // time -> cell -> intensity
    vector<vector<size_t>> areaVec; // time -> cell -> intensity
    vector<vector<vector<float>>> ctrPt; // time -> cell -> Point Center (y,x,z)
    vector<vector<vector<size_t>>> CoordinateRgs; // time -> cell -> boundary (min y, max y, min x, max x, min z, max z);

    vector<vector<size_t>> localId2GlobalId;
    vector<vector<size_t>> globalId2LocalId;

    // estimated null hypothesis
    vector<float> nullMtScPhat;
    vector<float> altMtScPhat;

    float intDiff_mu;
    float intDiff_sigma;
    float sizeChange_mu;
    float sizeChange_sigma;
    float sizeChangeRatio_mu;
    float sizeChangeRatio_sigma;

    ShapeParameters apprPara;
    MotionParameters posParaStatic;

    // intensity difference
    vector<vector<vector<float>>> dtctMtScVec;

public:
    ///////////////////////////////////////////////
    /// parameters for MCF
    ///////////////////////////////////////////////

    // the motion of static cells
    MotionParameters xStaticPara;
    // the motion of motile cells
    MotionParameters xRpMtPara;
    // the probability of this is a motile cell based on intensity difference
    vector<vector<vector<float>>> rpdMtRltdDtctPostP1;
    // the p value of intensity difference by static cell null hypothesis
    vector<vector<vector<float>>> dtctPvalVec;
    // other auxiliary parameters
    AuxiliaryParaForLink para;

    // networkflow solver
    NetworkFlowSolver<double,int> mSolver;

    // nodes number
    size_t nNode=0;

public:
    ///////////////////////////////////////////////
    /// result for MCF
    ///////////////////////////////////////////////
    vector<vector<QVector4D>> traces;

};

#endif // MICTRACKERLINKAGE_H
