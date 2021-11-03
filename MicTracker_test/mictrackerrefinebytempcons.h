#ifndef MICTRACKERREFINEBYTEMPCONS_H
#define MICTRACKERREFINEBYTEMPCONS_H
#include "mictrackermain.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include "mincostflow.h"
#include <string>
#include <chrono> // time elapsed
#include <fstream> // for file stream
#include <QTextStream>
#include <QLabel>
#include <QWidget>
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
#include <chrono> // time elapsed
#include <fstream> // for file stream
#include <QTextStream>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <QDebug>
#include <algorithm>    // std::max



class MicTrackerRefineByTempCons
{
public:

    MicTrackerRefineByTempCons(void *data_grayim4d, int _data_type, long bufSize[5],
    const QStringList &filelist);

    MicTrackerRefineByTempCons(MicTrackerMain micPerFrame,const QStringList &filelist);

    ///////////////////////////////////////////////
    /// learn parameter
    /// ///////////////////////////////////////////
    void paraLearn();
    // check valid z slice
    void checkZStacks();
    // build null hypothesis
    void getDtctMtScVec();
    void getFeature4linking();

    float gamfitLkeqn(float a, float cnst);
    float fzero_gamfitLkeqn(float lower, float upper, float cnst);
    void gamfit(vector<float> x, vector<float> &parmhat);
    void histogram(vector<float> x);
    static void CorrespondingCell(vector<size_t> vxLst,vector<size_t> &candIds,vector<size_t> &ovlpAreas);
    static float avg ( vector<float> v );
    static float stdev ( vector<float> v  , float &mean );
    static void TruncGaussFit(vector<float> x, float &muhat, float &sigmahat);
    static float gradiantTruncSigma(vector<float> x, float mu, float sigma, float alpha, float beta);
    static float gradiantTruncMu(vector<float> x, float mu, float sigma, float alpha, float beta);
    static float normal_pdf(float x, float m, float s);
    static float standardNormal_cdf(float x);
    static size_t findVecBoundary(vector<float> x_sorted,float boundary);

    ///////////////////////////////////////////////
    /// temporal refine
    /// ///////////////////////////////////////////
    void segRefineByTempCons();
    void refineInitial();
    // estimate the cell number of each detection
    void tempConsistSort();
    // refine the detection based cell number

    // get the link between detection
    void incldRelationship();
    void incldRelationshipInitial();
    void childDtctSet(size_t t0,size_t iSubId0,bool flagPre);

    static void getCommonElement(vector<size_t> vector1,vector<size_t> vector2,vector<size_t> &result);
    void getCommonElementSorted(vector<size_t> &vector1, vector<size_t> &vector2, vector<size_t> &result);
    static void getCommonElementFast(vector<size_t> vector1,vector<size_t> vector2,vector<size_t> &result);

    void segRefineOperations();
    void updateRefineResult();
    void updateDtctMtScVec();



public:

    ///////////////////////////////////////////
    /// input from in frame detection
    ///////////////////////////////////////////
    // data
    MicTrackerMain micTracker;

    std::vector<cv::Mat> cell_label_maps;
    std::vector<cv::Mat> principalCurv3d;
    std::vector<cv::Mat> segScore3d;
    vector<int> data_rows_cols_slices;
    long time_points = 0;
    Mat1b normalized_data4d;


    // cell info
    std::vector<std::size_t> number_cells; // time -> object num
    std::vector<std::vector<std::vector<std::size_t>>> VoxelIdxList; // time -> cell -> pixel->idx
    std::vector<std::vector<float>> avgItstyVec; // time -> cell -> intensity
    std::vector<std::vector<std::size_t>> areaVec; // time -> cell -> intensity
    std::vector<std::vector<std::vector<float>>> ctrPt; // time -> cell -> Point Center (y,x,z)
    std::vector<std::vector<std::vector<size_t>>> CoordinateRgs; // time -> cell -> boundary (min y, max y, min x, max x, min z, max z);

    ///////////////////////////////////////////
    /// learn parameters
    ///////////////////////////////////////////
    // effective z slice
    Mat effZinFmTwoEndIdx;

    // motion score of detections
    std::vector<std::vector<std::vector<float>>> dtctMtScVec; // time -> cell -> preScore,PostScore

    // null probability distribution with empirical distribution
    vector<float> nullMtScPhat;
    vector<float> altMtScPhat;
    float intDiff_mu=0;
    float intDiff_sigma=0;
    float sizeChange_mu=0;
    float sizeChange_sigma=0;
    float sizeChangeRatio_mu=0;
    float sizeChangeRatio_sigma=0;
    float DsplcmtMu3d[3]={0,0,0};
    float DsplcmtSigma3d[3]={0,0,0};
    vector<float> distPhat={0,0};
    float invOvrlpRtPhat=0;
    float nbIncldRtPhat=0;

    ///////////////////////////////////////////
    /// temporal refine
    ///////////////////////////////////////////

    // current cell info
    std::vector<std::size_t> number_cells_cur; // time -> object num
    std::vector<std::vector<std::vector<std::size_t>>> VoxelIdxList_cur; // time -> cell -> pixel->idx
    std::vector<std::vector<float>> avgItstyVec_cur; // time -> cell -> intensity
    std::vector<std::vector<std::size_t>> areaVec_cur; // time -> cell -> intensity
    std::vector<std::vector<std::vector<float>>> ctrPt_cur; // time -> cell -> Point Center (y,x,z)
    std::vector<std::vector<std::vector<size_t>>> CoordinateRgs_cur; // time -> cell -> boundary (min y, max y, min x, max x, min z, max z);

    std::vector<std::vector<std::vector<size_t>>> dtctIncldLst_cur; // {pre,post} -> detection -> included cell->id
    std::vector<std::size_t> numCellResVec_cur; // detection -> estimated cell number
    std::vector<std::vector<std::size_t>> cellIdxsLst; // time -> detection -> overall id
    bool isChanged=true; // if the temporal refine give change
    size_t nDtct;
    // threshold
    float pMtScThre=0.005;

    // graph solver
    MinCostFlow minCostSolver;

    ///////////////////////////////////////////
    /// for debug
    ///////////////////////////////////////////

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    float timer=((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000;
    float timer_c=0;

};

#endif // MICTRACKERREFINEBYTEMPCONS_H
