#ifndef MICTRACKERMAIN_H
#define MICTRACKERMAIN_H

#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include "cellsegmentation/synquant_simple.h"




class MicTrackerMain
{
public:
    MicTrackerMain(void *data_grayim4d, int _data_type, long bufSize[5]/*(x,y,z,c,t)*/);
    void processSingleFrameAndReturn(int curr_timePoint_in_canvas,const QString &fileName);
    void cellSegmentSingleFrame(Mat *data_grayim3d, size_t curr_frame);
    void regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat,
                              Mat *idMapIn /*int*/, int seed_num,
                              int curr_frame);
    int getMaxContrast(Mat *data_grayim3d);
    //void cropSeed(vector<size_t> idx_yxz, Mat *data_3d, int shift_yxz[3]);
    void cropSeed(int seed_id, vector<size_t> idx_yxz, Mat *data_grayim3d,
                                  Mat *idMap, int curr_frame, singleCellSeed &seed,
                                  segParameter p4segVol);
    // the key function to get cell territory
    void refineSeed2Region(singleCellSeed &seed, odStatsParameter p4odStats, segParameter p4segVol);
    void retrieve_seeds(Mat *dataVolFloat, Mat *label_map_1stRound, size_t cell_num_1stRound,
                        Mat *cellGapMap, Mat &idMap_2ndRound, int &seed_num_2ndRound);
    bool saveSegResults(const QString &fileName);


    Mat * extract3d(Mat *data_grayim4d,int curr_frame);
    Mat * extract2d(Mat *data_grayim3d,int curr_slice,int datatype);
    void showSlice(Mat *data_grayim3d,int curr_slice,int datatype);
    size_t sub2idx(int x,int y,int z,int y_size, int xy_size);

public:
    string debug_folder;
    float scale_cell_seed_shift_bndtest[6] = {6, 2, 3, 3, 2, 1.5};
    string default_name;
    int data_type;
    Mat1b normalized_data4d;
    vector<int> data_rows_cols_slices;
    vector<bool> time_points_processed;
    long time_points = 0;
    long curr_time_point;
    segParameter p4segVol;
    odStatsParameter p4odStats;
    std::vector<cv::Mat> cell_label_maps;
    std::vector<cv::Mat> threshold_maps;
    std::vector<cv::Mat> principalCurv2d;
    std::vector<cv::Mat> principalCurv3d;
    std::vector<cv::Mat> varMaps;
    std::vector<cv::Mat> stblizedVarMaps;
    std::vector<std::vector<float>> varTrends;
    std::vector<std::vector<float>> stblizedVarTrends;
    std::vector<float> variances;
    std::vector<std::size_t> number_cells;
    //std::vector<std::vector<std::size_t>> voxIdxList; // cell voxIdx list

public:
    void init_parameter(){
        debug_folder = "/home/ccw/Desktop/embryo_res_folder/";
        default_name = debug_folder + "test.tiff";
        p4segVol.min_intensity = -1.0;
        p4segVol.fdr = 0.0;
        p4segVol.min_cell_sz = 50;
        p4segVol.max_cell_sz = 10000;
        p4segVol.min_fill = 0.0001;
        p4segVol.max_WHRatio = 100;
        p4segVol.min_seed_size = 10;
        p4segVol.graph_cost_design[0] = ARITHMETIC_AVERAGE; //default 1, GEOMETRIC_AVERAGE = 2;
        p4segVol.graph_cost_design[1] = 2;
        p4segVol.growConnectInTest = 4;
        p4segVol.growConnectInRefine = 6;
        p4segVol.edgeConnect = 48;
        p4segVol.neiMap = 26;
        p4segVol.connect4fgGapRemoval = 26;
        p4segVol.shift_yxz[0] = round(20 / scale_cell_seed_shift_bndtest[2]);
        p4segVol.shift_yxz[1] = round(20 / scale_cell_seed_shift_bndtest[3]);
        p4segVol.shift_yxz[2] = round(4 / scale_cell_seed_shift_bndtest[4]);
        p4segVol.shrink_flag = true;
        p4segVol.shrink_scale_yxz[0] = 2;
        p4segVol.shrink_scale_yxz[1] = 2;
        p4segVol.shrink_scale_yxz[2] = 1;
        p4segVol.fgBoundaryHandle = LEAVEALONEFIRST;
        p4segVol.gapTestMinMaxRadius[0] = 2;
        p4segVol.gapTestMinMaxRadius[1] = 4;
        p4segVol.growSeedInTracking = false;

        p4odStats.gap4varTrendEst = 2;
        p4odStats.gap4fgbgCompare = 0; // the fg and bg should be adjacent
        p4odStats.roundNum4fgbgCompare = round(3 / scale_cell_seed_shift_bndtest[5]); // fg/bg both are width-3 band
        p4odStats.varAtRatio = 0.80;
        p4odStats.fgSignificanceTestWay = KSEC;
        p4odStats.minGapWithOtherCell_yxz[0] = 3;
        p4odStats.minGapWithOtherCell_yxz[1] = 3;
        p4odStats.minGapWithOtherCell_yxz[2] = 1;
        p4odStats.connectInSeedRefine = 6;
        p4odStats.gapTestMethod = GAP_LOCALORDERSTATS;
        p4odStats.gapTestSkippedBandWidth = 2;
        p4odStats.gapTestThreshold = 0.01;
    }
    void reset_shift(){
        p4segVol.shift_yxz[0] = round(20 / scale_cell_seed_shift_bndtest[2]);
        p4segVol.shift_yxz[1] = round(20 / scale_cell_seed_shift_bndtest[3]);
        p4segVol.shift_yxz[2] = round(4 / scale_cell_seed_shift_bndtest[4]);
    }

};

#endif // MICTRACKERMAIN_H
