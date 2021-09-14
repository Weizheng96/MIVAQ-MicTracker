#include "mictrackermain.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"
#include <string>
#include <chrono> // time elapsed
#include <fstream> // for file stream
#include <QTextStream>
#include <QLabel>
#include <QWidget>

#include <QDebug>

using namespace cv;
using namespace std;

MicTrackerMain::MicTrackerMain(void *data_grayim4d, int _data_type, long bufSize[5]/*(x,y,z,c,t)*/)
{
    init_parameter();
    data_rows_cols_slices.resize(3);

    if (bufSize[0] > INT_MAX || bufSize[1] > INT_MAX || bufSize[2] > INT_MAX){
        qDebug("Data is too large to process: %d!", INT_MAX);
    }
    data_rows_cols_slices[0] = bufSize[1]; // !!! V3D transpose the image
    data_rows_cols_slices[1] = bufSize[0];
    data_rows_cols_slices[2] = bufSize[2];
    if (bufSize[3] != 1){
        //qFatal("Input is not a gray image\n");
        Q_ASSERT(bufSize[3] != 1);
    }
    time_points = bufSize[4];
    time_points_processed.resize(time_points);
    fill(time_points_processed.begin(), time_points_processed.end(), false);
    data_type = _data_type;
    cell_label_maps.resize(time_points);
    threshold_maps.resize(time_points);
    principalCurv2d.resize(time_points);
    principalCurv3d.resize(time_points);
    varMaps.resize(time_points);
    stblizedVarMaps.resize(time_points);
    varTrends.resize(time_points);
    stblizedVarTrends.resize(time_points);
    variances.resize(time_points);
    number_cells.resize(time_points);

    Mat *data4d;
    int data_sz[4] = {data_rows_cols_slices[0], data_rows_cols_slices[1],
                      data_rows_cols_slices[2], (int)time_points};
    if (data_type == V3D_UINT16) {
        data4d = new Mat(4, data_sz, CV_16U, data_grayim4d);
    }
    else if(data_type == V3D_UINT8){
        data4d = new Mat(4, data_sz, CV_8U, data_grayim4d);
    }else{
        qFatal("Unsupported data type\n");
    }
    //normalize(*data4d, *data4d, 0, 255, NORM_MINMAX, CV_8U);

//    Mat tmp;
//    data4d->copyTo(tmp);
//    qInfo("%d - %d - %d - %d \n", data4d->size[0], data4d->size[1],
//            data4d->size[2], data4d->size[3]);
//    normalize(tmp, normalized_data4d, 255, 0, NORM_MINMAX, CV_8U);
    //Mat normalized_data4d;
//    normalize(*data4d, normalized_data4d, 255, 0, NORM_MINMAX, CV_8U);
    data4d->copyTo(normalized_data4d);
    //ccShowSlice3Dmat(data4d, CV_16U);
    assert(normalized_data4d.type() == CV_8U);
    delete data4d;
}
void MicTrackerMain::processSingleFrameAndReturn(int curr_timePoint_in_canvas){
    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];
    curr_time_point = curr_timePoint_in_canvas;

    unsigned char *ind = (unsigned char*)normalized_data4d.data + sz_single_frame*curr_time_point; // sub-matrix pointer
    Mat *single_frame = new Mat(3, normalized_data4d.size, CV_8U, ind);

    cell_label_maps[curr_time_point] = Mat::zeros(3, normalized_data4d.size, CV_32S); // int label
    threshold_maps[curr_time_point] = Mat::zeros(3, normalized_data4d.size, CV_8U);
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    cellSegmentSingleFrame(single_frame, curr_time_point);
    time_points_processed[curr_time_point] = true;
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    qInfo("----------------time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);

    if(number_cells[curr_time_point] > 50000){
        qDebug("!!!!!!!!!!!!!:We should not detect that many cells !");
    }
    qInfo("----------------totally: %ld cells are detected", number_cells[curr_time_point]);
    //ccShowSliceLabelMat(cell_label_maps[curr_time_point]);
    delete single_frame;
}

void MicTrackerMain::cellSegmentSingleFrame(Mat *data_grayim3d, size_t curr_frame)
{
    //data_grayim3d is uint8 0-255 datatype
    Mat *dataVolFloat = new Mat(data_grayim3d->dims, data_grayim3d->size, CV_32F);
    data_grayim3d->convertTo(*dataVolFloat, CV_32F);
    int thres=getMaxContrast(data_grayim3d);
//    Mat tmp = *dataVolFloat>25.5;
//    ccShowSlice3Dmat(tmp, CV_8U,20);
//    ccShowSlice3Dmat(data_grayim3d, CV_8U,20);
    /******** start to do cell segmentation *******/
//    float sigma2d[3] = {3.0, 3.0, 0.0};
//    principalCv2d(dataVolFloat, principalCurv2d[curr_frame], sigma2d, p4segVol.min_intensity);
    //ccShowSlice3Dmat(&principalCurv2d[curr_frame], CV_32F, 3);
    float sigma3d[3] = {3.0, 3.0, 3.0};
    principalCv3d(dataVolFloat, principalCurv3d[curr_frame], sigma3d, true,p4segVol.min_intensity);
//    Mat tmp = principalCurv3d[curr_frame]<=0;
//    ccShowSlice3Dmat(tmp, CV_8U,20);

    // synquant based on principalCurv3d


    Mat imgIn_temp, imgIn_temp2;
    sqrt(imgIn_temp,imgIn_temp2);
    Mat imgPCfloat=1-imgIn_temp2;
    Mat imgPCfloat255,imgPC8U;
    imgPCfloat255=imgPCfloat*255;

    variances[curr_frame] = calVarianceStablization(&imgPCfloat255, varMaps[curr_frame], varTrends[curr_frame],
                                                   p4odStats.varAtRatio, p4odStats.gap4varTrendEst);
    //ccShowSlice3Dmat(dataVolFloat, CV_32F, 3);
    ccShowSlice3Dmat(&varMaps[curr_frame], CV_32F, 3);

    //////////////////////////////////////////////////////////////////
    //           1. use synQuant for PC score
    //////////////////////////////////////////////////////////////////

    imgPCfloat255.convertTo(imgPC8U, CV_8U);
    synQuantSimple seeds_from_synQuant(&imgPC8U, variances[curr_frame], p4segVol, p4odStats);
//    for(int i=53;i<66;i++){
//        ccShowSliceLabelMat(seeds_from_synQuant.idMap,i);
//    }


    //////////////////////////////////////////////////////////////////
    //                 2. refine the seed regions                   //
    //////////////////////////////////////////////////////////////////
    regionWiseAnalysis4d(&imgPC8U, &imgPCfloat255,seeds_from_synQuant.idMap,
                         seeds_from_synQuant.cell_num,curr_frame);


}

int MicTrackerMain::getMaxContrast(Mat *data_grayim3d){
    Mat roiMap;
    Mat roiOuterEdgeMap;
    Mat dataFloat;
    Mat dataFloat2d;
    vector<vector<size_t>> voxList_ROI;
    vector<vector<size_t>> voxList_edge;
    float contrastVec[254]={0};

    data_grayim3d->convertTo(dataFloat,CV_32F);

    float iROI, iEdge;
    int nROI, nEdge;
    bool increaseFlag=false;
    int sz_single_frame=dataFloat.size[0]*dataFloat.size[1];
    Mat element=getStructuringElement(MORPH_ELLIPSE,Size(3,3));

    for(int thre=1;thre<255;thre++){

        iROI=0.0;
        iEdge=0.0;
        nROI=0;
        nEdge=0;

        float *ind = (float*)dataFloat.data; // sub-matrix pointer
        dataFloat2d= Mat(2, dataFloat.size, CV_32F, ind);
        for(int z=1;z<dataFloat.size[2];z++){
            float *ind = (float*)dataFloat.data + sz_single_frame*z; // sub-matrix pointer
            Mat temp= Mat(2, dataFloat.size, CV_32F, ind);
            max(temp,dataFloat2d,dataFloat2d);
        }


        roiMap=dataFloat2d>thre;
        dilate(roiMap,roiOuterEdgeMap,element);
        roiOuterEdgeMap=roiOuterEdgeMap-roiMap;
        for (size_t i = 0; i < roiMap.total(); i++){
            if(roiMap.at<unsigned char>(i)>0){
                iROI+=dataFloat2d.at<float>(i);
                nROI++;
            }else if(roiOuterEdgeMap.at<unsigned char>(i)>0){
                iEdge+=dataFloat2d.at<float>(i);
                nEdge++;
            }
        }

        if(nROI==0){
            break;
        }
        contrastVec[thre-1]=(iROI/nROI)-(iEdge/nEdge);
    }

    float contrastVecSmooth[254]={0};
    for(int thre=2;thre<254;thre++){
        contrastVecSmooth[thre-1]=0;
        for(int i=-1;i<2;i++){
            contrastVecSmooth[thre-1]+=contrastVec[thre-1+i];
        }

        if(contrastVecSmooth[thre-1]-contrastVecSmooth[thre-2]>0){
            increaseFlag=true;
        }else{
            if(increaseFlag==true){
                return thre;
            }else{
                increaseFlag=false;
            }
        }
    }
    return 0;
}

/**
 * @brief retrieve_seeds: retrieve the seeds of cells lost in the first round of synQuant
 * @param smoothed_VolFloat: 0-255 but in float format
 * @param label_map_1stRound
 * @param cell_num_1stRound
 * @param test_ids
 */
void MicTrackerMain::retrieve_seeds(Mat *smoothed_VolFloat, Mat *label_map_1stRound, size_t cell_num_1stRound,
                                     Mat *cellGapMap, Mat &idMap_2ndRound, int &seed_num_2ndRound){
    if(cell_num_1stRound == 0){
        idMap_2ndRound = Mat::zeros(label_map_1stRound->dims, label_map_1stRound->size, CV_32S);
        seed_num_2ndRound = 0;
        return;
    }
    // 1. extract the cell size and intensity
    vector<vector<size_t>> voxList_exist_cells;
    extractVoxIdxList(label_map_1stRound, voxList_exist_cells, cell_num_1stRound);
    vector<float> cell_mean_intensity;
    regionAvgIntensity(smoothed_VolFloat, voxList_exist_cells, cell_mean_intensity);
    vector<size_t> exist_cell_size (voxList_exist_cells.size());
    FOREACH_i(voxList_exist_cells) exist_cell_size[i] = voxList_exist_cells[i].size();

    size_t min_exist_cell_size = *min_element(exist_cell_size.begin(), exist_cell_size.end());
    size_t max_exist_cell_size = *max_element(exist_cell_size.begin(), exist_cell_size.end());
    float min_cell_intensity = *min_element(cell_mean_intensity.begin(), cell_mean_intensity.end());

    // 2. extract the new seed map if it satisifies the constraints from existing cells
    Mat seed_map;
    Mat exist_cell_territory, exist_cell_fg;
    exist_cell_fg = *label_map_1stRound > 0;
    int dilate_radius[] = {2, 2, 1};
    volumeDilate(&exist_cell_fg, exist_cell_territory, dilate_radius, MORPH_ELLIPSE);
    //setValMat(*smoothed_VolFloat, CV_32F, &exist_cell_territory, 0);
    //ccShowSlice3Dmat(smoothed_VolFloat, CV_32F);
    //volumeWrite(smoothed_VolFloat, "/home/ccw/Desktop/remain.tif");
    //ccShowSlice3Dmat(cellGapMap, CV_8U);
    //volumeWrite(&exist_cell_territory, "/home/ccw/Desktop/cell_territory.tif");
    bitwise_or(*cellGapMap, exist_cell_territory, exist_cell_territory);
    //ccShowSlice3Dmat(exist_cell_territory, CV_8U);
    //setValMat(*smoothed_VolFloat, CV_32F, &exist_cell_territory, 0);
    //volumeWrite(smoothed_VolFloat, "/home/ccw/Desktop/remain_further.tif");
    //ccShowSlice3Dmat(smoothed_VolFloat, CV_32F);
    //volumeWrite(cellGapMap, "/home/ccw/Desktop/gap.tif");
    //volumeWrite(&exist_cell_territory, "/home/ccw/Desktop/overlap.tif");
    seed_map = exist_cell_territory == 0;
    dilate_radius[0] = 1; dilate_radius[1] = 1; dilate_radius[2] = 0;
    Mat seed_map_dilate;
    volumeDilate(&seed_map, seed_map_dilate, dilate_radius, MORPH_ELLIPSE);
    Mat label_map;
    int numCC;
    numCC = connectedComponents3d(&seed_map_dilate, label_map, 26);
    setValMat(label_map, CV_32S, &exist_cell_territory, 0);
    //ccShowSliceLabelMat(label_map);
    vector<vector<size_t>> voxList_seeds;
    extractVoxIdxList(&label_map, voxList_seeds, numCC);
    vector<float> seeds_mean_intensity;
    regionAvgIntensity(smoothed_VolFloat, voxList_seeds, seeds_mean_intensity);
    idMap_2ndRound = Mat::zeros(label_map.dims, label_map.size, CV_32S);
    seed_num_2ndRound = 0;
    FOREACH_i(voxList_seeds){
        if(voxList_seeds[i].size() >= min_exist_cell_size &&
                voxList_seeds[i].size() <= max_exist_cell_size &&
                seeds_mean_intensity[i] > min_cell_intensity){
            seed_num_2ndRound ++;
            setValMat(idMap_2ndRound, CV_32S, voxList_seeds[i], (float)seed_num_2ndRound);
        }
    }
    //ccShowSliceLabelMat(idMap_2ndRound);
}

void MicTrackerMain::boundaryTouchedTest(Mat *label_map, Mat *fgMap, bool &xy_touched, bool &z_touched){
    int n_y[] = { -1, -1, -1,  1, 1, 1,  0, 0 };// 8 shifts to neighbors
    int n_x[] = { -1,  0,  1, -1, 0, 1, -1, 1 };// used in functions
    //n_z = {  0,  0,  0,  0, 0, 0,  0, 0 };
    int z_id, y, x;
    size_t remain;
    size_t page_sz = label_map->size[0] * label_map->size[1];
    FOREACH_i_ptrMAT(label_map){
        if(label_map->at<int>(i) > 0){
            z_id = i / page_sz;
            if(z_id == 0 || z_id == (label_map->size[2]-1)){
                z_touched = true;
                break;
            }
        }
    }
    int im_sz[] = {label_map->size[0], label_map->size[1]};
    FOREACH_i_ptrMAT(label_map){
        if(label_map->at<int>(i) > 0){
            z_id = i / page_sz;
            remain = i - z_id * page_sz;
            y = remain / label_map->size[1];
            x = remain - y * label_map->size[1];
            for(int i=0; i < 8; i++){
                if(inField(y + n_y[i], x + n_x[i], im_sz)
                        && fgMap->at<unsigned char>(y + n_y[i], x + n_x[i], z_id) == 0){
                    xy_touched = true;
                    return;
                }
            }
        }
    }
}

//void cellSegmentMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat, Mat * volStblizedFloat, Mat *idMap /*int*/, int seed_num, Mat *eigMap2d,
//                                           Mat *eigMap3d, Mat *varMap, Mat *stblizedVarMap, vector<int> test_ids){
void MicTrackerMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat,
                                          Mat *idMapIn /*int*/, int seed_num,int curr_frame){
    Mat idMap;
    idMapIn->copyTo(idMap);
    //1. sort the seeds based on intensity levels
//    vector<float> seed_intensity;
//    regionAvgIntensity(dataVolFloat, &idMap, seed_num, seed_intensity);
    //ccShowSlice3Dmat(dataVolFloat, CV_32F);
//    vector<size_t> seed_intensity_order;
//    seed_intensity_order = sort_indexes(seed_intensity, false); // false->descending
    //2. for each seed, refine it region
    vector<vector<size_t>> voxIdxList(seed_num);
    extractVoxIdxList(&idMap, voxIdxList, seed_num);

    int cell_cnt = 0;
    //int debug_cell_id = 2;
    FOREACH_i(voxIdxList){
        int seed_id = i + 1;
        qInfo("------------------#%ld, start %d seed process, size: %ld---------------------", i, seed_id,
              voxIdxList[seed_id-1].size());

        singleCellSeed seed;
        cropSeed(seed_id, voxIdxList[seed_id-1], data_grayim3d, &idMap,
                curr_frame, seed, p4segVol);
        refineSeed2Region(seed, p4odStats, p4segVol);


        qInfo("----------------%d cell output from this seed-------------------", seed.outCell_num);
        if(seed.outCell_num <= 0) continue;
        subVolReplace(cell_label_maps[curr_time_point], CV_32S, seed.outputIdMap, seed.crop_range_yxz, cell_cnt);
        cell_cnt += seed.outCell_num;
    }
    for(int z=0;z<100;z++){
        ccShowSliceLabelMat(cell_label_maps[curr_time_point],z);
    }
//    removeSmallCC(cell_label_maps[curr_time_point], cell_cnt, p4segVol.min_cell_sz, true);
    number_cells[curr_time_point] = cell_cnt;
}

void MicTrackerMain::cropSeed(int seed_id, vector<size_t> idx_yxz, Mat *data_grayim3d,
                              Mat *idMap, int curr_frame, singleCellSeed &seed,
                              segParameter p4segVol){
    seed.id = seed_id;
    seed.idx_yxz = idx_yxz; // deep copy
    seed.y.resize(idx_yxz.size());
    seed.x.resize(idx_yxz.size());
    seed.z.resize(idx_yxz.size());
    vec_ind2sub(idx_yxz, seed.y, seed.x, seed.z, idMap->size);// MatSize itself is a vector
//    ccShowSliceLabelMat(idMap);
    getRange(seed.y, p4segVol.shift_yxz[0], idMap->size[0], seed.crop_range_yxz[0]);
    getRange(seed.x, p4segVol.shift_yxz[1], idMap->size[1], seed.crop_range_yxz[1]);
    getRange(seed.z, p4segVol.shift_yxz[2], idMap->size[2], seed.crop_range_yxz[2]);


    subVolExtract(idMap, CV_32S, seed.idMap, seed.crop_range_yxz);
    seed.seedMap = seed.idMap == seed.id;

    subVolExtract(&principalCurv3d[curr_frame], CV_32F, seed.eigMap3d, seed.crop_range_yxz);// deep copy
    seed.gap3dMap = seed.eigMap3d > 0;
    subVolExtract(data_grayim3d, CV_8U, seed.volUint8, seed.crop_range_yxz);// deep copy
    seed.outputIdMap = Mat(seed.idMap.dims, seed.idMap.size, CV_32S, Scalar(0));
    //ccShowSliceLabelMat(idMap);
    vec_sub2ind(seed.idx_yxz_cropped, vec_Minus(seed.y, seed.crop_range_yxz[0].start),
            vec_Minus(seed.x, seed.crop_range_yxz[1].start),
            vec_Minus(seed.z, seed.crop_range_yxz[2].start), seed.idMap.size);// MatSize itself is a vector
    //ccShowSliceLabelMat(idMap);
    seed.otherIdMap = seed.idMap > 0;
    FOREACH_i(seed.idx_yxz_cropped){
        seed.otherIdMap.at<unsigned char>(seed.idx_yxz_cropped[i]) = 0;
    }
    //ccShowSliceLabelMat(idMap);
    //p4segVol.shift_yxz[2] = 0;
    //// NOTE: if there is a segmentation error, check first if the Mat format is wrongly
    //// assigned. This is a problem of using opencv, where data format is predefined everywhere.
}

void MicTrackerMain::refineSeed2Region(singleCellSeed &seed, odStatsParameter p4odStats, segParameter p4segVol){

//    double tmp_min, tmp_max;

    Mat PCseed_all= (seed.eigMap3d <= 0.0);
    PCseed_all.convertTo(PCseed_all,CV_8U);


    Mat PCseed_valid;
    bitwise_and(PCseed_all,seed.seedMap,PCseed_valid);

    // get each PC seed
    Mat PCseed_valid_CC;
    int PCseedNum=connectedComponents3d(&PCseed_valid,PCseed_valid_CC,26);

    Mat seed_valid_CC;
    seed.seedMap.copyTo(seed_valid_CC);
    seed_valid_CC/=255;

    removeSmallCC(PCseed_valid_CC, PCseedNum, p4segVol.min_seed_size, true);


//    minMaxIdx(PCseed_valid_CC, &tmp_min, &tmp_max);


    if(PCseedNum<0){
        seed.outCell_num=0;
        PCseed_valid_CC.convertTo(seed.outputIdMap,CV_32S);
    }else if(PCseedNum<2){
        seed.outCell_num=1;
        seed_valid_CC.convertTo(seed.outputIdMap,CV_32S);
    }else{
        // Find total markers
        vector<vector<Point> > contours;
        Mat img4watershed=255-seed.volUint8;
        Mat output;
        PCseed_valid_CC.convertTo(output,CV_32S);

        Mat in[] = {img4watershed, img4watershed, img4watershed};
        merge(in, 3, img4watershed);

        watershed(img4watershed, output);
        output.convertTo(seed.outputIdMap,CV_32S);

//        double tmp_min2, tmp_max2;
//        minMaxIdx(output, &tmp_min2, &tmp_max2);
//        if(tmp_max!=tmp_max2){

//            for(int z=0;z<(seed.crop_range_yxz[2].end-seed.crop_range_yxz[2].start);z++){
//                ccShowSliceLabelMat(PCseed_valid_CC,z);
//            }

//            for(int z=0;z<(seed.crop_range_yxz[2].end-seed.crop_range_yxz[2].start);z++){
//                ccShowSliceLabelMat(seed.outputIdMap,z);
//            }
//        }

    }

//    for(int z=0;z<(seed.crop_range_yxz[2].end-seed.crop_range_yxz[2].start);z++){
//        ccShowSliceLabelMat(seed.outputIdMap,z);
//    }


}



/**
 * @brief removeOtherSeedsInfgMap: if we use thresholding, the fgMap may contains more than one seed. This will
 * not happen if we consider othercellterritory as in cellTerritoryExtractFromSeed();
 * @param seed
 * @param p4segVol
 */
void MicTrackerMain::removeOtherSeedsInfgMap(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    //// 1. find if there is other cells in fgMap
    vector<size_t> fg_idx = fgMapIdx(&cellSegFromSynQuant.fgMap, CV_8U, 0);
    bool other_cell_exist = false;
    FOREACH_i(fg_idx){
        if(seed.idMap.at<int>(fg_idx[i]) > 0 &&
                seed.idMap.at<int>(fg_idx[i]) != seed.id){
            other_cell_exist = true;
            break;
        }
    }
    if(!other_cell_exist){
        return;
    }
    Mat seedMap = Mat::zeros(cellSegFromSynQuant.fgMap.dims, cellSegFromSynQuant.fgMap.size, CV_32S);
    FOREACH_i(fg_idx){
        if(seed.idMap.at<int>(fg_idx[i]) > 0){
            if(seed.idMap.at<int>(fg_idx[i]) == seed.id){
                seedMap.at<int>(fg_idx[i]) = 1;
            }else{
                seedMap.at<int>(fg_idx[i]) = 2;
            }
        }else if(isOnBoundary2d(&seed.validSearchAreaMap, fg_idx[i])){
            // if the boundary of init foreground is also contained, label as sink
            // this is only considered when we double seed.shift_yxz, since the area are
            // too large now. The boudnary of init foreground should not be part of target
            // cell. Before doubling seed.shift_yxz, it is possible.
            seedMap.at<int>(fg_idx[i]) = 2;
        }
    }
    Mat grownSeedMap2d, grownSeedMap3d;
    bool bg2sink = true;
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    //ccShowSliceLabelMat(seedMap);
    regionGrow(&seedMap, 2, grownSeedMap2d, &seed.score2d, &cellSegFromSynQuant.fgMap,
               p4segVol.growConnectInTest, p4segVol.graph_cost_design, bg2sink);
    //ccShowSliceLabelMat(grownSeedMap2d);
    bg2sink = false;
    regionGrow(&grownSeedMap2d, 2, grownSeedMap3d, &seed.scoreMap, &cellSegFromSynQuant.fgMap,
               p4segVol.growConnectInRefine, p4segVol.graph_cost_design, bg2sink);
    //ccShowSliceLabelMat(grownSeedMap3d);
    cellSegFromSynQuant.fgMap = grownSeedMap3d == 1;
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
}
void MicTrackerMain::cellShrinkTest(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    int extra_cell = 0;
    for(int i = 1; i <= cellSegFromSynQuant.cell_num; i++){
        Mat cur_cell = seed.outputIdMap == i;
        Mat shrinked_cell, label_shrinked_cell;
        volumeErode(&cur_cell, shrinked_cell, p4segVol.shrink_scale_yxz, MORPH_ELLIPSE);
//        if(seed.crop_range_yxz[0].end == 94 && seed.crop_range_yxz[1].end == 123){
//            ccShowSlice3Dmat(&shrinked_cell, CV_8U);
//        }
        int n = connectedComponents3d(&shrinked_cell, label_shrinked_cell, p4segVol.growConnectInRefine);

        if(n > 1){
            removeSmallCC(label_shrinked_cell, n, p4segVol.min_seed_size, true);
            if(n > 1){
                Mat grown_shrinked_cells;
                bool link_bg2sink = false;
                regionGrow(&label_shrinked_cell, n, grown_shrinked_cells, &seed.scoreMap, &cur_cell,
                           p4segVol.growConnectInRefine, p4segVol.graph_cost_design, link_bg2sink);
                setValMat(seed.outputIdMap, CV_32S, &cur_cell, 0.0);
                Mat sub_cell_mask = grown_shrinked_cells == 1;
                setValMat(seed.outputIdMap, CV_32S, &sub_cell_mask, (float)i);

                for(int j=2; j<=n; j++){
                    sub_cell_mask = grown_shrinked_cells == j;
                    extra_cell ++;
                    setValMat(seed.outputIdMap, CV_32S, &sub_cell_mask, (float)(extra_cell + cellSegFromSynQuant.cell_num));
                }
            }
        }
    }
    cellSegFromSynQuant.cell_num += extra_cell;
}

/**
 * @brief refineFgWithSeedRegion, extract the largest region if fg and remove those
 * voxels that has no z-direction neighbors
 * @param seed
 * @param p4segVol
 */
void MicTrackerMain::refineCellTerritoryWithSeedRegion(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    Mat labelMap;
    int n = connectedComponents3d(&cellSegFromSynQuant.fgMap, labelMap, 6);
    //ccShowSlice3Dmat(seed.seedMap, CV_8U);

    if (n > 1){
        int id = largestRegionIdExtract(&labelMap, n, &seed.seedMap);
        cellSegFromSynQuant.fgMap = labelMap == id;
        //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    }
    // remove pixels that happen only at one slice
    Mat tmp_map;
    cellSegFromSynQuant.fgMap.copyTo(tmp_map);
    size_t page_sz = cellSegFromSynQuant.fgMap.size[1] * cellSegFromSynQuant.fgMap.size[0];
    int width = cellSegFromSynQuant.fgMap.size[1];
    FOREACH_ijk_ptrMAT(cellSegFromSynQuant.idMap){
        size_t idx = vol_sub2ind(i,j,k, width, page_sz);
        if(tmp_map.at<unsigned char>(idx) > 0){
            if((k-1>=0 && tmp_map.at<unsigned char>(idx - page_sz) == 0) &&
                (k+1<tmp_map.size[2] && tmp_map.at<unsigned char>(idx + page_sz) == 0)){
                cellSegFromSynQuant.fgMap.at<unsigned char>(idx) = 0;
            }
        }
    }
    int radius[] = {1,1,0};
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
    volumeDilate(&cellSegFromSynQuant.fgMap, tmp_map, radius, MORPH_ELLIPSE);
    n = connectedComponents3d(&tmp_map, labelMap, 26);
    removeSmallCC(labelMap, n, p4segVol.min_cell_sz, false);
    //cellSegFromSynQuant.fgMap = labelMap > 0;
    bitwise_and(labelMap > 0, seed.otherIdMap==0, cellSegFromSynQuant.fgMap);
    //ccShowSlice3Dmat(cellSegFromSynQuant.fgMap, CV_8U);
}

/**
 * @brief fgGapRemoval: remove the gap defined by principal curvature. We jointly consider the 2d and 3d
 * principal curvature, since 3d may be too liberal to remove too much areas
 * @param seed
 * @param p4segVol
 */
void MicTrackerMain::fgGapRemoval(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol){
    bitwise_and(cellSegFromSynQuant.fgMap, seed.gap3dMap == 0, cellSegFromSynQuant.fgMapGapRemoved);

//    if (!isempty(&cellSegFromSynQuant.fgMapGapRemoved, CV_8U)){
//        Mat seedMapFrom2dMap, mapUnion, mapUnion_label;
//        bitwise_and(cellSegFromSynQuant.fgMap, seed.gap2dMap == 0, seedMapFrom2dMap);// or use -
//        bitwise_or(seedMapFrom2dMap, cellSegFromSynQuant.fgMapGapRemoved, mapUnion); // or use +
//        int numCC = connectedComponents3d(&mapUnion, mapUnion_label, p4segVol.connect4fgGapRemoval);
//        Mat newSeedMap; // CV_8U
//        bool found = findUnrelatedCC(&mapUnion_label, numCC, &cellSegFromSynQuant.fgMapGapRemoved, newSeedMap);
//        if(found){
//            cellSegFromSynQuant.fgMapGapRemoved = cellSegFromSynQuant.fgMapGapRemoved + newSeedMap; //CV_8U
//        }
//    }
}
/**
 * @brief gapBasedRegionSegment: use principal curvature to test if current fg contains > 1 cells.
 * @param seed
 * @param p4segVol
 * @param p4odStats
 */
void MicTrackerMain::gapBasedRegionSegment(synQuantSimple &cellSegFromSynQuant, singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats){
    // 1. remove gaps defined by principal curvature
    fgGapRemoval(cellSegFromSynQuant, seed, p4segVol);
    // 2. check if gaps divided the region into multiple valid seeds
    Mat label_map;
    int n = connectedComponents3d(&cellSegFromSynQuant.fgMapGapRemoved, label_map, p4segVol.neiMap);
    removeSmallCC(label_map, n, p4segVol.min_seed_size, REARRANGE_IDS);
    if(n <= 1){
        cellSegFromSynQuant.cell_num = 1;
        return;
    }
    // 3. gap test: test if the gap is true, if so, split the fg with respect to the multiple seeds
    gapTest2SplitCellTerritory(cellSegFromSynQuant, &label_map, n, seed, p4segVol, p4odStats);
}
/**
 * @brief gapTest2SplitCellTerritory: Given a seeds_map, test if the gap between any two seeds are true or false.
 * If true, keep them. Otherwise, merge the two seeds divided by the gap. If the gap is indeed dimmer than the
 * two seeds it split, it is true. Otherwise false.
 * @param seeds_Map:cv_32s
 * @param n
 * @param seed
 * @param p4segVol
 * @param p4odStats
 */
void MicTrackerMain::gapTest2SplitCellTerritory(synQuantSimple &cellSegFromSynQuant, Mat* seeds_Map /*CV_32S*/, int n,
                                                singleCellSeed &seed, segParameter &p4segVol,
                                                odStatsParameter &p4odStats){
    double maxId;
    minMaxIdx(*seeds_Map, nullptr, &maxId);
    // 1. first grow all seeds until they touch with each other
    Mat grown_seedMap;
    //Mat scoreMap = seed.score2d + seed.score3d;
    bool link_bg2sink = false; // this can force seeds to grow as much as they can
    regionGrow(seeds_Map, n, grown_seedMap, &seed.scoreMap, &cellSegFromSynQuant.fgMap,
               p4segVol.growConnectInRefine, p4segVol.graph_cost_design, link_bg2sink);
    //ccShowSliceLabelMat(seeds_Map);
    //ccShowSliceLabelMat(&grown_seedMap);


    // 2. for the connected seeds, test if the gap among them are true
    // the gap are choosen by dilating the seed regions
    vector<int> gap_tested_true(n*n, 0); // 0 not tested, -1 tested but merge, 1 tested split
    float p_treshold = p4odStats.gapTestThreshold / (p4segVol.gapTestMinMaxRadius[1] - p4segVol.gapTestMinMaxRadius[0]);
    vector<size_t> real_gap_idx;
    for(int r = p4segVol.gapTestMinMaxRadius[0]; r <= p4segVol.gapTestMinMaxRadius[1]; r++){
        vector<vector<size_t>> gap_idx_list;
        extractGapVoxel(&grown_seedMap, &cellSegFromSynQuant.fgMap, n, r, gap_idx_list, gap_tested_true);
        FOREACH_i(gap_idx_list){
            gap_tested_true[i] = 1;
        }
        gap_idx_list.clear();
    }
    if (real_gap_idx.size() > 0){
        vector<vector<int>> groups;
        FOREACH_i(gap_tested_true){
            if(gap_tested_true[i] == -1){
                int target0 = i / n + 1;
                int target1 = i % n + 1;
                groups.push_back({target0, target1});
            }
        }
        mergeIntersectGroups(groups);
        vector<int> label_map (n);
        FOREACH_i(label_map){
            label_map[i] = i + 1;
        }
        FOREACH_i(groups){
            //vec_unique(groups[i]); //sorted also
            for(size_t j = 1; j < groups[i].size(); j++){
                assert(label_map[groups[i][j] - 1] == groups[i][j] || label_map[groups[i][j] - 1] == groups[i][0]);// "Duplicated Regions in two groups");
                label_map[groups[i][j] - 1] = groups[i][0];
            }
        }

        FOREACH_i_MAT(grown_seedMap){
            if(grown_seedMap.at<int>(i) > 0){
                grown_seedMap.at<int>(i) = label_map[grown_seedMap.at<int>(i) - 1];
            }
        }
        //ccShowSliceLabelMat(&grown_seedMap);
        vector<size_t> non_used;
        cellSegFromSynQuant.cell_num = rearrangeIdMap(&grown_seedMap, *cellSegFromSynQuant.idMap, non_used);

        FOREACH_i(real_gap_idx){ // remove the real gaps
            cellSegFromSynQuant.idMap->at<int>(real_gap_idx[i]) = 0;
        }
    }else{
        if(cellSegFromSynQuant.idMap->empty()){
            cellSegFromSynQuant.idMap->create(cellSegFromSynQuant.fgMap.dims,
                                              cellSegFromSynQuant.fgMap.size,
                                              CV_32S);
        }
        FOREACH_i_MAT(cellSegFromSynQuant.fgMap){
            if (cellSegFromSynQuant.fgMap.at<unsigned char>(i) > 0){
                cellSegFromSynQuant.idMap->at<int>(i) = 1;
            }
        }

        cellSegFromSynQuant.cell_num = 1;
    }
}


Mat * MicTrackerMain::extract3d(Mat *data_grayim4d,int curr_frame){
    // test
    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];
    curr_time_point = curr_frame;
    unsigned char *ind = (unsigned char*)data_grayim4d->data + sz_single_frame*curr_time_point; // sub-matrix pointer
    Mat *data_grayim3d = new Mat(3, data_grayim4d->size, CV_8U, ind);
    return data_grayim3d;
}

Mat * MicTrackerMain::extract2d(Mat *data_grayim3d,int curr_slice,int datatype){



    int sz_single_slice =data_rows_cols_slices[0]*data_rows_cols_slices[1];
    Mat *single_slice;
    if(datatype == CV_8U){
        uchar * ind;
        ind = ((uchar*)data_grayim3d->data) + sz_single_slice*curr_slice; // sub-matrix pointer
        single_slice=new Mat(2, data_grayim3d->size, datatype, ind);
    }else if(datatype == CV_32F){
        float * ind;
        ind = ((float*)data_grayim3d->data) + sz_single_slice*curr_slice; // sub-matrix pointer
        single_slice=new Mat(2, data_grayim3d->size, datatype, ind);
    }
    return single_slice;
}

void MicTrackerMain::showSlice(Mat *data_grayim3d,int curr_slice,int datatype){



    int sz_single_slice =data_rows_cols_slices[0]*data_rows_cols_slices[1];
    Mat *single_slice;
    if(datatype == CV_8U){
        uchar * ind;
        ind = ((uchar*)data_grayim3d->data) + sz_single_slice*curr_slice; // sub-matrix pointer
        single_slice=new Mat(2, data_grayim3d->size, datatype, ind);
    }else if(datatype == CV_32F){
        float * ind;
        ind = ((float*)data_grayim3d->data) + sz_single_slice*curr_slice; // sub-matrix pointer
        single_slice=new Mat(2, data_grayim3d->size, datatype, ind);
    }


    QTextStream out(stdout);
    out<<"test:"<<Qt::endl;
    for(int i=0;i<data_rows_cols_slices[0];i++){
        for(int j=0;j<data_rows_cols_slices[1];j++){
            if(datatype == CV_8U){
                out<<single_slice->at<uchar>(i,j)<<"\t";
            }else if(datatype == CV_32F){
                out<<single_slice->at<float>(i,j)<<"\t";
            }

        }
        out<<Qt::endl;
    }
}

size_t MicTrackerMain::sub2idx(int x,int y,int z,int y_size, int xy_size){
    size_t idx = z * xy_size +  x * y_size + y;
    return idx;
}
