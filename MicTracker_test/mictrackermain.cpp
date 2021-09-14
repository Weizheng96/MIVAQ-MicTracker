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

    data4d->copyTo(normalized_data4d);
    assert(normalized_data4d.type() == CV_8U);
    delete data4d;
}
void MicTrackerMain::processSingleFrameAndReturn(int curr_timePoint_in_canvas,const QString &fileName){

    //////////////////////////////////////////////////////////////////
    //           1. get one frame
    //////////////////////////////////////////////////////////////////
    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];
    curr_time_point = curr_timePoint_in_canvas;
    unsigned char *ind = (unsigned char*)normalized_data4d.data + sz_single_frame*curr_time_point; // sub-matrix pointer
    Mat *single_frame = new Mat(3, normalized_data4d.size, CV_8U, ind);
    cell_label_maps[curr_time_point] = Mat::zeros(3, normalized_data4d.size, CV_32S); // int label
    threshold_maps[curr_time_point] = Mat::zeros(3, normalized_data4d.size, CV_8U);


    //////////////////////////////////////////////////////////////////
    //           2. segment
    //////////////////////////////////////////////////////////////////
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    cellSegmentSingleFrame(single_frame, curr_time_point);
    time_points_processed[curr_time_point] = true;
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    qInfo("----------------time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);


    //////////////////////////////////////////////////////////////////
    //           3. save
    //////////////////////////////////////////////////////////////////
    if(!saveSegResults(fileName)){
        qDebug("Failed to save the segmentation results to hard drive");
    }

    qInfo("----------------totally: %ld cells are detected", number_cells[curr_time_point]);
    delete single_frame;
}

void MicTrackerMain::cellSegmentSingleFrame(Mat *data_grayim3d, size_t curr_frame)
{
    //convert data to float
    Mat *dataVolFloat = new Mat(data_grayim3d->dims, data_grayim3d->size, CV_32F);
    data_grayim3d->convertTo(*dataVolFloat, CV_32F);


    //////////////////////////////////////////////////////////////////
    //           1. get PC score
    //////////////////////////////////////////////////////////////////
    float sigma3d[3] = {3.0, 3.0, 3.0};
    principalCv3d(dataVolFloat, principalCurv3d[curr_frame], sigma3d, true,p4segVol.min_intensity);

//    for(int z=0;z<100;z++){
//        ccShowSlice3Dmat(principalCurv3d[curr_frame],CV_32F,z);
//    }

    //////////////////////////////////////////////////////////////////
    //           2. use synQuant for PC score
    //////////////////////////////////////////////////////////////////
    double tmp_min, tmp_max;
    Mat imgIn_temp, imgIn_temp2;

    minMaxIdx(principalCurv3d[curr_frame], &tmp_min, &tmp_max);
    max(principalCurv3d[curr_frame]/tmp_max,0,imgIn_temp);
    sqrt(imgIn_temp,imgIn_temp2);
    Mat imgPCfloat=1-imgIn_temp2;
    Mat imgPCfloat255,imgPC8U;
    imgPCfloat255=imgPCfloat*255;

    variances[curr_frame] = calVarianceStablization(&imgPCfloat255, varMaps[curr_frame], varTrends[curr_frame],
                                                   p4odStats.varAtRatio, p4odStats.gap4varTrendEst);


    imgPCfloat255.convertTo(imgPC8U, CV_8U);
    synQuantSimple seeds_from_synQuant(&imgPC8U, variances[curr_frame], p4segVol, p4odStats);

    for(int z=0;z<100;z++){
        ccShowSliceLabelMat(seeds_from_synQuant.idMap,z);
    }

    //////////////////////////////////////////////////////////////////
    //                 3. refine the seed regions                   //
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


//void cellSegmentMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat, Mat * volStblizedFloat, Mat *idMap /*int*/, int seed_num, Mat *eigMap2d,
//                                           Mat *eigMap3d, Mat *varMap, Mat *stblizedVarMap, vector<int> test_ids){
void MicTrackerMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat,
                                          Mat *idMapIn /*int*/, int seed_num,int curr_frame){
    Mat idMap;
    idMapIn->copyTo(idMap);
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

    }
}

bool MicTrackerMain::saveSegResults(const QString &fileName){
    string fileNameNoExt = fileName.left(fileName.lastIndexOf('.')).toStdString();
    // save label data
    string label_file_name = fileNameNoExt + "_label_map_int32.bin";
    ofstream label_file(label_file_name, ios::binary);
    if (!label_file.is_open()) return false;
    // tmp is local variable, which will be released soon, so we need copyTo
    label_file.write((const char*)(cell_label_maps[curr_time_point].data),
                     cell_label_maps[curr_time_point].elemSize() * cell_label_maps[curr_time_point].total());
    label_file.close();

    return true;
}
