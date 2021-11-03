#include "mictrackermain.h"


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
    principalCurv3d.resize(time_points);
    segScore3d.resize(time_points);
    varMaps.resize(time_points);
    stblizedVarMaps.resize(time_points);
    varTrends.resize(time_points);
    stblizedVarTrends.resize(time_points);
    variances.resize(time_points);

    number_cells.resize(time_points);
    VoxelIdxList.resize(time_points);
    avgItstyVec.resize(time_points);
    areaVec.resize(time_points);
    ctrPt.resize(time_points);
    CoordinateRgs.resize(time_points);


    Mat *data4d;
    int data_sz[4] = {data_rows_cols_slices[0], data_rows_cols_slices[1],
                      data_rows_cols_slices[2], (int)time_points};
    if (data_type == V3D_UINT16) {
        data4d = new Mat(4, data_sz, CV_16U, data_grayim4d);
        *data4d /= 256;
        data4d->convertTo(*data4d,CV_8U);
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

MicTrackerMain::MicTrackerMain()
{
}

void MicTrackerMain::micInitialize(void *data_grayim4d, int _data_type, long bufSize[5]/*(x,y,z,c,t)*/)
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
    principalCurv3d.resize(time_points);
    segScore3d.resize(time_points);
    varMaps.resize(time_points);
    stblizedVarMaps.resize(time_points);
    varTrends.resize(time_points);
    stblizedVarTrends.resize(time_points);
    variances.resize(time_points);

    number_cells.resize(time_points);
    VoxelIdxList.resize(time_points);
    avgItstyVec.resize(time_points);
    areaVec.resize(time_points);
    ctrPt.resize(time_points);
    CoordinateRgs.resize(time_points);


    Mat *data4d;
    int data_sz[4] = {data_rows_cols_slices[0], data_rows_cols_slices[1],
                      data_rows_cols_slices[2], (int)time_points};
    if (data_type == V3D_UINT16) {
        data4d = new Mat(4, data_sz, CV_16U, data_grayim4d);
        *data4d /= 256;
        data4d->convertTo(*data4d,CV_8U);
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
    //           1. get one frame and intialize
    //////////////////////////////////////////////////////////////////
    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];
    curr_time_point = curr_timePoint_in_canvas;
    unsigned char *ind = (unsigned char*)normalized_data4d.data + sz_single_frame*curr_time_point; // sub-matrix pointer
    Mat *single_frame = new Mat(3, normalized_data4d.size, CV_8U, ind);
    cell_label_maps[curr_time_point] = Mat::zeros(3, normalized_data4d.size, CV_32S); // int label

    if(loadSegResults(fileName)){

        //////////////////////////////////////////////////////////////////
        //           2. try load
        //////////////////////////////////////////////////////////////////
        time_points_processed[curr_time_point] = true;
    }else{
        //////////////////////////////////////////////////////////////////
        //           3. segment
        //////////////////////////////////////////////////////////////////
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        cellSegmentSingleFrame(single_frame, curr_time_point);
        time_points_processed[curr_time_point] = true;
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        qInfo("----------------time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);


        //////////////////////////////////////////////////////////////////
        //           4. save
        //////////////////////////////////////////////////////////////////
        if(!saveSegResults(fileName)){
            qDebug("Failed to save the segmentation results to hard drive");
        }else{
            time_points_processed[curr_time_point] = true;
        }
    }

    //////////////////////////////////////////////////////////////////
    //           5. extract information for tracking
    //////////////////////////////////////////////////////////////////

    seedsInfo(single_frame,&cell_label_maps[curr_time_point],number_cells[curr_time_point]);

    qInfo("----------------totally: %ld cells are detected", number_cells[curr_time_point]);
    delete single_frame;
}

void MicTrackerMain::cellSegmentSingleFrame(Mat *data_grayim3d, size_t curr_frame)
{
    //////////////////////////////////////////////////////////////////
    //           1. get foreground
    //////////////////////////////////////////////////////////////////
    Mat *dataVolFloat = new Mat(data_grayim3d->dims, data_grayim3d->size, CV_32F);
    data_grayim3d->convertTo(*dataVolFloat, CV_32F);

    //////////////////////////////////////////////////////////////////
    //           2. get PC score
    //////////////////////////////////////////////////////////////////
    float sigma3d[3] = {3.0, 3.0, 3.0};
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    principalCv3d(dataVolFloat, principalCurv3d[curr_frame], sigma3d, true,p4segVol.min_intensity);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    qInfo("---------------- Principal curvature time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);


    //////////////////////////////////////////////////////////////////
    //           3. use synQuant for PC score
    //////////////////////////////////////////////////////////////////
    double tmp_min, tmp_max;
    Mat imgIn_temp, imgIn_temp2;

    minMaxIdx(principalCurv3d[curr_frame], &tmp_min, &tmp_max);
    max(principalCurv3d[curr_frame]/tmp_max,0,imgIn_temp);
    sqrt(imgIn_temp,imgIn_temp2);
    Mat imgPCfloat=1-imgIn_temp2;
    Mat imgPCfloat255,imgPC8U;
    imgPCfloat255=imgPCfloat*255;
    imgPCfloat.copyTo(segScore3d[curr_frame]);


    variances[curr_frame] = calVarianceStablization(&imgPCfloat255, varMaps[curr_frame], varTrends[curr_frame],
                                                   p4odStats.varAtRatio, p4odStats.gap4varTrendEst);
    imgPCfloat255.convertTo(imgPC8U, CV_8U);
    synQuantSimple seeds_from_synQuant(&imgPC8U, variances[curr_frame], p4segVol, p4odStats);

    fgThreshold=getMaxContrast(data_grayim3d);
//    Mat foreground;
//    getForeground(data_grayim3d, foreground);
//    Mat FgMask;
//    foreground.convertTo(FgMask,CV_32S);
//    FgMask/=255;
//    *seeds_from_synQuant.idMap=seeds_from_synQuant.idMap->mul(FgMask);
//    removeSmallCC(*seeds_from_synQuant.idMap, seeds_from_synQuant.cell_num, p4segVol.min_cell_sz, true);

    //////////////////////////////////////////////////////////////////
    //                 4. refine the seed regions                   //
    //////////////////////////////////////////////////////////////////
    regionWiseAnalysis4d(data_grayim3d, &imgPCfloat255,seeds_from_synQuant.idMap,
                         seeds_from_synQuant.cell_num,curr_frame);
}

void MicTrackerMain::getForeground(Mat * data_grayim3d, Mat &foreground){
    size_t thre=getMaxContrast(data_grayim3d);
    Mat foreGround_raw=*data_grayim3d>thre;
    foreGround_raw.convertTo(foreGround_raw,CV_32F);
    // 3d median filter
    averageSmooth3Ddata(foreGround_raw, 3);
    Mat foreGround_median=(foreGround_raw>127);
    // fill holes
    Mat foreGround_median_reverse=255-foreGround_median;
    Mat HoleCC;
    int holeNum=connectedComponents3d(&foreGround_median_reverse,HoleCC,26);
    vector<vector<size_t>> voxList_holes;
    extractVoxIdxList(&HoleCC, voxList_holes, holeNum);
    size_t sz[3]={0};
    sz[0]=foreGround_raw.size[1];
    sz[1]=foreGround_raw.size[0];
    sz[2]=foreGround_raw.size[2];
    for(size_t i=0;i<voxList_holes.size();i++){
        if(notTouchBoundary(voxList_holes[i],sz)){
            for(size_t j=0;j<voxList_holes[i].size();j++){
                foreGround_median.at<unsigned char>(voxList_holes[i][j])=255;
            }
        }
    }

    // imopen
    int elementSize[3]={7,7,7};
    Mat MatTemp,foreground_open;
    volumeErode(&foreGround_median, MatTemp, elementSize, MORPH_RECT);
    volumeDilate(&MatTemp, foreground_open, elementSize, MORPH_RECT);

    // imclose
    volumeDilate(&foreground_open, MatTemp, elementSize, MORPH_RECT);
    volumeErode(&MatTemp, foreground, elementSize, MORPH_RECT);
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
void MicTrackerMain::seedsInfo(Mat *single_frame, Mat *label_map, size_t cell_num){
    //////////////////////////////////////////////////////////////////
    //           1. extract information
    //////////////////////////////////////////////////////////////////
    Mat VolFloat;
    single_frame->convertTo(VolFloat, CV_32F);
    size_t sz[3]={0};
    sz[0]=VolFloat.size[1];
    sz[1]=VolFloat.size[0];
    sz[2]=VolFloat.size[2];

//    vector<float> ctr_temp(3);
//    vector<size_t> bg_temp(6);
    ctrPt[curr_time_point].resize(cell_num);
    CoordinateRgs[curr_time_point].resize(cell_num);

    // cell index list
    vector<vector<size_t>> voxList_exist_cells;
    extractVoxIdxList(label_map, voxList_exist_cells, cell_num);
    // cell intensity
    vector<float> cell_mean_intensity;
    regionAvgIntensity(&VolFloat, voxList_exist_cells, cell_mean_intensity);
    // for each cell
    vector<size_t> exist_cell_size (voxList_exist_cells.size());
    FOREACH_i(voxList_exist_cells) {
        // cell size
        exist_cell_size[i] = voxList_exist_cells[i].size();
        // cell center & cell boundary
        getCenterAndBoundary(voxList_exist_cells[i],sz,ctrPt[curr_time_point][i],CoordinateRgs[curr_time_point][i]);
//        ctrPt[curr_time_point][i]=ctr_temp;
//        CoordinateRgs[curr_time_point][i]=bg_temp;
    }

    //////////////////////////////////////////////////////////////////
    //           2. send information to class
    //////////////////////////////////////////////////////////////////
    VoxelIdxList[curr_time_point]=voxList_exist_cells; //cell index
    avgItstyVec[curr_time_point]=cell_mean_intensity; // cell intensity
    areaVec[curr_time_point]=exist_cell_size; // cell size
}

void MicTrackerMain::getCenterAndBoundary(vector<size_t> voxList, size_t sz[3], vector<float> &ctr, vector<size_t> &boundary){
    size_t y_bd[2]={sz[0],0};
    size_t x_bd[2]={sz[1],0};
    size_t z_bd[2]={sz[2],0};

    size_t y_sum=0;
    size_t x_sum=0;
    size_t z_sum=0;

    size_t y_temp=0;
    size_t x_temp=0;
    size_t z_temp=0;

    size_t sz_slice=sz[0]*sz[1];
    // size_t idx = z*xy_size + c + increment;
    for(size_t i=0;i<voxList.size();i++){
        z_temp=voxList[i]/sz_slice;
        x_temp=(voxList[i]%sz_slice)/sz[0];
        y_temp=voxList[i]%sz[0];

        z_sum+=z_temp;
        x_sum+=x_temp;
        y_sum+=y_temp;

        if(z_temp>z_bd[1]){z_bd[1]=z_temp;}
        if(z_temp<z_bd[0]){z_bd[0]=z_temp;}
        if(x_temp>x_bd[1]){x_bd[1]=x_temp;}
        if(x_temp<x_bd[0]){x_bd[0]=x_temp;}
        if(y_temp>y_bd[1]){y_bd[1]=y_temp;}
        if(y_temp<y_bd[0]){y_bd[0]=y_temp;}
    }


    ctr.resize(3);
    ctr[0]=((float) y_sum)/voxList.size();
    ctr[1]=((float) x_sum)/voxList.size();
    ctr[2]=((float) z_sum)/voxList.size();

    boundary.resize(6);
    boundary[0]=y_bd[0];
    boundary[1]=y_bd[1];
    boundary[2]=x_bd[0];
    boundary[3]=x_bd[1];
    boundary[4]=z_bd[0];
    boundary[5]=z_bd[1];
}

bool MicTrackerMain::notTouchBoundary(vector<size_t> voxList, size_t sz[3]){


    size_t y=0;
    size_t x=0;
    size_t z=0;

    size_t sz_slice=sz[0]*sz[1];
    // size_t idx = z*xy_size + c + increment;
    for(size_t i=0;i<voxList.size();i++){
        z=voxList[i]/sz_slice;
        if(z==0||z==sz[2]-1){return false;}
        x=(voxList[i]%sz_slice)/sz[0];
        if(x==0||x==sz[1]-1){return false;}
        y=voxList[i]%sz[0];
        if(y==0||y==sz[0]-1){return false;}
    }
    return true;
}


//void cellSegmentMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataVolFloat, Mat * volStblizedFloat, Mat *idMap /*int*/, int seed_num, Mat *eigMap2d,
//                                           Mat *eigMap3d, Mat *varMap, Mat *stblizedVarMap, vector<int> test_ids){
void MicTrackerMain::regionWiseAnalysis4d(Mat *data_grayim3d, Mat *dataPC255Float,
                                          Mat *idMapIn /*int*/, int seed_num,int curr_frame){
    Mat idMap;
    idMapIn->copyTo(idMap);
    vector<vector<size_t>> voxIdxList(seed_num);
    extractVoxIdxList(&idMap, voxIdxList, seed_num);

    int cell_cnt = 0;
    //int debug_cell_id = 2;
    bool test=true;

    FOREACH_i(voxIdxList){
        int seed_id = i + 1;
        qInfo("------------------#%ld, start %d seed process, size: %ld---------------------", i, seed_id,
              voxIdxList[seed_id-1].size());

        singleCellSeed seed;
        cropSeed(seed_id, voxIdxList[seed_id-1], data_grayim3d,dataPC255Float, &idMap,
                curr_frame, seed, p4segVol);
        refineSeed2Region(seed, p4odStats, p4segVol);

        qInfo("----------------%d cell output from this seed-------------------", seed.outCell_num);
        if(seed.outCell_num <= 0) continue;
        subVolReplace(cell_label_maps[curr_time_point], CV_32S, seed.outputIdMap, seed.crop_range_yxz, cell_cnt);
        cell_cnt += seed.outCell_num;
    }
    number_cells[curr_time_point] = cell_cnt;
}

void MicTrackerMain::cropSeed(int seed_id, vector<size_t> idx_yxz, Mat *data_grayim3d,Mat *dataPC255Float,
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
    seed.outputIdMap = Mat(seed.idMap.dims, seed.idMap.size, CV_32S, Scalar(-1));
    subVolExtract(dataPC255Float, CV_32F, seed.scoreMap, seed.crop_range_yxz);
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


    if(false){// for test
        size_t intensityAll=0;
        size_t sz=0;
        for(size_t i=0;i<seed.volUint8.total();i++){
            if(seed.seedMap.at<unsigned char>(i)>0){
                intensityAll+=seed.volUint8.at<unsigned char>(i);
                sz++;
            }
        }
        if(fgThreshold*sz>intensityAll){
            seed.outCell_num=0;
        }else{
            seed.outCell_num=1;
            Mat seedNew=seed.seedMap/255;
            seedNew.convertTo(seed.outputIdMap,CV_32S);
        }
    }else{
        Mat PCseed_all= (seed.eigMap3d <= 0.0);
        PCseed_all.convertTo(PCseed_all,CV_8U);
        Mat PCseed_valid;
        bitwise_and(PCseed_all,seed.seedMap,PCseed_valid);
        // get each PC seed
        Mat PCseed_valid_CC;
        int PCseedNum=connectedComponents3d(&PCseed_valid,PCseed_valid_CC,26);

        Mat seed_valid_CC;
        seed.seedMap.convertTo(seed_valid_CC,CV_32S);
        seed_valid_CC/=255;

        removeSmallCC(PCseed_valid_CC, PCseedNum,p4segVol.min_cell_sz, true);//p4segVol.min_cell_sz,p4segVol.min_seed_size

        vector<vector<int>> seedLst;
        extractVoxIdxList(&PCseed_valid_CC, seedLst, PCseedNum, false);

//        ccShowSliceLabelMat(seed_valid_CC,8);
//        ccShowSliceLabelMat(PCseed_valid_CC,8);

        // minMaxIdx(PCseed_valid_CC, &tmp_min, &tmp_max);


        if(PCseedNum<2){
            seed.outCell_num=1;
            seed_valid_CC.convertTo(seed.outputIdMap,CV_32S);
        }else{

            int width=seed.crop_range_yxz[1].end-seed.crop_range_yxz[1].start;
            int height=seed.crop_range_yxz[0].end-seed.crop_range_yxz[0].start;
            int thick=seed.crop_range_yxz[2].end-seed.crop_range_yxz[2].start;


            seed.scoreMap=255-seed.scoreMap;
            seed.scoreMap.convertTo(seed.scoreMap,CV_8U);
//            ccShowSlice3Dmat(seed.scoreMap, CV_8U,8);

            seed.outputIdMap=seed.outputIdMap.mul(seed_valid_CC);

            WaterShed_WZ::Watershed3D_WZ(
                seed.scoreMap,
                width,
                height,
                thick,
                &seed.outputIdMap,
                seedLst);

//            ccShowSliceLabelMat(seed.outputIdMap,8);

//            seed.outputIdMap=seed.outputIdMap.mul(seed_valid_CC);

//            ccShowSliceLabelMat(seed.outputIdMap,8);

            seed.outCell_num=PCseedNum;
        }
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

    // save 3d principal map
    string p3d_file_name = fileNameNoExt + "_principal3d_map_single.bin";
    ofstream p3d_file = ofstream(p3d_file_name, ios::binary);
    if (!p3d_file.is_open())  return false;
    p3d_file.write((const char*)(principalCurv3d[curr_time_point].data),
                     principalCurv3d[curr_time_point].elemSize() * principalCurv3d[curr_time_point].total());
    p3d_file.close();
    //ccShowSlice3Dmat(principalCurv3d[curr_time_point], CV_32F);

    // save 3d score map
    string sc3d_file_name = fileNameNoExt + "_segScore_map.bin";
    ofstream sc3d_file = ofstream(sc3d_file_name, ios::binary);
    if (!sc3d_file.is_open())  return false;
    sc3d_file.write((const char*)(segScore3d[curr_time_point].data),
                     segScore3d[curr_time_point].elemSize() * segScore3d[curr_time_point].total());
    sc3d_file.close();
    //ccShowSlice3Dmat(principalCurv3d[curr_time_point], CV_32F);

    return true;
}

bool MicTrackerMain::loadSegResults(const QString &fileName){
    QString fileNameNoExt = fileName.left(fileName.lastIndexOf('.'));
    // read label data
    QString label_file_name = fileNameNoExt + "_label_map_int32.bin";
    QFile label_file(label_file_name);
    if (!label_file.open(QIODevice::ReadOnly)){
        return false;
    }
    Mat(3, normalized_data4d.size, CV_32S, label_file.readAll().data()).copyTo(cell_label_maps[curr_time_point]);
    label_file.close();
    double tmp_maxVal;
    minMaxIdx(cell_label_maps[curr_time_point], nullptr, &tmp_maxVal);//minMaxLoc for 2d
    number_cells[curr_time_point] = round(tmp_maxVal);
    //ccShowSliceLabelMat(cell_label_maps[curr_time_point]);

    // read 3d principal map
    QString p3d_file_name = fileNameNoExt + "_principal3d_map_single.bin";
    QFile p3d_file(p3d_file_name);
    if (!p3d_file.open(QIODevice::ReadOnly))  return false;
    Mat(3, normalized_data4d.size, CV_32F, p3d_file.readAll().data()).copyTo(principalCurv3d[curr_time_point]);
    p3d_file.close();
    //ccShowSlice3Dmat(principalCurv3d[curr_time_point], CV_32F);

    // read 3d score map
    QString sc3d_file_name = fileNameNoExt + "_segScore_map.bin";
    QFile sc3d_file(sc3d_file_name);
    if (!sc3d_file.open(QIODevice::ReadOnly))  return false;
    Mat(3, normalized_data4d.size, CV_32F, sc3d_file.readAll().data()).copyTo(segScore3d[curr_time_point]);
    sc3d_file.close();

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////
/// refine part
/// ////////////////////////////////////////////////////////////////////////////////////

void MicTrackerMain::localIdMap2GlobalIdMap(vector<Mat> &localMap, vector<Mat> &globalMap, vector<size_t> cellNum)
{
    size_t N=localMap.size();
    size_t cell_cnt=0;
    Range range_allyxz[3];
    for(size_t i=0;i<3;i++){range_allyxz[i]=Range(0, localMap[0].size[i]);}
    globalMap.resize(N);
    for(size_t frameCnt=0;frameCnt<N;frameCnt++){
        globalMap[frameCnt]=Mat::zeros(localMap[0].dims,localMap[0].size,CV_32S);
        subVolReplace(globalMap[frameCnt], CV_32S, localMap[frameCnt], range_allyxz, cell_cnt);
        cell_cnt += cellNum[frameCnt];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief loop2_Initialization
/// 1. initialzie "cell_label_maps_old" based on "cell_label_maps" and "number_cells"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MicTrackerMain::loop2_Initialization(){

    // initialize
//    localIdMap2GlobalIdMap(cell_label_maps, cell_label_maps_globalId, number_cells);

    // initialize local to global id
    localId2GlobalId.clear();
    size_t cnt=0;
    localId2GlobalId.resize(cell_label_maps.size());
    for(size_t time=0;time<cell_label_maps.size();time++){
        localId2GlobalId[time].resize(number_cells[time]);
        for(size_t localId=0;localId<number_cells[time];localId++){
            localId2GlobalId[time][localId]=cnt;
            cnt++;
        }
    }

    // initialize global to local id
    globalId2LocalId.clear();
    globalId2LocalId.resize(cnt);
    cnt=0;
    for(size_t time=0;time<cell_label_maps.size();time++){
        for(size_t localId=0;localId<number_cells[time];localId++){
            globalId2LocalId[cnt].push_back(time);
            globalId2LocalId[cnt].push_back(localId);
            cnt++;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief MicTrackerMain::refineBytemporalInfo_loop2
/// find reference and split valid detection based on the reference, update result to "cell_label_maps_old"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MicTrackerMain::refineBytemporalInfo_loop2(vector<size_t> numCellResVec_cur,vector<vector<vector<size_t>>> dtctIncldLst_cur)
{
    ////////////////////////////////////////////////////////////////////////////////////
    /// 1. initialize
    ////////////////////////////////////////////////////////////////////////////////////
    loop2_Initialization();

    //////////////////////////////////////////////////////
    /// find valid detections: detections have more than one cells
    vector<size_t> dtcts2SplitId;
    for(size_t i=0;i<numCellResVec_cur.size();i++){
       if(numCellResVec_cur[i]>1){
           dtcts2SplitId.push_back(i);
       }
    }
    const size_t N=dtcts2SplitId.size();

    //////////////////////////////////////////////////////
    /// find valid detections with valid reference
    vector<vector<bool>> dtctsValidRefer;
    dtctsValidRefer.resize(2);
    size_t temp_overallId;
    bool temp_flag; // if need break
    size_t ii;
    for(size_t dirCnt=0;dirCnt<2;dirCnt++){
        dtctsValidRefer[dirCnt].resize(N);
        fill (dtctsValidRefer[dirCnt].begin(),dtctsValidRefer[dirCnt].end(),false);
        for(size_t i=0; i<N ;i++){
            // ii is glocal id;
            ii=dtcts2SplitId[i];
            // have at least two reference detections
            if(dtctIncldLst_cur[dirCnt][ii].size()<2){continue;}
            // each reference detection only have one cell
            temp_flag=false;
            for(size_t refCnt=0;refCnt<dtctIncldLst_cur[dirCnt][ii].size();refCnt++){
                temp_overallId=dtctIncldLst_cur[dirCnt][ii][refCnt];
                if(numCellResVec_cur[temp_overallId]!=1){
                    temp_flag=true;
                    break;
                }
            }
            if(temp_flag){continue;}
            dtctsValidRefer[dirCnt][i]=true;
        }
    }

    ////////////////////////////////////////////////////////
    /// If both directions are valid, find the better one
    size_t temp_diffPre;
    size_t temp_diffPost;
    size_t temp_refDtctGlobalId;
    vector<size_t> temp_refDtctLocalId;
    vector<size_t> temp_maxArea;
    vector<size_t> temp_minArea;
    size_t temp_Area;

    for(size_t i=0; i<N ;i++)
    {
        if(dtctsValidRefer[0][i]==true&&dtctsValidRefer[1][i]==true)
        {
            ii=dtcts2SplitId[i];
            // 1. If cell number is different, find the direction whose cell number is closest to estimate cell number
            temp_diffPre=abs(numCellResVec_cur[ii]-dtctIncldLst_cur[0][ii].size());
            temp_diffPost=abs(numCellResVec_cur[ii]-dtctIncldLst_cur[1][ii].size());
            if(temp_diffPre!=temp_diffPost)
            {
                if(temp_diffPre<temp_diffPost)
                {
                    dtctsValidRefer[1][i]=false; // choose pre frame
                }
                else
                {
                    dtctsValidRefer[0][i]=false; // chose post frame
                }
            }
            else
            {
                // 2. find the direction whose detection size is more stable
                temp_maxArea={0,0};
                temp_minArea={p4segVol.max_cell_sz,p4segVol.max_cell_sz};
                for(size_t dirCnt=0;dirCnt<2;dirCnt++)
                {
                    for(size_t refCnt=0;refCnt<dtctIncldLst_cur[dirCnt][ii].size();refCnt++)
                    {
                        temp_refDtctGlobalId=refCnt<dtctIncldLst_cur[dirCnt][ii][refCnt];
                        temp_refDtctLocalId=globalId2LocalId[temp_refDtctGlobalId];
                        temp_Area=areaVec[temp_refDtctLocalId[0]][temp_refDtctLocalId[1]];
                        if(temp_Area>temp_maxArea[dirCnt]){temp_maxArea[dirCnt]=temp_Area;}
                        if(temp_Area<temp_minArea[dirCnt]){temp_minArea[dirCnt]=temp_Area;}
                    }
                }
                if((temp_minArea[0]/temp_maxArea[0])>(temp_minArea[1]/temp_maxArea[1]))
                {
                    dtctsValidRefer[1][i]=false; // choose pre frame
                }
                else
                {
                    dtctsValidRefer[0][i]=false; // choose post frame
                }
            }

        }
    }


    ////////////////////////////////////////////////////////////////////////////////////
    /// 2. split detections with valid reference and update to cell_label_maps_cur
    ////////////////////////////////////////////////////////////////////////////////////
    vector<size_t> cell_cnt(cell_label_maps.size(),0);
    size_t NN=numCellResVec_cur.size();
    vector<bool> splitedDetection(NN,false);
    vector<size_t> currentLocalId;
    vector<size_t> currentIdx;
    vector<size_t> refLocalId;
    vector<vector<size_t>> refIdx;
    size_t globalId;

    for(size_t i=0; i<N ;i++){
        for(size_t dirCnt=0;dirCnt<2;dirCnt++){
            if(dtctsValidRefer[dirCnt][i]){
                ii=dtcts2SplitId[i];
                // find input for refineBytemporalInfo_single
                currentLocalId=globalId2LocalId[ii];
                currentIdx=VoxelIdxList[currentLocalId[0]][currentLocalId[1]];
                refIdx.clear();
                for(size_t refCnt=0;refCnt<dtctIncldLst_cur[dirCnt][ii].size();refCnt++)
                {
                    globalId=dtctIncldLst_cur[dirCnt][ii][refCnt];
                    refLocalId=globalId2LocalId[globalId];
                    refIdx.push_back(VoxelIdxList[refLocalId[0]][refLocalId[1]]);
                }
                // split and update cell_label_maps_cur
                singleCellSeed seed;
                refineBytemporalInfo_single(currentLocalId[1]+1,currentIdx,
                        currentLocalId[0],refIdx,cell_cnt,seed);
                // only one direction is valid
                splitedDetection[ii]=true;
                break;
            }
        }
    }

    ////////////////
    /// test
//    ccShowSliceLabelMat(cell_label_maps[0],12);
//    ccShowSliceLabelMat(cell_label_maps_cur[0],12);

    ////////////////////////////////////////////////////////////////////////////////////
    /// 3. update unsplitted detections to cell_label_maps_cur
    ////////////////////////////////////////////////////////////////////////////////////
    vector<size_t> temp_Idx;
    for(size_t cnt=0;cnt<NN;cnt++)
    {
        if(!splitedDetection[cnt]){
            currentLocalId=globalId2LocalId[cnt];
            temp_Idx=VoxelIdxList[currentLocalId[0]][currentLocalId[1]];
            cell_cnt[currentLocalId[0]]++;
            fillMat(cell_label_maps_cur[currentLocalId[0]], temp_Idx, cell_cnt[currentLocalId[0]]);
        }
    }

    ////////////////
    /// test
//    ccShowSliceLabelMat(cell_label_maps[0],12);
//    ccShowSliceLabelMat(cell_label_maps_cur[0],12);

    ////////////////////////////////////////////////////////////////////////////////////
    /// 4. update detections' infomation
    ////////////////////////////////////////////////////////////////////////////////////
    number_cells=cell_cnt;
    for(size_t i=0;i<(size_t)time_points;i++){
        cell_label_maps_cur[i].copyTo(cell_label_maps[i]);
    }

    ////////////////
    /// test
    ccShowSliceLabelMat(cell_label_maps[0],12);

    long sz_single_frame = data_rows_cols_slices[0]*data_rows_cols_slices[1]*data_rows_cols_slices[2];
    for(curr_time_point=0;curr_time_point<time_points;curr_time_point++){
        unsigned char *ind = (unsigned char*)normalized_data4d.data + sz_single_frame*curr_time_point; // sub-matrix pointer
        Mat *single_frame = new Mat(3, normalized_data4d.size, CV_8U, ind);
        seedsInfo(single_frame,&cell_label_maps[curr_time_point],number_cells[curr_time_point]);
    }
}

void MicTrackerMain::fillMat(Mat &dat, vector<size_t> idx, size_t value)
{
    for(size_t i=0;i<idx.size();i++)
    {
        dat.at<int>(idx[i])=value;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief MicTrackerMain::refineBytemporalInfo_single
/// This function is used to split one detection into several based on reference and save the result to "cell_label_maps_cur"
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MicTrackerMain::refineBytemporalInfo_single(int currentID,vector<size_t> currentIdx, size_t currentFrame,
                                                    vector<vector<size_t>> refIdx,vector<size_t> &cell_cnt,
                                                 singleCellSeed &seed)
{
//    singleCellSeed seed;
    Mat *scMat=&segScore3d[currentFrame];
    Mat *idMap;
    // if this is first time use this function, intialize "cell_label_maps_cur"
    if(accumulate(cell_cnt.begin(), cell_cnt.end(), 0)==0){
        cell_label_maps_cur.resize(time_points);
        for(size_t i=0;i<(size_t)time_points;i++){
            cell_label_maps_cur[i] = Mat::zeros(3, normalized_data4d.size, CV_32S); // int label
        }
    }

    // crop the detection and spilit based on the linkage
    idMap=&cell_label_maps[currentFrame];
    cropSeed(currentID,currentIdx,refIdx,scMat,idMap, seed, p4segVol);

    // update cell_label_maps_cur
    subVolReplace(cell_label_maps_cur[currentFrame], CV_32S, seed.outputIdMap, seed.crop_range_yxz, cell_cnt[currentFrame]);
    cell_cnt[currentFrame] += seed.outCell_num;

    ////////////////
    /// test
//    int z=50;
//    if(currentFrame==0){
//        if(seed.crop_range_yxz[2].start<=z){
//            if(seed.crop_range_yxz[2].end>z){
//                ccShowSliceLabelMat(cell_label_maps[0],z);
//                ccShowSliceLabelMat(seed.outputIdMap,z-seed.crop_range_yxz[2].start);
//                ccShowSliceLabelMat(cell_label_maps_cur[0],z);
//            }
//        }
//    }

}

///////////////////////////////////////////////////////////////////////////////////////////
/// \brief MicTrackerMain::cropSeed
/// This function is used to split one detection into several based on reference, used by "refineBytemporalInfo_single"
///////////////////////////////////////////////////////////////////////////////////////////

void MicTrackerMain::cropSeed(int currentID, vector<size_t> currentIdx,vector<vector<size_t>> refIdx,
                              Mat *scMat, Mat *idMap, singleCellSeed &seed, segParameter p4segVol)
{
    vector<vector<int>> seedLst;
    seedLst.resize(refIdx.size());

    size_t yc;
    size_t xc;
    size_t zc;
    size_t vxc;
    Mat mask;


    seed.id = currentID;
    seed.idx_yxz = currentIdx; // deep copy
    seed.y.resize(currentIdx.size());
    seed.x.resize(currentIdx.size());
    seed.z.resize(currentIdx.size());
    vec_ind2sub(currentIdx, seed.y, seed.x, seed.z, idMap->size);// MatSize itself is a vector

    getRange(seed.y, p4segVol.shift_yxz[0], idMap->size[0], seed.crop_range_yxz[0]);
    getRange(seed.x, p4segVol.shift_yxz[1], idMap->size[1], seed.crop_range_yxz[1]);
    getRange(seed.z, p4segVol.shift_yxz[2], idMap->size[2], seed.crop_range_yxz[2]);


    subVolExtract(idMap, CV_32S, seed.idMap, seed.crop_range_yxz);
    seed.seedMap = seed.idMap == seed.id;

    seed.seedMap.convertTo(mask,CV_32S);;
    mask=mask/255;



    seed.outputIdMap = Mat(seed.idMap.dims, seed.idMap.size, CV_32S, Scalar(-1));

    subVolExtract(scMat, CV_32F, seed.scoreMap, seed.crop_range_yxz);
    //ccShowSliceLabelMat(idMap);
    vec_sub2ind(seed.idx_yxz_cropped, vec_Minus(seed.y, seed.crop_range_yxz[0].start),
            vec_Minus(seed.x, seed.crop_range_yxz[1].start),
            vec_Minus(seed.z, seed.crop_range_yxz[2].start), seed.idMap.size);// MatSize itself is a vector
    //ccShowSliceLabelMat(idMap);

    ///////////////////////////////////////////////
    /// for referrence cells in the other frame
    ///////////////////////////////////////////////
    seed.otherIdMap=Mat::zeros(seed.idMap.dims,seed.idMap.size,CV_32S);
    bool isMember;
    for(size_t i=0;i<refIdx.size();i++){
//        vector<size_t> vxLst=refIdx[i];
        vector<int> y;
        vector<int> x;
        vector<int> z;
        y.resize(refIdx[i].size());
        x.resize(refIdx[i].size());
        z.resize(refIdx[i].size());
        vec_ind2sub(refIdx[i], y, x, z, idMap->size);
        yc=avg(y);
        xc=avg(x);
        zc=avg(z);
        vol_sub2ind(vxc, yc, xc, zc, idMap->size);

        isMember=(find(currentIdx.begin(), currentIdx.end(), vxc) !=currentIdx.end());
        if(isMember){
            // if the center of refer is in the detection
            // for method #1
            vol_sub2ind(vxc, yc-seed.crop_range_yxz[0].start, xc-seed.crop_range_yxz[1].start,
                    zc-seed.crop_range_yxz[2].start, seed.idMap.size);
            seed.otherIdMap.at<int>(vxc)=i+1;
            seedLst[i].push_back(vxc);
        }else{
            // if the center of refer is not in the detection
            for(size_t j=0;j<refIdx[i].size();j++){
                // for method #1
                yc=y[j];
                xc=x[j];
                zc=z[j];
                vol_sub2ind(vxc, yc-seed.crop_range_yxz[0].start, xc-seed.crop_range_yxz[1].start,
                        zc-seed.crop_range_yxz[2].start, seed.idMap.size);
                seed.otherIdMap.at<int>(vxc)=i+1;
                seedLst[i].push_back(vxc);
            }
        }
        y.clear();
        x.clear();
        z.clear();

    }
    seed.otherIdMap=seed.otherIdMap.mul(mask);

    seed.outCell_num=refIdx.size();


//    for(int i=0;i<(seed.crop_range_yxz[2].end-seed.crop_range_yxz[2].start);i++)
//    {
//        ccShowSliceLabelMat(seed.otherIdMap,i);
//    }

    // method #1

    int width=seed.crop_range_yxz[1].end-seed.crop_range_yxz[1].start;
    int height=seed.crop_range_yxz[0].end-seed.crop_range_yxz[0].start;
    int thick=seed.crop_range_yxz[2].end-seed.crop_range_yxz[2].start;


    seed.scoreMap=(1-seed.scoreMap)*255;
    seed.scoreMap.convertTo(seed.scoreMap,CV_8U);

    seed.outputIdMap=seed.outputIdMap.mul(mask);
    WaterShed_WZ::Watershed3D_WZ(
        seed.scoreMap,
        width,
        height,
        thick,
        &seed.outputIdMap,
        seedLst);
    seed.outputIdMap=seed.outputIdMap.mul(mask);
//    mask.convertTo(seed.outputIdMap,CV_32S);


//    for(int i=0;i<(seed.crop_range_yxz[2].end-seed.crop_range_yxz[2].start);i++)
//    {
//        ccShowSliceLabelMat(seed.outputIdMap,i);
//    }

//    vxLst.clear();


}

int MicTrackerMain::avg ( vector<int> v )
{
        float sum = accumulate(v.begin(), v.end(), 0.0);
        size_t mean =round( sum / v.size() );
        return mean;
}
