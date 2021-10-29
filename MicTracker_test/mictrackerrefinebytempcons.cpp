#include "mictrackerrefinebytempcons.h"


#ifndef Pi
#define Pi 3.141592653589793238462643
#endif

using namespace cv;
using namespace std;

MicTrackerRefineByTempCons::MicTrackerRefineByTempCons(void *data_grayim4d, int _data_type, long bufSize[5],
const QStringList &filelist){

    micTracker.micInitialize(data_grayim4d,_data_type,bufSize);
    for(size_t i=0;i<(size_t)micTracker.time_points;i++){
        micTracker.processSingleFrameAndReturn(i,filelist.at(i));
    }

    // info from micPerFrame (this is deep copy of mat!!)

    cell_label_maps.resize(micTracker.time_points);
    for(size_t i=0;i<(size_t)micTracker.time_points;i++){
        micTracker.cell_label_maps[i].copyTo(cell_label_maps[i]);
    }

    principalCurv3d.resize(micTracker.time_points);
    for(size_t i=0;i<(size_t)micTracker.time_points;i++){
        micTracker.principalCurv3d[i].copyTo(principalCurv3d[i]);
        principalCurv3d[i]/=255;
    }

    segScore3d.resize(micTracker.time_points);
    for(size_t i=0;i<(size_t)micTracker.time_points;i++){
        micTracker.segScore3d[i].copyTo(segScore3d[i]);
    }

    data_rows_cols_slices=micTracker.data_rows_cols_slices;
    time_points=micTracker.time_points;
    micTracker.normalized_data4d.copyTo(normalized_data4d);

    number_cells=micTracker.number_cells;
    VoxelIdxList=micTracker.VoxelIdxList;
    avgItstyVec=micTracker.avgItstyVec;
    for(size_t i=0;i<avgItstyVec.size();i++){
        for(size_t j=0;j<avgItstyVec[i].size();j++){
            avgItstyVec[i][j]/=255;
        }
    }
    areaVec=micTracker.areaVec;
    ctrPt=micTracker.ctrPt;
    CoordinateRgs=micTracker.CoordinateRgs;


    effZinFmTwoEndIdx=Mat::zeros(micTracker.time_points,2,CV_8U);
    dtctMtScVec.resize(micTracker.time_points);

}

// not used
MicTrackerRefineByTempCons::MicTrackerRefineByTempCons(MicTrackerMain micPerFrame,
                                                       const QStringList &filelist)
{
    for(size_t i=0;i<micPerFrame.time_points;i++){
        micPerFrame.processSingleFrameAndReturn(i,filelist.at(i));
    }

    // info from micPerFrame (this is deep copy of mat!!)

    cell_label_maps.resize(micPerFrame.time_points);
    for(size_t i=0;i<micPerFrame.time_points;i++){
        micPerFrame.cell_label_maps[i].copyTo(cell_label_maps[i]);
    }

    principalCurv3d.resize(micPerFrame.time_points);
    for(size_t i=0;i<micPerFrame.time_points;i++){
        micPerFrame.principalCurv3d[i].copyTo(principalCurv3d[i]);
        principalCurv3d[i]/=255;
    }

    segScore3d.resize(micPerFrame.time_points);
    for(size_t i=0;i<micPerFrame.time_points;i++){
        micPerFrame.segScore3d[i].copyTo(segScore3d[i]);
    }

    data_rows_cols_slices=micPerFrame.data_rows_cols_slices;
    time_points=micPerFrame.time_points;
    micPerFrame.normalized_data4d.copyTo(normalized_data4d);

    number_cells=micPerFrame.number_cells;
    VoxelIdxList=micPerFrame.VoxelIdxList;
    avgItstyVec=micPerFrame.avgItstyVec;
    for(size_t i=0;i<avgItstyVec.size();i++){
        for(size_t j=0;j<avgItstyVec[i].size();j++){
            avgItstyVec[i][j]/=255;
        }
    }
    areaVec=micPerFrame.areaVec;
    ctrPt=micPerFrame.ctrPt;
    CoordinateRgs=micPerFrame.CoordinateRgs;


    effZinFmTwoEndIdx=Mat::zeros(micPerFrame.time_points,2,CV_8U);
    dtctMtScVec.resize(micPerFrame.time_points);
}

void MicTrackerRefineByTempCons::paraLearn(){

    checkZStacks();
    getDtctMtScVec();
    getFeature4linking();
}

void MicTrackerRefineByTempCons::segRefineByTempCons(){
    ////////////////////////////////////
    /// intial parameters
    ////////////////////////////////////
    refineInitial();

    fill (numCellResVec_cur.begin(),numCellResVec_cur.end(), nDtct);
    size_t max_detection = *max_element(numCellResVec_cur.begin(), numCellResVec_cur.end());


    size_t iOutIter=0;
    QTextStream out(stdout);

    while((max_detection>1)&&isChanged){
        iOutIter++;
        out <<"Outer Loop #"<<iOutIter<<"..."<< Qt::endl;
        ///////////////////////////////////////////////
        /// estimate cell number of each detection
        /// ///////////////////////////////////////////

        tempConsistSort();

        ///////////////////////////////////////////////
        /// refine segmentation based on estimation
        /// ///////////////////////////////////////////
        segRefineOperations();

    }

}

void MicTrackerRefineByTempCons::refineInitial(){
    number_cells_cur=number_cells;
    VoxelIdxList_cur=VoxelIdxList;
    avgItstyVec_cur=avgItstyVec;
    areaVec_cur=areaVec;
    ctrPt_cur=ctrPt;
    CoordinateRgs_cur=CoordinateRgs;



    isChanged=true;
    nDtct =accumulate(number_cells_cur.begin(), number_cells_cur.end(), 0);
    numCellResVec_cur.resize(nDtct);
}

void MicTrackerRefineByTempCons::tempConsistSort(){


    incldRelationship();

    minCostSolver.buildGraph((int)nDtct,dtctIncldLst_cur);
    numCellResVec_cur=minCostSolver.MinCostFlow:: my_min_cost();
}

void MicTrackerRefineByTempCons::incldRelationship(){
    // initialize vector: cellIdxsLst, dtctIncldLst_cur
    incldRelationshipInitial();

    // find intersected detection id in pre & post frame
    for(size_t t0=0;t0<number_cells_cur.size();t0++){
        for(size_t iSubId0=0;iSubId0<number_cells_cur[t0];iSubId0++){
            // get overall idx
            if(t0>0){
                // get pre frame overlapping cell
                childDtctSet(t0,iSubId0,true);
            }
            if(t0<number_cells_cur.size()-1){
                childDtctSet(t0,iSubId0,false);
            }
        }
    }
}

void MicTrackerRefineByTempCons::incldRelationshipInitial(){
    nDtct =accumulate(number_cells_cur.begin(), number_cells_cur.end(), 0);

    cellIdxsLst.clear();
    size_t cnt=0;
    for(size_t t0=0;t0<number_cells_cur.size();t0++){
        vector<size_t> temp;
        for(size_t iSubId0=0;iSubId0<number_cells_cur[t0];iSubId0++){
            temp.push_back(cnt);
            cnt++;
        }
        cellIdxsLst.push_back(temp);
    }

    dtctIncldLst_cur.clear();
    dtctIncldLst_cur.resize(2);
    for(size_t i=0;i<2;i++){
        dtctIncldLst_cur[i].resize(nDtct);
    }
}

void MicTrackerRefineByTempCons::childDtctSet(size_t t0,size_t iSubId0,bool flagPre){



    ///////////////////////////////////////////////////////////
    /// initialize
    ///////////////////////////////////////////////////////////
    size_t iDtct0=cellIdxsLst[t0][iSubId0];
    vector<size_t> iDtct0_vxLst=VoxelIdxList_cur[t0][iSubId0];

    float mtSc;
    size_t t_tgt;
    size_t drctId;
    if(flagPre){
        t_tgt=t0-1;
        drctId=0;
    }else{
        t_tgt=t0+1;
        drctId=1;
    }

    mtSc=dtctMtScVec[t0][iSubId0][drctId];

    ///////////////////////////////////////////////////////////
    /// if is moving cell, no child set and return
    ///////////////////////////////////////////////////////////
    float p=boost::math::gamma_q(nullMtScPhat[0],mtSc/nullMtScPhat[1]);
    if(p<pMtScThre){
        dtctIncldLst_cur[drctId][iDtct0].clear();
        return;
    }

    ///////////////////////////////////////////////////////////
    /// find all overlapping cells as candidate
    ///////////////////////////////////////////////////////////



    vector<vector<size_t>> tgtIdxList=VoxelIdxList_cur[t_tgt];

    vector<size_t> temp_res;
    vector<size_t> candDtctIds;
    vector<size_t> candDtctVlm;
    vector<size_t> candDtctOvlp;

    for(size_t i=0;i<tgtIdxList.size();i++){

        getCommonElementSorted(iDtct0_vxLst,tgtIdxList[i],temp_res);

        if(temp_res.size()>0){
            candDtctIds.push_back(i);
            candDtctVlm.push_back(areaVec_cur[t_tgt][i]);
            candDtctOvlp.push_back(temp_res.size());
        }



    }


    ///////////////////////////////////////////////////////////
    /// iteratively remove poor candidates based on IOU score
    ///////////////////////////////////////////////////////////


    size_t N=candDtctIds.size();
    size_t N_remain=N;
    vector<bool> goodCandi(N,true);
    vector<bool> goodCandi_temp(N,true);
    vector<float> tempOvlpScVec(N,0);
    size_t ovlp_sum_temp=0;
    size_t all_sum_temp=0;
    size_t vlm0=iDtct0_vxLst.size();



    // get original IOU
    ovlp_sum_temp=0;
    all_sum_temp=0;
    for(size_t k=0;k<N;k++){
        ovlp_sum_temp+=candDtctOvlp[k];
        all_sum_temp+=candDtctVlm[k];
    }
    float ovlpSc=((float)ovlp_sum_temp)/(all_sum_temp+vlm0-ovlp_sum_temp);

    for(size_t i=0;i<N;i++){
        // calculate the IOU score after remove jth candidate
        for(size_t j=0;j<N;j++){
            if(goodCandi[j]){
                ovlp_sum_temp=0;
                all_sum_temp=0;

                goodCandi_temp=goodCandi;
                goodCandi_temp[j]=false;
                for(size_t k=0;k<N;k++){
                    if(goodCandi_temp[k]){
                        ovlp_sum_temp+=candDtctOvlp[k];
                        all_sum_temp+=candDtctVlm[k];
                    }
                }
                tempOvlpScVec[j]=((float)ovlp_sum_temp)/
                        (all_sum_temp+vlm0-ovlp_sum_temp);
            }else{
                tempOvlpScVec[j]=0;
            }

            int maxElementIndex = max_element(tempOvlpScVec.begin(),tempOvlpScVec.end()) - tempOvlpScVec.begin();
            float maxElement = *max_element(tempOvlpScVec.begin(), tempOvlpScVec.end());

            if((maxElement>ovlpSc)||(N_remain>numCellResVec_cur[iDtct0])){
                ovlpSc=maxElement;
                goodCandi[maxElementIndex]=false;
                N_remain--;
            }else{
                break;
            }
        }
    }


    // send result to dtctIncldLst_cur
    for(size_t i=0;i<N;i++){
        if(goodCandi[i]){
            size_t inFramId=candDtctIds[i];
            size_t overallId=cellIdxsLst[t_tgt][inFramId];
            dtctIncldLst_cur[drctId][iDtct0].push_back(overallId);
        }
    }


}


void MicTrackerRefineByTempCons::getCommonElement(vector<size_t> vector1,vector<size_t> vector2,vector<size_t> &result){

    // Sort the vector (slow, need modification)
    sort(vector1.begin(), vector1.end());
    sort(vector2.begin(), vector2.end());


    // Initialise a vector
    // to store the common values
    // and an iterator
    // to traverse this vector
    vector<size_t> v(vector1.size() + vector2.size());
    vector<size_t>::iterator it, st;

    it = set_intersection(vector1.begin(),
                          vector1.end(),
                          vector2.begin(),
                          vector2.end(),
                          v.begin());
    result.clear();
    for (st = v.begin(); st != it; ++st){
        result.push_back(*st);
    }
}

void MicTrackerRefineByTempCons::getCommonElementSorted(vector<size_t> &vector1,vector<size_t> &vector2,vector<size_t> &result){

    // Initialise a vector
    // to store the common values
    // and an iterator
    // to traverse this vector
    vector<size_t> v(vector1.size() + vector2.size());
    typename vector<size_t>::iterator it, st;

    it = set_intersection(vector1.begin(),
                          vector1.end(),
                          vector2.begin(),
                          vector2.end(),
                          v.begin());

    result.clear();
    for (st = v.begin(); st != it; ++st){
        result.push_back(*st);
    }
}

void MicTrackerRefineByTempCons::getCommonElementFast(vector<size_t> vector1,vector<size_t> vector2,vector<size_t> &result){

    // find max element
    size_t maxElement1 = *max_element(vector1.begin(), vector1.end());
    size_t maxElement2 = *max_element(vector2.begin(), vector2.end());
    size_t maxElement=max(maxElement1,maxElement2);
    size_t minElement1 = *min_element(vector1.begin(), vector1.end());
    size_t minElement2 = *min_element(vector2.begin(), vector2.end());
    size_t minElement=min(minElement1,minElement2);


    // build map
    size_t N=maxElement-minElement+1;
    vector<bool> ref(N,false);
//    bool ref[N]={false};
    for(size_t i=0;i<vector1.size();i++){
        ref[vector1[i]-minElement]=true;
    }
    for(size_t i=0;i<vector2.size();i++){
        if(ref[vector2[i]-minElement]){
            result.push_back(i);
        }
    }
}

void MicTrackerRefineByTempCons::checkZStacks(){
    double tmp_min, tmp_max;
    long sz_single_frame =
            data_rows_cols_slices[0]
            *data_rows_cols_slices[1]
            *data_rows_cols_slices[2];
    long sz_single_slice =
            data_rows_cols_slices[0]
            *data_rows_cols_slices[1];
    for(size_t tt=0;tt<time_points;tt++){
        for(size_t zz=0;zz<data_rows_cols_slices[2];zz+=3){
            unsigned char *ind = (unsigned char*)normalized_data4d.data
                    +sz_single_frame*tt+sz_single_slice*zz; // sub-matrix pointer
            Mat single_frame = Mat(2, normalized_data4d.size, CV_8U, ind);
            minMaxIdx(single_frame, &tmp_min, &tmp_max);
            if(tmp_max>tmp_min){
                effZinFmTwoEndIdx.at<unsigned char>(tt,0)=zz;
                break;
            }
        }
    }

    for(size_t tt=0;tt<time_points;tt++){
        for(size_t zz=data_rows_cols_slices[2]-1;zz<-1;zz-=3){
            unsigned char *ind = (unsigned char*)normalized_data4d.data
                    +sz_single_frame*tt+sz_single_slice*zz; // sub-matrix pointer
            Mat single_frame = Mat(2, normalized_data4d.size, CV_8U, ind);
            minMaxIdx(single_frame, &tmp_min, &tmp_max);
            if(tmp_max>tmp_min){
                effZinFmTwoEndIdx.at<unsigned char>(tt,1)=zz;
                break;
            }
        }
    }
}

void MicTrackerRefineByTempCons::getDtctMtScVec(){
    long sz_single_frame =
            data_rows_cols_slices[0]
            *data_rows_cols_slices[1]
            *data_rows_cols_slices[2];
    long sz_single_slice =
            data_rows_cols_slices[0]
            *data_rows_cols_slices[1];

    size_t effVxBnd1[2]={0,0};
    size_t effVxBnd2[2]={0,0};
    vector<std::size_t> vxVec;
    size_t validSz=0;
    float sumIntensityDiffSQ=0;

    for(size_t tt=0;tt<(size_t)time_points;tt++){
        dtctMtScVec[tt].resize(number_cells[tt]);
        for(size_t cellCnt=0;cellCnt<number_cells[tt];cellCnt++){
            dtctMtScVec[tt][cellCnt].resize(2);
        }
    }


    // collect data for null and estimate distribution
    vector<float> data4null;
    for(size_t tt=0;tt<(size_t)time_points-1;tt++){

        unsigned char *ind = (unsigned char*)normalized_data4d.data+sz_single_frame*tt;
        Mat img1=Mat(3, normalized_data4d.size, CV_8U, ind);
        ind = (unsigned char*)normalized_data4d.data+sz_single_frame*(tt+1);
        Mat img2=Mat(3, normalized_data4d.size, CV_8U, ind);

        img1.convertTo(img1,CV_32F);
        img2.convertTo(img2,CV_32F);
        img1/=255;
        img2/=255;


        effVxBnd1[0]=(effZinFmTwoEndIdx.at<unsigned char>(tt,0))*sz_single_slice;
        effVxBnd1[1]=(effZinFmTwoEndIdx.at<unsigned char>(tt,1)+1)*sz_single_slice-1;
        effVxBnd2[0]=(effZinFmTwoEndIdx.at<unsigned char>(tt+1,0))*sz_single_slice;
        effVxBnd2[1]=(effZinFmTwoEndIdx.at<unsigned char>(tt+1,1)+1)*sz_single_slice-1;

        for(size_t cellCnt=0;cellCnt<number_cells[tt];cellCnt++){
            vxVec=VoxelIdxList[tt][cellCnt];

            validSz=0;
            sumIntensityDiffSQ=0;
            for(size_t pxlCnt=0;pxlCnt<vxVec.size();pxlCnt++){
                if(vxVec[pxlCnt]>=effVxBnd1[0]&&vxVec[pxlCnt]<=effVxBnd1[1]){
                    validSz++;
                    sumIntensityDiffSQ+=pow(img1.at<float>(vxVec[pxlCnt])-img2.at<float>(vxVec[pxlCnt]),2);
                }
            }

            if(validSz*2>vxVec.size()){
                dtctMtScVec[tt][cellCnt][1]=sqrt((sumIntensityDiffSQ)/validSz);
                data4null.push_back(dtctMtScVec[tt][cellCnt][1]);
            }
        }

        for(size_t cellCnt=0;cellCnt<number_cells[tt+1];cellCnt++){
            vxVec=VoxelIdxList[tt+1][cellCnt];
            validSz=0;
            for(size_t pxlCnt=0;pxlCnt<vxVec.size();pxlCnt++){
                if(vxVec[pxlCnt]>=effVxBnd2[0]&&vxVec[pxlCnt]<=effVxBnd2[1]){
                    validSz++;
                    sumIntensityDiffSQ+=pow(img1.at<float>(vxVec[pxlCnt])-img2.at<float>(vxVec[pxlCnt]),2);
                }
            }


            if(validSz*2>vxVec.size()){
                dtctMtScVec[tt+1][cellCnt][0]=sqrt((sumIntensityDiffSQ)/validSz);
                data4null.push_back(dtctMtScVec[tt+1][cellCnt][0]);
            }
        }
    }

    gamfit(data4null, nullMtScPhat);

    // collect data for alt and estimate distribution
    vector<float> data4alt;
    for(size_t tt=0;tt<(size_t)time_points;tt++){
        for(size_t cellCnt=0;cellCnt<number_cells[tt];cellCnt++){
            data4alt.push_back(avgItstyVec[tt][cellCnt]);
        }
    }
    gamfit(data4alt, altMtScPhat);
}

void MicTrackerRefineByTempCons::getFeature4linking(){
    vector<float> OvrlpRate_m;
    vector<float> nbIncldRate_m;
    vector<float> intDiff;
    vector<vector<float>> dsplsmt;
    dsplsmt.resize(3);
    vector<float> dist;
    vector<float> sizeChange;
    vector<float> sizeChangeRate;

    Mat cellIdMap2;
    size_t areal;
    vector<size_t> unqVals;


    for(size_t tt=0;tt<time_points-1;tt++){

        cellIdMap2=cell_label_maps[tt+1];

        for(size_t cellCnt=0;cellCnt<number_cells[tt];cellCnt++){
            //////////////////////////////////////////////////////////////////
            //                 1. all overlapping cells in the next frame
            //////////////////////////////////////////////////////////////////
            size_t area=areaVec[tt][cellCnt];
            vector<size_t> candIds;
            vector<size_t> ovlpAreas;
            vector<size_t> idAll;
            for(size_t vxCnt=0;vxCnt<VoxelIdxList[tt][cellCnt].size();vxCnt++){
                idAll.push_back(cellIdMap2.at<int>(VoxelIdxList[tt][cellCnt][vxCnt]));
            }
            CorrespondingCell(idAll,candIds,ovlpAreas);
            if(candIds.size()==0){continue;}
            //////////////////////////////////////////////////////////////////
            //                 2. max overlapping ratio cell and max include ratio
            //////////////////////////////////////////////////////////////////
            vector<float> candOvrlpRates;
            vector<float> candIncldRates;
            for(size_t cellCnt2=0;cellCnt2<candIds.size();cellCnt2++){
                size_t area2=areaVec[tt+1][cellCnt2];
                candOvrlpRates.push_back(ovlpAreas[cellCnt2]/((float) area+area2-ovlpAreas[cellCnt2]));
                candIncldRates.push_back(ovlpAreas[cellCnt2]/((float) area2));
            }
            int clstNbId =max_element(candOvrlpRates.begin(),candOvrlpRates.end()) - candOvrlpRates.begin();
            if(candOvrlpRates[clstNbId]>0){
                OvrlpRate_m.push_back(candOvrlpRates[clstNbId]);
            }

            candIncldRates[clstNbId]=0;
            float maxElement = *max_element(candIncldRates.begin(), candIncldRates.end());
            if(maxElement>0){
                nbIncldRate_m.push_back(maxElement);
            }
            //////////////////////////////////////////////////////////////////
            //                 3. intensity diff, dist, change size
            //////////////////////////////////////////////////////////////////

            float intDiffNow=(avgItstyVec[tt+1][candIds[clstNbId]]-avgItstyVec[tt][cellCnt]);
            intDiff.push_back(intDiffNow);
            vector<float> dsplsmt_temp={0,0,0};
            float disp_temp=0;
            float temp=0;
            for(int i=0;i<3;i++){
                temp=ctrPt[tt+1][candIds[clstNbId]][i]-ctrPt[tt][cellCnt][i];
                dsplsmt[i].push_back(temp);
                disp_temp+=pow(temp,2);
            }
            dist.push_back(sqrt(disp_temp));
            float vlm1=areaVec[tt][cellCnt];
            float vlm2=areaVec[tt+1][candIds[clstNbId]];
            sizeChange.push_back(vlm2-vlm1);
            sizeChangeRate.push_back((vlm2-vlm1)/((vlm2+vlm1)/2));
        }
    }

    //////////////////////////////////////////////////////////////////
    //                 4. estimate distribution
    //////////////////////////////////////////////////////////////////
    intDiff_sigma=stdev ( intDiff , intDiff_mu);
    sizeChange_sigma=stdev ( sizeChange, sizeChange_mu);
    sizeChangeRatio_sigma=stdev ( sizeChangeRate , sizeChangeRatio_mu);
    gamfit(dist,distPhat);
    for(int i=0;i<3;i++){
        TruncGaussFit(dsplsmt[i],DsplcmtMu3d[i], DsplcmtSigma3d[i]);
    }
    invOvrlpRtPhat=avg(OvrlpRate_m);
    nbIncldRtPhat=avg(nbIncldRate_m);
}


void MicTrackerRefineByTempCons::CorrespondingCell(vector<size_t> vxLst,
                                                   vector<size_t> &candIds,vector<size_t> &ovlpAreas){
    sort(vxLst.begin(), vxLst.end());
    size_t lastVx=0;
    for(size_t cnt=0;cnt<vxLst.size();cnt++){
        if(vxLst[cnt]==0){continue;}
        if(vxLst[cnt]!=lastVx){
            lastVx=vxLst[cnt];
            candIds.push_back(lastVx-1);
            ovlpAreas.push_back(1);
        }
        else
        {
            ovlpAreas[ovlpAreas.size()-1]++;
        }
    }
}

void MicTrackerRefineByTempCons::TruncGaussFit(vector<float> x, float &muhat, float &sigmahat){
    float x_mu0;
    float x_sigma0=stdev (x, x_mu0);
    sort(x.begin(), x.end());
    float bd_low=x_mu0-x_sigma0;
    float bd_up=x_mu0+x_sigma0;

    size_t ind_low=findVecBoundary(x,bd_low);
    size_t ind_up=findVecBoundary(x,bd_up);
    while(ind_low==ind_up){
        bd_low-=x_sigma0;
        bd_up+=x_sigma0;
        ind_low=findVecBoundary(x,bd_low);
        ind_up=findVecBoundary(x,bd_up);
    }

    vector<float> x_trunc;
    for(size_t i=ind_low;i<ind_up;i++){
        x_trunc.push_back(x[i]);
    }

    sigmahat=stdev (x, muhat);

    float gMu=gradiantTruncMu(x_trunc,muhat, sigmahat, bd_up, bd_low);
    float gSig=gradiantTruncSigma(x_trunc,muhat,sigmahat, bd_up, bd_low);
    float thre=1e-5;

    while(abs(gMu)>abs(thre*muhat)||abs(gSig)>thre*sigmahat){
        muhat+=gMu;
        sigmahat+=gSig;
        gMu=gradiantTruncMu(x_trunc,muhat, sigmahat, bd_up, bd_low);
        gSig=gradiantTruncSigma(x_trunc,muhat, sigmahat, bd_up, bd_low);
    }
}

float MicTrackerRefineByTempCons::gradiantTruncSigma(vector<float> x, float mu, float sigma, float alpha, float beta){
    float temp1_1=0;
    for(int i=0;i<x.size();i++){
        temp1_1+=pow((x[i]-mu),2)/pow(sigma,3);
    }
    temp1_1/=x.size();
    float temp1=temp1_1-1/sigma;

    float alpha_n=(alpha-mu)/sigma;
    float beta_n=(beta-mu)/sigma;
    float temp2_1=normal_pdf(alpha_n,0,1)*(alpha-mu);
    float temp2_2=normal_pdf(beta_n,0,1)*(beta-mu);
    float temp2_3=standardNormal_cdf(alpha_n);
    float temp2_4=standardNormal_cdf(beta_n);

    float temp2=(temp2_1-temp2_2)/(temp2_3-temp2_4);
    float temp3=temp2/pow(sigma,2);

    return temp1+temp3;
}

float MicTrackerRefineByTempCons::gradiantTruncMu(vector<float> x, float mu, float sigma, float alpha, float beta){
    float temp1=0;
    for(int i=0;i<x.size();i++){
        temp1+=(x[i]-mu)/pow(sigma,2);
    }
    temp1/=x.size();

    float alpha_n=(alpha-mu)/sigma;
    float beta_n=(beta-mu)/sigma;
    float temp2_1=normal_pdf(alpha_n,0,1);
    float temp2_2=normal_pdf(beta_n,0,1);
    float temp2_3=standardNormal_cdf(alpha_n);
    float temp2_4=standardNormal_cdf(beta_n);

    float temp2=(temp2_1-temp2_2)/(temp2_3-temp2_4);
    float temp3=temp2/sigma;

    return temp1+temp3;
}

float MicTrackerRefineByTempCons::normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

float MicTrackerRefineByTempCons::standardNormal_cdf(float x)
{
//    double const Pi=3.141592653589793238462643;
//    double L, K, w ;
//  /* constants */
//    double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
//    double const a4 = -1.821255978, a5 = 1.330274429;

//    L = fabs(x);
//    K = 1.0 / (1.0 + 0.2316419 * L);
//    w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

//    if (x < 0 ){
//    w= 1.0 - w;
//    }
//    return (float) w;
    return 0.5 * erfc(-x * M_SQRT1_2);
}

size_t MicTrackerRefineByTempCons::findVecBoundary(vector<float> x_sorted,float boundary){
    size_t i=0;
    for(i=0;i<x_sorted.size();i++){
        if(x_sorted[i]>boundary){break;}
    }
    return i;
}

void MicTrackerRefineByTempCons::gamfit(vector<float> x, vector<float> &parmhat){
//    histogram(x);
    float sumSQ=0;
    float sumx=0;
    float sumlogx=0;
    float xbar=1;
    for(size_t i=0;i<x.size();i++){
        sumx+=x[i];
    }
    float scalex=sumx/x.size();

    for(size_t i=0;i<x.size();i++){
        x[i]=x[i]/scalex;
        sumSQ+=pow(x[i]-xbar,2);
        sumlogx+=log(x[i]);
    }

    float s2=sumSQ/(x.size()-1);
    float ahat=pow(xbar,2)/s2;
    float bhat=s2/xbar;
    float cnst=sumlogx/x.size();
    float lower;
    float upper;
    if(gamfitLkeqn(ahat,cnst)>0){
        lower=0.5*ahat;
        upper=2*lower;
        while(gamfitLkeqn(lower,cnst)>0){
            upper=lower;
            lower=0.5*upper;
        }
    }else{
        lower=ahat;
        upper=2*lower;
        while(gamfitLkeqn(upper,cnst)<0){
            lower=upper;
            upper=2*lower;
        }
    }
    ahat=fzero_gamfitLkeqn(lower,upper,cnst);
    parmhat.resize(2);
    parmhat[0]=ahat;
    parmhat[1]=xbar/ahat*scalex;
}

float MicTrackerRefineByTempCons::gamfitLkeqn(float a, float cnst){
    return -cnst-log(a)+boost::math::digamma(a);
}

float MicTrackerRefineByTempCons::fzero_gamfitLkeqn(float lower, float upper, float cnst){
    float in=(lower+upper)/2;
    float out=gamfitLkeqn(in,cnst);
    while(out!=0){
        out=gamfitLkeqn(in,cnst);
        if(out>0){upper=in;}
        if(out<0){lower=in;}
        in=(upper+lower)/2;
        if((upper-lower<1e-8)||(abs(out)<1e-6)){break;}
    }
    return in;
}

void MicTrackerRefineByTempCons::histogram(vector<float> x){
    int histSize = 50;
    float range[] = { 0, 0.1 }; //the upper boundary is exclusive
    const float* histRange[] = { range };
    bool uniform = true, accumulate = false;
    Mat hist;
    int n=x.size();
    float arr[n];
    for (int i = 0; i < n; i++) {
        arr[i] = x[i];
    }
    Mat m(1, n, CV_32FC1, arr);

    calcHist( &m, 1, 0, Mat(), hist, 1, &histSize, histRange, uniform, accumulate );

    hist*=20;

    int hist_w = 1000, hist_h = 1000;
    int bin_w = cvRound( (double) hist_w/histSize );
    Mat histImage( hist_h, hist_w, CV_8U, Scalar( 0,0,0) );


    for( int i = 1; i < histSize; i++ )
    {
        line( histImage, Point( bin_w*(i-1), hist_h - cvRound(hist.at<float>(i-1)) ),
              Point( bin_w*(i), hist_h - cvRound(hist.at<float>(i)) ),
              Scalar( 255, 0, 0), 2, 8, 0  );
    }
    imshow("calcHist Demo", histImage );
    waitKey();
}

float MicTrackerRefineByTempCons::avg ( vector<float> v )
{
        float sum = std::accumulate(v.begin(), v.end(), 0.0);
        float mean = sum / v.size();;
        return mean;
}

float MicTrackerRefineByTempCons::stdev ( vector<float> v , float &mean )
{
    mean = avg(v);

    std::vector<float> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   std::bind2nd(std::minus<float>(), mean));
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / v.size());
    return stdev;
}


void MicTrackerRefineByTempCons::segRefineOperations()
{
    micTracker.refineBytemporalInfo_loop2(numCellResVec_cur,dtctIncldLst_cur);
    updateRefineResult();
}

void MicTrackerRefineByTempCons::updateRefineResult()
{

    // info from micPerFrame (this is deep copy of mat!!)

    if(number_cells_cur==micTracker.number_cells){
        isChanged=false;
        return;
    }else{
        number_cells_cur=micTracker.number_cells;
    }


    cell_label_maps.resize(micTracker.time_points);
    for(size_t i=0;i<(size_t) micTracker.time_points;i++){
        micTracker.cell_label_maps[i].copyTo(cell_label_maps[i]);
    }

    VoxelIdxList_cur=micTracker.VoxelIdxList;
    avgItstyVec_cur=micTracker.avgItstyVec;
    for(size_t i=0;i<avgItstyVec_cur.size();i++){
        for(size_t j=0;j<avgItstyVec_cur[i].size();j++){
            avgItstyVec_cur[i][j]/=255;
        }
    }

    areaVec_cur=micTracker.areaVec;
    ctrPt_cur=micTracker.ctrPt;
    CoordinateRgs_cur=micTracker.CoordinateRgs;
    updateDtctMtScVec();
}

void MicTrackerRefineByTempCons::updateDtctMtScVec(){
    long sz_single_frame =
            data_rows_cols_slices[0]
            *data_rows_cols_slices[1]
            *data_rows_cols_slices[2];
    long sz_single_slice =
            data_rows_cols_slices[0]
            *data_rows_cols_slices[1];

    size_t effVxBnd1[2]={0,0};
    size_t effVxBnd2[2]={0,0};
    vector<std::size_t> vxVec;
    size_t validSz=0;
    float sumIntensityDiffSQ=0;

    for(size_t tt=0;tt<(size_t)time_points;tt++){
        dtctMtScVec[tt].resize(number_cells_cur[tt]);
        for(size_t cellCnt=0;cellCnt<number_cells_cur[tt];cellCnt++){
            dtctMtScVec[tt][cellCnt].resize(2);
        }
    }


    // collect data for null and estimate distribution
    for(size_t tt=0;tt<(size_t)time_points-1;tt++){

        unsigned char *ind = (unsigned char*)normalized_data4d.data+sz_single_frame*tt;
        Mat img1=Mat(3, normalized_data4d.size, CV_8U, ind);
        ind = (unsigned char*)normalized_data4d.data+sz_single_frame*(tt+1);
        Mat img2=Mat(3, normalized_data4d.size, CV_8U, ind);

        img1.convertTo(img1,CV_32F);
        img2.convertTo(img2,CV_32F);
        img1/=255;
        img2/=255;


        effVxBnd1[0]=(effZinFmTwoEndIdx.at<unsigned char>(tt,0))*sz_single_slice;
        effVxBnd1[1]=(effZinFmTwoEndIdx.at<unsigned char>(tt,1)+1)*sz_single_slice-1;
        effVxBnd2[0]=(effZinFmTwoEndIdx.at<unsigned char>(tt+1,0))*sz_single_slice;
        effVxBnd2[1]=(effZinFmTwoEndIdx.at<unsigned char>(tt+1,1)+1)*sz_single_slice-1;

        for(size_t cellCnt=0;cellCnt<number_cells_cur[tt];cellCnt++){
            vxVec=VoxelIdxList_cur[tt][cellCnt];

            validSz=0;
            sumIntensityDiffSQ=0;
            for(size_t pxlCnt=0;pxlCnt<vxVec.size();pxlCnt++){
                if(vxVec[pxlCnt]>=effVxBnd1[0]&&vxVec[pxlCnt]<=effVxBnd1[1]){
                    validSz++;
                    sumIntensityDiffSQ+=pow(img1.at<float>(vxVec[pxlCnt])-img2.at<float>(vxVec[pxlCnt]),2);
                }
            }

            if(validSz*2>vxVec.size()){
                dtctMtScVec[tt][cellCnt][1]=sqrt((sumIntensityDiffSQ)/validSz);
            }
        }

        for(size_t cellCnt=0;cellCnt<number_cells_cur[tt+1];cellCnt++){
            vxVec=VoxelIdxList_cur[tt+1][cellCnt];
            validSz=0;
            for(size_t pxlCnt=0;pxlCnt<vxVec.size();pxlCnt++){
                if(vxVec[pxlCnt]>=effVxBnd2[0]&&vxVec[pxlCnt]<=effVxBnd2[1]){
                    validSz++;
                    sumIntensityDiffSQ+=pow(img1.at<float>(vxVec[pxlCnt])-img2.at<float>(vxVec[pxlCnt]),2);
                }
            }


            if(validSz*2>vxVec.size()){
                dtctMtScVec[tt+1][cellCnt][0]=sqrt((sumIntensityDiffSQ)/validSz);
            }
        }
    }
}
