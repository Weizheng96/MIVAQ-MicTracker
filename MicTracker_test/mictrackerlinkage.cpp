#include "mictrackerlinkage.h"

/////////////////////////////////////
/// functions for test
/////////////////////////////////////
void MicTrackerLinkage::test()
{
    int a=0;
}

/////////////////////////////////////
/// functions about intitialization
/////////////////////////////////////

MicTrackerLinkage::MicTrackerLinkage()
{

}

MicTrackerLinkage::MicTrackerLinkage(MicTrackerRefineByTempCons refineResult)
{
    Initialize(refineResult);
    getIntegratedMotionModelMCFlinking();
    getTraces();
}

void MicTrackerLinkage::Initialize(MicTrackerRefineByTempCons refineResult)
{
    //////////////////////////////////////
    /// Get parameters from refine result
    //////////////////////////////////////

    // label map
    cell_label_maps=refineResult.cell_label_maps;

    // video size
    data_rows_cols_slices=refineResult.data_rows_cols_slices;
    time_points=refineResult.time_points;

    // current cell info
    number_cells=refineResult.number_cells_cur;
    VoxelIdxList=refineResult.VoxelIdxList_cur; // time -> cell -> pixel->idx
    avgItstyVec=refineResult.avgItstyVec_cur; // time -> cell -> intensity
    areaVec=refineResult.areaVec_cur; // time -> cell -> intensity
    ctrPt=refineResult.ctrPt_cur; // time -> cell -> Point Center (y,x,z)
    CoordinateRgs=refineResult.CoordinateRgs_cur; // time -> cell -> boundary (min y, max y, min x, max x, min z, max z);

    // set index LUT
    setIndLUT();

    // estimated null hypothesis
    nullMtScPhat=refineResult.nullMtScPhat;
    altMtScPhat=refineResult.altMtScPhat;

    apprPara.initialize(refineResult);
    posParaStatic.Initialize(refineResult);

    // score based on null hypothesis
    dtctMtScVec=refineResult.dtctMtScVec;

    /////////////////////////////////////////
    /// Set parameters for linkage
    /////////////////////////////////////////
    xStaticPara =posParaStatic;
    xRpMtPara=MotionParameters().static2motile(xStaticPara);
}

void MicTrackerLinkage::setIndLUT()
{
    localId2GlobalId.clear();
    size_t cnt=0;
    localId2GlobalId.resize(time_points);
    for(size_t time=0;time<(size_t)time_points;time++){
        localId2GlobalId[time].resize(number_cells[time]);
        for(size_t localId=0;localId<number_cells[time];localId++){
            localId2GlobalId[time][localId]=cnt;
            ++cnt;
        }
    }

    cellNumInAllFrames=cnt;

    // initialize global to local id
    globalId2LocalId.clear();
    globalId2LocalId.resize(cnt);
    cnt=0;
    for(size_t time=0;time<(size_t)time_points;time++){
        for(size_t localId=0;localId<number_cells[time];localId++){
            globalId2LocalId[cnt].push_back(time);
            globalId2LocalId[cnt].push_back(localId);
            ++cnt;
        }
    }
}

/////////////////////////////////////
/// functions for MCF
/////////////////////////////////////

void MicTrackerLinkage::getIntegratedMotionModelMCFlinking()
{
    // get posterior probability
    getPosteriorProbability();
    buildGraph();
    mSolver.edges2graph(mSolver.edges,nNode,mSolver.mGraph);
    mSolver.min_cost(mSolver.mGraph);

}

void MicTrackerLinkage::getPosteriorProbability()
{
    // data is vector<vector<vector<float>>> dtctMtScVec;
    // get rpdMtRltdDtctPostP1
    rpdMtRltdDtctPostP1=dtctMtScVec;
    float temp;
    float temp0;
    float temp1;
    for(vector<vector<float>> & v1:rpdMtRltdDtctPostP1)
    {
        for(vector<float> & v2:v1)
        {
            for(float &v3:v2){
                if(v3<=0)
                {
                    v3=0;
                }
                else
                {
                    temp=v3;
                    // temp1 is the pdf of intensity difference motile cell
                    temp0=gammapdf(temp,nullMtScPhat[0], nullMtScPhat[1]);
                    // temp2 is the pdf of intensity difference static cell
                    temp1=gammapdf(temp,altMtScPhat[0], altMtScPhat[1]);
                    // postPvec1
                    v3=temp1/(temp0+temp1);
                }
            }
        }
    }
    // get dtctPvalVec
    dtctPvalVec=dtctMtScVec;
    for(vector<vector<float>> & v1:dtctPvalVec)
    {
        for(vector<float> & v2:v1)
        {
            for(float &v3:v2){
                v3=gammacdf(v3,nullMtScPhat[0], nullMtScPhat[1]);
            }
        }
    }
    // auxiliary parameters
    para.initialize(data_rows_cols_slices,time_points,cellNumInAllFrames,xRpMtPara,xStaticPara,apprPara);
}



void MicTrackerLinkage::buildGraph()
{
    // get the dummy scores
    double pInvOvrlpRtRefStatic=1;
    double pNbIncldRtRefStatic=1;

    double pInvOvrlpRtStaticBase=1/expcdf<double>(1,xStaticPara.invOvrlpRtPhat,false);
    double pNbIncldRtStaticBase=1/expcdf<double>(1,xStaticPara.nbIncldRtPhat,false);

    nNode=para.nCellRep*2+2;

    auto &edges=mSolver.edges;
    edges.resize(para.nCellRep*3);

    size_t source=nNode-2;
    size_t sink=nNode-1;
    para.source=source;
    para.sink=sink;

    // p(yi)
    double c_yi=log(para.FPdtctRate/(1-para.FPdtctRate));

    for(size_t i=0;i<para.nCellRep;++i)
    {
        edges[i].from=i*2;
        edges[i].to=i*2+1;
        edges[i].cost=c_yi;
        edges[i].capacity=1;
    }

    // P_enter, P_exit
    double c_enter=-c_yi/2-10e-10;
    double c_exit=c_enter;

    for(int i=0;i<para.nCellRep;++i)
    {
        edges[i+para.nCellRep].from=source;
        edges[i+para.nCellRep].to=i*2;
        edges[i+para.nCellRep].cost=c_enter;
        edges[i+para.nCellRep].capacity=1;

        edges[i+2*para.nCellRep].from=i*2+1;
        edges[i+2*para.nCellRep].to=sink;
        edges[i+2*para.nCellRep].cost=c_exit;
        edges[i+2*para.nCellRep].capacity=1;
    }

    // P_xij
    int countEdges=0;
    for(size_t t1=0;t1<para.nnT;++t1)
    {
        qDebug()<<"t: "<<t1<<Qt::endl;
        size_t nCnt=number_cells[t1];

        for(size_t iCnt1=0;iCnt1<nCnt;++iCnt1)
        {
            vector<size_t> pxVec1=VoxelIdxList[t1][iCnt1];
            vector<float> pos1=ctrPt[t1][iCnt1];
            vector<size_t> coord1=CoordinateRgs[t1][iCnt1];
            double intst1=avgItstyVec[t1][iCnt1];
            int vlm1=areaVec[t1][iCnt1];
            double rpMtP1=rpdMtRltdDtctPostP1[t1][iCnt1][1];

            for(size_t kJump=0;kJump<=para.maxJump;++kJump)
            {
                size_t t2=t1+1+kJump;
                if(t2>=para.nnT)
                {
                    break;
                }

                /// get valid candidate

                // get candidate distance
                vector<double> candDists(ctrPt[t2].size(),0);
                for(size_t i=0;i<ctrPt[t2].size();++i) // t2 is the  time, i is the cell
                {
                    vector<float> &ctrTemp=ctrPt[t2][i];
                    double temp=0;
                    for(int j=0;j<3;++j)
                    {
                        temp+=pow(ctrTemp[j]-pos1[j],2);
                    }
                    candDists[i]=sqrt(temp);
                }

                // get p value of candidate distance by motile and static hypothesis
                vector<double> pCtrPos1=candDists;
                for(double &v:pCtrPos1)
                {
                    v=gammapdf<double>(v,xRpMtPara.distPhat[0],xRpMtPara.distPhat[1]);
                }
                vector<double> pCtrPos0=candDists;
                for(double &v:pCtrPos0)
                {
                    v=gammapdf<double>(v,xStaticPara.distPhat[0],xStaticPara.distPhat[1]);
                }

                // get the probability of belong to motile cell
                vector<double> rpMtEdgeP1(candDists.size(),0);
                for(size_t i=0;i<rpMtEdgeP1.size();++i)
                {
                    rpMtEdgeP1[i]=rpMtP1*rpdMtRltdDtctPostP1[t2][i][0];
                }
                // get center position overall pdf
                vector<double> pCtrPos(candDists.size(),0);
                for(size_t i=0;i<pCtrPos.size();++i)
                {
                    pCtrPos[i]=(1-rpMtEdgeP1[i])*pCtrPos0[i]+rpMtEdgeP1[i]*pCtrPos1[i];
                }
                // get the valid candidate
                vector<size_t> candIds;
                for(size_t i=0;i<pCtrPos.size();++i)
                {
                    if(pCtrPos[i]>=max(para.pLinkThr,((double)para.FPdtctRate)/(1-para.FPdtctRate)))
                    {
                        candIds.push_back(i);
                    }
                }

                // if no valid candidate, skip
                if(candIds.size()<=0)
                {
                    continue;
                }

                /// rapid motion ralated probability

                vector<double> pOvlpRt1(candIds.size(),0);
                vector<double> pNbIncldRt1(candIds.size(),0);

                if(kJump==0)
                {
                    // remove invalid candidate
                    vector<double> pCtrPos1_temp(candIds.size(),0);
                    for(size_t i=0;i<candIds.size();++i)
                    {
                        pCtrPos1_temp[i]=pCtrPos1[candIds[i]];
                    }
                    pCtrPos1.swap(pCtrPos1_temp);
                    fill(pOvlpRt1.begin(),pOvlpRt1.end(),pInvOvrlpRtRefStatic);
                    fill(pNbIncldRt1.begin(),pNbIncldRt1.end(),pNbIncldRtRefStatic);

                }
                else
                {
                    pCtrPos1.resize(candIds.size());
                    fill(pCtrPos1.begin(),pCtrPos1.end(),0);
                }

                /// rapid motion unrelated (static) probabilities
                // remove invalid candidate
                vector<double> pCtrPos0_temp(candIds.size(),0);
                for(size_t i=0;i<candIds.size();++i)
                {
                    pCtrPos0_temp[i]=pCtrPos0[candIds[i]];
                }
                pCtrPos0.swap(pCtrPos0_temp);

                // find candidate not overlapping to remove
                vector<bool> cands2rm0(candIds.size(),true);
                for(size_t i=0;i<candIds.size();++i)
                {
                   if(CoordinateRgs[t2][candIds[i]][1] <coord1[0]){continue;};
                   if(CoordinateRgs[t2][candIds[i]][0] >coord1[1]){continue;};
                   if(CoordinateRgs[t2][candIds[i]][3] <coord1[2]){continue;};
                   if(CoordinateRgs[t2][candIds[i]][2] >coord1[3]){continue;};
                   if(CoordinateRgs[t2][candIds[i]][5] <coord1[4]){continue;};
                   if(CoordinateRgs[t2][candIds[i]][4] >coord1[5]){continue;};
                   cands2rm0[i]=false;
                }

                // Key feature: overlapRate
                vector<double> candOvlpRts(candIds.size(),0);
                vector<double> candOvlps(candIds.size(),0);
                vector<double> pOvlpRt0(candIds.size(),0);

                for(size_t iiDtct2=0;iiDtct2<candIds.size();++iiDtct2)
                {
                    if(cands2rm0[iiDtct2]){continue;}

                    int dtctId2=candIds[iiDtct2];
                    vector<size_t> pxVec2=VoxelIdxList[t2][dtctId2];
                    vector<size_t> pxVec_temp;
                    MicTrackerRefineByTempCons::getCommonElementSorted(pxVec1,pxVec2,pxVec_temp);
                    size_t ovlp=pxVec_temp.size();
                    candOvlpRts[iiDtct2]=ovlp/(pxVec1.size()+pxVec2.size()-ovlp);
                    candOvlps[iiDtct2]=ovlp;
                    if(ovlp==0)
                    {
                        cands2rm0[iiDtct2]=true;
                    }
                    else
                    {
                        pOvlpRt0[iiDtct2]=exppdf<double>(1-candOvlpRts[iiDtct2],xStaticPara.invOvrlpRtPhat)*pInvOvrlpRtStaticBase;
                    }
                }

                // key feature: 2nd neighbor inclusion (both sides)
                vector<double> candVlms(candIds.size(),0);
                for(size_t i=0;i<candIds.size();++i)
                {
                    candVlms[i]=areaVec[t2][candIds[i]];
                }
                vector<double> candIncldRts(candIds.size(),0);
                vector<double> pNbIncldRt0(candIds.size(),0);

                for(size_t iiDtct2=0;iiDtct2<candIds.size();++iiDtct2)
                {
                    if(cands2rm0[iiDtct2]){continue;}
                    int iCnt2=candIds[iiDtct2];
                    vector<size_t> pxVec2=VoxelIdxList[t2][iCnt2];
                    vector<float> pos2=ctrPt[t2][iCnt2];
                    vector<size_t> coord2=CoordinateRgs[t2][iCnt2];

                    // check possible 2nd neighbor (positive time order)
                    vector<double> candIncldRtPosDir(candIds.size(),0);
                    for(size_t i=0;i<candIds.size();++i)
                    {
                        if(i==iiDtct2){continue;}
                        candIncldRtPosDir[i]=candOvlps[i]/candVlms[i];
                    }
                    // check possible 2nd neighbor (negative time order)
                    vector<double> candPctrPos1(ctrPt[t1].size(),1);
                    vector<size_t> candIds1;
                    for(size_t i=0;i<candPctrPos1.size();++i)
                    {
                        for(int j=0;j<3;++j)
                        {
                            candPctrPos1[i]*=normalPDF<double>(ctrPt[t1][i][j],pos2[j],xStaticPara.DsplcmtSigma3d[j]*(kJump+1));
                        }
                        if(candPctrPos1[i]>0)
                        {
                            if(i!=iCnt1)
                            {
                                if(CoordinateRgs[t1][i][1] <coord2[0]){continue;};
                                if(CoordinateRgs[t1][i][0] >coord2[1]){continue;};
                                if(CoordinateRgs[t1][i][3] <coord2[2]){continue;};
                                if(CoordinateRgs[t1][i][2] >coord2[3]){continue;};
                                if(CoordinateRgs[t1][i][5] <coord2[4]){continue;};
                                if(CoordinateRgs[t1][i][4] >coord2[5]){continue;};

                                candIds1.push_back(i);
                            }
                        }
                    }

                    vector<double> candIncldRtNegDir(candIds1.size(),0);
                    for(size_t iiDtct1=0;iiDtct1<candIds1.size();++iiDtct1)
                    {
                        int candDtctId1=candIds1[iiDtct1];
                        vector<size_t> pxVec_temp;
                        MicTrackerRefineByTempCons::getCommonElementSorted(VoxelIdxList[t1][candDtctId1],pxVec2,pxVec_temp);
                        candIncldRtNegDir[iiDtct1]=pxVec_temp.size()/VoxelIdxList[t1][candDtctId1].size();
                    }


                    // remove non-overlapping
                    vector<double> candIncldRtNegDir_temp;
                    vector<size_t> candIds1_temp;
                    for(size_t iiDtct1=0;iiDtct1<candIds1.size();++iiDtct1)
                    {
                        if(candIncldRtNegDir[iiDtct2]==0)
                        {
                            candIncldRtNegDir_temp.push_back(candIncldRtNegDir[iiDtct1]);
                            candIds1_temp.push_back(candIds1[iiDtct1]);
                        }
                    }
                    candIncldRtNegDir.swap(candIncldRtNegDir_temp);
                    candIds1.swap(candIds1_temp);
                    if(candIncldRtNegDir.size()>0)
                    {
                        candIncldRts[iiDtct2]=max(vec_max(candIncldRtNegDir),vec_max(candIncldRtPosDir));
                    }
                    else
                    {
                        candIncldRts[iiDtct2]=vec_max(candIncldRtPosDir);
                    }


                }

                for(size_t i=0;i<candIncldRts.size();++i)
                {
                    if(!cands2rm0[i])
                    {
                        pNbIncldRt0[i]=exppdf<double>(candIncldRts[i],xStaticPara.nbIncldRtPhat)*pNbIncldRtStaticBase;
                    }
                }

                /// combined motion model
                rpMtEdgeP1.clear();
                vector<double> pSpace(candIds.size(),0);
                for(size_t i=0;i<candIds.size();++i)
                {
                    rpMtEdgeP1[i]=rpMtP1*rpdMtRltdDtctPostP1[t2][candIds[i]][0];
                    double pSpace0=pCtrPos0[i]*pOvlpRt0[i]*pNbIncldRt0[i];
                    double pSpace1=pCtrPos1[i]*pOvlpRt1[i]*pNbIncldRt1[i];
                    pSpace[i]=(1-rpMtEdgeP1[i])*pSpace0+rpMtEdgeP1[i]*pSpace1;
                }

                // Appearance features: intensity, sizeChangeRate
                vector<double> candDtctIntst(candIds.size(),0);
                vector<double> pIntst(candIds.size(),0);
                vector<double> candDtctVlm(candIds.size(),0);
                vector<double> pSize(candIds.size(),0);
                for(size_t i=0;i<candIds.size();++i)
                {
                    candDtctIntst[i]=avgItstyVec[t2][candIds[i]];
                    pIntst[i]=normalPDF<double>(candDtctIntst[i],intst1,apprPara.intDiff_sigma);

                    candDtctVlm[i]=areaVec[t2][candIds[i]];
                    pSize[i]=normalPDF<double>((candDtctVlm[i]-vlm1)/((candDtctVlm[i]+vlm1)/2),apprPara.sizeChangeRatio_mu,apprPara.sizeChangeRatio_sigma );
                }

                double pT=exppdf<double>(1+kJump,(para.maxJump+1)/2);

                // overall score
                vector<double> pLink(candIds.size(),0);
                vector<bool> boolActvLink(candIds.size(),false);
                for(size_t i=0;i<candIds.size();++i)
                {
                    pLink[i]=pIntst[i]*pSize[i]*pSpace[i]*pT;
                    if(pLink[i]>=max(para.pLinkThr,(double)para.FPdtctRate/(1-para.FPdtctRate)))
                    {
                        boolActvLink[i]=true;

                        NetworkFlowSolver<double,int>::Edge e;
                        e.from=localId2GlobalId[t1][iCnt1]*2+1;
                        e.to=localId2GlobalId[t2][candIds[i]]*2;
                        e.cost=-log(pLink[i]);
                        e.capacity=1;
                        edges.push_back(e);
                        ++countEdges;
                    }
                }

            }
        }
    }




}

void MicTrackerLinkage::getTraces()
{
    traces.clear();

    auto &flowMap=mSolver.mGraph.flowMap;
    size_t s=mSolver.mGraph.s;
    size_t t=mSolver.mGraph.t;
    size_t N=mSolver.mGraph.N;

    for(size_t cnt=0;cnt<N;++cnt)
    {
        if(flowMap[s][cnt]>0)
        {
            // find a new trace
            vector<QVector4D> tempTrace;
            getFullTrace(flowMap,cnt,t,tempTrace);
            traces.push_back(tempTrace);
        }
    }


}

void MicTrackerLinkage::getFullTrace(const auto &flowMap,const size_t &cur,const size_t &t,vector<QVector4D> &tempTrace)
{
    if(cur==t)
    {
        return;
    }
    else
    {
        size_t globalIdx=cur/2;
        const vector<size_t> &localIdx=globalId2LocalId[globalIdx];
        const vector<float> &ctr=ctrPt[localIdx[0]][localIdx[1]];
        tempTrace.emplace_back(QVector4D(localIdx[0],ctr[0],ctr[1],ctr[2]));

        size_t post=cur+1;
        size_t nxt=0;
        for(nxt=0;nxt<flowMap[post].size();++nxt)
        {
            if(flowMap[post][nxt]>0)
            {
                break;
            }
        }
        getFullTrace(flowMap,nxt,t,tempTrace);
        return;
    }
}
