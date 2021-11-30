#include "mincostflow.h"

using namespace cv;
using namespace std;

MinCostFlow::MinCostFlow()
{

}



void MinCostFlow:: shortest_paths(int n, int v0, vector<int>& d, vector<int>& p) {
    d.assign(n, INF);
    d[v0] = 0;
    vector<bool> inq(n, false);
    queue<int> q;
    q.push(v0);
    p.assign(n, -1);

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        inq[u] = false;
        for (int v : adj[u]) {
            if (capacity[u][v] > 0 && d[v] > d[u] + cost[u][v]) {
                d[v] = d[u] + cost[u][v];
                p[v] = u;
                if (!inq[v]) {
                    inq[v] = true;
                    q.push(v);
                }
            }
        }
    }
}

int MinCostFlow:: min_cost_flow(int N, vector<Edge> edges, int K, int s, int t) {

    adj.assign(N, vector<int>());
    cost.assign(N, vector<int>(N, 0));
    capacity.assign(N, vector<int>(N, 0));

    for (Edge e : edges) {
        adj[e.from].push_back(e.to);
        adj[e.to].push_back(e.from);
        cost[e.from][e.to] = e.cost;
        cost[e.to][e.from] = -e.cost;
        capacity[e.from][e.to] = e.capacity;
    }

    int flow = 0;
    int cost_all = 0;
    vector<int> d, p;
    while (flow < K) {
        shortest_paths(N, s, d, p);
        if (d[t] == INF)
            break;

        // find max flow on that path
        int f = K - flow;
        int cur = t;
        while (cur != s) {
            f = min(f, capacity[p[cur]][cur]);
            cur = p[cur];
        }

        // apply flow
        flow += f;
        cost_all += f * d[t];
        cur = t;
        while (cur != s) {
            capacity[p[cur]][cur] -= f;
            capacity[cur][p[cur]] += f;
            cur = p[cur];
        }
    }

    if (flow < K)
        return -1;
    else
        return cost_all;
}

int MinCostFlow:: min_cost(int N, vector<Edge> edges, int s, int t) {

    adj.assign(N, vector<int>());
    cost.assign(N, vector<int>(N, 0));
    capacity.assign(N, vector<int>(N, 0));
    flowMap.assign(N, vector<int>(N, 0));
    for (Edge e : edges) {
        adj[e.from].push_back(e.to);
        adj[e.to].push_back(e.from);
        cost[e.from][e.to] = e.cost;
        cost[e.to][e.from] = -e.cost;
        capacity[e.from][e.to] = e.capacity;
    }
    capacity_org=capacity;

    int flow = 0;
    int cost_all = 0;
    vector<int> d, p;
    while (flow < N) {
        shortest_paths(N, s, d, p);
        if (d[t]>0){
            break;
        }


        // find max flow on that path
        int f = INF - flow;
        int cur = t;
        while (cur != s) {
            f = min(f, capacity[p[cur]][cur]);
            cur = p[cur];
        }

        // apply flow
        flow += f;
        cost_all += f * d[t];
        cur = t;
        while (cur != s) {
            capacity[p[cur]][cur] -= f;
            capacity[cur][p[cur]] += f;
            cur = p[cur];
        }
    }

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            flowMap[i][j]=capacity_org[i][j]-capacity[i][j];
        }
    }

    return cost_all;
}

vector<size_t> MinCostFlow:: my_min_cost(){
    int s=mN-2;
    int t=mN-1;
    min_cost(mN, mEdges, s, t);

    vector<size_t> numCellResVec;
    vector<int> tempVec;
    size_t temp;
    for(int i=0;i<mDetectionNumber;i++){
        tempVec=flowMap[i*4];
        temp=0;
        for(size_t j=0;j<tempVec.size();j++){
            if(tempVec[j]>0){
                temp+=tempVec[j];
            }
        }
        numCellResVec.push_back(temp);
    }
    return numCellResVec;
}

void MinCostFlow::buildGraph(int N, vector<vector<vector<size_t> > > &dtctIncldLst_cur){
    // N : detection number
    //
    mDetectionNumber=N;
    mN=4*N+2;
    int s=4*N; // source
    int t=4*N+1; // sink
    mEdges.clear();
    Edge e;
    int alphaFP=3;
    int alphaFC=1;
    int beta=4;

    /////////////////////
    /// pre->unit
    /////////////////////
    for(int i=0;i<(4*N);i=i+4){
        e.from=i;
        e.to=i+1;
        e.cost=0;
        e.capacity=1;
        mEdges.push_back(e);
    }

    /////////////////////
    /// unit->post
    /////////////////////
    for(int i=0;i<(4*N);i=i+4){
        e.from=i+1;
        e.to=i+3;
        e.cost=-alphaFP;
        e.capacity=1;
        mEdges.push_back(e);
    }

    /////////////////////
    /// pre->extra
    /////////////////////
    for(int i=0;i<(4*N);i=i+4){
        e.from=i;
        e.to=i+2;
        e.cost=0;
        e.capacity=INF;
        mEdges.push_back(e);
    }

    /////////////////////
    /// extra->post
    /////////////////////
    for(int i=0;i<(4*N);i=i+4){
        e.from=i+2;
        e.to=i+3;
        e.cost=alphaFC;
        e.capacity=INF;
        mEdges.push_back(e);
    }

    /////////////////////
    /// source->pre
    /////////////////////
    for(int i=0;i<N;i++){
        e.from=s;
        e.to=i*4;
        if(dtctIncldLst_cur[0][i].empty()){
            e.cost=0;
        }
        else{
            e.cost=beta;
        }
        e.capacity=INF;
        mEdges.push_back(e);
    }

    /////////////////////
    /// post->sink
    /////////////////////
    for(int i=0;i<N;i++){
        e.from=i*4+3;
        e.to=t;
        if(dtctIncldLst_cur[1][i].empty()){
            e.cost=0;
        }
        else{
            e.cost=beta;
        }
        e.capacity=INF;
        mEdges.push_back(e);
    }

    /////////////////////
    /// post->pre
    /////////////////////
    int Nedge=dtctIncldLst_cur[1].size();
    for(int i=0;i<Nedge;i++){
        e.from=i*4+3;
        e.cost=0;
        e.capacity=INF;
        for(int j=0;j< (int) dtctIncldLst_cur[1][i].size();j++){
            e.to=dtctIncldLst_cur[1][i][j]*4;
            mEdges.push_back(e);
        }
    }
}


//void MinCostFlow:: shortest_paths_floatCost(int n, int v0, vector<float>& d) {
//    d.assign(n, INF);
//    d[v0] = 0;
//    vector<bool> inq(n, false);
//    vector<int> p;
//    queue<int> q;
//    q.push(v0);
//    p.assign(n, -1);

//    while (!q.empty()) {
//        int u = q.front();
//        q.pop();
//        inq[u] = false;
//        for (int v : adj[u]) {
//            if (capacity[u][v] > 0 && d[v] > d[u] + cost_f[u][v]) {
//                d[v] = d[u] + cost_f[u][v];
//                p[v] = u;
//                if (!inq[v]) {
//                    inq[v] = true;
//                    q.push(v);
//                }
//            }
//        }
//    }
//}

//void MinCostFlow:: mincost_forwaterShed_buildGraph(Mat &scoreMap, vector<size_t> &seedVx,Mat &Mask) {

//    const float MINCOST=0.0000001;

//    size_t sz_y=scoreMap.size.p[0];
//    size_t sz_x=scoreMap.size.p[1];
//    size_t sz_z=scoreMap.size.p[2];
//    size_t sz_slice=sz_y*sz_x;
//    int N=sz_slice*sz_z;
//    size_t n;
//    size_t n2;
//    float c;

//    adj.clear();
//    cost_f.clear();
//    capacity.clear();

//    adj.assign(N, vector<int>());
//    cost_f.assign(N, vector<float>(N, 0));
//    capacity.assign(N, vector<int>(N, 0));

////    vol_sub2ind(vxc, yc, xc, zc, idMap->size);
//    for(size_t z=0;z<sz_z;z++){
//        for(size_t y=0;y<sz_y;y++){
//            for(size_t x=0;x<sz_x;x++){
//                vol_sub2ind(n, y, x, z, scoreMap.size);
//                if(x>0){
//                    vol_sub2ind(n2, y, x-1, z, scoreMap.size);
//                    if(Mask.at<int>(n2)>0){
//                        adj[n].push_back(n2);
//                        c=scoreMap.at<float>(n)-scoreMap.at<float>(n2);
//                        cost_f[n][n2] = max(c,MINCOST);
//                        capacity[n][n2] = 1;
//                    }
//                }
//                if(x<sz_x-1){
//                    vol_sub2ind(n2, y, x+1, z, scoreMap.size);
//                    if(Mask.at<int>(n2)>0){
//                        adj[n].push_back(n2);
//                        c=scoreMap.at<float>(n)-scoreMap.at<float>(n2);
//                        cost_f[n][n2] = max(c,MINCOST);
//                        capacity[n][n2] = 1;
//                    }
//                }
//                if(y>0){
//                    vol_sub2ind(n2, y-1, x, z, scoreMap.size);
//                    if(Mask.at<int>(n2)>0){
//                        adj[n].push_back(n2);
//                        c=scoreMap.at<float>(n)-scoreMap.at<float>(n2);
//                        cost_f[n][n2] = max(c,MINCOST);
//                        capacity[n][n2] = 1;
//                    }
//                }
//                if(y<sz_y-1){
//                    vol_sub2ind(n2, y+1, x, z, scoreMap.size);
//                    if(Mask.at<int>(n2)>0){
//                        adj[n].push_back(n2);
//                        c=scoreMap.at<float>(n)-scoreMap.at<float>(n2);
//                        cost_f[n][n2] = max(c,MINCOST);
//                        capacity[n][n2] = 1;
//                    }
//                }
//                if(z>0){
//                    vol_sub2ind(n2, y, x, z-1, scoreMap.size);
//                    if(Mask.at<int>(n2)>0){
//                        adj[n].push_back(n2);
//                        c=scoreMap.at<float>(n)-scoreMap.at<float>(n2);
//                        cost_f[n][n2] = max(c,MINCOST);
//                        capacity[n][n2] = 1;
//                    }
//                }
//                if(z<sz_z-1){
//                    vol_sub2ind(n2, y, x, z+1, scoreMap.size);
//                    if(Mask.at<int>(n2)>0){
//                        adj[n].push_back(n2);
//                        c=scoreMap.at<float>(n)-scoreMap.at<float>(n2);
//                        cost_f[n][n2] = max(c,MINCOST);
//                        capacity[n][n2] = 1;
//                    }
//                }
//            }
//        }
//    }
//    size_t id;
//    for(size_t i=0;i<seedVx.size();i++){
//        id=seedVx[i];
//        for(size_t j=0;j<(size_t)N;j++){
//            if(capacity[j][id]>MINCOST){
//                capacity[j][id]=MINCOST;
//            }
//        }
//    }

//}

//void MinCostFlow::mincost_forwaterShed(Mat &scoreMap, vector<vector<size_t>> seedLst,Mat Mask,Mat &resMap){

//    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

//    size_t sz_y=scoreMap.size.p[0];
//    size_t sz_x=scoreMap.size.p[1];
//    size_t sz_z=scoreMap.size.p[2];
//    size_t sz_slice=sz_y*sz_x;
//    int N=sz_slice*sz_z;
//    int SeedN=seedLst.size();

////    const float BoundaryValue=-10000;
////    for(size_t i=0;i<N;i++){
////        if(Mask.at<int>(i)==0){
////            scoreMap.at<float>(i)=BoundaryValue;
////        }
////    }

//    vector<size_t> sLst;
//    vector<size_t> seedVx;
//    for(size_t i=0;i<(size_t)SeedN;i++){
//        for(size_t j=0;j<seedLst[i].size();j++){
//            if(j==0){
//                sLst.push_back(seedLst[i][j]); // use the first voxel as source
//            }
//            seedVx.push_back(seedLst[i][j]); // all the seed voxels
//        }
//    }

//    // 1. build graph

//    mincost_forwaterShed_buildGraph(scoreMap,seedVx,Mask);

//    // 2. min cost
//    vector<vector<float>> d;
//    d.resize(sLst.size());
//    for(size_t i=0;i<d.size();i++){
//        shortest_paths_floatCost(N,sLst[i],d[i]);
//    }


//    // 3. give watershed result
//    resMap=Mat(scoreMap.dims, scoreMap.size, CV_32S, Scalar(0));
//    float minSc;

//    for(size_t i=0;i<(size_t)N;i++){
//        minSc=INF;
//        for(size_t j=0;j<sLst.size();j++){
//            if(d[j][i]<minSc){
//                minSc=d[j][i];
//                resMap.at<int>(i)=j+1;
//            }
//        }
//    }

//    chrono::steady_clock::time_point end = chrono::steady_clock::now();
//    qInfo("---------------- watershed time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);


//}
