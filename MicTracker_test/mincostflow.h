#ifndef MINCOSTFLOW_H
#define MINCOSTFLOW_H

#include <vector>
#include <queue>
#include <numeric>

#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include "cellsegmentation/synquant_simple.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"

#include <chrono> // time elapsed

using namespace std;
using namespace cv;

class MinCostFlow
{
public:
    struct Edge
    {
        int from, to, capacity, cost;
    };
    const int INF = 1e9;

public:
    MinCostFlow();

    void shortest_paths(int n, int v0, vector<int>& d, vector<int>& p);
    int min_cost_flow(int N, vector<Edge> edges, int K, int s, int t);
    int min_cost(int N, vector<Edge> edges, int s, int t);

    void buildGraph(int N, vector<vector<vector<size_t> > > &dtctIncldLst_cur);
    vector<size_t> my_min_cost();

    // for watershed
//    void shortest_paths_floatCost(int n, int v0, vector<float>& d);
//    void mincost_forwaterShed_buildGraph(Mat &scoreMap, vector<size_t> &seedVx,Mat &Mask) ;
//    void mincost_forwaterShed(Mat &scoreMap, vector<vector<size_t>> seedLst,Mat Mask,Mat &resMap);



public:
    vector<vector<int>> adj, cost, capacity,capacity_org,flowMap;
    vector<vector<float>> cost_f;
    int mN;
    size_t mDetectionNumber;
    vector<Edge> mEdges;



};

#endif // MINCOSTFLOW_H
