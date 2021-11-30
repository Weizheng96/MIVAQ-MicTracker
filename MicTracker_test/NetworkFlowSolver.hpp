#ifndef NETWORKFLOWSOLVER_HPP
#define NETWORKFLOWSOLVER_HPP

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

template <class Tcost, class Tflow>
class NetworkFlowSolver
{
public:
    NetworkFlowSolver()
    {
        // empty constructor
    }

public:
    struct Edge
    {
        size_t from, to;
        Tflow capacity;
        Tcost cost;
    };
    struct Graph
    {
        size_t s;
        size_t t;
        size_t N;

        vector<vector<size_t>> adj;
        vector<vector<Tcost>> cost;
        vector<vector<Tflow>> capacity,capacity_org,flowMap;
        size_t mN;
        size_t mDetectionNumber;
    };
    const Tcost INFcost = 1e9;
    const Tflow INFflow=1e9;
    const size_t INF=1e9;

public:
    Graph mGraph;
    vector<Edge> edges;

public:
    // build graph based on edges
    void edges2graph(const vector<Edge> &edges,size_t N, Graph &mGraph)
    {
        mGraph.s=N-2;
        mGraph.t=N-1;
        mGraph.N=N;

        mGraph.adj.assign(N, vector<size_t>());
        mGraph.cost.assign(N, vector<Tcost>(N, 0));
        mGraph.capacity.assign(N, vector<Tflow>(N, 0));
        mGraph.flowMap.assign(N, vector<Tflow>(N, 0));
        for (Edge e : edges) {
            mGraph.adj[e.from].push_back(e.to);
            mGraph.adj[e.to].push_back(e.from);
            mGraph.cost[e.from][e.to] = e.cost;
            mGraph.cost[e.to][e.from] = -e.cost;
            mGraph.capacity[e.from][e.to] = e.capacity;
        }
        mGraph.capacity_org=mGraph.capacity;
    }

    // shortest path
    void shortest_paths(const Graph &mGraph, vector<Tcost>& d, vector<size_t>& p) {

        size_t n=mGraph.N;
        size_t v0=mGraph.s;

        d.assign(n, INFcost);
        d[v0] = 0;
        vector<bool> inq(n, false);
        queue<int> q;
        q.push(v0);
        p.assign(n, -1);

        while (!q.empty()) {
            size_t u = q.front();
            q.pop();
            inq[u] = false;
            for (size_t v : mGraph.adj[u]) {
                if (mGraph.capacity[u][v] > 0 && d[v] > d[u] + mGraph.cost[u][v]) {
                    d[v] = d[u] + mGraph.cost[u][v];
                    p[v] = u;
                    if (!inq[v]) {
                        inq[v] = true;
                        q.push(v);
                    }
                }
            }
        }
    }

    // max flow by SSP
    int min_cost(Graph &mGraph) {

        auto &s=mGraph.s;
        auto &t=mGraph.t;
        auto &N=mGraph.N;
        auto &adj=mGraph.adj;
        auto &cost=mGraph.cost;
        auto &capacity=mGraph.capacity;
        auto &capacity_org=mGraph.capacity_org;
        auto &flowMap=mGraph.flowMap;

        Tflow flow = 0;
        Tcost cost_all = 0;
        vector<Tcost> d;
        vector<size_t> p;

        while (true) {
            shortest_paths(mGraph,d,p);
            if (d[t]>0){
                break;
            }


            // find max flow on that path
            Tflow f = INFflow - flow;
            auto cur = t;
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

        for(size_t i=0;i<N;i++){
            for(int j=0;j<N;j++){
                flowMap[i][j]=capacity_org[i][j]-capacity[i][j];
            }
        }

        return cost_all;
    }

};


#endif // NETWORKFLOWSOLVER_HPP
