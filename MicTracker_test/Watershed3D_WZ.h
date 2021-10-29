#ifndef WATERSHED3D_WZ_H
#define WATERSHED3D_WZ_H

#endif // WATERSHED3D_WZ_H

// Marker-based 3D watershed header

// By Wei Zheng
// Referred code (2D traditional watershed): https://ofstack.com/C++/22898/c++-watershed-algorithm-(watershed-algorithm).html

#define _USE_MATH_DEFINES

#include <vector>
#include <array>
#include <list>

#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display


using namespace std;
using namespace cv;

typedef unsigned char          GVByte;
typedef int32_t               GVInt32;
//typedef uint32_t               GVInt32U;

const GVByte              GV_BYTE_MAX = UCHAR_MAX;


array<GVInt32, 6> getNeighbor_WZ(const GVInt32 idx, const GVInt32 width, const GVInt32 height, const GVInt32 thick) {

    GVInt32 SLICE = width * height;
    GVInt32 z = idx / SLICE;
    GVInt32 y = idx / width;
    GVInt32 x = idx % width;

    array<GVInt32, 6> nIndex;
    nIndex[0] = (x == 0) ? -1 : (idx - 1);
    nIndex[1] = ((x + 1) == width) ? -1 : (idx + 1);
    nIndex[2] = (y == 0) ? -1 : (idx - width);
    nIndex[3] = ((y + 1) == height) ? -1 : (idx + width);
    nIndex[4] = (z == 0) ? -1 : (idx - SLICE);
    nIndex[5] = ((z + 1) == thick) ? -1 : (idx + SLICE);

    return nIndex;
}



void Watershed3D_WZ(
    const Mat im,
    const GVInt32 width,
    const GVInt32 height,
    const GVInt32 thick,
    Mat &L,
    const vector<vector<GVInt32>> marker)
{
    //<Parameter>
    //<image>  the image for watershed
    //<width>  the width of the image
    //<height> the height of the image
    //<thick>  the thick of the image
    //<label>  the map to save result. need to allocate memory before use watershed
    //<marker> the marker's index

    const GVByte* image=im.data;
    GVInt32* label=(GVInt32*)L.data;


    const GVInt32 SZ_slice = width * height;
    const GVInt32 SZ = SZ_slice * thick;
    const size_t markerNum = marker.size();

    // create labelMap. first is height, second is label
    vector<pair<GVByte, GVInt32>> labelBar;

    // create toSearchQueue. Saved pixel connected to labeled pixels and wait to search
    vector<list<GVInt32>> toSearchList;
    toSearchList.resize(markerNum);

    // set label to INIT (unsearched)
    ::memset(label, -1, sizeof(GVInt32) * SZ);

    // initialize
    array<GVInt32, 6> nIdx;
    labelBar.push_back(make_pair<GVByte, GVInt32>(0, 0)); // label for watershed (not used, just a place holder)
    for (size_t i = 0; i < markerNum; i++)
    {
        labelBar.push_back(make_pair<GVByte, GVInt32>(0, i + 1)); // labels for markers
        for (GVInt32 idx : marker[i])
        {
            // initialize label (which can be considered as a map of pointer to labelBar)
            label[idx] = i + 1;
            // initialize toSearchQueue
            nIdx = getNeighbor_WZ(idx, width, height, thick);
            for (GVInt32 newidx : nIdx)
            {
                if (newidx != -1)
                {
                    if (label[newidx] == -1) {
                        toSearchList[i].push_back(newidx);
                        label[newidx] = -2;
                    }
                }
            }
        }
    }

    //watershed
    GVByte h;
    for (int h_cnt = 0; h_cnt < (1+(int)GV_BYTE_MAX); h_cnt++) // water height
    {
        h = (GVByte)h_cnt;
        for (int cnt = 0; cnt < markerNum; cnt++) { // for each marker

            list<GVInt32>::iterator it = toSearchList[cnt].begin();
            while (!toSearchList[cnt].empty())
            {
                if (it == toSearchList[cnt].end()) { break; }
                // for each pixel connected to the cnt th labeled region
                GVInt32 idx = *it;
                // if this pixel is higher than water, ignore it
                if (image[idx] > h)
                {
                    it++;
                    continue;
                }
                // this pixel is lower than water, assign it
                label[idx] = cnt + 1;

                // add new neighbor
                nIdx = getNeighbor_WZ(idx, width, height, thick);
                for (GVInt32 newidx : nIdx)
                {
                    if (newidx != -1)
                    {
                        if (label[newidx] == -1) {
                            toSearchList[cnt].push_back(newidx);
                            label[newidx] = -2;
                        }
                    }
                }
                // erase searched pixel
                it = toSearchList[cnt].erase(it);
            }
        }
    }
}
