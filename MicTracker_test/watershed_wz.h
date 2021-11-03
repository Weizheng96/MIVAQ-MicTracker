#ifndef WATERSHED_WZ_H
#define WATERSHED_WZ_H

#define _USE_MATH_DEFINES

#include <vector>
#include <array>
#include <list>
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
#include <QDebug>
#include <QTextStream>

#include <chrono>

using namespace std;
using namespace cv;

typedef unsigned char          GVByte;
typedef int32_t               GVInt32;
//typedef uint32_t               GVInt32U;

const GVByte              GV_BYTE_MAX = UCHAR_MAX;


class WaterShed_WZ
{
public:
    WaterShed_WZ();
    static array<GVInt32, 6> getNeighbor_WZ(const GVInt32 idx, const GVInt32 width, const GVInt32 height, const GVInt32 thick);
    static void Watershed3D_WZ(Mat im,
            GVInt32 width,
            GVInt32 height,
            GVInt32 thick,
            Mat *L,
            vector<vector<GVInt32> > marker);
};

#endif // WATERSHED_WZ_H
