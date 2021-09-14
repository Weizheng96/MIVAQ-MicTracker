#ifndef PARAMETER4MICTRACKER_H
#define PARAMETER4MICTRACKER_H
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <string>
#include <vector>
#include <QVector3D>
#include <QRect>
#include <QMatrix4x4>


class Parameter4Mictracker
{
public:
    Parameter4Mictracker();

public:
    // the field of view of interest
    int yRg[2]={105,254};
    int xRg[2]={150,299};
    int zRg[2]={0,138};
    int tRg[2]={60,79};
//    int yRg[2]={0,5};
//    int xRg[2]={0,10};
//    int zRg[2]={0,5};
//    int tRg[2]={0,5};

    int nY=yRg[1]-yRg[0]+1;
    int nX=xRg[1]-xRg[0]+1;
    int nZ=zRg[1]-zRg[0]+1;
    int nT=tRg[1]-tRg[0]+1;

    // for segment
    double PCsigma=3;
    int minSz=50;
    int maxSz=10000;

};

#endif // PARAMETER4MICTRACKER_H
