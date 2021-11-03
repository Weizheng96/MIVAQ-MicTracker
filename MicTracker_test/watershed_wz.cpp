#include "watershed_wz.h"

WaterShed_WZ::WaterShed_WZ()
{

}


array<GVInt32, 6> WaterShed_WZ::getNeighbor_WZ(const GVInt32 idx, const GVInt32 width, const GVInt32 height, const GVInt32 thick) {

    GVInt32 SLICE = width * height;
    GVInt32 z = idx / SLICE;
    GVInt32 y = (idx%SLICE) / width;
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



//void WaterShed_WZ::Watershed3D_WZ(
//    const Mat im,
//    const GVInt32 width,
//    const GVInt32 height,
//    const GVInt32 thick,
//    GVInt32* label,
//    const vector<vector<GVInt32>> marker)
//{
//    //<Parameter>
//    //<image>  the image for watershed
//    //<width>  the width of the image
//    //<height> the height of the image
//    //<thick>  the thick of the image
//    //<label>  the map to save result. need to allocate memory before use watershed
//    //<marker> the marker's index

////    const GVByte* image=im.data;

//    auto t0 = chrono::high_resolution_clock::now();
//    QTextStream out(stdout);


////    const GVInt32 SZ_slice = width * height;
////    const GVInt32 SZ = SZ_slice * thick;
//    const GVInt32 markerNum = marker.size();

//    // create toSearchList. Saved pixel connected to labeled pixels and wait to search
//    vector<list<GVInt32>> toSearchList;
//    toSearchList.resize(markerNum);

//    // set label to INIT (unsearched)
////    ::memset(label, -1, sizeof(GVInt32) * SZ);

//    // initialize
//    array<GVInt32, 6> nIdx;
//    for (size_t i = 0; i < markerNum; i++)
//    {
//        for (GVInt32 idx : marker[i])
//        {
//            // initialize label (which can be considered as a map of pointer to labelBar)
//            label[idx] = i + 1;
//            nIdx = getNeighbor_WZ(idx, width, height, thick);
//            for (GVInt32 newidx : nIdx)
//            {
//                if (newidx != -1)
//                {
//                    if (label[newidx] == -1) {

//                        toSearchList[i].push_back(newidx);

//                        label[newidx] = -2;
//                    }
//                }
//            }
//        }
//    }

//    //watershed
//    GVByte h;
//    GVInt32 idx;
//    for (int h_cnt = 0; h_cnt < (1+(int)GV_BYTE_MAX); h_cnt++) // water height
//    {
//        h = (GVByte)h_cnt;
//        for (GVInt32 cnt = 0; cnt < markerNum; cnt++) { // for each marker

//            list<GVInt32>::iterator it = toSearchList[cnt].begin();
//            while (!toSearchList[cnt].empty())
//            {
//                // for each pixel connected to the cnt-th labeled region

//                idx = *it;
//                // if this pixel is higher than water, ignore it
//                if (im.at<unsigned char>(idx) > h)
//                {
//                    ++it;
//                    if(it == toSearchList[cnt].end())
//                    {
//                        break;
//                    }
//                    else
//                    {
//                        continue;
//                    }
//                }
//                // this pixel is lower than water, assign it
//                label[idx] = cnt + 1;
////                L.at<int>(idx)=cnt + 1;

//                // add new neighbor
//                nIdx = getNeighbor_WZ(idx, width, height, thick);
//                for (GVInt32 newidx : nIdx)
//                {
//                    if (newidx != -1)
//                    {
//                        if (label[newidx]== -1) {
//                            toSearchList[cnt].push_back(newidx);
//                            label[newidx] = -2;
////                            L.at<int>(newidx)=-2;
//                        }
//                    }
//                }
//                // erase searched pixel
//                it = toSearchList[cnt].erase(it);

//                if(it == toSearchList[cnt].end())
//                {
//                    break;
//                }
//                else
//                {
//                    continue;
//                }
//            }
//        }
//    }

//    auto t1 = chrono::high_resolution_clock::now();
//    auto dt = 1.e-9 * chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
//    out << "Watershed used " << dt << " seconds.\n\n" << Qt::endl;
//}

// slow version

//void WaterShed_WZ::Watershed3D_WZ(
//    Mat im,
//    GVInt32 width,
//    GVInt32 height,
//    GVInt32 thick,
//    Mat *L,//GVInt32* label,
//    vector<vector<GVInt32>> marker)
//{
//    //<Parameter>
//    //<image>  the image for watershed
//    //<width>  the width of the image
//    //<height> the height of the image
//    //<thick>  the thick of the image
//    //<label>  the map to save result. need to allocate memory before use watershed
//    //<marker> the marker's index

////    const GVByte* image=im.data;

//    chrono::steady_clock::time_point begin = chrono::steady_clock::now();


////    const GVInt32 SZ_slice = width * height;
////    const GVInt32 SZ = SZ_slice * thick;
//    const GVInt32 markerNum = marker.size();

//    // create toSearchList. Saved pixel connected to labeled pixels and wait to search
//    vector<vector<GVInt32>> toSearchList;
//    toSearchList.resize(markerNum);

//    // set label to INIT (unsearched)
////    ::memset(label, -1, sizeof(GVInt32) * SZ);

//    // initialize
//    array<GVInt32, 6> nIdx;
//    for (size_t i = 0; i < markerNum; i++)
//    {
//        for (GVInt32 idx : marker[i])
//        {
//            // initialize label (which can be considered as a map of pointer to labelBar)
////            label[idx] = i + 1;
//            L->at<int>(idx)=i+1;
//            nIdx = getNeighbor_WZ(idx, width, height, thick);
//            for (GVInt32 newidx : nIdx)
//            {
//                if (newidx != -1)
//                {
//                    if (L->at<int>(newidx)==-1) {

//                        toSearchList[i].push_back(newidx);

////                        label[newidx] = -2;
//                        L->at<int>(newidx)=-2;
//                    }
//                }
//            }
//        }
//    }

//    //watershed
//    GVByte h;
//    GVInt32 idx;
//    for (int h_cnt = 0; h_cnt < (1+(int)GV_BYTE_MAX); h_cnt++) // water height
//    {
//        h = (GVByte)h_cnt;
//        for (GVInt32 cnt = 0; cnt < markerNum; cnt++) { // for each marker
//            GVInt32 cnt_idx=0;
//            while (cnt_idx<toSearchList[cnt].size())
//            {
//                // for each pixel connected to the cnt-th labeled region

//                idx = toSearchList[cnt][cnt_idx];
//                // if this pixel is higher than water, ignore it
//                if (im.at<unsigned char>(idx) > h)
//                {
//                    ++cnt_idx;
//                    if(cnt_idx==toSearchList[cnt].size())
//                    {
//                        break;
//                    }
//                    else
//                    {
//                        continue;
//                    }
//                }
//                // this pixel is lower than water, assign it
////                label[idx] = cnt + 1;
//                L->at<int>(idx)=cnt+1;

//                // add new neighbor
//                nIdx = getNeighbor_WZ(idx, width, height, thick);
//                for (GVInt32 newidx : nIdx)
//                {
//                    if (newidx != -1)
//                    {
//                        if (L->at<int>(newidx)== -1) {
//                            toSearchList[cnt].push_back(newidx);
//                            L->at<int>(newidx) = -2;
//                        }
//                    }
//                }
//                // erase searched pixel
//                toSearchList[cnt].erase(toSearchList[cnt].begin() + cnt_idx);

//                if(cnt_idx==toSearchList[cnt].size())
//                {
//                    break;
//                }
//                else
//                {
//                    continue;
//                }
//            }
//        }
//    }

//    chrono::steady_clock::time_point end = chrono::steady_clock::now();
//    qInfo("Watershed used used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);
//}

void WaterShed_WZ::Watershed3D_WZ(
    Mat im,
    GVInt32 width,
    GVInt32 height,
    GVInt32 thick,
    Mat *L,//GVInt32* label,
    vector<vector<GVInt32>> marker)
{
    //<Parameter>
    //<image>  the image for watershed
    //<width>  the width of the image
    //<height> the height of the image
    //<thick>  the thick of the image
    //<label>  the map to save result. need to allocate memory before use watershed
    //<marker> the marker's index

//    const GVByte* image=im.data;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();


//    const GVInt32 SZ_slice = width * height;
//    const GVInt32 SZ = SZ_slice * thick;
    const GVInt32 markerNum = marker.size();

    // create toSearchList. Saved pixel connected to labeled pixels and wait to search
    vector<list<GVInt32>> toSearchList;
    toSearchList.resize(markerNum);

    // set label to INIT (unsearched)
//    ::memset(label, -1, sizeof(GVInt32) * SZ);

    // initialize
    array<GVInt32, 6> nIdx;
    for (size_t i = 0; i < markerNum; i++)
    {
        for (GVInt32 idx : marker[i])
        {
            // initialize label (which can be considered as a map of pointer to labelBar)
//            label[idx] = i + 1;
            L->at<int>(idx)=i+1;
            nIdx = getNeighbor_WZ(idx, width, height, thick);
            for (GVInt32 newidx : nIdx)
            {
                if (newidx != -1)
                {
                    if (L->at<int>(newidx)==-1) {

                        toSearchList[i].push_back(newidx);

//                        label[newidx] = -2;
                        L->at<int>(newidx)=-2;
                    }
                }
            }
        }
    }

    //watershed
    GVByte h;
    GVInt32 idx;
    for (int h_cnt = 0; h_cnt < (1+(int)GV_BYTE_MAX); h_cnt++) // water height
    {
        h = (GVByte)h_cnt;
        for (GVInt32 cnt = 0; cnt < markerNum; cnt++) { // for each marker
            list<GVInt32>::iterator it = toSearchList[cnt].begin();
            while (!toSearchList[cnt].empty())
            {
                // for each pixel connected to the cnt-th labeled region

                idx = *it;
                // if this pixel is higher than water, ignore it
                if (im.at<unsigned char>(idx) > h)
                {
                    ++it;
                    if(it == toSearchList[cnt].end())
                    {
                        break;
                    }
                    else
                    {
                        continue;
                    }
                }
                // this pixel is lower than water, assign it
//                label[idx] = cnt + 1;
                L->at<int>(idx)=cnt+1;

                // add new neighbor
                nIdx = getNeighbor_WZ(idx, width, height, thick);
                for (GVInt32 newidx : nIdx)
                {
                    if (newidx != -1)
                    {
                        if (L->at<int>(newidx)== -1) {
                            toSearchList[cnt].push_back(newidx);
                            L->at<int>(newidx) = -2;
                        }
                    }
                }
                // erase searched pixel
                it = toSearchList[cnt].erase(it);

                if(it == toSearchList[cnt].end())
                {
                    break;
                }
                else
                {
                    continue;
                }
            }
        }
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    qInfo("Watershed used used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);
}
