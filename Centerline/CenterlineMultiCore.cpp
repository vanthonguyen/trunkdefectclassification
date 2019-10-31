#include <thread>
#include <queue>          // std::priority_queue
#include <functional>     // std::greater

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <DGtal/math/Statistic.h>

#include "CenterlineMultiCore.h" 
#include "BSplines.h"
#include "Regression.h"
#include "../Common/IOHelper.h"
#include "../Common/Rosin.h"
#include "../Common/IOHelper.h"

using namespace DGtal;


std::vector<Z3i::RealPoint> 
CenterlineMultiCore::compute(){
    Z3i::RealPoint upPt = points[0];
    Z3i::RealPoint lowPt = points[0];
    for(Z3i::RealPoint p: points){
        lowPt = lowPt.inf(p);
        upPt = upPt.sup(p);
    }

    //get main centerline
    std::vector<Z3i::RealPoint> mainLine = getRawLine(lowPt, upPt);

//IOHelper::export2Text(mainLine, "raw.xyz");
std::cout<<"found " << mainLine.size()<<" points "<<std::endl;
    std::vector<Z3i::RealPoint> smoothLine;
    if(mainLine.size() < nbControlPoint + 2){
        std::cout<<"too few points for bsplines: " << mainLine.size()<<" points "<<std::endl;
        std::cout<<"CenterlineMultiCore::compute() -> Try a smaller value of confidenceThreshold or Inverse normal might work !!! "<<std::endl;
        return smoothLine;
    }

    PointOrderByZ pointOrderByZ;
    //Spline with 8 break
    std::sort(mainLine.begin(), mainLine.end(), pointOrderByZ);

    int nbKnots = nbControlPoint + 2;
/*
    double minZ = mainLine[0][2];
    double maxZ = mainLine[mainLine.size() - 1][2];
    //knots
    int segmentLengthForKnots = (maxZ - minZ)  / nbKnots;
    
    //O^2 algo
    std::vector<int> breakPointInds(nbKnots);
    std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > > dists(nbKnots);

    for(int i = 0; i < mainLine.size(); i++){
        int sId = (mainLine[i][2] - minZ)/segmentLengthForKnots;
        if(sId > nbKnots - 1){
            sId = nbKnots - 1;
        }
        for (int j = 0; j < mainLine.size(); j++){
            if(j != i){
                double d = mainLine[i] - mainLine 
                dists[sId].push();
            }
        }
    }
    */

//IOHelper::export2Text(mainLine, "rawsort.xyz");
    smoothLine = BSplines::bsplines(mainLine, nbKnots);
    /*
    //test different nb of control points
    trace.info()<< "test different nb of control points"<<std::endl;
    for(int i = 2; i < 10; i++){
        trace.info()<<"nb control :" << i <<std::endl;
        BSplines::bsplines(mainLine, i + 2);
    }
    */
    //smoothLine = BSplines::splines(pointsForSpline, 0.05);
    
    Z3i::RealPoint firstVect = smoothLine[0] - smoothLine[1];
    Z3i::RealPoint firstPoint = smoothLine[0];

    Z3i::RealPoint lastVect = smoothLine[smoothLine.size()-1] - smoothLine[smoothLine.size()-2];
    Z3i::RealPoint lastPoint = smoothLine[smoothLine.size()-1];

    std::vector<Z3i::RealPoint> frontVect;

    Z3i::RealPoint nextFront = firstPoint + firstVect;

    //extend to the 2 ends of trunk
    //mm -> m
    double zLow = lowPt[2] - 0.1;
    double zUp = upPt[2] + 0.1;

    while(nextFront[2] > zLow && nextFront[2] < zUp){
        frontVect.push_back(nextFront);
        nextFront = nextFront + firstVect;
    }

    std::vector<Z3i::RealPoint> backVect;
    Z3i::RealPoint nextBack = lastPoint + lastVect;
    while(nextBack[2] > zLow && nextBack[2] < zUp){
        backVect.push_back(nextBack);
        nextBack = nextBack + lastVect;
    }

    std::reverse(frontVect.begin(), frontVect.end());

    frontVect.insert(frontVect.end(), smoothLine.begin(), smoothLine.end());
    frontVect.insert(frontVect.end(), backVect.begin(), backVect.end());
    smoothLine = frontVect;


    return smoothLine;

}

std::vector<Z3i::RealPoint> 
CenterlineMultiCore::getRawLine(const Z3i::RealPoint &lowPt, const Z3i::RealPoint &upPt){

    //if voxel size is 3 -> 900m
    //double segLengthInitial = 300.0;
    //separe into multiple parts
    //zLength = 
    
    double zLength = upPt[2] -lowPt[2];
    //int nbSegment = zLength / segLengthInitial;

    double zMin = lowPt[2];

    //adjusted
    double segLength = zLength / nbSegment;
    double segLengthWithBuff = segLength*1.1;
    double segOverlapLength = segLength*0.1;

    std::vector<double> beginSegments(nbSegment);
    std::vector<double> endSegments(nbSegment);

    for(int i = 0; i < nbSegment; i++){
        beginSegments[i] = zMin + segLength * i;
        endSegments[i] = zMin + segLength * (i + 1) + segOverlapLength;
        //endSegments[i] = zMin + segLengthWithBuff * (i + 1);
    }

    std::vector< std::vector<Z3i::RealPoint> > segPoints(nbSegment);
    //for presentation
    std::vector< std::vector<int> > segPids(nbSegment);
    std::vector< std::vector<Z3i::RealPoint> > segNormals(nbSegment);
    Z3i::RealPoint zero(0,0,0);
    float dMax = std::numeric_limits<float>::max();
    std::vector<Z3i::RealPoint> segLowPts(nbSegment, Z3i::RealPoint(dMax, dMax, dMax));
    std::vector<Z3i::RealPoint> segUpPts(nbSegment, Z3i::RealPoint(-dMax, -dMax, -dMax));
    std::vector<Z3i::Domain> segDomains(nbSegment);

    for(int pId = 0; pId < points.size(); pId++){
        Z3i::RealPoint p = points[pId];
        Z3i::RealPoint n = normals[pId];
        for(int segId = 0; segId < nbSegment; segId++){
            if( p[2] > beginSegments[segId] && p[2] < endSegments[segId]) {
                //for presentation
                segPids[segId].push_back(pId);
                segPoints[segId].push_back(p);
                segNormals[segId].push_back(n);
                segLowPts[segId] = segLowPts[segId].inf(p);
                segUpPts[segId] = segUpPts[segId].sup(p);
            }
        }
    }

    std::vector<Z3i::RealPoint>  rawLine;
    /**
     * one core
    for(int segId = 0; segId < nbCores; segId++){
        Z3i::Domain sDomain = Z3i::Domain(segLowPts[segId] - DGtal::Z3i::RealPoint::diagonal(1), segUpPts[segId] + DGtal::Z3i::RealPoint::diagonal(1));
        trace.info()<<accRadius << "  "<< sDomain << " "<< segNormals[segId].size()<<"  "<<segPoints[segId].size()<< "  "<< nbControlPoint<<std::endl;
        if(segPoints[segId].size() == 0){
            continue;
        }
        //threshold is deprecated
        Centerline cen(accRadius, sDomain, segNormals[segId], segPoints[segId], 0.5, nbControlPoint);
        //futures[segId] = std::async(std::launch::async, [&cen](){ return cen.getRawLine();});
        
        std::vector<Z3i::RealPoint> line = cen.getRawLine();
        rawLine.insert(rawLine.end(), line.begin(), line.end());
    }
    **/


    std::vector<std::thread> ts(nbSegment);
    //std::vector<std::future<std::vector<Z3i::RealPoint> > > futures(nbSegment);
    std::vector<std::vector<Z3i::RealPoint> > raws(nbSegment);
    for(int segId = 0; segId < nbSegment; segId++){
        if(segPoints[segId].size() == 0){
            continue;
        }
        Z3i::Domain sDomain = Z3i::Domain(segLowPts[segId] - DGtal::Z3i::RealPoint::diagonal(1), segUpPts[segId] + DGtal::Z3i::RealPoint::diagonal(1));
        trace.info()<<accRadius << "  "<< sDomain << " "<< segNormals[segId].size()<<"  "<<segPoints[segId].size()<< "  "<< nbControlPoint<<std::endl;
        //threshold is deprecated
        Centerline cen(accRadius, sDomain, segNormals[segId], segPoints[segId], threshold, nbControlPoint);
        //ts[segId] = std::thread(&Centerline::computeRawLine, &cen, std::ref(raws[segId]));
        //std::vector<Z3i::RealPoint> raw;
        //ts[segId] = std::thread(&CenterlineMultiCore::computeRaw, this, accRadius, sDomain, segNormals[segId], segPoints[segId], nbControlPoint, raws[segId]);
        //ts[segId] = std::thread([&cen] { computeRawLine(&raws[segId]); });
        cen.computeRawLine(raws[segId]);

        //for presentation
        std::string segPointName = "segPoints" + std::to_string(segId) + ".xyz";
        IOHelper::export2Text(segPoints[segId], segPointName);

        std::string segPiDName = "segPid" + std::to_string(segId) + ".id";
        IOHelper::export2Text(segPids[segId], segPiDName);

        std::string rawLineName = "rawline" + std::to_string(segId) + ".xyz";
        IOHelper::export2Text(raws[segId], rawLineName);
    }
    //for(int i = 0; i < nbSegment; i++){
    //    ts[i].join();
    //}

    
    for(int i = 0; i < nbSegment; i++){
        rawLine.insert(rawLine.end(), raws[i].begin(), raws[i].end());
    }
    return rawLine; 
}


