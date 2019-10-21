#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//#include <DGtal/math/Statistic.h>

#include "Centerline.h" 
#include "BSplines.h"
#include "Regression.h"
#include "../Common/IOHelper.h"
#include "../Common/Rosin.h"
#include "../Common/Statistic.h"

using namespace DGtal;

/**
 * Depricated
 */
std::pair<Z3i::RealPoint, Z3i::RealPoint>
Centerline::computeLine(double &radius, RadiusEstimationType est){
    std::vector<double> radiis;
    std::vector<Z3i::RealPoint> mainLine;
    computeRawLine(mainLine);

    /*
    DGtal::Statistic<double> stat(true);
    for(int i = 0; i < radiis.size(); i++){
        stat.addValue(radiis[i]);
    }
    stat.terminate();
    
    switch(est){
        case mean   : radius = stat.mean(); break;
        case median : radius = stat.median(); break;
    }
    */

    radius = Statistic::getMode(radiis, 0.01);

    if(mainLine.size() == 0){
        Z3i::RealPoint p0(0,0,0);
        Z3i::RealPoint dir(1,0,0);
        std::pair<Z3i::RealPoint, Z3i::RealPoint> lineEq(p0, dir);
        return lineEq;
    }

    //Spline with 8 break
    //std::vector<Z3i::RealPoint> smoothLine;
    //std::sort(mainLine.begin(), mainLine.end(), pointOrderByZ);
        
    //3D parametric line regression
    //Xp = X0 + Vx*t
    //Yp = Y0 + Vy*t
    //Zp = Z0 + Vz*t
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> T;
    for(int i = 0; i < mainLine.size(); i++){
        X.push_back(mainLine[i][0]);
        Y.push_back(mainLine[i][1]);
        Z.push_back(mainLine[i][2]);
        T.push_back(i);
    }
    std::pair<double, double> xeq =  Regression::rmse(T,X);
    std::pair<double, double> yeq =  Regression::rmse(T,Y);
    std::pair<double, double> zeq =  Regression::rmse(T,Z);

    Z3i::RealPoint p0(xeq.second, yeq.second, zeq.second);
    Z3i::RealPoint dir(xeq.first, yeq.first, zeq.first);
    std::pair<Z3i::RealPoint, Z3i::RealPoint> lineEq(p0, dir);
    return lineEq;
}

std::pair<Z3i::RealPoint, Z3i::RealPoint>
Centerline::computeLine(){
    std::vector<Z3i::RealPoint> mainLine;
    computeRawLine(mainLine);

    if(mainLine.size() == 0){
        Z3i::RealPoint p0(0,0,0);
        Z3i::RealPoint dir(1,0,0);
        std::pair<Z3i::RealPoint, Z3i::RealPoint> lineEq(p0, dir);
        return lineEq;
    }

    //Spline with 8 break
    //std::vector<Z3i::RealPoint> smoothLine;
    //std::sort(mainLine.begin(), mainLine.end(), pointOrderByZ);
        
    //3D parametric line regression
    //Xp = X0 + Vx*t
    //Yp = Y0 + Vy*t
    //Zp = Z0 + Vz*t
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> T;
    for(int i = 0; i < mainLine.size(); i++){
        X.push_back(mainLine[i][0]);
        Y.push_back(mainLine[i][1]);
        Z.push_back(mainLine[i][2]);
        T.push_back(i);
    }
    std::pair<double, double> xeq =  Regression::rmse(T,X);
    std::pair<double, double> yeq =  Regression::rmse(T,Y);
    std::pair<double, double> zeq =  Regression::rmse(T,Z);

    Z3i::RealPoint p0(xeq.second, yeq.second, zeq.second);
    Z3i::RealPoint dir(xeq.first, yeq.first, zeq.first);
    std::pair<Z3i::RealPoint, Z3i::RealPoint> lineEq(p0, dir);
    return lineEq;
}
std::vector<Z3i::RealPoint> 
Centerline::compute(){

    //get main centerline
    std::vector<Z3i::RealPoint> mainLine;
    computeRawLine(mainLine);
//IOHelper::export2Text(mainLine, "raw.xyz");
//debug
//for(Z3i::RealPoint p:mainLine){
//trace.info()<<p<<std::endl;
//}
std::cout<<"found " << mainLine.size()<<" points "<<std::endl;
    std::vector<Z3i::RealPoint> smoothLine;
    if(mainLine.size() < nbControlPoint + 2){
std::cout<<"too few points for bsplines: " << mainLine.size()<<" points "<<std::endl;
std::cout<<"Try a smaller value of confidenceThreshold or Inverse normal might work !!! "<<std::endl;
        return smoothLine;
    }
    PointOrderByZ pointOrderByZ;
    //Spline with 8 break
    std::sort(mainLine.begin(), mainLine.end(), pointOrderByZ);

    smoothLine = BSplines::bsplines(mainLine, nbControlPoint + 2);
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

    while(myDomain.isInside(nextFront)){
        frontVect.push_back(nextFront);
        nextFront = nextFront + firstVect;
    }

    std::vector<Z3i::RealPoint> backVect;
    Z3i::RealPoint nextBack = lastPoint + lastVect;
    while(myDomain.isInside(nextBack)){
        backVect.push_back(nextBack);
        nextBack = nextBack + lastVect;
    }

    std::reverse(frontVect.begin(), frontVect.end());

    frontVect.insert(frontVect.end(), smoothLine.begin(), smoothLine.end());
    frontVect.insert(frontVect.end(), backVect.begin(), backVect.end());
    smoothLine = frontVect;


    return smoothLine;
}

/**
 * current version
 */
void  
Centerline::computeRawLine(std::vector<Z3i::RealPoint> &rawLine){
    computeAccumulation(true);
    computeConfidence();

    Image3DDouble imageConfidence = getConfidenceImage();
    Image3DDouble imageRadiusAcc = getRadiusImageAcc();

    std::vector<Z3i::RealPoint> rs;
    std::vector<double> radiis;
    Z3i::Domain aDomain = imageConfidence.domain();
    //DGtal::Statistic<double> stat(true);

    std::vector<double> confidences;
    for(auto &v: aDomain){
        if(imageConfidence(v) > 0){
            confidences.push_back(imageConfidence(v));
        }
    }
    /*
     * if we use a small value of binwidth, the histogram might not be perfect
     */
    //confidenceThreshold = Rosin::compute(confidences, 0.001);
    confidenceThreshold = Rosin::compute(confidences, 0.1);
    std::cout<<"conf Rosin: "<<confidenceThreshold<<std::endl;

    /* initialize random seed: */
    srand (time(NULL));

    // generate random number, hoping that it is not used
    std::string rIndex = std::to_string(rand());
    IOHelper::export2Text(confidences, "confRosin" + rIndex);
    //Z3i::RealPoint centroid(0,0,0);
    double minRa = std::numeric_limits<double>::max();
    double maxRa = -std::numeric_limits<double>::max();
    for(auto &v: aDomain){
        if (imageConfidence(v) > confidenceThreshold ){
            rs.push_back(v);
            //centroid +=v;
            radiis.push_back(imageRadiusAcc(v));
           // stat.addValue(imageRadiusAcc(v));
        }
    } 
    //stat.terminate();

    //centroid /= rs.size();

    std::cout<< "nbPoint: "<< radiis.size()<<std::endl;
    IOHelper::export2Text(radiis, "radiis");

    //double median = stat.median();

    double radiiMode = Statistic::getMode(radiis, 0.01);
    //get main centerline
std::cout<<"radii Mode:"<<radiiMode<<std::endl;
    std::vector<Z3i::RealPoint> goodRadiis;
    //double epsilon = 0.5;
    for(int i = 0; i < radiis.size(); i++){
        //if(radiis[i] > (1-epsilon)*radiiMode && radiis[i] < (1 + epsilon)*radiiMode){
        //if(true){
            goodRadiis.push_back(rs[i]);
        //}
    }
    assert(goodRadiis.size() > 0);
    IOHelper::export2Text(goodRadiis, "goodRadiis" + rIndex + ".xyz");
    IOHelper::export2Text(rs, "nofilter" + rIndex + ".xyz");
    // analyse by Z segment
    double zLength = 100;

    PointOrderByZ pointOrderByZ;
    std::sort(goodRadiis.begin(), goodRadiis.end(), pointOrderByZ);
    double minZ = goodRadiis.at(0)[2];
    double maxZ = goodRadiis.at(goodRadiis.size() - 1) [2];

    int nbSegment = (maxZ - minZ) / zLength + 1;
    
    //@TODO: finish filter by x, y coordinates
    std::vector<std::vector<int> > segments(nbSegment);
    std::vector<double> modeX(nbSegment);
    std::vector<double> modeY(nbSegment);
    std::vector<std::vector<double> > Xs(nbSegment);
    std::vector<std::vector<double> > Ys(nbSegment);
    for ( int i = 0; i < goodRadiis.size(); i++){
        Z3i::RealPoint zp = goodRadiis[i];
        int segId = (zp[2] - minZ) / zLength;
        Xs[segId].push_back(zp[0]);
        Ys[segId].push_back(zp[1]);
        segments[segId].push_back(i);
    }

    double radiiThreshold = radiiMode * 0.05;
    int nbPoint = (int)(maxZ - minZ +1);
    std::vector<Z3i::RealPoint> meansByZ(nbPoint, Z3i::RealPoint(0,0,0));
    std::vector<int> counts(nbPoint, 0);
//std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";
    std::vector<Z3i::RealPoint> nofilter;
    for ( int segId = 0; segId < nbSegment; segId++ ){
        if(segments.size() > 0){
            double mX = Statistic::getMode(Xs[segId], 0.01);
            double mY = Statistic::getMode(Ys[segId], 0.01);
            std::cout<<mX<<std::endl;
            std::cout<<mY<<std::endl;
            for(int pId : segments[segId]){
                if(std::abs(goodRadiis[pId][0] - mX) < radiiThreshold && std::abs(goodRadiis[pId][1] - mY) < radiiThreshold){
                    //rawLine.push_back(goodRadiis[pId]);

                //if(true){
                    meansByZ[(int)(goodRadiis[pId][2] - minZ)] += goodRadiis[pId];
                    counts[(int)(goodRadiis[pId][2] - minZ)]++;
                }
            }
        }
    }

    for (int i = 0; i < nbPoint; i++){
        if(counts[i] > 0){
            rawLine.push_back(meansByZ[i] / counts[i]);
        }
    }
    // clustering
}

