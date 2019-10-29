#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <thread>
//debug
#include <stdlib.h>
#include <time.h>


#include "../Common/Statistic.h"
#include "../Common/IOHelper.h"
#include "../Common/DefectType.h"
#include "../Common/MultiThreadHelper.h"
#include "../Feature/PredictRF.h"
#include "../Feature/FeatureExtractionXyz.h"


#include "DefectQuantification.h"


using namespace DGtal;


Defect
DefectQuantification::quantify(){
    Defect defectInfo;

    CylindricalCoordinateSystem ccs(centerline, seedPoint);

    Z3i::RealPoint defectLowPoint;
    Z3i::RealPoint defectHighPoint;

    //use points with depth < 20mm to compute center of defect
    double maxDepth = 10;

    double minR = std::numeric_limits<double>::max();

    bool first = true;
    std::vector<CylindricalPoint> cylPoints(pointCloud.size());

    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint xyzPoint = pointCloud[i];
        cylPoints[i] = ccs.xyz2Cylindrical(xyzPoint);

//std::cout<<"#####"<< cylPoints[i].radius << "#"<< cylPoints[i].angle<<"#"<<cylPoints[i].height<<std::endl;
        if(cylPoints[i].radius < minR){
            minR = cylPoints[i].radius;
        }
    }
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        CylindricalPoint cylPoint = cylPoints[i];

        if(cylPoint.radius - minR < maxDepth){
            Z3i::RealPoint xyzPoint = pointCloud[i];
            if(first){
                defectHighPoint = xyzPoint;
                defectLowPoint = xyzPoint;
                first = false;
            }else{
                defectLowPoint = defectLowPoint.inf(xyzPoint);
                defectHighPoint = defectHighPoint.sup(xyzPoint);
            }
        }
    }

    /*
    double width = feature.at(7);
    double height = width * feature.at(8);
    double depth = width * feature.at(9);
    Z3i::RealPoint dim(width, height, depth);
    */
    
    double width = sqrt( (defectHighPoint[0] - defectLowPoint[0])*(defectHighPoint[0] - defectLowPoint[0]) +  (defectHighPoint[1] - defectLowPoint[1])*(defectHighPoint[1] - defectLowPoint[1]) ) ;
    //double height = defectHighPoint[2] - defectLowPoint[2];

    Z3i::RealPoint defectCenter = (defectLowPoint + defectHighPoint)/2;
    CylindricalPoint cylDefectCenter = ccs.xyz2Cylindrical(defectCenter);
    double trunkRadiusAtDefect = segmentRadius[cylDefectCenter.segmentId];
    defectInfo.coordinateZ = defectCenter[2] - seedPoint[2] + padding;

    double defA = cylDefectCenter.angle;
    //arc length from defect to repere point
    if (defA > M_PI){
        defA = -2*M_PI + defA;
    }

    defectInfo.coordinateL = defA * trunkRadiusAtDefect;

    //CylindricalPoint cylDefLowPoint = ccs.xyz2Cylindrical(defectLowPoint);
    //CylindricalPoint cylDefHighPoint = ccs.xyz2Cylindrical(defectHighPoint);


    defectInfo.height = defectHighPoint[2] - defectLowPoint[2];
    //defectInfo.width = std::abs(cylDefHighPoint.angle - cylDefLowPoint.angle)*trunkRadiusAtDefect;
    defectInfo.width = width;
    

    if(pointCloud.size() < 4){
        defectInfo.type = BARK;
        return defectInfo;
    }

    PredictRF prf(rfFilePath);

    FeatureExtractionXyz ef(pointCloud, centerline, trunkRadiusAtDefect);
    std::vector<float> feature = ef.compute();

    double alpha = feature.at(8) / trunkRadiusAtDefect;
    if(alpha > 0.2){
        width = feature.at(8);
        defectInfo.width = width;
    }

    feature.insert(feature.begin(), species);

    std::vector<float> featureForRf = feature;
    std::cout<<"####Feature:"<<std::endl;
    for(unsigned int i = 0; i < feature.size() - 1; i++){
        std::cout<< feature[i]<<" ";
    }

    std::cout<< feature[feature.size() - 1]<<std::endl;
    //std::vector<float> featureForRf;

    //I want to test without moments
    //I want to test without axis
//    for(int i = 0; i < feature.size(); i++){
        /*
        switch(i){
            case 0:
            case 1:
            case 9 ... 15:
                featureForRf.push_back(feature[i]);
                break;
        }
        */
        
//        switch(i){
//            case 0 ... 12:
//                featureForRf.push_back(feature[i]);
//                break;
//        }
        //featureForSvm.push_back(feature[i]);
        //std::cout<< feature[i]<<" ";
//    }

    //debug
//    for(int i = 0; i < feature.size() - 1; i++){
//        std::cout<< feature[i]<<" ";
//    }
//    std::cout<< feature[feature.size() - 1]<<std::endl;
    //end debug
 

    int label = prf.predict(featureForRf);

    defectInfo.type = label;

    defectInfo.branchDiameter = 0;
    if(label == BRANCH || label == BURL ){
        defectInfo.branchDiameter = getBranchDiameter();
    }
    //return Defect();
    return defectInfo;
}


double
DefectQuantification::getBranchDiameter(){

    std::vector<Z3i::RealPoint> branches;
    Z3i::RealPoint seedPoint(0,1,0);
    CylindricalCoordinateSystem ccs(centerline, seedPoint);
    double minR = std::numeric_limits<double>::max();

    //@TODO: using parameters:
    double skip = 30;
    double maxLength = 200;
    double minNbPoints = 500;

    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint xyzPoint = pointCloud.at(i);
        CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);
        if(cylPoint.radius < minR){
            minR = cylPoint.radius;
        }
    }
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint xyzPoint = pointCloud.at(i);
        CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);
        if(cylPoint.radius > minR + skip && cylPoint.radius < minR + skip + maxLength){
            branches.push_back(xyzPoint);
        }
    }

    if(branches.size() < minNbPoints){
        return 0;
    }

std::cout<<"nb points in branch"<<branches.size()<<std::endl;
    //@TODO: using parameters:
    double accRadius = 100;
    double searchRadius = 20;

    std::vector<Z3i::RealPoint> branhceCenterline = CenterlineHelper::computeCenterline(branches, 1, searchRadius, accRadius, 2, 5, 0);
    
    if( branhceCenterline.size() < 2 ){
        return 0;
    }

    CylindricalCoordinateSystem branchCcs(branhceCenterline, seedPoint);
    std::vector<double> radiis;

    double radiiMin = std::numeric_limits<double>::max();
    double radiiMax = -std::numeric_limits<double>::max();

    std::vector<double> segmentRadiiMin(branhceCenterline.size() - 1, std::numeric_limits<double>::max());
    std::vector<double> segmentRadiiMax(branhceCenterline.size() - 1, -std::numeric_limits<double>::max());
    std::vector<std::vector<double> > radiiBySegment(branhceCenterline.size() - 1);
    for(unsigned int i = 0; i < branches.size(); i++){
        Z3i::RealPoint xyzPoint = branches.at(i);
        CylindricalPoint cylPoint = branchCcs.xyz2Cylindrical(xyzPoint);
        CylindricalPoint pointInTrunkCcs = ccs.xyz2Cylindrical(xyzPoint);
        if(pointInTrunkCcs.radius < minR + skip + maxLength){
            radiis.push_back(cylPoint.radius);
            if(cylPoint.radius < radiiMin){
                radiiMin = cylPoint.radius;
            }
            if(cylPoint.radius > radiiMax){
                radiiMax = cylPoint.radius;
            }
        }
        radiiBySegment[cylPoint.segmentId].push_back(cylPoint.radius);
        int segmentId = cylPoint.segmentId;
        //std::cout<<segmentId;
        if(cylPoint.radius < segmentRadiiMin[segmentId]){
            segmentRadiiMin[segmentId] = cylPoint.radius;
        }

        if(cylPoint.radius > segmentRadiiMax[segmentId]){
            segmentRadiiMax[segmentId] = cylPoint.radius;
        }
    }
    
std::cout<<"radiis nb:"<<radiis.size()<<std::endl;

    double branchRadiusEst = 2*Statistic::getMode(radiis, radiiMin, radiiMax, 1);

    return branchRadiusEst;
}
