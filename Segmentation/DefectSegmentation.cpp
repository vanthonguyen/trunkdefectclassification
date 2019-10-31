#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <thread>
//debug
#include <stdlib.h>
#include <time.h>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Householder>

#include "../Common/Regression.h"
#include "../Common/Statistic.h"
#include "../Common/IOHelper.h"
#include "../Common/MultiThreadHelper.h"

#include "DefectSegmentation.h"

using namespace DGtal;


const std::pair<double, double> DefectSegmentation::coeffZero = std::make_pair(0.0, 0.0);

void
DefectSegmentation::init(){
    lightInit();
    computeAngleLimit();

    /**fore debug*/
    isForExport.resize(pointCloud.size());
    /**
    int maxRand = pointCloud.size() / 100;

    srand (time(NULL));
    isForExport.resize(pointCloud.size());
    for (int i = 0; i < pointCloud.size(); i++){
        int aNum = rand() % maxRand;
std::cout<<"aNum:"<<aNum<<std::endl;
        if(aNum == 0){
            isForExport[i] = true;
            std::cout<<"to export:"<<i<<std::endl;
        }else{
            isForExport[i] = false;
        }
    }

    computeEquationsFullPcl();

    **/
    /*end debug*/
    computeEquations();
    //computeEquationsKdTree2D();

    computeDistances();
    //debug
}

void
DefectSegmentation::lightInit(){
    allocate();
    allocateExtra();

    computeBeginOfSegment();
    computeVectorMarks();
    computePlaneNormals();

	convertToCcs();
    //debug
}

void 
DefectSegmentation::allocateExtra(){
    //allocate
    coefficients.resize(pointCloud.size(), coeffZero);
}


void
DefectSegmentation::computeEquationsKdTree2D(){

}

void
DefectSegmentation::computeEquationsKdTree2DMultiThreads(int threadId, int nbThread,
            const pcl::KdTreeFLANN<pcl::PointXY> &qtree1,
            const pcl::KdTreeFLANN<pcl::PointXY> &qtree2,
            const std::vector<CylindricalPoint> &points1,
            const std::vector<CylindricalPoint> &points2,
            double minHeight, double maxHeight){}

void 
DefectSegmentation::computeEquations(){

    //w = patch width, wh = patch height
    //build kdtree using pcl 
    trace.info()<<"Begin Compute Equation"<<std::endl;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
    //get subsample cloud
    //
    std::vector<Z3i::RealPoint> subsampleCloud;

    std::vector<unsigned int> subsampleCloudIds;

    //estimate radius
    
    //@FIX? == voxelsize
    double resolution = 3;
    //hardcode, resolution = 5, radius = 180 is not a good idea!!
    SegmentationHelper::subSample(pointCloud, resolution, radii, subsampleCloudIds);

    IOHelper::export2Text(subsampleCloud, "subsample.xyz");

    Z3i::RealPoint upPt = pointCloud[subsampleCloudIds[0]];
    Z3i::RealPoint lowPt = pointCloud[subsampleCloudIds[0]];
    for(int i = 0; i < subsampleCloudIds.size(); i++){
        Z3i::RealPoint p = pointCloud.at(subsampleCloudIds.at(i));
        cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
        upPt = upPt.sup(p);
        lowPt = lowPt.inf(p);
    }

    upPt += Z3i::RealPoint(resolution,resolution,resolution);
    lowPt -= Z3i::RealPoint(resolution,resolution,resolution);

    std::map<unsigned int, std::vector<unsigned int> > indexMap;
    SegmentationHelper::voxelize(pointCloud, resolution, lowPt, upPt, indexMap);

    //convert map to vector
    std::vector<std::vector<unsigned int> > voxelPointMaps;
    for(auto it = indexMap.begin(); it != indexMap.end(); ++it ) {
        voxelPointMaps.push_back( it->second );
    }


    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud (cloudPcl);

    CylindricalPointOrder heightOrder;
    
    auto minMaxElem = std::minmax_element(myPoints.begin(), myPoints.end(), heightOrder);
    double minHeight = (*minMaxElem.first).height;
    double maxHeight = (*minMaxElem.second).height;

    int nbCores = getNumCores();
    std::vector<std::thread> ts;
    for(int i = 0; i < nbCores - 1; i++){
        ts.push_back(std::thread(&DefectSegmentation::computeEquationsMultiThread, this, i, nbCores, kdtree, minHeight, maxHeight, subsampleCloudIds, voxelPointMaps));
    }
    computeEquationsMultiThread(nbCores - 1, nbCores, kdtree, minHeight, maxHeight, subsampleCloudIds, voxelPointMaps);
    for(int i = 0; i < nbCores - 1; i++){
        ts[i].join();
    }
    trace.info()<<"finish eq"<<std::endl;
    IOHelper::export2Text(coefficients, "coeffs");
}

void 
DefectSegmentation::computeEquationsFullPcl(){
    //w = patch width, wh = patch height
    //build kdtree using pcl 
    trace.info()<<"Begin Compute Equation Full"<<std::endl;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
    //get subsample cloud

    std::vector<unsigned int> subsampleCloudIds;

    double resolution = 5;
    Z3i::RealPoint upPt = pointCloud[0];
    Z3i::RealPoint lowPt = pointCloud[0];
    for(int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint p = pointCloud.at(i);
        cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
        subsampleCloudIds.push_back(i);
        upPt = upPt.sup(p);
        lowPt = lowPt.inf(p);
    }

    //upPt += Z3i::RealPoint(resolution,resolution,resolution);
    //lowPt -= Z3i::RealPoint(resolution,resolution,resolution);

    std::map<unsigned int, std::vector<unsigned int> > indexMap;
    SegmentationHelper::voxelize(pointCloud, resolution, lowPt, upPt, indexMap);

    //convert map to vector
    std::vector<std::vector<unsigned int> > voxelPointMaps;
    for(auto it = indexMap.begin(); it != indexMap.end(); ++it ) {
        voxelPointMaps.push_back( it->second );
    }

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud (cloudPcl);

    CylindricalPointOrder heightOrder;
    
    auto minMaxElem = std::minmax_element(myPoints.begin(), myPoints.end(), heightOrder);
    double minHeight = (*minMaxElem.first).height;
    double maxHeight = (*minMaxElem.second).height;

    //for testing, build subsampleCloudIds
    int nbCores = getNumCores();
    std::vector<std::thread> ts;
    for(int i = 0; i < nbCores - 1; i++){
        ts.push_back(std::thread(&DefectSegmentation::computeEquationsMultiThread, this, i, nbCores, kdtree, minHeight, maxHeight, subsampleCloudIds, voxelPointMaps));
    }
    computeEquationsMultiThread(nbCores - 1, nbCores, kdtree, minHeight, maxHeight, subsampleCloudIds, voxelPointMaps);
    for(int i = 0; i < nbCores - 1; i++){
        ts[i].join();
    }
    trace.info()<<"finish eq"<<std::endl;
    IOHelper::export2Text(coefficients, "coeffs");
}



void 
DefectSegmentation::computeEquationsMultiThread(int threadId, int nbThread,
        const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, 
        double minHeight, double maxHeight,
        const std::vector<unsigned int> &subsampleCloudIds, const std::vector<std::vector<unsigned int> > &voxelPointMaps){

    /**debug*/
    std::string prefix = "";
    if(subsampleCloudIds.size() == pointCloud.size()){
        prefix = "f-";
    }
    /**end debug*/
    double neighborSearchRadius = arcLength / 5;
std::cout<<"search: "<<neighborSearchRadius<<std::endl;
    double patchAngle = arcLength / radii;
bool isExported = false;
std::vector<unsigned int> voxels;
//debug, random an id for export patch info ( 100 patches)
    for(unsigned int voxelId = threadId; voxelId < voxelPointMaps.size();voxelId+=nbThread){
        //if(isComputed(i)){
            //continue;
        //}
        std::vector<unsigned int> pointsInVoxel = voxelPointMaps[voxelId];
        if(pointsInVoxel.size() <= 0){
            continue;
        }
        unsigned int i = pointsInVoxel[0];
voxels.push_back(i);
        std::cout<<"nb points in vx "<<voxelId<< "  "<< pointsInVoxel.size()<<std::endl;
        if( pointsInVoxel.size() > 100){
            for(unsigned int dummyId: pointsInVoxel){
                std::cout<< pointCloud[dummyId][0]<< "  "<< pointCloud[dummyId][1]<< "  " << pointCloud[dummyId][2]<<"  "<< std::endl;
            }
        }
        Z3i::RealPoint currentPoint = pointCloud.at(i);
        CylindricalPoint mpCurrent = myPoints.at(i);
        pcl::PointXYZ searchPoint(currentPoint[0], currentPoint[1], currentPoint[2]);
        std::vector<int> pointIdx;
        std::vector<float> pointRadiusSquaredDistance;

        std::vector<double> radiusForEstimate;
        std::vector<double> lengthForEstimate;

        double searchRadius = patchHeight / 2 + 1;
        if(mpCurrent.height - minHeight < patchHeight /2){
            searchRadius += mpCurrent.height - minHeight;
        }else if(maxHeight - mpCurrent.height < patchHeight/2){
            searchRadius += maxHeight - mpCurrent.height;
        }

        //store neighbors in neighborSearchRadius
        std::vector<unsigned int> nearestNeighbors;
        if ( kdtree.radiusSearch (searchPoint, searchRadius, pointIdx, pointRadiusSquaredDistance) > 0 ){
            for (unsigned int idx = 0; idx < pointIdx.size (); ++idx){
                //index of dgtal and pcl is the same
                unsigned int foundedIndex = pointIdx.at(idx);
                unsigned int gIndex = subsampleCloudIds[foundedIndex];
                Z3i::RealPoint found = pointCloud.at(gIndex);
                CylindricalPoint mpFound = myPoints.at(gIndex);
                double angleDiff = std::abs(mpFound.angle - mpCurrent.angle);
                if(angleDiff > patchAngle/2 && 2*M_PI - angleDiff > patchAngle / 2){
                    continue;
                }
                /*
                if(pointRadiusSquaredDistance[idx] < neighborSearchRadius){
                    nearestNeighbors.push_back(foundedIndex);
                }
                */
                radiusForEstimate.push_back(mpFound.radius);
                lengthForEstimate.push_back(mpFound.height);
            }
        }
        //coefficients[i] = Regression::robustLinearOls(lengthForEstimate, radiusForEstimate);
        //coefficients2[i] = coefficients[i];// = Regression::robustLinearOls(lengthForEstimate, radiusForEstimate);
        //std::lock_guard<std::recursive_mutex> lock(coeffMutex);
        
        coefficients[i] = Regression::linearRegression(lengthForEstimate, radiusForEstimate);
        if( pointsInVoxel.size() > 1){
            for(unsigned int pId: pointsInVoxel){
                coefficients[pId] = coefficients[i];
            }
        }
//std::cout<< "nb Point:"<< lengthForEstimate.size()<< std::endl;
if(  lengthForEstimate.size() == 0 ){
    trace.info()<< searchPoint<<std::endl;
    trace.info()<< searchRadius<<std::endl;
}

/**debug
        std::pair<double,double> ran = ransacLine(lengthForEstimate, radiusForEstimate, 2, 1, 5, 30);
        std::pair<double,double> ran1 = ransacLine(lengthForEstimate, radiusForEstimate);
        std::pair<double, double> coli = Regression::ransac(lengthForEstimate, radiusForEstimate, 3, 100, 0.3, 1);
        coefficients[i] = ran;
        //std::pair<double, double> ran1 = ran;
        //std::pair<double, double> coli = ran;
//std::cout<<"skip "<<neaestNeighbors.size()<<std::endl;
std::cout<<"linear coeff:"<<coefficients[i].first<<"  "<<coefficients[i].second<<std::endl;
std::cout<<"my ran coeff:"<<coli.first<<"  "<<coli.second<<std::endl;
std::cout<<"ransac1 coeff:"<<ran.first<<"  "<<ran.second<<std::endl;
std::cout<<"ransac2 coeff:"<<ran1.first<<"  "<<ran1.second<<std::endl;
std::cout<< "est radii 1 ln: "<< coefficients[i].first*lengthForEstimate[0] + coefficients[i].second << " vs "<<radiusForEstimate[0] <<std::endl;
std::cout<< "est radii 1 r1: "<<ran.first*lengthForEstimate[0] + ran.second << " vs "<<radiusForEstimate[0] <<std::endl;
std::cout<< "est radii 1 r2: "<<coli.first*lengthForEstimate[0] + coli.second << " vs "<<radiusForEstimate[0] <<std::endl;
std::cout<< "est radii n ln: "<< coefficients[i].first*lengthForEstimate[lengthForEstimate.size() - 1] + coefficients[i].second <<" vs "<< radiusForEstimate[lengthForEstimate.size() - 1]<<std::endl;
std::cout<< "est radii n r1: "<< ran.first*lengthForEstimate[lengthForEstimate.size() - 1] + ran.second <<" vs "<< radiusForEstimate[lengthForEstimate.size() - 1]<<std::endl;
std::cout<< "est radii n r2: "<< coli.first*lengthForEstimate[lengthForEstimate.size() - 1] + coli.second <<" vs "<< radiusForEstimate[lengthForEstimate.size() - 1]<<std::endl;
        
//random number
        
        if(isForExport[i]){
            std::string xyFile = prefix + "xy-" + std::to_string(i);
            std::string coeffFile = prefix + "coeff-" + std::to_string(i);
            if (lengthForEstimate.size() == radiusForEstimate.size()){
                IOHelper::export2Text(lengthForEstimate, radiusForEstimate, xyFile);
                std::vector<std::pair<double, double> > coeff1;
                coeff1.push_back(ran);
                coeff1.push_back( ran1);
                coeff1.push_back( coli );
                IOHelper::export2Text(coeff1, coeffFile);

            }
        }
        */ //end debug 

        //  uncomment for skip
        /*
        for(unsigned int nId : nearestNeighbors){
            coefficients[nId] = coefficients[i];
        }*/
        

if(threadId == 0){
    trace.progressBar(voxelId , voxelPointMaps.size());
    trace.info()<< voxelId << "/" << voxelPointMaps.size()<<std::endl;
}
    }
IOHelper::export2Text(pointCloud, voxels, "voxels.xyz");
}


void DefectSegmentation::computeDistances(){

    for(unsigned int i = 0; i < myPoints.size(); i++){
        std::pair<double, double> coeffs = coefficients[i];
        if(coeffs.second == 0.0){
        distances[i] = 0;
        }else{
            double estimateRadii = myPoints[i].height * coeffs.first + coeffs.second; 
            distances[i] = myPoints[i].radius - estimateRadii;
        }
    }
    IOHelper::export2Text(distances, "distances");
}

std::pair<double, double> DefectSegmentation::getCoeffs(unsigned int pointId){
    return coefficients.at(pointId);
}

std::pair<double, double> DefectSegmentation::getCoeffs2(unsigned int pointId){
    return coefficients.at(pointId);
}

int DefectSegmentation::computePatch(unsigned int pointId, double searchRadius, const pcl::KdTreeFLANN<pcl::PointXY> &qtree1, 
        const pcl::KdTreeFLANN<pcl::PointXY> &qtree2, std::vector<int> &pointIdx, std::vector<float> &pointRadiusSquaredDistance){
    CylindricalPoint mpFound = myPoints.at(pointId);
    double angle = mpFound.angle;

    pcl::PointXY pxy;
    pxy.y = mpFound.height;
    pxy.x = getScaledX(angle);

    if(angle < angleLimit/2 || angle > 2*M_PI - angleLimit/2){
        return qtree2.radiusSearch (pxy, searchRadius, pointIdx, pointRadiusSquaredDistance);
    }
    return qtree1.radiusSearch (pxy, searchRadius, pointIdx, pointRadiusSquaredDistance);
}

bool DefectSegmentation::isComputed(int pointIndex){
    //std::lock_guard<std::recursive_mutex> lock(coeffMutex);
    bool computed = true;
//std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
//std::cout<< coefficients[pointIndex].first<< coefficients[pointIndex].second <<std::endl;
//std::cout<< coeffZero.first<<"  "<<coeffZero.second <<std::endl;
//std::cout<< (coefficients[pointIndex] == coeffZero ) <<std::endl;
//std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
    if(coefficients[pointIndex] == coeffZero ){
        computed = false;
    }
    return computed;
}

double DefectSegmentation::getScaledX(double angle){
    if(isInAlternativeTree(angle)){
        if( angle < angleLimit ){
            angle += angleLimit;
        }else{
            angle = 2*M_PI - angle;
        }
    }
    return angle * radii * yxRatio;
}

bool DefectSegmentation::isInAlternativeTree(double angle){
    return (angle < angleLimit || angle > 2*M_PI - angleLimit);
}


void DefectSegmentation::computeAngleLimit(){
    angleLimit = arcLength / radii;
    yxRatio = patchHeight / arcLength;
}



std::pair<double, double> DefectSegmentation::ransacLine(const std::vector<double> &xs, const std::vector<double> &ys){
    assert(xs.size() == ys.size());
    //point cloud with x,y, 0
    //
    pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud (new pcl::PointCloud<pcl::PointXYZ>);

    for (int i = 0; i < xs.size(); i++){
        pcl::PointXYZ p(xs[i], ys[i], 0.0);
        pointCloud->push_back(p);
    }

    pcl::SampleConsensusModelLine<pcl::PointXYZ>::Ptr lineModel(new pcl::SampleConsensusModelLine<pcl::PointXYZ>(pointCloud)); 

    //lineModel->setRadiusLimits(0.001,10000); 

    pcl::RandomSampleConsensus<pcl::PointXYZ> sac(lineModel, 3);   
    sac.setMaxIterations(100000);
    sac.setProbability(0.3);
    sac.computeModel();

    Eigen::VectorXf coeff; 
    sac.getModelCoefficients(coeff); 


    std::cout<< "######################################"<<std::endl;
    std::cout<< coeff[0]<< "  "<< coeff[1] << "  "<< coeff[2] << "  "<< coeff[3] <<"  "<< coeff[4] << "  "<<coeff[5] << std::endl;
    std::cout<< "######################################"<<std::endl;
    
    //a = y1/x1 = coeff[4] / coeff[3]
    double a =  coeff[4] / coeff[3];
    //b = y0 -x0y1/x1
    double b =  coeff[1] - coeff[0]*coeff[4]/coeff[3];
    std::pair<double, double> line2D = std::make_pair(a, b);
    return line2D;
}

std::pair<double, double> DefectSegmentation::ransacLine(const std::vector<double> &xs, const std::vector<double> &ys, double threshold, double resolution, double maxPoints, double ep){
    assert(xs.size() == ys.size());
    //point cloud with x,y, 0
    //
    // segment by x
    auto xMinMax = std::minmax_element(xs.begin(), xs.end());
    double xMin = xs.at( xMinMax.first - xs.begin() );
    double xMax = xs.at( xMinMax.second - xs.begin() );
    //double ep = arcLength / 5.0;
    int nbSegment = (xMax - xMin) / resolution + 1;

    //element with min radii
    std::vector<double> minElems(nbSegment, std::numeric_limits<double>::max());
    std::vector<double> minElemIds(nbSegment, -1);
    for(int i = 0; i < xs.size(); i++){
        int segId = (xs[i] - xMin) / resolution;
        if(ys[i] < minElems[segId]){
            minElems[segId] = ys[i];
            minElemIds[segId] = i;
        }
    }
    

    pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud (new pcl::PointCloud<pcl::PointXYZ>);

    std::vector<int> nbPoints(nbSegment, 0);
    for (int i = 0; i < xs.size(); i++){
        int segId = (xs[i] - xMin) / resolution;
        if(minElemIds[segId] >= 0 && nbPoints[segId] < maxPoints){
            if(ys[i] - minElems[segId] < ep){
                nbPoints[segId]++;
                pcl::PointXYZ p(xs[i], ys[i], 0.0);
                pointCloud->push_back(p);
            }
           // int pId = minElemIds[i];
            //pcl::PointXYZ p(xs[pId], ys[pId], 0.0);
            //pointCloud->push_back(p);
        }
    }
    
    /*
    for (int i = 0; i < minElemIds.size(); i++){
        //int segId = (xs[i] - xMin) / resolution;
        if(minElemIds[i] >= 0){
            int pId = minElemIds[i];
            pcl::PointXYZ p(xs[pId], ys[pId], 0.0);
            pointCloud->push_back(p);
        }
    }
    */

    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
    pcl::ModelCoefficients::Ptr coeffs(new pcl::ModelCoefficients);

    Eigen::Vector3f axis(0,0,1);
    pcl::SACSegmentation<pcl::PointXYZ> seg;

    seg.setAxis(axis);
    seg.setEpsAngle(0.5);
    seg.setOptimizeCoefficients (true);
    seg.setModelType (pcl::SACMODEL_LINE);
    seg.setMethodType (pcl::SAC_RANSAC);
    seg.setDistanceThreshold (threshold);
    //limit ?
    //seg.setRadiusLimits (radius * 0.7, radius * 1.3);

    seg.setInputCloud (pointCloud);
    seg.segment (*inliers, *coeffs);
    /*
    pcl::SampleConsensusModelLine<pcl::PointXYZ>::Ptr lineModel(new pcl::SampleConsensusModelLine<pcl::PointXYZ>(pointCloud)); 

    //lineModel->setRadiusLimits(0.001,10000); 

    pcl::RandomSampleConsensus<pcl::PointXYZ> sac(lineModel, 3);   
    sac.setMaxIterations(100000);
    sac.setProbability(0.3);
    sac.computeModel();

    Eigen::VectorXf coeff; 
    sac.getModelCoefficients(coeff); 
    */

    double x0 = coeffs->values[0];
    double y0 = coeffs->values[1];
    double z0 = coeffs->values[2];
    double x1 = coeffs->values[3];
    double y1 = coeffs->values[4];
    double z1 = coeffs->values[5];

    //a = y1/x1 = coeff[4] / coeff[3]
    double a =  y1 / x1;
    //b = y0 -x0y1/x1
    double b =  y0 - x0*y1/x1;
    std::pair<double, double> line2D = std::make_pair(a, b);
    return line2D;
}


