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
#include "../Common/Rosin.h"
#include "../Common/IOHelper.h"
#include "../Common/MultiThreadHelper.h"
#include "../Segmentation/EuclideCluster.h"

#include "DefectBranchSegmentation.h"

using namespace DGtal;

void
DefectBranchSegmentation::init(){
    allocate();
    allocateExtra();

    computeBeginOfSegment();
    computeVectorMarks();
    computePlaneNormals();

	convertToCcs();

    computeAngleLimit();

    computeEquations();

    computeDistances();
}

void
DefectBranchSegmentation::segment(std::vector<std::vector<unsigned int> > &clusters){

    std::vector<unsigned int> defectIds = getDefect();
    //clustering
    clusters = EuclideCluster::cluster(pointCloud, defectIds, clusterTolerance, minClusterSize);
}


std::vector<unsigned int>
DefectBranchSegmentation::getDefect(){
    std::vector<unsigned int> defects;

    for(unsigned int globalPid: trunkIds){
        if(distances[globalPid] > threshold){
            defects.push_back(globalPid);
        }
    }

    for(unsigned int globalPid: branchIds){
        defects.push_back(globalPid);
    }

    //IOHelper::export2Text(branchIds, "branchIds.id");
    return defects;

}


void DefectBranchSegmentation::computeDistances(){
    /**
     * if point is on branch, distance = R - radii
     * is this stable? distance always > threshold ?
     **/
    for(unsigned int globalPid : branchIds){
        distances[globalPid] = myPoints[globalPid].radius - radii;
    }

    std::numeric_limits<double>::epsilon();
    std::vector<double> distanceOfPointOnTrunk(trunkIds.size(), 0);
std::cout<<"DefectBranchSegmentation->trunkIds.size()"<< trunkIds.size()<<std::endl;
    for(unsigned int i = 0; i < trunkIds.size(); i++){
        unsigned int globalPid = trunkIds[i];
        std::pair<double, double> coeffs = coefficients[globalPid];
        if(coeffs.second == 0.0){
            distances[globalPid] = 0;
        }else{
            double estimateRadii = myPoints[globalPid].height * coeffs.first + coeffs.second; 
            distances[globalPid] = myPoints[globalPid].radius - estimateRadii;
            distanceOfPointOnTrunk[i] = distances[globalPid];
        }
    }
    IOHelper::export2Text(distances, "distances");
    IOHelper::export2Text(distanceOfPointOnTrunk, "distancesTrunk");
std::cout<<"####binWidth"<<binWidth<<std::endl;
    threshold = Rosin::compute(distanceOfPointOnTrunk, binWidth);
}


void 
DefectBranchSegmentation::computeEquations(){

std::cout<<"DefectBranchSegmentation->Begin compute equations"<<std::endl;
    SegmentationHelper::segmentBranchesUsingCenterline(pointCloud, sectorLength, centerline, trunkIds, branchIds);

    //w = patch width, wh = patch height
    //build kdtree using pcl 
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
    //get subsample cloud
    //
    std::vector<Z3i::RealPoint> trunkCloud;
    for(unsigned int pid : trunkIds){
        trunkCloud.push_back(pointCloud[pid]);
    }

    std::vector<unsigned int> subsampleCloudIds;

    //estimate radius
    
    //@FIX?
    double resolution = voxelSize;
    //hardcode, resolution = 5, radius = 180 is not a good idea!!
    SegmentationHelper::subSample(trunkCloud, resolution, radii, subsampleCloudIds);

    //IOHelper::export2Text(subsampleCloud, "subsample.xyz");

    Z3i::RealPoint upPt = trunkCloud[subsampleCloudIds[0]];
    Z3i::RealPoint lowPt = trunkCloud[subsampleCloudIds[0]];
    for(int i = 0; i < subsampleCloudIds.size(); i++){
        Z3i::RealPoint p = trunkCloud.at(subsampleCloudIds.at(i));
        cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
        upPt = upPt.sup(p);
        lowPt = lowPt.inf(p);
    }

    upPt += Z3i::RealPoint(resolution,resolution,resolution);
    lowPt -= Z3i::RealPoint(resolution,resolution,resolution);

    std::map<unsigned int, std::vector<unsigned int> > indexMap;
    SegmentationHelper::voxelize(trunkCloud, resolution, lowPt, upPt, indexMap);

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
        ts.push_back(std::thread(&DefectBranchSegmentation::computeEquationsMultiThread, this, i, nbCores, kdtree, minHeight, maxHeight, subsampleCloudIds, voxelPointMaps));
    }
    computeEquationsMultiThread(nbCores - 1, nbCores, kdtree, minHeight, maxHeight, subsampleCloudIds, voxelPointMaps);
    for(int i = 0; i < nbCores - 1; i++){
        ts[i].join();
    }
    trace.info()<<"finish eq"<<std::endl;
    IOHelper::export2Text(coefficients, "coeffs");
//    IOHelper::export2Text(trunkCloud, "trunk.xyz");
//    IOHelper::export2Text(trunkCloud, subsampleCloudIds, "subtrunk.xyz");
}



void 
DefectBranchSegmentation::computeEquationsMultiThread(int threadId, int nbThread,
        const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, 
        double minHeight, double maxHeight,
        const std::vector<unsigned int> &subsampleCloudIds, 
        const std::vector<std::vector<unsigned int> > &voxelPointMaps){

    double patchAngle = arcLength / radii;
    for(unsigned int voxelId = threadId; voxelId < voxelPointMaps.size();voxelId+=nbThread){
        std::vector<unsigned int> pointsInVoxel = voxelPointMaps[voxelId];
        if(pointsInVoxel.size() <= 0){
            continue;
        }
        unsigned int i = pointsInVoxel[0];
        //std::cout<<"nb points in vx "<<voxelId<< "  "<< pointsInVoxel.size()<<std::endl;
 
        int currentGlobalId = trunkIds.at(i);
        Z3i::RealPoint currentPoint = pointCloud.at(currentGlobalId);
        CylindricalPoint mpCurrent = myPoints.at(currentGlobalId);
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
                /**
                 * Something complicate here
                 **/
                // index on subsampled cloud of trunk without branch
                unsigned int foundedIndex = pointIdx.at(idx);

                // index on trunk
                unsigned int indexOntrunk = subsampleCloudIds[foundedIndex];
                // index on original point cloud (trunk + branch)
                unsigned int gIndex = trunkIds.at(indexOntrunk);
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
        
        coefficients[currentGlobalId] = Regression::linearRegression(lengthForEstimate, radiusForEstimate);
        if( pointsInVoxel.size() > 1){
            for(unsigned int pId: pointsInVoxel){
                unsigned int globalPid = trunkIds[pId];
                coefficients[globalPid] = coefficients[currentGlobalId];
            }
        }
        /*
        if(threadId == 0){
            trace.progressBar(voxelId , voxelPointMaps.size());
            trace.info()<< voxelId << "/" << voxelPointMaps.size()<<std::endl;
        }*/
    }
}



std::vector<unsigned int> 
DefectBranchSegmentation::getTrunkIds(){
    return trunkIds;
}

std::vector<unsigned int> 
DefectBranchSegmentation::getBranchIds(){
    return branchIds;
}


