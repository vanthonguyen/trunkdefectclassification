#ifndef EFS_FEATURES_H 
#define EFS_FEATURES_H

#include<vector>

#include <pcl/point_types.h>
#include <pcl/features/esf.h>

#include "../IOHelperPcl.h"



class EFSFeatures{
public:

    //EFSFeatures(std::string f): fileName(f){}
    void static compute(const std::string &filename, pcl::PointCloud<pcl::ESFSignature640>::Ptr &descriptor){
        pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud(new pcl::PointCloud<pcl::PointXYZ>);
        IOHelperPcl::loadXyz(filename, pointCloud);
        compute(pointCloud, descriptor);
    }

    void static compute(const pcl::PointCloud<pcl::PointXYZ>::Ptr &pointCloud, pcl::PointCloud<pcl::ESFSignature640>::Ptr &descriptor){
        pcl::ESFEstimation<pcl::PointXYZ, pcl::ESFSignature640> esf;
        esf.setInputCloud(pointCloud);
        esf.compute(*descriptor);
    }
    //
//    static double getMedian(std::vector<double> v);
//private:


//    void computeEstimatorESF(const pcl::PointCloud<pcl::PointXYZ>::Ptr &pointCloud, pcl::PointCloud<pcl::ESFSignature640>::Ptr &descriptor);
//    std::string fileName;

};

#endif // STATISTIC_H
