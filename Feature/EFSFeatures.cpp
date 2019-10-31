#include "EFSFeatures.h"

void compute(pcl::PointCloud<pcl::ESFSignature640>::Ptr &descriptor){
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    IOHelperPcl::loadXyz(fileName, cloud);
    computeEstimatorESF(cloud, descriptor);
}

void EFSFeatures::computeEstimatorESF(const pcl::PointCloud<pcl::PointXYZ>::Ptr &pointCloud, pcl::PointCloud<pcl::ESFSignature640>::Ptr &descriptor){
    // ESF estimation object.
    pcl::ESFEstimation<pcl::PointXYZ, pcl::ESFSignature640> esf;
    esf.setInputCloud(pointCloud);
    esf.compute(*descriptor);
}
