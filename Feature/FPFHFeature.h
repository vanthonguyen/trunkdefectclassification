#ifndef FPFH_FEATURES_H 
#define FPFH_FEATURES_H

#include<vector>

#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/fpfh.h>




class FPFHFeature{
public:

    void static 
        compute(const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, double normalRadiusSearch, 
                double fpfhRadiusSearch,
                pcl::PointCloud<pcl::FPFHSignature33>::Ptr &descriptors){
        // Object for storing the normals.
        pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);

        // Estimate the normals.
        pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
        normalEstimation.setInputCloud(cloud);
        normalEstimation.setRadiusSearch(normalRadiusSearch);
        pcl::search::KdTree<pcl::PointXYZ>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZ>);
        normalEstimation.setSearchMethod(kdtree);
        normalEstimation.compute(*normals);
        std::cout<<normals->points.size()<<std::endl;
        // FPFH estimation object.
        pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh;
        fpfh.setInputCloud(cloud);
        fpfh.setInputNormals(normals);
        fpfh.setSearchMethod(kdtree);
        // Search radius, to look for neighbors. Note: the value given here has to be
        // larger than the radius used to estimate the normals.
        fpfh.setRadiusSearch(fpfhRadiusSearch);

        fpfh.compute(*descriptors);
    }
    //
//    static double getMedian(std::vector<double> v);
};

#endif // STATISTIC_H
