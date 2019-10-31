#ifndef PCA_CYLINDRICAL_H
#define PCA_CYLINDRICAL_H

#include<vector>

#include <pcl/common/common.h>
#include <pcl/point_cloud.h>

#include <pcl/common/pca.h>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "../Common/RlzPoint.h"


using namespace DGtal;

class Pca{
public:

    static void compute(const std::vector<RlzPoint> & points, Eigen::Matrix3f &eigenVectors, Eigen::Vector3f &eigenValues){
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
        for(int i = 0; i < points.size(); i++){
            RlzPoint p = points[i];
            cloudPcl->points.push_back(pcl::PointXYZ(p.r, p.l, p.z));
        }
        pcl::PCA<pcl::PointXYZ> pca;
        pca.setInputCloud (cloudPcl);
        eigenVectors = pca.getEigenVectors ();
        eigenValues = pca.getEigenValues ();
    }

    static void compute(const std::vector<Z3i::RealPoint> & points, Eigen::Matrix3f &eigenVectors, Eigen::Vector3f &eigenValues){
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
        for(int i = 0; i < points.size(); i++){
            Z3i::RealPoint p = points[i];
            cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
        }
        pcl::PCA<pcl::PointXYZ> pca;
        pca.setInputCloud (cloudPcl);
        eigenVectors = pca.getEigenVectors ();
        eigenValues = pca.getEigenValues ();
    }

};

#endif // PCA_CYLINDRICAL_H
