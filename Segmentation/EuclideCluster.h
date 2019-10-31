#ifndef EUCLIDE_CLUSTER_H
#define EUCLIDE_CLUSTER_H

#include "pcl/common/common.h"
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>

#include <iostream>
#include <vector>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

using namespace DGtal;

#define MINIMUM_DEFAULT_SIZE 20
class EuclideCluster{
public:
    static std::vector<unsigned int> removeSmallCluster(const std::vector<Z3i::RealPoint> &aCloud, double tolerance = 8, int minClusterSize = 20){
        std::vector<unsigned int> newDefs;
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
        std::vector<pcl::PointIndices> clusterIndices;
        //Z3i to pcl
        for(int i = 0; i < aCloud.size(); i++){
            Z3i::RealPoint p = aCloud[i];
            cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
        }
        pcl::search::KdTree<pcl::PointXYZ>::Ptr kd (new pcl::search::KdTree<pcl::PointXYZ>);
        kd->setInputCloud (cloudPcl);

        pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
        ec.setClusterTolerance (tolerance);
        ec.setMinClusterSize (minClusterSize);
        ec.setMaxClusterSize (cloudPcl->points.size());
        ec.setSearchMethod (kd);
        ec.setInputCloud (cloudPcl);
        ec.extract (clusterIndices);
        //covert pcl clusterIndices to std::vector
        for (std::vector<pcl::PointIndices>::const_iterator it = clusterIndices.begin (); it != clusterIndices.end (); it++){
            //clusters are ordered by the number of points
            std::vector<unsigned int> aCluster;
            for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); pit++){
                newDefs.push_back(*pit);
            }
        }
        return newDefs;
    }

    static std::vector<unsigned int> removeSmallCluster(const std::vector<Z3i::RealPoint> &aCloud, const std::vector<unsigned int> &defectIndices, double  tolerance = 8, int minClusterSize = 20){
        std::vector<unsigned int> newDefs;
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointIndices::Ptr defectIndicesPcl (new pcl::PointIndices);
        std::vector<pcl::PointIndices> clusterIndices;
        //Z3i to pcl
        for(int i = 0; i < aCloud.size(); i++){
            Z3i::RealPoint p = aCloud[i];
            cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
        }
        //std vector to pcl indices
        for( size_t i = 0; i < defectIndices.size(); i++)
        {
            defectIndicesPcl->indices.push_back(defectIndices[i]);
        } 

        pcl::search::KdTree<pcl::PointXYZ>::Ptr kd (new pcl::search::KdTree<pcl::PointXYZ>);
        kd->setInputCloud (cloudPcl);

        pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
        ec.setClusterTolerance (tolerance);
        ec.setMinClusterSize (minClusterSize);
        ec.setMaxClusterSize (cloudPcl->points.size());
        ec.setSearchMethod (kd);
        ec.setInputCloud (cloudPcl);
        ec.setIndices(defectIndicesPcl);
        ec.extract (clusterIndices);
        //covert pcl clusterIndices to std::vector
        for (std::vector<pcl::PointIndices>::const_iterator it = clusterIndices.begin (); it != clusterIndices.end (); it++)
        {
            //clusters are ordered by the number of points
            if(it->indices.size() > MINIMUM_DEFAULT_SIZE)
            {
                std::vector<unsigned int> aCluster;
                for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); pit++)
                {
                    newDefs.push_back(*pit);
                }
            }
        }
        return newDefs;
    }

    static std::vector<std::vector<unsigned int> > 
    cluster(const std::vector<Z3i::RealPoint> &aCloud, const std::vector<unsigned int> &defectIndices, double tolerance = 8, int minClusterSize = 20){
        std::vector<std::vector<unsigned int> > clusters;
        trace.info()<<"EuclideCluster->begin cluster"<< std::endl;
        trace.info()<<std::endl;
        trace.info()<<"EuclideCluster->using PCL ....";
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointIndices::Ptr defectIndicesPcl (new pcl::PointIndices);
        std::vector<pcl::PointIndices> clusterIndices;
        //Z3i to pcl
        for(int i = 0; i < aCloud.size(); i++){
            Z3i::RealPoint p = aCloud[i];
            cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
        }
        //std vector to pcl indices
        for( size_t i = 0; i < defectIndices.size(); i++)
        {
            defectIndicesPcl->indices.push_back(defectIndices[i]);
        } 

        pcl::search::KdTree<pcl::PointXYZ>::Ptr kd (new pcl::search::KdTree<pcl::PointXYZ>);
        kd->setInputCloud (cloudPcl);

        pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
        ec.setClusterTolerance (tolerance);
        ec.setMinClusterSize (minClusterSize);
        ec.setMaxClusterSize (cloudPcl->points.size());
        ec.setSearchMethod (kd);
        ec.setInputCloud (cloudPcl);
        ec.setIndices(defectIndicesPcl);
        ec.extract (clusterIndices);
        //covert pcl clusterIndices to std::vector
        for (std::vector<pcl::PointIndices>::const_iterator it = clusterIndices.begin (); it != clusterIndices.end (); it++)
        {
            //clusters are ordered by the number of points
            std::vector<unsigned int> aCluster;
            for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); pit++)
            {
                aCluster.push_back(*pit);
            }
            clusters.push_back(aCluster);
        }
        clusters.push_back(defectIndices);
        return clusters;
    }

    static std::vector<std::vector<unsigned int> > 
    clusterPcl(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl, const std::vector<unsigned int> &defectIndices , double tolerance = 8, int minClusterSize = 20){
        std::vector<std::vector<unsigned int> > clusters;
        std::vector<pcl::PointIndices> clusterIndices;
        pcl::PointIndices::Ptr defectIndicesPcl (new pcl::PointIndices);

        //std vector to pcl indices
        for( size_t i = 0; i < defectIndices.size(); i++)
        {
            defectIndicesPcl->indices.push_back(defectIndices[i]);
        } 

        pcl::search::KdTree<pcl::PointXYZ>::Ptr kd (new pcl::search::KdTree<pcl::PointXYZ>);
        kd->setInputCloud (cloudPcl);

        pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
        ec.setClusterTolerance (tolerance);
        ec.setMinClusterSize (minClusterSize);
        ec.setMaxClusterSize (cloudPcl->points.size());
        ec.setSearchMethod (kd);
        ec.setInputCloud (cloudPcl);
        ec.setIndices(defectIndicesPcl);
        ec.extract (clusterIndices);
        //covert pcl clusterIndices to std::vector
        for (std::vector<pcl::PointIndices>::const_iterator it = clusterIndices.begin (); it != clusterIndices.end (); it++)
        {
            //clusters are ordered by the number of points
            std::vector<unsigned int> aCluster;
            for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); pit++)
            {
                aCluster.push_back(*pit);
            }
            clusters.push_back(aCluster);
        }
        clusters.push_back(defectIndices);
        return clusters;
    }

};
#endif //EUCLIDE_CLUSTER_H
