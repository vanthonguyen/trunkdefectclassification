#include <iostream>
#include <fstream>
#include <utility>
#include <math.h>
#include <cmath>
#include <stdlib.h>


#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/segmentation/extract_clusters.h>



#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/PointListReader.h"
#include <DGtal/math/Statistic.h>



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/shapes/Mesh.h"

#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

#include "../Common/IOHelper.h"
#include "../Common/Rosin.h"
#include "../Segmentation/EuclideCluster.h"
#include "../Centerline/Centerline.h"


using namespace DGtal;
namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;


int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input pointCloud.")
        ("searchRadius,R", po::value<double>()->default_value(10), "search for neighbor in radius.")
        ("searchRadiusSmall,r", po::value<double>()->default_value(3), "search for neighbor in a smaller radius.")
        ("nbPoint,n", po::value<int>()->default_value(1), "nb neighbor in radius.")
        ("tolerance", po::value<double>()->default_value(10), "min distance between two clusters")
        ("clusterSize", po::value<int>()->default_value(30), "min points in a cluster")
        ("output,o", po::value<std::string>()->default_value("output.xyz"), "output.xyz");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }

    po::notify(vm);
    if(vm.count("help") || argc<=1 || !parseOK || !vm.count("input") || !vm.count("searchRadius") ) {
        if(!vm.count("input")){
            trace.error()<<"the input mesh is required!"<<std::endl;
        }else if( !vm.count("searchRadius") ){
            trace.error()<<"the searchRadius is required!"<<std::endl;
        }
        trace.info()<< "Point cloud denoising" <<std::endl << "Options: "<<std::endl
            << general_opt << "\n";
        return 0;
    }

    if(!vm.count("nbPoint") ) {
        trace.warning()<<"nb of neighbors was not specified, used default value 1 !"<<std::endl;
    }

    if(!vm.count("output") ) {
        trace.warning()<<"ouput filename was not specified, used default value output.xyz !"<<std::endl;
    }

    
    
    int nbPoint = vm["nbPoint"].as<int>();
    double tolerance = vm["tolerance"].as<double>();
    int minClusterSize = vm["clusterSize"].as<int>();
    std::cout<<nbPoint<<std::endl;

    //double binWidth = vm["binWidth"].as<double>();
    double searchRadiusSmall = vm["searchRadiusSmall"].as<double>();
    double searchRadius = vm["searchRadius"].as<double>();

    std::vector<Z3i::RealPoint> points = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["input"].as<std::string>());

std::cout<<"begin"<<std::endl;
    //estimate normal for all points
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    for(Z3i::RealPoint &p: points){
        //scale down
        //
        //p /= voxelSize;
        pcl::PointXYZ pxyz(p[0], p[1], p[2]);
        cloud->points.push_back(pxyz);
    }
        
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud (cloud);
    // Create the normal estimation class, and pass the input dataset to it
     // Create the normal estimation class, and pass the input dataset to it
    std::vector<int> pointIdx;
    std::vector<float> pointRadiusSquaredDistance;
    //for debug
    //Search for neighbours
    //for(Z3i::RealPoint p: points){
    std::vector<Z3i::RealPoint> filtered;
    std::cout<<"begin filter"<<std::endl;
    std::vector<Z3i::RealPoint> noises;
    for (int pId = 0; pId < cloud->points.size(); pId++){
        //std::cout << "processing point" << pId <<"/" << cloud->points.size()<<std::endl;
        pcl::PointXYZ searchPoint = cloud->points[pId];
        trace.progressBar(pId, points.size());
        if ( kdtree.radiusSearch (searchPoint, searchRadius, pointIdx, pointRadiusSquaredDistance) >= nbPoint ) {
            filtered.push_back(points[pId]);
        }else{
            noises.push_back(points[pId]);
        }
    }

    //
    std::cout<<"remove "<< points.size() - filtered.size() << "  points by density" <<std::endl;  

    std::string fileName = vm["output"].as<std::string>();
    if(tolerance > 0){

        std::vector<unsigned int> indices = EuclideCluster::removeSmallCluster(filtered, tolerance, minClusterSize);

        std::cout<<"remove "<< filtered.size() - indices.size() << "  points by clustering" <<std::endl;  

        IOHelper::export2Text(filtered, indices, fileName);
    }else{
        IOHelper::export2Text(filtered, fileName);
    }
    //IOHelper::export2Text(filtered, fileName);
    IOHelper::export2Text(noises, "noises-" + fileName);

    return 0;
}
