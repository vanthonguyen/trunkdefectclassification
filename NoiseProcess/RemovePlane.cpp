#include <iostream>
#include <fstream>
#include <utility>
#include <math.h>
#include <cmath>
#include <stdlib.h>


#include <pcl/common/common.h>
#include <pcl/point_cloud.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/PointListReader.h"



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "../Common/IOHelper.h"


using namespace DGtal;
namespace po = boost::program_options;

int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input pointCloud.")
        ("arcLength,a", po::value<double>()->default_value(3), "arcLength")
        ("radius,r", po::value<double>()->default_value(200), "estimated radius of grume")
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
    if(vm.count("help") || argc<=1 || !parseOK || !vm.count("input") ) {
        if(!vm.count("input")){
            trace.error()<<"the input mesh is required!"<<std::endl;
        }
        trace.info()<< "Point cloud denoising" <<std::endl << "Options: "<<std::endl
            << general_opt << "\n";
        return 0;
    }

    if(!vm.count("output") ) {
        trace.warning()<<"ouput filename was not specified, used default value output.xyz !"<<std::endl;
    }

    //double binWidth = vm["binWidth"].as<double>();
    double length = vm["arcLength"].as<double>();
    double radius = vm["radius"].as<double>();

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


  pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
  pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
  // Create the segmentation object
  pcl::SACSegmentation<pcl::PointXYZ> seg;
  // Optional
  seg.setOptimizeCoefficients (true);
  // Mandatory
  seg.setModelType (pcl::SACMODEL_PLANE);
  seg.setMethodType (pcl::SAC_RANSAC);
  seg.setDistanceThreshold (5);

  seg.setInputCloud (cloud);
  seg.segment (*inliers, *coefficients);

  std::vector<unsigned int> plane;
  std::vector<unsigned int> compplane;
    std::vector<bool> onPlane(points.size(), false);

  for (int i = 0; i < inliers->indices.size(); i++){
    int ind = inliers->indices.at(i);
    onPlane[ind] = true;
    plane.push_back(ind);
  }
  for(int i = 0; i < onPlane.size(); i++){
      bool isOnPlane = onPlane[i];
      if( !isOnPlane){
          compplane.push_back(i);
      }
  }

    std::string fileName = vm["output"].as<std::string>();
    IOHelper::export2Text(points, plane, "plane-"+fileName);
    IOHelper::export2Text(points, compplane, "compplane-"+fileName);

    return 0;
}
