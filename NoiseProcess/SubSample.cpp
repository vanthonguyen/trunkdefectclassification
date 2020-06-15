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


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/shapes/Mesh.h"

#include "../Common/IOHelper.h"
#include "../Segmentation/SegmentationHelper.h"


using namespace DGtal;
namespace po = boost::program_options;


int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input pointCloud.")
        ("resolution,r", po::value<double>()->default_value(1), " in a smaller radius.")
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

    
    double resolution  = vm["resolution"].as<double>();

    std::vector<Z3i::RealPoint> points = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["input"].as<std::string>());

    std::vector<Z3i::RealPoint> sampledCloud;

    //radius???
    SegmentationHelper::simpleSubSample(points, resolution, sampledCloud);

    
    std::string outputFile = vm["output"].as<std::string>();
    //IOHelper::export2Text(filtered, fileName);
    IOHelper::export2Text(sampledCloud, outputFile); 

    return 0;
}
