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
//#include "../TrunkBranchSegmentation.h"
#include "../Common/CylindricalCoordinateSystem.h"


using namespace DGtal;
namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

int getNearestPoint(const Z3i::RealPoint &p, const std::vector<Z3i::RealPoint> &centerline){
    int id = 0;
    double minDist = 100000000000;
    for(int i = 1; i< centerline.size(); i++){
        Z3i::RealPoint vect = p - centerline[i];
        if( vect.norm() < minDist){
            minDist = vect.norm();
            id = i;
        }
    }
    return id;
}

int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input mesh.")
        ("centerline,l", po::value<std::string>(), "input centerline.")
        ("maxRadius,R", po::value<double>()->default_value(600), "max radius, remove all point whose radius > max radius")
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
    double maxRadius = vm["maxRadius"].as<double>();

    std::vector<Z3i::RealPoint> centerline = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["centerline"].as<std::string>());

    std::string inputMeshName = vm["input"].as<std::string>();
    DGtal::Mesh<Z3i::RealPoint> oriMesh(true);
    MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, oriMesh, false);

    Z3i::RealPoint seedPoint(0,1,0);
    CylindricalCoordinateSystem ccs(centerline, seedPoint);

    std::vector<unsigned int> facesToBeRemoved;
    for (unsigned int i = 0; i < oriMesh.nbFaces(); i++){

        Z3i::RealPoint xyzPoint = oriMesh.getFaceBarycenter(i);
        CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);

        if(cylPoint.radius > maxRadius){
            facesToBeRemoved.push_back(i);
        }
    }

    oriMesh.removeFaces(facesToBeRemoved);

    std::string fileName = vm["output"].as<std::string>();
    IOHelper::export2OFF(oriMesh, fileName);
//    IOHelper::export2Text(points, noiseIds, "noises-"+fileName);

    return 0;
}
