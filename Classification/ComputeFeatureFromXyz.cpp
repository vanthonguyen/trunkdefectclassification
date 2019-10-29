#include <iostream>
#include <fstream>
#include <utility>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/shapes/Mesh.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "../Common/IOHelper.h"
#include "../Feature/FeatureExtractionXyz.h"
#include "../Segmentation/SegmentationHelper.h"


using namespace DGtal;

namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input file name a point cloud in text xyz format.")
        ("resolution,r", po::value<double>()->default_value(3), "Resolution for subsample point cloud")
        ("segmentRadius,s", po::value<std::string>(), "radius of segments.")
        ("centerline,l", po::value<std::string>(), "centerline");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }
    //better if use custom reader ?
    std::vector<Z3i::RealPoint> pointCloud = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["input"].as<std::string>());
    
    double resolution = vm["resolution"].as<double>();
    //subsample
    std::vector<Z3i::RealPoint> sampledCloud;
    SegmentationHelper::simpleSubSample(pointCloud, resolution, sampledCloud);
    std::vector<Z3i::RealPoint> centerline = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["centerline"].as<std::string>());
    std::vector<double> segRadiis;
    IOHelper::readDistanceFromFile(vm["segmentRadius"].as<std::string>(), segRadiis);

    Z3i::RealPoint defectLowPoint;
    Z3i::RealPoint defectHighPoint;

    //use points with depth < 20mm to compute center of defect
    double maxDepth = 10;

    double minR = std::numeric_limits<double>::max();

    bool first = true;
    std::vector<CylindricalPoint> cylPoints(pointCloud.size());

    Z3i::RealPoint seedPoint(0,1,0);
    CylindricalCoordinateSystem ccs(centerline, seedPoint);
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint xyzPoint = pointCloud[i];
        cylPoints[i] = ccs.xyz2Cylindrical(xyzPoint);

//std::cout<<"#####"<< cylPoints[i].radius << "#"<< cylPoints[i].angle<<"#"<<cylPoints[i].height<<std::endl;
        if(cylPoints[i].radius < minR){
            minR = cylPoints[i].radius;
        }
    }
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        CylindricalPoint cylPoint = cylPoints[i];

        if(cylPoint.radius - minR < maxDepth){
            Z3i::RealPoint xyzPoint = pointCloud[i];
            if(first){
                defectHighPoint = xyzPoint;
                defectLowPoint = xyzPoint;
                first = false;
            }else{
                defectLowPoint = defectLowPoint.inf(xyzPoint);
                defectHighPoint = defectHighPoint.sup(xyzPoint);
            }
        }
    }

    Z3i::RealPoint defectCenter = (defectLowPoint + defectHighPoint)/2;
    CylindricalPoint cylDefectCenter = ccs.xyz2Cylindrical(defectCenter);
    double trunkRadiusAtDefect = segRadiis[cylDefectCenter.segmentId];

    FeatureExtractionXyz ef(sampledCloud, centerline, trunkRadiusAtDefect);
    std::vector<float> feature = ef.compute();
    for(unsigned int i = 0; i < feature.size() - 1; i++){
        std::cout<< feature[i]<<" ";
    }

    std::cout<< feature[feature.size() - 1]<<std::endl;
    return 0;
 }

