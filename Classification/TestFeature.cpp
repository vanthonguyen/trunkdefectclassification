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

#include "../IOHelper.h"
#include "../FeatureExtraction.h"


using namespace DGtal;

namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input file name a point cloud in text xyz format.");

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
    std::vector<CylindricalPoint> cpoints;
    for (Z3i::RealPoint p: pointCloud){
        CylindricalPoint cp(p[0], p[1], p[2]);
        cpoints.push_back(cp);
    }

    FeatureExtraction ef(cpoints);
    std::vector<float> feature = ef.compute();
    for(int i = 0; i < feature.size() - 1; i++){
        std::cout<< feature[i]<<" ";
    }

    std::cout<< feature[feature.size() - 1]<<std::endl;
    return 0;
}

