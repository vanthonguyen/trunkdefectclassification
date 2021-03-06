#include <iostream>
#include <fstream>
#include <utility>


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
#include "../Segmentation/EuclideCluster.h"


using namespace DGtal;
namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

/**
 * Clustering the mesh by euclidean distance
 * Remove all but the biggist cluster
 **/

int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input mesh.")
        ("tolerance", po::value<int>()->default_value(10), "min distance between two clusters")
        ("clusterSize", po::value<int>()->default_value(30), "min points in a cluster")
        ("output,o", po::value<std::string>()->default_value("output"), "output prefix: output-defect.off, output-def-faces-ids, ...");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }

    po::notify(vm);
    if(vm.count("help") || argc<=1 || !parseOK || !vm.count("input")){
        if(!vm.count("input")){
            trace.error()<<"the input mesh is required!"<<std::endl;
        }
        trace.info()<< "Remove clusters: " <<std::endl << "Options: "<<std::endl
            << general_opt << "\n";
        return 0;
    }

    DGtal::Mesh<Z3i::RealPoint> oriMesh(true);
    std::string inputMeshName = vm["input"].as<std::string>();
    MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, oriMesh, false);

    int tolerance = vm["tolerance"].as<int>();
    int minClusterSize = vm["clusterSize"].as<int>();


    std::vector<Z3i::RealPoint> cloud(oriMesh.nbFaces());

    for (unsigned int i = 0; i < oriMesh.nbFaces(); i++){
        cloud[i] = oriMesh.getFaceBarycenter(i);
    }

    std::cout<<"begin remove small clusters"<<std::endl;
    std::cout<< "cloud size: "<< cloud.size()<<std::endl;
    std::vector<unsigned int> indices = EuclideCluster::removeSmallCluster(cloud, tolerance, minClusterSize);
    std::cout<<"end remove small clusters"<<std::endl;

    std::vector<bool> toBeRemoveFlags(cloud.size(), true);
    for(unsigned int i = 0; i < indices.size();i++){
        toBeRemoveFlags[indices[i]] = false;
    }

    std::vector<unsigned int> facesToBeRemoved;
    for (unsigned int i = 0; i < oriMesh.nbFaces(); i++){
        if(toBeRemoveFlags[i]){
            facesToBeRemoved.push_back(i);
        }
    }


    oriMesh.removeFaces(facesToBeRemoved);

    std::cout<<std::endl<< "removed "<< facesToBeRemoved.size()<<std::endl;
    std::string fileName = vm["output"].as<std::string>();
    IOHelper::export2OFF(oriMesh, fileName);
    return 0;
}
