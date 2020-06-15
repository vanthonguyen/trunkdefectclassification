#include <iostream>
#include <fstream>
#include <utility>


#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>


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

#include <DGtal/io/writers/VolWriter.h>                                                                                 
#include <DGtal/io/writers/LongvolWriter.h>  


#include "../IOHelperPcl.h"

#include "../Common/IOHelper.h"
#include "../Common/Rosin.h"


using namespace DGtal;
namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

double computeAreaThreshold(const Mesh<Z3i::RealPoint> &mesh){
    std::vector<double> faceAreas;
    for (unsigned int i = 0; i < mesh.nbFaces(); i++){
        //trace.progressBar(i , mesh.nbFaces());
        Face aFace = mesh.getFace(i);
        if(aFace.size() != 3){
            //facesToBeRemoved.push_back(i);
            continue;
        }
        Z3i::RealPoint p0 = mesh.getVertex(aFace.at(0));
        Z3i::RealPoint p1 = mesh.getVertex(aFace.at(1));
        Z3i::RealPoint p2 = mesh.getVertex(aFace.at(2));
        
        double l1 = (p1 - p0).norm();
        double l2 = (p2 - p0).norm();
        double l3 = (p2 - p1).norm();
        //half of perimeter
        double s = (l1 + l2 + l3)/2;
        //Heron's formula
        double a = sqrt(s * ( s - l1 ) * (s - l2) * (s - l3));
        faceAreas.push_back(a);
    }
    return Rosin::compute(faceAreas, 0.1);
}

int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input mesh.")
        ("area,a", po::value<double>(), "max area of a face.")
        ("vertexLength,l", po::value<double>()->default_value(30), "max vertex length.")
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
        trace.info()<< "Remove faces: " <<std::endl << "Options: "<<std::endl
            << general_opt << "\n";
        return 0;
    }

    DGtal::Mesh<Z3i::RealPoint> oriMesh(true);
    std::string inputMeshName = vm["input"].as<std::string>();
    MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, oriMesh, false);

    double maxArea = 0;
    if(vm.count("area")){
        maxArea = vm["area"].as<double>();
    }else{
        maxArea = computeAreaThreshold(oriMesh);
    }
    double maxVertexLength = vm["vertexLength"].as<double>();

std::cout<<"max A: "<<maxArea<<std::endl;
    std::vector<unsigned int> facesToBeRemoved;
    for (unsigned int i = 0; i < oriMesh.nbFaces(); i++){
        trace.progressBar(i , oriMesh.nbFaces());
        Face aFace = oriMesh.getFace(i);
        if(aFace.size() != 3){
            facesToBeRemoved.push_back(i);
            continue;
        }
        Z3i::RealPoint p0 = oriMesh.getVertex(aFace.at(0));
        Z3i::RealPoint p1 = oriMesh.getVertex(aFace.at(1));
        Z3i::RealPoint p2 = oriMesh.getVertex(aFace.at(2));
        
        double l1 = (p1 - p0).norm();
        double l2 = (p2 - p0).norm();
        double l3 = (p2 - p1).norm();
        //half of perimeter
        double s = (l1 + l2 + l3)/2;
        //Heron's formula
        double a = sqrt(s * ( s - l1 ) * (s - l2) * (s - l3));
        if(l1 > maxVertexLength || l2 > maxVertexLength || l3 > maxVertexLength || a > maxArea){
            facesToBeRemoved.push_back(i);
        }
    }

    oriMesh.removeFaces(facesToBeRemoved);
    std::cout<< "removed "<< facesToBeRemoved.size()<<std::endl;
    std::string fileName = vm["output"].as<std::string>();
    IOHelper::export2OFF(oriMesh, fileName);
    return 0;
}
