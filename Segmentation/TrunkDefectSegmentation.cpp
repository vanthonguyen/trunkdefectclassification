#include <iostream>
#include <fstream>
#include <utility>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


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


#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/shapes/Mesh.h"

#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"


#include "../Common/IOHelper.h"
#include "../Common/Statistic.h"
#include "../Common/DefectType.h"
#include "../Common/Defect.h"
#include "../Common/CylindricalCoordinateSystem.h"

#include "../Centerline/CenterlineMultiCore.h"
#include "../Centerline/CenterlineHelper.h"

#include "../Segmentation/EuclideCluster.h"
#include "../Segmentation/SegmentationHelper.h"
#include "../Segmentation/DefectBranchSegmentation.h"

#include "../Feature/FeatureExtractionXyz.h"
#include "../Feature/PredictRF.h"

#include "DefectQuantification.h"

using namespace DGtal;

namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;



int
main(int argc,char **argv){
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input mesh.")
        ("accRadius,r", po::value<double>()->default_value(350), "accumulation radius.")
        ("searchRadius,R", po::value<double>()->default_value(30), "search for neighbor in radius.")
        ("binWidth,b", po::value<double>()->default_value(0.001), "bin width used to compute threshold")
        ("patchWidth,a", po::value<double>()->default_value(25), "Arc length/ width of patch")
        ("patchHeight,e", po::value<int>()->default_value(100), "Height of patch")
        ("nbControl", po::value<int>()->default_value(2), "Number of control points for bsplines")
        ("voxelSize", po::value<int>()->default_value(3), "Voxel size")
        ("minClusterSize", po::value<int>()->default_value(100), "nbPoints to be considered as cluster")
        ("sectorLength", po::value<double>()->default_value(50), "used to segment branch")
        ("clusterTolerance", po::value<int>()->default_value(20), "minimum distance between two clusters")
        ("centerline", po::value<std::string>(), "Centerline of log or trunk")
        //("seedPoint", po::value<std::string>()->default_value("seedPoint.xyz"), "SeedPoint to compute angle")
        ("species", po::value<int>()->default_value(0), "Species of defect, 1 for oak, 2: birch")
        ("rfFile", po::value<std::string>()->default_value("learnedForest"), "random forest serialized file")
        ("padding", po::value<double>()->default_value(1000), "distance of trunk to soil")
        ("nbSegment", po::value<int>()->default_value(10), "Number of segments to compute the centerline")
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
    if(vm.count("help") || argc<=1 || !parseOK || !vm.count("input") || !vm.count("binWidth") || 
            !vm.count("patchWidth") || !vm.count("patchHeight") || !vm.count("rfFile")){
        if(!vm.count("input")){
            trace.error()<<"the input mesh is required!"<<std::endl;
        }else if(!vm.count("binWidth")){
            trace.error()<<"binWidth is required!"<<std::endl;
        }else if(!vm.count("patchWidth")){
            trace.error()<<"patchWidth is required!"<<std::endl;
        }else if(!vm.count("patchHeight")){
            trace.error()<<"patchHeight is required!"<<std::endl;
        }else if(!vm.count("rfFile")){
            trace.error()<<"rfFile is required!"<<std::endl;
        }
        trace.info()<< "Segmentation log defects" <<std::endl << "Options: "<<std::endl
            << general_opt << "\n";
        return 0;
    }

    double binWidth = vm["binWidth"].as<double>();

    double patchWidth = vm["patchWidth"].as<double>(); 
    int patchHeight = vm["patchHeight"].as<int>(); 

    double sectorLength = vm["sectorLength"].as<double>(); 
    double padding = vm["padding"].as<double>(); 
    std::string rfFile = vm["rfFile"].as<std::string>();

    int species = vm["species"].as<int>();

    DGtal::Mesh<Z3i::RealPoint> oriMesh(true);
    std::string inputMeshName = vm["input"].as<std::string>();

    //MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, scaledMesh, false);
    MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, oriMesh, false);
    std::vector<Z3i::RealPoint> centerline;
    std::string outputPrefix = vm["output"].as<std::string>();

    int minClusterSize = vm["minClusterSize"].as<int>();
    int clusterTolerance = vm["clusterTolerance"].as<int>();

    int voxelSize = vm["voxelSize"].as<int>();
std::cout<<"voxelSize:"<<voxelSize<<std::endl;
    assert(voxelSize > 0);

    double searchRadius = vm["searchRadius"].as<double>();
    double accRadius = vm["accRadius"].as<double>();
    int nbControlPoint = vm["nbControl"].as<int>();
    int nbSegment = vm["nbSegment"].as<int>();

    //std::vector<Z3i::RealPoint> seeds = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["seedPoint"].as<std::string>());
    //assert(seeds.size() > 0);
    Z3i::RealPoint seedPoint(0,1,0);

    /**step 1: build point cloud from face centers**/
    std::vector<Z3i::RealPoint> pointCloud(oriMesh.nbFaces());
    for (int i = 0; i < oriMesh.nbFaces(); i++){
        pointCloud[i] = oriMesh.getFaceBarycenter(i);
        //set color
        oriMesh.setFaceColor(i, DEFECT_COLOR[0]);
    }



    std::vector<Z3i::RealPoint> sampledCloud;

    //radius???
    SegmentationHelper::simpleSubSample(pointCloud, voxelSize, sampledCloud);

    //centerline
    if ( vm.count("centerline") ){
        std::cout<<"reading centerline"<<std::endl;
        centerline = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["centerline"].as<std::string>());
    }else{
        std::cout<<"compute centerline"<<std::endl;

        int nbControlPoint = vm["nbControl"].as<int>();

        double searchRadius = vm["searchRadius"].as<double>();
        double accRadius = vm["accRadius"].as<double>();

        centerline = CenterlineHelper::computeCenterline(sampledCloud, voxelSize, searchRadius, accRadius, nbControlPoint, nbSegment);
        IOHelper::export2Text(centerline, "centerline.xyz");
    }
    if(centerline.size() < 1){
        std::cerr<<"Centerline with 0 point !!"<<std::endl;
        return 0;
    }


    DefectBranchSegmentation defectSegmentation( pointCloud, centerline, patchWidth, patchHeight, binWidth, sectorLength, minClusterSize, clusterTolerance, voxelSize);

    std::vector< std::vector<unsigned int> > defectClusters;
    defectSegmentation.segment(defectClusters);

    //compute radius at each trunk segment defined by centerline
    
    CylindricalCoordinateSystem ccs(centerline, seedPoint);
    
    std::vector<std::vector<double> > radiiBySegment(centerline.size() - 1);
    std::vector<double> segmentRadiiMode(centerline.size() - 1);
    std::vector<double> segmentRadiiMin(centerline.size() - 1, std::numeric_limits<double>::max());
    std::vector<double> segmentRadiiMax(centerline.size() - 1, -std::numeric_limits<double>::max());
    for(Z3i::RealPoint xyzPoint : pointCloud){

        CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);

        radiiBySegment[cylPoint.segmentId].push_back(cylPoint.radius);

        int segmentId = cylPoint.segmentId;
        //std::cout<<segmentId;
        if(cylPoint.radius < segmentRadiiMin[segmentId]){
            segmentRadiiMin[segmentId] = cylPoint.radius;
        }

        if(cylPoint.radius > segmentRadiiMax[segmentId]){
            segmentRadiiMax[segmentId] = cylPoint.radius;
        }
    }

    for(int i = 0; i < centerline.size() - 1; i++){
        std::cout<<i <<"/"<<centerline.size() - 1 << " segment points:" << radiiBySegment[i].size()<<std::endl;
        if(radiiBySegment[i].size() == 0){
            continue;
        }
        segmentRadiiMode[i] = Statistic::getMode(radiiBySegment[i], segmentRadiiMin[i], segmentRadiiMax[i], 1);
    }
   

    //IOHelper::export2Text(segmentRadiiMode, "segmentradii");
    int defId = 0;
    int defIdAllClass = 0;
    std::ofstream defStatsStream;
    std::ofstream defStatsAllClassStream;
    std::string dfile = outputPrefix + "-defstats";
    std::string dfileAllClass = outputPrefix + "-defstatsallclass";
    defStatsStream.open(dfile.c_str(), std::ofstream::out);
    defStatsAllClassStream.open(dfileAllClass.c_str(), std::ofstream::out);


    DGtal::Mesh<Z3i::RealPoint> rsMesh = oriMesh; 
    Mesh<Z3i::RealPoint> allDefMesh = oriMesh;
    Mesh<Z3i::RealPoint> defUniColorMesh = oriMesh;

    DGtal::Color defUniCol(89, 166, 43);
    std::vector<unsigned int> classified;
    for(int cId = 0; cId < defectClusters.size() - 1; cId++ ){
    //for(std::vector<unsigned int> defCluster: defectClusters){
        std::vector<unsigned int> defCluster = defectClusters[cId];
        std::vector<Z3i::RealPoint> defPointCloud(defCluster.size());
        for(int i = 0; i < defCluster.size(); i++){
            unsigned int pid = defCluster[i];
            defPointCloud[i] = pointCloud[pid];
        }
        DefectQuantification defquan(defPointCloud, centerline, seedPoint, segmentRadiiMode, species, rfFile);
        Defect def = defquan.quantify();
        std::cout<< def.toString()<<std::endl;
        //write stat file
        defStatsAllClassStream << defIdAllClass << " " << def.type << " " << def.height    << " " <<def.width << " "<< def.coordinateL <<" "<<def.coordinateZ<<std::endl;
        defIdAllClass++;
        if( def.type <= SMALLDEFECT ){
            defStatsStream << defId << " " << def.type << " " << def.height    << " " <<def.width << " "<< def.coordinateL <<" "<<def.coordinateZ<<std::endl;
            defId++;
            for(int i = 0; i < defCluster.size(); i++){
                unsigned int pid = defCluster[i];
                rsMesh.setFaceColor(pid, DEFECT_COLOR[def.type]);
                defUniColorMesh.setFaceColor(pid, defUniCol);
                classified.push_back(cluster[pId]);
            }
        }
        for(int i = 0; i < defCluster.size(); i++){
            unsigned int pid = defCluster[i];
            allDefMesh.setFaceColor(pid, defUniCol);
        }

    }

    defStatsStream.close();
    defStatsAllClassStream.close();

    
    std::string allDefectIdFile = outputPrefix + "-alldef.id";
    IOHelper::export2Text(defectClusters[defectClusters.size() - 1], allDefectIdFile);
    
    std::string classifiedDefFile = outputPrefix + "-classifed.id";
    IOHelper::export2Text(classified, classifiedDefFile);

    //write color file
    std::string defectFile = outputPrefix + "-defect.off";
    IOHelper::export2OFF(rsMesh, defectFile);

    //unicolor
    std::string defectUnicolFile = outputPrefix + "-defectunicol.off";
    IOHelper::export2OFF(defUniColorMesh, defectUnicolFile);

    std::string defectFileNoClass = outputPrefix + "-defectnoclass.off";
    IOHelper::export2OFF(allDefMesh, defectFileNoClass);

    //write error map file
    
    double minValue = -10;
    double maxValue = 20.0;
    DGtal::GradientColorMap<double, CMAP_JET>  gradientShade(minValue, maxValue); 
    DGtal::Mesh<Z3i::RealPoint> errorMesh = oriMesh;

    std::vector<double> distances = defectSegmentation.getDistances();
    for (int i = 0; i < errorMesh.nbFaces(); i++){

        double err = distances.at(i);

        //centroid
        errorMesh.setFaceColor(i, gradientShade(err));
    }

    std::string errorFile = outputPrefix + "-error.off";
    IOHelper::export2OFF(errorMesh,errorFile);

    return 0;
}
