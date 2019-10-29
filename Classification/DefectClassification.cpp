#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>
#include <stdio.h>


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

using namespace DGtal;

namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;


#define RESOLUTION_FOR_FEATURE 3

int
main(int argc,char **argv){
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("verbose,v", "write more meshes.")
        ("input,i", po::value<std::string>(), "input mesh.")
        ("accRadius,r", po::value<double>()->default_value(350), "accumulation radius.")
        ("searchRadius,R", po::value<double>()->default_value(30), "search for neighbor in radius.")
        ("binWidth,b", po::value<double>()->default_value(0.001), "bin width used to compute threshold")
        ("patchWidth,a", po::value<double>()->default_value(25), "Arc length/ width of patch")
        ("patchHeight,e", po::value<int>()->default_value(100), "Height of patch")
        ("nbControl", po::value<int>()->default_value(2), "Number of control points for bsplines")
        ("confThreshold,t", po::value<double>()->default_value(0.0), "threshold in the confidence estimation.")
        ("voxelSize", po::value<int>()->default_value(3), "Voxel size")
        ("minClusterSize", po::value<int>()->default_value(100), "nbPoints to be considered as cluster")
        ("sectorLength", po::value<double>()->default_value(50), "used to segment branch")
        ("clusterTolerance", po::value<int>()->default_value(20), "minimum distance between two clusters")
        ("centerline", po::value<std::string>(), "Centerline of log or trunk")
        ("seedPoint", po::value<std::string>(), "SeedPoint to compute angle")
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

    //bizzard
    bool isVerbose = true;
    if ( vm.count("verbose") ){
        isVerbose = true;
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

    //MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshNam/mnt/data/dataPointCloud/Defauts2017/ScanWithMesure/enMm/$trunkName/centerline.xyz e, scaledMesh, false);
    MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, oriMesh, false);
    std::vector<Z3i::RealPoint> centerline;
    std::string outputPrefix = vm["output"].as<std::string>();

    int minClusterSize = vm["minClusterSize"].as<int>();
    int clusterTolerance = vm["clusterTolerance"].as<int>();

    int voxelSize = vm["voxelSize"].as<int>();
    double confThreshold = vm["confThreshold"].as<double>();
std::cout<<"voxelSize:"<<voxelSize<<std::endl;
    assert(voxelSize > 0);

    double searchRadius = vm["searchRadius"].as<double>();
    double accRadius = vm["accRadius"].as<double>();
    int nbControlPoint = vm["nbControl"].as<int>();
    int nbSegment = vm["nbSegment"].as<int>();

    //default -> oy
    Z3i::RealPoint seedPoint(0,1,0);
    if ( vm.count("seedPoint") ){
        std::vector<Z3i::RealPoint> seeds = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["seedPoint"].as<std::string>());
        assert(seeds.size() > 0);
        seedPoint = seeds[0];
    }

    /**step 1: build point cloud from face centers**/
    std::vector<Z3i::RealPoint> pointCloud(oriMesh.nbFaces());
    for (int i = 0; i < oriMesh.nbFaces(); i++){
        pointCloud[i] = oriMesh.getFaceBarycenter(i);
        //set color
        oriMesh.setFaceColor(i, DEFECT_COLOR[0]);
    }



    std::vector<Z3i::RealPoint> sampledCloud;

    time_t begin,end; // time_t is a datatype to store time values.
    //radius???

    time (&begin); 
    SegmentationHelper::simpleSubSample(pointCloud, voxelSize, sampledCloud);
    time (&end); // note time after execution
    double difference = difftime (end,begin);
    printf ("time taken for subsample %.2lf seconds.\n", difference );
    //centerline
    if ( vm.count("centerline") ){
        std::cout<<"reading centerline"<<std::endl;
        centerline = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["centerline"].as<std::string>());
    }else{
        std::cout<<"compute centerline"<<std::endl;

        int nbControlPoint = vm["nbControl"].as<int>();

        double searchRadius = vm["searchRadius"].as<double>();
        double accRadius = vm["accRadius"].as<double>();
        time (&begin); 
        centerline = CenterlineHelper::computeCenterline(sampledCloud, voxelSize, searchRadius, accRadius, nbControlPoint, nbSegment, confThreshold);

        time (&end); // note time after execution
        difference = difftime (end,begin);
        printf ("time taken for centerline computation %.2lf seconds.\n", difference );
        IOHelper::export2Text(centerline, "centerline.xyz");
    }
    if(centerline.size() < 1){
        std::cerr<<"Centerline with 0 point !!"<<std::endl;
        return 0;
    }


    time (&begin); // note time after execution
    DefectBranchSegmentation defectSegmentation( pointCloud, centerline, patchWidth, patchHeight, binWidth, sectorLength, minClusterSize, clusterTolerance, voxelSize);

    time (&end); // note time after execution
    difference = difftime (end,begin);
    printf ("time taken for branch segmentation %.2lf seconds.\n", difference );
    std::vector< std::vector<unsigned int> > defectClusters;

    time (&begin); // note time after execution
    defectSegmentation.segment(defectClusters);
    time (&end); // note time after execution
    difference = difftime (end,begin);
    printf ("time taken for segmentation %.2lf seconds.\n", difference );

    //compute radius at each trunk segment defined by centerline
    
    CylindricalCoordinateSystem ccs(centerline, seedPoint);
    
    std::vector<std::vector<double> > radiiBySegment(centerline.size() - 1);
    std::vector<double> segmentRadiiMode(centerline.size() - 1);
    std::vector<double> segmentRadiiMin(centerline.size() - 1, std::numeric_limits<double>::max());
    std::vector<double> segmentRadiiMax(centerline.size() - 1, -std::numeric_limits<double>::max());
    std::vector<double> trunkRadiis;
    double trunkRadMin;
    double trunkRadMax;
    for(Z3i::RealPoint xyzPoint : pointCloud){

        CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);

        radiiBySegment[cylPoint.segmentId].push_back(cylPoint.radius);
        trunkRadiis.push_back(cylPoint.radius);

        int segmentId = cylPoint.segmentId;
        //std::cout<<segmentId;
        if(cylPoint.radius < segmentRadiiMin[segmentId]){
            segmentRadiiMin[segmentId] = cylPoint.radius;
        }

        if(cylPoint.radius > segmentRadiiMax[segmentId]){
            segmentRadiiMax[segmentId] = cylPoint.radius;
        }

        if(cylPoint.radius < trunkRadMin){
            trunkRadMin = cylPoint.radius;
        }

        if(cylPoint.radius > trunkRadMax){
            trunkRadMax = cylPoint.radius;
        }
    }

    double trunkRadiusEst = Statistic::getMode(trunkRadiis, trunkRadMin, trunkRadMax, 1);

    for(unsigned int i = 0; i < centerline.size() - 1; i++){
        std::vector<double> radiisWithVoisins = radiiBySegment[i];

        double radiiMin = segmentRadiiMin[i];
        double radiiMax = segmentRadiiMax[i];
        for (int vid = 1; vid < 5; vid++){
            int index1 = i - vid;
            int index2 = i + vid;
            if( index1 >= 0 ){
                radiisWithVoisins.insert( radiisWithVoisins.end(), radiiBySegment[index1].begin(), radiiBySegment[index1].end() );
                if(segmentRadiiMin[index1] < radiiMin ){
                    radiiMin = segmentRadiiMin[index1];
                }

                if(segmentRadiiMax[index1] > radiiMax ){
                    radiiMax = segmentRadiiMax[index1];
                }
            }

            if( index2 < centerline.size() - 1 ){
                radiisWithVoisins.insert( radiisWithVoisins.end(), radiiBySegment[index2].begin(), radiiBySegment[index2].end() );

                if(segmentRadiiMin[index2] < radiiMin ){
                    radiiMin = segmentRadiiMin[index2];
                }

                if(segmentRadiiMax[index2] > radiiMax ){
                    radiiMax = segmentRadiiMax[index2];
                }
            }

        }
        if(radiisWithVoisins.size() == 0){
            segmentRadiiMode[i] = 0;
        }else{
            segmentRadiiMode[i]= Statistic::getMode(radiisWithVoisins, radiiMin, radiiMax, 1);
        }

        if ( segmentRadiiMode[i] > 1.5*trunkRadiusEst || segmentRadiiMode[i] < 0.7*trunkRadiusEst){
            segmentRadiiMode[i] = trunkRadiusEst;
        }
        std::cout<<i <<"/"<<centerline.size() - 1 << " segment points:" << radiisWithVoisins.size()<<std::endl;
    }
   

    //IOHelper::export2Text(segmentRadiiMode, "segmentradii");
    int defId = 0;
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

    time (&begin); // note time after execution
    for(int cId = 0; cId < defectClusters.size() - 1; cId++ ){
    //for(std::vector<unsigned int> defCluster: defectClusters){
        std::vector<unsigned int> defCluster = defectClusters[cId];
        std::vector<Z3i::RealPoint> defPointCloud(defCluster.size());
        for(int i = 0; i < defCluster.size(); i++){
            unsigned int pid = defCluster[i];
            defPointCloud[i] = pointCloud[pid];
        }

        std::vector<Z3i::RealPoint> sampledDefCloud;
        SegmentationHelper::simpleSubSample(defPointCloud, RESOLUTION_FOR_FEATURE , sampledDefCloud);
        DefectQuantification defquan(sampledDefCloud, centerline, seedPoint, segmentRadiiMode, species, rfFile);
        Defect def = defquan.quantify();
        //std::cout<< def.toString()<<std::endl;
        //write stat file
        defStatsAllClassStream << cId << " " << def.type << " " << def.height    << " " <<def.width << " "<< def.coordinateL <<" "<<def.coordinateZ<< " "<< def.branchDiameter<< std::endl;
        //defIdAllClass++;
        if( def.type <= SMALLDEFECT ){
            defStatsStream << defId << " " << def.type << " " << def.height    << " " <<def.width << " "<< def.coordinateL <<" "<<def.coordinateZ<< " "<< def.branchDiameter<<std::endl;
            defId++;
            for(int i = 0; i < defCluster.size(); i++){
                unsigned int pid = defCluster[i];
                rsMesh.setFaceColor(pid, DEFECT_COLOR[def.type]);
                defUniColorMesh.setFaceColor(pid, defUniCol);
                classified.push_back(pid);
            }
        }
        IOHelper::export2Text(defPointCloud, "def" + std::to_string(cId) + ".xyz");
        for(int i = 0; i < defCluster.size(); i++){
            unsigned int pid = defCluster[i];
            allDefMesh.setFaceColor(pid, defUniCol);
        }

    }

    time(&end);
    difference = difftime (end,begin);
    printf ("time taken for classification %.2lf seconds.\n", difference );

    defStatsStream.close();
    defStatsAllClassStream.close();

    
    if(isVerbose){
        std::string allDefectIdFile = outputPrefix + "-alldef.id";
        IOHelper::export2Text(defectClusters[defectClusters.size() - 1], allDefectIdFile);

        std::string classifiedDefFile = outputPrefix + "-classifed.id";
        IOHelper::export2Text(classified, classifiedDefFile);
    }
    //write color file
    std::string defectFile = outputPrefix + "-defect.off";
    IOHelper::export2OFF(rsMesh, defectFile);

    //unicolor
    if(isVerbose){
        std::string defectUnicolFile = outputPrefix + "-defectunicol.off";
        IOHelper::export2OFF(defUniColorMesh, defectUnicolFile);

        std::string defectFileNoClass = outputPrefix + "-defectnoclass.off";
        IOHelper::export2OFF(allDefMesh, defectFileNoClass);
    }


    //write error map file
    
    if(isVerbose){
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
        
        std::vector<unsigned int> branchIds = defectSegmentation.getBranchIds();

        errorMesh.removeFaces(branchIds);
        IOHelper::export2OFF(errorMesh, "trunk-error.off");
    }

    return 0;
}
