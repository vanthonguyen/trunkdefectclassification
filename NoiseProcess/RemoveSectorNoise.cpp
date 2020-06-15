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

#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

#include "../Common/IOHelper.h"
#include "../Segmentation/SegmentationHelper.h"


using namespace DGtal;
namespace po = boost::program_options;

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

    std::vector<bool> notNoise(points.size(), true);

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

    double minZ = points[0][2];
    double maxZ = points[0][2];

    for( Z3i::RealPoint p : points ){
        if(p[2] < minZ) {
            minZ = p[2];
        }
        if(p[2] > maxZ){
            maxZ = p[2];
        }
    }
    int nbSlice = (maxZ - minZ) / length + 1;

    std::vector< std::vector<unsigned int> > slices(nbSlice);
    //std::vector<Z3i::RealPoint> centroids(nbSlice, Z3i::RealPoint(0,0,0));
    //std::vector<double> meanR(nbSlice, 0);

    for( int i = 0; i < points.size(); i++ ){
        int sliceId = (points[i][2] - minZ) / length;
        slices[sliceId].push_back(i);
        //centroids[sliceId] += points[i];
    }

   

    Z3i::RealPoint oz(0,0,1);
   
    double dummyR = 10000;
    std::vector<unsigned int> filtered;
    std::vector<bool> onTrunks(points.size(), false);
    for(int i = 0; i < nbSlice; i++){
        std::vector<unsigned int> slice = slices[i];
        if( slices[i].size() < 50 ){
            std::cout<< "slice with <  "<< 50<<" point" <<std::endl;

            for(int pind = 1; pind < slice.size(); pind++){
                unsigned int pId = slice[pind];
                filtered.push_back(pId);
                onTrunks[pId] = true;
            }
            continue;
        }
        std::cout<<"processing slide: "<<  i<< "/"<< nbSlice <<"with "<< slice.size()<<" points "<< std::endl;
        CylinderModel cylModel = SegmentationHelper::fitCylinderByCircle2D(points, slice, radius);
        std::cout<< cylModel.radius<<std::endl;
        if(cylModel.radius <= 0){
            std::cout<< "Invalide radius:  "<< cylModel.radius <<std::endl;
            for(int pind = 1; pind < slice.size(); pind++){
                unsigned int pId = slice[pind];
                filtered.push_back(pId);
                onTrunks[pId] = true;
            }
            continue;
        }
        //Z3i::RealPoint centroid = centroids[i];
        //double radius = meanR[i];
        double angleStep = length/cylModel.radius;
        //first point as marqueur
        Z3i::RealPoint ma = SegmentationHelper::getRadialVector(cylModel, points[slice[0]]);
        int nbSector = 2*M_PI / angleStep + 1;
        std::vector<unsigned int> mIndex(nbSector);
        std::vector<double> minR(nbSector, dummyR);

        for(int pind = 1; pind < slice.size(); pind++){
            unsigned int pId = slice[pind];
           // Z3i::RealPoint vectRadial = points[pId] - centroid;
            Z3i::RealPoint vectRadial = SegmentationHelper::getRadialVector(cylModel, points[pId]);

            //angle bt this point and ma

            double sca = vectRadial.dot(ma)/ ma.norm()/vectRadial.norm();
            //avoid sca > 1 due to floating point arithmatic?
            if(sca > 1){
                sca = 1;
            }
            double angle = acos(sca);
            //Z3i::RealPoint crossProduct = vectMarks[segmentId].crossProduct(vectRadial);
            Z3i::RealPoint u = ma.crossProduct(oz);
            if (u.dot(vectRadial) < 0){
                angle = 2 * M_PI - angle;
            }
            
            //current sector
            int sectId = angle / angleStep;
            double curR = vectRadial.norm();
            if(minR[sectId] > curR){
                minR[sectId] = curR;
                mIndex[sectId] = pId;
            }
        }
        for(int mInd = 0; mInd < mIndex.size(); mInd++){
            if(minR[mInd] < dummyR){
                filtered.push_back(mIndex[mInd]);
                onTrunks[mIndex[mInd]] = true;
            }
        }
    }

    //extending
   
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud (cloud);

    std::vector<int> pointIdx;
    std::vector<float> pointRadiusSquaredDistance;
    std::vector<int> extendPoints;
    //length * sqrt(2)
    double searchRadius = length*1.5;
    for (int pId = 0; pId < filtered.size(); pId++){
        pcl::PointXYZ searchPoint = cloud->points[filtered[pId]];
        if ( kdtree.radiusSearch (searchPoint, searchRadius, pointIdx, pointRadiusSquaredDistance) > 0 && pointIdx.size() > 0) {
            for(int foundId: pointIdx){
                if(!onTrunks[foundId]){
                    //branchIdsExtend.push_back(foundId);
                    onTrunks[foundId] = true;
                }
            }
        }
        //std::cout << "extending processing point " << pId <<"/" << branchIds.size()<<std::endl;
    }
    //output
    std::vector<unsigned int> trunkIds;
    std::vector<unsigned int> noiseIds;
    for (int i = 0; i < points.size(); i++){
        if(onTrunks[i]){
            trunkIds.push_back(i);
        }else{
            noiseIds.push_back(i);
        }
    }
    std::string fileName = vm["output"].as<std::string>();
    IOHelper::export2Text(points, trunkIds, "trunk-"+fileName);
    IOHelper::export2Text(points, noiseIds, "noises-"+fileName);

    return 0;
}
