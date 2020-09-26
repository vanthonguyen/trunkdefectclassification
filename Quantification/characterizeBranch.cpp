#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <limits>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/kernel/sets/DigitalSetFromMap.h>
#include <DGtal/kernel/sets/DigitalSetBySTLSet.h>
#include "DGtal/geometry/volumes/distance/FMM.h"
#include "DGtal/shapes/Mesh.h"

#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

//#include "SimpleNormalAccumulator.h"
//
//

#include "../Common/IOHelper.h"
#include "../Common/Statistic.h"
#include "../Common/CylindricalCoordinateSystem.h"
#include "../Segmentation/SegmentationHelper.h"
#include "../Segmentation/EuclideCluster.h"
#include "../Centerline/CenterlineMultiCore.h"
#include "../Centerline/CenterlineHelper.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;



typedef DigitalSetBySTLSet<Z3i::Domain> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;



using namespace DGtal::Z3i;
using namespace DGtal;
typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;


std::pair<Z3i::RealPoint, Z3i::RealPoint> 
getExtremityOfCylinder(const std::vector<Z3i::RealPoint> &pointCloud, 
                       Z3i::RealPoint center1, Z3i::RealPoint center2){
    Z3i::RealPoint vectDir = (center2 - center1).getNormalized();
    double minScalar = std::numeric_limits<double>::max();
    double maxScalar = -std::numeric_limits<double>::max();
    Z3i::RealPoint minPoint(0,0,0);
    Z3i::RealPoint maxPoint(0,0,0);
    for(Z3i::RealPoint p : pointCloud){

        Z3i::RealPoint aVect = p - center1;
        double scalar = aVect.dot(vectDir);
        if(scalar < minScalar ){
            minScalar = scalar;
            minPoint = center1 + scalar*vectDir;
        }
        if(scalar > maxScalar){
            maxScalar = scalar;
            maxPoint = center1 + scalar*vectDir;
        }
    }
    std::pair<Z3i::RealPoint, Z3i::RealPoint> extre(minPoint, maxPoint);
    return  extre;
}

/*
 * If p1p2 and p3p4 are screw, the volume of tetrahedron is > 0
 */
bool testSkewness(Z3i::RealPoint p1, Z3i::RealPoint p2, Z3::RealPoint p3, Z3i::RealPoint p4){

//trace.info()<<p1<<std::endl;
//trace.info()<<p2<<std::endl;
//trace.info()<<p3<<std::endl;
//trace.info()<<p4<<std::endl;
    //V = 1/6 * (p1 - p4). ((p2 - p4)x(p3-p4))
    Z3i::RealPoint p14 = p1 - p4;
    Z3i::RealPoint p24 = p2 - p4;
    Z3i::RealPoint p34 = p3 - p4;
    //volume 
    double V = 1.0/6 * std::abs(p14.dot(p24.crossProduct(p34)));
//    std::cout<<p14<<std::endl;
//    std::cout<<p24<<std::endl;
//    std::cout<<p34<<std::endl;
//    std::cout<<p24.crossProduct(p34)<<std::endl;
//    std::cout<<"V"<<V<<std::endl;
    return V > 0; 
}

std::pair<Z3i::RealPoint, Z3i::RealPoint> 
getNearestPoints(Z3i::RealPoint p1, Z3i::RealPoint d1, Z3::RealPoint p2, Z3i::RealPoint d2){
    Z3i::RealPoint n = d1.crossProduct(d2);
    //The plane formed by the translations of Line 2 along n contains the point p2 and is perpendicular to n2 = d2 × n
    Z3i::RealPoint n2 = d2.crossProduct(n);
    Z3i::RealPoint n1 = d1.crossProduct(n);
    //nearest 1 = (p1 + (p2 - p1).n2 / d1.n2) d1 
    Z3i::RealPoint nearest1 = p1 + (( p2 - p1 ).dot(n2)/(d1.dot(n2))) * d1;
    Z3i::RealPoint nearest2 = p2 + (( p1 - p2 ).dot(n1)/(d2.dot(n1))) * d2;
std::cout<<"n1#"<<nearest1<<std::endl;
std::cout<<"p2#"<<p2<<std::endl;

    std::pair<Z3i::RealPoint, Z3i::RealPoint> nearests(nearest1, nearest2);
    return nearests;
}

bool
getIntersectionSegment(const std::vector<Z3i::RealPoint> &lines, Z3i::RealPoint p1, Z3i::RealPoint d1, int &segmentId, Z3i::RealPoint &nearest){
    for(int i = 0; i < lines.size() -1; i++){
std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
        //test for skewness
        Z3i::RealPoint p12 = p1 + d1;
        Z3i::RealPoint p2 = lines[i];
        Z3i::RealPoint p22 = lines[i + 1];
        bool isSkew = testSkewness(p1, p12, p2, p22);
        Z3i::RealPoint d2 = (p22 - p2).getNormalized();
        Z3i::RealPoint nearestPoint(std::numeric_limits<double>::max(), 0, 0);
        if( isSkew ){
            //nearest point on line segment i
            std::pair<Z3i::RealPoint, Z3i::RealPoint> nearests = getNearestPoints(p2, d2, p1, d1);
            nearestPoint = nearests.first;
std::cout<<"@nr"<<nearestPoint<<std::endl;
        }else{
            // c = p1 + td1 = p2 + td2
            //d1[0]*t - d2[0]*v =  -p1[0] + p2[0; 
            //d1[1] *t - d2[1]*v = - p1[1] + p2[1];
            
            double a = d1[0];
            double b = -d2[0];
            double c = p1[0] - p2[0];

            double d = d1[1];
            double e = -d2[1];
            double f = p1[1] - p2[1];

            double det = a* d - b * e;
            if(det == 0){
                trace.info()<<"getIntersectionSegment : det = 0 !!!!";
                continue;
            }
            double t = (e *d - b*f)/det;
            double v = (a * f - c * e) / det;
            nearestPoint = p1 + t * d1;
//trace.info()<<"intersection"<<std::endl;
        }
        //double scalar = (p2 - nearestPoint).dot(segmentPoint2 - nearestPoint);
        //std::cout<<"scalar "<< scalar<<std::endl;
        double segmentLength = (p2-p22).norm();
        double nearestP2 = (nearestPoint - p2).norm();
        double nearestSegmentPoint2 = (nearestPoint - p22).norm();
        std::cout<<"nr"<< nearestPoint<<std::endl;
        std::cout<<"ps2"<< p22<<std::endl;
        std::cout<<"s length:" <<segmentLength<<std::endl;
        std::cout<<"p2 :" <<nearestP2<<std::endl;
        std::cout<<"p3:" <<nearestSegmentPoint2<<std::endl;
        if(nearestP2 < segmentLength && nearestSegmentPoint2 < segmentLength){
            segmentId = i;
            nearest = nearestPoint;
trace.info()<<nearest<<"xxx"<<std::endl;
            return true;
        }
        std::cout<<"###################"<<std::endl;
    }
    return false;
}

Z3i::RealPoint getLogVectorAtBranchPosition2(const std::vector<Z3i::RealPoint> &lines, const int segmentId, const double length){
    Z3i::RealPoint upperSum(0,0,0);
    Z3i::RealPoint lowerSum(0,0,0);
    int upIndex = segmentId;
    double currentLength = 0.0;
    while(upIndex < lines.size() && currentLength < length / 2){
        upperSum += lines.at(upIndex);
        currentLength = (lines[upIndex] - lines[segmentId]).norm();
        upIndex++;
    }
    int nbUpperSegment = upIndex - segmentId;

    int downIndex = segmentId - 1;
    currentLength = 0.0;

    while(downIndex > 0 && currentLength < length / 2){
        lowerSum += lines.at(downIndex);
        currentLength = (lines[downIndex] - lines[segmentId - 1]).norm();
        downIndex--;
    }
    int nbLowerSegment = segmentId - 1 - downIndex;

    return upperSum / nbUpperSegment - lowerSum / nbLowerSegment;
}

std::pair<Z3i::RealPoint, Z3i::RealPoint> 
getLogVectorAtBranchPosition(const std::vector<Z3i::RealPoint> &lines, const int segmentId, const double length){
    Z3i::RealPoint upperSum(0,0,0);
    Z3i::RealPoint lowerSum(0,0,0);
    int upIndex = segmentId;
    double currentLength = 0.0;
    while(upIndex < lines.size() && currentLength < length / 2){
        upperSum += lines.at(upIndex);
        currentLength = (lines[upIndex] - lines[segmentId]).norm();
        upIndex++;
    }

    int nbUpperSegment = upIndex - segmentId;

    int downIndex = segmentId - 1;
    currentLength = 0.0;

    while(downIndex > 0 && currentLength < length / 2){
        lowerSum += lines.at(downIndex);
        currentLength = (lines[downIndex] - lines[segmentId - 1]).norm();
        downIndex--;
    }
    int nbLowerSegment = segmentId - 1 - downIndex;
    std::pair<Z3i::RealPoint, Z3i::RealPoint> vects(lowerSum/ nbLowerSegment, upperSum/ nbUpperSegment);
    return vects;

}


#define MAX_LENGTH_FOR_EST 200

/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv){
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input defect branch.")
        ("accRadius,r", po::value<double>()->default_value(50), "accumulation radius.")
        ("searchRadius,R", po::value<double>()->default_value(20), "search for neighbor in radius.")
        ("centerline,l", po::value<std::string>(), "input trunk centerline branch.")
        ("maxLength,m", po::value<double>()->default_value(200), "skip")
        ("skip,s", po::value<double>()->default_value(30), "skip")
        ("output,o", po::value<std::string>()->default_value("branch"), "branch, ...");


    bool parseOK = true;
    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }
    catch (const std::exception &ex)
    {
        trace.info() << "Error checking program options: " << ex.what() << std::endl;
        parseOK = false;
    }
    po::notify(vm);
    if ( !parseOK || vm.count("help") || argc <= 1 || !vm.count("input") )
    {
        trace.info() << "Characterize branches of standing tree" << std::endl
            << "Options: " << std::endl
            << general_opt << std::endl;
        return 0;
    }


    // Reading parameters:
    //double distanceSearch = vm["distanceSearch"].as<double>();
    double skip = vm["skip"].as<double>();
    double maxLength = vm["maxLength"].as<double>();

    std::vector<Z3i::RealPoint> pointCloud = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["input"].as<std::string>());
    std::vector<Z3i::RealPoint> trunkCenterline = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["centerline"].as<std::string>());
std::cout<<"nb points"<<pointCloud.size()<<std::endl;

    double searchRadius = vm["searchRadius"].as<double>();
    double accRadius = vm["accRadius"].as<double>();
    std::vector<Z3i::RealPoint> branches;
    Z3i::RealPoint seedPoint(0,1,0);
    CylindricalCoordinateSystem ccs(trunkCenterline, seedPoint);
    double minR = std::numeric_limits<double>::max();
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint xyzPoint = pointCloud.at(i);
        CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);
        if(cylPoint.radius < minR){
            minR = cylPoint.radius;
        }
    }
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint xyzPoint = pointCloud.at(i);
        CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);
        if(cylPoint.radius > minR + skip && cylPoint.radius < minR + skip + MAX_LENGTH_FOR_EST){
            branches.push_back(xyzPoint);
        }
    }
std::cout<<"nb points in branch"<<branches.size()<<std::endl;
    std::vector<Z3i::RealPoint> centerline = CenterlineHelper::computeCenterline(branches, 1, searchRadius, accRadius , 2, 1, 0);
    
    CylindricalCoordinateSystem branchCcs(centerline, seedPoint);
    std::vector<double> radiis;

    double radiiMin = std::numeric_limits<double>::max();
    double radiiMax = -std::numeric_limits<double>::max();
    std::vector<double> segmentRadiiMin(centerline.size() - 1, std::numeric_limits<double>::max());
    std::vector<double> segmentRadiiMax(centerline.size() - 1, -std::numeric_limits<double>::max());
    std::vector<std::vector<double> > radiiBySegment(centerline.size() - 1);
    for(unsigned int i = 0; i < branches.size(); i++){
        Z3i::RealPoint xyzPoint = branches.at(i);
        CylindricalPoint cylPoint = branchCcs.xyz2Cylindrical(xyzPoint);
        CylindricalPoint pointInTrunkCcs = ccs.xyz2Cylindrical(xyzPoint);
        if(pointInTrunkCcs.radius < minR + skip + maxLength){
            radiis.push_back(cylPoint.radius);
            if(cylPoint.radius < radiiMin){
                radiiMin = cylPoint.radius;
            }
            if(cylPoint.radius > radiiMax){
                radiiMax = cylPoint.radius;
            }
        }
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
    
std::cout<<"radiis nb:"<<radiis.size()<<std::endl;

    double branchRadiusEst = Statistic::getMode(radiis, radiiMin, radiiMax, 1);


    std::vector<double> segmentRadiiMode(centerline.size() - 1);

    for(unsigned int i = 0; i < centerline.size() - 1; i++){
        std::vector<double> radiisWithVoisins = radiiBySegment[i];

        double radiiMin = segmentRadiiMin[i];
        double radiiMax = segmentRadiiMax[i];
        for (int vid = 1; vid < 2; vid++){
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
            std::cout<<"Est"<< segmentRadiiMode[i]<<std::endl;
        }

        if ( segmentRadiiMode[i] > 1.5*branchRadiusEst || segmentRadiiMode[i] < 0.7*branchRadiusEst){
            //segmentRadiiMode[i] = branchRadiusEst;
            std::cout<<"radii is too small or large"<<std::endl;
        }
        std::cout<<segmentRadiiMode[i]<<std::endl;
    }


    std::cout<<"Centerline method"<<std::endl;
    std::cout<< "estimated radius: "<<branchRadiusEst<<std::endl;

    std::string outputPrefix = vm["output"].as<std::string>();
    IOHelper::export2Text(branches, outputPrefix + ".xyz");

    CylinderModel cylModel = SegmentationHelper::fitCylinder(branches, accRadius);

    Z3i::RealPoint center1 = cylModel.point;
    Z3i::RealPoint center2 = center1 + cylModel.direction;
    Z3i::RealPoint dir = cylModel.direction;

    std::pair<Z3i::RealPoint, Z3i::RealPoint> exs = getExtremityOfCylinder(branches, center1, center2);

    std::vector<Z3i::RealPoint> cylCenterline;
    cylCenterline.push_back(exs.first);
    cylCenterline.push_back(exs.second);

    std::cout<<"Cylinder method"<<std::endl;
    std::cout<< "estimated radius: "<<cylModel.radius<<std::endl;


    Mesh<Z3i::RealPoint> transMesh(true);
    //Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, cylCenterline, cylModel.radius, 0.1, DGtal::Color::Red);
    std::vector<Z3i::RealPoint> segs;
    bool first = true;
    for(unsigned int i = 0; i < centerline.size() - 1;i++){
        if(segmentRadiiMode[i] > 0){
            if(first){
                segs.push_back(centerline[i]);
                first = false;
            }
            segs.push_back(centerline[i + 1]);
        }
    }
    Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, segs, branchRadiusEst, 0.1, DGtal::Color::Green);
    //Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, cylCenterline, cylModel.radius, 0.1, DGtal::Color::Red);
   /* for(unsigned int i = 0; i < centerline.size() - 1;i++){
        if(segmentRadiiMode[i] > 0){
            std::vector<Z3i::RealPoint> seg;
            seg.push_back(centerline[i]);
            seg.push_back(centerline[i + 1]);
            Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, seg, segmentRadiiMode[i], 0.1, DGtal::Color::Green);
        }
    }*/

    IOHelper::export2OFF(transMesh, outputPrefix + "-cen.off");

    Mesh<Z3i::RealPoint> transMesh2(true);
    Mesh<Z3i::RealPoint>::createTubularMesh(transMesh2, cylCenterline, cylModel.radius, 0.1, DGtal::Color::Red);
    //Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, centerline, branchRadiusEst, 0.1, DGtal::Color::Green);

    IOHelper::export2OFF(transMesh2, outputPrefix + "-cyl.off");

  return 0;
}
