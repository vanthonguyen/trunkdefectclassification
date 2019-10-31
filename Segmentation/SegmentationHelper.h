#ifndef SEGMENTATION_HELPER_H 
#define SEGMENTATION_HELPER_H 

#include <set>
#include <map>

#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>

#include <pcl/filters/extract_indices.h>
#include <pcl/features/normal_3d.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_circle.h>

#include "DGtal/shapes/Mesh.h"

#include "../Common/IOHelper.h"
#include "../Common/Statistic.h"
#include "../Common/CylindricalCoordinateSystem.h"

#include "../Centerline/Centerline.h"
#include "../Centerline/CenterlineHelper.h"


#define NB_POINTS_MIN 30
#define BRANCH_SEG_SEARCH_COEFF 2


using namespace DGtal;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

typedef pcl::PointXYZ PointT;

struct CylinderModel{
    Z3i::RealPoint point;
    Z3i::RealPoint direction;
    double radius;
};

struct Circle2DModel{
    double x;
    double y;
    double radius;
};


class SegmentationHelper{
public:

    static void
    segmentBranchesUsingCenterline(const DGtal::Mesh<Z3i::RealPoint> &mesh, const int &voxelSize,
            const double &searchRadius, const double &accRadius, const int &nbSegment, const int &nbControlPoint,
            const double &threshold,
            const double &length,
            const std::vector<Z3i::RealPoint> &centerline,
            std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch){

        std::vector<Z3i::RealPoint> points(mesh.nbFaces());
        for (int i = 0; i < mesh.nbFaces(); i++){
            Z3i::RealPoint p = mesh.getFaceBarycenter(i);
            points[i] = p;
        }
        segmentBranchesUsingCenterline(points, voxelSize, searchRadius, accRadius, nbControlPoint, nbSegment, threshold,length, centerline, trunk, branch);
    }


    static void
    segmentBranchesUsingCenterline(const DGtal::Mesh<Z3i::RealPoint> &mesh, const int &voxelSize,
            const double &searchRadius, const double &accRadius, const int &nbSegment, const int &nbControlPoint,
            const double &threshold,
            const double &length,
            std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch){

        std::vector<Z3i::RealPoint> points(mesh.nbFaces());
        for (int i = 0; i < mesh.nbFaces(); i++){
            Z3i::RealPoint p = mesh.getFaceBarycenter(i);
            points[i] = p;
        }
        segmentBranchesUsingCenterline(points, voxelSize, searchRadius, accRadius, nbControlPoint, nbSegment, threshold,length, trunk, branch);
    }

    static void
    segmentBranchesUsingCenterline(const std::vector<Z3i::RealPoint> &points, const int &voxelSize,
            const double &searchRadius, const double &accRadius, const int &nbSegment, const int &nbControlPoint,
            const double &threshold,
            const double &length,
            std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch){
        assert(points.size() > 0);

        std::vector<Z3i::RealPoint> centerline = CenterlineHelper::computeCenterline(points, voxelSize, searchRadius, accRadius, nbControlPoint, nbSegment);

        segmentBranchesUsingCenterline(points, length, centerline, trunk, branch);
    }


    static void
    segmentBranchesUsingCenterline(const std::vector<Z3i::RealPoint> &points,
            const double &length,
            const std::vector<Z3i::RealPoint> &centerline,
            std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch){
        assert(points.size() > 0);

        Z3i::RealPoint p0 = points.at(0);
        Z3i::RealPoint ptLow = p0;
        Z3i::RealPoint ptUp = p0;


        for(unsigned int i = 0; i < points.size(); i++){
            Z3i::RealPoint pi = points.at(i);
            ptLow = ptLow.inf(pi);
            ptUp = ptUp.sup(pi);
        }

        double minZ = ptLow[2];
        double maxZ = ptUp[2];

        int nbSlice = (maxZ - minZ) / length + 1;

        std::vector< std::vector<unsigned int> > slices(nbSlice);

        CylindricalCoordinateSystem ccs(centerline, Z3i::RealPoint(0.0,0.0,0.0));
        //store the distance to centerline for each segment
        std::vector<std::vector<double> > radiiBySlice(nbSlice);
        std::vector<double> sliceRadiiMode(nbSlice);
        std::vector<double> sliceRadiiMin(nbSlice, std::numeric_limits<double>::max());
        std::vector<double> sliceRadiiMax(nbSlice, -std::numeric_limits<double>::max());
        std::vector<CylindricalPoint> cylPoints(points.size());

        std::vector<double> radiis(points.size());
        double minRadii = std::numeric_limits<double>::max();
        double maxRadii = -std::numeric_limits<double>::max();

        for(int i = 0; i < points.size();i++){
            Z3i::RealPoint xyzPoint = points.at(i);
            int sliceId = (points[i][2] - minZ) / length;
            slices[sliceId].push_back(i);

            CylindricalPoint cylPoint = ccs.xyz2Cylindrical(xyzPoint);

            radiiBySlice[sliceId].push_back(cylPoint.radius);

            //std::cout<<sliceId;
            if(cylPoint.radius < sliceRadiiMin[sliceId]){
                sliceRadiiMin[sliceId] = cylPoint.radius;
            }

            if(cylPoint.radius > sliceRadiiMax[sliceId]){
                sliceRadiiMax[sliceId] = cylPoint.radius;
            }
            cylPoints[i] = cylPoint;

            radiis[i] = cylPoint.radius;
            //global
            if(cylPoint.radius < minRadii){
                minRadii = cylPoint.radius;
            }

            if(cylPoint.radius > maxRadii){
                maxRadii = cylPoint.radius;
            }
        }

        double radiiMode = Statistic::getMode(radiis, minRadii, maxRadii, 1);
        for(int i = 0; i < slices.size() - 1; i++){
            std::cout<<"CenterlineHelper::segmentBranchesUsingCenterline: "<<i <<"/"<<nbSlice << " segment points:" << radiiBySlice[i].size()<<std::endl;
            if(radiiBySlice[i].size() == 0){
                continue;
            }
            //using a binwidth=5 for the histogram
            //sliceRadiiMode[i] = Statistic::getMode(radiiBySlice[i], sliceRadiiMin[i], sliceRadiiMax[i], 5);
            //global radii is better?
            sliceRadiiMode[i] = radiiMode;
        }

        segmentBranches(points, cylPoints, slices, sliceRadiiMode, length, centerline, trunk, branch);


       
        /*
        for(int i = 0; i < cylPoints.size(); i++){
            
            CylindricalPoint cylPoint = cylPoints.at(i);
            int segmentId = cylPoint.segmentId;
            if(cylPoint.radius > segmentRadiiMode[segmentId] + threshold){
                branch.push_back(i);
            }else{
                trunk.push_back(i);
            }
        }
        */

        std::cout<<"SegmentationHelper::segmentBranchesUsingCenterline: finish with"<< trunk.size() << " points on trunk and " << branch.size()<<" points on branches"<<std::endl;
 
    }
   

    static void
    segmentBranches(const std::vector<Z3i::RealPoint> &points, double length,
        std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch){

        std::cout<<"begin segment branches"<<std::endl;

        assert(points.size() > 0);

        Z3i::RealPoint p0 = points.at(0);
        Z3i::RealPoint ptLow = p0;
        Z3i::RealPoint ptUp = p0;


        for(unsigned int i = 0; i < points.size(); i++){
            Z3i::RealPoint pi = points.at(i);
            ptLow = ptLow.inf(pi);
            ptUp = ptUp.sup(pi);
        }


        double minZ = ptLow[2];
        double maxZ = ptUp[2];

        int nbSlice = (maxZ - minZ) / length + 1;
std::cout<<"maxZ: "<<maxZ<<std::endl;
std::cout<<"minZ: "<<minZ<<std::endl;
std::cout<<"nbSlice: "<<nbSlice<<std::endl;
        std::vector< std::vector<unsigned int> > slices(nbSlice);
        std::vector< std::vector<unsigned int> > slicesReduce(nbSlice);
        //used to compute centroid
        //if we use center of faces to compute centroid
        //it might not be precise
        //std::vector< std::set<Z3i::RealPoint> > slicePoints(nbSlice);
        //use median?

        std::map<unsigned int, unsigned int> selectedPoints;

        getHomogenityCloud(points, 5, ptLow, ptUp, selectedPoints);
std::cout<<"selected Point size:" << selectedPoints.size()<<std::endl;
std::cout<<"points size:" << points.size()<<std::endl;
        for (std::map<unsigned int, unsigned int>::iterator it=selectedPoints.begin(); it!=selectedPoints.end(); ++it){
            std::cout << it->first << " => " << it->second << '\n';
            unsigned int index = it->second;
            int sliceId = (points[index][2] - minZ) / length;
            slicesReduce[sliceId].push_back(index);
        }

        for(int i = 0; i < points.size(); i++){
            int sliceId = (points[i][2] - minZ) / length;
            slices[sliceId].push_back(i);
        }

std::cout<<"slices size:" <<slices.size()<<std::endl;
std::cout<<"slices reduced size:" <<slicesReduce.size()<<std::endl;
        segmentBranches(points, slices, slicesReduce, length, trunk, branch);
    }
    ////end segment for point cloud///

    static void
    segmentBranches(const DGtal::Mesh<Z3i::RealPoint> &mesh, double length,
            std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch){
        Face f0 = mesh.getFace(0);
        double minZ = mesh.getFaceBarycenter(0)[2];
        double maxZ = mesh.getFaceBarycenter(0)[2];
        std::vector<Z3i::RealPoint> centers(mesh.nbFaces());
        for( int i = 0; i < mesh.nbFaces(); i++){
            Z3i::RealPoint p = mesh.getFaceBarycenter(i);
            centers[i] = p;
            if(p[2] < minZ) {
                minZ = p[2];
            }
            if(p[2] > maxZ){
                maxZ = p[2];
            }
        }
        int nbSlice = (maxZ - minZ) / length + 1;

        std::vector< std::vector<unsigned int> > slices(nbSlice);
        std::vector< std::vector<unsigned int> > slicesReduce(nbSlice);
        //used to compute centroid
        //if we use center of faces to compute centroid
        //it might not be precise 
        //std::vector< std::set<Z3i::RealPoint> > slicePoints(nbSlice);
        //use median?

        std::map<unsigned int, unsigned int> selectedFaces;
        getHomogenityCloud(mesh, 10, selectedFaces);

        for (std::map<unsigned int, unsigned int>::iterator it=selectedFaces.begin(); it!=selectedFaces.end(); ++it){
            //std::cout << it->first << " => " << it->second << '\n';
            unsigned int index = it->second;
            int sliceId = (centers[index][2] - minZ) / length;
            slicesReduce[sliceId].push_back(index);
        }
        for(int i = 0; i < centers.size(); i++){
            int sliceId = (centers[i][2] - minZ) / length;
            slices[sliceId].push_back(i);
        }

        segmentBranches(centers, slices, slicesReduce, length, trunk, branch);
    }

    /////////////////////////////////////////////////////////////////////////

    static void simpleSubSample(const std::vector<Z3i::RealPoint> &pointCloud, const double &resolution, std::vector<unsigned int> &subsampleCloud){
        assert(pointCloud.size() > 0);
        Z3i::RealPoint p0 = pointCloud.at(0);
        Z3i::RealPoint ptLow = p0;
        Z3i::RealPoint ptUp = p0;

        
        for(unsigned int i = 0; i < pointCloud.size(); i++){
            Z3i::RealPoint pi = pointCloud.at(i);
            ptLow = ptLow.inf(pi);
            ptUp = ptUp.sup(pi);
        }

        double minZ = ptLow[2];
        double maxZ = ptUp[2];

        std::map<unsigned int, unsigned int> selectedFaces;
        getHomogenityCloud(pointCloud, resolution, ptLow, ptUp, selectedFaces);

        for (std::map<unsigned int, unsigned int>::iterator it=selectedFaces.begin(); it!=selectedFaces.end(); ++it){
            //std::cout << it->first << " => " << it->second << '\n';
            unsigned int index = it->second;
            subsampleCloud.push_back(index);
        }

    }
    
    /////////////////////////////////////////////////////////////////////////
    static void simpleSubSample(const std::vector<Z3i::RealPoint> &pointCloud, const double &resolution, std::vector<Z3i::RealPoint> &subsampleCloud){
        assert(pointCloud.size() > 0);
        Z3i::RealPoint p0 = pointCloud.at(0);
        Z3i::RealPoint ptLow = p0;
        Z3i::RealPoint ptUp = p0;

        
        for(unsigned int i = 0; i < pointCloud.size(); i++){
            Z3i::RealPoint pi = pointCloud.at(i);
            ptLow = ptLow.inf(pi);
            ptUp = ptUp.sup(pi);
        }

        double minZ = ptLow[2];
        double maxZ = ptUp[2];

        std::map<unsigned int, unsigned int> selectedFaces;
        getHomogenityCloud(pointCloud, resolution, ptLow, ptUp, selectedFaces);

        for (std::map<unsigned int, unsigned int>::iterator it=selectedFaces.begin(); it!=selectedFaces.end(); ++it){
            //std::cout << it->first << " => " << it->second << '\n';
            unsigned int index = it->second;
            subsampleCloud.push_back(pointCloud[index]);
        }
    }

    //
    static void subSample(const std::vector<Z3i::RealPoint> &pointCloud, const double &resolution, const double &radius, std::vector<unsigned int> &subsampleCloud){
        assert(pointCloud.size() > 0);
        Z3i::RealPoint p0 = pointCloud.at(0);
        Z3i::RealPoint ptLow = p0;
        Z3i::RealPoint ptUp = p0;

        
        for(unsigned int i = 0; i < pointCloud.size(); i++){
            Z3i::RealPoint pi = pointCloud.at(i);
            ptLow = ptLow.inf(pi);
            ptUp = ptUp.sup(pi);
        }


        double minZ = ptLow[2];
        double maxZ = ptUp[2];

        //std::cout<< minZ << " vs "<< maxZ<<std::endl;
        int nbSlice = (maxZ - minZ) / resolution + 1;

        std::vector< std::vector<unsigned int> > slices(nbSlice);
        std::vector< std::vector<unsigned int> > slicesReduce(nbSlice);
        //used to compute centroid
        //if we use center of faces to compute centroid
        //it might not be precise 
        //std::vector< std::set<Z3i::RealPoint> > slicePoints(nbSlice);
        //use median?

        std::map<unsigned int, unsigned int> selectedFaces;
        getHomogenityCloud(pointCloud, 5, ptLow, ptUp, selectedFaces);

        for (std::map<unsigned int, unsigned int>::iterator it=selectedFaces.begin(); it!=selectedFaces.end(); ++it){
            //std::cout << it->first << " => " << it->second << '\n';
            unsigned int index = it->second;
            int sliceId = (pointCloud[index][2] - minZ) / resolution;
            slicesReduce[sliceId].push_back(index);
        }

        for(int i = 0; i < pointCloud.size(); i++){
            int sliceId = (pointCloud[i][2] - minZ) / resolution;
            slices[sliceId].push_back(i);
        }

        Z3i::RealPoint oz(0,0,1);

        double maxDotProduct = 1 - std::numeric_limits<double>::epsilon(); 
        double dummyR = std::numeric_limits<double>::max();

        for(int i = 0; i < nbSlice; i++){
            if( slices[i].size() == 0 || slicesReduce[i].size() < NB_POINTS_MIN ){
                std::cout<< "SegmentationHelper.h::subSample : insuficient points: slice.size = "<< slices[i].size()<<"  slicesReduce.size="<< slicesReduce[i].size() <<std::endl;
                continue;
            }
            std::vector<unsigned int> slice = slices[i];

std::cout<<"processing slide: "<<  i<< "/"<< nbSlice <<"with "<< slicesReduce[i].size()<<" points "<< std::endl;
            CylinderModel cylModel = fitCylinderByCircle2D(pointCloud, slicesReduce[i], radius);
            std::cout<< "estimate cyl radii:"<< cylModel.radius<<std::endl;
            if(cylModel.radius <= 0){
                continue;
            }
            double angleStep = resolution/cylModel.radius;

            //first point as marqueur
            Z3i::RealPoint ma = getRadialVector(cylModel, pointCloud[slice[0]]);
            int nbSector = 2*M_PI / angleStep + 1;

            std::vector<unsigned int> mIndex(nbSector);

            std::vector<double> minR(nbSector, dummyR);
            for(int pind = 1; pind < slice.size(); pind++){
                unsigned int pId = slice[pind];
                Z3i::RealPoint vectRadial = getRadialVector(cylModel, pointCloud[pId]);
                //angle bt this point and ma

                /*
                std::cout<<"culModel: "<< cylModel.radius <<std::endl;
                std::cout<<"vectRadial: "<< vectRadial[0] << "  "<< vectRadial[1] << "  "<< vectRadial[2]<<std::endl;
                */
                double dotProduct = vectRadial.dot(ma)/ ma.norm()/vectRadial.norm();
                //prevent dot product > 1 due to floating point operations
                if(dotProduct > maxDotProduct){
                    dotProduct = maxDotProduct;
                }
                double angle = acos(dotProduct);
                //Z3i::RealPoint crossProduct = vectMarks[segmentId].crossProduct(vectRadial);
                Z3i::RealPoint u = ma.crossProduct(oz);
                if (u.dot(vectRadial) < 0){
                    angle = 2 * M_PI - angle;
                }

                //current sector
                int sectId = angle / angleStep;
                double curR = vectRadial.norm();
                /**
                  std::cout<<"angle:" << angle<<" step "<<angleStep<<std::endl;
                  std::cout<<"sectId:" << sectId<<" nb "<<nbSector<<std::endl;
                  std::cout.precision(std::numeric_limits<double>::max_digits10);
                  std::cout<<"acos:" << vectRadial.dot(ma)/ ma.norm()/vectRadial.norm()<<std::endl; 
                  */

                assert(sectId < nbSector);
                if (sectId < 0 && sectId >= nbSector){
                    //continue;
                }
                if(minR[sectId] > curR){
                    minR[sectId] = curR;
                    mIndex[sectId] = pId;
                }
            }
            for(int mInd = 0; mInd < mIndex.size(); mInd++){
                //if(minR[mInd] < dummyR){
                    subsampleCloud.push_back(mIndex[mInd]);
                //}
            }
        }
    }

    /*
    static void subSample(const DGtal::Mesh<Z3i::RealPoint> &mesh, const double &resolution, const double &radius,  std::vector<Z3i::RealPoint> &subsampleCloud){
        Face f0 = mesh.getFace(0);
        double minZ = mesh.getFaceBarycenter(0)[2];
        double maxZ = mesh.getFaceBarycenter(0)[2];
        std::vector<Z3i::RealPoint> centers(mesh.nbFaces());
        for( int i = 0; i < mesh.nbFaces(); i++){
            Z3i::RealPoint p = mesh.getFaceBarycenter(i);
            centers[i] = p;
            if(p[2] < minZ) {
                minZ = p[2];
            }
            if(p[2] > maxZ){
                maxZ = p[2];
            }
        }
        int nbSlice = (maxZ - minZ) / resolution + 1;

        std::vector< std::vector<unsigned int> > slices(nbSlice);
        std::vector< std::vector<unsigned int> > slicesReduce(nbSlice);
        //used to compute centroid
        //if we use center of faces to compute centroid
        //it might not be precise 
        //std::vector< std::set<Z3i::RealPoint> > slicePoints(nbSlice);
        //use median?

        std::map<unsigned int, unsigned int> selectedFaces;
        getHomogenityCloud(mesh, 10, selectedFaces);

        for (std::map<unsigned int, unsigned int>::iterator it=selectedFaces.begin(); it!=selectedFaces.end(); ++it){
            //std::cout << it->first << " => " << it->second << '\n';
            unsigned int index = it->second;
            int sliceId = (centers[index][2] - minZ) / resolution;
            slicesReduce[sliceId].push_back(index);
        }
        for(int i = 0; i < centers.size(); i++){
            int sliceId = (centers[i][2] - minZ) / resolution;
            slices[sliceId].push_back(i);
        }

        Z3i::RealPoint oz(0,0,1);

        double dummyR = 10000;
        //std::vector<unsigned int> filtered;
        //std::vector<bool> onTrunks(centers.size(), false);
        for(int i = 0; i < nbSlice; i++){
            if( slices[i].size() == 0 || slicesReduce[i].size() <100){
                std::cout<< "empty slice: "<< i <<std::endl;
                continue;
            }
            std::vector<unsigned int> slice = slices[i];

std::cout<<"processing slide: "<<  i<< "/"<< nbSlice <<"with "<< slicesReduce[i].size()<<" points "<< std::endl;
            CylinderModel cylModel = fitCylinderByCircle2D(mesh, slicesReduce[i], radius);
            std::cout<< cylModel.radius<<std::endl;
            if(cylModel.radius <= 0){
                continue;
            }
            double angleStep = resolution/cylModel.radius;

            //first point as marqueur
            Z3i::RealPoint ma = getRadialVector(cylModel, centers[slice[0]]);
            int nbSector = 2*M_PI / angleStep + 1;

            std::cout<< nbSector<< "  "<< angleStep << std::endl;

            std::vector<unsigned int> mIndex(nbSector);

            std::vector<double> minR(nbSector, dummyR);
            for(int pind = 1; pind < slice.size(); pind++){
                unsigned int pId = slice[pind];
                Z3i::RealPoint vectRadial = getRadialVector(cylModel, centers[pId]);
                //angle bt this point and ma

                double angle = acos(vectRadial.dot(ma)/ ma.norm()/vectRadial.norm());
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
                    subsampleCloud.push_back(mesh.getFaceBarycenter(mIndex[mInd]));
                    //filtered.push_back(mIndex[mInd]);
                    //onTrunks[mIndex[mInd]] = true;
                }
            }
        }
    }*/

    static CylinderModel fitCylinder(const std::vector<Z3i::RealPoint> &points, const double &radius){
        pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);
        pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
        for(unsigned int pointId = 0; pointId < points.size(); pointId++){
            Z3i::RealPoint p = points.at(pointId);
            PointT pt(p[0], p[1], p[2]);
            cloud->points.push_back(pt);
        } 

        pcl::NormalEstimation<PointT, pcl::Normal> ne;
        pcl::SACSegmentationFromNormals<PointT, pcl::Normal> seg; 

        pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT> ());
        pcl::PointCloud<pcl::Normal>::Ptr cloudNormals (new pcl::PointCloud<pcl::Normal>);
        pcl::ModelCoefficients::Ptr coeffs(new pcl::ModelCoefficients);
        pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

        ne.setSearchMethod (tree);
        ne.setInputCloud (cloud);
        ne.setKSearch (30);
        ne.compute (*cloudNormals);

        Eigen::Vector3f axis(0,0,1);
        seg.setOptimizeCoefficients (true);
        seg.setModelType (pcl::SACMODEL_CYLINDER);
        seg.setMethodType (pcl::SAC_RANSAC);

        seg.setAxis(axis);
        seg.setEpsAngle(0.5);
        seg.setProbability(0.3);
        seg.setNormalDistanceWeight (0.2);
        seg.setMaxIterations (100000);
        seg.setDistanceThreshold (5);
        seg.setRadiusLimits (radius * 0.3, radius);
        seg.setInputCloud (cloud);
        seg.setInputNormals (cloudNormals);
        seg.segment (*inliers, *coeffs);

        double x1 = coeffs->values[0];
        double y1 = coeffs->values[1];
        double z1 = coeffs->values[2];
        double x2 = x1 + coeffs->values[3];
        double y2 = y1 + coeffs->values[4];
        double z2 = z1 + coeffs->values[5];
        double ra = coeffs->values[6];

        Z3i::RealPoint aPointOnCenter(x1, y1, z1);
        Z3i::RealPoint aDirection(coeffs->values[3], coeffs->values[4], coeffs->values[5]);

        CylinderModel model;
        model.point = aPointOnCenter;
        model.direction = aDirection;
        model.radius = ra;

        return model;
    }



    static CylinderModel fitCylinder(const std::vector<Z3i::RealPoint> &points, const std::vector<unsigned int> &slice, const double &radius){
        pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);
        pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
        for(unsigned int j = 0; j < slice.size(); j++){
            unsigned int pIndex = slice[j];
            Z3i::RealPoint p = points.at(pIndex);
            PointT pt(p[0], p[1], p[2]);
            cloud->points.push_back(pt);
        } 
        pcl::NormalEstimation<PointT, pcl::Normal> ne;
        pcl::SACSegmentationFromNormals<PointT, pcl::Normal> seg; 

        pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT> ());
        pcl::PointCloud<pcl::Normal>::Ptr cloudNormals (new pcl::PointCloud<pcl::Normal>);
        pcl::ModelCoefficients::Ptr coeffs(new pcl::ModelCoefficients);
        pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

        ne.setSearchMethod (tree);
        ne.setInputCloud (cloud);
        ne.setKSearch (30);
        ne.compute (*cloudNormals);

        Eigen::Vector3f axis(0,0,1);
        seg.setOptimizeCoefficients (true);
        seg.setModelType (pcl::SACMODEL_CYLINDER);
        seg.setMethodType (pcl::SAC_RANSAC);

        seg.setAxis(axis);
        seg.setEpsAngle(0.5);
        seg.setProbability(0.3);
        seg.setNormalDistanceWeight (0.2);
        seg.setMaxIterations (100000);
        seg.setDistanceThreshold (5);
        seg.setRadiusLimits (radius * 0.7, radius * 1.3);
        seg.setInputCloud (cloud);
        seg.setInputNormals (cloudNormals);
        seg.segment (*inliers, *coeffs);

        double x1 = coeffs->values[0];
        double y1 = coeffs->values[1];
        double z1 = coeffs->values[2];
        double x2 = x1 + coeffs->values[3];
        double y2 = y1 + coeffs->values[4];
        double z2 = z1 + coeffs->values[5];
        double ra = coeffs->values[6];

        Z3i::RealPoint aPointOnCenter(x1, y1, z1);
        Z3i::RealPoint aDirection(coeffs->values[3], coeffs->values[4], coeffs->values[5]);

        CylinderModel model;
        model.point = aPointOnCenter;
            model.direction = aDirection;
            model.radius = ra;

            return model;
    }

    static CylinderModel fitCylinderByCircle2D(const std::vector<Z3i::RealPoint> &points, const std::vector<unsigned int> &slice, double radius){
            pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);
            pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
            for(unsigned int j = 0; j < slice.size(); j++){
                unsigned int pIndex = slice[j];
                Z3i::RealPoint p = points.at(pIndex);
                PointT pt(p[0], p[1], p[2]);
                cloud->points.push_back(pt);
            } 

            //pcl::SampleConsensusModelCircle2D<PointT> sac(cloud);


            pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
            pcl::ModelCoefficients::Ptr coeffs(new pcl::ModelCoefficients);

            pcl::SACSegmentation<pcl::PointXYZ> seg;

            seg.setOptimizeCoefficients (true);
            seg.setModelType (pcl::SACMODEL_CIRCLE2D);
            seg.setMethodType (pcl::SAC_RANSAC);
            seg.setDistanceThreshold (5);
            //limit ?
            seg.setRadiusLimits (radius * 0.7, radius * 1.3);

            seg.setInputCloud (cloud);
            seg.segment (*inliers, *coeffs);


            double x0 = coeffs->values[0];
            double y0 = coeffs->values[1];
            double ra = coeffs->values[2];

            CylinderModel model;

            Z3i::RealPoint aPointOnCenter(x0, y0, 0);
            Z3i::RealPoint aDirection(0,0,1);
            model.point = aPointOnCenter;
            model.direction = aDirection;
            model.radius = ra;

            return model;
    }

    static CylinderModel fitCylinder(const DGtal::Mesh<Z3i::RealPoint> &mesh, const std::vector<unsigned int> &slice){
            pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);
            pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
            for(unsigned int j = 0; j < slice.size(); j++){
                unsigned int pIndex = slice[j];
                Z3i::RealPoint p = mesh.getFaceBarycenter(pIndex);
                PointT pt(p[0], p[1], p[2]);
                cloud->points.push_back(pt);
            } 
            pcl::NormalEstimation<PointT, pcl::Normal> ne;
            pcl::SACSegmentationFromNormals<PointT, pcl::Normal> seg; 

            pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT> ());
            pcl::PointCloud<pcl::Normal>::Ptr cloudNormals (new pcl::PointCloud<pcl::Normal>);
            pcl::ModelCoefficients::Ptr coeffs(new pcl::ModelCoefficients);
            pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

            ne.setSearchMethod (tree);
            ne.setInputCloud (cloud);
            ne.setKSearch (30);
            ne.compute (*cloudNormals);

            Eigen::Vector3f axis(0,0,1);
            seg.setOptimizeCoefficients (true);
            seg.setModelType (pcl::SACMODEL_CYLINDER);
            seg.setMethodType (pcl::SAC_RANSAC);

            seg.setAxis(axis);
            seg.setEpsAngle(0.5);
            seg.setProbability(0.3);
            seg.setNormalDistanceWeight (0.2);
            seg.setMaxIterations (100000);
            seg.setDistanceThreshold (5);
            seg.setRadiusLimits (130, 200);
            seg.setInputCloud (cloud);
            seg.setInputNormals (cloudNormals);
            seg.segment (*inliers, *coeffs);

            double x1 = coeffs->values[0];
            double y1 = coeffs->values[1];
            double z1 = coeffs->values[2];
            double x2 = x1 + coeffs->values[3];
            double y2 = y1 + coeffs->values[4];
            double z2 = z1 + coeffs->values[5];
            double ra = coeffs->values[6];

            Z3i::RealPoint aPointOnCenter(x1, y1, z1);
            Z3i::RealPoint aDirection(coeffs->values[3], coeffs->values[4], coeffs->values[5]);

            CylinderModel model;
            model.point = aPointOnCenter;
            model.direction = aDirection;
            model.radius = ra;

            return model;
    }

    static CylinderModel fitCylinderByCircle2D(const DGtal::Mesh<Z3i::RealPoint> &mesh, const std::vector<unsigned int> &slice, double radius){
            pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);
            pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
            for(unsigned int j = 0; j < slice.size(); j++){
                unsigned int pIndex = slice[j];
                Z3i::RealPoint p = mesh.getFaceBarycenter(pIndex);
                PointT pt(p[0], p[1], p[2]);
                cloud->points.push_back(pt);
            } 

            //pcl::SampleConsensusModelCircle2D<PointT> sac(cloud);

            
            pcl::SampleConsensusModelCircle2D<PointT>::Ptr sacModel (new pcl::SampleConsensusModelCircle2D<PointT> (cloud));

            // Create the RANSAC object
            pcl::RandomSampleConsensus<PointT> sac (sacModel, 5);

            sac.computeModel();
            Eigen::VectorXf coeff;
            sac.getModelCoefficients (coeff);

            CylinderModel model;


            Z3i::RealPoint aPointOnCenter(coeff[0], coeff[1], 0);
            Z3i::RealPoint aDirection(0,0,1);
            model.point = aPointOnCenter;
            model.direction = aDirection;
            model.radius = coeff[2];

            return model;
    }

    static Z3i::RealPoint getRadialVector(const CylinderModel &cylModel, const Z3i::RealPoint &p){
        Z3i::RealPoint aVect = p - cylModel.point;
        double dist = aVect.dot(cylModel.direction);
        Z3i::RealPoint prj = cylModel.point + dist * cylModel.direction;
        return p -  prj;
    }

    static void getHomogenityCloud(const DGtal::Mesh<Z3i::RealPoint> &mesh, const double &resolution, std::map<unsigned int, unsigned int> &indexMap){
        std::vector<bool> selected(mesh.nbFaces(), false);

        std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> boudingBox = mesh.getBoundingBox();
        Z3i::RealPoint ptLow = boudingBox.first;
        Z3i::RealPoint ptUp = boudingBox.second;

        //save only one point in a voxel 10mm
        double minx = ptLow[0];
        double miny = ptLow[1];
        double minz = ptLow[2];

        double dimx = ceil((ptUp[0] - minx)/resolution);
        double dimy = ceil((ptUp[1] - miny)/resolution);
        double dimz = ceil((ptUp[2] - minz)/resolution);

        for (unsigned int i = 0; i < mesh.nbFaces(); i++){
            Z3i::RealPoint point = mesh.getFaceBarycenter(i);
            unsigned int colx = (unsigned int) floor((point[0] - minx) / resolution);
            unsigned int liny = (unsigned int) floor((point[1] - miny) / resolution);
            unsigned int levz = (unsigned int) floor((point[2] - minz) / resolution);

            unsigned int grdIndex = levz * dimx * dimy + liny * dimx + colx;
            

            if( indexMap.find(grdIndex) != indexMap.end()){
                indexMap[grdIndex] = i;
            }else{
                unsigned int previousPointIndex  = indexMap[grdIndex];

                double gridx = minx + colx*resolution + resolution/2;
                double gridy = miny + liny*resolution + resolution/2;
                double gridz = minz + levz*resolution + resolution/2;

                double distance = pow(point[0] - gridx, 2) + pow(point[1] - gridy, 2) + pow(point[2] - gridz, 2);

                Z3i::RealPoint previousPoint = mesh.getFaceBarycenter(previousPointIndex);

                double previousDistance = pow(previousPoint[0] - gridx, 2) + pow(previousPoint[1] - gridy, 2) + pow(previousPoint[2] - gridz, 2);

                if (distance < previousDistance){
                    indexMap[grdIndex] = i;
                }

            }
        }
    }

    static void getHomogenityCloud(const std::vector<Z3i::RealPoint> &pointCloud, const double &resolution, 
            const Z3i::RealPoint &ptLow, const Z3i::RealPoint &ptUp, std::map<unsigned int, unsigned int> &indexMap){
//std::cout<<"Get Homogene"<<std::endl;
        //std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> boudingBox = mesh.getBoundingBox();
        //Z3i::RealPoint ptLow = boudingBox.first;
        //Z3i::RealPoint ptUp = boudingBox.second;
        //save only one point in a voxel 10mm
        double minx = ptLow[0];
        double miny = ptLow[1];
        double minz = ptLow[2];
        double maxx = ptUp[0];
        double maxy = ptUp[1];
        double maxz = ptUp[2];

        double dimx = ceil((maxx - minx)/resolution);
        double dimy = ceil((maxy - miny)/resolution);
        double dimz = ceil((maxz - minz)/resolution);

        for (unsigned int i = 0; i < pointCloud.size(); i++){
            Z3i::RealPoint point = pointCloud.at(i);
            unsigned int colx = (unsigned int) floor((point[0] - minx) / resolution);
            unsigned int liny = (unsigned int) floor((point[1] - miny) / resolution);
            unsigned int levz = (unsigned int) floor((point[2] - minz) / resolution);

            unsigned int grdIndex = levz * dimx * dimy + liny * dimx + colx;
            

            if( indexMap.find(grdIndex) != indexMap.end()){
                indexMap[grdIndex] = i;
            }else{
                unsigned int previousPointIndex  = indexMap[grdIndex];

                double gridx = minx + colx*resolution + resolution/2;
                double gridy = miny + liny*resolution + resolution/2;
                double gridz = minz + levz*resolution + resolution/2;

                double distance = pow(point[0] - gridx, 2) + pow(point[1] - gridy, 2) + pow(point[2] - gridz, 2);

                Z3i::RealPoint previousPoint = pointCloud.at(previousPointIndex);

                double previousDistance = pow(previousPoint[0] - gridx, 2) + pow(previousPoint[1] - gridy, 2) + pow(previousPoint[2] - gridz, 2);

                if (distance < previousDistance){
                    indexMap[grdIndex] = i;
                }

            }
        }
        //std::cout<<"map size:"<<indexMap.size()<<std::endl;
    }


    static void voxelize(const std::vector<Z3i::RealPoint> &pointCloud, const double &resolution, 
            const Z3i::RealPoint &ptLow, const Z3i::RealPoint &ptUp, std::map<unsigned int, std::vector<unsigned int> > &indexMap){
std::cout<<"Get Homogene"<<std::endl;
        //std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> boudingBox = mesh.getBoundingBox();
        //Z3i::RealPoint ptLow = boudingBox.first;
        //Z3i::RealPoint ptUp = boudingBox.second;
        //save only one point in a voxel 10mm
        double minx = ptLow[0];
        double miny = ptLow[1];
        double minz = ptLow[2];
        double maxx = ptUp[0];
        double maxy = ptUp[1];
        double maxz = ptUp[2];

        double dimx = ceil((maxx - minx)/resolution);
        double dimy = ceil((maxy - miny)/resolution);
        double dimz = ceil((maxz - minz)/resolution);

        for (unsigned int i = 0; i < pointCloud.size(); i++){
            Z3i::RealPoint point = pointCloud.at(i);
            unsigned int colx = (unsigned int) floor((point[0] - minx) / resolution);
            unsigned int liny = (unsigned int) floor((point[1] - miny) / resolution);
            unsigned int levz = (unsigned int) floor((point[2] - minz) / resolution);

            unsigned int grdIndex = levz * dimx * dimy + liny * dimx + colx;
            
            std::vector<unsigned int> indexes;
            const auto result = indexMap.emplace(grdIndex, indexes);
            result.first->second.push_back(i);
        }
        std::cout<<"voxels size:"<<indexMap.size()<<std::endl;
    }

private:
    static
    void segmentBranches(const std::vector<Z3i::RealPoint> &points,const std::vector<CylindricalPoint> &cylPoints, const std::vector<std::vector<unsigned int> > &slices,
            const std::vector<double> &sliceRadius,
            const double &length,
            const std::vector<Z3i::RealPoint> &centerline,
            std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch
    ){

        std::cout<<"SegmentationHelper::segmentBranches:begin process"<<std::endl;
        int nbSlice = slices.size();
        std::cout<<"nbSlice:"<<nbSlice<<std::endl;
        Z3i::RealPoint oz(0,0,1);

        double dummyR = std::numeric_limits<double>::max();
        std::vector<unsigned int> filtered;
        std::vector<bool> onTrunks(cylPoints.size(), false);

        for(int i = 0; i < nbSlice; i++){
            if( slices[i].size() == 0 ){
                std::cout<< "empty slice: "<< i <<std::endl;
                continue;
            }
            std::vector<unsigned int> slice = slices[i];

std::cout<<"processing slide: "<<  i<< "/"<< nbSlice <<"with "<< slices[i].size()<<" points "<< std::endl;

            double angleStep = length/sliceRadius[i];

            //first point as marqueur
            //Z3i::RealPoint ma = getRadialVector(cylModel, points[slice[0]]);
            int nbSector = 2*M_PI / angleStep + 1;

            std::cout<< nbSector<< "  "<< angleStep << std::endl;

            std::vector<unsigned int> mIndex(nbSector);

            std::vector<double> minR(nbSector, dummyR);
            for(int pind = 1; pind < slice.size(); pind++){
                unsigned int pId = slice[pind];
                CylindricalPoint cp = cylPoints[pId];

                //current sector
                int sectId = cp.angle / angleStep;
                if(minR[sectId] > cp.radius){
                    minR[sectId] = cp.radius;
                    mIndex[sectId] = pId;
                }
            }
            for(int mInd = 0; mInd < mIndex.size(); mInd++){
                if(minR[mInd] < 2*sliceRadius[i]){
                    filtered.push_back(mIndex[mInd]);
                    onTrunks[mIndex[mInd]] = true;
                }
            }
        }

        //Now push back the nearest point in a sphere!!!
        std::cout<<"begin extending"<<std::endl;
        //std::vector<unsigned int> filteredExtends;
        //Search for neighbours
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        for(Z3i::RealPoint p: points){
            //p /= voxelSize;
            pcl::PointXYZ pxyz(p[0], p[1], p[2]);
            cloud->points.push_back(pxyz);
        }


        pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
        kdtree.setInputCloud (cloud);

        std::vector<int> pointIdx;
        std::vector<float> pointRadiusSquaredDistance;

        //search for points in the trunk set in a sphere with a radius of patch length * sqrt(2) or 2 for sure
        double searchRadius = length*BRANCH_SEG_SEARCH_COEFF;
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
        for (int i = 0; i < points.size(); i++){
            if(onTrunks[i]){
                trunk.push_back(i);
            }else{
                branch.push_back(i);
            }
        }
    }


    static
    void segmentBranches(const std::vector<Z3i::RealPoint> &points, const std::vector<std::vector<unsigned int> > &slices,
            const std::vector<std::vector<unsigned int> > &slicesReduce, double length,
            std::vector<unsigned int> &trunk, std::vector<unsigned int> &branch
    ){

        int nbSlice = slices.size();
        Z3i::RealPoint oz(0,0,1);

        double dummyR = std::numeric_limits<double>::max();
        std::vector<unsigned int> filtered;
        std::vector<bool> onTrunks(points.size(), false);

        for(int i = 0; i < nbSlice; i++){
            if( slices[i].size() == 0 || slicesReduce[i].size() <100){
                std::cout<< "empty slice: "<< i <<std::endl;
                std::cout<< "slice reduce size:" <<slicesReduce[i].size()<<std::endl;
                continue;
            }
            std::vector<unsigned int> slice = slices[i];

std::cout<<"processing slide: "<<  i<< "/"<< nbSlice <<"with "<< slicesReduce[i].size()<<" points "<< std::endl;
            CylinderModel cylModel = fitCylinderByCircle2D(points, slicesReduce[i], 200);
            std::cout<< cylModel.radius<<std::endl;
            if(cylModel.radius <= 0){
                continue;
            }

            double angleStep = length/cylModel.radius;

            //first point as marqueur
            Z3i::RealPoint ma = getRadialVector(cylModel, points[slice[0]]);
            int nbSector = 2*M_PI / angleStep + 1;

            std::cout<< nbSector<< "  "<< angleStep << std::endl;

            std::vector<unsigned int> mIndex(nbSector);

            std::vector<double> minR(nbSector, dummyR);
            for(int pind = 1; pind < slice.size(); pind++){
                unsigned int pId = slice[pind];
                Z3i::RealPoint vectRadial = getRadialVector(cylModel, points[pId]);
                //angle bt this point and ma

                double angle = acos(vectRadial.dot(ma)/ ma.norm()/vectRadial.norm());
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

        //Now push back the nearest point in a sphere!!!
        std::cout<<"begin extending"<<std::endl;
        //std::vector<unsigned int> filteredExtends;
        //Search for neighbours
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        for(Z3i::RealPoint p: points){
            //p /= voxelSize;
            pcl::PointXYZ pxyz(p[0], p[1], p[2]);
            cloud->points.push_back(pxyz);
        }


        pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
        kdtree.setInputCloud (cloud);

        std::vector<int> pointIdx;
        std::vector<float> pointRadiusSquaredDistance;
        double searchRadius = length + 10;
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
        for (int i = 0; i < points.size(); i++){
            if(onTrunks[i]){
                trunk.push_back(i);
            }else{
                branch.push_back(i);
            }
        }
    }
};
#endif

