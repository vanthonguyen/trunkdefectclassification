#ifndef CENTERLINE_H
#define CENTERLINE_H 

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/writers/VolWriter.h"
#include <DGtal/base/BasicFunctors.h>
#include <DGtal/shapes/implicit/ImplicitBall.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/kernel/BasicPointFunctors.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/shapes/EuclideanShapesDecorator.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>

#include "DGtal/shapes/Mesh.h"


#include "SimpleNormalAccumulator.h"

///////////////////////////////////////////////////////////////////////////////
// class Centerline
/**
 * Description of class 'Centerline' <p>
 * 
 * @brief Class to compute centerline of log using volumic accumulation from normal vector field
 * and Splines to smoothing
 *
 *
 */


using namespace DGtal;

// types of image containers:
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned int> Image3D;
typedef ImageContainerBySTLVector<DGtal::Z3i::Domain, double> Image3DDouble;
typedef ImageContainerBySTLVector<Z3i::Domain,  Z3i::RealPoint> ImageVector;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

typedef ConstImageAdapter<Image3D, Z2i::Domain, functors::Point2DEmbedderIn3D<Z3i::Domain>,
                                 unsigned int, functors::Identity >  ImageAdapterExtractor;

enum RadiusEstimationType { mean, median };
/**
 * compare two Point by z value
 **/
struct PointOrderByZ {
      bool operator() (const Z3i::RealPoint &p1, const Z3i::RealPoint &p2) { return p1[2] < p2[2];}
};


/**
 * Private class to calculate points and domain to fit the input of SimpleNormalAccumulator
 * it might not be a optimise solution
 **/
/**class NormalComputation {
  // friend class the_class;
public: 
    NormalComputation(){}

    NormalComputation(const DGtal::Mesh<Z3i::RealPoint> &aMesh){

        Z3i::RealPoint lowPt = aMesh.getVertex(0);
        Z3i::RealPoint upPt = aMesh.getVertex(0);

        for ( int i = 0; i < aMesh.nbFaces(); i++ ){
            Face aFace = aMesh.getFace(i);

            //centroid
            Z3i::RealPoint cp = aMesh.getFaceBarycenter(i);
            startingPoints.push_back(cp);

            lowPt = lowPt.inf(cp);
            upPt = upPt.sup(cp);

            Z3i::RealPoint vectorNormal = ((p1-p0).crossProduct(p2 - p0)).getNormalized();
            pointsNormal.push_back(vectorNormal);
        }
        aDomain = Z3i::Domain(lowPt - DGtal::Z3i::RealPoint::diagonal(3*accRadius), upPt + DGtal::Z3i::RealPoint::diagonal(3*accRadius));

        // however that calculates x, y, and z from str I wouldn't know
    }

    const Z3i::Domain aDomain,
    const std::vector<Z3i::RealP oint> &pointsNormal,
    const std::vector<Z3i::RealPoint> &startingPoints,

}*/

//class Centerline : private NormalComputation, public SimpleNormalAccumulator{
class Centerline : public SimpleNormalAccumulator{

// ----------------------- Standard methods ------------------------------
public:
    //using SimpleNormalAccumulator::SimpleNormalAccumulator;
    Centerline(const double aRadius,
               const Z3i::Domain aDomain,
               const std::vector<Z3i::RealPoint> &pointsNormal,
               const std::vector<Z3i::RealPoint> &startingPoints,
               const double threshold,
               const int nbCp
            ): SimpleNormalAccumulator(aRadius, aDomain, pointsNormal, startingPoints), confidenceThreshold(threshold), nbControlPoint(nbCp){}
/*
    Centerline(const DGtal::Mesh<Z3i::RealPoint> &oriMesh,
               const double aRadius,
               const double threshold,
               const int nbCp
            ): NormalComputation(oriMesh), SimpleNormalAccumulator(aRadius, aDomain, pointsNormal, startingPoints), confidenceThreshold(threshold), nbControlPoint(nbCp){}

*/

    //return smoothed centerline
    std::vector<Z3i::RealPoint> compute();
    //std::vector<Z3i::RealPoint> interpolate()
    std::pair<Z3i::RealPoint, Z3i::RealPoint> computeLine();
    void computeRawLine(std::vector<Z3i::RealPoint> & rawLine);
    void computeRawLineWithLinearRegression(std::vector<Z3i::RealPoint> &rawLine);
    std::vector<Z3i::RealPoint> getRawLine();
    std::pair<Z3i::RealPoint, Z3i::RealPoint> computeLine(double &radius, RadiusEstimationType est);

//protected functions
protected:
//    std::vector<Z3i::RealPoint> getRawLine(std::vector<double> &radiis);
    double confidenceThreshold;
    int nbControlPoint;
};

#endif //end CENTERLINE_H
