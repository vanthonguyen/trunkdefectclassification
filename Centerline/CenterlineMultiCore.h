#ifndef CENTERLINE_MULTICORE_H
#define CENTERLINE_MULTICORE_H 

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


#include "Centerline.h"
#include "../Common/MultiThreadHelper.h"

///////////////////////////////////////////////////////////////////////////////
// class CenterlineMulticore
/**
 * Description of class 'CenterlineMulticore' <p>
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

class CenterlineMultiCore{

// ----------------------- Standard methods ------------------------------
public:
    CenterlineMultiCore(const double aRadius,
               const std::vector<Z3i::RealPoint> &pointsNormal,
               const std::vector<Z3i::RealPoint> &startingPoints,
               const int nbCp,
               int nbSeg = 1,
               double th = 0.0
            ): accRadius(aRadius), normals(pointsNormal), points(startingPoints), nbControlPoint(nbCp), nbSegment(nbSeg), threshold(th){}

    std::vector<Z3i::RealPoint> compute();

//protected functions
protected:

    std::vector<Z3i::RealPoint> getRawLine(const Z3i::RealPoint &lowPt, const Z3i::RealPoint &upPt);
    const std::vector<Z3i::RealPoint> &points;
    const std::vector<Z3i::RealPoint> &normals;
    double accRadius;
    int nbControlPoint;
    int nbSegment;
    int threshold;
};

#endif //end CENTERLINE_MULTICORE_H
