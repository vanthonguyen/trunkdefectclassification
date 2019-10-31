
#ifndef FEATURE_EXTRACTION_XYZ_H
#define FEATURE_EXTRACTION_XYZ_H 

#include <utility>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
//for eigen?
#include <pcl/common/common.h>

#include "MomentCylindrical.h"
#include "Pca.h"
#include "../Common/CylindricalCoordinateSystem.h"
#include "../Common/RlzPoint.h"

///////////////////////////////////////////////////////////////////////////////
// class CenterlineByZ
/**
 * Description of class 'CenterlineByZ' <p>
 * 
 * @brief Class to compute centerline of log using volumic accumulation from normal vector field
 * and Splines to smoothing
 *
 *
 */
// types of image containers:

using namespace DGtal;

class FeatureExtractionXyz{

// ----------------------- Standard methods ------------------------------
public:
    //using SimpleNormalAccumulator::SimpleNormalAccumulator;
    FeatureExtractionXyz(const std::vector<Z3i::RealPoint> &points, const std::vector<Z3i::RealPoint> &cl, const double &radii = 0): 
        xyzPoints(points), centerline(cl), radius(radii) {
        init();
    }

    /**
     * Compute the feature vector from : 
     * dimension x,y,z
     * moments ...
     * 
     **/
    std::vector<float> compute();


//protected functions
protected:
    void init();

    Z3i::RealPoint trunkAxis;
    std::vector<Z3i::RealPoint> xyzPoints;
    std::vector<Z3i::RealPoint> centerline;
    std::vector<CylindricalPoint> cylPoints;
    CylindricalPoint centroid;
    Z3i::RealPoint centroidXyz;
    Z3i::RealPoint pointOnSegment;
    int segmentId;
    double radius;

    bool testSkewness(Z3i::RealPoint p1, Z3i::RealPoint p2, Z3i::RealPoint p3, Z3i::RealPoint p4);

    double
    distanceBetween2Lines(Z3i::RealPoint p1, Z3i::RealPoint d1, Z3i::RealPoint p2, Z3i::RealPoint d2);

    //Pca related
    Z3i::RealPoint longitudinalEigenVector;
    double longitudinalEigenValue;

    Z3i::RealPoint tangentialEigenVector;
    double tangentialEigenValue;

    Z3i::RealPoint radialEigenVector;
    double radialEigenValue;
};

#endif //end CENTERLINE_H
