
#ifndef FEATURE_EXTRACTION_H
#define FEATURE_EXTRACTION_H 

//for eigen?
#include <pcl/common/common.h>

#include "MomentCylindrical.h"
#include "Pca.h"
#include "../Common/CylindricalPoint.h"
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

class FeatureExtraction{

// ----------------------- Standard methods ------------------------------
public:
    //using SimpleNormalAccumulator::SimpleNormalAccumulator;
    FeatureExtraction(const std::vector<CylindricalPoint> &ps ): points(ps) {}

    /**
     * Compute the feature vector from : 
     * dimension x,y,z
     * moments ...
     * 
     **/
    std::vector<float> compute();


//protected functions
protected:
    std::vector<CylindricalPoint> points;
};

#endif //end CENTERLINE_H
