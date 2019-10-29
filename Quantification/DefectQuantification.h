#ifndef DEFECT_QUANTIFICATION_H
#define DEFECT_QUANTIFICATION_H

#include <utility>
#include <mutex>
#include <map>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/sample_consensus/sac_model_line.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>


#include "../Segmentation/DefectSegmentation.h"
#include "../Segmentation/SegmentationHelper.h"
#include "../Common/CylindricalPoint.h"
#include "../Common/Defect.h"

//#define BIN_SIZE 0.01
//#define eps 0.0000000001

using namespace DGtal;



/**
 * Quantify a defect segmented by precedence step
 *
 **/

class DefectQuantification {
    public:

        DefectQuantification(std::vector<Z3i::RealPoint> &aCloud, std::vector<Z3i::RealPoint> &aFib, Z3i::RealPoint seed,
                const std::vector<double> &segmentRadii,
                int spec, std::string rfFile, double pad = 1000):
            pointCloud(aCloud), centerline(aFib), seedPoint(seed), segmentRadius(segmentRadii),
            species(spec), rfFilePath(rfFile), padding(pad) {
        }

        /*@return: def type, dimension at insersection with the trunk, diameter if it is branch, ...*/
        Defect quantify();

    protected:
        //std::pair<double, double> getDimensionAtInsersection();
        double getBranchDiameter();
        std::vector<Z3i::RealPoint> &pointCloud;
        std::vector<Z3i::RealPoint> &centerline;
        //radius of trunk segment defined by centerline
        std::vector<double> segmentRadius;
        //double trunkRadiusAtDefect;

        Z3i::RealPoint seedPoint;
        //Classification parameter
        int species;
        double padding;
        std::string rfFilePath;
 };
#endif
