#ifndef CYLINDRDIRCAL_CONVERSION_H 
#define CYLINDRDIRCAL_CONVERSION_H

#include <utility>
#include <mutex>


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

#include "../Segmentation/SegmentationAbstract.h"
#include "../Common/CylindricalPoint.h"

//#define BIN_SIZE 0.01
#define eps 0.0000000001

using namespace DGtal;




class CylindricalConversion : public SegmentationAbstract {
    public:
        //CylindricalConversion(std::vector<Z3i::RealPoint> &aCloud, std::vector<Z3i::RealPoint> &aFib, double arcLen, int wh):
        //               pointCloud(aCloud), fiber(aFib), arcLength(arcLen), patchHeight(wh){
        //}
        using SegmentationAbstract::SegmentationAbstract;

        void init() override;

        std::vector<CylindricalPoint> getPointsInCylindric();

    protected:
        /** Brief 
         ** Allocate (resize) memory for array
         **/
        void allocateExtra() override;
        void computeEquations() override;
        void computeDistances() override;

 };
#endif
