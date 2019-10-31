#ifndef DEFECT_BRANCH_SEGMENTATION_H
#define DEFECT_BRANCH_SEGMENTATION_H

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


#include "DefectSegmentation.h"
#include "SegmentationHelper.h"
#include "../Common/CylindricalPoint.h"

//#define BIN_SIZE 0.01
//#define eps 0.0000000001

using namespace DGtal;




class DefectBranchSegmentation : public DefectSegmentation {
    public:

        DefectBranchSegmentation(std::vector<Z3i::RealPoint> &aCloud, std::vector<Z3i::RealPoint> &aFib, 
            double arcLen, int wh, double bw, double secLength, int minClusterNb, double clusterTol, int voz = 3):
            DefectSegmentation(aCloud, aFib, arcLen, wh, bw), 
            sectorLength(secLength), minClusterSize(minClusterNb), clusterTolerance(clusterTol), voxelSize(voz) {
                init();
        }

        void segment(std::vector<std::vector<unsigned int> > &clusters);

        std::vector<unsigned int> getDefect() override;
        void init() override;
        std::vector<unsigned int> getTrunkIds();
        std::vector<unsigned int> getBranchIds();

    protected:

        //also threshold
        void computeDistances() override;

        void computeEquations() override;

        void computeEquationsMultiThread(int threadId, int nbThread, const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, 
                             double minHeight, double maxHeight,
                             const std::vector<unsigned int> &subsampleCloudIds,
                             const std::vector<std::vector<unsigned int> > &voxelPointMaps);

        //segment the trunk in z coordinate to reduce the volume of image 3D when doing accumulation
        //int nbSegment; 
        //resolution using for branch segmentation 
        double sectorLength;
        int minClusterSize;
        double clusterTolerance;

        std::vector<unsigned int> trunkIds;
        std::vector<unsigned int> branchIds;
        double threshold;
        int voxelSize;
 };
#endif
