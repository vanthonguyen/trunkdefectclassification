#ifndef DEFECT_SEGMENTATION_H
#define DEFECT_SEGMENTATION_H

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


#include "SegmentationAbstract.h"
#include "SegmentationHelper.h"
#include "../Common/CylindricalPoint.h"

//#define BIN_SIZE 0.01
//#define eps 0.0000000001

using namespace DGtal;




class DefectSegmentation : public SegmentationAbstract {
    public:
        //DefectSegmentation(std::vector<Z3i::RealPoint> &aCloud, std::vector<Z3i::RealPoint> &aFib, double arcLen, int wh):
        //               pointCloud(aCloud), centerline(aFib), arcLength(arcLen), patchHeight(wh){
        //}
        using SegmentationAbstract::SegmentationAbstract;

        void init() override;
        void lightInit();

        std::pair<double, double> getCoeffs(unsigned int index);
        std::pair<double, double> getCoeffs2(unsigned int index);

        static const std::pair<double, double> coeffZero;
    protected:
        /** Brief 
         ** Allocate (resize) memory for array
         **/
        void allocateExtra() override;
        /** Brief 
         * Compute the equation of the relation between distance to center line and z
         * of patches associated to points
         */
        void computeEquations() override;
        void computeEquationsFullPcl();
        void computeEquationsKdTree2D();
        void computeEquationsKdTree2DMultiThreads(int threadId, int nbThread,
            const pcl::KdTreeFLANN<pcl::PointXY> &qtree1,
            const pcl::KdTreeFLANN<pcl::PointXY> &qtree2,
            const std::vector<CylindricalPoint> &points1,
            const std::vector<CylindricalPoint> &points2,
            double minHeight, double maxHeight);
        /** Brief 
         * Call by computeEquations
         */
        void computeEquationsMultiThread(int threadId, int nbThread, const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, 
                             double minHeight, double maxHeight,
                             const std::vector<unsigned int> &subsampleCloudIds,
                             const std::vector<std::vector<unsigned int> > &voxelPointMaps);

        void computeDistances() override;

        void computeAngleLimit();
        //return coefficients a, b of y = ax + b
        std::pair<double, double> ransacLine(const std::vector<double> &xs, const std::vector<double> &ys);
        std::pair<double, double> ransacLine(const std::vector<double> &xs, const std::vector<double> &ys, double threshold, double resolution = 1, double maxPoints = 5, double ep=10);
        
        int computePatch(unsigned int pointId, double searchRadius, const pcl::KdTreeFLANN<pcl::PointXY> &qtree1, 
            const pcl::KdTreeFLANN<pcl::PointXY> &qtree2, std::vector<int> &pointIdx, std::vector<float> &pointRadiusSquaredDistance);


        //coefficients of regressed lines, one line for each windows
        //a window = some bands consecutives
        std::vector<std::pair<double, double> > coefficients;

        bool isComputed(int pointIndex);


        bool isInAlternativeTree(double angle);

        //convert angle (theta) in cylindrical coordinate to arc length and scale by ration height/arclength
        double getScaledX(double angle);
        
        //used for quadtree search / to overcome the circulation problem.
        double angleLimit;
        double yxRatio;

        //debug
        //std::map<int, bool> forExport;
        std::vector<bool> isForExport;
        int pid;
        std::recursive_mutex coeffMutex;
 };
#endif
