#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <algorithm>    // std::min

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
#include "../Common/RlzPoint.h"
#include "../Feature/MomentCylindrical.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;



typedef DigitalSetBySTLSet<Z3i::Domain> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;



using namespace DGtal::Z3i;
using namespace DGtal;
typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;



#define MIN_DIFF 0.4

double estimateRadius(std::vector<CylindricalPoint> cylPoints){
    float minL = std::numeric_limits<float>::max();
    float minZ = std::numeric_limits<float>::max();
    float minR = std::numeric_limits<float>::max();
    float maxL = -std::numeric_limits<float>::max();
    float maxZ = -std::numeric_limits<float>::max();
    float maxR = -std::numeric_limits<float>::max();

   
    //c'est pas corect
    std::vector<double> radiis(cylPoints.size());
    for (unsigned int i = 0; i < cylPoints.size(); i++){
        CylindricalPoint p = cylPoints[i];
        radiis[i] = cylPoints[i].radius;
        if(p.radius > maxR){
            maxR = p.radius;
        }
        if(p.angle > maxL){
            maxL = p.angle;
        }
        if(p.height > maxZ){
            maxZ = p.height;
        }

        if(p.radius < minR){
            minR = p.radius;
        }

        if(p.angle < minL){
            minL = p.angle;
        }

        if(p.height < minZ){
            minZ = p.height;
        }

    }

    auto minmax = std::minmax_element(radiis.begin(), radiis.end());
    return (float)Statistic::getMode(radiis, radiis.at(minmax.first - radiis.begin()), radiis.at(minmax.second - radiis.begin()), 0.001);

}

std::vector<int> getExtremas(const std::vector<double> &hist){
    std::vector<int> extres;
//biminmax.push_back(0.0);
    //std::cout<<"xxxxxx"<<std::endl;
    for(unsigned int i = 1; i < hist.size() -1 ; i++){
        double di = hist[i] - hist[i - 1];
        double di1 = hist[i + 1] - hist[i];
//std::cout<<di<<"  "<<di1<<std::endl;
        //no
        if(std::abs(di) <= std::numeric_limits<double>::epsilon() || std::abs(di1) <= std::numeric_limits<double>::epsilon()){
            continue;
        }
        if(di * di1 > 0){
            continue;
        }
        extres.push_back(i);
        //std::cout<<hist[i]<<std::endl;
    }
    return extres;
}

std::vector<double> smooth(const std::vector<double> &heights){
    int winSize = 7;
    int winMidSize = winSize / 2;

    std::vector<double> hist(heights.size());
    for(unsigned int i = winMidSize; i < heights.size() - winMidSize; ++i){
        double mean = 0;
        for(unsigned int j = i - winMidSize; j <= (i + winMidSize); ++j){
            mean += heights[j];
        }

        hist[i] = mean / winSize;
    }
    return hist;
}
std::vector<double> getHeights(const std::vector<Z3i::RealPoint> &points, const std::vector<Z3i::RealPoint> &centerline, double res){
    Z3i::RealPoint seedPoint(0,1,0);
    CylindricalCoordinateSystem tmpCcs(centerline, seedPoint);
    CylindricalCoordinateSystem ccs(centerline, seedPoint);
    std::vector<CylindricalPoint> cylPoints(points.size());

    for (unsigned int i = 0; i < points.size(); i++){
        Z3i::RealPoint xyzPoint = points[i];
        cylPoints[i] = ccs.xyz2Cylindrical(xyzPoint);
    }

    double radius = estimateRadius(cylPoints);
    MomentCylindrical mc(cylPoints, radius);

    std::vector<RlzPoint> rlzPoints = mc.getRlzPointCloud();

    double minL = std::numeric_limits<double>::max();
    double maxL = -std::numeric_limits<double>::max();
    for(RlzPoint rlz: rlzPoints){
//std::cout<< rlz.l << "    ";
        minL = std::min(rlz.l, minL);
        maxL = std::max(rlz.l, maxL);
    }
    int nbBands = (maxL - minL)/res + 1;
    std::vector<std::vector<double> > zByBand(nbBands);
    for(RlzPoint rlz: rlzPoints){
        int bandId = (rlz.l - minL) / res;
        zByBand[bandId].push_back(rlz.z);
    }

    std::vector<double> heightByBand(nbBands);
    for (int i = 0; i < nbBands; i++){
        if(zByBand[i].size() == 0){
            heightByBand[i] = 0;
            continue;
        }
        double minZ = std::numeric_limits<double>::max();
        double maxZ = -std::numeric_limits<double>::max();
        for (double z: zByBand[i]){
            minZ = std::min(z, minZ);
            maxZ = std::max(z, maxZ);
        }
        heightByBand[i] = maxZ - minZ;
        //std::cout<<heightByBand[i]<<std::endl;
    }
    return heightByBand;
}
/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv){
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input defect scar in xyz.")
        ("resolution,r", po::value<double>()->default_value(1), "resolution.")
        ("centerline,l", po::value<std::string>(), "input trunk centerline branch.")
        ("output,o", po::value<std::string>()->default_value("scar"), "branch, ...");


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

    std::vector<Z3i::RealPoint> pointCloud = PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["input"].as<std::string>());
    std::vector<Z3i::RealPoint>  centerline= PointListReader<Z3i::RealPoint>::getPointsFromFile(vm["centerline"].as<std::string>());
//std::cout<<"nb points"<<pointCloud.size()<<std::endl;

    double resolution = vm["resolution"].as<double>();

    std::string outputPrefix = vm["output"].as<std::string>();
    std::vector<double> hs = getHeights(pointCloud, centerline, resolution);
    std::vector<double> histtemp = smooth(hs);
    std::vector<double> hist = smooth(histtemp);
    //std::cout<<hist.size()<<std::endl;
    std::vector<int> extremas = getExtremas(hist);
    double maxHeight = hist[0];
    int maxHeightId = 0;
    
    int maximaId = 0;

    for(int m = 0; m < extremas.size(); m++){
        int i = extremas[m];
        if(maxHeight < hist[i]){
            maxHeight = hist[i];
            maxHeightId = i;
            maximaId = m;
        }
        maxHeight = std::max(maxHeight, hist[i]);

//std::cout<<i<<" "<<hist[i]<<std::endl;
    }

    int leftId = 0;
    int rightId = hist.size() - 1;

//std::cout<<"xxx:"<< maxHeightId;
//std::cout<<"maxHeight:"<< maxHeight;

    if (maximaId > 0){
        for(int i = maximaId - 1; i >= 0; i -= 2 ){
            double minHeight = hist[extremas[i]];

//std::cout<<"minHeight:"<< minHeight<< std::endl;
            if(minHeight < (1 - MIN_DIFF) * maxHeight){
                leftId = extremas[i];
                break;
            }
        }
    }

    if (maximaId < extremas.size() - 1){
        for(int i = maximaId + 1; i < extremas.size(); i += 2 ){
            double minHeight = hist[extremas[i]];
            if(minHeight < (1 - MIN_DIFF)*maxHeight){
                rightId = extremas[i];
                break;
            }
        }
    }


    std::cout<< std::endl<< rightId <<"  "<< leftId << " "<< hist.size() <<std::endl;
//std::cout<< std::endl<< (rightId - leftId)*resolution << " "<< maxHeight <<std::endl;
    std::cout<< (rightId - leftId)*resolution << " "<< maxHeight <<std::endl;
   //find leftMin

    //std::cout<<"done ext"<<std::endl;
    IOHelper::export2Text(hist, outputPrefix);

  return 0;
}
