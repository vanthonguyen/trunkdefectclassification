#include <algorithm>
#include "FeatureExtraction.h"
#include "../Common/Statistic.h"


std::vector<float> 
FeatureExtraction::compute(){
    //estimate radius
    //double std::vector
    //volume
    double minL = std::numeric_limits<double>::max();
    double minZ = std::numeric_limits<double>::max();
    double minR = std::numeric_limits<double>::max();
    double maxL = -std::numeric_limits<double>::max();
    double maxZ = -std::numeric_limits<double>::max();
    double maxR = -std::numeric_limits<double>::max();

   

    std::vector<double> radiis(points.size());
    for (int i = 0; i < points.size(); i++){
        CylindricalPoint p = points[i];
        radiis[i] = points[i].radius;
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
    double radiiEst = Statistic::getMode(radiis, radiis.at(minmax.first - radiis.begin()), radiis.at(minmax.second - radiis.begin()), 0.001);

    float meanRadius = (float)Statistic::getMean(radiis);
    float sdRadius = (float)Statistic::standardDeviation(radiis, meanRadius);

    
    double boundingboxVolume = (maxZ - minZ)*(maxL - minL)*(maxR - minR);
    double volumeRatio = points.size() / boundingboxVolume;

    MomentCylindrical mc(points, radiiEst);
    std::vector<float> feature{volumeRatio, mc.nu11, mc.nu02, mc.nu20, mc.nu12, mc.nu21, mc.nu03, mc.nu30};

    //add dimension ration
    //width/height
    feature.push_back(mc.width);
    feature.push_back(mc.width / mc.height);
    feature.push_back(mc.width / mc.depth);
    feature.push_back(meanRadius);
    feature.push_back(sdRadius);

    std::vector<RlzPoint> rlzPoints = mc.getRlzPointCloud();
    Eigen::Matrix3f evt;
    Eigen::Vector3f evv;
    Pca::compute(rlzPoints, evt, evv);
    
    for(int row = 0; row < evt.rows(); row++){
        for(int col = 0; col < evt.cols(); col++){
            feature.push_back(evt(row,col));
        }
    }

    feature.push_back(evv(0)/evv(2));
    feature.push_back(evv(1)/evv(2));
    //feature.push_back(evv(2));

    return feature;
}
