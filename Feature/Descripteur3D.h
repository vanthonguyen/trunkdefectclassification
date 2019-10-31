#ifndef MOMENTS_3D 
#define MOMENTS_3D

#include <iostream>
#include <vector>
#include <array>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

using namespace DGtal;

class Descripteur3D{
public:
    Descripteur3D(const std::vector<Z3i::RealPoint> &pc, const std::vector<unsigned int> &ind):
        pointCloud(pc), indices(ind){}
    /**
     * compute normalized moments invariant of point cloud j1, j2, j3
     * j1=u002 + u020 + u 002
     * ....
     */
    std::array<double, 3> compute();
    std::array<double, 3> compute(double &u200, double &u020, double &u002);
private:
    const std::vector<Z3i::RealPoint> &pointCloud;
    const std::vector<unsigned int> &indices;
    Z3i::RealPoint computeCentroid();
    double boundingboxVolume();
};
#endif //MOMENTS_3D
