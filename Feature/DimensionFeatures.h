#ifndef MOMENT_CYLINDRICAL_H
#define MOMENT_CYLINDRICAL_H

#include<vector>
#include<limits>
#include <math.h>       // pow



struct CylindricalPointThetaOrder {
    bool operator() (Z3i::RealPoint p1, Z3i::RealPoint p2){
        return p1[1] < p2[1];
    }
};

class DimentionFeatures{
public:

    DimentionFeatures(const std::vector<CylindricalPoint> &cpoints){
    }
    //
//    static double getMedian(std::vector<double> v);
public: 
    double width;
    double height;
    double profond;

};

#endif // STATISTIC_H
