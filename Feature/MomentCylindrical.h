#ifndef MOMENT_CYLINDRICAL_H
#define MOMENT_CYLINDRICAL_H

#include<vector>
#include<limits>
#include <math.h>       // pow


#include "../Common/CylindricalPoint.h"
#include "../Common/RlzPoint.h"


class MomentCylindrical{
public:

    MomentCylindrical(const std::vector<CylindricalPoint> &cpoints, double r): cylindricalPoints(cpoints), radius(r){
        m00 = m10 = m01 = m20 = m11 = m02 = m30 = m21 = m12 = m03 =
            mu20 = mu11 = mu02 = mu30 = mu21 = mu12 = mu03 =
            nu20 = nu11 = nu02 = nu30 = nu21 = nu12 = nu03 = 0.;
        //convert from (r, t, z) to (r l z)
        nomalizedPointCloud();
        compute();
    }
    //
//    static double getMedian(std::vector<double> v);
public: 

    std::vector<RlzPoint> getRlzPointCloud();

    double  m00, m10, m01, m20, m11, m02, m30, m21, m12, m03;
    // central moments
    double  mu00, mu10, mu01, mu20, mu11, mu02, mu30, mu21, mu12, mu03;
    // central normalized moments
    double  nu00, nu10, nu01, nu20, nu11, nu02, nu30, nu21, nu12, nu03;

    double hu[7]; 

    double centerX;
    double centerY;

    double width;
    double height;
    double depth;

    //double meanRadius;
    //standard deviation
    //double sdRadius;

    //double minTheta;
    //double maxTheta;

    //double minZ;
    //double maxZ
private:
    void compute();
    void nomalizedPointCloud();
    
    double compute2DMoment(int p, int q );
    double compute2DCenterMoment(int p, int q);

    std::vector<CylindricalPoint> cylindricalPoints;
    std::vector<RlzPoint> rlzPoints;

    double radius;

};

#endif // STATISTIC_H
