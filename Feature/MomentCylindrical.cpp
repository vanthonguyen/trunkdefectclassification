#include<algorithm>
#include<iostream>

#include "../Common/Statistic.h"

#include "MomentCylindrical.h"

void MomentCylindrical::compute(){
    //compute centroid M00
    //double sumZ = 0;
    //spatial moments
    m00 = compute2DMoment(0, 0);
    m10 = compute2DMoment(1, 0);
    m01 = compute2DMoment(0, 1);
    m11 = compute2DMoment(1, 1);
    m20 = compute2DMoment(2, 0);
    m02 = compute2DMoment(0, 2);
    m21 = compute2DMoment(2, 1);
    m12 = compute2DMoment(1, 2);
    m03 = compute2DMoment(0, 3);
    m30 = compute2DMoment(3, 0);

    if( m00 > std::numeric_limits<double>::epsilon() ){
        centerX = m10 / m00;
        centerY = m01 / m00;
    }

    mu00 = m00;
    mu10 = 0;
    mu01 = 0;
    mu11 = compute2DCenterMoment(1, 1);
    mu20 = compute2DCenterMoment(2, 0);
    mu02 = compute2DCenterMoment(0, 2);
    mu21 = compute2DCenterMoment(2, 1);
    mu12 = compute2DCenterMoment(1, 2);
    mu03 = compute2DCenterMoment(0, 3);
    mu30 = compute2DCenterMoment(3, 0);

    //std::cout<<"m00"<< m00 <<std::endl;
    //nuij=muij/m00^((i+j)/2+1)
    nu00 = 1;
    nu01 = 0;
    nu10 = 0;
    if( m00 > std::numeric_limits<double>::epsilon() ){
        double divisor2 = pow(m00, 2);
        double divisor3 = pow(m00, 5/2);
        nu11 = mu11 / divisor2;
        nu20 = mu20 / divisor2;
        nu02 = mu02 / divisor2;
        nu12 = mu12 / divisor3;
        nu21 = mu21 / divisor3;
        nu03 = mu03 / divisor3;
        nu30 = mu30 / divisor3;
    }

    //hu
    hu[0] = nu20 + nu02;
    hu[1] = (nu20 - nu02)*(nu20 - nu02) + 4*nu11*nu11;
    hu[2] = (nu30 - 3*nu12) * (nu30 - 3*nu12) + (nu21 - 3*nu30)* (nu21 - 3*nu30);
    hu[3] = (nu30 + nu12)*(nu30 + nu12) + (nu21 + nu03)*(nu21 + nu03);
    hu[4] = (nu30 -3*nu12)*(nu30 + nu12)*((nu30 + nu12)*(nu30 + nu12) - 3*(nu21 + nu03)*(nu21 + nu03)) + (3*nu21 - nu03)*(nu21 + nu03)*(3*(nu30 + nu12)*(nu30 + nu12) - (nu21 + nu03)*(nu21 + nu03));
    hu[5] = (nu20 - nu02) * ( (nu30 + nu12)*(nu30 + nu12) - (nu21 + nu03)*(nu21 + nu03)) + 4*nu11*(nu30 + nu12)*(nu21 + nu03);
    hu[6] = (3*nu21 - nu03)*(nu21 + nu03)*(3*(nu30 + nu12)*(nu30 + nu12) - (nu21 + nu03)*(nu21 + nu03)) - (nu30 - 3*nu12)*(nu21 + nu03)*(3*(nu30 + nu12)*(nu30 + nu12) - (nu21 + nu03)*(nu21 + nu03));

    //mean and sd
    //std::vector<double> radiis(cylindricalPoints.size()); 
    //for (int i = 0; i < cylindricalPoints.size(); i++){
    //    radiis[i] = cylindricalPoints[i].radius;
    //}
    //meanRadius = Statistic::getMean(radiis);
    //sdRadius = Statistic::standardDeviation(radiis, meanRadius);
}
/**
 * p: order of arclength
 * q: order of height
 **/
double MomentCylindrical::compute2DMoment(int p, int q ){
    double Mpq = 0;
    for(RlzPoint rlz: rlzPoints){
        //r is intensity
        Mpq += pow(rlz.l, p) * pow(rlz.z, q)*rlz.r;
        //Mpq += pow(rlz.l, p) * pow(rlz.z, q) * rlz.r;
    }
    return Mpq;
}

double MomentCylindrical::compute2DCenterMoment(int p, int q){
    double Mupq = 0;
    for(RlzPoint rlz: rlzPoints){
        Mupq += pow(rlz.l - centerX, p) * pow(rlz.z - centerY, q) * rlz.r;
        //std::cout<< rlz.r<< "  "<< rlz.l<< "  "<< rlz.z <<std::endl;
    }
    return Mupq;
}


void MomentCylindrical::nomalizedPointCloud(){
    CylindricalPointThetaOrder order;
    std::sort(cylindricalPoints.begin(), cylindricalPoints.end(), order);

    bool isCounterClockWise = false;

    //min
    unsigned int border1 = 0;
    //max distance between two consecutive points
    double maxThetaDist = 0;
    for (unsigned int index = 0; index < cylindricalPoints.size() - 1; index++){
        double deltaTheta = cylindricalPoints[index + 1].angle - cylindricalPoints[index].angle;
        if(deltaTheta > maxThetaDist){
            maxThetaDist = deltaTheta;
            border1 = index;
        }
    }

    unsigned int border2 = border1 + 1;

    //resolve the circular problem, compare first and last element by clockwise distance
    if(cylindricalPoints[0].angle + 2*M_PI - cylindricalPoints[cylindricalPoints.size() - 1].angle > maxThetaDist){
        border1 = 0;
        border2 = cylindricalPoints.size() - 1;
    }else{
        isCounterClockWise = true;
    }

    double minZ = std::numeric_limits<double>::max();
    double maxZ = -std::numeric_limits<double>::max();

    double minR = std::numeric_limits<double>::max();
    double maxR = -std::numeric_limits<double>::max();

    double minTheta = std::numeric_limits<double>::max();
    double maxTheta = -std::numeric_limits<double>::max();

    for (CylindricalPoint cp: cylindricalPoints){
        if(cp.height < minZ){
            minZ = cp.height;
        }
        if(cp.height > maxZ){
            maxZ = cp.height;
        }
        if(cp.radius < minR){
            minR = cp.radius;
        }
        if(cp.radius > maxR){
            maxR = cp.radius;
        }
    }


    double borderTheta1 = cylindricalPoints[border1].angle;
    double borderTheta2 = cylindricalPoints[border2].angle;
//std::cout<<borderTheta1<<"@"<<borderTheta2<<std::endl;
    double thetaRange = borderTheta2 - borderTheta1;

    if( isCounterClockWise ){
        thetaRange = borderTheta1 + 2*M_PI - borderTheta2;
    }
    
    height = maxZ - minZ;
    depth = maxR - minR;
    width = thetaRange * radius;

    rlzPoints.resize(cylindricalPoints.size());
    for(unsigned int i = 0; i < cylindricalPoints.size(); i++){
        CylindricalPoint cp = cylindricalPoints[i];
        double th = cp.angle; 
        double arcLength = (th - borderTheta1)* radius;
        if(isCounterClockWise){
            double diff = th - borderTheta2;
            if (diff < 0){
                diff += 2*M_PI;
            }
            arcLength = diff * radius;
        }
        RlzPoint rlz( cp.radius, arcLength, cp.height);

        rlzPoints[i] = rlz;
        //cp.angle = arcLength;
        //cp.height = cp.height - minZ;
        //cylindricalPoints[i] = cp;
    }
}


std::vector<RlzPoint> MomentCylindrical::getRlzPointCloud(){
    return rlzPoints;
}
