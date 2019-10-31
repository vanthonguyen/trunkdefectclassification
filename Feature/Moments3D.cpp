#include "Moments3D.h"

using namespace DGtal;

std::array<double, 3> Moments3D::compute(double &u200, double &u020, double &u002){
    double j1, j2, j3;
    Z3i::RealPoint centroid = computeCentroid();
   // Initalize the centralized moments
    double mu200 = 0, mu020 = 0, mu002 = 0, mu110 = 0, mu101 = 0, mu011  = 0;
   
   // Iterate over the nearest neighbors set
    for(unsigned int index : indices){
        // Demean the points
        Z3i::RealPoint tmpPoint = pointCloud.at(index) - centroid; 
        mu200 += tmpPoint[0] * tmpPoint[0];
        mu020 += tmpPoint[1] * tmpPoint[1];
        mu002 += tmpPoint[2] * tmpPoint[2];
        mu110 += tmpPoint[0] * tmpPoint[1];
        mu101 += tmpPoint[0] * tmpPoint[2];
        mu011 += tmpPoint[1] * tmpPoint[2];
    }
    // Normalize
    int m00 = indices.size();
    mu200 /= pow(m00,(2+0+0+3)/3.0);
    mu020 /= pow(m00, (0+2+0+3)/3.0);
    mu002 /= pow(m00, (0+0+2+3)/3.0);
    mu110 /= pow(m00, (1+1+0+3)/3.0);
    mu101 /= pow(m00, (1+0+1+3)/3.0);
    mu011 /= pow(m00, (0+1+1+3)/3.0);

std::cout<<"mu200:"<< mu200<<std::endl;
std::cout<<"mu020:"<< mu020<<std::endl;
std::cout<<"mu002:"<< mu002<<std::endl;
    // moment invariants
    j1 = mu200 + mu020 + mu002;
    j2 = mu200*mu020 + mu200*mu002 + mu020*mu002 - mu110*mu110 - mu101*mu101 - mu011*mu011;
    j3 = mu200*mu020*mu002 + 2*mu110*mu101*mu011 - mu002*mu110*mu110 - mu020*mu101*mu101 - mu200*mu011*mu011;
    std::array<double, 3> mm{ {j1, j2, j3} };
    u200 = mu200;
    u020 = mu020;
    u002 = mu002;
    return mm;
}

std::array<double, 3> Moments3D::compute(){
    double u200, u020, u002;
    return compute(u200, u020, u002);
}

Z3i::RealPoint Moments3D::computeCentroid(){
    Z3i::RealPoint centroid(0,0,0);
    if (indices.size() == 0){
        return centroid;
    }
    for(unsigned int index : indices){
        Z3i::RealPoint p = pointCloud.at(index);
        centroid += p;
    }
    return centroid/indices.size();
}
