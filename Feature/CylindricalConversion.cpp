#include "CylindricalConversion.h"

void
CylindricalConversion::init(){
    allocate();
    allocateExtra();

    computeBeginOfSegment();
    computeVectorMarks();
    computePlaneNormals();

    convertToCcs();
}
std::vector<CylindricalPoint>
CylindricalConversion::getPointsInCylindric(){
    return myPoints;
}



void 
CylindricalConversion::allocateExtra(){
}


void 
CylindricalConversion::computeEquations(){
}

void
CylindricalConversion::computeDistances(){
}
