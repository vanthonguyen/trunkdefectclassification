#ifndef RLZ_POINT_H
#define RLZ_POINT_H

//struct cylindrical point
class RlzPoint{
    public:
        RlzPoint(double radius, double arcLength, double height): r(radius), l(arcLength), z(height){}
        RlzPoint(){}
        double r;
        double l;
        double z;
};

#endif
