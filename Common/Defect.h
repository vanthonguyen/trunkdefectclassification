#ifndef DEFECT_H
#define DEFECT_H
#include <string>
//#include "DefectType.h"

class Defect {
    public:

        Defect(const double &w, const double &h, const double &z, const double &l):
            width(w), height(h), coordinateZ(z), coordinateL(l){}
        Defect(){}
        /*
        std::string getName(){
            return defectName;
        }
        */
        std::string toString(){
            std::string str = "type: " + std::to_string(type) + " coordinateZ: " + std::to_string(coordinateZ) +
                " coordinateL " + std::to_string(coordinateL) + " width (horizontal diameter): " + std::to_string(width) +
                " height (vertical diameter) :" + std::to_string(height);
            return str;
        }


        int species;
        double width; //horizontal dimension 
        double height; // vectical dimension
        double coordinateZ;
        double coordinateL;

        //branch
        double branchDiameter;

        //chinese moustache
        double moustacheHeight;
        double moustacheWidth;

        int type;
           protected:

        //Classification parameter
        std::string defectName;
};


#endif
