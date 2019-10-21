#ifndef DEFECT_TYPE_H
#define DEFECT_TYPE_H
#include <string>

#include "DGtal/base/Common.h"

//const static std::string[] DEFECT_TYPE={"Unknown", "chêne", "hêtre", "sapin", "pin", "merisier", "épicea"}
//const static std::string DEFECT_TYPE[] = {"Unknown", "Amas", "Branche", "Broussin", "Cicatrice", "Gourmand", "Sphéroblast", "Pico", "Unknown", "Ecorce"};
/*
const static std::string DEFECT_TYPE[] = {  "Inconnu", 
                                            "Branche", 
                                            "Cicatrice", 
                                            "Broussin", 
                                            "Gourmand", 
                                            "ASP", 
                                            "Ecorce", 
                                            "Etiquette", 
                                            "Inconnu 1", 
                                            "Inconnu 2"};
*/ 

//enum DefectType {UNKNOWN = 0, BRANCH, SCAR, MOUSTACHE, BURL, EPISHOOT, SMALLDEFECT, BARK, STICKER};
//                  0           1       2       3       4       5           6       7       8       9
enum DefectType {UNKNOWN = 0, BRANCH, SCAR, MOUSTACHE, BURL, EPISHOOT, SMALLDEFECT, BARK, STICKER, NOISE, NOTPROCESS};

const static DGtal::Color DEFECT_COLOR[] = {DGtal::Color(200, 200, 200), //0
                                            DGtal::Color(0,81,255),      //1
                                            DGtal::Color(0,255,4),       //2
                                            DGtal::Color(0,255,255),     //3
                                            DGtal::Color(227,255,0),     //4
                                            DGtal::Color(255,141,0),     //5
                                            DGtal::Color(255,26,0),      //6
                                            DGtal::Color(200,200,200),   //7
                                            DGtal::Color(200,200,200),   //8
                                            DGtal::Color(200,200,200),   //9
                                            DGtal::Color(4,4,4),}; //10


#endif
