#include <iostream>
#include <set>
#include <map>
#include <utility>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "IOHelper.h"

using namespace DGtal;

void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        Z3i::RealPoint p = v[i];
        outStream << p[0] << " "<< p[1] << " "<< p[2]<<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &pointCloud, 
        const std::vector<unsigned int> &indices, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < indices.size();i++){
        Z3i::RealPoint p = pointCloud.at(indices.at(i));
        outStream << p[0] << " "<< p[1] << " "<< p[2]<<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<CylindricalPoint> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        CylindricalPoint p = v[i];
        outStream << p.radius << " "<< p.angle << " "<< p.height <<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<std::pair<double, double> > &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        std::pair<double, double> pa = v.at(i);
        outStream << pa.first << " "<< pa.second <<std::endl;
    }
    outStream.close();
}

/*
void IOHelper::export2Text(const std::vector<double> &xs, const std::vector<double> &ys, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    if( xs.size() == ys.size() ){
        for(unsigned int i = 0; i < xs.size();i++){
            outStream << xs.at(i) << " "<< ys.at(i) <<std::endl;
        }
    }
    outStream.close();

}
*/

void IOHelper::readDistanceFromFile(const std::string &fileName, std::vector<double> &vectDistances){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
      if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
          vectDistances.push_back(std::stod(str));
      }
      getline(infile, str);
    }
  }


void IOHelper::export2OFF(const Mesh<Z3i::RealPoint> &mesh, std::string fileName){
    std::ofstream offMesh (fileName.c_str());
    DGtal::MeshWriter<Z3i::RealPoint>::export2OFF(offMesh, mesh);
    offMesh.close();
}


void IOHelper::export2OFF(const Mesh<Z3i::RealPoint> &mesh, const std::vector<unsigned int> &faceIds, const std::string &fileName){
    //slow
    std::set<unsigned int> vertexes;
    std::vector<unsigned int> vertexList;

    //map between old vertex id and new one
    std::map<unsigned int, unsigned int> vertexIdMap;
    unsigned int currentId = 0;
    for(unsigned int fId: faceIds){
        std::vector<unsigned int>  aFace = mesh.getFace(fId);
        for(unsigned int indexVertex : aFace){
            std::pair<std::set<unsigned int>::iterator, bool> ret = vertexes.insert(indexVertex);
            //new item
            if( ret.second == true ){
                vertexIdMap[indexVertex] = currentId;
                vertexList.push_back(indexVertex);
                currentId++;
            }
        }
    }
    assert(vertexIdMap.size() == vertexList.size());
    std::ofstream out (fileName.c_str());
    out << "OFF"<< std::endl;
    out << vertexes.size()  << " " << faceIds.size() << " " << 0 << " " << std::endl;

    for(unsigned int vId: vertexList){
        out << mesh.getVertex(vId)[0] << " " << mesh.getVertex(vId)[1] << " "<< mesh.getVertex(vId)[2] << std::endl;	
    }

    //print faces!!
    for(unsigned int fId: faceIds){
        std::vector<unsigned int>  aFace = mesh.getFace(fId);
        out << aFace.size() << " " ;
        for(unsigned int indexVertex : aFace){
            out << vertexIdMap[indexVertex] << " " ;
        }
        DGtal::Color col = mesh.getFaceColor(fId);
        if( mesh.isStoringFaceColors() ){
            out << " ";
            out << ((double) col.red())/255.0 << " "
                << ((double) col.green())/255.0 << " "<< ((double) col.blue())/255.0 
                << " " << ((double) col.alpha())/255.0 ;
        }  
        out << std::endl;
    }
    out.close();

}

void IOHelper::readIntsFromFile(const std::string &fileName, std::vector<int> &rs){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
        if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
            rs.push_back(std::stoi(str));
        }
        getline(infile, str);
    }
}

void IOHelper::readStringsFromFile(const std::string &fileName, std::vector<std::string> &rs){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
        if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
            rs.push_back(str);
        }
        getline(infile, str);
    }
}




void IOHelper::readFeatures(const std::string &fileName, std::vector<std::vector<float> > &rs){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
        if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
            std::vector<float> feature;
            std::istringstream s(str);
            float f;
            while (s >> f) {
                feature.push_back(f);
            }
            rs.push_back(feature);
        }
        getline(infile, str);
    }
}


