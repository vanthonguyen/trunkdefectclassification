#include <algorithm>
#include "../Common/Statistic.h"

#include "FeatureExtractionXyz.h"
#include "DGtal/shapes/Mesh.h"
#include "../Common/IOHelper.h"

std::vector<float> 
FeatureExtractionXyz::compute(){
    //estimate radius
    //double std::vector
    //volume
    float minL = std::numeric_limits<float>::max();
    float minZ = std::numeric_limits<float>::max();
    float minR = std::numeric_limits<float>::max();
    float maxL = -std::numeric_limits<float>::max();
    float maxZ = -std::numeric_limits<float>::max();
    float maxR = -std::numeric_limits<float>::max();

   
    //c'est pas corect
    std::vector<double> radiis(cylPoints.size());
    for (int i = 0; i < cylPoints.size(); i++){
        CylindricalPoint p = cylPoints[i];
        radiis[i] = cylPoints[i].radius;
        if(p.radius > maxR){
            maxR = p.radius;
        }
        if(p.angle > maxL){
            maxL = p.angle;
        }
        if(p.height > maxZ){
            maxZ = p.height;
        }

        if(p.radius < minR){
            minR = p.radius;
        }

        if(p.angle < minL){
            minL = p.angle;
        }

        if(p.height < minZ){
            minZ = p.height;
        }

    }

    auto minmax = std::minmax_element(radiis.begin(), radiis.end());
    double radiiEst = (float)Statistic::getMode(radiis, radiis.at(minmax.first - radiis.begin()), radiis.at(minmax.second - radiis.begin()), 0.001);

    float meanRadius = (float)Statistic::getMean(radiis);
    float sdRadius = (float)Statistic::standardDeviation(radiis, meanRadius);

    
    float boundingboxVolume = (maxZ - minZ)*(maxL - minL)*(maxR - minR);
    float volumeRatio = cylPoints.size() / boundingboxVolume;

    if(radius == 0){
        radius = radiiEst;
//std::cout<<"RadiiEst:"<<radiiEst<<std::endl;
    }
//std::cout<<"rrrr"<<radius<<"@"<<radiiEst<<std::endl;
    MomentCylindrical mc(cylPoints, radius);
    std::vector<float> feature{volumeRatio, mc.nu11, mc.nu02, mc.nu20, mc.nu12, mc.nu21, mc.nu03, mc.nu30};

    //add dimension ration
    //width/height
    feature.push_back(mc.width);
    feature.push_back(mc.width / mc.height);
    feature.push_back(mc.width / mc.depth);
    feature.push_back(meanRadius);
    feature.push_back(sdRadius);

    std::vector<RlzPoint> rlzPoints = mc.getRlzPointCloud();
    Eigen::Matrix3f evt;
    Eigen::Vector3f evv;
    Pca::compute(xyzPoints, evt, evv);

    std::vector<double> dists(evt.cols(),0);
    std::vector<double> angles(evt.cols(),0);
    std::vector<double> dotProducts(evt.cols(),0);
    std::vector<Z3i::RealPoint> evtP(evt.cols());
    double minDist = std::numeric_limits<double>::max();
    int minIndex = 0;
    for(int col = 0; col < evt.cols(); col++){
        Z3i::RealPoint vect;
        for(int row = 0; row < evt.rows(); row++){
            vect[row] = evt(row, col);
        }
        bool isSkew = testSkewness(centroidXyz, centroidXyz + vect, pointOnSegment, pointOnSegment + trunkAxis);
        if( isSkew){
            dists[col] = distanceBetween2Lines(centroidXyz, vect, pointOnSegment, trunkAxis);
        }
        if(dists[col] < minDist){
            minDist = dists[col];
            minIndex = col;
        }
        //vect and trunkAxis were normalized
        angles[col] = acos(vect.dot(trunkAxis));
        dotProducts[col] = std::abs(vect.dot(trunkAxis));
        evtP[col] = vect;
    }


    int longitudinalIndex = 0;
    int tangentialIndex = 0;
    int radialIndex = minIndex;
    bool longitudinalIndexSet = false;

    for(int i = 0; i < evt.cols(); i++){
        if(i != radialIndex){
            if( !longitudinalIndexSet ){
                longitudinalIndex = i;
                longitudinalIndexSet = true;
            }else{
                tangentialIndex = i;
            }
        }
        //vector radial is the one with the minimum distance and ...
    }

    if(dotProducts[tangentialIndex] > dotProducts[longitudinalIndex]){
        int tmp = tangentialIndex;
        tangentialIndex = longitudinalIndex;
        longitudinalIndex = tmp;
    }

/*
    for(int row = 0; row < evt.rows(); row++){
        std::cout<<evt.row(row)<< "   val:"<< evv(row)<<std::endl;
        std::cout<<"dist:"<<dists[row]<<std::endl; 
        std::cout<<"angle:"<<angles[row]<<std::endl; 
    }

*/


//debug
/*
    std::cout<< "trunk axis: "<< trunkAxis[0] << " "<<trunkAxis[1] << "  "<< trunkAxis[2] <<std::endl;
    std::cout<< "longitudinal:"<< evt.row(longitudinalIndex) << "   val:"<< evv(longitudinalIndex)<<std::endl;
    std::cout<<"dist:"<<dists[longitudinalIndex]<<std::endl; 
    std::cout<<"angle:"<<angles[longitudinalIndex]<<std::endl; 
    std::cout<<"dot product"<<dotProducts[longitudinalIndex]<<std::endl; 

    std::cout<< "radial:"<< evt.row(radialIndex) << "   val:"<< evv(radialIndex)<<std::endl;
    std::cout<<"dist:"<<dists[radialIndex]<<std::endl; 
    std::cout<<"angle:"<<angles[radialIndex]<<std::endl; 
    std::cout<<"dot product"<<dotProducts[radialIndex]<<std::endl; 

    std::cout<< "tangentialIndex: "<< evt.row (tangentialIndex) << "   val:"<< evv(tangentialIndex)<<std::endl;
    std::cout<<"dist:"<<dists[tangentialIndex]<<std::endl; 
    std::cout<<"angle:"<<angles[tangentialIndex]<<std::endl; 
    std::cout<<"dot product"<<dotProducts[tangentialIndex]<<std::endl; 

    DGtal::Mesh<Z3i::RealPoint> mesh(true);
    Mesh<Z3i::RealPoint>::createTubularMesh(mesh, centerline, 5, 0.1, DGtal::Color::Cyan);
    //longitudinal
    std::vector<Z3i::RealPoint> longi{centroidXyz};
    Z3i::RealPoint longiEnd = centroidXyz + 100.0 * evtP[longitudinalIndex];
    longi.push_back(longiEnd);
    Mesh<Z3i::RealPoint>::createTubularMesh(mesh, longi, 3, 0.1, DGtal::Color::Red);

    std::vector<Z3i::RealPoint> radi{centroidXyz};
    Z3i::RealPoint radiEnd = centroidXyz + (double)(100*evv(radialIndex)/evv(longitudinalIndex))*evtP[radialIndex];
    radi.push_back(radiEnd);
    Mesh<Z3i::RealPoint>::createTubularMesh(mesh, radi, 3, 0.1, DGtal::Color::Green);

    std::vector<Z3i::RealPoint> tange{centroidXyz};
    Z3i::RealPoint tangeEnd = centroidXyz + (double)(100*evv(tangentialIndex)/evv(longitudinalIndex))*evtP[tangentialIndex];
    tange.push_back(tangeEnd);
    Mesh<Z3i::RealPoint>::createTubularMesh(mesh, tange, 3, 0.1, DGtal::Color::Blue);

    IOHelper::export2OFF(mesh, "dir.off");
*/
//end debug
    feature.push_back(evv(radialIndex)/evv(longitudinalIndex));
    feature.push_back(evv(tangentialIndex)/evv(longitudinalIndex));
    feature.push_back(angles[radialIndex]);

    /*
    //feature.push_back(evv(2));
    //histogram
    int nbInterval = 5;

    //??? 
    double binWidth = (maxR - minR)*1.00 / nbInterval;
    std::vector<int> histogram(nbInterval, 0);
    for(unsigned int i = 0; i < radiis.size(); i++){
        int index = (radiis[i] - minR)/binWidth;
//trace.info()<<"i: "<<index<<"  s: "<< histogram.size()<<std::endl;
        if(index >= histogram.size()){
            index = histogram.size() - 1;
        }

        histogram[index]++;
    }

    for(unsigned int i = 1; i < histogram.size(); i++){
        feature.push_back(1.0*histogram[i]/histogram[0]);
    }

    //histogram
    //??? 
    binWidth = (maxZ - minZ)*1.00 / nbInterval;
    std::vector<int> histZ(nbInterval, 0);
    for(unsigned int i = 0; i < radiis.size(); i++){
        int index = (cylPoints[i].height - minZ)/binWidth;
//trace.info()<<"i: "<<index<<"  s: "<< histogram.size()<<std::endl;
        if(index >= histZ.size()){
            index = histZ.size() - 1;
        }

        histZ[index]++;
    }

    for(unsigned int i = 1; i < histZ.size(); i++){
        feature.push_back(1.0*histZ[i]/histZ[0]);
    }
    */


    return feature;
}


void FeatureExtractionXyz::init(){
    //convert to cylindrical
    int nbPoint = xyzPoints.size();
    cylPoints.resize(nbPoint);
    CylindricalCoordinateSystem ccs(centerline, Z3i::RealPoint(0,1,0));
    centroidXyz = Z3i::RealPoint(0.0, 0.0, 0.0);
    for (int i = 0; i < xyzPoints.size(); i++){
        cylPoints[i] = ccs.xyz2Cylindrical(xyzPoints[i]);
        //std::cout<<"#####"<< cylPoints[i].radius << "#"<< cylPoints[i].angle<<"#"<<cylPoints[i].height<<std::endl;
        centroidXyz += xyzPoints[i];
    }
    centroidXyz /= xyzPoints.size();
    segmentId = ccs.getSegment(centroidXyz);
    pointOnSegment = centerline[segmentId];
    /*
    trunkAxis += ccs.getDirectionVector(segmentId);
    //sum of 3 vector
    if(segmentId - 1 >= 0){
        trunkAxis += ccs.getSegment(segmentId - 1);
    }
    if(segmentId + 1 < centerline.size()){
        trunkAxis += ccs.getSegment(segmentId + 1);
    }
    */
    trunkAxis += centerline[segmentId + 1] - centerline[segmentId];
    if(segmentId - 1 >= 0){
        trunkAxis += centerline[segmentId] - centerline[segmentId - 1];
    }
    if(segmentId + 1 < centerline.size() - 1){
        trunkAxis += centerline[segmentId + 2] - centerline[segmentId + 1];
    }

    trunkAxis = trunkAxis.getNormalized();
}



double
FeatureExtractionXyz::distanceBetween2Lines(Z3i::RealPoint p1, Z3i::RealPoint d1, Z3i::RealPoint p2, Z3i::RealPoint d2){
    // n the normal vector of the plan paralel to 2 lines
    Z3i::RealPoint n = d1.crossProduct(d2);
    n = n.getNormalized();

    //vector p1p2
    Z3i::RealPoint p1p2 = p2 - p1;

    //distance = project of p1p2 to n
    return std::abs(p1p2.dot(n));
}

bool FeatureExtractionXyz::testSkewness(Z3i::RealPoint p1, Z3i::RealPoint p2, Z3i::RealPoint p3, Z3i::RealPoint p4){

//trace.info()<<p1<<std::endl;
//trace.info()<<p2<<std::endl;
//trace.info()<<p3<<std::endl;
//trace.info()<<p4<<std::endl;
    //V = 1/6 * (p1 - p4). ((p2 - p4)x(p3-p4))
    Z3i::RealPoint p14 = p1 - p4;
    Z3i::RealPoint p24 = p2 - p4;
    Z3i::RealPoint p34 = p3 - p4;
    //volume 
    double V = 1.0/6 * std::abs(p14.dot(p24.crossProduct(p34)));
//    std::cout<<p14<<std::endl;
//    std::cout<<p24<<std::endl;
//    std::cout<<p34<<std::endl;
//    std::cout<<p24.crossProduct(p34)<<std::endl;
//    std::cout<<"V"<<V<<std::endl;
    return V > 0; 
}
