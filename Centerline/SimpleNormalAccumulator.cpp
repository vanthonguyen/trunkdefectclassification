#include <DGtal/math/Statistic.h>

#include "SimpleNormalAccumulator.h"
//#include "AccumulatorHelper.h"

#include "DGtal/io/writers/VolWriter.h"

#if defined USE_PCL
#include <pcl/features/integral_image_normal.h>
#include <pcl/features/normal_3d.h>
#include <pcl/point_types.h>
#endif


using namespace DGtal;






void
SimpleNormalAccumulator::computeAccumulation(bool retainVertexAsso, bool verbose){
    assert(myIsInitialized);
    if(verbose)
        DGtal::trace.progressBar(0, directions.size());
    trace.info()<<"domain@:"<<myDomain<<std::endl;
    for(auto &v: myDomain){
        myImageMinRadius.setValue(v, myRadius);
        myImageMaxRadius.setValue(v, 0);
    }

    for (unsigned int numNorm = 0; numNorm < directions.size(); numNorm++){
        //    if(verbose)
        DGtal::trace.progressBar(numNorm, directions.size());
        Z3i::RealPoint aNormal  = directions[numNorm];    
        Z3i::RealPoint centerPoint = points[numNorm];
        Z3i::RealPoint currentPoint = centerPoint;
        Z3i::RealPoint previousPoint;

        while((currentPoint - centerPoint).norm()<myRadius){
            if(myDomain.isInside(currentPoint) && previousPoint != currentPoint){
                // retain for each voxel its associated faces contributing to the accumulation
                if(retainVertexAsso){
                    std::vector<unsigned int> v = myAssociationImage(currentPoint);
                    v.push_back(numNorm);        
                    myAssociationImage.setValue(currentPoint, v);
                }
                myAccumulationImage.setValue(currentPoint, myAccumulationImage(currentPoint)+1);
                myRadiusImageAcc.setValue(currentPoint, myRadiusImageAcc(currentPoint) + (currentPoint - centerPoint).norm());
                previousPoint = currentPoint;
                double dist = (currentPoint - centerPoint).norm(); 
                if( dist < myImageMinRadius(currentPoint) ){
                    myImageMinRadius.setValue(currentPoint, dist);
                }

                if( dist > myImageMaxRadius(currentPoint) ) {
                    myImageMaxRadius.setValue(currentPoint, dist);
                }
                // update the max of accumulation
                if( myAccumulationImage(currentPoint) > myMaxAccumulation ) {
                    myMaxAccumulation = myAccumulationImage(currentPoint);
                    myMaxAccumulationPoint = currentPoint;
                }
            }
            previousPoint = currentPoint;
            currentPoint += aNormal;
        }
    }

    //other direction
    for (unsigned int numNorm = 0; numNorm< directions.size(); numNorm++){
        //    if(verbose)
        //      DGtal::trace.progressBar(numNorm, directions.size());
        Z3i::RealPoint aNormal = -directions[numNorm];    
        Z3i::RealPoint centerPoint = points[numNorm];
        Z3i::RealPoint currentPoint = centerPoint;
        Z3i::RealPoint previousPoint;

        while((currentPoint - centerPoint).norm()<myRadius){
            if(myDomain.isInside(currentPoint) && previousPoint != currentPoint){
                // retain for each voxel its associated faces contributing to the accumulation
                if(retainVertexAsso){
                    std::vector<unsigned int> v = myAssociationImage(currentPoint);
                    v.push_back(numNorm);        
                    myAssociationImage.setValue(currentPoint, v);
                }
                myAccumulationImage.setValue(currentPoint, myAccumulationImage(currentPoint)+1);
                myRadiusImageAcc.setValue(currentPoint, myRadiusImageAcc(currentPoint) + (currentPoint - centerPoint).norm());
                previousPoint = currentPoint;
                double dist = (currentPoint - centerPoint).norm(); 
                if( dist < myImageMinRadius(currentPoint) ){
                    myImageMinRadius.setValue(currentPoint, dist);
                }

                if( dist > myImageMaxRadius(currentPoint) ){
                    myImageMaxRadius.setValue(currentPoint, dist);
                }
                // update the max of accumulation
                if( myAccumulationImage(currentPoint)>myMaxAccumulation ) {
                    myMaxAccumulation = myAccumulationImage(currentPoint);
                    myMaxAccumulationPoint = currentPoint;
                }
            }
            previousPoint = currentPoint;
            currentPoint += aNormal;
        }
  }

  for(auto &v: myDomain){
    if ( myAccumulationImage(v) > 1 )
      myRadiusImageAcc.setValue(v, myRadiusImageAcc(v)/(double)myAccumulationImage(v));
  }

//  if(verbose)
//    DGtal::trace.progressBar(directions.size(), directions.size());
  
  myIsAccumulationComputed = true;
}





void
SimpleNormalAccumulator::computeRadiusFromOrigins() {
    // Step 1: ensure that accumulation is computed.
    assert(myIsAccumulationComputed);

    // Step 2:  Recover for each voxel the list of point for which the normal contribute and compute radius from stats
    myMaxRadius = 0;
    for(auto &v: myDomain){
        myRadiusImage.setValue(v, 0);
    }
    double r;
    for(auto &v: myDomain){
        if(myAssociationImage(v).size() >= 2){
            DGtal::Statistic<double> stat(true);
            for(unsigned int ptId :myAssociationImage(v)){
                stat.addValue((points[ptId] -v).norm());
            }
            stat.terminate();
            switch (myRadiusStatEstim)
            {
                case StatRadiusDef::min:
                    r = stat.min();
                    break;
                case StatRadiusDef::mean:
                    r = stat.mean();  
                    break;
                case  StatRadiusDef::max:
                    r = stat.max();  
                    break;
                case  StatRadiusDef::median:
                    r = stat.median();  
                    break;
            }

            myRadiusImage.setValue(v, r);
            if ( r > myMaxRadius ) 
            {
                myMaxRadiusPoint = v; 
                myMaxRadius = r;
            }
        }
    }
    myIsRadiusComputed = true;
}


void
SimpleNormalAccumulator::computeRadiusFromConfidence(){
    if(!myIsAccumulationComputed){
        computeAccumulation();
    }
    if(!myIsAssociationCompFromConfidence){
        computeConfidence(true);
    }
    computeRadiusFromOrigins();
}




SimpleNormalAccumulator::PointContainer
SimpleNormalAccumulator::getAssociatedPoints(const DGtal::Z3i::Point &aVoxel){
    PointContainer result; 
    for(unsigned int ptId :myAssociationImage(aVoxel)){
        result.push_back(points[ptId]);
    }
    return result;
}




void
SimpleNormalAccumulator::computeConfidence(bool updateVertexAsso){
    // Step 1: ensure that accumulation is computed and clean association image if needed.
    assert(myIsAccumulationComputed);
    if (updateVertexAsso){
        for (auto &p: myDomain){
            myAssociationImage.setValue(p, std::vector<unsigned int>());
        }
    }  

    // Step 2:  Apply second face scan to add 1 to the voxel which is the maximal value of the current scan.
    // stored in scoreConfidance.  
    unsigned int maxAcc = 0;
    unsigned int maxAccVerified = 0;
    Image3DDouble scoreConfidance (myDomain);
    for (unsigned int numNorm = 0; numNorm< directions.size(); numNorm++){
        Z3i::RealPoint aNormal  = directions[numNorm];    
        DGtal::Z3i::RealPoint centerPoint = points[numNorm];
        DGtal::Z3i::RealPoint currentPoint = centerPoint;
        DGtal::Z3i::RealPoint previousPoint;

        DGtal::Z3i::RealPoint aPtMaxAcc = currentPoint;
        unsigned int aMaxAcc = 0;

        while((currentPoint - centerPoint).norm()<myRadius){
            if(myDomain.isInside(currentPoint) && previousPoint != currentPoint){        
                unsigned int valAcc = myAccumulationImage(currentPoint);
                if( valAcc > aMaxAcc){
                    aMaxAcc = valAcc;
                    aPtMaxAcc = currentPoint;
                }        
                previousPoint = currentPoint;
            }
            previousPoint = currentPoint;
            currentPoint += aNormal;
        }
        //other direction
        aNormal  = -directions[numNorm];
        centerPoint = points[numNorm];
        currentPoint = centerPoint;
        while((currentPoint - centerPoint).norm()<myRadius){
            if(myDomain.isInside(currentPoint) && previousPoint != currentPoint){        
                unsigned int valAcc = myAccumulationImage(currentPoint);
                if( valAcc > aMaxAcc){
                    aMaxAcc = valAcc;
                    aPtMaxAcc = currentPoint;
                }        
                previousPoint = currentPoint;
            }
            previousPoint = currentPoint;
            currentPoint += aNormal;
        }

        if( aMaxAcc > maxAcc ){
            maxAcc = aMaxAcc;
        }

        if(aPtMaxAcc != centerPoint){
            scoreConfidance.setValue(aPtMaxAcc, scoreConfidance(aPtMaxAcc) + 1);
            if (scoreConfidance(aPtMaxAcc) > maxAccVerified) {
                maxAccVerified = scoreConfidance(aPtMaxAcc);
            }
            if (updateVertexAsso){
                std::vector<unsigned int> v = myAssociationImage(aPtMaxAcc);
                v.push_back(numNorm);        
                myAssociationImage.setValue(aPtMaxAcc, v);
            }
        }

    }
std::cout<<"maxAccVerified#"<< maxAccVerified<<std::endl;
std::cout<<"maxAcc#        "<< maxAcc<<std::endl;
    // Step 3: Compute confidance image indicating the rate between accIsMax/acc
    double minScore = 1;
    double maxScore = 0;
    for(auto &v: myDomain){
        if ( myAccumulationImage(v) > 1 ){
            //double scoreValue = scoreConfidance(v)/maxAccVerified*(double)myAccumulationImage(v)/maxAcc;
            double scoreValue = scoreConfidance(v)/(double)myAccumulationImage(v);
            //double scoreValue = (double)myAccumulationImage(v)/maxAcc;
            //double scoreValue = (double)scoreConfidance(v)/maxAccVerified;
            //double scoreValue = (double)myAccumulationImage(v)/maxAcc;
//std::cout<<"xxx:"<<scoreValue<<std::endl;
            myConfidenceImage.setValue(v, scoreValue);
            if(scoreValue > maxScore){
                maxScore = scoreValue;
            }
            if(scoreValue < minScore){
                minScore = scoreValue;
            }
            //myConfidenceImage.setValue(v, scoreConfidance(v)/(double)myAccumulationImage(v) * myImageMinRadius(v) / myImageMaxRadius(v));
            //myConfidenceImage.setValue(v, myImageMinRadius(v) / myImageMaxRadius(v));
        }
    }
    //normalized

    double range = maxScore - minScore;
    for(auto &v: myDomain){
        if ( myAccumulationImage(v) > 1 ){
            myConfidenceImage.setValue(v, ((double)myConfidenceImage(v) - minScore)/range);
        }
    }
    myIsAssociationCompFromConfidence = updateVertexAsso;
    myIsConfidenceComputed = true;
}



unsigned int 
SimpleNormalAccumulator::getMaxAccumulation() const {
    return myMaxAccumulation;
}

double
SimpleNormalAccumulator::getMaxRadius() const {
  return myMaxRadius;
}

DGtal::Z3i::Point
SimpleNormalAccumulator::getMaxRadiusPoint() const {
  return myMaxRadiusPoint;
}


DGtal::Z3i::Point
SimpleNormalAccumulator::getMaxAccumulationPoint() const{
  return myMaxAccumulationPoint;
}


SimpleNormalAccumulator::Image3D &
SimpleNormalAccumulator::getAccumulationImage() {
  return myAccumulationImage;
}


SimpleNormalAccumulator::Image3DDouble & 
SimpleNormalAccumulator::getConfidenceImage(){
  return myConfidenceImage;
}


SimpleNormalAccumulator::Image3DDouble & 
SimpleNormalAccumulator::getRadiusImage(){
  return myRadiusImage;
}

SimpleNormalAccumulator::Image3DDouble & 
SimpleNormalAccumulator::getRadiusImageAcc(){
  return myRadiusImageAcc;
}

const 
SimpleNormalAccumulator::VectorContainer& 
SimpleNormalAccumulator::getNormalField () const{
  return directions;
}

const 
SimpleNormalAccumulator::VectorContainer& 
SimpleNormalAccumulator::getNormalOrigins () const{
  return points;
}

Z3i::Domain
SimpleNormalAccumulator::getDomain() const
{
  return myDomain;
}

/**
 * Overloads 'operator<<' for displaying objects of class 'NormalAccumulator'.
 * @param out the output stream where the object is written.
 * @param aColor the object of class 'NormalAccumulator' to write.
 * @return the output stream after the writing.
 */
void
SimpleNormalAccumulator::selfDisplay( std::ostream & out ) const {
  out << "----" << std::endl;
  DGtal::Z3i::Point dim = myDomain.upperBound()-myDomain.lowerBound();
  out << "NormalAccumulator: \ninitialized with:\n " << " \t# normals: " << directions.size()
      << "\n\t# normal base: " << points.size() << std::endl;    
  out << "Domain: \n\t" << "size:" << myDomain.size() << " (" << dim[0] << " " << dim[1] << " " << dim[2] << ")\n\t";
  out << "bounds: "<< myDomain.lowerBound() << " " << myDomain.upperBound() << std::endl;
  out << "\tmax accumulation: " << myMaxAccumulation  << " (" << myMaxAccumulationPoint << ")" << std::endl;
  out << "\t radius estim stat: " << myRadiusStatEstim   << std::endl;
  out << "----" << std::endl;
}


std::ostream&
operator<< ( std::ostream & out,
             const SimpleNormalAccumulator & aNormalAcc )
{
    aNormalAcc.selfDisplay ( out );
    return out;
}

