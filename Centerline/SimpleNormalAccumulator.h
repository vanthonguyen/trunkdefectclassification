#ifndef SIMPLE_NORMAL_ACCUMULATOR_H
#define SIMPLE_NORMAL_ACCUMULATOR_H



#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Mesh.h"
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/images/ImageContainerBySTLVector.h>


///////////////////////////////////////////////////////////////////////////////
// class AccumulateNormal
/**
 * Description of class 'AcumulateNormal' <p>
 * 
 * @brief Class to compute volumic accumulation from normal vector field.
 *
 *
 */


using namespace DGtal;

class SimpleNormalAccumulator{


  // types of main object
public:
  enum StatRadiusDef {mean, median, min, max};
  typedef typename DGtal::Z3i::RealPoint Point; 
  typedef typename DGtal::Z3i::RealPoint Vector; 
  typedef std::vector<Point> PointContainer;
  typedef std::vector<Vector> VectorContainer;  
  
  // types of image containers:
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, DGtal::uint64_t> Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, double> Image3DDouble;
  // to recover the origin point which has contributed to a particular voxel.
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  std::vector<unsigned int> > ImagePointAssociation;

  

  // ----------------------- Standard services ------------------------------
  
public:
  
  SimpleNormalAccumulator(  const double aRadius,
                            const Z3i::Domain aDomain,
                            const std::vector<Z3i::RealPoint> &pointsNormal,
                            const std::vector<Z3i::RealPoint> &startingPoints
                         ): myRadius(aRadius), 
                            directions(pointsNormal),
                            points(startingPoints),
                            myDomain(aDomain),
                            myAccumulationImage(aDomain),
                            myConfidenceImage(aDomain),
                            myRadiusImage(aDomain),
                            myRadiusImageAcc(aDomain),
                            myImageMinRadius(aDomain),
                            myImageMaxRadius(aDomain),
                            myAssociationImage(aDomain){

    myIsAccumulationComputed = false;
    myIsRadiusComputed = false;
    myIsConfidenceComputed = false;
    myIsAssociationCompFromConfidence = false;
    myMaxAccumulation = 0;
    myMaxRadius = 0;    
    myIsInitialized = true;
}


  
  /**
   * @todo import normal computation from DGtal surface
   * 
   **/
  //void initFromDGtalSurface(const DGtal::Surface &aSurface);



  /**
   * Compute accumulation  from normal vectors and update the maximum accumulation value.
   *
   **/
  
  void computeAccumulation(bool retainVertexAsso = true, bool verbose = true);

  
  /**
   * Compute the confidence image
   * 
   * @param updateVertexAsso if true update the image association with only the confident vertex. 
   * 
   **/
  void computeConfidence(bool updateVertexAsso = true);
  
  
  /**
   * Compute the radius from all contributing normal origins
   *
   **/
  void computeRadiusFromOrigins();



  /**
   * Compute the radius from all contributing normal with confidance rate.
   *
   **/
  void computeRadiusFromConfidence();
  


  /**
   * Get the associated input points associated to a voxel. It can be
   * all the points with normal contributing to the accumulation or
   * the confidence (if the point association is updated whe when
   * calling computeConfidence and passing true on argument
   * updateVertexAsso).
   *
   * @param aVoxel the input voxel
   * @return the container with all point contributing to the condidence.
   **/

  PointContainer getAssociatedPoints(const DGtal::Z3i::Point &aVoxel); 
  
  
  /**
   * return the maximum accumulation value. 
   *
   **/  
  unsigned int getMaxAccumulation() const;


  /**
   * return the maximum radius value. 
   *
   **/  
  double getMaxRadius() const;


  /**
   * return the maximum radius point. 
   *
   **/  
  DGtal::Z3i::Point getMaxRadiusPoint() const;
  

  /**
   * Return the point of maximum accumulation value.
   *
   **/  
  DGtal::Z3i::Point getMaxAccumulationPoint() const;

  const VectorContainer& getNormalField ()const;
 
  const VectorContainer& getNormalOrigins ()const;
  
  
  Image3D & getAccumulationImage() ;
  

  Image3DDouble & getConfidenceImage() ;


  Image3DDouble & getRadiusImageAcc() ;


  Image3DDouble & getRadiusImage() ;

  Z3i::Domain getDomain() const;
  
  /**
   * self display method.
   **/
  void  selfDisplay( std::ostream & out) const;

  
  
  // protected attributes:
protected:

  std::vector<Z3i::RealPoint> directions;
  std::vector<Z3i::RealPoint> points;
  double myRadius; // the maximal radius of accumulation
  Z3i::Domain myDomain; // the myDomain of the source objet.
  

  Image3D myAccumulationImage;
  Image3DDouble myConfidenceImage;
  Image3DDouble myRadiusImage;    
  Image3DDouble myRadiusImageAcc;

  Image3DDouble myImageMinRadius;
  Image3DDouble myImageMaxRadius;
  ImagePointAssociation myAssociationImage; // to recover all normal origins contributing to an accumulation point.


private:
  bool myIsVerbose = true;
  bool myIsInitialized = false;
  bool myIsAccumulationComputed = false;
  bool myIsConfidenceComputed = false;
  bool myIsRadiusComputed = false;
  bool myIsAssociationCompFromConfidence = false;
  DGtal::Z3i::Point myMaxAccumulationPoint;
  unsigned int myMaxAccumulation = 0;
  double myMaxRadius = 0;
  DGtal::Z3i::Point myMaxRadiusPoint; 
  StatRadiusDef myRadiusStatEstim = StatRadiusDef::median;  
};



/**
 * Overloads 'operator<<' for displaying objects of class 'NormalAccumulator'.
 * @param out the output stream where the object is written.
 * @param aColor the object of class 'NormalAccumulator' to write.
 * @return the output stream after the writing.
 */
std::ostream&
operator<< ( std::ostream & out, const SimpleNormalAccumulator & aNormalAcc );


#endif //end SIMPLE_NORMAL_ACCUMULATOR

