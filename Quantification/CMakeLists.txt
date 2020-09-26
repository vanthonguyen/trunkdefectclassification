PROJECT(SubQuantification)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_CXX_FLAGS "-std=gnu++11 -pthread -O2")

FIND_PACKAGE(Threads)

FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
FIND_PROGRAM(CMAKE_PKGCONFIG_EXECUTABLE pkg-config)

if(CMAKE_PKGCONFIG_EXECUTABLE)
  find_package( PkgConfig )    
  pkg_check_modules( EIGEN3 REQUIRED eigen3 )
ELSE(CMAKE_PKGCONFIG_EXECUTABLE)
        FIND_PACKAGE(EIGEN3 REQUIRED)  
ENDIF(CMAKE_PKGCONFIG_EXECUTABLE)

include_directories( ${EIGEN3_INCLUDE_DIRS} )


FIND_PACKAGE(GSL REQUIRED) 
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
LINK_DIRECTORIES(${GSL_LIBRARIES})

find_package(PCL 1.7 REQUIRED COMPONENTS common segmentation features)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
      SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND) 

SET(PCLLib ${PCL_SEGMENTATION_LIBRARIES} ${PCL_FEATURES_LIBRARIES})


FIND_PACKAGE ( OpenCV REQUIRED )


MESSAGE( STATUS "DGTAL_INCLUDE_DIRS:         " ${DGTAL_INCLUDE_DIRS} )
MESSAGE( STATUS "DGTAL_LIB_DIRS:         " ${DGTAL_LIBRARY_DIRS} )


ADD_EXECUTABLE(trunkQuan TrunkQuantification DefectQuantification
    ../Common/IOHelper ../Centerline/SimpleNormalAccumulator 
    ../Common/CylindricalCoordinateSystem ../Common/Statistic 
    ../Centerline/Centerline ../Centerline/CenterlineMultiCore
    ../Segmentation/DefectBranchSegmentation ../Segmentation/SegmentationAbstract ../Segmentation/DefectSegmentation 
    ../Feature/FeatureExtractionXyz ../Feature/PredictRF ../Feature/MomentCylindrical)
TARGET_LINK_LIBRARIES(trunkQuan ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}  ${OpenCV_LIBS} )

ADD_EXECUTABLE(characterizeBranch characterizeBranch 
    ../Common/IOHelper ../Centerline/SimpleNormalAccumulator 
    ../Common/CylindricalCoordinateSystem ../Common/Statistic 
    ../Centerline/Centerline ../Centerline/CenterlineMultiCore
    )
TARGET_LINK_LIBRARIES(characterizeBranch ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}  ${OpenCV_LIBS} )
#scar
ADD_EXECUTABLE(characterizeScar characterizeScar
    ../Common/IOHelper
    ../Common/CylindricalCoordinateSystem
    ../Feature/MomentCylindrical ../Common/Statistic 
    )
TARGET_LINK_LIBRARIES(characterizeScar ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}  ${OpenCV_LIBS} )


