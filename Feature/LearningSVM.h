#ifndef LEARNING_SVM_H
#define LEARNING_SVM_H

#include <iostream>
#include <vector>
#include <array>

#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/ml.hpp"
#include "opencv2/highgui.hpp"

#include "Learning.h"
/**
 * Input is a list of feature vectors
 * and associate with class
 */
class LearningSVM : Learning{
public:
    using Learning::Learning;
    void learn();
};
#endif //MOMENTS_3D
