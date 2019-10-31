#ifndef LEARNING_RF_H
#define LEARNING_RF_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <array>

#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/ml.hpp"
#include "opencv2/highgui.hpp"

/**
 * Input is a list of feature vectors
 * and associate with class
 */
class LearningRF{
public:

    LearningRF(const std::vector<std::vector<float> > &fs, const std::vector<int> &lbls, const int &nt, const int &nv, double oob = 0.01): 
        features(fs), labels(lbls), nbTree(nt), nbVar(nv), oobErr(oob){}
    void learn();

    static const std::string RF_FILE;
private:
    cv::Ptr<cv::ml::TrainData> createTrainData();
    std::vector<std::vector<float> > features;
    std::vector<int> labels;
    double oobErr;
    int nbTree;
    int nbVar;
};
#endif //MOMENTS_3D
