#ifndef PREDICT_RF_H
#define PREDICT_RF_H

#include <iostream>
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
class PredictRF{
public:
    PredictRF(const std::string &treeFile): serialRFFile(treeFile){
        rtrees = cv::ml::StatModel::load<cv::ml::RTrees> (serialRFFile);
    }
    std::vector<int> predicts(const std::vector<std::vector<float>> &features);
    int predict(const std::vector<float> &feature);
    int predict(const std::vector<double> &feature);
private:
    std::string serialRFFile;
    cv::Ptr<cv::ml::RTrees> rtrees;
    //std::vector<std::vector<float> > features;
};
#endif //MOMENTS_3D
