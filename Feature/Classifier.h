#ifndef LEARNING_H
#define LEARNING_H

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
class Learning{

protected:
    //virtual std::string getFile();
    cv::Ptr<cv::ml::TrainData> createTrainData(){
        int nbSample = features.size();
        int nbAttribute = features.at(0).size();
        cv::Mat samples(nbSample, nbAttribute, CV_32FC1);
        for( int featureId = 0; featureId < features.size(); featureId++ ){
            std::vector<float> feature = features.at(featureId);
            for(int attrId = 0; attrId < feature.size(); attrId++){
                samples.at<float>(featureId, attrId) = feature.at(attrId);
            }

        }
        return cv::ml::TrainData::create(samples, cv::ml::ROW_SAMPLE, cv::Mat(labels));
    }
    std::vector<std::vector<float> > features;
    std::vector<int> labels;
};
#endif //Learning
