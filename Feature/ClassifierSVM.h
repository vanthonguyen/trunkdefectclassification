#ifndef CLASSIFIER_SVM_H
#define CLASSIFIER_SVM_H
#include <vector>
#include <iostream>

#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/ml.hpp"
#include "opencv2/highgui.hpp"





#define debug(str)\
    std::cout<<str<<std::endl;
class ClassifierSVM{
    public:
        ClassifierSVM(){
            svmAvail = false;
        }
        //train(int dictionarySize);
        void train(const std::vector<std::vector<float> > &fs, const std::vector<int> &lbls, double c = 300.0, double gamma = 0.5);
        void trainAuto(const std::vector<std::vector<float> > &fs, const std::vector<int> &lbls, int kernelType);
        float predict(const std::vector<float> &feature);
        std::vector<float> predict(const std::vector<std::vector<float>> &features);

        static const std::string SVM_FILE;
    private:
        cv::Ptr<cv::ml::TrainData> createTrainData(const std::vector<std::vector<float>> &features, const std::vector<int> &labels);

        void loadFromFile();

        bool svmAvail;
        cv::Ptr<cv::ml::SVM> svm;
};
#endif
