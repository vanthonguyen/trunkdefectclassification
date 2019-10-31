#include "ClassifierSVM.h"

const std::string ClassifierSVM::SVM_FILE = "svmFile";

void 
ClassifierSVM::trainAuto(const std::vector<std::vector<float> > &fs, const std::vector<int> &lbls, int kernelType){
    cv::Ptr<cv::ml::SVM> svm = cv::ml::SVM::create();
    svm->setType(cv::ml::SVM::C_SVC);
    switch (kernelType){
        case 1:
            svm->setKernel(cv::ml::SVM::LINEAR); //SVM::LINEAR;
            break;
        case 2:
            svm->setKernel(cv::ml::SVM::RBF); //SVM::LINEAR;
            break;
        case 3:
            svm->setKernel(cv::ml::SVM::CHI2); //SVM::LINEAR;
            break;
        case 4:
            svm->setKernel(cv::ml::SVM::INTER); //SVM::LINEAR;
            break;
        default:
            svm->setKernel(cv::ml::SVM::POLY); //SVM::LINEAR;
            break;
    }

    //svm->setKernel(cv::ml::SVM::INTER); //SVM::LINEAR;
    svm->trainAuto(createTrainData(fs, lbls)); 
    svm->save(SVM_FILE);
    predict(fs[0]);
    svmAvail = true;
}


void 
ClassifierSVM::train(const std::vector<std::vector<float> > &fs, const std::vector<int> &lbls, double c, double gamma){
    cv::Ptr<cv::ml::SVM> svm = cv::ml::SVM::create();
    svm->setType(cv::ml::SVM::C_SVC);
    svm->setKernel(cv::ml::SVM::RBF); //SVM::LINEAR;
    //svm->setKernel(cv::ml::SVM::INTER); //SVM::LINEAR;
    svm->setDegree(0.5);
    svm->setGamma(gamma);
    svm->setCoef0(1);
    svm->setNu(0.5);
    svm->setP(0);
    svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER+cv::TermCriteria::EPS, 1000, 0.0001));
    svm->setC(c);
    svm->train(createTrainData(fs, lbls)); 
    svm->save(SVM_FILE);
    predict(fs[0]);
    svmAvail = true;
}

float
ClassifierSVM::predict(const std::vector<float> &feature){
    if(!svmAvail){
        loadFromFile();
    }
    return svm->predict( feature );
}
std::vector<float> 
ClassifierSVM::predict(const std::vector<std::vector<float>> &features){
    std::vector<float> predictLabels;
    if(!svmAvail){
        loadFromFile();
    }

    for(std::vector<float> feature: features){

        //debug(feature[0]);
        float response = svm->predict( feature );
        predictLabels.push_back(response);
    }
    return predictLabels;
}


void
ClassifierSVM::loadFromFile(){
    svm = cv::ml::StatModel::load<cv::ml::SVM> (SVM_FILE);
    svmAvail = true;
}

cv::Ptr<cv::ml::TrainData> 
ClassifierSVM::createTrainData(const std::vector<std::vector<float>> &features, const std::vector<int> &labels){
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

