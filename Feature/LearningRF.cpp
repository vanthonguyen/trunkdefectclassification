#include "LearningRF.h"


const std::string LearningRF::RF_FILE = "learnedForest";

void 
LearningRF::learn(){

    cv::Ptr<cv::ml::RTrees> rtrees = cv::ml::RTrees::create();

    cv::Ptr<cv::ml::TrainData> data = createTrainData();

    rtrees->setMaxDepth(256);
    rtrees->setMinSampleCount(2);
    //rtrees->setRegressionAccuracy(0.f);
    rtrees->setUseSurrogates(false);
    rtrees->setMaxCategories(4);
    rtrees->setPriors(cv::Mat());
    rtrees->setCalculateVarImportance(true);
    rtrees->setActiveVarCount(nbVar);
    //rtrees->setActiveVarCount(2);
    rtrees->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS, nbTree, oobErr));
    rtrees->train(data);
    //save trees for predict
    rtrees->save(RF_FILE);

    printf("test with ooErr: %f\n", oobErr);
    std::cout<<rtrees->getVarImportance()<<std::endl;
    //std::cout<<"calc error: "<< rtrees->calc_error()<<std::endl;
    //std::cout<<"train error: "<< rtrees->get_train_error()<<std::endl;
    std::cout<<"Active Var Count:"<<rtrees->getActiveVarCount()<<std::endl;

    printf( "train error: %f\n", rtrees->calcError(data, false, cv::noArray()) );
    printf( "test error: %f\n\n", rtrees->calcError(data, true, cv::noArray()) );
}


cv::Ptr<cv::ml::TrainData> 
LearningRF::createTrainData(){
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
