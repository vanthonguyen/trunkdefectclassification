#include "PredictRF.h"

std::vector<int> PredictRF::predicts(const std::vector<std::vector<float>> &features){
    std::vector<int> predictLabels;
    for(std::vector<float> feature: features){
        //cv::Ptr<cv::ml::RTrees> rtrees = cv::ml::StatModel::load<cv::ml::RTrees> (serialRFFile);
            //cv::ml::RTrees::load(serialRFFile);

        //Ptr<RTrees> model_read = StatModel::load<RTrees>( filename_model );
        float response = rtrees->predict( feature );
//std::cerr<<"Class:"<< response<<std::endl;
        predictLabels.push_back((int)response);
    }
    return predictLabels;
}



int PredictRF::predict(const std::vector<float> &feature){
    float response = rtrees->predict( feature );
std::cerr<<"Class:"<< response<<std::endl;
    return (int)response;
}

int PredictRF::predict(const std::vector<double> &feature){
    std::vector<float> featuresInFloat;
    featuresInFloat.insert(featuresInFloat.end(),feature.begin(),feature.end());
    return predict(featuresInFloat);
}
