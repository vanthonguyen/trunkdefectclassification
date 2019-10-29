#include <iostream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <ctime>



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/ml.hpp"
#include "opencv2/highgui.hpp"



#include "../Common/IOHelper.h"
#include "../Feature/PredictRF.h"

namespace po = boost::program_options;

const std::string RED("\033[0;31m");
const std::string DEFAULT("\033[0m");

#define debug(str)\
    std::cout<<str<<std::endl;
int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input descripteur.")
        ("flag,f", po::value<int>()->default_value(0), "flags.")
        ("rffile,r", po::value<std::string>(), "input random forest file.");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }


    std::string featureFile = vm["input"].as<std::string>();
    std::string rfFile = vm["rffile"].as<std::string>();
    int flags = vm["flag"].as<int>();
    cv::Ptr<cv::ml::RTrees> rtrees = cv::ml::StatModel::load<cv::ml::RTrees> (rfFile);



    std::vector<std::vector<float> > features;
    IOHelper::readFeatures(featureFile, features);

    std::vector<float> feature = features.at(0);

    cv::Mat results; 
    float label = rtrees->predict(feature, results, cv::ml::StatModel::RAW_OUTPUT);

    //std::cout<< label <<std::endl;
    //std::cout<< results <<std::endl;

    cv::Mat votes;
    
    rtrees->getVotes(feature, votes, flags);

    //std::cout<<flags<< " votes"<<std::endl;
    std::cout<< votes <<std::endl;

    PredictRF prf(rfFile);


    int label2 = prf.predict(feature);
    std::cout<< "label2:"<< label2<<std::endl;
    return 0;
}

