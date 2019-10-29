#include <iostream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <ctime>



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "../IOHelper.h"
#include "../Feature/LearningRF.h"
#include "../Feature/PredictRF.h"


using namespace DGtal;
namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;
const std::string RED("\033[0;31m");
const std::string DEFAULT("\033[0m");

int
main(int argc,char **argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input learning features file.")
        ("label,l", po::value<std::string>(), "input learning labels file.")
        ("test,t", po::value<std::string>(), "input test features file.")
        ("labelTest", po::value<std::string>(), "input test labels file, for verification")
        ("defnameFileTest", po::value<std::string>()->default_value("testNames"), "defectNameDefFile")
        ("output,o", po::value<std::string>()->default_value("defect"), "output defect clusters into file defect0.xyz, defect1.xyz, ...");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }


    std::string learningFeatureFile = vm["input"].as<std::string>();
    std::string learningLabelFile = vm["label"].as<std::string>();
    std::string testFeatureFile = vm["test"].as<std::string>();
    std::string testLabelFile = vm["labelTest"].as<std::string>();
    std::string testNameFile = vm["defnameFileTest"].as<std::string>();


    std::vector<std::vector<float> > learningFeatures;
    IOHelper::readFeatures(learningFeatureFile, learningFeatures);
    
    std::vector<int> labels;
    IOHelper::readIntsFromFile(learningLabelFile, labels);

    std::vector<std::vector<float> > testFeatures;
    IOHelper::readFeatures(testFeatureFile, testFeatures);

    std::vector<int> testLabels;
    IOHelper::readIntsFromFile(testLabelFile, testLabels);

    std::vector<std::string> testNames;
    IOHelper::readStringsFromFile(testNameFile, testNames);

    LearningRF lrf(learningFeatures, labels);
    lrf.learn();


    PredictRF prf(LearningRF::RF_FILE);
    std::vector<int> predictResults = prf.predicts(testFeatures);
    assert(testLabels.size() ==  predictResults.size());

    std::srand(std::time(0)); 

    std::cout<<std::setw(24)<< "def"
             <<std::setw(4)<< "pre"
             <<std::setw(4)<< "gt"
             <<std::setw(4)<<"random"<<std::endl;
    int tp = 0;
    int randTp = 0;
    for ( int i = 0; i < testLabels.size(); i ++ ){
        int r = std::rand() % labels.size() + 1; 
        int randomLabel = labels[r];
//        std::cout<< predictResults[i] << "  "<< testLabels[i]<< "  "<< randomLabel<<std::endl;
        if(predictResults[i] != testLabels[i]){
            std::cout<<RED;
        }
        std::cout<<std::setw(24)<< testNames[i]
                 <<std::setw(4)<< predictResults[i]
                 <<std::setw(4)<< testLabels[i]
                 <<std::setw(4)<<randomLabel<<DEFAULT<<std::endl;
        if(predictResults[i] == testLabels[i]){
            tp++;
        }
        if(randomLabel == testLabels[i]){
            randTp++;
        }
    }
    std::cout<< "correct / total = "<< tp << "/" << predictResults.size()<<"="<<tp*1.0/predictResults.size()<<std::endl;
    std::cout<< "random: correct / total = "<< randTp << "/" << predictResults.size()<<"="<<randTp*1.0/predictResults.size()<<std::endl;

    return 0;
}

