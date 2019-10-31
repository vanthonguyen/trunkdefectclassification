#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <utility>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "MomentsBin.h"


namespace po = boost::program_options;
template <typename T>
void export2Text(const std::vector<T> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        outStream << v[i]<<std::endl;
    }
    outStream.close();
}

/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input point cloud.")
        ("output,o", po::value<std::string>(), "output bmp file name.");


    bool parseOK = true;
    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }
    catch (const std::exception &ex)
    {
        std::cout << "Error checking program options: " << ex.what() << std::endl;
        parseOK = false;
    }
    po::notify(vm);
    if ( !parseOK || vm.count("help") || argc <= 1 || !vm.count("input") )
    {
        std::cout << "Test:" << std::endl
            << "Options: " << std::endl
            << general_opt << std::endl;
        return 0;
    }

    std::string inputImg = vm["input"].as<std::string>();
    std::string outname = vm["output"].as<std::string>();

    //read a binary image
    cv::Mat srcImg = cv::imread(inputImg);


	MomentsBin mo(inputImg);
	cv::Moments mu = mo.getMoment();
	std::vector<double> hu = mo.getHuMoments();
	std::vector<double> nu;
	//nu20, nu11, nu02, nu30, nu21, nu12, nu03
	nu.push_back(mu.nu20);
	nu.push_back(mu.nu11);
    nu.push_back(mu.nu02); 
	nu.push_back(mu.nu30); 
	nu.push_back(mu.nu21); 
	nu.push_back(mu.nu12); 
	nu.push_back(mu.nu03);


	for(double x : nu){
		std::cout<< x << "  ";
	}
	std::cout<<std::endl;

	for(double x : hu){
		std::cout<< x << "  ";
	}
	std::cout<<std::endl;

    export2Text(nu, outname);
    export2Text(hu, outname + "hu");

    return 0;
}
