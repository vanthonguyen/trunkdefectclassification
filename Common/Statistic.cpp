#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <assert.h>

#include "Statistic.h"
double Statistic::standardDeviation(const std::vector<double> &v, const double mean){
    double s = 0.0;
    for (size_t i = 0; i < v.size(); i++)
    {
        s += (v[i] - mean) * (v[i] - mean);
    }
    return sqrt(s/(v.size() - 1));
}

double Statistic::getMean(const std::vector<double> &v){
    //calculate mean
    double sum = 0.0;
    for (unsigned int i = 0; i < v.size(); i++)
    {
        sum += v[i];
    }

    return sum / v.size();
}

//double Statistic::getMedian(std::vector<double> v){
//    return 0.0;
//}


double Statistic::getMode(const std::vector<double> &v, double min, double max, double binWidth){
    int nbInterval = (max - min) / binWidth + 1;
//std::cout<<"nbIn:"<< nbInterval<<std::endl;
//std::cout<<"min:"<< min<<std::endl;
//std::cout<<"max:"<< max<<std::endl;
    std::vector<int> histogram(nbInterval, 0);
    for(unsigned int i = 0; i < v.size(); i++){
        int index = (v.at(i) - min)/binWidth;
//std::cout<<v.at(i)<<std::endl;
//std::cout<<index<<std::endl;
        assert(index < histogram.size());
        histogram[index]++;
    }
    std::vector<int>::iterator maxFreq = std::max_element(histogram.begin(), histogram.end());
    return min + std::distance(histogram.begin(), maxFreq) * binWidth;
}

double Statistic::getMode(const std::vector<double> &v, double binWidth){
    auto xMinMax = std::minmax_element(v.begin(), v.end());
    double min = *xMinMax.first;
    double max = *xMinMax.second;
    return getMode(v, min, max, binWidth);
}
