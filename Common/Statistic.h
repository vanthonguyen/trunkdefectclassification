#ifndef STATISTIC_H
#define STATISTIC_H

#include<vector>

class Statistic
{
public:
    Statistic();
    static double standardDeviation(const std::vector<double> &v, double mean);
    static double getMean(const std::vector<double> &v);
    static double getMode(const std::vector<double> &v, double min, double max, double binWidth);
    static double getMode(const std::vector<double> &v, double binWidth);
//    static double getMedian(std::vector<double> v);
};

#endif // STATISTIC_H
