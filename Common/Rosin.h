#ifndef ROSIN_H
#define ROSIN_H

#include <utility>
#include <vector>
#include <cassert>
#include <iostream>

class Rosin
{
public:
    Rosin(){}
    static double compute(const std::vector<double> &data, const double &binWidth){

        //trace.error()<<"binWidth:"<<binWidth<<std::endl;
        double maxValue = *std::max_element(data.begin(), data.end());
        double minValue = *std::min_element(data.begin(), data.end());
        double range = maxValue - minValue;
        int nbInterval = range / binWidth + 1;

        std::vector<int> histogram(nbInterval, 0);
        for(unsigned int i = 0; i < data.size(); i++){
            int index = (data.at(i) - minValue)/binWidth;
            histogram[index]++;
        }
        std::vector<int>::iterator maxFreq = std::max_element(histogram.begin(), histogram.end());

        int maxFreqValue = *maxFreq;
        int maxFreqIndex = std::distance(histogram.begin(), maxFreq);
        assert(maxFreqValue == histogram.at(maxFreqIndex));

        unsigned int lastIndex = histogram.size() - 1;
        int lastValue = histogram.at(lastIndex);

        for(unsigned int i = maxFreqIndex; i < histogram.size(); i++){
            if(histogram.at(i) == 0){
                lastIndex = i;
                lastValue = 0;
                break;
            }
        }

        double valueDiff = lastValue - maxFreqValue;
        double valueDiff2 = valueDiff *valueDiff;
        double indexDiff = lastIndex - maxFreqIndex;
        double indexDiff2 = indexDiff * indexDiff;
        double bestThresIndex = maxFreqIndex;
        double bestDist = 0;

        //line between maxFreq and last element of historgram
        double a = (lastValue - maxFreqValue)*1.0/(lastIndex - maxFreqIndex);
        double b = maxFreqValue - a * maxFreqIndex;

        //for (std::vector<int>::iterator it = maxFreq ; it != histogram.end(); ++it){
        for (unsigned int i = maxFreqIndex; i < lastIndex; i++){
            //trace.info()<<valueDiff *  i - indexDiff*histogram.at(i) + maxFreqIndex*lastValue - maxFreqValue * lastIndex<<std::endl;
            double dist = std::abs(valueDiff *  i - indexDiff*histogram.at(i) + maxFreqValue*lastIndex - maxFreqIndex * lastValue)/
                sqrt(valueDiff2 + indexDiff2 );
            if(dist > bestDist){
                bestDist = dist;
                bestThresIndex = i;
            }
        }

        /**
         * perpendicular line: -1/a + c passe through bestThresIndex
         * -1/ax1 + c = y1
         *  c = y1 + 1/ax1
         *  ax2 + b = y2
         * -1/ax2 + c = y2
         *  -1/ax2 + y1 + 1/ax1 = y2
         *  (a + 1/a)x2 = y1 + 1/ax1 -b
         */
        int bestVal = histogram.at(bestThresIndex);
        double x2 = (bestVal + bestThresIndex/a - b)/(a + 1/a);
        double y2 = b + a*x2;

        std::pair<double, double> maxPoint(maxFreqIndex * binWidth+ minValue, maxFreqValue);
        std::pair<double, double> lastPoint(lastIndex * binWidth+ minValue, lastValue);
        std::pair<double, double> bestPoint(bestThresIndex * binWidth+ minValue, bestVal);
        std::pair<double, double> projBestPoint(x2* binWidth+ minValue, y2);

        //histogram
        /*
        std::vector<std::pair<double, double>> histForPlot;
        for(unsigned int i = 0; i< histogram.size(); i++){
            std::pair<double, double> aBin(i* binWidth+ minValue, histogram.at(i));
            histForPlot.push_back(aBin);
        }
        */

        std::cout<<"threshold: "<< bestThresIndex*binWidth+ minValue<<std::endl;
        return bestThresIndex*binWidth+ minValue;
    }

};
#endif // ROSIN
