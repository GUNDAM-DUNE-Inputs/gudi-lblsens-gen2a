#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>

struct HistogramConfig
{
    std::vector<std::string> dialNames;
    std::vector<std::string> sampleNames;
    std::vector<double> variationValues;
    std::vector<double> binEdges;
};

#endif
