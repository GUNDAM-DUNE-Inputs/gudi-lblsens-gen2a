#ifndef HISTOGRAMMANAGER_H
#define HISTOGRAMMANAGER_H

#include <vector>
#include <string>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include "Config.h"

class HistogramManager
{
public:
    HistogramManager(const HistogramConfig& config, std::string filepath);
    ~HistogramManager();

    void Fill(int dialIndex, int variationIndex, int sampleIndex, double value, double weight);

    void WriteAndClose();

private:
    std::vector<std::string> fDialNames;
    std::vector<std::string> fSampleNames;
    std::vector<double>      fVariationValues;
    std::vector<double>      fBinEdges;

    int fNdials;
    int fNvariations;
    int fNsamples;

    std::vector<TFile*> fFiles;
    std::vector<std::vector<TDirectory*>> fDirectories;
    std::vector<std::vector<std::vector<TH1D*>>> fHistograms;
};

#endif
