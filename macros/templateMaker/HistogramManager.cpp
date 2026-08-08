#include "HistogramManager.h"
#include <iostream>

HistogramManager::HistogramManager(const HistogramConfig& cfg, std::string Outfilelocation)
{
    // Disable ROOT global histogram ownership
    TH1::AddDirectory(false);

    fDialNames       = cfg.dialNames;
    fSampleNames     = cfg.sampleNames;
    fVariationValues = cfg.variationValues;
    fBinEdges        = cfg.binEdges;

    fNdials      = fDialNames.size();
    fNvariations = fVariationValues.size();
    fNsamples    = fSampleNames.size();

    fFiles.resize(fNdials);
    fDirectories.resize(fNdials);
    fHistograms.resize(fNdials);

    for(int ip = 0; ip < fNdials; ++ip)
    {
        std::string fileName = fDialNames[ip] + "_var.root";

        fFiles[ip] = new TFile((Outfilelocation+fileName).c_str(), "RECREATE");

        fDirectories[ip].resize(fNvariations);
        fHistograms[ip].resize(fNvariations);

        for(int ivar = 0; ivar < fNvariations; ++ivar)
        {
            std::string dirName = Form("variation_%g", fVariationValues[ivar]);

            fDirectories[ip][ivar] = fFiles[ip]->mkdir(dirName.c_str());

            fDirectories[ip][ivar]->cd();

            fHistograms[ip][ivar].resize(fNsamples);

            for(int isample = 0; isample < fNsamples; ++isample)
            {
                std::string histName = Form("h_%s", fSampleNames[isample].c_str());

                fHistograms[ip][ivar][isample] = new TH1D(histName.c_str(), histName.c_str(), fBinEdges.size() - 1, fBinEdges.data());
            }
        }
    }
}

void HistogramManager::Fill(int dialIndex, int variationIndex, int sampleIndex, double value, double weight)
{
    if(dialIndex >= fNdials || variationIndex >= fNvariations || sampleIndex >= fNsamples)
    {
        std::cerr << "HistogramManager::Fill index out of bounds.\n";
        return;
    }

    fHistograms[dialIndex][variationIndex][sampleIndex]->Fill(value, weight);

    // // // // // Overflow handling // // // // //

    int lastBin = fHistograms[dialIndex][variationIndex][sampleIndex]->GetNbinsX();

    double lastBinContent = fHistograms[dialIndex][variationIndex][sampleIndex]->GetBinContent(lastBin);
    double overflowBinContent = fHistograms[dialIndex][variationIndex][sampleIndex]->GetBinContent(lastBin+1);

    fHistograms[dialIndex][variationIndex][sampleIndex]->SetBinContent(lastBin, lastBinContent+overflowBinContent);
    fHistograms[dialIndex][variationIndex][sampleIndex]->SetBinContent(lastBin+1, 0.0);
}

void HistogramManager::WriteAndClose()
{
    for(int ip = 0; ip < fNdials; ++ip)
    {
        fFiles[ip]->cd();

        for(int ivar = 0; ivar < fNvariations; ++ivar)
        {
            fDirectories[ip][ivar]->cd();

            for(int isample = 0; isample < fNsamples; ++isample)
            {
                fHistograms[ip][ivar][isample]->Write();
            }
        }

        fFiles[ip]->Write();
        fFiles[ip]->Close();
    }
}

HistogramManager::~HistogramManager()
{
    // ROOT handles cleanup after file close
}
