#include <TFile.h>
#include <TDirectoryFile.h>
#include <TKey.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <map>
#include <vector>
#include <string>

// A small struct to hold interaction-mode info
struct ModeInfo {
    int modeId;
    std::string modeName;
    int color;
};

// Build a lookup table for your SIMB_Mode enumerations
std::map<int, ModeInfo> buildModeMap()
{
    std::map<int, ModeInfo> modeMap;
    modeMap[-1]  = { -1, "Unknown",        kGray+2   };
    modeMap[0]   = {  0, "QE",             kOrange   };
    modeMap[1]   = {  1, "Resonant",       kGreen+2  };
    modeMap[2]   = {  2, "DIS",            kBlue     };
    modeMap[3]   = {  3, "Coherent",       kMagenta+1};
    modeMap[4]   = {  4, "CohElastic",     kOrange+3 };
    modeMap[5]   = {  5, "e-Scattering",   kCyan+1   };
    modeMap[6]   = {  6, "IMD",            kSpring+6 };
    modeMap[7]   = {  7, "Inv Beta Decay", kYellow+2 };
    modeMap[8]   = {  8, "GlashowRes",     kViolet+2 };
    modeMap[9]   = {  9, "AMNuGamma",      kAzure+7  };
    modeMap[10]  = { 10, "MEC (2p2h)",     kRed      };
    modeMap[11]  = { 11, "Diffractive",    kPink+2   };
    modeMap[12]  = { 12, "EM",             kYellow+3 };
    modeMap[13]  = { 13, "WeakMix",        kGray     };
    return modeMap;
}

// Helper to get a subdirectory from a parent TDirectoryFile
// Returns nullptr if not found.
TDirectoryFile* GetSubdir(TDirectoryFile* parent, const char* name)
{
    if (!parent) return nullptr;
    return dynamic_cast<TDirectoryFile*>( parent->Get(name) );
}

// A struct to hold the result of building one sample's stack
struct SamplePlot {
    THStack* stack;
    TLegend* legend;
    std::vector<TH1D*> histList; // keep histogram clones alive
};

// Helper function to build a THStack for one "merged sample"
SamplePlot buildStackForSample(const std::string &mergedName,
                               const std::vector<std::string> &subDirs,
                               TDirectoryFile* histDir,
                               const std::map<int, ModeInfo> &modeMap)
{
    SamplePlot result;
    result.stack = nullptr;
    result.legend = nullptr;

    // Decide if we want Ev_reco_numu or Ev_reco_nue
    std::string evRecoName;
    if (mergedName.find("#nu_{#mu}") != std::string::npos) {
        evRecoName = "Ev_reco_numu";
        // evRecoName = "Raw";
    }
    else if (mergedName.find("#nu_{e}") != std::string::npos) {
        evRecoName = "Ev_reco_nue";
        // evRecoName = "Raw";
    }
    else {
        std::cerr << "Cannot determine Ev_reco folder for " << mergedName << std::endl;
        return result; 
    }

    // Create a THStack
    THStack* hs = new THStack(
        Form("stack_%s", mergedName.c_str()),
        Form("%s;Reconstructed E (GeV);Events / bin", mergedName.c_str())
    );

    // Create a TLegend
    TLegend* leg = new TLegend(0.7, 0.5, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // We'll keep a list of hist clones in memory
    std::vector<TH1D*> histClones;

    // Loop over each subdirectory (CC, NC)
    for (auto &subDirName : subDirs) {

        // Retrieve the subdirectory from histDir
        TDirectoryFile* sampleDir = GetSubdir(histDir, subDirName.c_str());
        if (!sampleDir) {
            std::cerr << "Warning: subdirectory not found: " << subDirName << std::endl;
            continue;
        }

        // Just Raw folder
        // std::string evRecoName = "Raw";
        // Then get Ev_reco_numu or Ev_reco_nue
        TDirectoryFile* evRecoDir = GetSubdir(sampleDir, evRecoName.c_str());
        if (!evRecoDir) {
            std::cerr << "Warning: " << evRecoName << " not found in "
                      << subDirName << std::endl;
            continue;
        }

        // Then "mode"
        TDirectoryFile* modeDir = GetSubdir(evRecoDir, "mode");
        if (!modeDir) {
            std::cerr << "Warning: 'mode' not found in "
                      << subDirName << "/" << evRecoName << std::endl;
            continue;
        }

        // Now loop over all possible mode subdirectories
        TIter nextMode(modeDir->GetListOfKeys());
        TKey* keyMode = nullptr;

        while ((keyMode = (TKey*)nextMode())) {
            if (std::string(keyMode->GetClassName()) != "TDirectoryFile") continue;

            // e.g. "0" => modeVal = 0, "2" => modeVal = 2, etc.
            const char* modeDirName = keyMode->GetName();
            int modeVal = atoi(modeDirName);

            TDirectoryFile* thisModeDir = GetSubdir(modeDir, modeDirName);
            if (!thisModeDir) continue;

            // Retrieve the MC_TH1D
            TH1D* h = dynamic_cast<TH1D*>( thisModeDir->Get("MC_TH1D") );
            if (!h) continue;

            // Clone so we can set style
            TH1D* cloneHist = (TH1D*)h->Clone(
                Form("clone_%s_%s_mode%s", mergedName.c_str(),
                     subDirName.c_str(), modeDirName)
            );

            // Base color from the mode
            auto it = modeMap.find(modeVal);
            int baseColor = kBlack;
            std::string baseName = Form("Mode%d", modeVal);
            if (it != modeMap.end()) {
                baseColor = it->second.color;
                baseName = it->second.modeName;
            }

            // Distinguish CC vs NC by offset color
            bool isNC = (subDirName.rfind(" NC ") != std::string::npos);
            int finalColor = baseColor;
            if (isNC) {
                finalColor += 2; // offset to visually differentiate NC
            }

            cloneHist->SetFillColor(finalColor);
            cloneHist->SetLineColor(finalColor);

            // Build legend label like "CC DIS" or "NC QE"
            std::string label;
            if (!isNC) {
                // subDirName ends in _CC
                label = "CC " + baseName;
            } else {
                // subDirName ends in _NC
                label = "NC " + baseName;
            }

            // ------------- PRINT THE INTEGRAL -------------
            // integrate only over 0–10 GeV (x‑axis)
            // int binLow  = h->GetXaxis()->FindBin(0.0);
            // int binHigh = h->GetXaxis()->FindBin(10.0 - 1e-6); // upper edge just below 10
            // double integral_0_10 = h->Integral(binLow, binHigh);
            double integral = h->Integral("width");

            std::cout << "Directory: " << subDirName
                      << " , Mode: " << modeVal
                      << ", " << label
                      << ", Integral: " << integral
                      << std::endl;
            // ----------------------------------------------

            // Add to stack
            hs->Add(cloneHist, "hist");
            histClones.push_back(cloneHist);

            // Add to legend
            leg->AddEntry(cloneHist, label.c_str(), "f");
        } // end while keys in "mode"
    } // end loop over subDirs (CC & NC)

    // fill in our result
    result.stack = hs;
    result.legend = leg;
    result.histList = histClones; 
    return result;
}

// Main function: 4 subplots in one canvas, with CC/NC color offset,
// and also printing integrals to stdout
void eventspectra(const char* inputFileName = "outputs/gundamFitter_config_DUNE_Asimov_DryRun.root")
{
    // 1) Open your ROOT file
    TFile* f = TFile::Open(inputFileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file: " << inputFileName << std::endl;
        return;
    }

    // 2) Go to FitterEngine/preFit/plots/histograms
    TDirectoryFile* histDir = dynamic_cast<TDirectoryFile*>(
        f->Get("FitterEngine/preFit/plots/histograms")
    );
    if (!histDir) {
        std::cerr << "Directory not found: FitterEngine/preFit/plots/histograms\n";
        f->Close();
        return;
    }

    // 3) Define the four merged categories:
    //    Each merges 2 subdirectories (CC, NC).
    std::map<std::string, std::vector<std::string>> mergedSamples = {
        { "FHC #nu_{#mu} like", { "FHC #nu_{#mu}_like CC ", "FHC #nu_{#mu}_like NC " } },
        { "FHC #nu_{e} like",  { "FHC #nu_{e}_like CC ",  "FHC #nu_{e}_like NC "  } },
        { "RHC #nu_{#mu} like", { "RHC #nu_{#mu}_like CC ", "RHC #nu_{#mu}_like NC " } },
        { "RHC #nu_{e} like",  { "RHC #nu_{e}_like CC ",  "RHC #nu_{e}_like NC "  } }
    };

    // 4) Build mode -> (modeName, baseColor) map
    std::map<int, ModeInfo> modeMap = buildModeMap();

    // Optionally hide stats boxes
    gStyle->SetOptStat(0);

    // We'll store the results (stacks & legends) for the 4 samples
    // in the order we want them on the canvas
    std::vector<std::string> sampleOrder = {
        "FHC #nu_{#mu} like",
        "FHC #nu_{e} like",
        "RHC #nu_{#mu} like",
        "RHC #nu_{e} like"
    };

    std::map<std::string, SamplePlot> plots;

    // 5) Build each sample's stack
    for (auto &sampName : sampleOrder) {
        auto iter = mergedSamples.find(sampName);
        if (iter == mergedSamples.end()) {
            std::cerr << "No merged subdirs for " << sampName << std::endl;
            continue;
        }
        const auto &subdirs = iter->second;

        SamplePlot sp = buildStackForSample(sampName, subdirs, histDir, modeMap);
        if (!sp.stack) {
            std::cerr << "Skipping sample " << sampName << " (no stack)\n";
            continue;
        }
        plots[sampName] = sp;
    }

    // 6) Create one canvas with 2x2 subpads
    TCanvas* cAll = new TCanvas("cAll", "Four Merged Samples in One Canvas", 1400, 1000);
    cAll->Divide(2,2);

    // 7) Loop over each sample in sampleOrder, draw in subpad i
    int padIndex = 1;
    for (auto &sampName : sampleOrder) {
        // If we didn't build it, skip
        if (plots.find(sampName) == plots.end()) continue;

        // Switch to next subpad
        cAll->cd(padIndex);
        padIndex++;

        // Retrieve the THStack and TLegend
        THStack* hs = plots[sampName].stack;
        TLegend* leg = plots[sampName].legend;

        if (hs) {
            hs->Draw("hist");
            // limit visible x axis to 0–10 GeV
            // hs->GetHistogram()->GetXaxis()->SetRangeUser(0.0, 10.0);    

            hs->GetXaxis()->SetTitle("Reconstructed E (GeV)");
            hs->GetYaxis()->SetTitle("Events / bin");
        }
        if (leg) {
            leg->Draw();
        }
    }

    // 8) Save as one PDF (one page)
    cAll->SaveAs("outputs/plots/merged.pdf");

    // 9) Close the file
    f->Close();
}