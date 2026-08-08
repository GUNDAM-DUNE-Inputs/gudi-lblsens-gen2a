#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1D.h>
#include <TClass.h>
#include <iostream>

void getIntegrals(const char* inputFileName = "default.root")
{   
    // A lookup table from interaction mode number to descriptive text
    // (You can rename or adjust titles to match exactly what you want printed)
    static std::map<int, std::string> modeTitle = {
        {-1, "UNDEFINED"},
        {0,  "QE"},
        {1,  "Resonant"},
        {2,  "DIS"},
        {3,  "Coherent"},
        {4,  "Coherent Elastic"},
        {5,  "Electron Scattering"},
        {6,  "Inverse Muon Decay Annihliation"},
        {7,  "Inverse Beta Decay"},
        {8,  "Glasgow Resonance"},
        {9,  "Atmospheric Muon Nu Gamma"},
        {10, "MEC aka 2p2h"},
        {11, "Diffractive"},
        {12, "kEM"},
        {13, "kWeakMix"}
    };

    // Read Files
    TFile* f = TFile::Open(inputFileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Navigate to the directory: "FitterEngine/preFit/plots/histograms"
    TDirectory* histDir = dynamic_cast<TDirectory*>(
        f->Get("FitterEngine/preFit/plots/histograms")
    );
    if (!histDir) {
        std::cerr << "Could not find directory: FitterEngine/preFit/plots/histograms" << std::endl;
        f->Close();
        return;
    }

    // Loop over the eight subfolders (e.g., "FHC #nu_{#mu}_like CC", "FHC #nu_{#mu}_like NC", etc.)
    TIter nextSubfolder(histDir->GetListOfKeys());
    TKey* keySubfolder;
    while ((keySubfolder = (TKey*)nextSubfolder())) {

        // Only process directories
        if (strcmp(keySubfolder->GetClassName(), "TDirectoryFile") != 0) {
            continue;
        }

        // subDir is something like "FHC #nu_{#mu}_like CC"
        TDirectory* subDir = dynamic_cast<TDirectory*>(histDir->Get(keySubfolder->GetName()));
        if (!subDir) continue;

        // Inside subDir, we expect "Raw/mode"
        TDirectory* rawDir = dynamic_cast<TDirectory*>(subDir->Get("Raw"));
        if (!rawDir) {
            // Not all subfolders might have this structure
            continue;
        }

        TDirectory* modeDir = dynamic_cast<TDirectory*>(rawDir->Get("mode"));
        if (!modeDir) {
            continue;
        }

        // Identify if this subdirectory name has "CC" or "NC" in it (for printing)
        std::string subDirName = subDir->GetName();
        std::string ccOrNc;
        if (subDirName.find("CC") != std::string::npos) {
            ccOrNc = "CC"; 
        } 
        else if (subDirName.find("NC") != std::string::npos) {
            ccOrNc = "NC";
        } 
        // If neither CC nor NC is found, ccOrNc remains empty (which is OK)

        // Now loop over the numeric subfolders ("0", "1", "2", etc.) in "mode"
        TIter nextMode(modeDir->GetListOfKeys());
        TKey* keyMode;
        while ((keyMode = (TKey*)nextMode())) {

            if (strcmp(keyMode->GetClassName(), "TDirectoryFile") != 0) {
                continue;
            }

            TDirectory* modeSubDir = dynamic_cast<TDirectory*>(modeDir->Get(keyMode->GetName()));
            if (!modeSubDir) continue;

            // Convert the subdirectory name to an integer to match the mode
            int modeVal = std::atoi(modeSubDir->GetName());

            // Retrieve the TH1D named "MC_TH1D"
            TH1D* h = dynamic_cast<TH1D*>(modeSubDir->Get("MC_TH1D"));
            if (!h) {
                // Some modes might not have this histogram
                continue;
            }

            // Get the integral (sum of bin contents)
            double integral = h->Integral();

            // Find the descriptive title in our map (or "Unknown" if not found)
            std::string title = "UnknownMode";
            if (modeTitle.find(modeVal) != modeTitle.end()) {
                title = modeTitle[modeVal];
            }

            // If ccOrNc is non-empty, prepend it to the mode title (e.g. "CC MEC aka 2p2h")
            std::string outModeTitle = title;
            if (!ccOrNc.empty()) {
                outModeTitle = ccOrNc + " " + title;
            }

            // Print the result in the desired format
            // Example:
            //   Directory: RHC #nu_{#mu}_like CC , Mode: 10, CC MEC aka 2p2h, Integral: 1657.63
            std::cout << "Directory: " << subDirName
                      << " , Mode: " << modeVal
                      << ", " << outModeTitle
                      << ", Integral: " << integral
                      << std::endl;
        }
    }

    f->Close();   
}
