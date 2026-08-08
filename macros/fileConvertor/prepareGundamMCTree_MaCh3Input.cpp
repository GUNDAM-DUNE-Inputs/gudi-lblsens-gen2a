#include <iostream>
#include <string>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <vector>
#include <map>



#include "TSystem.h"

using namespace TMath;

const double DELMSQ_31 = 2.515e-3; // In eV^2
const double LOSC = 1300.; // In km

const double THETA_23 = 0.859;
const double THETA_13 = 0.150;

// OA_INPUTS is /storage/shared/DUNE/OA-inputs
std::vector<TString> FHCnonswap_MaCh3 = {
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nuebar_x_nuebar_nueselec.root",      //file_idx=1
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nuebar_x_nuebar_numuselec.root",     //file_idx=2
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nue_x_nue_nueselec.root",        //file_idx=3
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nue_x_nue_numuselec.root",       //file_idx=4
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numubar_x_numubar_nueselec.root",        //file_idx=5
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numubar_x_numubar_numuselec.root",       //file_idx=6
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numu_x_numu_nueselec.root",      //file_idx=7
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numu_x_numu_numuselec.root"      //file_idx=8
};
std::vector<TString> FHCnueswap_MaCh3 = {
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nuebar_x_nutaubar_nueselec.root",        //file_idx=9
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nuebar_x_nutaubar_numuselec.root",       //file_idx=10
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nue_x_nutau_nueselec.root",      //file_idx=11
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nue_x_nutau_numuselec.root",     //file_idx=12
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numubar_x_nuebar_nueselec.root",     //file_idx=13
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numubar_x_nuebar_numuselec.root",        //file_idx=14
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numu_x_nue_nueselec.root",       //file_idx=15
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numu_x_nue_numuselec.root"       //file_idx=16
};
std::vector<TString> FHCtauswap_MaCh3 = {
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nuebar_x_numubar_nueselec.root",     //file_idx=17   
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nuebar_x_numubar_numuselec.root",        //file_idx=18
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nue_x_numu_nueselec.root",       //file_idx=19
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_nue_x_numu_numuselec.root",      //file_idx=20
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numubar_x_nutaubar_nueselec.root",       //file_idx=21
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numubar_x_nutaubar_numuselec.root",      //file_idx=22
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numu_x_nutau_nueselec.root",     //file_idx=23
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_FHC_ger_numu_x_nutau_numuselec.root"     //file_idx=24
};
std::vector<TString> RHCnonswap_MaCh3 = {
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nuebar_x_nuebar_nueselec.root",       //file_idx=25
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nuebar_x_nuebar_numuselec.root",       //file_idx=26
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nue_x_nue_nueselec.root",       //file_idx=27
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nue_x_nue_numuselec.root",       //file_idx=28
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numubar_x_numubar_nueselec.root",       //file_idx=29
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numubar_x_numubar_numuselec.root",       //file_idx=30
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numu_x_numu_nueselec.root",       //file_idx=31
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numu_x_numu_numuselec.root"       //file_idx=32
};
std::vector<TString> RHCnueswap_MaCh3 = {
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nuebar_x_nutaubar_nueselec.root",       //file_idx=33
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nuebar_x_nutaubar_numuselec.root",       //file_idx=34
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nue_x_nutau_nueselec.root",       //file_idx=35
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nue_x_nutau_numuselec.root",       //file_idx=36
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numubar_x_nuebar_nueselec.root",       //file_idx=37
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numubar_x_nuebar_numuselec.root",       //file_idx=38
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numu_x_nue_nueselec.root",       //file_idx=39
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numu_x_nue_numuselec.root"       //file_idx=40
};
std::vector<TString> RHCtauswap_MaCh3 = {
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nuebar_x_numubar_nueselec.root",       //file_idx=41
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nuebar_x_numubar_numuselec.root",       //file_idx=42
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nue_x_numu_nueselec.root",       //file_idx=43
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_nue_x_numu_numuselec.root",       //file_idx=44
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numubar_x_nutaubar_nueselec.root",       //file_idx=45
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numubar_x_nutaubar_numuselec.root",       //file_idx=46
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numu_x_nutau_nueselec.root",       //file_idx=47
    "${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_RHC_ger_numu_x_nutau_numuselec.root"       //file_idx=48
};


double POTperYear = 1.1e21;
double FHCnonswapPOT = 1.62824e24;
double FHCnueswapPOT = 1.62063e24;
double FHCtauswapPOT = 5.18551e24;

double RHCnonswapPOT = 3.27608e24;
double RHCnueswapPOT = 3.24713e24;
double RHCtauswapPOT = 8.58955e24;

double scaleCorrection = 40. / 1.13; //1.13 is kt
// For the 7-year nominal exposure is 40/1.13 * 1.1e21 * 3.5 which is the same
// for FHC and RHC since we're splitting half/half
double POT_desired = scaleCorrection * POTperYear * (3.5);


void fillArrayWithWeights(int size, std::vector<double>& array, int centralIndex, double centralValue, double shiftStep) {
    for (int i = 0; i < size; i++) {
        double value = centralValue + (i - centralIndex) * shiftStep;
        array.push_back(value);
    }
}

void expandAndAdd(std::vector<TString>& output, const std::vector<TString>& input) {
    for (TString name : input) {
        gSystem->ExpandPathName(name);
        output.emplace_back(name);
    }
}

TString createOutputFileName(const TString& inputFile) {
    TString outputDir = gSystem->Getenv("OA_OUTPUTS");  // Expand env variable
    gSystem->Exec(Form("mkdir -p %s", outputDir.Data())); // Ensure the output directory exists

    TString baseFileName = gSystem->BaseName(inputFile);  // Extract file name
    TString outputFileName = outputDir + "/" + baseFileName;
    outputFileName.ReplaceAll(".root", "_gudi.root");

    return outputFileName;
}

void CopyMCTree(const std::string& inputFile) {

    std::vector<TString> all_files;
    expandAndAdd(all_files, FHCnonswap_MaCh3);   // index 1–8
    expandAndAdd(all_files, FHCnueswap_MaCh3);   // index 9–16
    expandAndAdd(all_files, FHCtauswap_MaCh3);   // index 17–24
    expandAndAdd(all_files, RHCnonswap_MaCh3);   // index 25–32
    expandAndAdd(all_files, RHCnueswap_MaCh3);   // index 33–40
    expandAndAdd(all_files, RHCtauswap_MaCh3);   // index 41–48

    bool useRHC = false;
    if (inputFile.find("RHC") != std::string::npos) useRHC = true;
    else if (inputFile.find("FHC") != std::string::npos) useRHC = false;
    else {
        std::cerr << "Input file not part of the FHC and RHC input files" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    int sample_idx = -1;
    if (!useRHC) {
        if (inputFile.find("numuselec") != std::string::npos) sample_idx = 1;
        else if (inputFile.find("nueselec") != std::string::npos) sample_idx = 2;
        else {
            std::cerr << "Input file not part of the nue and numu selection input files" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    } else {

        if (inputFile.find("numuselec") != std::string::npos) sample_idx = 3;
        else if (inputFile.find("nueselec") != std::string::npos) sample_idx = 4;
        else {
            std::cerr << "Input file not part of the nue and numu selection input files" << std::endl;
            std::exit(EXIT_FAILURE);
        }

    }

    std::vector<TString> nonswap;
    std::vector<TString> nueswap;
    std::vector<TString> tauswap;

    // Print information for debugging purposes
    std::cout << "Processing " << (useRHC ? "RHC" : "FHC") << " file " << inputFile << std::endl;

    if (!useRHC) {
        expandAndAdd(nonswap,FHCnonswap_MaCh3);
        expandAndAdd(nueswap,FHCnueswap_MaCh3);
        expandAndAdd(tauswap,FHCtauswap_MaCh3);
    } else {
        expandAndAdd(nonswap,RHCnonswap_MaCh3);
        expandAndAdd(nueswap,RHCnueswap_MaCh3);
        expandAndAdd(tauswap,RHCtauswap_MaCh3);
    }

    TString filename = "";

    // Fill the chain to be copied.
    std::unique_ptr<TChain> ch1 = std::make_unique<TChain>("caf");
    ch1->Add(inputFile.c_str());

    // Create the output file with a unique suffix
    //std::string baseName = gSystem->BaseName(inputFile.c_str());
    //std::string outputFileName = baseName.substr(0,baseName.find(".root"));
    //outputFileName = "${OA_OUTPUTS}" + outputFileName + "_gudi.root";

    TString outputFileName = createOutputFileName(inputFile);
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    int isCC = 0;
    double POTScaledOscweight = 0.0;
    double POT_generated = 0.0;
    int isNC = 0;
    int isRHC = 0;

    int SwapType = -1; // 1: nonswap files, 2: nueswap files, 3: tauswap files, -1: unknown file (an error)

    ch1->SetBranchStatus("*", true);
    ch1->SetBranchAddress("isCC", &isCC);

    std::vector<std::string> branchPrefixes = {
        "FSILikeEAvailSmearing",
        "SPPLowQ2Suppression",
        "nuenuebar_xsec_ratio",
        "nuenumu_xsec_ratio",
        "C12ToAr40_2p2hScaling_nubar",
        "C12ToAr40_2p2hScaling_nu",
        "EbFSLepMomShift",
        "BeRPA_E",
        "BeRPA_D",
        "BeRPA_B",
        "BeRPA_A",
        "NR_nubar_p_NC_3Pi",
        "NR_nubar_p_NC_2Pi",
        "NR_nubar_p_NC_1Pi",
        "NR_nubar_n_NC_3Pi",
        "NR_nubar_n_NC_2Pi",
        "NR_nubar_n_NC_1Pi",
        "NR_nubar_p_CC_3Pi",
        "NR_nubar_p_CC_2Pi",
        "NR_nubar_p_CC_1Pi",
        "NR_nubar_n_CC_3Pi",
        "NR_nubar_n_CC_2Pi",
        "NR_nubar_n_CC_1Pi",
        "NR_nu_p_NC_3Pi",
        "NR_nu_p_NC_2Pi",
        "NR_nu_p_NC_1Pi",
        "NR_nu_n_NC_3Pi",
        "NR_nu_n_NC_2Pi",
        "NR_nu_n_NC_1Pi",
        "NR_nu_np_CC_1Pi",
        "NR_nu_p_CC_3Pi",
        "NR_nu_p_CC_2Pi",
        "NR_nu_n_CC_3Pi",
        "NR_nu_n_CC_2Pi",
        "E2p2h_B_nubar",
        "E2p2h_A_nubar",
        "E2p2h_B_nu",
        "E2p2h_A_nu",
        "MKSPP_ReWeight",
        "Mnv2p2hGaussEnhancement",
        "CCQEPauliSupViaKF",
        "FrPiProd_N",
        "FrAbs_N",
        "FrInel_N",
        "FrElas_N",
        "FrCEx_N",
        "MFP_N",
        "FrPiProd_pi",
        "FrAbs_pi",
        "FrInel_pi",
        "FrElas_pi",
        "FrCEx_pi",
        "MFP_pi",
        "FormZone",
        "CV2uBY",
        "CV1uBY",
        "BhtBY",
        "AhtBY",
        "Theta_Delta2Npi",
        "RDecBR1eta",
        "RDecBR1gamma",
        "MvNCRES",
        "MaNCRES",
        "MvCCRES",
        "MaCCRES",
        "EtaNCEL",
        "MaNCEL",
        "VecFFCCQEshape",
        "MaCCQE"
    };

    std::map<std::string, int> nshiftsMap;
    std::map<std::string, double> cvwgtMap;
    std::map<std::string, Double_t*> wgtMap;  // Pointer for array of weights
    std::map<std::string, std::vector<Double_t>> Xarray;
    std::map<std::string, double> priorMap;
    std::map<std::string, double> sigmaMap;
    double prior = 0;

    for (const auto& prefix : branchPrefixes) {


        std::string prefix_centr = prefix + "_cvwgt";
        std::string prefix_wgt = "wgt_"+prefix;

        if (prefix == "EbFSLepMomShift"){

            prefix_centr = prefix + "_cvvar";
            prefix_wgt = "var_"+prefix;

        }

        ch1->SetBranchStatus((prefix + "_nshifts").c_str(), 1);
        ch1->SetBranchStatus((prefix_centr).c_str(), 1);
        ch1->SetBranchStatus((prefix_wgt).c_str(), 1);

        nshiftsMap[prefix] = 7;
        cvwgtMap[prefix] = 0.0;
        wgtMap[prefix] = nullptr;  // Initialize to nullptr
        priorMap[prefix] = 0;
        sigmaMap[prefix] = 1;
        if (prefix == "MaCCRES"){ priorMap[prefix] = -0.36; sigmaMap[prefix] = 0.0900056; }
        if (prefix == "NR_nu_np_CC_1Pi"){ priorMap[prefix] = -1.14; sigmaMap[prefix] = 0.200002; }        
        if (prefix == "nuenuebar_xsec_ratio") priorMap[prefix] = 1;
        if (prefix == "nuenumu_xsec_ratio") priorMap[prefix] = 1;

        ch1->SetBranchAddress((prefix + "_nshifts").c_str(), &nshiftsMap[prefix]);
        ch1->SetBranchAddress((prefix_centr).c_str(), &cvwgtMap[prefix]);

        // Make sure to clean up old memory before allocating new array
        if (wgtMap[prefix] != nullptr) {
            delete[] wgtMap[prefix];  // Free previously allocated memory
            wgtMap[prefix] = nullptr;  // Reset to nullptr to prevent double deletion
        }

        // Only allocate if there are valid shifts
        if (nshiftsMap[prefix] > 0) {
            wgtMap[prefix] = new Double_t[nshiftsMap[prefix]];  // Allocate new memory for the array
            ch1->SetBranchAddress((prefix_wgt).c_str(), wgtMap[prefix]);  // Set the new branch address
        } else {
            std::cerr << "Warning: nshiftsMap[" << prefix << "] is zero or negative." << std::endl;
        }
    }

    TTree* event_tree = ch1->CloneTree(0);  // Cloning structure only (no events)

    event_tree->SetBranchStatus("*", true);

    event_tree->Branch("POTScaledWeight", &POTScaledOscweight);
    event_tree->Branch("POT_generated",&POT_generated);
    event_tree->Branch("SwapType", &SwapType);
    event_tree->Branch("isNC", &isNC);
    event_tree->Branch("isRHC", &isRHC);
    event_tree->Branch("sample_idx", &sample_idx);

    // add file_idx for each file
    int file_idx = -1;
    event_tree->Branch("file_idx", &file_idx, "file_idx/I");
        


    std::vector<TClonesArray *> response_splines;
    std::vector<TBranch *> branch_vector;

    for (const auto& prefix : branchPrefixes) {

        response_splines.push_back(new TClonesArray("TGraph", 1));
        branch_vector.push_back(event_tree->Branch((prefix+"_TGraph").c_str(), "TClonesArray", response_splines.data()[(int)(response_splines.size() -1)], 256000, 0));

    }

    std::cout<<"response_splines size: " << response_splines.size() << std::endl;
    std::cout<<"branch_vector size:  " << branch_vector.size() << std::endl;

    // Get the number of entries in the chain
    int nentries = ch1->GetEntries();
    std::cout << "Copy " << nentries << " entries" << std::endl;
    for (int i = 0; i < nentries; ++i) {
        ch1->GetEntry(i);
        // Reset variables for each entry
        SwapType = -1;
        isNC = 0;
        isRHC = 0;
        POTScaledOscweight = 0.0;
        POT_generated = 0.0;

        filename = ch1->GetCurrentFile()->GetName();
        for (size_t i = 0; i < all_files.size(); ++i) {
            if (filename == all_files[i]) {
                file_idx = i + 1;  // make it 1-based
                break;
            }
        }

        if (file_idx <= 0) {
            std::cerr << "Error: couldn't find input file in known list" << std::endl;
            std::exit(EXIT_FAILURE);
        }



        if (std::find(nonswap.begin(), nonswap.end(), filename) != nonswap.end()) {
            SwapType = 1;
        }else if (std::find(nueswap.begin(), nueswap.end(), filename) != nueswap.end()) {
            SwapType = 2;
        }else if (std::find(tauswap.begin(), tauswap.end(), filename) != tauswap.end()) {
            SwapType = 3;
        }else {
            std::cout << "Input file not in the list of known files" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (isCC==0 && std::find(nonswap.begin(), nonswap.end(), filename) != nonswap.end()) isNC = 1;

        if(useRHC) isRHC = 1;

        if (SwapType == 1 && isRHC == 0) {
            POT_generated = FHCnonswapPOT;
        }
        else if (SwapType == 2 && isRHC == 0) {
            POT_generated = FHCnueswapPOT;
        }
        else if (SwapType == 3 && isRHC == 0) {
            POT_generated = FHCtauswapPOT;
        }
        else if (SwapType == 1 && isRHC == 1) {
            POT_generated = RHCnonswapPOT;
        }
        else if (SwapType == 2 && isRHC == 1) {
            POT_generated = RHCnueswapPOT;
        }
        else if (SwapType == 3 && isRHC == 1) {
            POT_generated = RHCtauswapPOT;
        }

        if (POT_generated>0){
            POTScaledOscweight = POT_desired * 1.0/ POT_generated;
        } else {
            POTScaledOscweight = 0.;
        }

        int parCount = -1;
        for (const auto& prefix : branchPrefixes) {
            parCount++;

            response_splines[parCount]->Clear();  // Clear the TClonesArray for the new event
            Xarray[prefix].clear();

            // Check for valid shifts before processing
            if (nshiftsMap[prefix] > 0) {
                fillArrayWithWeights(nshiftsMap[prefix], Xarray[prefix], nshiftsMap[prefix] / 2, priorMap[prefix], sigmaMap[prefix]);

                TClonesArray& arr_tmp = *response_splines[parCount];
                new(arr_tmp[0]) TGraph();  // Construct new TGraph in TClonesArray
                TGraph* graph = (TGraph*)(arr_tmp[0]);

                if (nshiftsMap[prefix] != 7) std::cout<<prefix<<": "<<nshiftsMap[prefix]<<std::endl;

                // Fill the graph with points
                for (size_t j = 0; j < (size_t)nshiftsMap[prefix]; j++) {
                    double x = Xarray[prefix][j];
                    double y;
                    
                    if (prefix == "BeRPA_A") y = wgtMap[prefix][j];
                    else y = wgtMap[prefix][j]*cvwgtMap[prefix];  // Now properly initialized

                    graph->SetPoint(j, x, y);
                }
                graph->SetName(prefix.c_str());
                graph->Sort();
            } else {
                std::cerr << "Skipping prefix: " << prefix << " due to zero or negative shifts." << std::endl;
            }
        }

        event_tree->Fill();
    }

    std::cout<<"Cleaning!"<<std::endl;

    outputFile->cd();
    outputFile->WriteObject(event_tree, "event_tree");
    outputFile->Close();

}

int main(int argc, char* argv[]) {
    // Check if the correct number of arguments is provided
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input-file>" << std::endl;
        std::cerr << "    Arguments: " << argc << std::endl;
        std::cerr << "    input-file -- The file to be converted." << std::endl;
        return 1;
    }

    // Parse the command-line arguments
    std::string inputFile{argv[1]};            // Get the file to convert.

    // Call the main function to process the ROOT file
    CopyMCTree(inputFile);

    return 0;
}
