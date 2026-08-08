#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

#include <yaml-cpp/yaml.h>

#include "TFile.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TVectorD.h"
#include "TObjArray.h"
#include "TObjString.h"

double readDouble(const YAML::Node& node, const std::string& label, int i)
{
    if (!node) {
        std::cerr << "Missing node: " << label
                  << " for parameter index " << i << std::endl;
        throw std::runtime_error("Missing YAML node");
    }

    try {
        return node.as<double>();
    }
    catch (const std::exception&) {
        std::cerr << "Bad conversion at " << label
                  << " for parameter index " << i
                  << ". Node content: " << node << std::endl;
        throw;
    }
}

std::string readString(const YAML::Node& node, const std::string& label, int i)
{
    if (!node) {
        std::cerr << "Missing node: " << label
                  << " for parameter index " << i << std::endl;
        throw std::runtime_error("Missing YAML node");
    }

    try {
        return node.as<std::string>();
    }
    catch (const std::exception&) {
        std::cerr << "Bad conversion at " << label
                  << " for parameter index " << i
                  << ". Node content: " << node << std::endl;
        throw;
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.yaml" << std::endl;
        return 1;
    }

    const std::string yamlFile = argv[1];

    YAML::Node root = YAML::LoadFile(yamlFile);
    YAML::Node systematicsHdr = root["Systematics"];

    if (!systematicsHdr || !systematicsHdr.IsSequence()) {
        std::cerr << "Error: 'Systematics' node is missing or is not a sequence." << std::endl;
        return 1;
    }

    const int N = static_cast<int>(systematicsHdr.size());
    std::cout << "Found " << N << " systematics." << std::endl;

    // Output object names
    const std::string parnameStr = "xsec_param_names";
    const std::string nupdgStr   = "xsec_norm_nupdg";

    // ROOT containers
    TObjArray* parnameArray = new TObjArray();
    parnameArray->SetOwner(kTRUE);
    parnameArray->SetName(parnameStr.c_str());

    TObjArray* nupdgArray = new TObjArray();
    nupdgArray->SetOwner(kTRUE);
    nupdgArray->SetName(nupdgStr.c_str());

    TVectorT<Double_t> priorvalVec(N);
    TVectorT<Double_t> nomvalVec(N);
    TVectorT<Double_t> paramvalLbVec(N);
    TVectorT<Double_t> paramvalUbVec(N);
    TVectorT<Double_t> kinematicLbVec(N);
    TVectorT<Double_t> kinematicUbVec(N);
    TVectorT<Double_t> stepscaleVec(N);

    TVectorD* nupdgVec = new TVectorD(N);

    TMatrixT<Double_t> covmat(N, N);
    covmat.Zero();

    // Metadata
    std::vector<std::string> paramNames(N);
    std::vector<double> errorVec(N);
    std::unordered_map<std::string, int> nameToIndex;

    // ---------------- First pass: names and errors ----------------
    for (int i = 0; i < N; ++i) {
        YAML::Node systematic = systematicsHdr[i]["Systematic"];
        if (!systematic) {
            std::cerr << "Missing 'Systematic' block at index " << i << std::endl;
            return 1;
        }

        std::string paramName =
            readString(systematic["Names"]["ParameterName"], "Names/ParameterName", i);
        double errval =
            readDouble(systematic["Error"], "Error", i);

        paramNames[i] = paramName;
        errorVec[i] = errval;

        if (nameToIndex.find(paramName) != nameToIndex.end()) {
            std::cerr << "Duplicate parameter name found: " << paramName << std::endl;
            return 1;
        }

        nameToIndex[paramName] = i;
    }

    // ---------------- Second pass: vectors and covariance ----------------
    for (int i = 0; i < N; ++i) {
        YAML::Node systematic = systematicsHdr[i]["Systematic"];
        const std::string& paramName = paramNames[i];

        std::cout << "Processing parameter " << i
                  << " : " << paramName << std::endl;

        // Store parameter name
        TObjString* s = new TObjString();
        s->SetString(paramName.c_str());
        parnameArray->Add(s);

        // Optional neutrino flavour
        if (systematic["NeutrinoFlavourUnosc"] &&
            systematic["NeutrinoFlavourUnosc"].IsSequence() &&
            systematic["NeutrinoFlavourUnosc"].size() > 0) {
            (*nupdgVec)[i] = systematic["NeutrinoFlavourUnosc"][0].as<double>();
        } else {
            (*nupdgVec)[i] = -999.0;
        }

        // Required values
        priorvalVec[i]   = readDouble(systematic["ParameterValues"]["PreFitValue"], "ParameterValues/PreFitValue", i);
        nomvalVec[i]     = readDouble(systematic["ParameterValues"]["Generated"],   "ParameterValues/Generated", i);
        paramvalLbVec[i] = readDouble(systematic["ParameterBounds"][0], "ParameterBounds[0]", i);
        paramvalUbVec[i] = readDouble(systematic["ParameterBounds"][1], "ParameterBounds[1]", i);

        // Optional kinematic cuts
        if (systematic["KinematicCuts"] &&
            systematic["KinematicCuts"].IsSequence() &&
            systematic["KinematicCuts"].size() > 0 &&
            systematic["KinematicCuts"][0]["TrueNeutrinoEnergy"] &&
            systematic["KinematicCuts"][0]["TrueNeutrinoEnergy"].IsSequence() &&
            systematic["KinematicCuts"][0]["TrueNeutrinoEnergy"].size() > 1) {
            kinematicLbVec[i] = systematic["KinematicCuts"][0]["TrueNeutrinoEnergy"][0].as<double>();
            kinematicUbVec[i] = systematic["KinematicCuts"][0]["TrueNeutrinoEnergy"][1].as<double>();
        } else {
            kinematicLbVec[i] = -999.0;
            kinematicUbVec[i] = -999.0;
        }

        // Optional step scale
        if (systematic["StepScale"] && systematic["StepScale"]["MCMC"]) {
            stepscaleVec[i] = systematic["StepScale"]["MCMC"].as<double>();
        } else {
            stepscaleVec[i] = 1.0;
        }

        // Diagonal = sigma^2
        covmat(i, i) = errorVec[i] * errorVec[i];

        // Off-diagonals from named correlation entries
        YAML::Node corrList = systematic["Correlations"];
        if (!corrList || !corrList.IsSequence()) {
            std::cerr << "Missing or invalid Correlations list for parameter "
                      << i << " (" << paramName << ")" << std::endl;
            return 1;
        }

        for (auto entry : corrList) {
            if (!entry.IsMap()) {
                std::cerr << "Correlation entry is not a map for parameter "
                          << i << " (" << paramName << ")" << std::endl;
                return 1;
            }

            for (auto it = entry.begin(); it != entry.end(); ++it) {
                std::string corrName = it->first.as<std::string>();
                double corrval = it->second.as<double>();

                auto found = nameToIndex.find(corrName);
                if (found == nameToIndex.end()) {
                    std::cerr << "Unknown correlation target '" << corrName
                              << "' in row " << paramName << std::endl;
                    return 1;
                }

                int j = found->second;

                // Keep diagonal from Error^2 even if self-entry appears
                if (j == i) {
                    continue;
                }

                covmat(i, j) = corrval * errorVec[i] * errorVec[j];
            }
        }
    }

    // ---------------- Symmetrize only where one side is missing ----------------
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            const double a = covmat(i, j);
            const double b = covmat(j, i);

            if (a == 0.0 && b != 0.0) {
                covmat(i, j) = b;
            }
            else if (b == 0.0 && a != 0.0) {
                covmat(j, i) = a;
            }
        }
    }

    // ---------------- Diagnostics ----------------
    std::cout << "\nChecking diagonal elements...\n";
    for (int i = 0; i < N; ++i) {
        if (covmat(i, i) <= 0.0) {
            std::cerr << "WARNING: non-positive diagonal at i=" << i
                      << " name=" << paramNames[i]
                      << " value=" << covmat(i, i) << std::endl;
        }
    }

    std::cout << "Checking symmetry...\n";
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double diff = std::abs(covmat(i, j) - covmat(j, i));
            if (diff > 1e-10) {
                std::cerr << "WARNING: covariance matrix is not symmetric at ("
                          << i << "," << j << ") : "
                          << covmat(i, j) << " vs " << covmat(j, i) << std::endl;
            }
        }
    }

    nupdgArray->Add(nupdgVec);

    // ---------------- Write output ----------------
    TFile outFile("flux_xsec_covariance_DUNE_inclNDcontr.root", "RECREATE");

    outFile.WriteObject(parnameArray, parnameStr.c_str());
    outFile.WriteObject(nupdgArray, nupdgStr.c_str());

    priorvalVec.Write("xsec_param_prior");
    nomvalVec.Write("xsec_param_nom");
    paramvalLbVec.Write("xsec_param_lb");
    paramvalUbVec.Write("xsec_param_ub");
    kinematicLbVec.Write("xsec_norm_kinematic_lb");
    kinematicUbVec.Write("xsec_norm_kinematic_ub");
    stepscaleVec.Write("xsec_stepscale");
    covmat.Write("xsec_cov");

    outFile.Close();

    std::cout << "\nWrote output file: flux_xsec_covariance_DUNE_inclNDcontr" << std::endl;
    return 0;
}