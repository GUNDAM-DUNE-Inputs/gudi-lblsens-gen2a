//  -------------------------------------------------------
//  Run, e.g.:
//      root -l -b -q makeOscCov.cpp

#include <TFile.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMatrixT.h>
#include <iostream>

void makeOscCov(const char* outfile = "../../inputs/parameters/oscProb/oscCovTDR_NO.root")
{
    // ------------------------------------------------------------------ //
    // 1. Parameter names
    // ------------------------------------------------------------------ //
    const char* names[] = {
        "PMNS_SIN_SQUARED_12",
        "PMNS_SIN_SQUARED_13",
        "PMNS_SIN_SQUARED_23",
        "PMNS_DELTA_MASS_SQUARED_21",
        "PMNS_DELTA_MASS_SQUARED_32",
        "PMNS_DELTA_CP",
        "PMNS_SIGN_MASS_SQUARED_32"
    };
    constexpr int N = sizeof(names) / sizeof(names[0]);

    auto* nameArray = new TObjArray(N);
    nameArray->SetName("osc_param_names");
    for (int i = 0; i < N; ++i) {
        nameArray->Add(new TObjString(names[i]));
    }

    // ------------------------------------------------------------------ //
    // 2. Define the parameter prior values and uncertainties
    // ------------------------------------------------------------------ //

    // PDG: http://www.nu-fit.org/sites/default/files/v40.tbl-parameters.pdf
    double parVals[N];
    double parSigs[N];

    //PMNS_SIN_SQUARED_12
    parVals[0] = 0.310;
    parSigs[0] = 0.013;

    //NO: PMNS_SIN_SQUARED_13  PDG 2018: 0.02240 +/- 0.00065/0.00066;
    //IO: PMNS_SIN_SQUARED_13  PDG 2018: 0.02263 +/- 0.00065/0.00066;
    parVals[1] = 0.0224;  // normal
    parSigs[1] = 0.0007;  // normal

    // PMNS_SIN_SQUARED_23  PDG 2018: 0.582 +0.015/-0.019  (normal ordering)
    // PMNS_SIN_SQUARED_23  PDG 2018: 0.582 +0.015/-0.018  (inverted ordering)
     parVals[2] = 0.582;  // normal
     parSigs[2] = 0.021;  // normal


    // PMNS_DELTA_MASS_SQUARED_21  PDG 2018: 7.39E-5 +/- 0.21E-5/0.20E-5

    parVals[3] = 7.39E-5;
    parSigs[3] = 1.8E-6;

    // PMNS_DELTA_MASS_SQUARED_32 should be free in any fit.
    //
    // PMNS_DELTA_MASS_SQUARED_32  PDG 2018:  2.525E-3 +/- 0.033E-3/0.031E-3 (normal)
    // PMNS_DELTA_MASS_SQUARED_32  PDG 2018: -2.512E-3 +/- 0.034E-3/0.031E-3 (inverted)

    parVals[4] = 2.525E-3; // Normal
    parSigs[4] = 3.4E-5; // Normal

    // PMNS_DELTA_CP PDG 2018: -2.496 +/- 0.698/0.489 (normal)
    // PMNS_DELTA_CP PDG 2018: -1.396 +/- 0.489/0.436 (inverted)
    parVals[5] = -2.498;
    parSigs[5] = 1.17;

    // PMNS_SIGN_MASS_SQUARED_32
    parVals[6] = 0.5;
    parSigs[6] = 0.1;  // Mostly unconstrained

    // ------------------------------------------------------------------ //
    // 3. Create the covariance matrix and prior vector
    // ------------------------------------------------------------------ //
    auto* cov = new TMatrixT<double>(N, N);
    cov->Zero();
    for (int i = 0; i < N; ++i){
        (*cov)(i, i) = parSigs[i]*parSigs[i];
    }
    auto* prior = new TVectorT<double>(N,parVals);

    // ------------------------------------------------------------------ //
    // 4. Write to file
    // ------------------------------------------------------------------ //
    TFile f(outfile, "RECREATE");

    //  Write the names array as ONE key only:
    //     kSingleKey = do NOT write the six sub-objects as individual keys
    nameArray->Write("osc_param_names", TObject::kSingleKey);

    // Write the priors
    f.WriteObject(prior, "osc_param_priors");

    // Write the matrix
    f.WriteObject(cov, "osc_param_cov");

    f.Close();

    std::cout << "Wrote " << outfile << " with:" << std::endl
              << "  • osc_param_names (TObjArray, single key)"  << std::endl
              << "  • osc_param_priors (TVectorT<double>)"  << std::endl
              << "  • osc_param_cov   (TMatrixT<double>)" << std::endl;
}
