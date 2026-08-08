#include <TFile.h>
#include <TMatrixT.h>
#include <TDecompChol.h>
#include <TKey.h>
#include <TObjString.h>
#include <iostream>
#include <cmath>
#include <map>

void decompose_with_regularization() {
    // Open the new ROOT file (a copy of the original file)
    TFile *file = TFile::Open("../inputs/parameters/xsec/regularize_xsec_covariance_DUNE_systs_2023.root", "UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve the TMatrixT<double> object
    TMatrixT<double> *xsec_cov = nullptr;
    file->GetObject("xsec_cov", xsec_cov);
    if (!xsec_cov) {
        std::cerr << "Error: could not retrieve TMatrixT<double> object!" << std::endl;
        file->Close();
        return;
    }

    // Print out elements greater than 0.99
    std::cout << "Elements greater than 0.99 in the original matrix:" << std::endl;
    for (int i = 0; i < xsec_cov->GetNrows(); ++i) {
        for (int j = 0; j < xsec_cov->GetNcols(); ++j) {
            if ((*xsec_cov)(i, j) > 0.99) {
                std::cout << "Element at (" << i << ", " << j << "): " << (*xsec_cov)(i, j) << std::endl;
            }
        }
    }
    
    // Initialize regularization parameter
    double lambda = 1e-6;  // Start with a slightly larger regularization value
    double increment = 1e-6;
    bool decomposable = false;

    // Regularize until decomposable
    TMatrixT<double> regularized_matrix = *xsec_cov;
    while (!decomposable) {
        for (int i = 0; i < regularized_matrix.GetNrows(); ++i) {
            regularized_matrix(i, i) += lambda;
        }

        // Attempt Cholesky decomposition
        TDecompChol decomp(regularized_matrix);
        decomposable = decomp.Decompose();

        if (!decomposable) {
            lambda += increment;
            if (lambda > 1.0) {
                std::cerr << "Error: Matrix could not be decomposed even after significant regularization." << std::endl;
                file->Close();
                delete file;
                return;
            }
        }
    }

    // Print the final regularization value
    std::cout << "The final regularization value needed to make the matrix decomposable: " << lambda << std::endl;

    // Save the regularized matrix and a new TString into the same file
    file->cd();
    regularized_matrix.Write("regularized_xsec_cov");

    // Clean up
    file->Close();
    delete file;
}
