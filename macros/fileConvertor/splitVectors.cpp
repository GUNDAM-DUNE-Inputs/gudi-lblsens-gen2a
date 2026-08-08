// Include necessary headers
#include <TFile.h>
#include <TKey.h>
#include <TVectorD.h>
#include <iostream>

void splitVectors() {
    // Open the input ROOT file
    // TFile* input_file = TFile::Open("../inputs/parameters/xsec/xsec_covariance_DUNE_systs_2023.root", "READ");
    // TFile* input_file = TFile::Open("../inputs/parameters/xsec/regularize_xsec_covariance_DUNE_systs_2023.root", "READ");
    TFile* input_file = TFile::Open("../inputs/parameters/xsec/xsec_covariance_DUNE_systs_2022a_FD_v3.root", "READ");
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "Error: Cannot open input.root" << std::endl;
        return;
    }

    // Create output ROOT files
    // TFile* output_file_first = TFile::Open("../inputs/parameters/xsec/xsec_covariance_only_DUNE_systs_2023.root", "RECREATE");
    // TFile* output_file_first = TFile::Open("../inputs/parameters/xsec/regularize_xsec_covariance_only_DUNE_systs_2023.root", "RECREATE");
    TFile* output_file_first = TFile::Open("../inputs/parameters/xsec/xsec_covariance_only_DUNE_systs_2022a_FD_v3.root", "RECREATE");
    // TFile* output_file_second = TFile::Open("../inputs/parameters/flux/flux_covariance_DUNE_systs_2023.root", "RECREATE");
    // TFile* output_file_second = TFile::Open("../inputs/parameters/flux/regularize_flux_covariance_DUNE_systs_2023.root", "RECREATE");
    TFile* output_file_second = TFile::Open("../inputs/parameters/flux/flux_covariance_DUNE_systs_2022a_FD_v3.root", "RECREATE");

    // Define the split point
    // int split_point = 99;   // for xsec_covariance_only_DUNE_systs_2023
    int split_point = 56;   // for flux_covariance_DUNE_systs_2022a_FD_v3
    

    // Get the list of keys in the input file
    TList* key_list = input_file->GetListOfKeys();
    TIter next(key_list);
    TKey* key;

    // Loop over each object in the input file
    while ((key = (TKey*)next())) {
        std::string obj_name = key->GetName();
        TObject* obj = input_file->Get(obj_name.c_str());

        // Process TVectorD objects
        if (auto vector = dynamic_cast<TVectorD*>(obj)) {
            int n = vector->GetNrows();
            int split_n = std::min(split_point, n);

            // Create two new TVectorD objects
            TVectorD first_part(split_n);
            TVectorD rest_part(n - split_n);

            // Copy the first 100 elements
            for (int i = 0; i < split_n; ++i) {
                first_part[i] = (*vector)[i];
            }

            // Copy the rest of the elements
            for (int i = split_n; i < n; ++i) {
                rest_part[i - split_n] = (*vector)[i];
            }

            // Write to the respective output files
            output_file_first->cd();
            first_part.Write(obj_name.c_str());

            output_file_second->cd();
            rest_part.Write(obj_name.c_str());
        }
        // Process TMatrixD objects
        else if (auto matrix = dynamic_cast<TMatrixD*>(obj)) {
            int n_rows = matrix->GetNrows();
            int n_cols = matrix->GetNcols();
            int split_rows = std::min(split_point, n_rows);
            int split_cols = std::min(split_point, n_cols);

            // Ensure the matrix is square for covariance matrices
            if (n_rows != n_cols) {
                std::cerr << "Warning: " << obj_name << " is not a square matrix (" 
                          << n_rows << "x" << n_cols << "). Skipping." << std::endl;
                continue;
            }

            // Create two new TMatrixD objects
            TMatrixD first_part(split_rows, split_cols);
            TMatrixD rest_part(n_rows - split_rows, n_cols - split_cols);

            // Fill the first 100x100 submatrix
            for (int i = 0; i < split_rows; ++i) {
                for (int j = 0; j < split_cols; ++j) {
                    first_part(i, j) = (*matrix)(i, j);
                }
            }

            // Fill the rest of the submatrix
            for (int i = split_rows; i < n_rows; ++i) {
                for (int j = split_cols; j < n_cols; ++j) {
                    rest_part(i - split_rows, j - split_cols) = (*matrix)(i, j);
                }
            }

            // Write to the respective output files using the original object name
            output_file_first->cd();
            first_part.Write(obj_name.c_str());

            output_file_second->cd();
            rest_part.Write(obj_name.c_str());
        }
        // Process TObjArray objects
        else if (auto array = dynamic_cast<TObjArray*>(obj)) {
            int n = array->GetEntries();
            int split_n = std::min(split_point, n);

            // Create two new TObjArray objects
            TObjArray* first_part = new TObjArray(split_n);
            TObjArray* rest_part = new TObjArray(n - split_n);

            // Copy the first 100 entries
            for (int i = 0; i < split_n; ++i) {
                TObject* item = array->At(i);
                first_part->Add(item);
            }

            // Copy the rest of the entries
            for (int i = split_n; i < n; ++i) {
                TObject* item = array->At(i);
                rest_part->Add(item);
            }

            // Write to the respective output files
            output_file_first->cd();
            first_part->Write(obj_name.c_str(), TObject::kSingleKey);

            output_file_second->cd();
            rest_part->Write(obj_name.c_str(), TObject::kSingleKey);

            // Clean up
            delete first_part;
            delete rest_part;
        }
        else {
            std::cout << "Skipping " << obj_name << ": unsupported object type" << std::endl;
            continue;
        }
    }

    // Close all files
    input_file->Close();
    output_file_first->Close();
    output_file_second->Close();

    std::cout << "Splitting completed successfully!" << std::endl;
}

