#include <cstdlib>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <set>
#include <list>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TStopwatch.h"
#include "TString.h"

#include "HistogramManager.h"
#include "ConfigLoader.h"


int main(int argc, char** argv)
{
    //argv[1] = yaml file, argv[2] = output file location

    HistogramConfig cfg = LoadHistogramConfig(argv[1]);
    HistogramManager histMgr(cfg, argv[2]);

    const int Ndials      = cfg.dialNames.size();
    const int Nvariations = cfg.variationValues.size();
    const int Nsamples    = cfg.sampleNames.size();

    std::vector<std::string> MakeFDSamples();
    std::vector<std::string> sampleInput = MakeFDSamples();

    // Important!--- Nsamples == sampleInput.size() && (same samples order)

    std::string Infile_location = "/Users/sjoshi/Work/DUNE/FD detector uncertainties/Energy scale systematics/";      // File location
    std::string Infile_name = "gundamFitter_config_DUNE_Asimov_DryRun.root";
    TFile *input_file = new TFile((Infile_location+Infile_name).c_str());

    double eRecProxy, eRecProxy_shift, Ehadron_sum;
    eRecProxy = eRecProxy_shift = Ehadron_sum = 0;

    double es_totalshift, es_hadronshift, es_muonshift, es_neutronshift, es_emshift, resolution_shift;
    es_totalshift = es_hadronshift = es_muonshift = es_neutronshift = es_emshift = resolution_shift = 0;

    std::vector<double> es_totalp0 = {0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_totalp1 = {0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_totalp2 = {0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_hadronp0 = {0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_hadronp1 = {0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_hadronp2 = {0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_muonp0 = {0, 0, 0, 0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_muonp1 = {0, 0, 0, 0, 0, 0, 0, 0.005, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_muonp2 = {0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_neutronp0 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2, 0, 0, 0, 0, 0, 0};
    std::vector<double> es_neutronp1 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0};
    std::vector<double> es_neutronp2 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0, 0, 0, 0};
    std::vector<double> es_emp0 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.025, 0, 0, 0};
    std::vector<double> es_emp1 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.025, 0, 0};
    std::vector<double> es_emp2 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.025, 0};
    std::vector<double> det_resolution = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1};


    for(int idial=0; idial<Ndials; idial++)                  //Ndials
    {
        for(int isample=0; isample<Nsamples; isample++)                               //Nsamples
        {
            TTree *input_tree = input_file->Get<TTree>(Form("FitterEngine/preFit/model/%s/events_TTree", sampleInput[isample].c_str()));

            std::cout << sampleInput[isample].c_str() << " sample to be processed..." << std::endl;

            for(int ientry=0; ientry<input_tree->GetEntries(); ientry++)                       //input_tree->GetEntries()
            {
                input_tree->GetEntry(ientry);

                if (ientry % 5000 == 0)
                {
                    std::cout <<  "Progress " << 100. * ((double)ientry) / input_tree->GetEntries() << "%" << std::endl;
                }

                double LepE = input_tree->GetLeaf("Leaves/LepE")->GetValue();
                double eP = input_tree->GetLeaf("Leaves/eP")->GetValue();
                double eN = input_tree->GetLeaf("Leaves/eN")->GetValue();
                double ePip = input_tree->GetLeaf("Leaves/ePip")->GetValue();
                double ePim = input_tree->GetLeaf("Leaves/ePim")->GetValue();
                double ePi0 = input_tree->GetLeaf("Leaves/ePi0")->GetValue();
                double eOther = input_tree->GetLeaf("Leaves/eOther")->GetValue();
                int nipi0 = input_tree->GetLeaf("Leaves/nipi0")->GetValue();
                double Ev = input_tree->GetLeaf("Leaves/Ev")->GetValue();
                int isCC = input_tree->GetLeaf("Leaves/isCC")->GetValue();
                int nuPDG = input_tree->GetLeaf("Leaves/nuPDG")->GetValue();

                double eventWeight = input_tree->GetLeaf("Event/eventWeight")->GetValue();

                for(int ivar=0; ivar<Nvariations; ivar++)             //Nvariations
                {
                    double sigma = cfg.variationValues[ivar];
                    eRecProxy = LepE + eP + ePip + ePim + ePi0 + (0.135 * nipi0) + eOther;
                    Ehadron_sum = eP + ePip + ePim;

                    es_totalshift = 0;
                    es_muonshift = 0;

                    if(eRecProxy >= 10)
                        continue;

                    if(abs(nuPDG) == 14)              //Include only (anti-)nue events
                        continue;

                    // Total energy scale uncertainty

                    if(isCC == 1)
                    {
                        if(abs(nuPDG) == 14)
                        {
                            es_totalshift = ((eRecProxy - LepE) * es_totalp0[idial] * sigma)
                            + ((eRecProxy - LepE) * es_totalp1[idial] * sigma * sqrt(eRecProxy - LepE))
                            + ((eRecProxy - LepE) * es_totalp2[idial] * sigma * (1/sqrt(eRecProxy - LepE + 0.1)));
                        }
                        else if(abs(nuPDG) == 12)
                        {
                            es_totalshift = (eRecProxy * es_totalp0[idial] * sigma)
                            + (eRecProxy * es_totalp1[idial] * sigma * sqrt(eRecProxy))
                            + (eRecProxy * es_totalp2[idial] * sigma * (1/sqrt(eRecProxy + 0.1)));
                        }
                    }

                    // EM particle response energy scale uncertainty

                    es_emshift = (ePi0 * es_emp0[idial] * sigma)
                                 + (ePi0 * es_emp1[idial] * sigma * sqrt(ePi0))
                                 + (ePi0 * es_emp2[idial] * sigma * (1/sqrt(ePi0 + 0.1)));

                    if(isCC == 1 && abs(nuPDG) == 12)
                    {
                        es_emshift += (LepE * es_emp0[idial] * sigma)
                        + (LepE * es_emp1[idial] * sigma * sqrt(LepE))
                        + (LepE * es_emp2[idial] * sigma * (1/sqrt(LepE + 0.1)));

                    }

                    // Charged hadrons particle response energy scale uncertainty

                    es_hadronshift = (Ehadron_sum * es_hadronp0[idial] * sigma)
                                     + (Ehadron_sum * es_hadronp1[idial] * sigma * sqrt(Ehadron_sum))
                                     + (Ehadron_sum * es_hadronp2[idial] * sigma * (1/sqrt(Ehadron_sum + 0.1)));

                    // Neutrons particle response energy scale uncertainty

                    es_neutronshift = (eN * es_neutronp0[idial] * sigma)
                                      + (eN * es_neutronp1[idial] * sigma * sqrt(eN))
                                      + (eN * es_neutronp2[idial] * sigma * (1/sqrt(eN + 0.1)));

                    // Muon particle response energy scale uncertainty

                    if(isCC == 1 && abs(nuPDG) == 14)
                    {
                        es_muonshift = (LepE * es_muonp0[idial] * sigma)
                        + (LepE * es_muonp1[idial] * sigma * sqrt(LepE))
                            + (LepE * es_muonp2[idial] * sigma * (1/sqrt(LepE + 0.1)));
                    }

                    // Effect of resolution on Erec

                    resolution_shift = (Ev - eRecProxy) * det_resolution[idial] * sigma;

                    // // // // // // // // // // // // // // // // // // // // // // // //

                    eRecProxy_shift = es_totalshift + es_emshift + es_hadronshift + es_neutronshift + es_muonshift + resolution_shift;       //+ eRecProxy for true 'shifted'

                    eRecProxy += eRecProxy_shift;

                    //    std::cout << i << " | sigma: " << sigma << " | particle energies LepE, eN, Ehadron_sum: "
                    //              << LepE << " " << eN << " " << Ehadron_sum << "   |  ERP shift: "
                    //              << eRecProxy_shift << "    ERP final: " << eRecProxy << "   | Bin index: " << BinIndex << std::endl;

                    histMgr.Fill(idial, ivar, isample, eRecProxy, eventWeight);

                }
            }
        }
    }

    histMgr.WriteAndClose();

    return 0;         

}

std::vector<std::string> MakeFDSamples()
{
    std::vector<std::string> FDSampleNames;

    FDSampleNames.push_back("FHC #nu_{#mu}_like CC ");
    FDSampleNames.push_back("FHC #nu_{e}_like CC ");
    FDSampleNames.push_back("RHC #nu_{#mu}_like CC ");
    FDSampleNames.push_back("RHC #nu_{e}_like CC ");

    return FDSampleNames;
}
