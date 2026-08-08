#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>

#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TFrame.h"

#include "TH1.h"
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <algorithm>

#include "TTree.h"
#include "TFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "Riostream.h"

#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TKey.h"

#include <yaml-cpp/yaml.h>

using namespace std;

string intToString(int t)
{
  std::string ch;
  ostringstream outs;
  outs << t;                                                  // Convert int into a string.
  ch = outs.str();

  return ch;
}

std::vector<std::string> SortAlphabetically(const std::vector<std::string>& input)
{
    std::vector<std::string> result = input;
    std::sort(result.begin(), result.end());
    return result;
}

int main(int argc, char* argv[])
{
    const std::string yamlFile = argv[1];

    YAML::Node root = YAML::LoadFile(yamlFile);
    YAML::Node InfoHdr = root["DialInformation"];

    std::vector<std::string> dial_names;
    for(const auto& n : InfoHdr["DialNames"])
    {
        dial_names.push_back(n.as<std::string>());
    }
    dial_names = SortAlphabetically(dial_names);
    int Ndials = dial_names.size();

    std::vector<std::string> sample_names;
    for(const auto& n : InfoHdr["SampleNames"])
    {
        sample_names.push_back(n.as<std::string>());
    }
    int Nsamples = sample_names.size();

    std::vector<double> varvalues;
    for(const auto& n : InfoHdr["Variations"])
    {
        varvalues.push_back(n.as<double>());
    }

    int Nbins = InfoHdr["Nbins"].as<int>();

    std::string Infile_location = "/Users/sjoshi/Work/DUNE/FD detector uncertainties/Energy scale systematics"
                                  "/Templates/Advanced_templates/";

    ifstream file;
    std::string filetxt = Infile_location + "filenames.txt";
    file.open(filetxt.c_str());

    std::vector<std::string> Infile_names;
    std::string Infile_name;

    while(getline(file,Infile_name))
    {
        Infile_names.push_back(Infile_name);
    }

    std::string Outfile_location = "/Users/sjoshi/Work/DUNE/FD detector uncertainties/Energy scale systematics/Splines/";
    std::string Outfile_name = "FDdetector_EScale_dials_v3.root";

    TFile *Outfile = new TFile((Outfile_location+Outfile_name).c_str(), "RECREATE");

    TObjArray *dialarr[Ndials];
    int Ngraphs = Nsamples * Nbins;
    TGraph *g1_sample[Ndials][Ngraphs];

    std::string g1_name, g1_title;

    Outfile->cd();

    for(int j = 0; j < Ndials; j++)
    {
        dialarr[j] = new TObjArray();
        dialarr[j]->SetName((dial_names[j]+"_graph").c_str());
        int c0 = 0;

        for(int k = 0; k < Nsamples; k++)
        {
            for(int l = 0; l < Nbins; l++)
            {
                g1_name = "Bin" + intToString(l+1) + "_Sample" + intToString(k+1);
                g1_title = sample_names[k] + " | " + dial_names[j] + " | Bin " + intToString(l+1) + "; #sigma; weight";

                g1_sample[j][c0] = new TGraph();
                g1_sample[j][c0]->SetName(g1_name.c_str());
                g1_sample[j][c0]->SetTitle(g1_title.c_str());
                g1_sample[j][c0]->SetMarkerStyle(8);
                g1_sample[j][c0]->SetMarkerSize(1.3);
                g1_sample[j][c0]->SetMarkerColor(kBlue);
                g1_sample[j][c0]->SetLineWidth(2);

                dialarr[j]->Add(g1_sample[j][c0]);
                c0++;
            }
        }
    }

    // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

    std::vector<std::string> tempdir_names;
    int Ndir = 0;
    std::string control_dirname = "variation_0";

    Outfile->cd();

    for(int j = 0; j < Ndials; j++)      //Ndials
    {
        TFile *template_file = TFile::Open((Infile_location+Infile_names[j]).c_str());
        int c1 = 0;

        if(j == 0)
        {
            TIter next(template_file->GetListOfKeys());
            TKey* key;

            while((key = (TKey*)next()))
            {
                if(std::string(key->GetClassName()) == "TDirectoryFile")
                {
                    tempdir_names.push_back(key->GetName());
                  //  cout << key->GetName() << endl;
                }
            }

            Ndir = tempdir_names.size();
        }

        for(int k = 0; k < Nsamples; k++)        //Nsamples
        {
            for(int l = 0; l < Nbins; l++)        //Nbins
            {
                cout << sample_names[k] << " " << dial_names[j] << " " << l << endl;

                TH1D *control_hist = template_file->Get<TH1D>(Form("%s/%s", control_dirname.c_str(), ("h_"+sample_names[k]).c_str()));
                double control_binvalue = control_hist->GetBinContent(l+1);

                for(int m = 0; m < Ndir; m++)           //Ndir
                {
                    TH1D *varhist = template_file->Get<TH1D>(Form("%s/%s", tempdir_names[m].c_str(), ("h_"+sample_names[k]).c_str()));
                    double varbinvalue = varhist->GetBinContent(l+1);

                    double weight;

                    if(control_binvalue == 0)
                    {
                        weight = 1;
                    }
                    else
                    {
                        weight = varbinvalue/control_binvalue;
                    }

                    g1_sample[j][c1]->SetPoint(m, varvalues[m], weight);

                    //  cout << c1 << " " << varbinvalue << " " << control_binvalue << " " << weight << endl;
                }

                c1++;
            }
        }

        template_file->Close();
    }

    for(int i = 0; i < Ndials; i++)
    {
        Outfile->WriteObject(dialarr[i], (dial_names[i]+"_graph").c_str());
    }
    Outfile->Close();

    return 0;

}
