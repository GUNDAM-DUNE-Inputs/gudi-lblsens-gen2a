#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TColor.h>
#include <TAxis.h>
#include <TPaveText.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <cctype>

namespace {

    // =========================================================
    // Configuration
    // =========================================================
    static const int NMODES = 14;

    // Set to true to use custom reco-energy binning separate for numu and nue selections.
    static const bool kUseCustomRecoBinning = true;

    enum class SelectionType {
        NuMu,
        NuE
    };

    enum class BeamMode {
        Unknown,
        FHC,
        RHC
    };

    // ---------------------------------------------------------
    // Custom reco-energy bin edges for numu selection.
    // NOTE: first edge must be > 0 if log-x is used.
    // ---------------------------------------------------------
    const std::vector<double> kRecoBinEdgesNuMu = {
        0.0, 0.5, 1.0, 1.25, 1.5,
        1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75,
        4.0, 5.0, 6.0, 10.0
    };

    // ---------------------------------------------------------
    // Custom reco-energy bin edges for nue selection.
    // NOTE: first edge must be > 0 if log-x is used.
    // ---------------------------------------------------------
    const std::vector<double> kRecoBinEdgesNuE = {
        0.0, 0.5, 1.0, 1.25, 1.5,
        1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75,
        4.0, 5.0, 6.0, 10.0
    };

    // Fallback uniform reco-energy binning if kUseCustomRecoBinning = false.
    static const int    kRecoNBinsUniform = 50;
    static const double kRecoMinUniform   = 0.0;
    static const double kRecoMaxUniform   = 10.0;

    // Truth-energy binning.
    static const int    kTrueNBins = 50;
    static const double kTrueMin   = 0.0;
    static const double kTrueMax   = 10.0;

    const char* modeTitle[NMODES] = {
        "QE", "RES", "DIS", "COH", "COHEL", "NuEl", "IMDAnn",
        "IBD", "GlRES", "AnuGam", "MEC", "DIFF", "kEM", "kWeakMix"
    };

    int modeColor[NMODES] = {
        kAzure+1, kOrange+7, kGreen+2, kMagenta+1, kCyan+1, kYellow+2,
        kRed+1, kBlue-7, kSpring+5, kViolet+1, kPink+7, kTeal+3,
        kGray+2, kOrange-3
    };

    // =========================================================
    // Category container
    // =========================================================
    struct Category {
        std::string name;
        std::string label;

        SelectionType selectionType;

        // MC 1D histograms stacked by interaction mode.
        std::vector<TH1D*> recoE_mc;
        std::vector<TH1D*> trueE_mc;

        // Data 1D histograms, inclusive.
        TH1D* recoE_data = nullptr;
        TH1D* trueE_data = nullptr;

        // Inclusive 2D histograms: x = reco energy, y = true energy.
        TH2D* recoTrue2D_mc   = nullptr;
        TH2D* recoTrue2D_data = nullptr;

        // Inclusive weighted event rates.
        double rate_mc   = 0.0;
        double rate_data = 0.0;

        // Inclusive raw entry counts.
        Long64_t n_mc   = 0;
        Long64_t n_data = 0;

        // Per-mode weighted event rates.
        std::vector<double> rate_mc_mode;
        std::vector<double> rate_data_mode;

        // Per-mode raw entry counts.
        std::vector<Long64_t> n_mc_mode;
        std::vector<Long64_t> n_data_mode;
    };

    // =========================================================
    // String helpers
    // =========================================================
    std::string toLower(std::string s)
    {
        std::transform(
            s.begin(),
            s.end(),
            s.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); }
        );

        return s;
    }

    bool contains(const std::string& s, const std::string& token)
    {
        return s.find(token) != std::string::npos;
    }

    BeamMode getBeamModeFromDirName(const std::string& dirName)
    {
        const std::string s = toLower(dirName);

        // Common labels.
        if (contains(s, "fhc") ||
            contains(s, "forward horn current") ||
            contains(s, "forward_horn_current") ||
            contains(s, "numode") ||
            contains(s, "nu-mode") ||
            contains(s, "nu_mode")) {
            return BeamMode::FHC;
        }

        if (contains(s, "rhc") ||
            contains(s, "reverse horn current") ||
            contains(s, "reverse_horn_current") ||
            contains(s, "antinumode") ||
            contains(s, "antinu-mode") ||
            contains(s, "antinu_mode")) {
            return BeamMode::RHC;
        }

        return BeamMode::Unknown;
    }

    bool isNuMuSampleName(const std::string& dirName)
    {
        const std::string s = toLower(dirName);

        return contains(s, "#nu_{#mu}") ||
               contains(s, "#bar{#nu}_{#mu}") ||
               contains(s, "numu") ||
               contains(s, "nu_mu") ||
               contains(s, "nubar_mu") ||
               contains(s, "anti_numu") ||
               contains(s, "antinumu");
    }

    bool isNuESampleName(const std::string& dirName)
    {
        const std::string s = toLower(dirName);

        return contains(s, "#nu_{e}") ||
               contains(s, "#bar{#nu}_{e}") ||
               contains(s, "nue") ||
               contains(s, "nu_e") ||
               contains(s, "nubar_e") ||
               contains(s, "anti_nue") ||
               contains(s, "antinue");
    }

    // =========================================================
    // Histogram helpers
    // =========================================================
    void styleData(TH1D* h)
    {
        if (!h) return;

        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.9);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetLineWidth(2);
    }

    const std::vector<double>& getRecoBinEdges(SelectionType selectionType)
    {
        if (selectionType == SelectionType::NuMu) {
            return kRecoBinEdgesNuMu;
        }

        return kRecoBinEdgesNuE;
    }

    bool hasCustomRecoBinning(SelectionType selectionType)
    {
        const auto& edges = getRecoBinEdges(selectionType);
        return kUseCustomRecoBinning && edges.size() >= 2;
    }

    TH1D* makeRecoHist1D(const TString& name,
                         const TString& title,
                         SelectionType selectionType)
    {
        if (hasCustomRecoBinning(selectionType)) {
            const auto& edges = getRecoBinEdges(selectionType);

            return new TH1D(
                name,
                title,
                static_cast<int>(edges.size()) - 1,
                edges.data()
            );
        }

        return new TH1D(
            name,
            title,
            kRecoNBinsUniform,
            kRecoMinUniform,
            kRecoMaxUniform
        );
    }

    TH1D* makeTrueHist1D(const TString& name, const TString& title)
    {
        return new TH1D(
            name,
            title,
            kTrueNBins,
            kTrueMin,
            kTrueMax
        );
    }

    TH2D* makeRecoTrueHist2D(const TString& name,
                             const TString& title,
                             SelectionType selectionType)
    {
        if (hasCustomRecoBinning(selectionType)) {
            const auto& recoEdges = getRecoBinEdges(selectionType);

            return new TH2D(
                name,
                title,
                static_cast<int>(recoEdges.size()) - 1,
                recoEdges.data(),
                kTrueNBins,
                kTrueMin,
                kTrueMax
            );
        }

        return new TH2D(
            name,
            title,
            kRecoNBinsUniform,
            kRecoMinUniform,
            kRecoMaxUniform,
            kTrueNBins,
            kTrueMin,
            kTrueMax
        );
    }

    void initCategory(Category& c,
                      const std::string& name,
                      const std::string& label,
                      SelectionType selectionType)
    {
        c.name = name;
        c.label = label;
        c.selectionType = selectionType;

        c.recoE_mc.resize(NMODES, nullptr);
        c.trueE_mc.resize(NMODES, nullptr);

        c.rate_mc_mode.assign(NMODES, 0.0);
        c.rate_data_mode.assign(NMODES, 0.0);

        c.n_mc_mode.assign(NMODES, 0);
        c.n_data_mode.assign(NMODES, 0);

        for (int m = 0; m < NMODES; ++m) {
            c.recoE_mc[m] = makeRecoHist1D(
                Form("hRecoE_%s_mc_mode%d", name.c_str(), m),
                Form("Reco E_{#nu} by mode (%s);Reco E_{#nu} [GeV];Events", label.c_str()),
                c.selectionType
            );

            c.trueE_mc[m] = makeTrueHist1D(
                Form("hTrueE_%s_mc_mode%d", name.c_str(), m),
                Form("True E_{#nu} by mode (%s);True E_{#nu} [GeV];Events", label.c_str())
            );

            c.recoE_mc[m]->Sumw2();
            c.trueE_mc[m]->Sumw2();
        }

        c.recoE_data = makeRecoHist1D(
            Form("hRecoE_%s_data", name.c_str()),
            Form("Reco E_{#nu} (%s);Reco E_{#nu} [GeV];Events", label.c_str()),
            c.selectionType
        );

        c.trueE_data = makeTrueHist1D(
            Form("hTrueE_%s_data", name.c_str()),
            Form("True E_{#nu} (%s);True E_{#nu} [GeV];Events", label.c_str())
        );

        c.recoE_data->Sumw2();
        c.trueE_data->Sumw2();

        styleData(c.recoE_data);
        styleData(c.trueE_data);

        c.recoTrue2D_mc = makeRecoTrueHist2D(
            Form("hRecoTrue2D_%s_mc", name.c_str()),
            Form("Model: True E_{#nu} vs Reco E_{#nu} (%s);Reco E_{#nu} [GeV];True E_{#nu} [GeV]", label.c_str()),
            c.selectionType
        );

        c.recoTrue2D_data = makeRecoTrueHist2D(
            Form("hRecoTrue2D_%s_data", name.c_str()),
            Form("Data: True E_{#nu} vs Reco E_{#nu} (%s);Reco E_{#nu} [GeV];True E_{#nu} [GeV]", label.c_str()),
            c.selectionType
        );

        c.recoTrue2D_mc->Sumw2();
        c.recoTrue2D_data->Sumw2();
    }

    THStack* makeStack(const std::vector<TH1D*>& hs,
                       const char* name,
                       TLegend*& leg,
                       int& nAdded)
    {
        THStack* st = new THStack(name, "");

        leg = new TLegend(0.70, 0.38, 0.96, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);

        nAdded = 0;

        for (int m = 0; m < NMODES; ++m) {
            if (!hs[m]) continue;

            const double integral = hs[m]->Integral();

            if (integral <= 0.0) continue;

            hs[m]->SetFillColor(modeColor[m]);
            hs[m]->SetLineColor(kBlack);
            hs[m]->SetLineWidth(1);

            st->Add(hs[m]);
            leg->AddEntry(hs[m], modeTitle[m], "f");

            ++nAdded;
        }

        return st;
    }

    void drawEmptyPage(TCanvas& can,
                       const std::string& pdf,
                       const TString& title,
                       const TString& message)
    {
        can.Clear();
        can.SetLogx(false);
        can.SetLogz(false);
        can.SetRightMargin(0.10);

        TH1D frame("empty_frame", title, 1, 0.0, 1.0);
        frame.SetMinimum(0.0);
        frame.SetMaximum(1.0);
        frame.Draw();

        TPaveText box(0.18, 0.42, 0.82, 0.58, "NDC");
        box.SetFillColor(0);
        box.SetBorderSize(1);
        box.AddText(message);
        box.Draw();

        can.Print(pdf.c_str());
    }

    void drawStackWithData(TCanvas& can,
                           const std::string& pdf,
                           const std::vector<TH1D*>& mcHists,
                           TH1D* dataHist,
                           const TString& title,
                           bool drawData = true,
                           bool logX = false)
    {
        TLegend* leg = nullptr;
        int nAdded = 0;

        THStack* st = makeStack(
            mcHists,
            Form("st_%s", dataHist ? dataHist->GetName() : "nullData"),
            leg,
            nAdded
        );

        const bool hasData =
            drawData &&
            dataHist &&
            dataHist->Integral() > 0.0;

        if (nAdded == 0 && !hasData) {
            drawEmptyPage(
                can,
                pdf,
                title,
                "No entries passed the selection for this plot."
            );

            delete st;
            delete leg;
            return;
        }

        can.Clear();
        can.SetLogz(false);
        can.SetLogx(logX && kUseCustomRecoBinning);

        st->SetTitle(title);
        st->Draw("hist");

        if (logX && kUseCustomRecoBinning && st->GetXaxis()) {
            st->GetXaxis()->SetMoreLogLabels(true);
            st->GetXaxis()->SetNoExponent(true);
        }

        double maxY = st->GetMaximum();

        if (hasData && dataHist->GetMaximum() > maxY) {
            maxY = dataHist->GetMaximum();
        }

        if (maxY > 0.0) {
            st->SetMaximum(1.25 * maxY);
        }

        if (hasData) {
            dataHist->Draw("E1 same");
            leg->AddEntry(dataHist, "Data", "lep");
        }

        leg->Draw();
        can.Print(pdf.c_str());

        can.SetLogx(false);
        can.SetLogz(false);

        delete st;
        delete leg;
    }

    void draw2DPage(TCanvas& can,
                    const std::string& pdf,
                    TH2D* h2,
                    const TString& title,
                    bool logZ = false,
                    bool logX = false)
    {
        if (!h2) return;

        if (h2->Integral() <= 0.0) {
            drawEmptyPage(
                can,
                pdf,
                title,
                "No entries passed the selection for this 2D plot."
            );
            return;
        }

        can.Clear();
        can.SetRightMargin(0.14);
        can.SetLogz(logZ);
        can.SetLogx(logX && kUseCustomRecoBinning);

        h2->SetTitle(title);
        h2->Draw("COLZ");

        if (logX && kUseCustomRecoBinning && h2->GetXaxis()) {
            h2->GetXaxis()->SetMoreLogLabels(true);
            h2->GetXaxis()->SetNoExponent(true);
        }

        can.Print(pdf.c_str());

        can.SetLogx(false);
        can.SetLogz(false);
        can.SetRightMargin(0.10);
    }

    // =========================================================
    // Filling
    // =========================================================
    void fillOneCategory(Category& cat,
                         bool isData,
                         int mode,
                         double w,
                         double Ereco,
                         double Etrue)
    {
        // Inclusive event-rate counters.
        if (!isData) {
            cat.rate_mc += w;
            cat.n_mc += 1;
        } else {
            cat.rate_data += w;
            cat.n_data += 1;
        }

        // Per-mode event-rate counters.
        if (mode >= 0 && mode < NMODES) {
            if (!isData) {
                cat.rate_mc_mode[mode] += w;
                cat.n_mc_mode[mode] += 1;
            } else {
                cat.rate_data_mode[mode] += w;
                cat.n_data_mode[mode] += 1;
            }
        }

        // Fill histograms.
        if (!isData) {
            // 1D MC histograms are stacked by interaction mode.
            if (mode >= 0 && mode < NMODES) {
                if (std::isfinite(Ereco)) {
                    cat.recoE_mc[mode]->Fill(Ereco, w);
                }

                if (std::isfinite(Etrue)) {
                    cat.trueE_mc[mode]->Fill(Etrue, w);
                }
            }

            // 2D MC histogram is inclusive.
            if (std::isfinite(Ereco) && std::isfinite(Etrue)) {
                cat.recoTrue2D_mc->Fill(Ereco, Etrue, w);
            }
        }
        else {
            if (std::isfinite(Ereco)) {
                cat.recoE_data->Fill(Ereco, w);
            }

            if (std::isfinite(Etrue)) {
                cat.trueE_data->Fill(Etrue, w);
            }

            if (std::isfinite(Ereco) && std::isfinite(Etrue)) {
                cat.recoTrue2D_data->Fill(Ereco, Etrue, w);
            }
        }
    }

    void fillFromBaseDir(TDirectory* baseDir,
                         Category& fhc_numu,
                         Category& fhc_nue,
                         Category& rhc_numu,
                         Category& rhc_nue,
                         bool isData)
    {
        if (!baseDir) return;

        std::cout << "\nScanning base directory: " << baseDir->GetName()
                  << "   isData = " << isData << "\n";

        TIter nextSubfolder(baseDir->GetListOfKeys());
        TKey* keySubfolder = nullptr;

        while ((keySubfolder = static_cast<TKey*>(nextSubfolder()))) {
            if (std::strcmp(keySubfolder->GetClassName(), "TDirectoryFile") != 0) {
                continue;
            }

            TDirectory* subDir =
                dynamic_cast<TDirectory*>(baseDir->Get(keySubfolder->GetName()));

            if (!subDir) continue;

            std::string subDirName = subDir->GetName();

            std::cout << "  Sample directory: [" << subDirName << "]\n";

            const bool isNuMu = isNuMuSampleName(subDirName);
            const bool isNuE  = isNuESampleName(subDirName);

            if (!isNuMu && !isNuE) {
                std::cout << "    skipped: cannot identify numu/nue selection\n";
                continue;
            }

            const BeamMode beamMode = getBeamModeFromDirName(subDirName);

            if (beamMode == BeamMode::Unknown) {
                std::cout << "    skipped: cannot identify FHC/RHC beam mode\n";
                continue;
            }

            TTree* tree = dynamic_cast<TTree*>(subDir->Get("events_TTree"));
            if (!tree) {
                std::cout << "    skipped: no events_TTree\n";
                continue;
            }

            TLeaf* leafMode      = tree->GetLeaf("mode");
            TLeaf* leafGenWeight = tree->GetLeaf("eventWeight");

            TLeaf* leafEnuRecoEle = tree->GetLeaf("Ev_reco_nue");
            TLeaf* leafEnuRecoMu  = tree->GetLeaf("Ev_reco_numu");
            TLeaf* leafEnuTruth   = tree->GetLeaf("Ev");

            if (!leafGenWeight) {
                std::cerr << "    skipped: eventWeight leaf is missing\n";
                continue;
            }

            std::cout << "    classified as: "
                      << (beamMode == BeamMode::FHC ? "FHC" : "RHC")
                      << " "
                      << (isNuMu ? "numu" : "")
                      << (isNuE  ? "nue"  : "")
                      << "\n";

            std::cout << "    leaves: "
                      << " mode=" << (leafMode ? "yes" : "no")
                      << " eventWeight=" << (leafGenWeight ? "yes" : "no")
                      << " Ev_reco_numu=" << (leafEnuRecoMu ? "yes" : "no")
                      << " Ev_reco_nue=" << (leafEnuRecoEle ? "yes" : "no")
                      << " Ev=" << (leafEnuTruth ? "yes" : "no")
                      << "\n";

            const Long64_t nEntries = tree->GetEntries();

            std::cout << "    entries = " << nEntries << "\n";

            Category* cat = nullptr;
            SelectionType selectionType = SelectionType::NuMu;

            if (beamMode == BeamMode::FHC && isNuMu) {
                cat = &fhc_numu;
                selectionType = SelectionType::NuMu;
            }
            else if (beamMode == BeamMode::FHC && isNuE) {
                cat = &fhc_nue;
                selectionType = SelectionType::NuE;
            }
            else if (beamMode == BeamMode::RHC && isNuMu) {
                cat = &rhc_numu;
                selectionType = SelectionType::NuMu;
            }
            else if (beamMode == BeamMode::RHC && isNuE) {
                cat = &rhc_nue;
                selectionType = SelectionType::NuE;
            }

            if (!cat) {
                std::cout << "    skipped: no matching category\n";
                continue;
            }

            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);

                int mode = -1;

                if (leafMode) {
                    mode = static_cast<int>(leafMode->GetValue(0));
                }

                const double w =
                    static_cast<double>(leafGenWeight->GetValue(0));

                double Etrue = NAN;

                if (leafEnuTruth) {
                    Etrue = static_cast<double>(leafEnuTruth->GetValue(0));
                }

                double Ereco = NAN;

                if (selectionType == SelectionType::NuMu) {
                    if (leafEnuRecoMu) {
                        Ereco = static_cast<double>(leafEnuRecoMu->GetValue(0));
                    }
                }
                else if (selectionType == SelectionType::NuE) {
                    if (leafEnuRecoEle) {
                        Ereco = static_cast<double>(leafEnuRecoEle->GetValue(0));
                    }
                }

                fillOneCategory(
                    *cat,
                    isData,
                    mode,
                    w,
                    Ereco,
                    Etrue
                );
            }
        }
    }

    // =========================================================
    // Printing and rate summaries
    // =========================================================
    void printCategoryPDF(const Category& c)
    {
        TCanvas can("can", "can", 1000, 700);

        std::string pdf = c.name + ".pdf";

        can.Print((pdf + "[").c_str());

        // 1D reco-energy stack by interaction mode, with data overlay.
        drawStackWithData(
            can,
            pdf,
            c.recoE_mc,
            c.recoE_data,
            Form("Reco E_{#nu} (%s);Reco E_{#nu} [GeV];Events", c.label.c_str()),
            true,
            false
        );

        // 1D true-energy stack by interaction mode, with data overlay if available.
        bool hasTrueData = (c.trueE_data && c.trueE_data->Integral() > 0.0);

        drawStackWithData(
            can,
            pdf,
            c.trueE_mc,
            c.trueE_data,
            Form("True E_{#nu} (%s);True E_{#nu} [GeV];Events", c.label.c_str()),
            hasTrueData,
            false
        );

        // 2D reco-vs-true energy plots.
        draw2DPage(
            can,
            pdf,
            c.recoTrue2D_mc,
            Form("Model: True E_{#nu} vs Reco E_{#nu} (%s);Reco E_{#nu} [GeV];True E_{#nu} [GeV]", c.label.c_str()),
            false,
            false
        );

        can.Print((pdf + "]").c_str());

        std::cout << "Wrote " << pdf << "\n";
    }

    void printEventRatesPerMode(const Category& c)
    {
        std::cout << "\n------------------------------------------------------------\n";
        std::cout << "Event rates per interaction mode for " << c.name
                  << " [" << c.label << "]\n";
        std::cout << "------------------------------------------------------------\n";

        std::cout << std::fixed << std::setprecision(6);

        std::cout << std::setw(14) << std::left << "Mode"
                  << " | " << std::setw(16) << "Model rate"
                  << " | " << std::setw(16) << "Data rate"
                  << " | " << std::setw(12) << "MC entries"
                  << " | " << std::setw(12) << "Data entries"
                  << "\n";

        std::cout << "------------------------------------------------------------\n";

        double sumModel = 0.0;
        double sumData = 0.0;

        Long64_t sumModelEntries = 0;
        Long64_t sumDataEntries = 0;

        for (int m = 0; m < NMODES; ++m) {
            const double rMC = c.rate_mc_mode[m];
            const double rData = c.rate_data_mode[m];

            const Long64_t nMC = c.n_mc_mode[m];
            const Long64_t nData = c.n_data_mode[m];

            if (rMC == 0.0 && rData == 0.0 && nMC == 0 && nData == 0) {
                continue;
            }

            sumModel += rMC;
            sumData += rData;

            sumModelEntries += nMC;
            sumDataEntries += nData;

            std::cout << std::setw(14) << std::left << modeTitle[m]
                      << " | " << std::setw(16) << rMC
                      << " | " << std::setw(16) << rData
                      << " | " << std::setw(12) << nMC
                      << " | " << std::setw(12) << nData
                      << "\n";
        }

        std::cout << "------------------------------------------------------------\n";
        std::cout << std::setw(14) << std::left << "TOTAL MODES"
                  << " | " << std::setw(16) << sumModel
                  << " | " << std::setw(16) << sumData
                  << " | " << std::setw(12) << sumModelEntries
                  << " | " << std::setw(12) << sumDataEntries
                  << "\n";

        std::cout << std::setw(14) << std::left << "INCLUSIVE"
                  << " | " << std::setw(16) << c.rate_mc
                  << " | " << std::setw(16) << c.rate_data
                  << " | " << std::setw(12) << c.n_mc
                  << " | " << std::setw(12) << c.n_data
                  << "\n";

        if (std::fabs(sumModel - c.rate_mc) > 1e-8 ||
            std::fabs(sumData - c.rate_data) > 1e-8) {
            std::cout << "  Note: TOTAL MODES may differ from INCLUSIVE if some entries "
                      << "have missing or invalid mode.\n";
        }
    }

    void printBeamCategoryRates(const Category& fhc_numu,
                                const Category& fhc_nue,
                                const Category& rhc_numu,
                                const Category& rhc_nue)
    {
        auto printOne = [](const Category& c) {
            std::cout << std::fixed << std::setprecision(6)
                      << "  " << std::setw(24) << std::left << c.name
                      << " | model rate = " << std::setw(16) << c.rate_mc
                      << " | data rate = " << std::setw(16) << c.rate_data
                      << " | model entries = " << std::setw(10) << c.n_mc
                      << " | data entries = " << c.n_data
                      << "\n";
        };

        std::cout << "\n============================================================\n";
        std::cout << "Beam-mode category event rates\n";
        std::cout << "============================================================\n";

        printOne(fhc_numu);
        printOne(fhc_nue);
        printOne(rhc_numu);
        printOne(rhc_nue);

        std::cout << "\n============================================================\n";
        std::cout << "Per-interaction-mode event rates\n";
        std::cout << "============================================================\n";

        printEventRatesPerMode(fhc_numu);
        printEventRatesPerMode(fhc_nue);
        printEventRatesPerMode(rhc_numu);
        printEventRatesPerMode(rhc_nue);

        std::cout << "\n============================================================\n\n";
    }

    void printRecoBinning()
    {
        auto printEdges = [](const std::string& label,
                             SelectionType selectionType)
        {
            std::cout << "  " << label << ": ";

            if (hasCustomRecoBinning(selectionType)) {
                const auto& edges = getRecoBinEdges(selectionType);

                std::cout << "custom binning with edges:\n    ";

                for (size_t i = 0; i < edges.size(); ++i) {
                    std::cout << edges[i];

                    if (i + 1 < edges.size()) {
                        std::cout << ", ";
                    }
                }

                std::cout << "\n";
            } else {
                std::cout << "uniform binning: "
                          << kRecoNBinsUniform << " bins from "
                          << kRecoMinUniform << " to "
                          << kRecoMaxUniform << "\n";
            }
        };

        std::cout << "\nReco binning configuration:\n";

        printEdges("numu selection", SelectionType::NuMu);
        printEdges("nue  selection", SelectionType::NuE);

        std::cout << "\n";
    }

} // namespace

double histIntegral(const TH1D* h, bool includeUnderOverflow = false)
{
    if (!h) return 0.0;

    if (includeUnderOverflow) {
        return h->Integral(0, h->GetNbinsX() + 1);
    }

    return h->Integral(1, h->GetNbinsX());
}


void print1DPlotIntegrals(const Category& c)
{
    std::cout << "\n============================================================\n";
    std::cout << "1D plot integrals for " << c.name
              << " [" << c.label << "]\n";
    std::cout << "============================================================\n";

    std::cout << std::fixed << std::setprecision(6);

    double recoMCTotal = 0.0;
    double trueMCTotal = 0.0;

    std::cout << "\nReco-energy MC integrals by interaction mode:\n";
    std::cout << "------------------------------------------------------------\n";

    for (int m = 0; m < NMODES; ++m) {
        const double integral = histIntegral(c.recoE_mc[m]);

        if (integral == 0.0) continue;

        recoMCTotal += integral;

        std::cout << "  " << std::setw(14) << std::left << modeTitle[m]
                  << " : " << integral << "\n";
    }

    std::cout << "------------------------------------------------------------\n";
    std::cout << "  " << std::setw(14) << std::left << "TOTAL MC"
              << " : " << recoMCTotal << "\n";

    std::cout << "  " << std::setw(14) << std::left << "DATA"
              << " : " << histIntegral(c.recoE_data) << "\n";


    std::cout << "\nTrue-energy MC integrals by interaction mode:\n";
    std::cout << "------------------------------------------------------------\n";

    for (int m = 0; m < NMODES; ++m) {
        const double integral = histIntegral(c.trueE_mc[m]);

        if (integral == 0.0) continue;

        trueMCTotal += integral;

        std::cout << "  " << std::setw(14) << std::left << modeTitle[m]
                  << " : " << integral << "\n";
    }

    std::cout << "------------------------------------------------------------\n";
    std::cout << "  " << std::setw(14) << std::left << "TOTAL MC"
              << " : " << trueMCTotal << "\n";

    std::cout << "  " << std::setw(14) << std::left << "DATA"
              << " : " << histIntegral(c.trueE_data) << "\n";

    std::cout << "============================================================\n\n";
}

// =========================================================
// Main entry point
// =========================================================
void getERPlots()
{
    gStyle->SetOptStat(0);

    const char* inputFileName =
        "/gpfs/projects/McGrewGroup/uyevarou/Work/lbl/lbl_profiling/gudi-lbl-tdr25a/outputs/gundamFitter_config_DUNE_Asimov_DryRun.root";

    printRecoBinning();

    TFile* f = TFile::Open(inputFileName, "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file: " << inputFileName << std::endl;
        return;
    }

    // MC model.
    TDirectory* modelDir =
        dynamic_cast<TDirectory*>(f->Get("FitterEngine/preFit/model"));

    if (!modelDir) {
        std::cerr << "Could not find FitterEngine/preFit/model" << std::endl;
        f->Close();
        return;
    }

    // Data / Asimov data.
    TDirectory* dataDir =
        dynamic_cast<TDirectory*>(f->Get("FitterEngine/preFit/data"));

    if (!dataDir) {
        std::cerr << "Could not find FitterEngine/preFit/data" << std::endl;
        f->Close();
        return;
    }

    // Four requested beam-mode categories.
    Category fhc_numu;
    Category fhc_nue;
    Category rhc_numu;
    Category rhc_nue;

    initCategory(
        fhc_numu,
        "modeStacks_FHC_numu",
        "FHC #nu_{#mu} selection, CC+NC",
        SelectionType::NuMu
    );

    initCategory(
        fhc_nue,
        "modeStacks_FHC_nue",
        "FHC #nu_{e} selection, CC+NC",
        SelectionType::NuE
    );

    initCategory(
        rhc_numu,
        "modeStacks_RHC_numu",
        "RHC #nu_{#mu} selection, CC+NC",
        SelectionType::NuMu
    );

    initCategory(
        rhc_nue,
        "modeStacks_RHC_nue",
        "RHC #nu_{e} selection, CC+NC",
        SelectionType::NuE
    );

    // Fill model.
    fillFromBaseDir(
        modelDir,
        fhc_numu,
        fhc_nue,
        rhc_numu,
        rhc_nue,
        false
    );

    // Fill data / Asimov.
    fillFromBaseDir(
        dataDir,
        fhc_numu,
        fhc_nue,
        rhc_numu,
        rhc_nue,
        true
    );

    printBeamCategoryRates(
        fhc_numu,
        fhc_nue,
        rhc_numu,
        rhc_nue
    );

    // Print one PDF per beam-mode category.
    printCategoryPDF(fhc_numu);
    printCategoryPDF(fhc_nue);
    printCategoryPDF(rhc_numu);
    printCategoryPDF(rhc_nue);

    //print1DPlotIntegrals(fhc_numu);
    //print1DPlotIntegrals(fhc_nue);
    //print1DPlotIntegrals(rhc_numu);
    //print1DPlotIntegrals(rhc_nue);

    f->Close();
}