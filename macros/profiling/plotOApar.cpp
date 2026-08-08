#include <filesystem>
#include <regex>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>   // sort, min_element, max_element
#include <cctype>      // std::isalnum
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>

namespace fs = std::filesystem;
#include <limits>
static inline bool isZero(double x, double eps=1e-12) { return std::abs(x) < eps; }


// -------------------- Helpers --------------------

static std::string xTitleFor(const std::string& param) {
  if (param == "dcp" || param == "delta" || param == "deltaCP" || param == "dcp_rad")
    return "#delta_{CP} [rad]";
  if (param == "m21" || param == "dm21" || param == "dm2_21")
    return "#Delta m^{2}_{21} [eV^{2}]";
  if (param == "m32" || param == "dm31" || param == "dm32" || param == "dm2_31" || param == "dm2_32")
    return "#Delta m^{2}_{32} [eV^{2}]";
  if (param == "sin12")
    return "sin^{2}#theta_{12}";
  if (param == "sin13")
    return "sin^{2}#theta_{13}";
  if (param == "sin23")
    return "sin^{2}#theta_{23}";
  if (param == "MO")
    return "MO";
  return param; // fallback
}

// extract the parameter name, e.g. "sin12"
static std::string extractParamName(const std::string& fn) {
  const std::string marker = "_With_";
  auto pos = fn.find(marker);
  if (pos == std::string::npos) {
    throw std::runtime_error("No '_With_' in filename: " + fn);
  }
  size_t start = pos + marker.size();
  size_t end = start;
  while (end < fn.size() && std::isalnum(static_cast<unsigned char>(fn[end]))) {
    ++end;
  }
  if (end == start) {
    throw std::runtime_error("Empty parameter name in filename: " + fn);
  }
  return fn.substr(start, end - start);
}

static double extractPrior(const std::string& fn) {
  // Find “…_With_<token>(underscores)(encoded-number)_…”
  // Examples matched:
  //   _With_dcp__1_256_...
  //   _With_m21_7_41e_05_...
  //   _With_sin12_0_303_...
  static const std::regex re(R"(_With_[A-Za-z0-9]+(_+)([0-9eE_]+)_)");
  std::smatch m;
  if (!std::regex_search(fn, m, re) || m.size() < 3) {
    throw std::runtime_error("Bad filename for prior: " + fn);
  }

  std::string delim = m[1];   // one or more underscores just after the token
  std::string core  = m[2];   // the underscore-encoded number (mantissa or mantissa+exp)

  const bool negative = (delim.size() > 1);  // “__” means negative

  // Handle scientific notation like "7_41e_05" (-> 7.41e-05)
  size_t posE = core.find_first_of("eE");
  double val;
  if (posE != std::string::npos) {
    std::string mant = core.substr(0, posE);
    std::string exp  = core.substr(posE + 1);

    // mantissa: underscores become dots   (7_41 -> 7.41)
    std::replace(mant.begin(), mant.end(), '_', '.');

    // exponent: leading '_' means negative sign; keep remaining digits as-is
    if (!exp.empty() && exp[0] == '_') {
      exp[0] = '-';
    }
    // (exp like "-05" or "05" is fine for std::stod)

    val = std::stod(mant + 'e' + exp);
  } else {
    // simple decimal: underscores become dots  (1_256 -> 1.256)
    std::replace(core.begin(), core.end(), '_', '.');
    val = std::stod(core);
  }

  return negative ? -val : val;
}


// -------------------- Main --------------------

int plotOApar() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // === Global axis style (bigger labels/titles everywhere) ===
  const double AXIS_TITLE_SIZE   = 0.055;
  const double AXIS_LABEL_SIZE   = 0.045;
  const double AXIS_TITLE_OFF_X  = 1.10;
  const double AXIS_TITLE_OFF_Y  = 1.2;

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(AXIS_TITLE_SIZE, "XYZ");
  gStyle->SetLabelSize(AXIS_LABEL_SIZE, "XYZ");
  gStyle->SetTitleOffset(AXIS_TITLE_OFF_X, "X");
  gStyle->SetTitleOffset(AXIS_TITLE_OFF_Y, "Y");

  // Give the bigger labels some breathing room
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.08);

  const std::string rootDir = "/gpfs/scratch/uyevarouskay/gundam_results/DUNE/lbl/profiling/output_frozen_TDR_inclND/";
  const std::string pdfOut = "OApar_plots_allFits_inclND.pdf";
  bool opened = false;

  const bool kEnableExtrapolation = false;   // OFF: do not extrapolate
  const bool kShowHighGuidesWhenPossible = true; // keep 5/10 lines only when axis already reaches them


  // param -> vectors of x(prior) and y(chi2 or NLL)
  std::map<std::string, std::vector<double>> chi2map, priormap;

  // 1) scan directory
  for (auto& entry : fs::directory_iterator(rootDir)) {
    if (!entry.is_regular_file()) continue;
    const std::string fn = entry.path().filename().string();
    if (entry.path().extension() != ".root") continue;

    std::string param, fullpath = rootDir + fn;
    double prior = 0.0;

    try {
      param = extractParamName(fn);   // e.g. "sin12"
      prior = extractPrior(fn);       // e.g. 0.00241

      if (param == "MO" && isZero(prior)) {
        std::cerr << "Skip file (param=MO, value=0): " << fn << "\n";
        continue;
      }
    } catch (const std::exception& e) {
      std::cerr << "Skip (name parse): " << e.what() << "\n";
      continue;
    }


    // open root file
    std::unique_ptr<TFile> f(TFile::Open(fullpath.c_str(), "READ"));
    if (!f || f->IsZombie()) { continue; }

    auto* tFit = dynamic_cast<TTree*>(f->Get("FitterEngine/postFit/bestFitStats"));
    if (!tFit) { continue; }

    if (tFit->GetEntries() < 1) { continue; }
    tFit->GetEntry(0);

    auto* leafFit = tFit->GetLeaf("fitStatusCode");
    auto* leafCov = tFit->GetLeaf("covStatusCode");
    auto* leafChi = tFit->GetLeaf("totalLikelihoodAtBestFit"); // may be NLL; Δ works fine

    if (!leafFit || !leafCov || !leafChi) { continue; }

    int fitOK = static_cast<int>(leafFit->GetValue(0));
    int covOK = static_cast<int>(leafCov->GetValue(0));
    //if (fitOK != 0 || covOK != 3) { continue; }

    double chi2 = leafChi->GetValue(0);
    priormap[param].push_back(prior);
    chi2map[param].push_back(chi2);

    std::cout << param << "  prior: " << prior << "  chi2: " << chi2 << "\n";
  }

  // 2) plot per parameter
  for (const auto& kv : priormap) {
    const std::string& param = kv.first;
    auto priors = kv.second;               // mutable copy
    auto vals   = chi2map.at(param);       // mutable copy
    const int n0 = static_cast<int>(priors.size());
    if (n0 == 0) continue;

    // Zip, sort by prior, unzip
    std::vector<std::pair<double,double>> data;
    data.reserve(n0);
    for (int i = 0; i < n0; ++i) data.emplace_back(priors[i], vals[i]);
    std::sort(data.begin(), data.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });
    for (int i = 0; i < n0; ++i) {
      priors[i] = data[i].first;
      vals[i]   = data[i].second;
    }

    // ---------- Build Δχ² ----------
    std::vector<double> dchi2s = vals;
    const double ymin = *std::min_element(dchi2s.begin(), dchi2s.end());
    for (auto& v : dchi2s) v -= ymin;

    // ---------- Parabolic (Gaussian) fit & conditional extrapolation ----------
    // Fit Δχ²(x) ≈ ((x - μ)/σ)^2 over current data span
    TGraph gfit(static_cast<int>(priors.size()), priors.data(), dchi2s.data());
    TF1 quad("quad", "pow((x-[0])/[1], 2)", priors.front(), priors.back());

    // Seed μ at min Δχ², σ ~ span/4
    int iMin = static_cast<int>(std::distance(
                 dchi2s.begin(), std::min_element(dchi2s.begin(), dchi2s.end())));
    double mu0    = priors[iMin];
    double sigma0 = std::max(1e-12, (priors.back() - priors.front()) / 4.0);
    quad.SetParameters(mu0, sigma0);

    // Quiet fit (no draw). If it fails, we just won't extrapolate.
    gfit.Fit(&quad, "Q0"); // Q=quiet, 0=no drawing
    double mu    = quad.GetParameter(0);
    double sigma = std::abs(quad.GetParameter(1));

    const double yTarget = 2.0; // target level for 1-dof ~ 1.414σ
    double yMaxData = *std::max_element(dchi2s.begin(), dchi2s.end());

    if (kEnableExtrapolation && yMaxData < yTarget && priors.size() >= 2 && sigma > 0 && std::isfinite(mu) && std::isfinite(sigma)) {
      // choose side consistent with last x
      double x_last = priors.back();
      double xNew = (x_last >= mu) ? (mu + sigma * std::sqrt(2.0))
                                   : (mu - sigma * std::sqrt(2.0));
      // append (xNew, 2). Keep x monotonic (if needed, re-sort).
      priors.push_back(xNew);
      dchi2s.push_back(yTarget);

      if (xNew < priors[priors.size()-2]) {
        std::vector<std::pair<double,double>> tmp;
        tmp.reserve(priors.size());
        for (size_t i=0;i<priors.size();++i) tmp.emplace_back(priors[i], dchi2s[i]);
        std::sort(tmp.begin(), tmp.end(), [](auto& a, auto& b){ return a.first < b.first; });
        for (size_t i=0;i<tmp.size();++i) { priors[i]=tmp[i].first; dchi2s[i]=tmp[i].second; }
      }
      yMaxData = std::max(yMaxData, yTarget);
    }

    // ---------- Plot (line only) ----------
    const int n = static_cast<int>(priors.size());
    TGraph graph(n, priors.data(), dchi2s.data());
    graph.SetTitle("");
    graph.GetYaxis()->SetTitle("#Delta #chi^{2}");
    graph.GetXaxis()->SetTitle(xTitleFor(param).c_str());

    graph.SetMarkerStyle(1);   // no markers
    graph.SetMarkerSize(0);
    graph.SetLineWidth(2);
    // graph.SetLineColor(kBlue+1); // optional styling

    const double xDataMin = *std::min_element(priors.begin(), priors.end());
    const double xDataMax = *std::max_element(priors.begin(), priors.end());
    const double yDataMin = *std::min_element(dchi2s.begin(), dchi2s.end());
    const double yDataMax = *std::max_element(dchi2s.begin(), dchi2s.end());

    const double yMinAxis = std::min(0.0, yDataMin);
    double yMaxAxis;

    // If the curve never reaches 2, keep the axis to the data (no forcing up to 2/5/10)
    if (yDataMax < 2.0) {
    yMaxAxis = yDataMax * 1.05;  // tiny headroom
    } else {
    // Otherwise, you can keep your previous behavior (e.g., allow 5/10 guides to show)
    yMaxAxis = (kShowHighGuidesWhenPossible ? std::max(10.0, yDataMax) : yDataMax) * 1.05;
    }


    // Decide upper y-axis so dashed guides are visible
    const bool showGuides510 = true;                 // set to false to only require 2
    //const double yMaxDesired = showGuides510 ? 10.0 : 2.0;

    //const double yMinAxis = std::min(0.0, yDataMin);
    //const double yMaxAxis = std::max(yMaxDesired, yDataMax) * 1.05;

    TCanvas canvas(("c_" + param).c_str(), param.c_str(), 900, 650);

    // Special range for δ_CP (if applicable)
    if (param == "dcp" || param == "delta" || param == "deltaCP" || param == "dcp_rad") {
      graph.GetXaxis()->SetLimits(-3.15, 3.15);
    } else {
      graph.GetXaxis()->SetLimits(xDataMin, xDataMax);
    }


    graph.GetYaxis()->SetRangeUser(yMinAxis, yMaxAxis);

    // Draw axes + line
    graph.Draw("AL");

    // Optional: show the fitted parabola (light dashed), comment out if not wanted
    // quad.SetRange(graph.GetXaxis()->GetXmin(), graph.GetXaxis()->GetXmax());
    // quad.SetLineStyle(3);
    // quad.SetLineWidth(1);
    // quad.SetLineColor(kGray+2);
    // quad.Draw("SAME");

    // Horizontal dashed lines at Δχ² = 1, 2, 5, 10
    double xMinAxis = graph.GetXaxis()->GetXmin();
    double xMaxAxis = graph.GetXaxis()->GetXmax();
    /*auto addHLine = [&](double y){
      if (y < yMinAxis || y > yMaxAxis) return;
      TLine* ln = new TLine(xMinAxis, y, xMaxAxis, y);
      ln->SetLineStyle(2);  // dashed
      ln->SetLineWidth(1);
      ln->SetLineColor(kGray+2);
      ln->Draw("SAME");
    };
    addHLine(1.0);
    addHLine(2.0);
    if (showGuides510) { addHLine(5.0); addHLine(10.0); }*/

auto addHLineWithLabel = [&](double y, const char* txt){
    if (y < yMinAxis || y > yMaxAxis) return;
    //double xL = isDCP ? -3.15 : xMinAll;
    //double xR = isDCP ?  3.15 : xMaxAll;
    if (y < yMinAxis || y > yMaxAxis) return;
    TLine* ln = new TLine(xMinAxis, y, xMaxAxis, y);
    // dashed line
    //auto* ln = new TLine(xL, y, xR, y);
    ln->SetLineStyle(2);
    ln->SetLineWidth(1);
    ln->SetLineColor(kGray+2);
    ln->Draw("SAME");

    // label near the right axis, slightly inset
    double dx = 0.01 * (xMaxAxis - xMinAxis);  // 1% inset
    auto* lab = new TLatex(xMaxAxis - dx, y, txt);
    lab->SetTextAlign(31);         // right-justified, vertically centered
    lab->SetTextSize(0.03);        // adjust to taste
    lab->SetTextColor(kGray+2);
    lab->Draw("SAME");
  };


    //addHLine(1.0);
    //addHLine(2.0);
    addHLineWithLabel(1.0, "1#sigma");
    //addHLineWithLabel(4.0, "2#sigma");
    //addHLineWithLabel(9.0, "3#sigma");
    if (yMaxAxis < 1.0) yMaxAxis = 1.05;

    //if (kShowHighGuidesWhenPossible) { addHLine(4.0); addHLine(9.0); }
    if (showGuides510) { addHLineWithLabel(4.0, "2#sigma"); addHLineWithLabel(9.0, "3#sigma"); }

    //canvas.SaveAs((param + ".png").c_str());
    if (!opened) {                     // open the PDF on the very first page
      canvas.Print((pdfOut + "[").c_str());
      opened = true;
    }
    canvas.Print(pdfOut.c_str());      // append this canvas as a page

  }

    // close the PDF after the loop
  if (opened) {
    TCanvas closer;
    closer.Print((pdfOut + "]").c_str());
  }

  return 0;
}

// If running as a ROOT macro: .x plotOApar.cxx+
//void plotOApar_wrapper() { plotOApar(); }
