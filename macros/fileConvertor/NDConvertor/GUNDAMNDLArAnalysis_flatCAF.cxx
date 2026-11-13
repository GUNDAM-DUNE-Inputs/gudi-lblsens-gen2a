// Build:
//   g++ -O2 -std=c++11 GUNDAMNDLArAnalysis_flatCAF.cpp \
//       $(root-config --cflags --libs) -lduneanaobj_StandardRecord \
//       -o GUNDAMNDLArAnalysis_flatCAF
//
// The duneanaobj headers and libraries must be available in the active DUNE
// environment.  SRGlobal is used to read the x coordinates from globalTree.
//
// Run:
//   ./GUNDAMNDLArAnalysis_flatCAF input1.root [input2.root ...]
//
// This program reads a *flat* CAF tree (cafTree) and writes one output entry
// for every selected reconstructed common.ixn.dlp interaction.

#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMatrixT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TVector3.h"
#include "TVectorT.h"

#include "duneanaobj/StandardRecord/SRGlobal.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>


namespace {

// These are deliberately generous for the MiniProdN5 flat CAF files.  The
// program checks the reported per-event sizes before using the buffers.
const int kMaxInteractions = 10000;
const int kMaxParticles = 1000000;
const int kMaxTracks = 1000000;
const int kMaxShowers = 1000000;

// Maximum flattened sizes for one cafTree entry.  These must be checked for
// every production before processing; the MiniProdN5 example is well below
// these limits (about 3.5e5 xsec weights in one entry).
const int kMaxMCInteractions = 10000;
const int kMaxTruthLinks = 1000000;
const int kMaxSystematicDials = 1000000;
const int kMaxSystematicWeights = 2000000;

// Diagnostics are deliberately bounded: metadata is printed once at startup,
// while detailed event/graph information is printed for the first few output
// interactions only.  Set either value to zero to silence that diagnostic.
const bool kPrintSystematicMetadata = false;
const Long64_t kDebugOutputEntries = 0;

constexpr double NDLArXLo = -346.9;
constexpr double NDLArXHi = 346.9;
constexpr double NDLArYLo = -215.5;
constexpr double NDLArYHi = 81.7;
constexpr double NDLArZLo = 418.2;
constexpr double NDLArZHi = 913.3;

enum Sample { kRejected = 0, kNuMuCCLikeContained, kNuECCLikeContained };

struct Particle {
  Int_t pdg;
  Bool_t primary;
  Bool_t contained;
  Float_t sx, sy, sz;
  Float_t ex, ey, ez;
  Float_t px, py, pz;
};

bool InFV(const Float_t x, const Float_t y, const Float_t z)
{
  return x > NDLArXLo + 25. && x < NDLArXHi - 25. &&
         y > NDLArYLo + 25. && y < NDLArYHi - 25. &&
         z > NDLArZLo + 25. && z < NDLArZHi - 25.;
}

double ParticleLength(const Particle& p)
{
  // Exact flat-array counterpart of ana::ParticleLength in RecoTools.h.
  return std::hypot(p.ex - p.sx, p.ey - p.sy, p.ez - p.sz);
}

const Particle* LongestParticle(const std::vector<Particle>& particles)
{
  const Particle* result = 0;
  double longest = -std::numeric_limits<double>::max();
  for (std::vector<Particle>::const_iterator it = particles.begin(); it != particles.end(); ++it) {
    // Exact counterpart of ana::GetLongestParticle: only primary particles
    // participate, and an infinite reconstructed length is not suppressed.
    if (it->primary != 1) continue;
    const double length = ParticleLength(*it);
    if (length > longest) {
      longest = length;
      result = &(*it);
    }
  }
  return result;
}

// Used only for the diagnostic theta-Z branches; it is not part of the
// selection.  The flat CAF stores -inf end points for shower-like DLP
// particles, so angles must use their reconstructed momentum, not end-start.
const Particle* LongestPrimaryParticleWithPID(const std::vector<Particle>& particles,
                                              const Int_t pdg)
{
  const Particle* result = 0;
  double longest = -std::numeric_limits<double>::max();
  for (std::vector<Particle>::const_iterator it = particles.begin(); it != particles.end(); ++it) {
    if (it->primary != 1 || it->pdg != pdg) continue;
    const double length = ParticleLength(*it);
    if (!std::isfinite(length)) continue;
    if (length > longest) {
      longest = length;
      result = &(*it);
    }
  }
  return result;
}

double MomentumThetaZ(const Particle* p)
{
    if (!p)
        return std::numeric_limits<double>::quiet_NaN();

    if (!std::isfinite(p->px) ||
        !std::isfinite(p->py) ||
        !std::isfinite(p->pz)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    TVector3 momentumVector(p->px, p->py, p->pz);

    if (momentumVector.Mag2() <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    momentumVector = momentumVector.Unit();

    const double cosTheta = -momentumVector.Y();

    return cosTheta;
}

bool AllPrimaryParticlesContained(const std::vector<Particle>& particles)
{
  for (std::vector<Particle>::const_iterator it = particles.begin(); it != particles.end(); ++it) {
    if (it->primary == 1 && it->contained != 1) return false;
  }
  return true;
}

bool HasPrimaryPID(const std::vector<Particle>& particles, const Int_t pdg)
{
  for (std::vector<Particle>::const_iterator it = particles.begin(); it != particles.end(); ++it) {
    if (it->primary == 1 && it->pdg == pdg) return true;
  }
  return false;
}

Sample ApplySelectionV0(const Float_t x, const Float_t y, const Float_t z,
                        const std::vector<Particle>& particles)
{
  if (!InFV(x, y, z)) return kRejected;

  const Particle* longest = LongestParticle(particles);
  if (AllPrimaryParticlesContained(particles) && longest && longest->pdg == 13) {
    return kNuMuCCLikeContained;
  }
  if (AllPrimaryParticlesContained(particles) && HasPrimaryPID(particles, 11)) {
    return kNuECCLikeContained;
  }
  return kRejected;
}

Sample ApplySelectionNuELikeStrict(const Float_t x, const Float_t y,
                                   const Float_t z,
                                   const std::vector<Particle>& particles)
{
  if (!InFV(x, y, z)) return kRejected;
  if (!AllPrimaryParticlesContained(particles)) return kRejected;

  // Keep the standard muon-like category unchanged.
  const Particle* longest = LongestParticle(particles);
  if (longest && longest->pdg == 13) {
    return kNuMuCCLikeContained;
  }

  const bool hasPrimaryElectron = HasPrimaryPID(particles, 11);
  const bool hasPrimaryMuon = HasPrimaryPID(particles, 13);

  if (hasPrimaryElectron && !hasPrimaryMuon) {
    return kNuECCLikeContained;
  }

  return kRejected;
}

// Flat-CAF counterpart of duneanasel ApplySelectionV0p1.  In contrast to V0,
// its muon category accepts a reconstructed primary muon that exits downstream
// towards the TMS, provided every other reconstructed particle is contained.
bool IntHasParticles(const std::vector<Particle>& particles)
{
  return !particles.empty();
}

//bool LongestPrimTrackLengthCut(const std::vector<Particle>& particles)
//{
//  const Particle* longest = LongestParticle(particles);
//  return longest && ParticleLength(*longest) > 0.1;
//}

bool LongestPrimTrackLengthCut(const std::vector<Particle>& particles)
{
    for (const Particle& p : particles) {

        if (!p.primary) continue;

        const double length = ParticleLength(p);

        if (!std::isfinite(length)) continue;

        if (length > 0.1) return true;
    }

    return false;
}

bool ParticlesNDLArContainedMuonsEscapeDownstream(
    const std::vector<Particle>& particles)
{
  for (std::vector<Particle>::const_iterator it = particles.begin();
       it != particles.end(); ++it) {
    const Particle& p = *it;
    const bool contained = (p.contained == 1);
    const bool isPrimaryMuon = (p.pdg == 13 && p.primary == 1);
    bool escapesToTMS = false;

    if (isPrimaryMuon) {
      escapesToTMS =
          (p.ex < (NDLArXHi - 25.)) && (p.ex > (NDLArXLo + 25.)) &&
          (p.ey < (NDLArYHi - 25.)) && (p.ey > (NDLArYLo + 25.)) &&
          (p.ez > (NDLArZLo + 25.));
    }
    if (!contained && !escapesToTMS) return false;
  }
  return true;
}

Sample ApplySelectionV0p1(const Float_t x, const Float_t y, const Float_t z,
                          const std::vector<Particle>& particles)
{
  if (!IntHasParticles(particles) ||
      !LongestPrimTrackLengthCut(particles) ||
      !InFV(x, y, z) ||
      !ParticlesNDLArContainedMuonsEscapeDownstream(particles)) {
    return kRejected;
  }

  // This is deliberately HasPrimaryPID rather than the longest-particle PID:
  // it reproduces duneanasel::numode::NuMuCCLike in ApplySelectionV0p1.
  if (HasPrimaryPID(particles, 13)) return kNuMuCCLikeContained;

  // This reproduces NuECCLikeContained.  It intentionally does not apply the
  // additional primary-muon veto used by ApplySelectionNuELikeStrict.
  if (AllPrimaryParticlesContained(particles) && HasPrimaryPID(particles, 11)) {
    return kNuECCLikeContained;
  }
  return kRejected;
}

enum SelectionMode {
  kSelectionV0,
  kSelectionV0p1,
  kSelectionNuELikeStrict
};

// Change this one line to choose the selection.  The strict nue mode is kept
// as the default so that the current output remains unchanged.
const SelectionMode kSelectionMode = kSelectionV0; //kSelectionV0 //kSelectionV0p1 //kSelectionNuELikeStrict

const char* SelectionModeName(const SelectionMode mode)
{
  switch (mode) {
  case kSelectionV0:
    return "V0";
  case kSelectionV0p1:
    return "V0p1";
  case kSelectionNuELikeStrict:
    return "NueStrict";
  }
  return "UnknownSelection";
}

Sample ApplyConfiguredSelection(const Float_t x, const Float_t y,
                                const Float_t z,
                                const std::vector<Particle>& particles)
{
  switch (kSelectionMode) {
  case kSelectionV0:
    return ApplySelectionV0(x, y, z, particles);
  case kSelectionV0p1:
    return ApplySelectionV0p1(x, y, z, particles);
  case kSelectionNuELikeStrict:
    return ApplySelectionNuELikeStrict(x, y, z, particles);
  }
  return kRejected;
}

bool Bind(TChain& chain, const char* branch, void* address)
{
  if (!chain.GetBranch(branch)) {
    std::cerr << "Missing branch: " << branch << std::endl;
    return false;
  }
  if (chain.SetBranchAddress(branch, address) < 0) {
    std::cerr << "Could not bind branch: " << branch << std::endl;
    return false;
  }
  return true;
}

// Each reconstructed interaction can be associated with more than one true
// interaction.  Store systematics for the true interaction with the largest
// overlap, rather than assuming that the reco and MC indices are identical.
Int_t BestMatchedTruth(const size_t iInt,
                       const Int_t* truthLength,
                       const Int_t* truthIdx,
                       const ULong_t* truth,
                       const Float_t* truthOverlap,
                       const Int_t nMC,
                       const Int_t totalTruthLinks,
                       Float_t& bestOverlap)
{
  bestOverlap = -1.f;
  const Int_t first = truthIdx[iInt];
  const Int_t count = truthLength[iInt];
  if (first < 0 || count < 0 || first > totalTruthLinks ||
      count > totalTruthLinks - first) return -1;

  Int_t best = -1;
  for (Int_t j = first; j < first + count; ++j) {
    const ULong_t truthIndex = truth[j];
    if (truthIndex >= static_cast<ULong_t>(nMC)) continue;
    if (truthOverlap[j] > bestOverlap) {
      bestOverlap = truthOverlap[j];
      best = static_cast<Int_t>(truthIndex);
    }
  }
  return best;
}

// Metadata for one systematic dial.  The order in this vector is the CAF's
// outer-vector order.  "id" is retained for diagnostics only: the flattened
// event arrays are indexed by that outer-vector order, not by id.
struct DialMetadata {
  Int_t id;
  std::string name;
  std::vector<double> x;
};

bool NearlyEqual(const double a, const double b)
{
  const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
  return std::abs(a - b) <= 1e-6 * scale;
}

template <class HeaderVector>
void CopyDialMetadata(const HeaderVector& source,
                      std::vector<DialMetadata>& destination)
{
  destination.clear();
  destination.reserve(source.size());
  for (size_t d = 0; d < source.size(); ++d) {
    DialMetadata metadata;
    metadata.id = source[d].id;
    metadata.name = source[d].name;
    metadata.x.reserve(source[d].vals.size());
    for (size_t w = 0; w < source[d].vals.size(); ++w) {
      metadata.x.push_back(static_cast<double>(source[d].vals[w]));
    }
    destination.push_back(metadata);
  }
}

void PrintDialMetadata(const char* family,
                       const std::vector<DialMetadata>& metadata)
{
  std::cout << "[DEBUG] globalTree " << family << " metadata: "
            << metadata.size() << " dials\n";
  for (size_t d = 0; d < metadata.size(); ++d) {
    std::cout << "[DEBUG]   " << family << "[" << d << "] id="
              << metadata[d].id << " name='" << metadata[d].name
              << "' vals:";
    for (size_t w = 0; w < metadata[d].x.size(); ++w) {
      std::cout << " " << metadata[d].x[w];
    }
    std::cout << "\n";
  }
}

// ROOT branch names must not contain punctuation such as '.', ':', '/', or
// whitespace.  Preserve metadata letters, digits, and underscores, replacing
// every other character with '_'.  The full, unmodified parameter names are
// still retained in FluxSysCov_ND.root and XsecSysCov_ND.root.
std::string MakeGraphBranchName(const std::string& parameterName)
{
  std::string result;

  for (size_t i = 0; i < parameterName.size(); ++i) {
    const unsigned char character =
        static_cast<unsigned char>(parameterName[i]);

    result += (std::isalnum(character) || character == '_')
        ? static_cast<char>(character) : '_';
  }

  if (result.empty()) result = "unnamed";
  return result + "_TGraph";
}

bool BuildGraphBranchNames(const std::vector<DialMetadata>& metadata,
                           const char* family,
                           std::vector<std::string>& graphNames,
                           std::vector<std::string>& allGraphNames)
{
  graphNames.clear();
  graphNames.reserve(metadata.size());
  for (size_t d = 0; d < metadata.size(); ++d) {
    const std::string graphName = MakeGraphBranchName(metadata[d].name);
    if (std::find(allGraphNames.begin(), allGraphNames.end(), graphName) !=
        allGraphNames.end()) {
      std::cerr << "Cannot create a unique TTree branch for " << family
                << " dial " << d << " ('" << metadata[d].name
                << "'): sanitized graph name '" << graphName
                << "' conflicts with another dial\n";
      return false;
    }
    graphNames.push_back(graphName);
    allGraphNames.push_back(graphName);
  }
  return true;
}

void PrintGraphBranchNames(const char* family,
                           const std::vector<DialMetadata>& metadata,
                           const std::vector<std::string>& graphNames)
{
  for (size_t d = 0; d < graphNames.size(); ++d) {
    std::cout << "[DEBUG] " << family << " dial " << d << " ('"
              << metadata[d].name << "') -> branch " << graphNames[d]
              << "\n";
  }
}

// Write an identity prior covariance in the object layout expected by the
// supplied FluxSysCov_ND.root and XsecSysCov_ND.root examples.  Parameter names are
// taken from globalTree, so their order is exactly the order used by the
// output flux_* and xsec_* graph branches.
//
// The default is prior=0 and sigma=1.  The explicitly configured xsec
// exceptions below use their required parameter coordinates.  This is still
// not a complete physics covariance model: correlations remain zero and must
// be replaced with the approved flux/xsec covariance when one is available.
void GetCovariancePriorAndSigma(const std::string& family,
                                const std::string& parameterName,
                                double& prior, double& sigma)
{
  prior = 0.0;
  sigma = 1.0;
  if (family != "xsec") return;

  if (parameterName.find("WSReweight_DUNE_v0_nucleus_radius") !=
      std::string::npos) {
    prior = 3.53;
    sigma = 0.5;
    return;
  }
  if (parameterName.find("WSReweight_DUNE_v0_surface_thickness") !=
      std::string::npos) {
    prior = 0.54;
    sigma = 0.05;
    return;
  }

  //if (parameterName.find("CCQETemplateReweight") != std::string::npos ||
  //    parameterName.find("CCQEXSecCorr") != std::string::npos ||
  //    parameterName.find("MECq0q3InterpWeighting") != std::string::npos ||
  //    parameterName.find("QEInterference") != std::string::npos ||
  //    parameterName.find("PionAbsWeighter_DUNE_v0_QuasiDeuteronFraction") !=
  //        std::string::npos) {
  //  prior = 1.0;
  //}
}

bool WriteSystematicCovariance(const char* outputFile,
                               const std::vector<DialMetadata>& metadata,
                               const char* family)
{
  if (metadata.empty()) {
    std::cerr << "Cannot write " << outputFile << ": no " << family
              << " metadata was found in globalTree\n";
    return false;
  }

  const Int_t nParameters = static_cast<Int_t>(metadata.size());
  const std::string prefix = std::string(family) + "_param_";

  TObjArray names(nParameters);
  names.SetName((prefix + "names").c_str());
  names.SetOwner(kTRUE);

  TVectorT<double> priors(nParameters);
  TMatrixT<double> covariance(nParameters, nParameters);
  covariance.Zero();

  for (Int_t d = 0; d < nParameters; ++d) {
    // Do not use a shortened branch name such as xsec_17 here.  The full
    // metadata name is necessary to identify the physical dial unambiguously.
    names.Add(new TObjString(metadata[d].name.c_str()));
    double prior = 0.0;
    double sigma = 1.0;
    GetCovariancePriorAndSigma(family, metadata[d].name, prior, sigma);
    priors[d] = prior;
    covariance(d, d) = sigma * sigma;

    if (std::string(family) == "xsec" &&
        (prior != 0.0 || sigma != 1.0)) {
      std::cout << "[DEBUG] Xsec covariance: '" << metadata[d].name
                << "' prior=" << prior << " sigma=" << sigma << "\n";
    }
  }

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Could not create covariance file " << outputFile << "\n";
    return false;
  }

  names.Write((prefix + "names").c_str(), TObject::kSingleKey);
  output.WriteObject(&priors, (prefix + "priors").c_str());
  output.WriteObject(&covariance, (prefix + "cov").c_str());
  output.Close();

  std::cout << "Wrote " << outputFile << " with:\n"
            << "  - " << nParameters << " parameter names from globalTree\n"
            << "  - Default prior: 0.0; default sigma: 1.0\n";
  if (std::string(family) == "xsec") {
    std::cout << "  - Applied configured xsec prior/sigma overrides (see DEBUG)\n";
  }
  std::cout << "  - Covariance matrix: " << nParameters << " x "
            << nParameters << " (diagonal)" << std::endl;
  return true;
}

bool SameDialMetadata(const std::vector<DialMetadata>& reference,
                      const std::vector<DialMetadata>& candidate,
                      const char* family, const char* filename)
{
  if (reference.size() != candidate.size()) {
    std::cerr << filename << ": globalTree " << family
              << " dial count differs from the first input file\n";
    return false;
  }
  for (size_t d = 0; d < reference.size(); ++d) {
    if (reference[d].id != candidate[d].id ||
        reference[d].name != candidate[d].name ||
        reference[d].x.size() != candidate[d].x.size()) {
      std::cerr << filename << ": globalTree " << family
                << " metadata differs at dial " << d << "\n";
      return false;
    }
    for (size_t w = 0; w < reference[d].x.size(); ++w) {
      if (!NearlyEqual(reference[d].x[w], candidate[d].x[w])) {
        std::cerr << filename << ": globalTree " << family
                  << " x coordinate differs at dial " << d
                  << ", knot " << w << "\n";
        return false;
      }
    }
  }
  return true;
}

bool LoadGlobalMetadata(const char* filename,
                        std::vector<DialMetadata>& xsecMetadata,
                        std::vector<DialMetadata>& fluxMetadata)
{
  TFile input(filename, "READ");
  if (input.IsZombie()) {
    std::cerr << "Could not open input file " << filename << "\n";
    return false;
  }

  TTree* globalTree = dynamic_cast<TTree*>(input.Get("globalTree"));
  if (!globalTree || globalTree->GetEntries() < 1) {
    std::cerr << filename << ": no populated globalTree; cannot build "
              << "metadata-based systematic graphs\n";
    return false;
  }

  caf::SRGlobal* global = 0;
  if (globalTree->SetBranchAddress("global", &global) < 0 ||
      globalTree->GetEntry(0) <= 0 || !global) {
    std::cerr << filename << ": could not read globalTree.global\n";
    return false;
  }

  CopyDialMetadata(global->wgts.xsec_params, xsecMetadata);
  CopyDialMetadata(global->wgts.flux_params, fluxMetadata);
  globalTree->ResetBranchAddresses();
  return true;
}

// Each TGraph is the raw NuSystematics response for exactly one dial:
//
//   x = global.wgts.<family>_params[d].vals[w]
//   y = rec.mc.nu.<family>_systs[d].weights[w]
//
// xsec_cvwgt is deliberately not used here.  It is copied unchanged to the
// output xsec_cvwgt branch for provenance, but this converter does not apply
// it to a nominal event weight or to any systematic-response graph.

bool FillSeparateDialBranches(
    std::vector<TClonesArray*>& outputs,
    const std::vector<DialMetadata>& metadata,
    const std::vector<std::string>& graphNames,
    const Int_t nDials,
    const Int_t firstDial,
    const Int_t totalDials,
    const Int_t* weightLength,
    const Int_t* weightIdx,
    const Int_t totalWeights,
    const Float_t* weights,
    const char* family,
    const bool printDebug)
{
  for (size_t d = 0; d < outputs.size(); ++d) outputs[d]->Clear("C");

  if (nDials != static_cast<Int_t>(outputs.size()) ||
      metadata.size() != outputs.size() ||
      graphNames.size() != outputs.size()) {
    std::cerr << family << ": inconsistent number of dials\n";
    return false;
  }

  if (firstDial < 0 || nDials < 0 ||
      firstDial > totalDials || nDials > totalDials - firstDial) {
    return false;
  }

  for (Int_t d = 0; d < nDials; ++d) {
    const Int_t dial = firstDial + d;
    const Int_t nWeights = weightLength[dial];
    const Int_t firstWeight = weightIdx[dial];

    if (nWeights < 1 || firstWeight < 0 ||
        firstWeight > totalWeights ||
        nWeights > totalWeights - firstWeight ||
        nWeights != static_cast<Int_t>(metadata[d].x.size())) {
      std::cerr << family << "_" << d << " (id=" << metadata[d].id
                << ", name='" << metadata[d].name
                << "'): event weight count does not match globalTree\n";
      return false;
    }

    std::vector<std::pair<double, double> > points;
    points.reserve(nWeights);

    for (Int_t w = 0; w < nWeights; ++w) {
      const double x = metadata[d].x[w];
      double y = weights[firstWeight + w];

      // A non-finite x value means invalid globalTree metadata: stop.
      if (!std::isfinite(x)) {
        std::cerr << family << "_" << d
                  << ": non-finite x coordinate\n";
        return false;
      }

      // For a one-knot dial, use the neutral response.
      if (!std::isfinite(y) && nWeights == 1) {
        std::cerr << "[WARNING] " << family << "_" << d
                  << " (name='" << metadata[d].name
                  << "'): one-knot non-finite weight replaced by 1.0\n";
        y = 1.0;
      }

      // For a multi-knot dial, omit only the bad knot.
      if (!std::isfinite(y)) {
        std::cerr << "[WARNING] " << family << "_" << d
                  << " (name='" << metadata[d].name
                  << "', variation=" << w
                  << "): non-finite knot omitted from TGraph\n";
        continue;
      }

      points.push_back(std::make_pair(x, y));
    }

    // Preserve only CAF-provided knots.  In particular, do not infer a
    // central x coordinate from the event weights and do not add a synthetic
    // knot for physical-coordinate grids such as the Woods-Saxon dials.
    bool hasZeroKnot = false;
    for (size_t p = 0; p < points.size(); ++p) {
      if (std::abs(points[p].first) < 1e-8) {
        hasZeroKnot = true;
        break;
      }
    }

    if (printDebug && !hasZeroKnot) {
      std::cout << "[DEBUG] " << family << "[" << d << "] id="
                << metadata[d].id << " has no x=0 knot; preserving its "
                << "globalTree coordinates without adding a knot\n";
    }

    std::sort(points.begin(), points.end(),
              [](const std::pair<double, double>& a,
                 const std::pair<double, double>& b) {
                return a.first < b.first;
              });

    TGraph* graph = new ((*outputs[d])[0]) TGraph();
    for (size_t p = 0; p < points.size(); ++p) {
      graph->SetPoint(p, points[p].first, points[p].second);
    }

    graph->SetName(graphNames[d].c_str());

    if (printDebug) {
      std::cout << "[DEBUG] " << family << " graph " << graphNames[d]
                << " (metadata id=" << metadata[d].id << ", flattened "
                << "dial=" << dial << ", firstWeight=" << firstWeight
                << ") knots:";
      for (size_t p = 0; p < points.size(); ++p) {
        std::cout << " (" << points[p].first << "," << points[p].second
                  << ")";
      }
      std::cout << "\n";
    }
  }

  return true;
}

} // namespace

int main(int argc, char const* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " input1.flat.root [input2.flat.root ...]" << std::endl;
    return 1;
  }

  // The flat cafTree holds the event weights, while globalTree holds the
  // corresponding dial names and x coordinates.  Load metadata from every
  // file before chaining, and refuse to combine files with different axes.
  std::vector<DialMetadata> xsecMetadata;
  std::vector<DialMetadata> fluxMetadata;
  if (!LoadGlobalMetadata(argv[1], xsecMetadata, fluxMetadata)) return 2;

  for (int inputIndex = 2; inputIndex < argc; ++inputIndex) {
    std::vector<DialMetadata> fileXsecMetadata;
    std::vector<DialMetadata> fileFluxMetadata;
    if (!LoadGlobalMetadata(argv[inputIndex], fileXsecMetadata,
                            fileFluxMetadata) ||
        !SameDialMetadata(xsecMetadata, fileXsecMetadata, "xsec",
                          argv[inputIndex]) ||
        !SameDialMetadata(fluxMetadata, fileFluxMetadata, "flux",
                          argv[inputIndex])) {
      return 2;
    }
  }

  // The branch and the TGraph stored in it share the same metadata-derived
  // name: TGraph_<parameter-name>.  Check both families together because all
  // branches live in the same output TTree namespace.
  std::vector<std::string> allGraphNames;
  std::vector<std::string> fluxGraphNames;
  std::vector<std::string> xsecGraphNames;
  if (!BuildGraphBranchNames(fluxMetadata, "flux", fluxGraphNames,
                             allGraphNames) ||
      !BuildGraphBranchNames(xsecMetadata, "xsec", xsecGraphNames,
                             allGraphNames)) {
    return 2;
  }

  if (kPrintSystematicMetadata) {
    PrintDialMetadata("xsec", xsecMetadata);
    PrintDialMetadata("flux", fluxMetadata);
    PrintGraphBranchNames("flux", fluxMetadata, fluxGraphNames);
    PrintGraphBranchNames("xsec", xsecMetadata, xsecGraphNames);
  }

  // These files are written once per converter run.  Their parameter ordering
  // and full names are guaranteed to match the graph branches below.
  if (!WriteSystematicCovariance("FluxSysCov_ND.root", fluxMetadata, "flux") ||
      !WriteSystematicCovariance("XsecSysCov_ND.root", xsecMetadata, "xsec")) {
    return 2;
  }

  TChain ch("cafTree");
  for (int i = 1; i < argc; ++i) ch.Add(argv[i]);

  // These buffers directly replace SR->common.ixn.dlp and
  // SR->common.ixn.dlp[i].part.dlp in the standard CAF program.
  static Int_t nCommon = 0;
  static Int_t nLar = 0;
  static Float_t vtxX[kMaxInteractions], vtxY[kMaxInteractions], vtxZ[kMaxInteractions];
  static Int_t partLength[kMaxInteractions], partIdx[kMaxInteractions];
  static Int_t partPDG[kMaxParticles];
  static Bool_t partPrimary[kMaxParticles], partContained[kMaxParticles];
  static Float_t partStartX[kMaxParticles], partStartY[kMaxParticles], partStartZ[kMaxParticles];
  static Float_t partEndX[kMaxParticles], partEndY[kMaxParticles], partEndZ[kMaxParticles];
  static Float_t partMomX[kMaxParticles], partMomY[kMaxParticles], partMomZ[kMaxParticles];

  // ND-LAr track and shower summaries, indexed by their parent reconstructed
  // interaction through the flat-array length and idx branches.
  static Int_t trackTotal = 0, showerTotal = 0;
  static Int_t trackLength[kMaxInteractions], trackIdx[kMaxInteractions];
  static Float_t trackE[kMaxTracks], trackEvis[kMaxTracks], trackLen[kMaxTracks];
  static Int_t showerLength[kMaxInteractions], showerIdx[kMaxInteractions];
  static Float_t showerEvis[kMaxShowers];

  // Common-interaction energy and CVN quantities.  In the current SPINE-only
  // MiniProdN5 files these are generally NaN, but they are copied if present.
  static Float_t enuCalo[kMaxInteractions], enuECalo[kMaxInteractions], enuLepCalo[kMaxInteractions];
  static Float_t cvnNue[kMaxInteractions], cvnNumu[kMaxInteractions];

  // MC truth and reconstructed-interaction to truth associations.
  static Int_t nMC = 0;
  static Int_t truthTotal = 0;
  static Int_t recoTruthLength[kMaxInteractions], recoTruthIdx[kMaxInteractions];
  static ULong_t recoTruth[kMaxTruthLinks];
  static Float_t recoTruthOverlap[kMaxTruthLinks];
  static Float_t mcXSec[kMaxMCInteractions], mcXSecCVWgt[kMaxMCInteractions];
  static Float_t mcE[kMaxMCInteractions], mcGenWeight[kMaxMCInteractions];
  static Float_t mcMomX[kMaxMCInteractions], mcMomY[kMaxMCInteractions], mcMomZ[kMaxMCInteractions];
  static Int_t mcPDG[kMaxMCInteractions], mcPDGOrig[kMaxMCInteractions], mcMode[kMaxMCInteractions];
  static Bool_t mcIsCC[kMaxMCInteractions];

  // Nested flat-CAF systematic arrays.  The top level is indexed by MC
  // interaction, then by dial, then by universe/variation weight.
  static Int_t fluxDialTotal = 0, fluxWeightTotal = 0;
  static Int_t fluxDialLength[kMaxMCInteractions], fluxDialIdx[kMaxMCInteractions];
  static Int_t fluxWeightLength[kMaxSystematicDials], fluxWeightIdx[kMaxSystematicDials];
  static Float_t fluxWeights[kMaxSystematicWeights];

  static Int_t xsecDialTotal = 0, xsecWeightTotal = 0;
  static Int_t xsecDialLength[kMaxMCInteractions], xsecDialIdx[kMaxMCInteractions];
  static Int_t xsecWeightLength[kMaxSystematicDials], xsecWeightIdx[kMaxSystematicDials];
  static Float_t xsecWeights[kMaxSystematicWeights];

  if (!Bind(ch, "rec.common.ixn.dlp..length", &nCommon) ||
      !Bind(ch, "rec.nd.lar.dlp..length", &nLar) ||
      !Bind(ch, "rec.common.ixn.dlp.vtx.x", vtxX) ||
      !Bind(ch, "rec.common.ixn.dlp.vtx.y", vtxY) ||
      !Bind(ch, "rec.common.ixn.dlp.vtx.z", vtxZ) ||
      !Bind(ch, "rec.common.ixn.dlp.Enu.calo", enuCalo) ||
      !Bind(ch, "rec.common.ixn.dlp.Enu.e_calo", enuECalo) ||
      !Bind(ch, "rec.common.ixn.dlp.Enu.lep_calo", enuLepCalo) ||
      !Bind(ch, "rec.common.ixn.dlp.nuhyp.cvn.nue", cvnNue) ||
      !Bind(ch, "rec.common.ixn.dlp.nuhyp.cvn.numu", cvnNumu) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp..length", partLength) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp..idx", partIdx) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.pdg", partPDG) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.primary", partPrimary) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.contained", partContained) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.start.x", partStartX) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.start.y", partStartY) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.start.z", partStartZ) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.end.x", partEndX) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.end.y", partEndY) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.end.z", partEndZ) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.p.x", partMomX) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.p.y", partMomY) ||
      !Bind(ch, "rec.common.ixn.dlp.part.dlp.p.z", partMomZ) ||

      !Bind(ch, "rec.nd.lar.dlp.tracks..totarraysize", &trackTotal) ||
      !Bind(ch, "rec.nd.lar.dlp.tracks..length", trackLength) ||
      !Bind(ch, "rec.nd.lar.dlp.tracks..idx", trackIdx) ||
      !Bind(ch, "rec.nd.lar.dlp.tracks.E", trackE) ||
      !Bind(ch, "rec.nd.lar.dlp.tracks.Evis", trackEvis) ||
      !Bind(ch, "rec.nd.lar.dlp.tracks.len_cm", trackLen) ||
      !Bind(ch, "rec.nd.lar.dlp.showers..totarraysize", &showerTotal) ||
      !Bind(ch, "rec.nd.lar.dlp.showers..length", showerLength) ||
      !Bind(ch, "rec.nd.lar.dlp.showers..idx", showerIdx) ||
      !Bind(ch, "rec.nd.lar.dlp.showers.Evis", showerEvis) ||

      !Bind(ch, "rec.mc.nu..length", &nMC) ||
      !Bind(ch, "rec.common.ixn.dlp.truth..totarraysize", &truthTotal) ||
      !Bind(ch, "rec.common.ixn.dlp.truth..length", recoTruthLength) ||
      !Bind(ch, "rec.common.ixn.dlp.truth..idx", recoTruthIdx) ||
      !Bind(ch, "rec.common.ixn.dlp.truth", recoTruth) ||
      !Bind(ch, "rec.common.ixn.dlp.truthOverlap", recoTruthOverlap) ||
      !Bind(ch, "rec.mc.nu.xsec", mcXSec) ||
      !Bind(ch, "rec.mc.nu.xsec_cvwgt", mcXSecCVWgt) ||
      !Bind(ch, "rec.mc.nu.E", mcE) ||
      !Bind(ch, "rec.mc.nu.genweight", mcGenWeight) ||
      !Bind(ch, "rec.mc.nu.momentum.x", mcMomX) ||
      !Bind(ch, "rec.mc.nu.momentum.y", mcMomY) ||
      !Bind(ch, "rec.mc.nu.momentum.z", mcMomZ) ||
      !Bind(ch, "rec.mc.nu.pdg", mcPDG) ||
      !Bind(ch, "rec.mc.nu.pdgorig", mcPDGOrig) ||
      !Bind(ch, "rec.mc.nu.mode", mcMode) ||
      !Bind(ch, "rec.mc.nu.iscc", mcIsCC) ||

      !Bind(ch, "rec.mc.nu.flux_systs..totarraysize", &fluxDialTotal) ||
      !Bind(ch, "rec.mc.nu.flux_systs..length", fluxDialLength) ||
      !Bind(ch, "rec.mc.nu.flux_systs..idx", fluxDialIdx) ||
      !Bind(ch, "rec.mc.nu.flux_systs.weights..totarraysize", &fluxWeightTotal) ||
      !Bind(ch, "rec.mc.nu.flux_systs.weights..length", fluxWeightLength) ||
      !Bind(ch, "rec.mc.nu.flux_systs.weights..idx", fluxWeightIdx) ||
      !Bind(ch, "rec.mc.nu.flux_systs.weights", fluxWeights) ||

      !Bind(ch, "rec.mc.nu.xsec_systs..totarraysize", &xsecDialTotal) ||
      !Bind(ch, "rec.mc.nu.xsec_systs..length", xsecDialLength) ||
      !Bind(ch, "rec.mc.nu.xsec_systs..idx", xsecDialIdx) ||
      !Bind(ch, "rec.mc.nu.xsec_systs.weights..totarraysize", &xsecWeightTotal) ||
      !Bind(ch, "rec.mc.nu.xsec_systs.weights..length", xsecWeightLength) ||
      !Bind(ch, "rec.mc.nu.xsec_systs.weights..idx", xsecWeightIdx) ||
      !Bind(ch, "rec.mc.nu.xsec_systs.weights", xsecWeights)) {
    return 2;
  }

  const std::string outputName = std::string("ndlar_selected_flat_") +
                                 SelectionModeName(kSelectionMode) + "_GlobalAxes.root";
  TFile output(outputName.c_str(), "RECREATE");
  TTree outtree("event_tree", "Selected ND-LAr DLP interactions from flat CAF");

  Int_t event = -1, interaction = -1, sample = kRejected;
  Int_t passFV = 0, passContained = 0, hasLongestMuon = 0, hasElectronPID = 0;
  Int_t longestPDG = 0;
  Float_t longestLen = -999.f;
  Bool_t nue_CC = kFALSE, numu_CC = kFALSE, anue_CC = kFALSE, anumu_CC = kFALSE;
  Float_t recoEnu = -999.f, recoEnuVis = -999.f, recoLepE = -999.f, yRec = -999.f;
  Int_t hasLarMatch = 0, nTracks = 0, nShowers = 0;
  Float_t longestTrackLength = -999.f;
  Float_t sumTrackE = 0.f, sumTrackEvis = 0.f, sumShowerEvis = 0.f;
  Float_t recoEnuFromTracksOnly = -999.f, recoEnuFromTracksPlusShowers = -999.f;
  Float_t recoEnuCalo = -999.f, recoLepEFromTracks = -999.f, recoLepEFromShowers = -999.f;
  Float_t caloVsTracksRelDiff = -999.f, caloVsTracksShowersRelDiff = -999.f;
  Int_t usedRecoEnuFallback = 0, usedRecoLepEFallback = 0;
  Double_t vtx_x = -999., vtx_y = -999., vtx_z = -999.;

  Int_t matchedMC = -1;
  Float_t matchedTruthOverlap = -1.f;
  // xsec_cvwgt is copied to the output tree for each matched event.  It is
  // informational only in this converter: it is not applied to a nominal
  // event weight and is not folded into any dial-response TGraph.
  Float_t xsec = -999.f, xsec_cvwgt = 1.f;
  Double_t Ev = -999., mcEnu = -999.;
  Int_t mcModeOut = -999, isCC = -1, isNC = -1, nuPDG = -1, nuPDGunosc = -1;
  Int_t file_idx = -1, iniNuFlavor = -1, interNuFlavor = -1, detNuFlavor = -1;
  Int_t ATM = 0, LBL = 1, isFD = 0, isRHC = 0;
  Float_t mcGenWeightOut = -999.f, mcthetaZ = -999.f;
  Double_t POTScaledOscweight = 1.312500e+17;
  Double_t recoEnu_numu = -999., recoEnu_nue = -999.;
  Double_t cvnnue = std::numeric_limits<double>::quiet_NaN();
  Double_t cvnnumu = std::numeric_limits<double>::quiet_NaN();
  Double_t recoEle_thetaZ = std::numeric_limits<double>::quiet_NaN();
  Double_t recoMu_thetaZ = std::numeric_limits<double>::quiet_NaN();

  // GUNDAM-style branch layout: one TClonesArray branch per dial.  Both the
  // number of branches and every graph x coordinate come from globalTree.
  std::vector<TClonesArray*> flux_dials(fluxMetadata.size(), nullptr);
  std::vector<TClonesArray*> xsec_dials(xsecMetadata.size(), nullptr);

  outtree.Branch("event", &event, "event/I");
  outtree.Branch("interaction", &interaction, "interaction/I");
  outtree.Branch("sample", &sample, "sample/I");
  outtree.Branch("passFV", &passFV, "passFV/I");
  outtree.Branch("passContained", &passContained, "passContained/I");
  outtree.Branch("hasLongestMuon", &hasLongestMuon, "hasLongestMuon/I");
  outtree.Branch("hasElectronPID", &hasElectronPID, "hasElectronPID/I");
  outtree.Branch("longestPDG", &longestPDG, "longestPDG/I");
  outtree.Branch("longestLen", &longestLen, "longestLen/F");
  outtree.Branch("nue_CC", &nue_CC, "nue_CC/O");
  outtree.Branch("numu_CC", &numu_CC, "numu_CC/O");
  outtree.Branch("anue_CC", &anue_CC, "anue_CC/O");
  outtree.Branch("anumu_CC", &anumu_CC, "anumu_CC/O");
  outtree.Branch("recoEnu", &recoEnu, "recoEnu/F");
  outtree.Branch("recoEnuVis", &recoEnuVis, "recoEnuVis/F");
  outtree.Branch("recoLepE", &recoLepE, "recoLepE/F");
  outtree.Branch("yRec", &yRec, "yRec/F");
  outtree.Branch("hasLarMatch", &hasLarMatch, "hasLarMatch/I");
  outtree.Branch("nTracks", &nTracks, "nTracks/I");
  outtree.Branch("nShowers", &nShowers, "nShowers/I");
  outtree.Branch("longestTrackLength", &longestTrackLength, "longestTrackLength/F");
  outtree.Branch("sumTrackE", &sumTrackE, "sumTrackE/F");
  outtree.Branch("sumTrackEvis", &sumTrackEvis, "sumTrackEvis/F");
  outtree.Branch("sumShowerEvis", &sumShowerEvis, "sumShowerEvis/F");
  outtree.Branch("recoEnuFromTracksOnly", &recoEnuFromTracksOnly, "recoEnuFromTracksOnly/F");
  outtree.Branch("recoEnuFromTracksPlusShowers", &recoEnuFromTracksPlusShowers, "recoEnuFromTracksPlusShowers/F");
  outtree.Branch("recoEnuCalo", &recoEnuCalo, "recoEnuCalo/F");
  outtree.Branch("recoLepEFromTracks", &recoLepEFromTracks, "recoLepEFromTracks/F");
  outtree.Branch("recoLepEFromShowers", &recoLepEFromShowers, "recoLepEFromShowers/F");
  outtree.Branch("caloVsTracksRelDiff", &caloVsTracksRelDiff, "caloVsTracksRelDiff/F");
  outtree.Branch("caloVsTracksShowersRelDiff", &caloVsTracksShowersRelDiff, "caloVsTracksShowersRelDiff/F");
  outtree.Branch("usedRecoEnuFallback", &usedRecoEnuFallback, "usedRecoEnuFallback/I");
  outtree.Branch("usedRecoLepEFallback", &usedRecoLepEFallback, "usedRecoLepEFallback/I");
  outtree.Branch("vtx_x", &vtx_x, "vtx_x/D");
  outtree.Branch("vtx_y", &vtx_y, "vtx_y/D");
  outtree.Branch("vtx_z", &vtx_z, "vtx_z/D");
  outtree.Branch("matched_mc", &matchedMC, "matched_mc/I");
  outtree.Branch("truth_overlap", &matchedTruthOverlap, "truth_overlap/F");
  outtree.Branch("xsec", &xsec, "xsec/F");
  outtree.Branch("xsec_cvwgt", &xsec_cvwgt, "xsec_cvwgt/F");
  outtree.Branch("Ev", &Ev, "Ev/D");
  outtree.Branch("mcEnu", &mcEnu, "mcEnu/D");
  outtree.Branch("mcMode", &mcModeOut, "mcMode/I");
  outtree.Branch("isCC", &isCC, "isCC/I");
  outtree.Branch("isNC", &isNC, "isNC/I");
  outtree.Branch("nuPDG", &nuPDG, "nuPDG/I");
  outtree.Branch("nuPDGunosc", &nuPDGunosc, "nuPDGunosc/I");
  outtree.Branch("file_idx", &file_idx, "file_idx/I");
  outtree.Branch("iniNuFlavor", &iniNuFlavor, "iniNuFlavor/I");
  outtree.Branch("interNuFlavor", &interNuFlavor, "interNuFlavor/I");
  outtree.Branch("detNuFlavor", &detNuFlavor, "detNuFlavor/I");
  outtree.Branch("ATM", &ATM, "ATM/I");
  outtree.Branch("LBL", &LBL, "LBL/I");
  outtree.Branch("isFD", &isFD, "isFD/I");
  outtree.Branch("isRHC", &isRHC, "isRHC/I");
  outtree.Branch("mcGenWeight", &mcGenWeightOut, "mcGenWeight/F");
  outtree.Branch("POTScaledOscweight", &POTScaledOscweight, "POTScaledOscweight/D");
  outtree.Branch("mcthetaZ", &mcthetaZ, "mcthetaZ/F");
  outtree.Branch("recoEnu_numu", &recoEnu_numu, "recoEnu_numu/D");
  outtree.Branch("recoEnu_nue", &recoEnu_nue, "recoEnu_nue/D");
  outtree.Branch("cvnnue", &cvnnue, "cvnnue/D");
  outtree.Branch("cvnnumu", &cvnnumu, "cvnnumu/D");
  outtree.Branch("recoEle_thetaZ", &recoEle_thetaZ, "recoEle_thetaZ/D");
  outtree.Branch("recoMu_thetaZ", &recoMu_thetaZ, "recoMu_thetaZ/D");

  for (size_t dial = 0; dial < flux_dials.size(); ++dial) {
    flux_dials[dial] = new TClonesArray("TGraph", 1);
    outtree.Branch(fluxGraphNames[dial].c_str(), &flux_dials[dial], 32000, 0);
  }
  for (size_t dial = 0; dial < xsec_dials.size(); ++dial) {
    xsec_dials[dial] = new TClonesArray("TGraph", 1);
    outtree.Branch(xsecGraphNames[dial].c_str(), &xsec_dials[dial], 32000, 0);
  }

  const Long64_t entries = ch.GetEntries();
  std::cout << "Input tree has " << entries << " entries." << std::endl;
  std::cout << "[DEBUG] Output will contain " << flux_dials.size()
            << " flux graphs and " << xsec_dials.size()
            << " xsec graphs per matched interaction." << std::endl;

  Long64_t nMu = 0, nE = 0;

  for (Long64_t i = 0; i < entries; ++i) {
    ch.GetEntry(i);
    //std::cout << "Event " << i << " nCommon=" << nCommon << " nLar=" << nLar << std::endl;

    if (nCommon < 0 || nLar < 0 || nMC < 0 ||
        nCommon > kMaxInteractions || nLar > kMaxInteractions ||
        nMC > kMaxMCInteractions || truthTotal < 0 ||
        truthTotal > kMaxTruthLinks || fluxDialTotal < 0 ||
        fluxDialTotal > kMaxSystematicDials || fluxWeightTotal < 0 ||
        fluxWeightTotal > kMaxSystematicWeights || xsecDialTotal < 0 ||
        xsecDialTotal > kMaxSystematicDials || xsecWeightTotal < 0 ||
        xsecWeightTotal > kMaxSystematicWeights || trackTotal < 0 ||
        trackTotal > kMaxTracks || showerTotal < 0 || showerTotal > kMaxShowers) {
      std::cerr << "Event " << i << " exceeds a configured input buffer." << std::endl;
      return 3;
    }
    const size_t nReco = std::min(static_cast<size_t>(nCommon), static_cast<size_t>(nLar));

    for (size_t iInt = 0; iInt < nReco; ++iInt) {
      const Int_t first = partIdx[iInt];
      const Int_t count = partLength[iInt];
      if (first < 0 || count < 0 || first > kMaxParticles || count > kMaxParticles - first) {
        std::cerr << "Event " << i << ", interaction " << iInt
                  << " exceeds particle buffer." << std::endl;
        return 3;
      }

      std::vector<Particle> particles;
      particles.reserve(count);
      for (Int_t j = first; j < first + count; ++j) {
        Particle p;
        p.pdg = partPDG[j]; p.primary = partPrimary[j]; p.contained = partContained[j];
        p.sx = partStartX[j]; p.sy = partStartY[j]; p.sz = partStartZ[j];
        p.ex = partEndX[j]; p.ey = partEndY[j]; p.ez = partEndZ[j];
        p.px = partMomX[j]; p.py = partMomY[j]; p.pz = partMomZ[j];
        particles.push_back(p);
      }

      event = static_cast<Int_t>(i);
      interaction = static_cast<Int_t>(iInt);
      file_idx = ch.GetTreeNumber();
      vtx_x = vtxX[iInt]; vtx_y = vtxY[iInt]; vtx_z = vtxZ[iInt];
      recoEnu = recoEnuVis = recoLepE = yRec = -999.f;
      hasLarMatch = 1;
      nTracks = nShowers = 0;
      longestTrackLength = -999.f;
      sumTrackE = sumTrackEvis = sumShowerEvis = 0.f;
      recoEnuFromTracksOnly = recoEnuFromTracksPlusShowers = recoEnuCalo = -999.f;
      recoLepEFromTracks = recoLepEFromShowers = -999.f;
      caloVsTracksRelDiff = caloVsTracksShowersRelDiff = -999.f;
      usedRecoEnuFallback = usedRecoLepEFallback = 0;
      recoEnu_numu = recoEnu_nue = -999.;
      cvnnue = static_cast<Double_t>(cvnNue[iInt]);
      cvnnumu = static_cast<Double_t>(cvnNumu[iInt]);
      recoEle_thetaZ = std::numeric_limits<double>::quiet_NaN();
      recoMu_thetaZ = std::numeric_limits<double>::quiet_NaN();
      matchedMC = -1;
      matchedTruthOverlap = -1.f;
      Ev = mcEnu = -999.;
      mcModeOut = -999;
      isCC = isNC = -1;
      nuPDG = nuPDGunosc = iniNuFlavor = interNuFlavor = detNuFlavor = -1;
      mcGenWeightOut = -999.f;
      mcthetaZ = -999.f;
      anue_CC = anumu_CC = kFALSE;
      passFV = InFV(vtxX[iInt], vtxY[iInt], vtxZ[iInt]) ? 1 : 0;
      passContained = AllPrimaryParticlesContained(particles) ? 1 : 0;
      hasElectronPID = HasPrimaryPID(particles, 11) ? 1 : 0;

      const Particle* longest = LongestParticle(particles);
      longestPDG = longest ? longest->pdg : 0;
      longestLen = longest ? static_cast<Float_t>(ParticleLength(*longest)) : -999.f;
      hasLongestMuon = (longest && longest->pdg == 13) ? 1 : 0;
      recoEle_thetaZ = MomentumThetaZ(LongestPrimaryParticleWithPID(particles, 11));
      recoMu_thetaZ = MomentumThetaZ(LongestPrimaryParticleWithPID(particles, 13));

      sample = ApplyConfiguredSelection(vtxX[iInt], vtxY[iInt],
                                        vtxZ[iInt], particles);
      numu_CC = (sample == kNuMuCCLikeContained);
      nue_CC = (sample == kNuECCLikeContained);
      detNuFlavor = numu_CC ? 14 : (nue_CC ? 12 : 0);

      // A common interaction is paired to the same-index ND-LAr interaction
      // because nReco is the minimum of the two collection sizes.
      const Int_t trackFirst = trackIdx[iInt];
      const Int_t trackCount = trackLength[iInt];
      const Int_t showerFirst = showerIdx[iInt];
      const Int_t showerCount = showerLength[iInt];
      if (trackFirst < 0 || trackCount < 0 || trackFirst > trackTotal ||
          trackCount > trackTotal - trackFirst || showerFirst < 0 ||
          showerCount < 0 || showerFirst > showerTotal ||
          showerCount > showerTotal - showerFirst) {
        std::cerr << "Event " << i << ", interaction " << iInt
                  << " has an invalid track or shower range." << std::endl;
        return 3;
      }

      nTracks = trackCount;
      nShowers = showerCount;
      for (Int_t j = trackFirst; j < trackFirst + trackCount; ++j) {
        if (std::isfinite(trackE[j])) sumTrackE += trackE[j];
        if (std::isfinite(trackEvis[j])) sumTrackEvis += trackEvis[j];
        if (std::isfinite(trackLen[j]) && trackLen[j] > longestTrackLength) {
          longestTrackLength = trackLen[j];
        }
      }
      for (Int_t j = showerFirst; j < showerFirst + showerCount; ++j) {
        if (std::isfinite(showerEvis[j])) sumShowerEvis += showerEvis[j];
      }

      recoEnuFromTracksOnly = sumTrackE;
      recoEnuFromTracksPlusShowers = sumTrackE + sumShowerEvis;
      recoEnuVis = sumTrackEvis + sumShowerEvis;
      recoLepEFromTracks = sumTrackE;
      recoLepEFromShowers = sumShowerEvis;

      // Use the CAF calorimetric estimator when it is available.  The present
      // SPINE-only production stores NaN here, so the documented track/shower
      // fallback provides the usable reconstructed quantities.
      if (std::isfinite(enuCalo[iInt]) && enuCalo[iInt] > 0.f) {
        recoEnuCalo = enuCalo[iInt];
        recoEnu = recoEnuCalo;
        caloVsTracksRelDiff = (recoEnuFromTracksOnly - recoEnuCalo) / recoEnuCalo;
        caloVsTracksShowersRelDiff =
          (recoEnuFromTracksPlusShowers - recoEnuCalo) / recoEnuCalo;
      }
      else {
        recoEnu = numu_CC ? recoEnuFromTracksOnly : recoEnuFromTracksPlusShowers;
        usedRecoEnuFallback = 1;
      }

      if (nue_CC && std::isfinite(enuECalo[iInt]) && enuECalo[iInt] >= 0.f) {
        recoLepE = enuECalo[iInt];
      }
      else if (numu_CC && std::isfinite(enuLepCalo[iInt]) && enuLepCalo[iInt] >= 0.f) {
        recoLepE = enuLepCalo[iInt];
      }
      else {
        recoLepE = nue_CC ? recoLepEFromShowers : recoLepEFromTracks;
        usedRecoLepEFallback = 1;
      }
      if (recoEnu > 0.f && recoLepE >= 0.f) yRec = 1.f - recoLepE / recoEnu;
      if (numu_CC) recoEnu_numu = recoEnuFromTracksOnly;
      if (nue_CC) recoEnu_nue = recoEnuFromTracksPlusShowers;

      if (sample == kRejected) continue;

      if (numu_CC) ++nMu;
      if (nue_CC) ++nE;

      matchedMC = BestMatchedTruth(iInt,
                                   recoTruthLength, recoTruthIdx,
                                   recoTruth, recoTruthOverlap,
                                   nMC, truthTotal,
                                   matchedTruthOverlap);

      // Clear the branches for every output entry.  An unmatched
      // reconstructed interaction is retained with empty dial arrays and
      // matched_mc = -1.
      for (size_t dial = 0; dial < flux_dials.size(); ++dial) {
        flux_dials[dial]->Clear("C");
      }
      for (size_t dial = 0; dial < xsec_dials.size(); ++dial) {
        xsec_dials[dial]->Clear("C");
      }
      xsec = -999.f;
      // Neutral value for retained reco interactions without a truth match.
      // matched_mc == -1 identifies those entries downstream.
      xsec_cvwgt = 1.f;
      const bool printThisOutput =
          kDebugOutputEntries > 0 && outtree.GetEntries() < kDebugOutputEntries;

      if (matchedMC >= 0) {
        xsec = mcXSec[matchedMC];
        xsec_cvwgt = mcXSecCVWgt[matchedMC];
        if (!std::isfinite(xsec_cvwgt)) {
          std::cerr << "Non-finite xsec_cvwgt in event " << i
                    << ", interaction " << iInt << "." << std::endl;
          return 3;
        }
        mcEnu = Ev = mcE[matchedMC];
        mcModeOut = mcMode[matchedMC];
        isCC = mcIsCC[matchedMC] ? 1 : 0;
        isNC = isCC ? 0 : 1;
        nuPDG = mcPDG[matchedMC];
        nuPDGunosc = mcPDGOrig[matchedMC];
        iniNuFlavor = mcPDGOrig[matchedMC];
        interNuFlavor = mcPDG[matchedMC];
        detNuFlavor = numu_CC ? 14 : (nue_CC ? 12 : 0);
        mcGenWeightOut = mcGenWeight[matchedMC];
        const double px = mcMomX[matchedMC];
        const double py = mcMomY[matchedMC];
        const double pz = mcMomZ[matchedMC];
        const double momentum = std::sqrt(px * px + py * py + pz * pz);
        if (momentum > 0.) {
          const double cosTheta = std::max(-1., std::min(1., pz / momentum));
          mcthetaZ = static_cast<Float_t>(std::acos(cosTheta));
        }

        // Truth-level antineutrino labels; the nue_CC and numu_CC branches
        // remain the reconstructed selection categories.
        anue_CC = isCC && nuPDG == -12;
        anumu_CC = isCC && nuPDG == -14;

        if (printThisOutput) {
          std::cout << "[DEBUG] output=" << outtree.GetEntries()
                    << " inputEvent=" << i
                    << " recoInteraction=" << iInt
                    << " matchedMC=" << matchedMC
                    << " overlap=" << matchedTruthOverlap
                    << " xsec_cvwgt=" << xsec_cvwgt
                    << " (stored only; not applied to nominal or graph y values)\n";
        }

        const bool fluxOK = FillSeparateDialBranches(
            flux_dials,
            fluxMetadata,
            fluxGraphNames,
            fluxDialLength[matchedMC],
            fluxDialIdx[matchedMC],
            fluxDialTotal,
            fluxWeightLength,
            fluxWeightIdx,
            fluxWeightTotal,
            fluxWeights,
            "flux",
            printThisOutput
        );

        const bool xsecOK = FillSeparateDialBranches(
            xsec_dials,
            xsecMetadata,
            xsecGraphNames,
            xsecDialLength[matchedMC],
            xsecDialIdx[matchedMC],
            xsecDialTotal,
            xsecWeightLength,
            xsecWeightIdx,
            xsecWeightTotal,
            xsecWeights,
            "xsec",
            printThisOutput
        );

        if (!fluxOK || !xsecOK) {
          std::cerr << "Invalid systematic-array range in event " << i
                    << ", interaction " << iInt << "." << std::endl;
          return 3;
        }
      }
      else if (printThisOutput) {
        std::cout << "[DEBUG] output=" << outtree.GetEntries()
                  << " inputEvent=" << i
                  << " recoInteraction=" << iInt
                  << " has no matched truth interaction; storing "
                  << "xsec_cvwgt=1 and empty systematic graphs\n";
      }

      outtree.Fill();
    }
  }

  output.cd();
  outtree.Write();
  output.Close();
  for (size_t dial = 0; dial < flux_dials.size(); ++dial) delete flux_dials[dial];
  for (size_t dial = 0; dial < xsec_dials.size(); ++dial) delete xsec_dials[dial];
  std::cout << "Selected " << outtree.GetEntries() << " interactions: "
            << nMu << " numu-like, " << nE << " nue-like." << std::endl;
  return 0;
}
