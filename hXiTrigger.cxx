// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
// O2 includes

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "AnalysisDataModel/PID/PIDResponse.h"
#include "AnalysisDataModel/TrackSelectionTables.h"
#include "AnalysisDataModel/Centrality.h"
#include "AnalysisDataModel/EventSelection.h"
#include "AnalysisDataModel/StrangenessTables.h"
#include "AnalysisCore/TrackSelection.h"
#include "AnalysisCore/RecoDecay.h"
#include "AnalysisCore/trackUtilities.h"
#include "filterTables.h"
#include "Framework/HistogramRegistry.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

namespace
{
  float rapidity(float pt, float eta, float m)
  {
    return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
  }

  //  static constexpr std::array<float, nNuclei> masses{
  //  constants::physics::MassDeuteron, constants::physics::MassTriton,
  //  constants::physics::MassHelium3, constants::physics::MassAlpha};
  //  static constexpr std::array<int, nNuclei> charges{1, 1, 2, 2};
  //  static const std::vector<std::string> nucleiNames{"H2", "H3", "He3", "He4"};

} // namespace

struct hXiTrigger {

  Produces<aod::hXiTriggers> tags;

  Configurable<float> yMin{"yMin", -0.8, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.8, "Minimum rapidity"};
  Configurable<float> yBeam{"yBeam", 0., "Beam rapidity"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutMinPt{"cfgCutEta", 8.0f, "Min pt for trigger particles"};

  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.9, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> minpt{"minpt", 1.0, "minpt"};


  //OLD  Configurable<LabeledArray<float>> cfgCutsPID{"nucleiCutsPID", {cutsPID[0], nNuclei, nCutsPID, nucleiNames, cutsNames}, "Nuclei PID selections"};

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    //OLD    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5.};
    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 40., 50., 70., 100.};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {centBinning, "V0M (%)"};

    spectra.add("fCollZpos", "collision z position", HistType::kTH1F, {{600, -20., +20., "z position (cm)"}});
    spectra.add("fProcessedEvents", "Selected events", HistType::kTH1F, {{4, -0.5, 3.5, "Event counter"}});
    spectra.add("PtTrigger", "PtTrigger",  HistType::kTH1F, {{300, 0, 30, "Pt of trigger particle"}});
    spectra.add("hMassXi", "hMassXi",  HistType::kTH1F, {{100, 1.28, 1.36, "Inv mass Xi candidates"}});
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  //other collision selections: no pile up events, events with multiplicity defined

  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::isGlobalTrack == static_cast<uint8_t>(1u)) && (aod::track::pt>cfgCutMinPt);
  //other selections for trigger particles: quality criteria, |DCAz| < , |DCAxy| < 

  //old  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl, aod::pidTOFDe, aod::pidTOFTr, aod::pidTOFHe, aod::pidTOFAl, aod::TrackSelection>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Cents>>::iterator const& collision, aod::BCsWithTimestamps const&, TrackCandidates const& tracks)
  {
    // collision process loop
    bool keepEvent{false};

    spectra.fill(HIST("fCollZpos"), collision.posZ());
    //    spectra.fill(HIST("fCollMultiplicity"), collision.Mult());

    for (auto track : tracks) { // start loop over tracks

      //      const float nSigmaTPC[nNuclei]{
      //  track.tpcNSigmaDe(), track.tpcNSigmaTr(), track.tpcNSigmaHe(), track.tpcNSigmaAl()};
      //const float nSigmaTOF[nNuclei]{
      //  track.tofNSigmaDe(), track.tofNSigmaTr(), track.tofNSigmaHe(), track.tofNSigmaAl()};

      //      if ( ) continue;
      keepEvent = true;

      //
      // fill QA histograms
      //
      spectra.fill(HIST("PtTrigger"), track.pt());

    } // end loop over tracks

    /*
    if (keepEvent) {
      keepEvent = kFALSE;
      //loop over Cascade candidates

      if () continue;
      keepEvent = true;
    }
    */

    if (keepEvent) {
      spectra.fill(HIST("fProcessedEvents"), 1);
    }

    tags(keepEvent);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<hXiTrigger.cxx>(cfg)};
}
