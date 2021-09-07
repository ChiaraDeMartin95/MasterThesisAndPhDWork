// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Example V0 tutorial task 1: a simple looper
//
// Step 1: add extra histogram for lambda and antilambda mass
// Step 2: add extra 3D histogram with masses, momenta and centralities
// Step 3: add configurable TPC dE/dx selection
// Step 4: vary selections
// Step 5: try with finder instead of builder
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "AnalysisCore/RecoDecay.h"
#include "AnalysisCore/trackUtilities.h"
#include "AnalysisDataModel/StrangenessTables.h"
#include "AnalysisCore/TrackSelection.h"
#include "AnalysisDataModel/TrackSelectionTables.h"
#include "AnalysisDataModel/EventSelection.h"
#include "AnalysisDataModel/Centrality.h"
#include "AnalysisDataModel/PID/PIDResponse.h"

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct lfstrangetut1 {

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {

  AxisSpec vtxZAxis = {100, -20, 20};
  AxisSpec centAxis = {100, 0, 100, "V0M (%)"};
  AxisSpec massK0sAxis = {100, 0.450f, 0.550f};
  AxisSpec massLambdaAxis = {200, 1.08f, 1.18f};
  AxisSpec ptAxis = {100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
  //  std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
  //  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
  //  std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
  //AxisSpec centAxis = {centBinning, "V0M (%)"};

  histos.add("NumberV0", "Number of selected V0", kTH1F, {{20, 0.5, 20.5}});
  histos.add("Centrality", "Centrality distribution (V0M)", kTH1F, {centAxis});
  histos.add("VtxZAfterSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});
  histos.add("hMassK0ShortBefSel", "hMassK0ShortBefSel", kTH1F, {massK0sAxis});
  histos.add("hMassK0ShortAfterSel", "hMassK0ShortAfterSel", kTH1F, {massK0sAxis});
  histos.add("hMassK0ShortAfterSel2D", "hMassK0ShortAfterSel2D", kTH2F, {massK0sAxis, ptAxis});
  histos.add("hMassLambdaBefSel", "hMassLambdaBefSel", kTH1F, {massLambdaAxis});
  histos.add("hMassLambdaAfterSel", "hMassLambdaAfterSel", kTH1F, {massLambdaAxis});
  histos.add("hMassLambdaAfterSel2D", "hMassLambdaAfterSel2D", kTH2F, {massLambdaAxis, ptAxis});
  histos.add("hv0cospa", "hv0cospa", kTH1F, {{100, 0.9, 1}});
  histos.add("hv0radius", "hv0radius", kTH1F, {{1000, 0, 50}});

  }

  //Selection criteria
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.9, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> minpt{"minpt", 0.1, "minpt"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  using CollisionCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Cents>>::iterator;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;   

  void process(CollisionCandidates const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, DaughterTracks &tracks )
  {
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }

    histos.fill(HIST("VtxZAfterSel"), collision.posZ());
    histos.fill(HIST("Centrality"), collision.centV0M());

    for (auto& v0 : fullV0s) {

      histos.fill(HIST("NumberV0"), 1);
      histos.fill(HIST("hMassK0ShortBefSel"), v0.mK0Short());
      histos.fill(HIST("hMassLambdaBefSel"), v0.mLambda());
      if (v0.v0radius() < v0radius) continue;
      histos.fill(HIST("hv0radius"), v0.v0radius());
      histos.fill(HIST("NumberV0"), 2);
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) continue;
      histos.fill(HIST("hv0cospa"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      histos.fill(HIST("NumberV0"), 3);

      if (TMath::Abs(v0.yK0Short()) < rapidity && fabs(v0.posTrack_as<DaughterTracks>().tpcNSigmaPi()) < 3.0 && fabs(v0.negTrack_as<DaughterTracks>().tpcNSigmaPi()) <3.0) {
	histos.fill(HIST("NumberV0"), 4);
	histos.fill(HIST("hMassK0ShortAfterSel"), v0.mK0Short());
	histos.fill(HIST("hMassK0ShortAfterSel2D"), v0.mK0Short(), v0.pt());
      }
      if (TMath::Abs(v0.yLambda()) < rapidity && fabs(v0.posTrack_as<DaughterTracks>().tpcNSigmaPr()) < 3.0 && fabs(v0.negTrack_as<DaughterTracks>().tpcNSigmaPi()) <3.0) {
	histos.fill(HIST("NumberV0"), 5);
	histos.fill(HIST("hMassLambdaAfterSel"), v0.mLambda());
	histos.fill(HIST("hMassLambdaAfterSel2D"), v0.mLambda(), v0.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lfstrangetut1>(cfgc, TaskName{"lf-lfstrangetut1"})};
}
