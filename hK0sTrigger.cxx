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

const float ctauK0s = 2.6844;   //from PDG                                                                                                    
const float massK0s = 0.497611; //from PDG    

struct hK0sTrigger {

  //Produces<aod::hK0sTriggers> tags;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histostracks{"HistosTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosK0s{"HistosK0s", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutMinPt{"cfgCutMinPt", 1.0f, "Min pt for trigger particles"};

  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"}; //
  Configurable<float> dcanegtopv{"dcanegtopv", .06, "DCA Neg To PV"}; //
  Configurable<float> dcapostopv{"dcapostopv", .06, "DCA Pos To PV"}; //
  Configurable<float> dcav0topv{"dcav0topv", .5, "DCA V0 To PV"}; //
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> v0radius{"v0radius", 0.9, "v0radius"}; //
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> eta{"eta", 0.8, "eta"}; //
  Configurable<float> ctau{"ctau", 20, "ctau"};
  Configurable<float> LRej{"LRej", 0.005, "LRej"}; //
  Configurable<float> minpt{"minpt", 0, "minpt"}; //
  Configurable<int> properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> nsigmatpc{"nsigmatpc", 6, "N Sigmas TPC"};

  void init(o2::framework::InitContext&)
  {

    AxisSpec vtxZAxis = {100, -20, 20};
    //    AxisSpec centAxis = {100, 0, 100, "V0M (%)"};
    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 40., 50., 70., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec massK0sAxis = {100, 0.450f, 0.550f};
    AxisSpec ptAxis = {100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("Centrality", "Centrality distribution (V0M)", kTH1F, {centAxis});
    histos.add("CentralityTriggered", "Centrality distribution of events with trigger particle h (V0M)", kTH1F, {centAxis});
    histos.add("CentralityFullTrigger", "Centrality distribution of events with h + Xi (V0M)", kTH1F, {centAxis});
    histos.add("VtxZAfterSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});
    histos.add("fProcessedEvents", "Selected events", HistType::kTH1F, {{4, 0.5, 4.5, "Event counter"}});

    histostracks.add("hTriggeredParticles", "Selected triggered particles", HistType::kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
    histostracks.add("hLengthvsCrossedRows", "hLengthvsCrossedRows", kTH2F, {{300, 0, 300}, {150, 0, 150}});
    histostracks.add("hTrackLength", "hTrackLength", kTH1F, {{1000, 0, 1000}});
    histostracks.add("PtTrigger", "PtTrigger",  HistType::kTH1F, {{300, 0, 30, "Pt of trigger particle"}});

    histosK0s.add("hMassK0ShortBefSel", "hMassK0ShortBefSel", kTH1F, {massK0sAxis});
    histosK0s.add("hMassK0ShortAfterSel", "hMassK0ShortAfterSel", kTH1F, {massK0sAxis});
    histosK0s.add("hptK0s", "hptK0s", kTH1F, {{100, 0, 10}});
    histosK0s.add("hctau", "hctau", kTH1F, {{250, 0, 25}});
    histosK0s.add("hdcav0dau", "hdcav0dau", kTH1F, {{20, 0, 2}});
    histosK0s.add("hdcanegtopv", "hdcanegtopv", kTH1F, {{10, 0, 0.1}});
    histosK0s.add("hdcapostopv", "hdcapostopv", kTH1F, {{10, 0, 0.1}}); 
    histosK0s.add("hdcav0topv", "hdcav0topv", kTH1F, {{50, 0, 1}});
    histosK0s.add("hv0cospa", "hv0cospa", kTH1F, {{50,0.95, 1}});
    histosK0s.add("hv0radius", "hv0radius", kTH1F, {{100, 0, 5}});
    histosK0s.add("hrapidity", "hrapidity", kTH1F, {{100, -1, 1}});
    histosK0s.add("heta", "heta", kTH1F, {{100, -1, 1}});

  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  //other collision selections: no pile up events
  
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::isGlobalTrack == static_cast<uint8_t>(1u)) && (aod::track::pt>cfgCutMinPt);

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv   && aod::v0data::dcaV0daughters < dcav0dau;

  using CollisionCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Cents>>::iterator;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>>;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;

  void process(CollisionCandidates const& collision, TrackCandidates const& tracks, soa::Filtered<aod::V0Datas> const& fullV0s, DaughterTracks &dtracks )
  {
    // collision process loop
    bool keepEvent{false};
    float ctauK0s{0};

    histos.fill(HIST("fProcessedEvents"), 1);

    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }

    histos.fill(HIST("VtxZAfterSel"), collision.posZ());
    histos.fill(HIST("Centrality"), collision.centV0M());
    histos.fill(HIST("fProcessedEvents"), 2);

    for (auto track : tracks) { // start loop over tracks
      histostracks.fill(HIST("hTriggeredParticles"), 1);
      /*
      if (track.tpcChi2NCl()>4) continue; 
      histostracks.fill(HIST("hTriggeredParticles"), 2);
      if (track.tpcNClsCrossedRows() < 80) continue; //~no changes
      histostracks.fill(HIST("hTriggeredParticles"), 3);
      if (track.tpcCrossedRowsOverFindableCls() < 0.8) continue; //no changes
      histostracks.fill(HIST("hTriggeredParticles"), 4);
      histostracks.fill(HIST("hTrackLength"), track.length());
      histostracks.fill(HIST("hLengthvsCrossedRows"), track.tpcNClsCrossedRows(), track.length());
      //if (track.length() < 90) continue; //some decrease observed, is the method correct? I get lenghts =0 or 370 < l < 550...
      histostracks.fill(HIST("hTriggeredParticles"), 5);
      //if (float(track.tpcNClsCrossedRows())/track.length() < 0.8) continue; //no track passes this selection!
      histostracks.fill(HIST("hTriggeredParticles"), 6);
      */
      histostracks.fill(HIST("PtTrigger"), track.pt());
      keepEvent = true; // a trigger particle is found
    } // end loop over tracks

    if (!keepEvent) return; 
    histos.fill(HIST("fProcessedEvents"), 3);
    histos.fill(HIST("CentralityTriggered"), collision.centV0M());

    keepEvent=false;
    float v0pos = -1.;
    float K0sproperlifetime = -1.;
    float v0ptotmom = -1.;

    for (auto& v0 : fullV0s) {

      //Position                                                                     
      v0pos = TMath::Sqrt(TMath::Power(v0.x() - collision.posX(), 2) + TMath::Power(v0.y() - collision.posY(), 2) + TMath::Power(v0.z() - collision.posZ(), 2));
      //Total momentum                                                                                                                           
      v0ptotmom = TMath::Sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
      //Proper lifetime                                                                                                                          
      K0sproperlifetime = massK0s * v0pos / (v0ptotmom + 1e-13);

      histosK0s.fill(HIST("hMassK0ShortBefSel"), v0.mK0Short());

      if (TMath::Abs(v0.posTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc) continue;
      if (TMath::Abs(v0.negTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc) continue; 
      //-----------------------
      // TOPOLOGICAL - KINEMATIC SELECTIONS
      //-----------------------
      if (v0.pt() < minpt) continue;
      if (TMath::Abs(v0.eta()) > eta) continue;
      if (v0.v0radius() < v0radius) continue;
      if (v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) continue;
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) continue;
      if (TMath::Abs(v0.mLambda() - constants::physics::MassLambda) < LRej) continue;
      if (K0sproperlifetime > properlifetimefactor * ctauK0s) continue;
      //OOB pile up: tofSignal or HasTOF

      histosK0s.fill(HIST("hptK0s"), v0.pt());
      histosK0s.fill(HIST("hctau"), K0sproperlifetime);
      histosK0s.fill(HIST("hMassK0ShortAfterSel"), v0.mK0Short());
      histosK0s.fill(HIST("hdcav0dau"), v0.dcaV0daughters());
      histosK0s.fill(HIST("hdcanegtopv"), v0.dcanegtopv());
      histosK0s.fill(HIST("hdcapostopv"), v0.dcapostopv());
      histosK0s.fill(HIST("hdcav0topv"), v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
      histosK0s.fill(HIST("hv0cospa"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      histosK0s.fill(HIST("hv0radius"), v0.v0radius());
      histosK0s.fill(HIST("hrapidity"), v0.yK0Short());
      histosK0s.fill(HIST("heta"), v0.eta());

      keepEvent = true;
    }

    if (keepEvent) {
      histos.fill(HIST("CentralityFullTrigger"), collision.centV0M());
      histos.fill(HIST("fProcessedEvents"), 4);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<hK0sTrigger>(cfg, TaskName{"hK0s-trigger"})};
}
