// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief A filter task for strangeness filter
//  usage:
/*  o2-analysis-timestamp --aod-file /data/dataalice/cdemart/O2/AO2D_15o_t180.root -b | \
    o2-analysis-multiplicity-table -b | \
    o2-analysis-centrality-table -b | \
    o2-analysis-event-selection -b | \
    o2-analysis-trackextension -b | \
    o2-analysis-pid-tpc | \
    o2-analysis-pid-tof | \
    o2-analysis-weak-decay-indices -b | \
    o2-analysis-lambdakzerobuilder --d_bz 5 -b | \
    o2-analysis-strangeness-filter-K0s -b 
*/
///
///
/// \author Chiara De Martin (chiara.de.martin@cern.ch) and Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since June 1, 2021

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
#include "filterTables.h"

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

struct strangenessFilterK0s {

  //Recall the output table
  Produces<aod::StrangenessFiltersK0s> strgtableK0s;

  //Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> hProcessedEvents{TH1F("hProcessedEvents", "Strangeness - event filtered; Event counter; Number of events", 3, 0., 3.)};

  //Selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //is it with respect to Xi decay vertex?
  Configurable<float> dcav0dau{"dcav0dau", 1, "DCA V0 Daughters"};       //is it in sigmas?
  Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 0.5, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 0.9, "V0 Radius"};
  Configurable<float> v0radiusupperlimit{"v0radiusupperlimit", 34, "V0 Radius Upper Limit"};
  Configurable<float> LRej{"LRej", 0.005, "LRej"}; 
  Configurable<float> rapidity{"rapidity", 2, "rapidity"};
  Configurable<float> eta{"eta", 2, "Eta"};
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.8, "EtaDaughters"};
  Configurable<int> properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> nsigmatpc{"nsigmatpc", 6, "N Sigmas TPC"};
  //missing selections: OOB pileup?
  //eta and y selections: loose enough?
  //eta selections of daughters

  void init(o2::framework::InitContext&)
  {

    //  std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    //  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    //  std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    //  AxisSpec centAxis = {centBinning, "V0M (%)"};

    QAHistos.add("Centrality", "Centrality distribution (V0M)", HistType::kTH1F, {{100, 0, 100, "V0M (%)"}});
    QAHistos.add("VtxZAfterSel", "Vertex distribution in Z;Z (cm)", HistType::kTH1F, {{100, -20, 20}});
    QAHistos.add("hMassK0sBefSel", "hMassK0sBefSel", HistType::kTH1F, {{100, 0.49f, 0.506f}});
    QAHistos.add("hMassK0sAfterSel", "hMassK0sAfterSel", HistType::kTH1F, {{100, 0.49f, 0.506f}});
    //QAHistos.add("hMassK0sAfterSelvsPt", "hMassK0sAfterSelvsPt", HistType::kTH2F, {{100, 1.30f, 1.34f}}, {{100, 0f, 10f, "#it{p}_{T} (GeV/#it{c})"}}); //add topological variables histos
    QAHistos.add("hTriggeredParticles", "Selected triggered particles", HistType::kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
    QAHistos.add("PtTrigger", "PtTrigger", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particle"}});

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "K0s-K0s");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "high-#it{p}_{T} hadron - K0s");
  }

  //Filters
  //  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::isGlobalTrack == static_cast<uint8_t>(1u)) && (aod::track::pt>cfgCutMinPt);
  Filter collisionFilter = (nabs(aod::collision::posZ) < cutzvertex);
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  //Tables
  
  using CollisionCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Cents>>::iterator;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;

  void process(CollisionCandidates const& collision,  TrackCandidates const& tracks, soa::Filtered<aod::V0Datas> const& fullV0, DaughterTracks& dtracks)
 
  {
    
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) { //what's that?
      return;
    }

    QAHistos.fill(HIST("VtxZAfterSel"), collision.posZ());
    QAHistos.fill(HIST("Centrality"), collision.centV0M());
    hProcessedEvents->Fill(0.5);

    //Is event good? [0] = DoubleK0s, [1] = high-pT hadron + K0s
    bool keepEvent[2]{false};
    //
    float v0pos = -1.;
    float K0sproperlifetime = -1.;
    float v0ptotmom = -1.;
    int K0scounter = 0;
    const float ctauK0s = 2.6844;   //from PDG
    const float massK0s = 0.497611; //from PDG

    
    for (auto& v0 : fullV0) { //loop over V0s

      QAHistos.fill(HIST("hMassK0sBefSel"), v0.mK0Short());

      //Position
      v0pos = TMath::Sqrt(TMath::Power(v0.x() - collision.posX(), 2) + TMath::Power(v0.y() - collision.posY(), 2) + TMath::Power(v0.z() - collision.posZ(), 2));
      //Total momentum
      v0ptotmom = TMath::Sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
      //Proper lifetime
      K0sproperlifetime = massK0s * v0pos / (v0ptotmom + 1e-13);

      if (TMath::Abs(v0.posTrack_as<DaughterTracks>().tpcNSigmaPi()) > 3.0) continue;
      if (TMath::Abs(v0.negTrack_as<DaughterTracks>().tpcNSigmaPi()) > 3.0) continue;
      //-----------------------                                                                                                                  
      // TOPOLOGICAL - KINEMATIC SELECTIONS                                                                                                      
      //-----------------------                                                                                                                  
      if (TMath::Abs(v0.eta()) > eta) continue;
      if (v0.v0radius() < v0radius) continue;
      if (v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) continue;
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) continue;
      if (TMath::Abs(v0.mLambda() - constants::physics::MassLambda) < LRej) continue;
      if (K0sproperlifetime > properlifetimefactor * ctauK0s) continue;
      
      QAHistos.fill(HIST("hMassK0sAfterSel"), v0.mK0Short());

      //Count number of K0s candidates
      K0scounter++;

    } //end loop over V0s
    
    //Double K0s trigger definition
    if (K0scounter > 1) {
      keepEvent[0] = true;
    }

    
    //High-pT hadron + Xi trigger definition
    if (K0scounter > 0) {
      for (auto track : tracks) { // start loop over tracks
        QAHistos.fill(HIST("hTriggeredParticles"), 1);
        if (track.itsChi2NCl() > 4)
          continue; //some decrease observed
        QAHistos.fill(HIST("hTriggeredParticles"), 2);
        if (track.tpcNClsCrossedRows() < 80)
          continue; //~no changes
        QAHistos.fill(HIST("hTriggeredParticles"), 3);
        if (track.tpcCrossedRowsOverFindableCls() < 0.8)
          continue; //no changes
        QAHistos.fill(HIST("hTriggeredParticles"), 4);
        if (track.length() < 90)
          continue; //some decrease observed
        QAHistos.fill(HIST("hTriggeredParticles"), 5);
	//        if (float(track.tpcNClsCrossedRows()) / track.length() < 0.8)
	// continue;
        QAHistos.fill(HIST("hTriggeredParticles"), 6);
        QAHistos.fill(HIST("PtTrigger"), track.pt());

        keepEvent[1] = true;
      } // end loop over tracks
    }
    
    if (keepEvent[0])
      hProcessedEvents->Fill(1.5);
    if (keepEvent[1])
      hProcessedEvents->Fill(2.5);
    
    //Filling the table
    strgtableK0s(keepEvent[0],keepEvent[1]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilterK0s>(cfgc, TaskName{"strangeness-filter-K0s"})};
}
