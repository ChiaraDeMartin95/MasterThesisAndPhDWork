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
/// \brief A filter task for strangeness filter
//  usage:
/*  o2-analysis-timestamp --aod-file ~/Scaricati/AO2D_PbPb.root -b | \
    o2-analysis-multiplicity-table -b | \
    o2-analysis-centrality-table -b | \
    o2-analysis-event-selection -b | \
    o2-analysis-trackextension -b | \
    o2-analysis-pid-tpc | \
    o2-analysis-pid-tof | \
    o2-analysis-weak-decay-indices -b | \
    o2-analysis-lambdakzerobuilder --d_bz 5 -b | \
    o2-analysis-cascadebuilder --d_bz 5 -b | \
    o2-analysis-strangeness-filter -b 
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

struct strangenessFilter {

  //Recall the output table
  Produces<aod::StrangenessFilters> strgtable;

  //Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> hProcessedEvents{TH1F("hProcessedEvents", "Strangeness - event filtered; Event counter; Number of events", 3, 0., 3.)};

  //Selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> v0cospa{"v0cospa", 0.97, "V0 CosPA"}; //is it with respect to Xi decay vertex?
  Configurable<float> casccospa{"casccospa", 0.995, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};       //is it in sigmas?
  Configurable<float> dcacascdau{"dcacascdau", 0.8, "DCA Casc Daughters"}; //is it in sigmas?
  Configurable<float> dcamesontopv{"dcamesontopv", 0.04, "DCA Meson To PV"};
  Configurable<float> dcabaryontopv{"dcabaryontopv", 0.03, "DCA Baryon To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.04, "DCA Bach To PV"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.02, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.02, "DCA Pos To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 1.2, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 1.2, "V0 Radius"};
  Configurable<float> v0radiusupperlimit{"v0radiusupperlimit", 34, "V0 Radius Upper Limit"};
  Configurable<float> cascradius{"cascradius", 0.6, "cascradius"};
  Configurable<float> cascradiusupperlimit{"cascradiusupperlimit", 34, "Casc Radius Upper Limit"};
  Configurable<float> rapidity{"rapidity", 2, "rapidity"};
  Configurable<float> eta{"eta", 2, "Eta"};
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.8, "EtaDaughters"};
  Configurable<float> masslambdalimit{"masslambdalimit", 0.01, "masslambdalimit"}; //0.006 Chiara
  Configurable<float> omegarej{"omegarej", 0.005, "omegarej"};
  Configurable<float> ximasswindow{"ximasswindow", 0.075, "Xi Mass Window"};
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
    QAHistos.add("hMassXiBefSel", "hMassXiBefSel", HistType::kTH1F, {{100, 1.30f, 1.34f}});
    QAHistos.add("hMassXiAfterSel", "hMassXiAfterSel", HistType::kTH1F, {{100, 1.30f, 1.34f}});
    //QAHistos.add("hMassXiAfterSelvsPt", "hMassXiAfterSelvsPt", HistType::kTH2F, {{100, 1.30f, 1.34f}}, {{100, 0f, 10f, "#it{p}_{T} (GeV/#it{c})"}}); //add topological variables histos
    QAHistos.add("hTriggeredParticles", "Selected triggered particles", HistType::kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
    QAHistos.add("PtTrigger", "PtTrigger", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particle"}});

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "#Xi-#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "high-#it{p}_{T} hadron - #Xi");
  }

  //Filters
  Filter collisionFilter = (nabs(aod::collision::posZ) < cutzvertex);
  Filter preFilterCasc = nabs(aod::cascdata::dcapostopv) > dcapostopv&& nabs(aod::cascdata::dcanegtopv) > dcanegtopv&& aod::cascdata::dcaV0daughters < dcav0dau&& aod::cascdata::dcacascdaughters < dcacascdau;

  //Tables
  using CollisionCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Cents>>::iterator;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;

  void process(CollisionCandidates const& collision, aod::CascDataExt const& fullCasc, DaughterTracks& tracks)
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

    //Is event good? [0] = DoubleXi, [1] = high-pT hadron + Xi
    bool keepEvent[2]{false};
    //
    float xipos = -1.;
    float xiproperlifetime = -1.;
    float xiptotmom = -1.;
    int xicounter = 0;
    const float ctauxi = 4.91;    //from PDG
    const float massxi = 1.32171; //from PDG

    for (auto& casc : fullCasc) { //loop over cascades

      QAHistos.fill(HIST("hMassXiBefSel"), casc.mXi());

      //Position
      xipos = TMath::Sqrt(TMath::Power(casc.x() - collision.posX(), 2) + TMath::Power(casc.y() - collision.posY(), 2) + TMath::Power(casc.z() - collision.posZ(), 2));
      //Total momentum
      xiptotmom = TMath::Sqrt(casc.px() * casc.px() + casc.py() * casc.py() + casc.pz() * casc.pz());
      //Proper lifetime
      xiproperlifetime = massxi * xipos / (xiptotmom + 1e-13);

      //Selections
      /*if (casc.sign() == 1) {
         if (TMath::Abs(casc.posTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc)
          continue;
        if (TMath::Abs(casc.negTrack_as<DaughterTracks>().tpcNSigmaPr()) > nsigmatpc)
          continue;
      } else {
        if (TMath::Abs(casc.posTrack_as<DaughterTracks>().tpcNSigmaPr()) > nsigmatpc)
          continue;
        if (TMath::Abs(casc.negTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc)
          continue;
      }
      if (TMath::Abs(casc.bachTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc)
        continue; //?

      if (TMath::Abs(casc.bachTrack_as<DaughterTracks>().eta()) > etadau)
        continue; //?
      if (TMath::Abs(casc.posTrack_as<DaughterTracks>().eta()) > etadau)
        continue; //?
      if (TMath::Abs(casc.negTrack_as<DaughterTracks>().eta()) > etadau)
        continue; //?
      */
      if (casc.sign() == 1) {
        if (TMath::Abs(casc.dcapostopv()) < dcabaryontopv)
          continue;
        if (TMath::Abs(casc.dcanegtopv()) < dcamesontopv)
          continue;
      } else {
        if (TMath::Abs(casc.dcanegtopv()) < dcabaryontopv)
          continue;
        if (TMath::Abs(casc.dcapostopv()) < dcamesontopv)
          continue;
      }
      if (TMath::Abs(casc.dcabachtopv()) < dcabachtopv)
        continue;
      if (casc.v0radius() > v0radiusupperlimit || casc.v0radius() < v0radius)
        continue;
      if (casc.cascradius() > cascradiusupperlimit || casc.cascradius() < cascradius)
        continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa)
        continue; //is it calculated with respect to the Xi decay vertex?
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < casccospa)
        continue;
      if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < dcav0topv)
        continue;
      if (casc.dcaV0daughters() > dcav0dau)
        continue;
      if (casc.dcacascdaughters() > dcacascdau)
        continue;
      if (TMath::Abs(casc.mLambda() - constants::physics::MassLambda) > masslambdalimit)
        continue;
      if (TMath::Abs(casc.mXi() - massxi) < ximasswindow) //if (TMath::Abs(casc.mXi() - constants::physics::Massxi) < ximasswindow) continue;
        continue;
      if (TMath::Abs(casc.mOmega() - 1.67245) < omegarej)
        continue; //To be used constants::physics::MassOmega
      if (xiproperlifetime > properlifetimefactor * ctauxi)
        continue;
      if (TMath::Abs(casc.yXi()) > rapidity)
        continue;
      if (TMath::Abs(casc.eta()) > eta)
        continue;

      QAHistos.fill(HIST("hMassXiAfterSel"), casc.mXi());
      QAHistos.fill(HIST("hMassXiAfterSel2D"), casc.mXi(), casc.pt());

      //Count number of Xi candidates
      xicounter++;

    } //end loop over cascades

    //Double Xi trigger definition
    if (xicounter > 1) {
      keepEvent[0] = true;
    }

    //High-pT hadron + Xi trigger definition
    if (xicounter > 0) {
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
        if (float(track.tpcNClsCrossedRows()) / track.length() < 0.8)
          continue;
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
    strgtable(keepEvent[0],keepEvent[1]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilter>(cfgc, TaskName{"strangeness-filter"})};
}
