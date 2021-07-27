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

  Produces<aod::StrangenessFilters> strgtable;
  HistogramRegistry histos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> hProcessedEvents{TH1F("hProcessedEvents", "Strangeness - event filtered; Event counter; Number of events", 3, 0., 3.)};

  void init(o2::framework::InitContext&)
  {

  AxisSpec vtxZAxis = {100, -20, 20};
  AxisSpec centAxis = {100, 0, 100, "V0M (%)"};
  AxisSpec massXiAxis = {100, 1.30f, 1.34f};
  AxisSpec ptAxis = {100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
  //  std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
  //  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
  //  std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
  //AxisSpec centAxis = {centBinning, "V0M (%)"};

  histos.add("Centrality", "Centrality distribution (V0M)", kTH1F, {centAxis});
  histos.add("VtxZAfterSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});

  histos.add("hMassXiBefSel", "hMassXiBefSel", kTH1F, {massXiAxis});
  histos.add("hMassXiAfterSel", "hMassXiAfterSel", kTH1F, {massXiAxis});
  histos.add("hMassXiAfterSelvsPt", "hMassXiAfterSelvsPt", kTH1F, {massXiAxis}, {ptAxis});
  //add topological variables histos

  histos.add("hTriggeredParticles", "Selected triggered particles", kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
  histos.add("PtTrigger", "PtTrigger", kTH1F, {{300, 0, 30, "Pt of trigger particle"}});

  hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
  hProcessedEvents->GetXaxis()->SetBinLabel(2, "DoubleXi");
  hProcessedEvents->GetXaxis()->SetBinLabel(2, "high pT hadron - Xi");
  }

  //Selection criteria
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};

  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"}; 
  //is it with respect to Xi decay vertex?
  Configurable<double> casccospa{"casccospa", 0.995, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};   //is it in sigmas?
  Configurable<float> dcaCascdau{"dcaCascdau", 0.8, "DCA Casc Daughters"};   //is it in sigmas?
  Configurable<float> dcaMesontopv{"dcaMesontopv", 0.04, "DCA Meson To PV"};
  Configurable<float> dcaBaryontopv{"dcaBaryontopv", 0.03, "DCA Baryon To PV"};
  Configurable<float> dcaBachtopv{"dcaBachtopv", 0.04, "DCA Bach To PV"};
  Configurable<float> dcaV0topv{"dcaV0topv", 1.2, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 1.2, "v0radius"};
  Configurable<float> v0radiusUpperLimit{"v0radiusUpperLimit", 34, "v0radius upper limit"};
  Configurable<float> cascradius{"cascradius", 0.6, "cascradius"};
  Configurable<float> cascradiusUpperLimit{"cascradiusUpperLimit", 34, "cascradius upper limit"};
  Configurable<float> rapidity{"rapidity", 2, "rapidity"}; 
  Configurable<float> eta{"eta", 2, "Eta"}; 
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.8, "EtaDaughters"}; 
  Configurable<float> massLambdaLimit{"masslambdaLimit", 0.008, "masslambdaLimit"}; //0.006 Chiara
  Configurable<float> omegarej{"omegarej", 0.005, "omegarej"};
  Configurable<float> ximassWindow{"ximassWindow", 0.075, "Xi Mass Window"};
  Configurable<int>   properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> nsigmatpc{"nsigmatpc", 6, "N Sigmas TPC"};
  //missing selection: OOB pileup?
  //eta and y selections: loose enough?
  //eta selections of daughters

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);

  //  Filter preFilterCasc = nabs(aod::cascdata::dcapostopv) > dcapostopv&& nabs(aod::cascdata::dcanegtopv) > dcanegtopv&& aod::cascdata::dcaV0daughters < dcav0dau && aod::cascdata::dcacascdaughters < dcacascdau;

  using CollisionCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Cents>>::iterator;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;   

  void process(CollisionCandidates const& collision, aod::CascDataExt const& fullCasc,  DaughterTracks &tracks )
  {
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) { //what's that?
      return;
    }

    histos.fill(HIST("VtxZAfterSel"), collision.posZ());
    histos.fill(HIST("Centrality"), collision.centV0M());
    hProcessedEvents->Fill(0.5);

    bool keepEvent[2]{false};

    float distxi = -1.;
    float properlifetimexi = -1.;
    float ptotmomxi = -1.;
    int   xicounter=0;
    const float ctauxi = 4.91;
    const float Massxi = 1.322;

    for (auto& casc : fullCasc) {
      histos.fill(HIST("hMassXiBefSel"), casc.mXi());

      distxi = TMath::Sqrt(
			     TMath::Power( casc.x() - collision.posX() , 2) +
			     TMath::Power( casc.y() - collision.posY() , 2) +
			     TMath::Power( casc.z() - collision.posZ() , 2)
			     );
      ptotmomxi = TMath::Sqrt( 
                              casc.px()*casc.px()+
                              casc.py()*casc.py()+
                              casc.pz()*casc.pz()
			      );
      properlifetimexi = constants::physics::MassXi*distxi/(ptotmomxi+1e-13);

      if (casc.sign() == 1){
	if (TMath::Abs(casc.posTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc) continue;
	if (TMath::Abs(casc.negTrack_as<DaughterTracks>().tpcNSigmaPr()) > nsigmatpc) continue;
      }
      else{
	if (TMath::Abs(casc.posTrack_as<DaughterTracks>().tpcNSigmaPr()) > nsigmatpc) continue;
	if (TMath::Abs(casc.negTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc) continue;
      }
      if (TMath::Abs(casc.bachTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc) continue; //?

      if (TMath::Abs(casc.bachTrack_as<DaughterTracks>().eta()) > etadau) continue; //?
      if (TMath::Abs(casc.posTrack_as<DaughterTracks>().eta()) > etadau) continue; //?
      if (TMath::Abs(casc.negTrack_as<DaughterTracks>().eta()) > etadau) continue; //?

      if (casc.sign() == 1){
	if (TMath::Abs(casc.dcapostopv()) < dcaBaryontopv ) continue;
	if (TMath::Abs(casc.dcanegtopv()) < dcaMesontopv ) continue;
      }
      else {
	if (TMath::Abs(casc.dcanegtopv()) < dcaBaryontopv ) continue;
	if (TMath::Abs(casc.dcapostopv()) < dcaMesontopv ) continue;
      }
      if (TMath::Abs(casc.dcabachtopv()) < dcaBachtopv) continue;
      if (casc.v0radius() > v0radiusUpperLimit || casc.v0radius() < v0radius) continue;
      if (casc.cascradius() > cascradiusUpperLimit || casc.cascradius() < cascradius) continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) continue; //is it calculated with respect to the Xi decay vertex?
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < casccospa) continue;
      if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < dcaV0topv) continue;
      if (casc.dcaV0daughters() > dcav0dau) continue;
      if (casc.dcacascdaughters() > dcaCascdau) continue;
      if (TMath::Abs(casc.mLambda() - constants::physics::MassLambda) > massLambdaLimit) continue;
      //      if (TMath::Abs(casc.mXi() - constants::physics::MassXi) < ximassWindow) continue;
      if (TMath::Abs(casc.mXi() - Massxi) < ximassWindow) continue;
      //      if (TMath::Abs(casc.mOmega() - constants::physics::MassOmega) < omegarej) continue;
      if (properlifetimexi > properlifetimefactor*ctauxi) continue;

      if (TMath::Abs(casc.yXi()) > rapidity ) continue;
      if (TMath::Abs(casc.eta()) > eta ) continue;

      histos.fill(HIST("hMassXiAfterSel"), casc.mXi());
      histos.fill(HIST("hMassXiAfterSel2D"), casc.mXi(), casc.pt());
      xicounter++;
    } //end loop casc

    if (xicounter > 1) keepEvent[0] = true;

    if (xicounter>0){
      for (auto track : tracks) { // start loop over tracks                                                  
	histos.fill(HIST("hTriggeredParticles"), 1);
	if (track.itsChi2NCl()>4) continue; //some decrease observed                                         
	histos.fill(HIST("hTriggeredParticles"), 2);
	if (track.tpcNClsCrossedRows() < 80) continue; //~no changes                                          
	histos.fill(HIST("hTriggeredParticles"), 3);
	if (track.tpcCrossedRowsOverFindableCls() < 0.8) continue; //no changes                               
	histos.fill(HIST("hTriggeredParticles"), 4);
	if (track.length() < 90) continue; //some decrease observed                                           
	histos.fill(HIST("hTriggeredParticles"), 5);
	if (float(track.tpcNClsCrossedRows())/track.length() < 0.8) continue;
	histos.fill(HIST("hTriggeredParticles"), 6);
	histos.fill(HIST("PtTrigger"), track.pt());

	keepEvent[1] = true;                                                  
      } // end loop over tracks                    
    }

    if (keepEvent[0])    hProcessedEvents->Fill(1.5);
    if (keepEvent[1])    hProcessedEvents->Fill(2.5);
    //Filling the table
    strgtable(keepEvent);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilter>(cfgc, TaskName{"strangeness-filter"})};
}
