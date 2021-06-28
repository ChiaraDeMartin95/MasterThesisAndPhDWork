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
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

struct MyAwesomeTask { // This is a task
  // The histogram registry is the container of the output histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {

    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Histogram is added to the ouput registry
    histos.add("VtxZBeforeSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});
    histos.add("VtxZAfterSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});
    histos.add("p", "Momentum distribution;#it{p} (GeV/#it{c})", kTH1F, {{100, 0, 20}});
    histos.add("pt", "Transverse momentum distribution", kTH1F, {ptAxis});
  }

  // Equivalent of the AliRoot task UserExec
  void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  {
    // Performing the event selection
    histos.fill(HIST("VtxZBeforeSel"), coll.posZ());
    if (abs(coll.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("VtxZAfterSel"), coll.posZ());

    for (auto track : inputTracks) { // Loop over tracks
      histos.fill(HIST("p"), track.p());
      histos.fill(HIST("pt"), track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) // This puts your task in the DPL workflow
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<MyAwesomeTask>(cfgc)};
  return workflow;
}
