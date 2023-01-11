#include <Riostream.h>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TCutG.h>
#include "TFitResult.h"
#include "TLegend.h"

void StyleCanvas(TCanvas *canvas, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  // gStyle->SetPalette(55, 0);
}

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX,
                TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
}

Bool_t reject;
Double_t fparab(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[3] == 0)
  {
    LimInf = 0.474;
    LimSup = 0.520;
  }
  else if (par[3] == 4 || par[3] == 5 || par[3] == 8)
  {
    LimInf = 1.310;
    LimSup = 1.335;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

Double_t fretta(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[2] == 0)
  {
    LimInf = 0.47;
    LimSup = 0.530;
  }
  else if (par[2] == 4 || par[2] == 5 || par[2] == 8)
  {
    LimInf = 1.310;
    LimSup = 1.335;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0];
}

const Int_t numPart = 7;
TString TitleInvMass[numPart] = {"(#pi^{+}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(p, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#bar{p}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#Lambda, #pi^{-}) invariant mass (GeV/#it{c}^{2})"};
TString namehisto[numPart] = {"h3dMassK0Short", "", "", "h2dMassXiMinus", "h2dMassXiPlus", "h2dMassOmegaMinus", "h2dMassOmegaPlus"};
Float_t LowLimitMass[numPart] = {0.42, 1.09, 1.09, 1.29, 1.29, 1.62, 1.62}; // 0.44
Float_t UpLimitMass[numPart] = {0.57, 1.14, 1.14, 1.35, 1.35, 1.72, 1.72};  // 0.55
Float_t LowMassRange[numPart] = {0.48, 1.09, 1.09, 1.31};
Float_t UpMassRange[numPart] = {0.51, 1.14, 1.14, 1.33};

Float_t min_range_signal[numPart] = {0.46, 1.105, 1.105, 1.31, 1.31, 1.66, 1.66}; // estremi region fit segnale (gaussiane)
Float_t max_range_signal[numPart] = {0.535, 1.125, 1.125, 1.334, 1.334, 1.685, 1.685};
Float_t min_histo[numPart] = {0.42, 1.09, 1.09, 1.30, 1.30, 1.62, 1.62}; // estremi del range degli istogrammi
Float_t max_histo[numPart] = {0.57, 1.14, 1.14, 1.342, 1.342, 1.72, 1.72};
Float_t liminf[numPart] = {0.45, 1.1153, 1.1153, 1.30, 1.30, 1.66, 1.66}; // estremi regione fit del bkg e total
Float_t limsup[numPart] = {0.545, 1.1168, 1.1168, 1.342, 1.342, 1.685, 1.685};

Float_t lim_inf_mean[numPart] = {0.495, 1.1153, 1.1153, 1.31, 1.31, 1.66, 1.66};
Float_t lim_sup_mean[numPart] = {0.500, 1.1168, 1.1168, 1.33, 1.33, 1.685, 1.685};
Float_t lim_inf_sigma[numPart] = {0};
Float_t lim_sup_sigma[numPart] = {0.008, 0.002, 0.002, 0.008, 0.008, 0.008, 0.008};
Float_t lim_inf_errmean[numPart] = {0};
Float_t lim_sup_errmean[numPart] = {10, 0.0006, 0.0006, 10, 10, 10, 10}; // loooooose
Float_t lim_inf_errsigma[numPart] = {0};
Float_t lim_sup_errsigma[numPart] = {10, 0.0004, 0.0004, 10, 10, 10, 10}; // loose

const Float_t massParticle[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
TString Spart[numPart] = {"K0s", "Lambda", "AntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus"};

void PostProcessing_Filters(TString year = "LHC22m_pass2",
                            TString SPathIn = "../TriggerForRun3/AnalysisResults_22mpass2_New_5_NoGlobalTrack.root",
                            Int_t part = 0,
                            Bool_t UseTwoGauss = 1,
                            Bool_t isBkgParab = 0,
                            Bool_t isMeanFixedPDG = 0,
                            Float_t sigmacentral = 4)
{

  TFile *filein = new TFile(SPathIn, "");
  if (!filein)
  {
    cout << "FileIn not available" << endl;
    return;
  }
  TDirectoryFile *dir;
  dir = (TFile *)filein->Get("lf-strangeness-filter");
  if (!dir)
  {
    cout << "dir not available" << endl;
    return;
  }

  TH1F *hEvents = (TH1F *)dir->Get("hProcessedEvents");
  if (!hEvents)
  {
    cout << "hProcessedEvents not available " << endl;
    return;
  }

  Int_t NEvents = hEvents->GetBinContent(1);
  hEvents->Scale(1. / NEvents);
  hEvents->GetYaxis()->SetRangeUser(10e-8, 1);

  TCanvas *canvasNEvents = new TCanvas("canvasNEvents", "canvasNEvents", 800, 500);
  gPad->SetLogy();
  hEvents->Draw("");

  TH1F *hCascCandidates = (TH1F *)dir->Get("hCandidate");
  if (!hCascCandidates)
  {
    cout << "hCandidate not available " << endl;
    return;
  }

  Int_t NInitialCasc = hCascCandidates->GetBinContent(1);
  hCascCandidates->Scale(1. / NInitialCasc);
  TCanvas *canvasCasc = new TCanvas("canvasCasc", "canvasCasc", 800, 500);
  hCascCandidates->Draw("text");

  // Track QA
  TDirectoryFile *dirQATrack;
  dirQATrack = (TDirectoryFile *)dir->Get("QAHistosTriggerParticles");
  if (!dirQATrack)
  {
    cout << "Directory QAHistosTriggerParticles not available" << endl;
    return;
  }

  TCanvas *canvasTrackQA = new TCanvas("canvasTrackQA", "canvasTrackQA", 800, 500);
  canvasTrackQA->Divide(3, 2);

  const Int_t numQATrackHistos = 5;
  TH1F *hTrackQAPt = (TH1F *)dirQATrack->Get("hPtTriggerAllEv");
  canvasTrackQA->cd(1);
  hTrackQAPt->Draw("");

  TH2F *hTrackQA[numQATrackHistos - 1];
  TH1F *hTrackQA1D[numQATrackHistos - 1];
  TString hSTrackQA[numQATrackHistos - 1] = {"hEtaTriggerAllEv", "hPhiTriggerAllEv", "hDCAxyTriggerAllEv", "hDCAzTriggerAllEv"};

  const Int_t numPtTrigg = 7;
  Float_t binptTrigg[numPtTrigg + 1] = {3, 4, 5, 6, 7, 8, 10, 15};
  Int_t ColorPtTrigg[numPtTrigg + 1] = {634, 628, 797, 815, 418, 429, 867};
  for (Int_t i = 0; i < numQATrackHistos - 1; i++)
  {
    hTrackQA[i] = (TH2F *)dirQATrack->Get(hSTrackQA[i]);
    if (!hTrackQA[i])
    {
      cout << hSTrackQA[i] << " not available" << endl;
      return;
    }

    canvasTrackQA->cd(i + 2);
    for (Int_t pt = 0; pt < numPtTrigg; pt++)
    {
      if (binptTrigg[pt] > 6)
        continue;
      hTrackQA1D[i] = (TH1F *)hTrackQA[i]->ProjectionX(Form("hTrackQA_var%i_pt%i", i, pt), hTrackQA[i]->GetYaxis()->FindBin(binptTrigg[pt] + 0.001), hTrackQA[i]->GetYaxis()->FindBin(binptTrigg[pt + 1] - 0.001));
      hTrackQA1D[i]->SetLineColor(ColorPtTrigg[pt]);
      hTrackQA1D[i]->SetMarkerColor(ColorPtTrigg[pt]);
      hTrackQA1D[i]->Scale(1. / hTrackQA1D[i]->GetEntries());

      if (hSTrackQA[i] == "hPhiTriggerAllEv")
        hTrackQA1D[i]->Rebin(2);
      hTrackQA1D[i]->Rebin(2);
      if (hSTrackQA[i] == "hDCAxyTriggerAllEv")
        hTrackQA1D[i]->GetXaxis()->SetRangeUser(-0.1, 0.1);
      hTrackQA1D[i]->Draw("same");
    }
  }

  TH1F *hTriggerParticles = (TH1F *)dirQATrack->Get("hTriggeredParticlesAllEv");
  if (!hTriggerParticles)
  {
    return;
  }
  TCanvas *canvasTrigger = new TCanvas("canvasTrigger", "canvasTrigger", 800, 500);
  hTriggerParticles->Draw("");

  // Cascade topological variables (after all topo sleections, for the time being)
  TDirectoryFile *dirCascVar;
  dirCascVar = (TDirectoryFile *)dir->Get("QAHistosTopologicalVariables");
  if (!dirCascVar)
  {
    cout << "Directory QAHistosTopologicalVariables not available" << endl;
    return;
  }

  const int NTopCascVar = 11;

  TH1F *hTopVar[NTopCascVar];
  TCanvas *canvasTopology[2];
  canvasTopology[0] = new TCanvas("canvasTopology1", "canvasTopology1", 1800, 1400);
  canvasTopology[0]->Divide(3, 2);
  canvasTopology[1] = new TCanvas("canvasTopology2", "canvasTopology2", 1800, 1400);
  canvasTopology[1]->Divide(3, 2);

  TString TopVarCascInput[NTopCascVar] = {"CascCosPA", "V0CosPA", "CascRadius", "V0Radius",
                                          "InvMassLambda", "DCACascDaughters", "DCAV0Daughters", "DCABachToPV",
                                          "DCAV0ToPV", "DCAPosToPV", "DCANegToPV"};
  TString TopVarCasc[NTopCascVar] = {"Casc #it{cos}#theta_{PA}", "V0 #it{cos}#theta_{PA}", "Casc #it{R}", "V0 #it{R}",
                                     "#it{m}_{inv} #Lambda Daughter", "DCA Casc Daughters", "DCA V0 Daughters", "DCA Bach. To PV",
                                     "DCA V0 To PV", "DCA Pos. To PV", "DCA Neg. To PV"};
  TString TopVarCascUnit[NTopCascVar] = {"", "", "(cm)", "(cm)",
                                         "(GeV/#it{c}^2)", "(cm)", "(cm)", "(cm)",
                                         "(cm)", "(cm)", "(cm)"};

  for (Int_t var = 0; var < NTopCascVar; var++)
  {
    hTopVar[var] = (TH1F *)dirCascVar->Get(TopVarCascInput[var]);
    if (!hTopVar[var])
    {
      cout << TopVarCascInput[var] << " not available" << endl;
      return;
    }
    hTopVar[var]->Scale(1. / NEvents);
    hTopVar[var]->GetYaxis()->SetRangeUser(0.1 * hTopVar[var]->GetMinimum(1.e-10), 10 * hTopVar[var]->GetMaximum());
    hTopVar[var]->SetTitle(TopVarCasc[var]);
    hTopVar[var]->GetXaxis()->SetTitle(TopVarCasc[var] + " " + TopVarCascUnit[var]);
    hTopVar[var]->GetYaxis()->SetTitle("1/N_{ev} Counts");
    if (var < 6)
      canvasTopology[0]->cd(var + 1);
    else
      canvasTopology[1]->cd(var + 1 - 6);
    hTopVar[var]->DrawCopy("hist");
  }

  // inv mass distributions of Xi (all selected xis) and Omegas (all selected omegas) vs pT
  TDirectoryFile *dirCasc;
  dirCasc = (TDirectoryFile *)dir->Get("QAHistos");
  if (!dirCasc)
  {
    cout << "Directory QAHistos not available" << endl;
    return;
  }
  TH2F *hMassXivsPt = (TH2F *)dirCasc->Get("hMassXiAfterSelvsPt");
  if (!hMassXivsPt)
  {
    cout << "hMassXiAfterSelvsPt not available" << endl;
    return;
  }

  const Int_t numPt = 6; // six pt intervals
  // Xi
  Float_t binpt[numPt + 1] = {1.0, 1.5, 1.8, 2.1, 2.4, 2.6, 6.0};

  // Omega
  // Float_t binpt[numPt + 1] = {1.0, 1.5, 2.0, 3.0, 6.0};

  TString SPt[numPt] = {""};
  TH1F *hInvMass[numPt];

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1800, 1400);
  canvas->Divide(numPt / 2, 2);
  StyleCanvas(canvas, 0.15, 0.05, 0.05, 0.15);

  TH1F *histoCountsPerEvent = new TH1F("histoCountsPerEvent", "histoCountsPerEvent", numPt, binpt);
  TH1F *histoYield = new TH1F("histoYield", "histoYield", numPt, binpt);

  Float_t counts = 0;
  Float_t errcount = 0;

  for (Int_t pt = 0; pt < numPt; pt++)
  {
    SPt[pt] = Form("%.1f < p_{T} < %.1f", binpt[pt], binpt[pt + 1]);
    // cout << binpt[pt] << endl;

    hInvMass[pt] = (TH1F *)hMassXivsPt->ProjectionY(Form("hInvMass_pt%i", pt), hMassXivsPt->GetXaxis()->FindBin(binpt[pt] + 0.001), hMassXivsPt->GetXaxis()->FindBin(binpt[pt + 1] - 0.001));
    if (part < 3)
      hInvMass[pt]->Rebin(4);
    else
      hInvMass[pt]->Rebin(2);
    StyleHisto(hInvMass[pt], 0, 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), 1, 20, TitleInvMass[part], "Counts", SPt[pt], 1, LowLimitMass[part], UpLimitMass[part], 1.4, 1.4, 1.2);
    canvas->cd(pt + 1);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.2);
    hInvMass[pt]->Draw("e same");

    counts = 0;
    errcount = 0;
    for (Int_t bmass = hInvMass[pt]->GetXaxis()->FindBin(LowMassRange[part]); bmass <= hInvMass[pt]->GetXaxis()->FindBin(UpMassRange[part]); bmass++)
    {
      counts += hInvMass[pt]->GetBinContent(bmass);
    }
    errcount = sqrt(counts);
    histoCountsPerEvent->SetBinContent(pt + 1, counts / NEvents / histoCountsPerEvent->GetBinWidth(pt + 1));
    histoCountsPerEvent->SetBinError(pt + 1, errcount / NEvents / histoCountsPerEvent->GetBinWidth(pt + 1));
  }

  // h-Xi events vs pt,thr

  TString Soutputfile = "../TriggerForRun3/" + year;
  canvasNEvents->SaveAs(Soutputfile + ".pdf(");
  canvasTrackQA->SaveAs(Soutputfile + ".pdf");
  canvasTrigger->SaveAs(Soutputfile + ".pdf");
  canvasCasc->SaveAs(Soutputfile + ".pdf");
  canvasTopology[0]->SaveAs(Soutputfile + ".pdf");
  canvasTopology[1]->SaveAs(Soutputfile + ".pdf");
  canvas->SaveAs(Soutputfile + ".pdf)");

  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  outputfile->WriteTObject(hEvents);
  outputfile->WriteTObject(hCascCandidates);
  outputfile->WriteTObject(hTriggerParticles);
  outputfile->WriteTObject(canvasTrackQA);
  outputfile->Close();
  cout << "Ho creato il file: " << Soutputfile << " (.pdf and .root)" << endl;

  cout << "\nTotal number of processed events " << NEvents << endl;
}
