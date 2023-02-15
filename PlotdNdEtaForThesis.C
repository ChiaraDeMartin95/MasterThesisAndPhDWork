#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TSpline.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include </data/dataalice/cdemart/ALICE_analysis_tutorial/Macros/ErrRatioCorr.C>
#include "Macros/constants.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
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

void SetFont(TH1F *histo){
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1F *histo, Float_t TickLengthX, Float_t TickLengthY){
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

void SetHistoTextSize(TH1F *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize,  Float_t YOffset, Float_t YLabelOffset){
  histo->GetXaxis()->SetTitleSize(XSize);
  histo->GetXaxis()->SetLabelSize(XLabelSize);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetXaxis()->SetLabelOffset(XLabelOffset);
  histo->GetYaxis()->SetTitleSize(YSize);
  histo->GetYaxis()->SetLabelSize(YLabelSize);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetYaxis()->SetLabelOffset(YLabelOffset);
}

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin){
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}
void StyleTGraphErrors(TGraphAsymmErrors *tgraph, Int_t color, Int_t style, Float_t mSize, Int_t linestyle){
  tgraph->SetLineColor(color);
  tgraph->SetLineWidth(3);
  tgraph->SetMarkerColor(color);
  tgraph->SetMarkerStyle(style);
  tgraph->SetMarkerSize(mSize);
  tgraph->SetLineStyle(linestyle);
}

//take spectra in input (MB + HM 13 TeV)
//produces ratio of spectra wrt 0-100% multiplciity class

//takes Xi/K0s ratios in input (MB + HM 13 TeV) 
//produces doubles ratios wrt 0-100% multiplciity class

Int_t nummoltMax = 5;
const Int_t numtipo=10;
const Int_t numregions=3;

TString TitleYPtRatio = "N_{#Xi}/N_{K^{0}_{S}}";
TString TitleYPtDRatio = "N_{#Xi}/N_{K^{0}_{S}}_{mult} / N_{#Xi}/N_{K^{0}_{S}}_{0-100%}_{0-100%}";
TString titleX=  "#it{p}_{T} (GeV/#it{c})";
TString titleY=  "1/(#Delta#it{#eta} #Delta#it{#varphi}) 1/#it{N}_{trigg} d#it{N}/d#it{p}_{T} [(GeV/#it{c})^{-1}]";
TString title = "Multiplicity class "; 
TString titleYSpectraRatio = "(d#it{N}^{#Xi}/d#it{p}_{T}) / (d#it{N}^{K^{0}_{S}}/d#it{p}_{T})";
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};
TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString SRegionTypeBis[3] = {"Jet", "OOj", "Full"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};

TString sRegion1K0s[3]={"|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 1.1", "0.86 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1Xi[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.1", "0.75 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1K0sGen[3]={"|#Delta#it{#eta}| < 1.1, |#Delta#it{#varphi}| < 1.09", "0.86 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1XiGen[3]={"|#Delta#it{#eta}| < 1.1, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1[3][2] ={""};

Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
TString tipo[numtipo]={"K0s", "Xi"};
TString Stipo[numtipo]={"K^{0}_{S}", "#Xi"};
TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
TString Srap[2] = {"_Eta0.8", "_y0.5"};
TString SSkipAssoc[2] = {"_AllAssoc", ""};

TString titledNdetaTrigg="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5, #it{p}_{T,trigg}>3 GeV/#it{c}}";
TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";

void PlotdNdEtaForThesis(){

  gStyle->SetOptStat(0);

  TCanvas * canvas = new TCanvas("canvas", "canvas", 1100, 900); //1100, 900
  canvas->SetFillColor(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  StyleCanvas(canvas, 0.02, 0.15, 0.1, 0.02);

  Float_t UpperValueX = 40;

  Float_t xTitle = 35;
  //  Float_t xOffset = 4.3;
  Float_t xOffset = 1.5;

  Float_t yTitle = 35;
  Float_t yOffset = 1.2; 

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.02;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.02;
  Float_t tickY = 0.02;
  Float_t tickXL = 0.03;
  Float_t tickYL = 0.03;

  Float_t DummyValue[8] = {0};
  Float_t DummyValueBis[8] = {0};
  Float_t DummyValueError[8] = {0};
  TGraphAsymmErrors *ghistoEvts;
  TGraphAsymmErrors *ghistoEvtsWTrig;
  TGraphAsymmErrors *ghistoEvtsMult[8];
  TGraphAsymmErrors *ghistoEvtsWTrigMult[8];
  Int_t Color[2] = {kOrange+2, kRed+2};
  //  Int_t ColorMult[8] = {634, 628, 797,815,418, 429, 867,601};
  Int_t ColorMult[8] = {601, 867, 429, 418, 815, 797, 628, 634};
  Int_t LineStyle[2] = {1, 1};
  //  TString SmoltLegend[8]={"50-100%", "30-50%", "10-30%", "5-10%", "0-5%", "0.05-0.1%", "0.01-0.05%", "0-0.01%"};
  TString SmoltLegend[8]={"#color[601]{50-100%}", "#color[867]{30-50%}", "#color[429]{10-30%}", "#color[418]{5-10%}", "#color[815]{0-5%}", "#color[797]{0.05-0.1%}", "#color[628]{0.01-0.05%}", "#color[634]{0-0.01%}"};

  TLegend *legendMultPaper;
  legendMultPaper  = new TLegend(0.15, 0.72, 0.87, 0.93);
  legendMultPaper->SetHeader("V0M Multiplicity Percentile");
  legendMultPaper->SetNColumns(2);
  legendMultPaper->SetNColumns(3);
  legendMultPaper->SetFillStyle(0);
  legendMultPaper->SetTextSize(0.04);
  TLegendEntry *lheaderMultPaper = (TLegendEntry*)legendMultPaper->GetListOfPrimitives()->First();
  lheaderMultPaper-> SetTextSize(0.05);

  TLegend * legend= new  TLegend(0.15, 0.22, 0.5, 0.33);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.04);


  TH1F* hdummy = new TH1F ("hdummy", "hdummy",  1000, 0, UpperValueX);
  StyleHisto(hdummy, 0.2, 2, 1, 1, titledNdeta, "Arbitrary value", "" , 1, 0, UpperValueX, xOffset, yOffset, 1);
  SetFont(hdummy);
  SetHistoTextSize(hdummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hdummy, tickXL, tickYL);
  hdummy->GetYaxis()->SetDecimals(kTRUE);
  //  hdummy->GetYaxis()->SetNdivisions(405, "false");

  Double_t  dNdEtaTabulated13TeVMult[1] ={0};
  Double_t  dNdEtaFinal13TeVMult[1] ={0};
  Double_t  dNdEtaTabulated13TeVMult_ErrorL[1] ={0};
  Double_t  dNdEtaFinal13TeVMult_ErrorL[1] ={0};
  Double_t  dNdEtaTabulated13TeVMult_ErrorR[1] ={0};
  Double_t  dNdEtaFinal13TeVMult_ErrorR[1] ={0};
  Double_t DummyValueMult[1] = {0};
  Double_t DummyValueMultBis[1] = {0};
  Double_t DummyValueErrorMult[1] = {0};

  for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
    if (Coll==1) continue; //for the time being, only 13 TeV
    for (Int_t m=0; m< 8; m++){
      DummyValue[m] = 1;
      DummyValueError[m] = 0.1;

      DummyValueMult[0] = 1.15;
      DummyValueMultBis[0] = 0.85;
      DummyValueErrorMult[0] = 0.1;
      dNdEtaTabulated13TeVMult[0] = dNdEtaTabulated13TeV[m];
      dNdEtaFinal13TeVMult[0] = dNdEtaFinal13TeV[m];
      dNdEtaTabulated13TeVMult_ErrorL[0] = dNdEtaTabulated13TeV_ErrorL[m];
      dNdEtaFinal13TeVMult_ErrorL[0] = dNdEtaFinal13TeV_ErrorL[m];
      dNdEtaTabulated13TeVMult_ErrorR[0] = dNdEtaTabulated13TeV_ErrorR[m];
      dNdEtaFinal13TeVMult_ErrorR[0] = dNdEtaFinal13TeV_ErrorR[m];

      cout << dNdEtaTabulated13TeV[m] << " + " << dNdEtaTabulated13TeV_ErrorR[m] << " - " << dNdEtaTabulated13TeV_ErrorL[m] << endl;

      ghistoEvtsMult[m] = new TGraphAsymmErrors(1, dNdEtaTabulated13TeVMult,DummyValueMult,dNdEtaTabulated13TeVMult_ErrorL,dNdEtaTabulated13TeVMult_ErrorR,DummyValueErrorMult,DummyValueErrorMult);
      ghistoEvtsMult[m]->SetName(Form("ghistoEvts_m%i", m));
      StyleTGraphErrors(ghistoEvtsMult[m], ColorMult[m], 24, 2, LineStyle[0]);
      ghistoEvtsMult[m]->SetFillStyle(0);
      ghistoEvtsMult[m]->SetFillColor(ColorMult[m]);
      ghistoEvtsMult[m]->SetLineWidth(0);
      ghistoEvtsWTrigMult[m] = new TGraphAsymmErrors(1, dNdEtaFinal13TeVMult,DummyValueMultBis,dNdEtaFinal13TeVMult_ErrorL,dNdEtaFinal13TeVMult_ErrorR,DummyValueErrorMult,DummyValueErrorMult);
      ghistoEvtsWTrigMult[m]->SetName(Form("ghistoEvtsWTrig_m%i", m));
      StyleTGraphErrors(ghistoEvtsWTrigMult[m], ColorMult[m], 33, 3, LineStyle[1]);
      ghistoEvtsWTrigMult[m]->SetFillStyle(0);
      ghistoEvtsWTrigMult[m]->SetFillColor(ColorMult[m]);
      ghistoEvtsWTrigMult[m]->SetLineWidth(0);

      legendMultPaper->AddEntry(ghistoEvtsWTrigMult[m], SmoltLegend[m], "p");

    }
    ghistoEvts = new TGraphAsymmErrors(8, dNdEtaTabulated13TeV,DummyValue,dNdEtaTabulated13TeV_ErrorL,dNdEtaTabulated13TeV_ErrorR,DummyValueError,DummyValueError);
    ghistoEvts->SetName("ghistoEvts");
    StyleTGraphErrors(ghistoEvts, Color[0], 33, 2, LineStyle[0]);
    ghistoEvts->SetFillStyle(0);
    ghistoEvts->SetFillColor(Color[0]);
    ghistoEvts->SetLineWidth(0);
    ghistoEvtsWTrig = new TGraphAsymmErrors(8, dNdEtaFinal13TeV,DummyValue,dNdEtaFinal13TeV_ErrorL,dNdEtaFinal13TeV_ErrorR,DummyValueError,DummyValueError);
    ghistoEvtsWTrig->SetName("ghistoEvtsWTrig");
    StyleTGraphErrors(ghistoEvtsWTrig, Color[1], 33, 2, LineStyle[1]);
    ghistoEvtsWTrig->SetFillStyle(0);
    ghistoEvtsWTrig->SetFillColor(Color[1]);
    ghistoEvtsWTrig->SetLineWidth(0);

  }

  hdummy->Draw("");
  //ghistoEvts->Draw("same p2");
  //ghistoEvtsWTrig->Draw("same p2");
  for (Int_t m=0; m< 8; m++){
    ghistoEvtsMult[m]->DrawClone("same p2");
    ghistoEvtsWTrigMult[m]->DrawClone("same p2");
  }
  legendMultPaper->Draw("");
  ghistoEvtsMult[0]->SetMarkerColor(kGray+2);
  ghistoEvtsWTrigMult[0]->SetMarkerColor(kGray+2);
  legend->AddEntry(ghistoEvtsMult[0], "All events", "p");
  legend->AddEntry(ghistoEvtsWTrigMult[0], "Events with #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "p");
  legend->Draw("");
  canvas->SaveAs("PlotdNdEta_TabulatedVsRecomputed.pdf");
}

