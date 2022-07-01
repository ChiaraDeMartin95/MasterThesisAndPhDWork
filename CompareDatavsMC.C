#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include <TSpline.h>
#include <TGraphErrors.h>
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
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include <Macros/constants.h>

TSpline3 *sp3;
Double_t spline(Double_t *x, Double_t* p) {
  Double_t xx = x[0];
  return sp3->Eval(xx);
}

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

void StyleTGraphErrors(TGraphErrors *tgraph, Int_t color, Int_t style, Float_t mSize, Int_t linestyle){
  tgraph->SetLineColor(color);
  tgraph->SetLineWidth(2);
  tgraph->SetMarkerColor(color);
  tgraph->SetMarkerStyle(style);
  tgraph->SetMarkerSize(mSize);
  tgraph->SetLineStyle(linestyle);
}

TString TitleYieldRatio="#Xi/K^{0}_{S} yield ratio vs multiplicity";
TString titleRelUnc = "Relative uncertainty";
TString titleMult = "Multiplicity class ";
TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
TString titleYield[2]={"K^{0}_{S}", "#Xi"};
TString TitleYYieldRatio="#it{N}_{#Xi}/#it{N}_{K^{0}_{S}}";

TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})", "#it{N}_{#Xi}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})"};
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString RegionTypeBis[3] = {"Jet", "OOj", "Full"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
//TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
//TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
TString sRegion[3]={"#color[634]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
//TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};
//TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
TString sRegion1[3]={"|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 0.85", "0.86 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 2.0", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
TString titleYToOOJ = "Toward / Transverse";
TString titleYMCToData = "Model / Data";

const Int_t NumberOfMCs=2;
const Int_t numColls =2; //pp13 TeV (MB + HM), pp5TeV
const Int_t numRegions =3;
const Int_t numDataorMC = NumberOfMCs+1;
Int_t nummoltMax = nummolt;
Int_t nummoltMaxDataorMC = 0;

void CompareDatavsMC( Int_t PlotType =0, Int_t ChosenRegion = -1,  Float_t ScalingFactorXiK0s = 1/*0.8458/1.08747*/, Bool_t isChangesIncluded=1, Bool_t isFit=0, Bool_t isdNdEtaTriggered=0){

  Int_t nummoltDataorMC =0;

  Int_t MonashTune =0;
  cout << "Should be analyse Ropes (=1) or MonashDefault (=2) or both? (=3)? " << endl;
  cin >> MonashTune;

  Bool_t isWingsCorrectionApplied=0;
  cout <<"Do you want to analyse the files with the wings correction applied FOR THE MC? Type 1 if you DO want" << endl;
  cin >> isWingsCorrectionApplied;

  //PlotType = 0: Xi/K0s ratio
  //PlotType = 1: K0s yield vs mult 
  //PlotType = 2: Xi yield vs mult 
  //PlotType = 3: K0s pt vs mult 
  //PlotType = 4: Xi pt vs mult 

  if (PlotType==1 || PlotType ==3) {
    //    sRegion1[0] = "|#Delta#it{#eta}| < 0.85, |#Delta#it{#varphi}| < 1.09";
    //    sRegion1[1] = "0.85 < |#Delta#it{#eta}| < 1.2, 1.09 < #Delta#it{#varphi} <1.8";
  }

  //Set titles
  TString titleY;
  if (PlotType==0) titleY = TitleYYieldRatio;
  else if (PlotType==1) titleY = titleYieldYType[0];
  else if (PlotType==2) titleY = titleYieldYType[1];
  else if (PlotType==3) titleY = titlePtvsMultYType[0];
  else if (PlotType==4) titleY = titlePtvsMultYType[1];

  gStyle->SetOptStat(0);

  //STYLE 
  TString titleX = titledNdeta;
  Float_t xOffset =1.2;
  Float_t yOffset =1.25;
  Float_t MarkerSize[3] ={2, 2, 3};
  Int_t MarkerType[3] = {20, 21, 33};
  Int_t Color[numRegions] = {628,418,601};
  Int_t ColorEnergy[numColls] = {922,920};
  Int_t ColorDiff[numRegions][numColls] = {{634,628}, {418,829} , {601, 867}};
  Int_t ColorDiffMC[numRegions][numColls] = {{634,628}, {418,829} , {601, 867}};
  Int_t ColorDiffMCRopes[numRegions][numColls] = {{634,628}, {418,829} , {601, 867}};
  //  Int_t ColorDiffMC[numRegions][numColls] = {{907, 907}, {812,812} , {870, 870}};
  //  Int_t ColorDiffMCRopes[numRegions][numColls] = {{807, 907}, {411,812} , {598, 870}};
  Int_t LineStyle[numDataorMC] = {1, 8, 1};

  //ChosenRegion == -1 //all regions (JET, OOJ, FULL) are plotted
  if (ChosenRegion>2) {cout << "Choose a valid region! " << endl; return;}

  TString CollisionsComp[2] = {"_HMMultBinning1_vsHM", "_vs5TeV"};
  TString SPlotType[5] = {"Ratio", "K0s", "Xi", "K0spt", "Xipt"};
  TString SPlotTypeBis[5] = {"", "K0s", "Xi", "K0s", "Xi"};
  TString SisFit[2] = {"", "_isFit"};
  TString Region[numRegions] = {"Jet", "Bulk", "All"};
  TString RegionBis[numRegions] = {"Jet", "Bulk", "Inclusive"};

  Float_t Up= 0.14-10e-6;
  Float_t Low = 0.02+10e-6;
  if (PlotType==0){
    //    Low = 0.014 + 10e-6;
    Low = 0.012 + 10e-6;
    //    Up = 0.12-10e-6;
    Up = 0.5-10e-6;
  }
  if (PlotType==1) {
    Low = 10e-6; Up = 0.50-10e-6; //0.45
    if (ChosenRegion==0) {Low = 0.015+10e-6; Up = 0.08-10e-6;} //0.035 
  }
  if (PlotType==2) {
    Low = 10e-6; Up = 0.04-10e-6;
    if (ChosenRegion==0) {Up = 0.0045-10e-6;}
  }
  if (PlotType==3) {
    Low = 0+10e-4; Up = 5.-10e-4;
    if (ChosenRegion==0) {Low = 1.2+10e-4; Up = 2.7-10e-4;}
    else if (ChosenRegion>0) {Low = 0.8+10e-4; Up = 1.3-10e-4;}
  }
  if (PlotType==4) {
    Low = 0+10e-4; Up = 6.999;
    if (ChosenRegion==0) {Low = 1.5; Up = 5;}
    else if (ChosenRegion>0) {Low = 1.0; Up = 2.2;}
  }

  Float_t LowRatio = 0.5+10e-4;
  Float_t UpRatio = 2.5-10e-4;
  if (PlotType==0) {
    LowRatio = 0.+10e-4;
    UpRatio = 1.5-10e-4;
  }
  else if (PlotType==1) {
    LowRatio = 0.5+10e-4;
    UpRatio = 2.5-10e-4;
  }
  else if (PlotType==2) {
    LowRatio = 0.05+10e-4;
    UpRatio = 2.5-10e-4;
  }
  else if (PlotType==3){
    LowRatio = 0.7+10e-4;
    UpRatio = 1.3-10e-4;
  }
  else if (PlotType ==4) {
    LowRatio = 0.7+10e-4;
    UpRatio = 1.3-10e-4;
  }

  TF1 * fitToRatioToOOJ = new TF1("fitToRatioToOOJ", "pol0", 0, 45);
  fitToRatioToOOJ->SetLineColor(1);
  fitToRatioToOOJ->SetLineStyle(2);
  fitToRatioToOOJ->SetLineWidth(2);

  TH1F *histoYield[numRegions][numColls][numDataorMC];
  TH1F *histoYieldSist[numRegions][numColls][numDataorMC];
  TH1F *histoYieldRatio[numRegions][numColls][numDataorMC];
  TH1F *hMCDataRatio[numRegions][numColls][numDataorMC];
  TSpline3 * splineY[numRegions][numColls][numDataorMC];
  TF1* fsplineY[numRegions][numColls][numDataorMC];
  TF1* fitY[numRegions][numColls][numDataorMC];
  TF1* fitRatios[numRegions][numColls][numDataorMC];

  Float_t Yields[numRegions][numColls][numDataorMC][nummoltMax] = {0};
  Float_t YieldsErrors[numRegions][numColls][numDataorMC][nummoltMax] = {0};
  Float_t YieldRatios[numRegions][numColls][numDataorMC][nummoltMax] = {0};
  Float_t YieldRatiosErrors[numRegions][numColls][numDataorMC][nummoltMax] = {0};
  TGraphErrors*	ghistoYield[numRegions][numColls][numDataorMC];
  TGraphErrors*	ghistoYieldGrey[numRegions][numColls][numDataorMC];
  TGraphErrors*	ghistoYieldRatio[numRegions][numColls][numDataorMC];

  TH1F*  fHistYieldStatBlack;
  TH1F*  fHistYieldSistBlack;
  TH1F*  fHistYieldStatGrey[numColls];

  TString NameHisto13TeV="histoYieldComparison";
  TString NameHistoSist13TeV="histoYieldSistComparison";
  TString NameHisto5TeV="fHistYield_pp5TeV";
  TString NameHistoSist5TeV="fHistYield_pp5TeV_Sist";
  TString NameHistoFinal[numRegions][numColls][numDataorMC];

  TString MCType[NumberOfMCs+1]= {"Data", "Monash default", "Monash ropes"};
  TString MCTypeBis[3]= {"PythiaRopes", "PythiaMonash", "PythiaAll"};
  
  Float_t dNdEta[numDataorMC][nummoltMax] = {0};
  Float_t dNdEtaError[numDataorMC][nummoltMax] = {0};
  /* MonashRopes, events with trigger particle */
  if (isdNdEtaTriggered){
    //Monash
    dNdEta[2][0] = 8.45;
    dNdEta[2][1] = 12.42;
    dNdEta[2][2] = 15.11;
    dNdEta[2][3] = 17.15;
    dNdEta[2][4] = 19.12;
    dNdEta[2][5] = 20.99;
    dNdEta[2][6] = 22.75;
    dNdEta[2][7] = 24.84;
    dNdEta[2][8] = 27.36;
    dNdEta[2][9] = 30.35;
    dNdEta[2][10] = 16.5;
    //Ropes
    dNdEta[1][0] = 8.27;
    dNdEta[1][1] = 12.48;
    dNdEta[1][2] = 15.23;
    dNdEta[1][3] = 17.26;
    dNdEta[1][4] = 19.22;
    dNdEta[1][5] = 21.1;
    dNdEta[1][6] = 22.9;
    dNdEta[1][7] = 25.1;
    dNdEta[1][8] = 27.81;
    dNdEta[1][8] = 27.81;
    dNdEta[1][9] = 31.26;
    dNdEta[1][10] = 16.78;
  }
  else {    
    //DATA
    dNdEta[0][0] = 3.33;
    dNdEta[0][1] = 7.14;
    dNdEta[0][2] = 11.46;
    dNdEta[0][3] = 16.17;
    dNdEta[0][4] = 21.2;
    dNdEta[0][5] = 30.43;
    dNdEta[0][6] = 32.57;
    dNdEta[0][7] = 36.29;
    dNdEta[0][8] = 0;
    dNdEta[0][9] = 0;
    dNdEta[0][10] = 0;   

    dNdEtaError[0][0] = 0;
    dNdEtaError[0][1] = 0;
    dNdEtaError[0][2] = 0;
    dNdEtaError[0][3] = 0;
    dNdEtaError[0][4] = 0;
    dNdEtaError[0][5] = 0;
    dNdEtaError[0][6] = 0;
    dNdEtaError[0][7] = 0;
    dNdEtaError[0][8] = 0;
    dNdEtaError[0][9] = 0;
    dNdEtaError[0][10] = 0;

    //Monash
    dNdEta[1][0] = 4.01;
    dNdEta[1][1] = 9.04;
    dNdEta[1][2] = 12.1;
    dNdEta[1][3] = 14.29;
    dNdEta[1][4] = 16.4;
    dNdEta[1][5] = 18.46;
    dNdEta[1][6] = 20.37;
    dNdEta[1][7] = 22.68;
    dNdEta[1][8] = 25.58;
    dNdEta[1][9] = 29.09;
    dNdEta[1][10] = 7.00;   

    dNdEtaError[1][0] = 0;
    dNdEtaError[1][1] = 0;
    dNdEtaError[1][2] = 0;
    dNdEtaError[1][3] = 0;
    dNdEtaError[1][4] = 0;
    dNdEtaError[1][5] = 0;
    dNdEtaError[1][6] = 0;
    dNdEtaError[1][7] = 0;
    dNdEtaError[1][8] = 0;
    dNdEtaError[1][9] = 0;
    dNdEtaError[1][10] = 0;

    //Ropes (now it's Monash...)
    dNdEta[2][0] = 4.01;
    dNdEta[2][1] = 9.04;
    dNdEta[2][2] = 12.1;
    dNdEta[2][3] = 14.29;
    dNdEta[2][4] = 16.4;
    dNdEta[2][5] = 18.46;
    dNdEta[2][6] = 20.37;
    dNdEta[2][7] = 22.68;
    dNdEta[2][8] = 25.58;
    dNdEta[2][9] = 29.09;
    dNdEta[2][10] = 7.00;   

    dNdEtaError[2][0] = 0;
    dNdEtaError[2][1] = 0;
    dNdEtaError[2][2] = 0;
    dNdEtaError[2][3] = 0;
    dNdEtaError[2][4] = 0;
    dNdEtaError[2][5] = 0;
    dNdEtaError[2][6] = 0;
    dNdEtaError[2][7] = 0;
    dNdEtaError[2][8] = 0;
    dNdEtaError[2][9] = 0;
    dNdEtaError[2][10] = 0;
  }

  TString pathin[NumberOfMCs+2] ={""};
  TString YieldOrAvgPt[5] = {"Yield", "Yield", "Yield", "AvgPt","AvgPt"};

  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  canvas->SetFillColor(0);
  TPad* pad1 = new TPad( "pad1" ,"pad1" ,0 ,0.36 ,1 ,1);
  TPad* pad2 = new TPad( "pad2" ,"pad2" ,0 ,0.01 ,1 ,0.35);

  if (PlotType==0)  pad1->SetLogy();
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.005);
  pad1->SetLeftMargin(0.13);
  pad1->SetRightMargin(0.05);

  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.32);
  pad2->SetLeftMargin(0.13);
  pad2->SetRightMargin(0.05);

  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad2->SetTickx(1);
  pad2->SetTicky(1);

  //legend
  TLegend *legendFit = new TLegend(0.62,0.72,0.9,0.92);
  TLegend *legendFitRatioToOOJ = new TLegend(0.7,0.8,0.8,0.96);
  TLegendEntry * legendFitRatio1;

  //  TLegend *LegendRatio=new TLegend(0.62,0.72,0.9,0.92);
  TLegend *LegendRatio=new TLegend(0.62,0.77,0.9,0.92);
  LegendRatio->SetFillStyle(0);
  //  TLegendEntry* E1Ratio =      LegendRatio->AddEntry("", "#bf{ALICE Preliminary}", "");
  TLegendEntry* E1Ratio =      LegendRatio->AddEntry("", "#bf{#color[0]{ALICE Preliminary}}", "");
  TLegendEntry* E3Ratio =           LegendRatio->AddEntry(""/*(TObject*)0*/, "#it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  E3Ratio->SetTextAlign(32);

  TLegend *LegendColor=new TLegend(0.16,0.80,0.5,0.92); 
  LegendColor->SetMargin(0);
  //  LegendColor->AddEntry("", "#bf{ALICE Preliminary}", "");
  LegendColor->AddEntry("", "#color[0]{#bf{ALICE Preliminary}}", "");
  //LegendColor->AddEntry("", "", "");

  TLegend *LegendYields=new TLegend(0.16,0.75,0.5,0.93);
  LegendYields->SetMargin(0);
  //  LegendYields->AddEntry("", "#bf{ALICE Preliminary}", "");
  LegendYields->AddEntry("", "#color[0]{#bf{ALICE Preliminary}}", "");
  //  LegendYields->AddEntry("", "", "");
  LegendYields->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  if (PlotType == 1 || PlotType ==3)     LegendYields->AddEntry(""/*(TObject*)0*/, "h#minusK_{S}^{0} correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  else if (PlotType == 2 || PlotType == 4) LegendYields->AddEntry(""/*(TObject*)0*/, "h#minus#Xi correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", ""); 


  TLegend *legendOneRegion=new TLegend(0.62, 0.8, 0.92, 0.9);

  TLegend *legendRegionJet=new TLegend(0.15, 0.82, 0.61, 0.905);
  legendRegionJet->SetFillStyle(0);
  legendRegionJet->SetMargin(0.1);
  TLegend *legendRegionBulk=new TLegend(0.15, 0.72, 0.61, 0.805);
  legendRegionBulk->SetFillStyle(0);
  legendRegionBulk->SetMargin(0.1);
  TLegend *legendRegionFull=new TLegend(0.15, 0.62, 0.61, 0.705);
  legendRegionFull->SetFillStyle(0);
  legendRegionFull->SetMargin(0.1);
  TLegendEntry * lReAll1Bis[3];
  TLegendEntry *lReAll2Bis[3];

  TLegend *legendRegionAll=new TLegend(0.15, 0.43, 0.59, 0.71);
  legendRegionAll->SetFillStyle(0);
  legendRegionAll->SetMargin(0.07);
  TLegendEntry * lReAll1[3];
  TLegendEntry *lReAll2[3];

  TLegend *legendEnergyBoxColor=new TLegend(0.67, 0.64, 0.90, 0.76);
  TLegend *legendEnergyBox1=new TLegend(0.68, 0.71, 0.90, 0.78);
  TLegend *legendEnergyBox2=new TLegend(0.656, 0.64, 0.90, 0.71);

  TLegend *legendEnergyBoxBis1=new TLegend(0.68, 0.69, 0.89, 0.76);
  TLegend *legendEnergyBoxBis2=new TLegend(0.68, 0.62, 0.89, 0.69);

  TLegend *legendStatBox=new TLegend(0.73, 0.80, 0.93, 0.92);
  TLegend *legendStatBoxBis=new TLegend(0.2, 0.48, 0.35, 0.58);
  TLegend *legendStatBoxColor=new TLegend(0.18, 0.56, 0.35, 0.71);

  TLegend *legendMCTypes = new TLegend(0.66, 0.65, 0.90, 0.76);

  Int_t NLoopRegion =-1;
  Int_t NTypeMC =-1;
  Int_t MCindexForPlotting = 0;
  if (MonashTune==1) MCindexForPlotting = 2;
  else if (MonashTune==2) MCindexForPlotting = 1;
  else if (MonashTune==3) MCindexForPlotting = 2;

  cout << "\n\e[35mPlot type:  " << SPlotType[PlotType] << "\e[39m" << endl;
  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    NLoopRegion++;
    cout << "\nRegion: " << SRegionType[ireg] << endl;
    for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
      if (Coll==1) continue; //for the time being, only 13 TeV
      cout << "\nCollisions: ";
      if (Coll==0) cout << " 13 TeV " << endl;
      else cout << " 5 TeV " << endl;
      NTypeMC=-1;
      for (Int_t isMC = 0; isMC <= NumberOfMCs; isMC++){
	if (MonashTune==1 && isMC==1) continue; //only Ropes 
	if (MonashTune==2 && isMC==2) continue; //only Pythia default
	NTypeMC++;
	cout << "\nisMC? " << MCType[isMC] << endl;
	//Input file 
	if (isMC==0) {
	  pathin[isMC] = "RatiosXiK0s3Systems_" + SPlotType[PlotType];
	  if (ChosenRegion>=0) pathin[isMC] += "_" + RegionType[ChosenRegion];
	  pathin[isMC] += ".root";
	}
	else {
	  pathin[isMC] = "Results_pp13TeVMB_FastMCPrediction"; //temporary solution
	  if (isdNdEtaTriggered)  pathin[isMC]+= "_isdNdEtaTriggered";
	  if (isWingsCorrectionApplied) pathin[isMC] += "_isWingsCorrectionAppliedNew";
	  if (MonashTune==2) pathin[isMC] += "_PythiaMonash";
	  else if (MonashTune==1) pathin[isMC] += "_PythiaRopes";
	  else if (MonashTune==3){
	    if (isMC==1) pathin[isMC] += "_PythiaMonash";
	    else if (isMC==2) pathin[isMC] += "_PythiaRopes";
	  }
	  pathin[isMC]+= ".root";
	}

	cout << "\n\e[35mPathin: " << pathin[isMC] << "\e[39m"<< endl;
	TFile *filein = new TFile(pathin[isMC], "");
	if (!filein) {cout << "Input file not available " << endl; return;}

	//input histo
	if (isMC ==0) {
	  NameHisto13TeV = SPlotType[PlotType] + "_" + RegionType[ireg] + "_13TeV_Stat";
	  NameHisto5TeV = SPlotType[PlotType] + "_" + RegionType[ireg] + "_5TeV_Stat";
	  NameHistoSist13TeV = SPlotType[PlotType] + "_" + RegionType[ireg] + "_13TeV_Sist";
	  NameHistoSist5TeV = SPlotType[PlotType] + "_" + RegionType[ireg] + "_5TeV_Sist";
	}
	else {
	  if (PlotType!=0){
	    if (PlotType<=2)   NameHisto13TeV = "Yield_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_StatErr";
	    else  {
	      if (ireg==0)   NameHisto13TeV = "PtvsMult_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_StatErr";
	      else NameHisto13TeV = "PtvsMult_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_StatErrFS"; 
	    }
	    if (PlotType<=2)   NameHistoSist13TeV = "Yield_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_SistErr";
	    else {
	      if (ireg==0)  NameHistoSist13TeV = "PtvsMult_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_SistErr";	
	      else  NameHistoSist13TeV = "PtvsMult_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_SistErrFS";	
	    }
	  }
	  else {
	    NameHisto13TeV = Form("fHistYieldStatRatio_iregion%i", ireg);
	    NameHistoSist13TeV = Form("fHistYieldSistRatio_iregion%i", ireg);
	  }
	}

	cout << "NameInputHisto " << NameHisto13TeV << endl;
	cout << "NameInputHistoSist " << NameHistoSist13TeV << endl;
	if (Coll==0)	histoYield[ireg][Coll][isMC]= (TH1F*) filein->Get(NameHisto13TeV);
	else 	histoYield[ireg][Coll][isMC]= (TH1F*) filein->Get(NameHisto5TeV);
	if (!histoYield[ireg][Coll][isMC]) {cout <<"no histo " << endl; return;}

	if (Coll==0)	histoYieldSist[ireg][Coll][isMC]= (TH1F*) filein->Get(NameHistoSist13TeV);
	else histoYieldSist[ireg][Coll][isMC]= (TH1F*) filein->Get(NameHistoSist5TeV);
	if (!histoYieldSist[ireg][Coll][isMC]) {cout <<"no histo " << endl; return;}

	cout << "Region: " << Region[ireg]<< endl;
	if (Coll ==0 ) cout << "Got histos for  13 TeV (MB + HM) " << endl;
	else  cout << "Got histos for 5 TeV " << endl;
	for (Int_t b=1; b<=  histoYield[ireg][Coll][isMC]->GetNbinsX(); b++){
	  if (histoYield[ireg][Coll][isMC]->GetBinContent(b) != 0)	 {
	    cout << histoYield[ireg][Coll][isMC]->GetBinContent(b) << " +- " << histoYield[ireg][Coll][isMC]->GetBinError(b) << " (stat.) +-  " <<  histoYieldSist[ireg][Coll][isMC]->GetBinError(b) << " (syst.) " << endl;
	  }
	}

	/* TGRAPH creation */
	if (isMC==0) nummoltDataorMC = 8;
	else nummoltDataorMC = nummolt;
	for (Int_t m=0; m< nummoltDataorMC; m++){
	  Yields[ireg][Coll][isMC][m] = histoYield[ireg][Coll][isMC]->GetBinContent(histoYield[ireg][Coll][isMC]->FindBin(dNdEta[isMC][m]));
	  YieldsErrors[ireg][Coll][isMC][m] = sqrt(pow(histoYield[ireg][Coll][isMC]->GetBinError(histoYield[ireg][Coll][isMC]->FindBin(dNdEta[isMC][m])), 2) +  pow(histoYieldSist[ireg][Coll][isMC]->GetBinError(histoYield[ireg][Coll][isMC]->FindBin(dNdEta[isMC][m])),2) );
	}
	ghistoYield[ireg][Coll][isMC] = new TGraphErrors(nummoltDataorMC,dNdEta[isMC],Yields[ireg][Coll][isMC],dNdEtaError[isMC],YieldsErrors[ireg][Coll][isMC]);
	ghistoYield[ireg][Coll][isMC]->SetName(Form("ghistoYield_reg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	StyleTGraphErrors(ghistoYield[ireg][Coll][isMC], ColorDiff[ireg][Coll], 1, 1, LineStyle[isMC]);

	//Splines to TGraphs
	splineY[ireg][Coll][isMC] = new TSpline3(Form("Spline_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), ghistoYield[ireg][Coll][isMC]);
	sp3= (TSpline3*) splineY[ireg][Coll][isMC] ->Clone(Form("sp3Spline_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	fsplineY[ireg][Coll][isMC]= new TF1(Form("fSpline_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), spline, 0, 40);
	fsplineY[ireg][Coll][isMC]->SetLineColor(ColorDiff[ireg][Coll]);
	fsplineY[ireg][Coll][isMC]->SetLineStyle(1);

	//Fit to yields
	if (ireg!=0){
	  fitY[ireg][Coll][isMC]= new TF1(Form("fitY_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), "pol1", 0, 40);
	  fitY[ireg][Coll][isMC]->SetLineColor(ColorDiff[ireg][Coll]);
	  fitY[ireg][Coll][isMC]->SetLineStyle(7);
	  histoYield[ireg][Coll][isMC]->Fit(fitY[ireg][Coll][isMC], "R0");
	}

	//create histo in order to perform division of TF1s
	hMCDataRatio[ireg][Coll][isMC]= new TH1F (Form("hMCDataRatio_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), Form("fitY_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), 1000, 0, 40);
	//	if (ireg!=0)	hMCDataRatio[ireg][Coll][isMC]->Add(fitY[ireg][Coll][isMC]);
	if (ireg!=0)	hMCDataRatio[ireg][Coll][isMC]->Add(fsplineY[ireg][Coll][isMC]);
	else hMCDataRatio[ireg][Coll][isMC]->Add(fsplineY[ireg][Coll][isMC]);
      
	canvas->cd();
	gStyle->SetLegendBorderSize(0);
        gStyle->SetLegendFillColor(0);
        gStyle->SetLegendFont(42);
        pad1->Draw();
        pad1->cd();

	if (isMC==1)	  ColorDiff[ireg][Coll] =  ColorDiffMC[ireg][Coll];
	else if (isMC==2) 	  ColorDiff[ireg][Coll] =  ColorDiffMCRopes[ireg][Coll];

	StyleHisto(histoYield[ireg][Coll][isMC], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	StyleHisto(histoYieldSist[ireg][Coll][isMC], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);

	if (PlotType==0) {
	  histoYield[ireg][Coll][isMC]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldSist[ireg][Coll][isMC]->GetYaxis()->SetTitleSize(0.07);
	  histoYield[ireg][Coll][isMC]->GetYaxis()->SetTitleOffset(0.7);
	  histoYieldSist[ireg][Coll][isMC]->GetYaxis()->SetTitleOffset(0.7);
	  histoYield[ireg][Coll][isMC]->GetYaxis()->SetTickLength(0.02);
	  histoYieldSist[ireg][Coll][isMC]->GetYaxis()->SetTickLength(0.02);
	}

	//legend
	if (NLoopRegion==0 && Coll==0 && isMC!=0){
	  ghistoYieldGrey[ireg][Coll][isMC] = (TGraphErrors*)ghistoYield[ireg][Coll][isMC]->Clone(Form("ghistoYieldGrey_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	  if (ChosenRegion==-1) 	  ghistoYieldGrey[ireg][Coll][isMC]->SetLineColor(kGray+3);
	  else  {
	    ghistoYieldGrey[ireg][Coll][isMC]->SetLineColor(ColorDiff[ireg][Coll]);
	    ghistoYieldGrey[ireg][Coll][isMC]->SetFillColor(ColorDiff[ireg][Coll]);
	  }
	  //	  if (/*PlotType==0 ||*/ ChosenRegion==0){
	  if (PlotType==0 || ChosenRegion==0){
	    if (isMC==1) ghistoYieldGrey[ireg][Coll][isMC]->SetFillStyle(3001);
	    else ghistoYieldGrey[ireg][Coll][isMC]->SetFillStyle(3002);
	    legendMCTypes->AddEntry(ghistoYieldGrey[ireg][Coll][isMC], MCType[isMC], "f");
	  }
	  else {
	    legendMCTypes->AddEntry(ghistoYieldGrey[ireg][Coll][isMC], MCType[isMC], "l");
	  }
	}

	if (NLoopRegion==0 && Coll==0 && NTypeMC==0){
	  fHistYieldStatGrey[Coll]= (TH1F*) histoYield[ireg][Coll][isMC]->Clone("fHistYieldStatBlack");
	  fHistYieldStatGrey[Coll]->SetLineColor(ColorEnergy[Coll]);
	  fHistYieldStatGrey[Coll]->SetMarkerColor(ColorEnergy[Coll]);
	  fHistYieldStatGrey[Coll]->SetMarkerStyle(20);

	  if (Coll==0){
	    TLegendEntry* L1 = legendEnergyBox1->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 13 TeV", "pe");
	    L1->SetTextAlign(32);
	    L1->SetTextSize(0.042);
	    TLegendEntry* L1Bis = legendEnergyBoxBis1->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 13 TeV", "pe");
	    L1Bis->SetTextAlign(22);
	    L1Bis->SetTextSize(0.044);
	  }
	  else {
	    TLegendEntry* L2 = legendEnergyBox2->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
	    L2->SetTextAlign(32);
	    L2->SetTextSize(0.042);
	    TLegendEntry* L2Bis = legendEnergyBoxBis2->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
	    L2Bis->SetTextAlign(22);
	    L2Bis->SetTextSize(0.044);
	  }
	}

	if (Coll==0 && NTypeMC ==0){

	  if (ireg==0) lReAll1Bis[ireg] = legendRegionJet->AddEntry(histoYield[ireg][Coll][isMC], sRegion[ireg], "p");
	  else if (ireg==1)lReAll1Bis[ireg] = legendRegionBulk->AddEntry(histoYield[ireg][Coll][isMC], sRegion[ireg], "p");
	  else lReAll1Bis[ireg] = legendRegionFull->AddEntry(histoYield[ireg][Coll][isMC], sRegion[ireg], "p");

	  lReAll1Bis[ireg]->SetTextSize(0.045);
	  //      lReAll1Bis[ireg]->SetTextAlign(32);

	  if (ireg==0) lReAll2Bis[ireg]=      legendRegionJet->AddEntry("", sRegion1[ireg], "");
	  else    if (ireg==1) lReAll2Bis[ireg]=      legendRegionBulk->AddEntry("", sRegion1[ireg], "");
	  else    lReAll2Bis[ireg]=      legendRegionFull->AddEntry("", sRegion1[ireg], "");
	  lReAll2Bis[ireg]->SetTextSize(0.035);

	  lReAll1[ireg] = legendRegionAll->AddEntry(histoYield[ireg][Coll][isMC], sRegion[ireg], "p");
	  lReAll1[ireg]->SetTextSize(0.048);
	  lReAll2[ireg]= legendRegionAll->AddEntry("", sRegion1[ireg], "");
	  lReAll2[ireg]->SetTextSize(0.038);

	  if (NLoopRegion==0){
	    fHistYieldStatBlack= (TH1F*)      histoYield[ireg][Coll][isMC]->Clone("fHistYieldStatBlack");
	    fHistYieldSistBlack= (TH1F*)      histoYieldSist[ireg][Coll][isMC]->Clone("fHistYieldSistBlack");
	    fHistYieldStatBlack->SetLineColor(1);
	    fHistYieldStatBlack->SetMarkerColor(1);
	    fHistYieldStatBlack->SetMarkerStyle(20);
	    fHistYieldSistBlack->SetLineColor(1);
	    fHistYieldSistBlack->SetMarkerColor(1);
	    fHistYieldSistBlack->SetMarkerStyle(20);

	    legendStatBoxBis->AddEntry(fHistYieldSistBlack, "syst. error", "ef");
	    legendStatBoxBis->AddEntry(fHistYieldStatBlack, "stat. error", "pe");
	    legendStatBox->AddEntry(fHistYieldSistBlack, "syst. error", "ef");
	    legendStatBox->AddEntry(fHistYieldStatBlack, "stat. error", "pe");

	  }
	}

	if (Coll==0 && NTypeMC==0){
	  legendStatBoxColor->AddEntry(histoYield[ireg][Coll][isMC], "stat. error", "pe");
	  legendStatBoxColor->AddEntry(histoYieldSist[ireg][Coll][isMC], "syst. error", "ef");

	  TLegendEntry * lRe1 = legendOneRegion->AddEntry("", sRegion[ireg], "");
	  lRe1->SetTextSize(0.06);
	  lRe1->SetTextAlign(32);
	  TLegendEntry *leR2=      legendOneRegion->AddEntry("", sRegion1[ireg], "");
	  leR2->SetTextSize(0.04);
	  leR2->SetTextAlign(32);
	}
	//	if (Coll==0)	LegendColor->AddEntry("", NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");

	if (Coll==1)	legendEnergyBoxColor->AddEntry(histoYield[ireg][Coll][isMC], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
	else 	legendEnergyBoxColor->AddEntry(histoYield[ireg][Coll][isMC], "pp, #sqrt{#it{s}} = 13 TeV", "pe");

	if (isMC ==1 || isMC==2){
	  //	  if (/*PlotType==0 ||*/ ChosenRegion==0){
	  if (PlotType==0 || ChosenRegion==0){
	    ghistoYield[ireg][Coll][isMC]->SetFillColor(ColorDiff[ireg][Coll]);
	    if (isMC==1)	  ghistoYield[ireg][Coll][isMC]->SetFillStyle(3001);
	    else 	  ghistoYield[ireg][Coll][isMC]->SetFillStyle(3002); //NO: 3003
	    ghistoYield[ireg][Coll][isMC]->Draw("same 3");
	    //	    legendMCTypes->AddEntry(ghistoYield[ireg][Coll][isMC], MCType[isMC], "f");
	  }
	  else {
	    ghistoYield[ireg][Coll][isMC]->Draw("same");
	  }
	}
	else {
	  histoYield[ireg][Coll][isMC]->Draw("same e0x0");
	  histoYieldSist[ireg][Coll][isMC]->SetFillStyle(0);
	  histoYieldSist[ireg][Coll][isMC]->Draw("same e2");
	  //ghistoYield[ireg][Coll][isMC]->Draw("same");
	}
	//if (ireg!=0)	fitY[ireg][Coll][isMC]->Draw("same");
	//else if (ireg==0)	fsplineY[ireg][Coll][isMC]->Draw("same");

	if (isMC == MCindexForPlotting){
	  if (PlotType==0){
	    LegendRatio->Draw("");
	    legendRegionJet->Draw("");
	    legendRegionBulk->Draw("");
	    legendRegionFull->Draw("");
	    legendStatBoxBis->Draw("");
	  }
	  else {
	    if (ChosenRegion == -1){
	      LegendYields->Draw("");
	      //legendEnergyBox1->Draw("");
	      //legendEnergyBox2->Draw("");
	      legendStatBox->Draw("");
	      legendRegionAll->Draw("");
	    }
	    else {
	      LegendYields->Draw("");
	      //LegendColor->Draw("");
	      legendStatBoxColor->Draw("");
	      legendOneRegion->Draw("");
	      //legendEnergyBoxColor->Draw("");
	    }
 	  }
	  legendMCTypes->Draw("");
	}

	//MC / DATA Ratios
	if (isMC >0){
	  //Ratios between splines (jet) or fits (OOJ, full) -> to be used to fill a TGraph
	  hMCDataRatio[ireg][Coll][isMC]->Divide(hMCDataRatio[ireg][Coll][0]);
	}
      }
    }
  }

  canvas->cd();
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  pad2->Draw();
  pad2->cd();

  TH1F* hdummy = new TH1F ("hdummy", "hdummy",  1000, 0, 40); 
  TF1 * lineAt1 = new TF1 ("lineAt1", "pol0", 0, 40);
  lineAt1->SetParameter(0,1);
  lineAt1->SetLineColor(1);
  lineAt1->SetLineStyle(2);

  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
      if (Coll==1) continue; //for the time being, only 13 TeV
      for (Int_t isMC = 0; isMC <= NumberOfMCs; isMC++){
	if (isMC==0) continue;
	if (MonashTune==1 && isMC!=2) continue; //only Ropes 
	if (MonashTune==2 && isMC!=1) continue; //only Pythia default

	StyleHisto(hMCDataRatio[ireg][Coll][isMC], LowRatio, UpRatio, ColorDiff[ireg][Coll], 1, titleX, titleYMCToData, "" , 1, 0, 40, xOffset, yOffset, 1);

	hMCDataRatio[ireg][Coll][isMC]->GetYaxis()->SetTitleSize(0.07);
	hMCDataRatio[ireg][Coll][isMC]->GetYaxis()->SetTitleOffset(0.7);

	hMCDataRatio[ireg][Coll][isMC]->GetYaxis()->SetLabelSize(0.10);
	hMCDataRatio[ireg][Coll][isMC]->GetYaxis()->SetLabelOffset(0.009);

	hMCDataRatio[ireg][Coll][isMC]->GetXaxis()->SetTitleSize(0.10);
	hMCDataRatio[ireg][Coll][isMC]->GetXaxis()->SetTitleOffset(1.3);

	hMCDataRatio[ireg][Coll][isMC]->GetXaxis()->SetLabelSize(0.10);
	hMCDataRatio[ireg][Coll][isMC]->GetXaxis()->SetLabelOffset(0.02);

	hMCDataRatio[ireg][Coll][isMC]->GetYaxis()->SetNdivisions(305);
	hMCDataRatio[ireg][Coll][isMC]->GetYaxis()->SetTickLength(0.02);

	//hMCDataRatio[ireg][Coll][isMC]->Draw("same");

	StyleHisto(hdummy, LowRatio, UpRatio, ColorDiff[ireg][Coll], 1, titleX, titleYMCToData, "" , 1, 0, 40, xOffset, yOffset, 1);

	hdummy->GetYaxis()->SetTitleSize(0.1);
	hdummy->GetYaxis()->SetTitleOffset(0.7);

	hdummy->GetYaxis()->SetLabelSize(0.10);
	hdummy->GetYaxis()->SetLabelOffset(0.009);

	hdummy->GetXaxis()->SetTitleSize(0.10);
	hdummy->GetXaxis()->SetTitleOffset(1.3);

	hdummy->GetXaxis()->SetLabelSize(0.10);
	hdummy->GetXaxis()->SetLabelOffset(0.02);

	hdummy->GetYaxis()->SetNdivisions(305);
	hdummy->GetYaxis()->SetTickLength(0.02);

	hdummy->Draw("same");

	for (Int_t m=0; m< nummolt; m++){
	  YieldRatios[ireg][Coll][isMC][m] = hMCDataRatio[ireg][Coll][isMC]->GetBinContent(hMCDataRatio[ireg][Coll][isMC]->FindBin(dNdEta[isMC][m]));
	  YieldRatiosErrors[ireg][Coll][isMC][m] = 0; 
	  //	  cout << "m " << dNdEta[isMC][m] << " " << YieldRatios[ireg][Coll][isMC][m] << endl;
	}
	ghistoYieldRatio[ireg][Coll][isMC] = new TGraphErrors(nummolt,dNdEta[isMC],YieldRatios[ireg][Coll][isMC],dNdEtaError[isMC],YieldRatiosErrors[ireg][Coll][isMC]);
	ghistoYieldRatio[ireg][Coll][isMC]->SetName(Form("ghistoYieldRatio_reg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	StyleTGraphErrors(ghistoYieldRatio[ireg][Coll][isMC], ColorDiff[ireg][Coll], 1, 1, LineStyle[isMC]);

	ghistoYieldRatio[ireg][Coll][isMC]->Draw("same");
	lineAt1->Draw("same");
      }
    }
  }
	
  /*
    if (PlotType==0){
    LegendRatio->Draw("");
    legendRegionJet->Draw("");
    legendRegionBulk->Draw("");
    legendRegionFull->Draw("");
    legendStatBoxBis->Draw("");
    if (!isFit){
    legendEnergyBoxBis1->Draw("");
    legendEnergyBoxBis2->Draw("");
    }
    if (isFit){
    legendFit->Draw("");
    }
    }
    else if (PlotType==1 || PlotType==2){
    if (ChosenRegion == -1){
    LegendYields->Draw("");
    legendEnergyBox1->Draw("");
    legendEnergyBox2->Draw("");
    legendStatBox->Draw("");
    legendRegionAll->Draw("");
    }
    else {
    LegendColor->Draw("");
    legendStatBoxColor->Draw("");
    legendOneRegion->Draw("");
    legendEnergyBoxColor->Draw("");
    }
    }
  */

  if (ChosenRegion<0){
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType] + MCTypeBis[MonashTune-1]+".pdf");
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType] + MCTypeBis[MonashTune-1]+".eps");
  }
  else {
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType]+"_" + RegionType[ChosenRegion]+ MCTypeBis[MonashTune-1]+".pdf");
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType]+"_" + RegionType[ChosenRegion]+ MCTypeBis[MonashTune-1]+".eps");
  }

  TFile *fileout = new TFile("CompareDatavsMC.root", "RECREATE");
  fileout->WriteTObject(canvas);
  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
      if (Coll==1) continue; //for the time being, only 13 TeV
      for (Int_t isMC = 0; isMC <= NumberOfMCs; isMC++){
	if (isMC==0) continue;
	if (MonashTune==1 && isMC!=2) continue; //only Ropes 
	if (MonashTune==2 && isMC!=1) continue; //only Pythia default
	ghistoYieldRatio[ireg][Coll][isMC]->Write();
      }
    }
  }
  fileout->Close();
}
