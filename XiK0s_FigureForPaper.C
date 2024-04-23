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
#include <TGraphAsymmErrors.h>
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
#include "Macros/constants.h"

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

void StyleTGraphErrors(TGraphAsymmErrors *tgraph, Int_t color, Int_t style, Float_t mSize, Int_t linestyle){
  tgraph->SetLineColor(color);
  tgraph->SetLineWidth(3);
  tgraph->SetMarkerColor(color);
  tgraph->SetMarkerStyle(style);
  tgraph->SetMarkerSize(mSize);
  tgraph->SetLineStyle(linestyle);
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

TString TitleYieldRatio="#Xi/K^{0}_{S} yield ratio vs multiplicity";
TString titleRelUnc = "Relative uncertainty";
TString titleMult = "Multiplicity class ";
//TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5, #it{p}_{T,trigg}>3 GeV/#it{c}}";
TString titleYield[2]={"K^{0}_{S}", "#Xi"};
TString TitleYYieldRatio="#it{N}_{#Xi}/#it{N}_{K^{0}_{S}}";

//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})", "#it{N}_{#Xi}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})"};
TString titleYieldYType[2]={"#it{N}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})", "#it{N}_{#Xi}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})"};
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString RegionTypeBis[3] = {"Jet", "OOj", "Full"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
//TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
//TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
TString sRegion[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};
//TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};
//TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
TString sRegion1[3]={"|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 1.1", "0.86 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
//TString sRegion1[3]={"|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 1.1", "0.86 < |#Delta#it{#eta}| < 1.2, 1.1 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};

TString titleYToOOJ = "Tow. / Transv.";
TString titleYMCToData = "Model / Data";

const Int_t NumberOfMCs=3;
const Int_t numColls =2; //pp13 TeV (MB + HM), pp5TeV
const Int_t numRegions =3;
const Int_t numDataorMC = NumberOfMCs+1;
Int_t nummoltMax = nummolt;
Int_t nummoltMaxDataorMC = 0;

void XiK0s_FigureForPaper(Int_t ChosenRegion = 0,  Float_t ScalingFactorXiK0s = 1/*0.8458/1.08747*/, Bool_t isChangesIncluded=1, Bool_t isFit=0, Bool_t isdNdEtaTriggered=1, Bool_t MaterialBudgetCorr=1){

  if (ChosenRegion < 0 || ChosenRegion > 2) {
    cout << "Chosen Region must be 0, 1 or 2 " << endl;
    return;
  }
  Double_t ScalingXi = 29.0692; //30
  Float_t UpperValueX = 40;

  Int_t nummoltDataorMC =0;

  Int_t MonashTune =0;
  cout << "Should be analyse Ropes (=1) or MonashDefault (=2) or EPOSLHC (=3) or only Pythia (=4) or all of them (=5)? " << endl;
  cin >> MonashTune;

  Bool_t isWingsCorrectionApplied=0;
  cout <<"Do you want to analyse the files with the wings correction applied (not for Xi DATA)? Type 1 if you DO want" << endl;
  cin >> isWingsCorrectionApplied;

  //Set titles
  TString titleY;
  titleY = titleYieldYType[0];

  gStyle->SetOptStat(0);

  //STYLE 
  const Int_t numParticles = 2;
  TString titleX = titledNdeta;
  Float_t MarkerSize[3] ={2, 2, 3};
  //Int_t MarkerType[3][numColls] = {{24,20} , {25, 21}, {27, 33}}; //thinner
  Int_t MarkerType[3][numColls] = {{53,20} , {54, 21}, {56, 33}}; //thicker marker (requires ROOT6)
  //Int_t MarkerType[3][numColls] = {{71,20} , {72, 21}, {74, 33}}; //thicker marker (requires ROOT6)
  //  Int_t MarkerType[3][numColls] = {{89,20} , {90, 21}, {92, 33}}; //thicker marker (requires ROOT6)
  Int_t Color[numRegions] = {628,418,601};
  Int_t ColorEnergy[numColls] = {1, 921};// {922,920};
  //  Int_t ColorDiff[numParticles][numColls] = {{634,628}, {418,829}};
  Int_t XiColor = kAzure+7;
  Int_t XiColor5TeV = kAzure + 8;
  Int_t Orange = 808;
  Int_t OrangePale = 806;
  Int_t VioletPale = 871;
  Int_t XiColorPale =   kAzure - 9;
  Int_t K0sColorPale =   kPink + 1;
  Int_t K0sColor = kPink + 10;
  Int_t K0sColor5TeV = kPink + 6;
  Int_t K0sColorLine = kPink + 8;
  Int_t XiColorLine = kAzure - 3;
  Int_t ColorDiff[numParticles][numColls] = {{K0sColor, K0sColor5TeV}, {XiColor,XiColor5TeV}};
  Int_t ColorDiffPale[numParticles][numColls] = {{K0sColorPale, K0sColor5TeV}, {XiColorPale,XiColor5TeV}};
  //Int_t ColorDiffMC[numParticles][numColls] = {{634,628}, {418,829}};
  Int_t ColorDiffMC[numParticles][numColls] = {{K0sColor, K0sColor5TeV}, {XiColor,XiColor5TeV}};
  Int_t ColorDiffLine[numParticles][numColls] = {{K0sColorLine, K0sColorLine}, {XiColorLine,XiColorLine}};
  //Int_t ColorDiffMCRopes[numParticles][numColls] = {{634,628}, {418,829}};
  Int_t ColorDiffMCRopes[numParticles][numColls] = {{K0sColor, K0sColor5TeV}, {XiColor,XiColor5TeV}};
  Int_t LineStyle[numDataorMC] = {1, 1, 8, 3};
  Float_t AlphaColor[numDataorMC] = {1, 0.9, 0.5, 1};
  Float_t FillStyle[numDataorMC] = {1, 1001, 1001, 3001}; //3001

  //ChosenRegion == -1 //all regions (JET, OOJ, FULL) are plotted
  if (ChosenRegion>2) {cout << "Choose a valid region! " << endl; return;}

  TString CollisionsComp[2] = {"_HMMultBinning1_vsHM", "_vs5TeV"};
  TString SPlotType[2] = {"K0s", "Xi"};
  TString SPlotTypeFinal[2] = {"K0s", "Xi"};
  TString SPlotTypeBis[2] = {"K0s", "Xi"};
  TString SisFit[2] = {"", "_isFit"};
  TString Region[numRegions] = {"Jet", "Bulk", "All"};
  TString RegionBis[numRegions] = {"Jet", "Bulk", "Inclusive"};

  Float_t Up= 0.14-10e-6;
  Float_t Low = 0.02+10e-6;
  Low = 10e-6; Up = 0.50-10e-6; //0.45
  if (ChosenRegion==0) {Low = 0.015+10e-6; Up = 0.10-10e-6;} //0.035 

  Float_t LowRatio = 0.5+10e-4;
  Float_t UpRatio = 2.5-10e-4;
  LowRatio = 0.3+10e-4;
  UpRatio = 2.2-10e-4;

  TF1 * fitToRatioToOOJ = new TF1("fitToRatioToOOJ", "pol0", 0, UpperValueX);
  fitToRatioToOOJ->SetLineColor(1);
  fitToRatioToOOJ->SetLineStyle(2);
  fitToRatioToOOJ->SetLineWidth(2);

  TF1* fitPol0[numDataorMC];
  TF1* fitPol1[numDataorMC];
  TF1* fitPol2[numDataorMC];

  TH1F *histoYield[numRegions][numColls][numDataorMC][numParticles];
  TH1F *histoYieldSist[numRegions][numColls][numDataorMC][numParticles];
  TH1F *histoYieldSistMultUnCorr[numRegions][numColls][numDataorMC][numParticles];
  TH1F *histoYieldStatSistMultUnCorr[numRegions][numColls][numDataorMC][numParticles];
  TH1F *histoYieldRatio[numRegions][numColls][numDataorMC][numParticles];
  Float_t MinValueYield[numRegions][numColls][numDataorMC][numParticles];
  Float_t MaxValueYield[numRegions][numColls][numDataorMC][numParticles];
  Float_t MinValueError[numRegions][numColls][numDataorMC][numParticles];
  Float_t MaxValueError[numRegions][numColls][numDataorMC][numParticles];
  TH1F *hMCDataRatio[numRegions][numColls][numDataorMC][numParticles];
  TSpline3 * splineY[numRegions][numColls][numDataorMC][numParticles];
  TF1* fsplineY[numRegions][numColls][numDataorMC][numParticles];
  TF1* fitY[numRegions][numColls][numDataorMC][numParticles];
  TF1* fitRatios[numRegions][numColls][numDataorMC][numParticles];

  Float_t Yields[numRegions][numColls][numDataorMC][numParticles][nummoltMax];
  Float_t YieldsErrors[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};
  Float_t YieldsErrorsStat[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};
  Float_t YieldsErrorsSist[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};
  Float_t YieldsErrorsSistMultUnCorr[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};

  Float_t YieldRatios[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};
  Float_t YieldRatiosErrors[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};
  Float_t YieldRatiosDATA[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};
  Float_t YieldRatiosErrorsDATA[numRegions][numColls][numDataorMC][numParticles][nummoltMax];// = {0};

  TGraphAsymmErrors*	ghistoYield[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldLine[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldSist[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldSistMultUnCorr[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldGrey[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldGreyBis[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldRed[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldRatio[numRegions][numColls][numDataorMC][numParticles];
  TGraphAsymmErrors*	ghistoYieldRatioDATA[numRegions][numColls][numDataorMC][numParticles];

  TH1F*  fHistYieldStatBlack;
  TH1F*  fHistYieldSistBlack;
  TH1F*  fHistYieldSistMultUnCorrBlack;
  TH1F*  fHistYieldStatGrey[numColls];
  TH1F*  fHistYieldSistGrey[numColls];
  TH1F*  fHistYieldSistMultUnCorrGrey[numColls];

  TString NameHisto13TeV="";
  TString NameHistoSist13TeV="";
  TString NameHistoSistMultUnCorr13TeV="";
  TString NameHisto5TeV="";
  TString NameHistoSist5TeV="";
  TString NameHistoSistMultUnCorr5TeV="";
  TString NameHistoFinal[numRegions][numColls][numDataorMC];

  TString MCType[NumberOfMCs+1]= {"Data", "PYTHIA8 Monash", "PYTHIA8 Ropes", "EPOS LHC"};
  TString MCTypeBis[5]= {"PythiaRopes", "PythiaMonash", "EPOSLHC", "AllPythia", "AllMC"};
  
  Float_t dNdEta[numDataorMC][numColls][nummoltMax];// = {0};
  Float_t dNdEtaErrorL[numDataorMC][numColls][nummoltMax];// = {0};
  Float_t dNdEtaErrorR[numDataorMC][numColls][nummoltMax];// = {0};
  /* MonashRopes, events with trigger particle */
  if (isdNdEtaTriggered){
    //DATA 13 TeV
    for (Int_t m=0; m<nummoltMax; m++){ 
      dNdEta[0][0][m] =  dNdEtaFinal13TeV[m];
      dNdEtaErrorL[0][0][m] =  dNdEtaFinal13TeV_ErrorL[m];
      dNdEtaErrorR[0][0][m] =  dNdEtaFinal13TeV_ErrorR[m];
      cout <<  "dNdeta values from macro " << dNdEta[0][0][m] << " + " << dNdEtaErrorR[0][0][m] << " - " <<  dNdEtaErrorL[0][0][m] << endl;
    }
    //DATA 5 TeV    
    for (Int_t m=0; m<nummoltMax; m++){ 
      dNdEta[0][1][m] =  dNdEtaFinal5TeV[m];
      dNdEtaErrorL[0][1][m] =  dNdEtaFinal5TeV_ErrorL[m];
      dNdEtaErrorR[0][1][m] =  dNdEtaFinal5TeV_ErrorR[m];
    }
    //Monash
    dNdEta[2][0][0] = 8.45;
    dNdEta[2][0][1] = 12.42;
    dNdEta[2][0][2] = 15.11;
    dNdEta[2][0][3] = 17.15;
    dNdEta[2][0][4] = 19.12;
    dNdEta[2][0][5] = 20.99;
    dNdEta[2][0][6] = 22.75;
    dNdEta[2][0][7] = 24.84;
    dNdEta[2][0][8] = 27.36;
    dNdEta[2][0][9] = 30.35;

    dNdEtaErrorL[2][0][0] = 0;
    dNdEtaErrorL[2][0][1] = 0;
    dNdEtaErrorL[2][0][2] = 0;
    dNdEtaErrorL[2][0][3] = 0;
    dNdEtaErrorL[2][0][4] = 0;
    dNdEtaErrorL[2][0][5] = 0;
    dNdEtaErrorL[2][0][6] = 0;
    dNdEtaErrorL[2][0][7] = 0;
    dNdEtaErrorL[2][0][8] = 0;
    dNdEtaErrorL[2][0][9] = 0;

    dNdEtaErrorR[2][0][0] = 0;
    dNdEtaErrorR[2][0][1] = 0;
    dNdEtaErrorR[2][0][2] = 0;
    dNdEtaErrorR[2][0][3] = 0;
    dNdEtaErrorR[2][0][4] = 0;
    dNdEtaErrorR[2][0][5] = 0;
    dNdEtaErrorR[2][0][6] = 0;
    dNdEtaErrorR[2][0][7] = 0;
    dNdEtaErrorR[2][0][8] = 0;
    dNdEtaErrorR[2][0][9] = 0;

    //Ropes
    dNdEta[1][0][0] = 8.27;
    dNdEta[1][0][1] = 12.48;
    dNdEta[1][0][2] = 15.23;
    dNdEta[1][0][3] = 17.26;
    dNdEta[1][0][4] = 19.22;
    dNdEta[1][0][5] = 21.1;
    dNdEta[1][0][6] = 22.9;
    dNdEta[1][0][7] = 25.1;
    dNdEta[1][0][8] = 27.81;
    dNdEta[1][0][9] = 31.26;

    dNdEtaErrorL[1][0][0] = 0;
    dNdEtaErrorL[1][0][1] = 0;
    dNdEtaErrorL[1][0][2] = 0;
    dNdEtaErrorL[1][0][3] = 0;
    dNdEtaErrorL[1][0][4] = 0;
    dNdEtaErrorL[1][0][5] = 0;
    dNdEtaErrorL[1][0][6] = 0;
    dNdEtaErrorL[1][0][7] = 0;
    dNdEtaErrorL[1][0][8] = 0;
    dNdEtaErrorL[1][0][9] = 0;

    dNdEtaErrorR[1][0][0] = 0;
    dNdEtaErrorR[1][0][1] = 0;
    dNdEtaErrorR[1][0][2] = 0;
    dNdEtaErrorR[1][0][3] = 0;
    dNdEtaErrorR[1][0][4] = 0;
    dNdEtaErrorR[1][0][5] = 0;
    dNdEtaErrorR[1][0][6] = 0;
    dNdEtaErrorR[1][0][7] = 0;
    dNdEtaErrorR[1][0][8] = 0;
    dNdEtaErrorR[1][0][9] = 0;

    //EPOS LHC
    dNdEta[3][0][0] = 7.97;
    dNdEta[3][0][1] = 11.63;
    dNdEta[3][0][2] = 14.26;
    dNdEta[3][0][3] = 16.29;
    dNdEta[3][0][4] = 18.29;
    dNdEta[3][0][5] = 20.25;
    dNdEta[3][0][6] = 22.16;
    dNdEta[3][0][7] = 24.55;
    dNdEta[3][0][8] = 27.49;
    dNdEta[3][0][9] = 31.61;

    dNdEtaErrorL[3][0][0] = 0;
    dNdEtaErrorL[3][0][1] = 0;
    dNdEtaErrorL[3][0][2] = 0;
    dNdEtaErrorL[3][0][3] = 0;
    dNdEtaErrorL[3][0][4] = 0;
    dNdEtaErrorL[3][0][5] = 0;
    dNdEtaErrorL[3][0][6] = 0;
    dNdEtaErrorL[3][0][7] = 0;
    dNdEtaErrorL[3][0][8] = 0;
    dNdEtaErrorL[3][0][9] = 0;

    dNdEtaErrorR[3][0][0] = 0;
    dNdEtaErrorR[3][0][1] = 0;
    dNdEtaErrorR[3][0][2] = 0;
    dNdEtaErrorR[3][0][3] = 0;
    dNdEtaErrorR[3][0][4] = 0;
    dNdEtaErrorR[3][0][5] = 0;
    dNdEtaErrorR[3][0][6] = 0;
    dNdEtaErrorR[3][0][7] = 0;
    dNdEtaErrorR[3][0][8] = 0;
    dNdEtaErrorR[3][0][9] = 0;

  }
  else {    
    //DATA
    dNdEta[0][0][0] = 3.33;
    dNdEta[0][0][1] = 7.14;
    dNdEta[0][0][2] = 11.46;
    dNdEta[0][0][3] = 16.17;
    dNdEta[0][0][4] = 21.2;
    dNdEta[0][0][5] = 30.43;
    dNdEta[0][0][6] = 32.57;
    dNdEta[0][0][7] = 36.29;
    dNdEta[0][0][8] = 0;
    dNdEta[0][0][9] = 0;

    //Monash
    dNdEta[1][0][0] = 4.01;
    dNdEta[1][0][1] = 9.04;
    dNdEta[1][0][2] = 12.1;
    dNdEta[1][0][3] = 14.29;
    dNdEta[1][0][4] = 16.4;
    dNdEta[1][0][5] = 18.46;
    dNdEta[1][0][6] = 20.37;
    dNdEta[1][0][7] = 22.68;
    dNdEta[1][0][8] = 25.58;
    dNdEta[1][0][9] = 29.09;

    //Ropes (now it's Monash...)
    dNdEta[2][0][0] = 4.01;
    dNdEta[2][0][1] = 9.04;
    dNdEta[2][0][2] = 12.1;
    dNdEta[2][0][3] = 14.29;
    dNdEta[2][0][4] = 16.4;
    dNdEta[2][0][5] = 18.46;
    dNdEta[2][0][6] = 20.37;
    dNdEta[2][0][7] = 22.68;
    dNdEta[2][0][8] = 25.58;
    dNdEta[2][0][9] = 29.09;

    for (Int_t iMC=0; iMC<numDataorMC; iMC++){
      for (Int_t jcoll=0; jcoll<numColls; jcoll++){
	for (Int_t mmult=0; mmult<nummoltMax; mmult++){
	  dNdEtaErrorL[iMC][jcoll][mmult] = 0;
	  dNdEtaErrorR[iMC][jcoll][mmult] = 0;
	}
      }
    }
  }

  TString pathin[NumberOfMCs+2] ={""};
  TString YieldOrAvgPt[5] = {"Yield", "Yield", "Yield", "AvgPt","AvgPt"};

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1100, 1000); //1100, 900
  canvas->SetFillColor(0);
  Float_t UpperPadValue = 0.36; 
  TPad* pad1 = new TPad( "pad1" ,"pad1" ,0 , UpperPadValue ,1 ,1);
  TPad* pad2 = new TPad( "pad2" ,"pad2" ,0 ,0.0 ,1 , UpperPadValue);

  StylePad(pad1, 0.11, 0.02, 0.03, 0.005); //L, R, T, B
  StylePad(pad2, 0.11, 0.02, 0.02, 0.33); //L, R, T, B

  //legend
  TLegend *legendFit = new TLegend(0.62,0.72,0.9,0.92);
  TLegend *legendFitRatioToOOJ = new TLegend(0.7,0.8,0.8,0.96);
  TLegendEntry * legendFitRatio1;

  //  TLegend *LegendRatio=new TLegend(0.62,0.72,0.9,0.92);
  TLegend *LegendRatio=new TLegend(0.65,0.77,0.93,0.92);
  LegendRatio->SetFillStyle(0);
  //  TLegendEntry* E1Ratio =      LegendRatio->AddEntry("", "#bf{ALICE Preliminary}", "");
  TLegendEntry* E1Ratio =  LegendRatio->AddEntry("", "ALICE, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  //  TLegendEntry* E3Ratio =  LegendRatio->AddEntry(""/*(TObject*)0*/, "#it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  E1Ratio->SetTextSize(0.05);
  //E3Ratio->SetTextSize(0.05);
  E1Ratio->SetTextAlign(32);
  //E3Ratio->SetTextAlign(32);

  TLegend *LegendColor=new TLegend(0.16,0.80,0.5,0.92); 
  LegendColor->SetMargin(0);
  //  LegendColor->AddEntry("", "#bf{ALICE Preliminary}", "");
  //LegendColor->AddEntry("", "#color[0]{#bf{ALICE Preliminary}}", "");
  LegendColor->AddEntry("", "#bf{This work}", "");
  //LegendColor->AddEntry("", "", "");

  //  TLegend *LegendYields=new TLegend(0.16,0.75,0.5,0.93);
  TLegend *LegendYields=new TLegend(0.14,0.85,0.5,0.93);
  LegendYields->SetMargin(0);
  LegendYields->SetTextSize(0.05);
  //LegendYields->AddEntry("", "ALICE pp, #sqrt{#it{s}} = 13 TeV", "");
  LegendYields->AddEntry("", "ALICE, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  //  LegendYields->AddEntry(""/*(TObject*)0*/, "h#minusK_{S}^{0} correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  //  TLegend *LegendParticle=new TLegend(0.14,0.7,0.4,0.85);
  TLegend *LegendParticle=new TLegend(0.14,0.7,0.35,0.85);
  LegendParticle->SetTextSize(0.042);

  TLegend *legendOneRegion=new TLegend(0.62, 0.8, 0.92, 0.9);

  //OLD before EB: 0.15, 0.84, 0.61, 0.915
  //TLegend *legendRegionJet=new TLegend(0.3, 0.8, 0.66, 0.875);
  TLegend *legendRegionJet=new TLegend(0.15, 0.84, 0.61, 0.915);
  legendRegionJet->SetFillStyle(0);
  //EB legendRegionJet->SetMargin(0.1);
  legendRegionJet->SetMargin(0);
  //  TLegend *legendRegionBulk=new TLegend(0.3, 0.71, 0.66, 0.785);
  TLegend *legendRegionBulk=new TLegend(0.15, 0.75, 0.61, 0.825);
  legendRegionBulk->SetFillStyle(0);
  //EB legendRegionBulk->SetMargin(0.1);
  legendRegionBulk->SetMargin(0);
  //  TLegend *legendRegionFull=new TLegend(0.3, 0.62, 0.66, 0.695);
  TLegend *legendRegionFull=new TLegend(0.15, 0.66, 0.61, 0.735);
  legendRegionFull->SetFillStyle(0);
  //EB legendRegionFull->SetMargin(0.1);
  legendRegionFull->SetMargin(0);
  TLegendEntry * lReAll1Bis[3];
  TLegendEntry *lReAll2Bis[3];

  TLegend *legendRegionAllNew=new TLegend(0.15, 0.765, 0.61, 0.84);
  legendRegionAllNew->SetFillStyle(0);
  legendRegionAllNew->SetMargin(0.1);

  //TLegend *legendRegionAll=new TLegend(0.15, 0.35, 0.59, 0.68);
  TLegend *legendRegionAll=new TLegend(0.15, 0.42, 0.59, 0.74);
  legendRegionAll->SetFillStyle(0);
  legendRegionAll->SetMargin(0.07);
  TLegendEntry * lReAll1[3];
  TLegendEntry *lReAll2[3];

  TLegend *legendEnergyBoxColor=new TLegend(0.16, 0.54, 0.31, 0.7);
  //TLegend *legendEnergyBoxColor1=new TLegend(0.15, 0.53, 0.31, 0.7);
  //TLegend *legendEnergyBoxColor3=new TLegend(0.18, 0.53, 0.35, 0.7);
  TLegend *legendEnergyBoxColor1=new TLegend(0.67, 0.59, 0.87, 0.76);
  legendEnergyBoxColor1->SetMargin(0.2);
  TLegend *legendEnergyBoxColor3=new TLegend(0.7, 0.59, 0.87, 0.76);
  legendEnergyBoxColor3->SetTextSize(0.042);

  TLegend *legendEnergyBoxColorBis1=new TLegend(0.65, 0.77, 0.81, 0.89);
  TLegend *legendEnergyBoxColorBis2=new TLegend(0.68, 0.77, 0.84, 0.89);
  TLegend *legendEnergyBoxColorBis3=new TLegend(0.71, 0.77, 0.87, 0.89);
  legendEnergyBoxColorBis3->SetTextSize(0.042);

  TLegend *legendEnergyBox1=new TLegend(0.15, 0.53, 0.31, 0.64);
  TLegend *legendEnergyBox2=new TLegend(0.65, 0.77, 0.81, 0.89);

  TLegend *legendEnergyBoxHelloUpper=new TLegend(0.15, 0.875, 0.28, 0.95);
  legendEnergyBoxHelloUpper->SetHeader("pp");
  legendEnergyBoxHelloUpper->SetNColumns(2);
  legendEnergyBoxHelloUpper->SetMargin(0);
  legendEnergyBoxHelloUpper->SetFillStyle(0);
  //legendEnergyBoxHelloUpper->AddEntry("", "#sqrt{#it{s}} = 13 TeV", "p");
  legendEnergyBoxHelloUpper->AddEntry("", "13 TeV", "p");
  //legendEnergyBoxHelloUpper->AddEntry("", "#sqrt{#it{s}} = 5.02 TeV", "p");
  legendEnergyBoxHelloUpper->AddEntry("", "5 TeV", "p");
  legendEnergyBoxHelloUpper->SetTextSize(0.042);
  TLegend *legendEnergyBoxHelloLower=new TLegend(0.15, 0.62, 0.28, 0.875);
  legendEnergyBoxHelloLower->SetTextSize(0.042);
  legendEnergyBoxHelloLower->SetNColumns(2);
  legendEnergyBoxHelloLower->SetMargin(0);

  //  TLegend *legendStatBox=new TLegend(0.45, 0.54, 0.6, 0.7);
  TLegend *legendStatBox=new TLegend(0.45, 0.5, 0.6, 0.66);
  legendStatBox->SetTextSize(0.04);
  //TLegend *legendStatBoxBis=new TLegend(0.19, 0.34, 0.34, 0.50);
  TLegend *legendStatBoxBis=new TLegend(0.18, 0.32, 0.3, 0.48);
  legendStatBoxBis->SetTextSize(0.03);
  TLegend *legendStatBoxTer=new TLegend(0.16, 0.35, 0.31, 0.54);
  TLegend *legendStatBoxColor=new TLegend(0.45, 0.54, 0.6, 0.7);
  legendStatBoxColor->SetTextSize(0.04);

  TLegend *legendMCTypes;
  //if (PlotType==0 && ChosenRegion==-1)  legendMCTypes = new TLegend(0.69, 0.6, 0.93, 0.73);
  //legendMCTypes = new TLegend(0.69, 0.54, 0.91, 0.71);
  legendMCTypes = new TLegend(0.16, 0.46, 0.41, 0.63);
  legendMCTypes->SetHeader("pp, #sqrt{#it{s}} = 13 TeV");
  legendMCTypes->SetTextSize(0.035);

  TLegend *legendMCTypesBis = new TLegend(0.69, 0.6, 0.93, 0.73);

  TLegend *smalllegendMCTypes;
  smalllegendMCTypes = new TLegend(0.14, 0.68, 0.27, 0.92);
  smalllegendMCTypes->SetTextSize(0.05);
  smalllegendMCTypes->SetMargin(0.2);

  TLegend *smalllegendMCTypesBis= new TLegend(0.15, 0.68, 0.28, 0.92);
  smalllegendMCTypesBis->SetTextSize(0.05);
  smalllegendMCTypesBis->SetMargin(0.2);

  Int_t NLoopRegion =-1;
  Int_t NTypeMC =-1;
  Int_t MCindexForPlotting = 0;
  if (MonashTune==1) MCindexForPlotting = 2;
  else if (MonashTune==2) MCindexForPlotting = 1;
  else if (MonashTune==3) MCindexForPlotting = 3;
  else if (MonashTune==4) MCindexForPlotting = 2;
  else if (MonashTune==5) MCindexForPlotting = 3;

  //DUMMY HISTO FOR AXES
  Float_t xTitle = 33;
  Float_t xOffset = 4; //1.6

  Float_t yTitle = 30; 
  Float_t yOffset = 1.8; 

  Float_t xLabel = 27;
  Float_t yLabel = 27;
  Float_t xLabelOffset = 0.03;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.03;
  Float_t tickXL = 0.06;
  Float_t tickYL = 0.045;

  TH1F*histoYieldDummy= new TH1F("histoYieldDummy", "histoYieldDummy", 100, 0, UpperValueX);
  StyleHisto(histoYieldDummy, Low, Up, 1, 1, titleX, titleY, "" , 1,0, UpperValueX, xOffset, yOffset, 1);
  SetFont(histoYieldDummy);
  SetHistoTextSize(histoYieldDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  histoYieldDummy->GetYaxis()->SetDecimals(kTRUE);
  SetTickLength(histoYieldDummy, tickX, tickY);

  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    NLoopRegion++;
    cout << "\nRegion: " << SRegionType[ireg] << endl;
    for (Int_t ParticleType=0; ParticleType < numParticles; ParticleType++){ //K0s and Xi
      cout << "Particle: " << ParticleType << endl;
      for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
	cout << "\nCollisions: ";
	if (Coll==0) cout << " 13 TeV " << endl;
	else cout << " 5 TeV " << endl;
	NTypeMC=-1;
	for (Int_t isMC = 0; isMC <= NumberOfMCs; isMC++){
	  if (Coll==1 && isMC!=0) continue; //5 TeV only for DATA
	  if (MonashTune==1 && isMC!=2 && isMC!=0) continue; //only Ropes 
	  if (MonashTune==2 && isMC!=1 && isMC!=0) continue; //only Pythia default
	  if (MonashTune==3 && isMC!=3 && isMC!=0) continue; //only EPOSLHC
	  if (MonashTune==4 && isMC==3) continue; //only Pythia
	  NTypeMC++;
	  cout << "\n" << MCType[isMC] << endl;
	  //Input file 
	  if (isMC==0) {
	    pathin[isMC] = "RatiosXiK0s3Systems_" + SPlotType[ParticleType];
	    if (ChosenRegion>=0) pathin[isMC] += "_" + RegionType[ChosenRegion];
	    if (isdNdEtaTriggered)  pathin[isMC]+= "_isdNdEtaTriggered";
	    if (ParticleType==0 && isWingsCorrectionApplied) pathin[isMC] += "_WingsCorrApplied";
	    if (MaterialBudgetCorr) pathin[isMC] += "_MatBudgetCorrFAST";
	    pathin[isMC] += "_MultCorrSistEval";
	    pathin[isMC] += ".root";
	  }
	  else {
	    pathin[isMC] = "Results_pp13TeVMB_FastMCPrediction"; //temporary solution
	    if (isdNdEtaTriggered)  pathin[isMC]+= "_isdNdEtaTriggered";
	    if (isWingsCorrectionApplied) pathin[isMC] += "_isWingsCorrectionAppliedNew";
	    if (MonashTune==2) pathin[isMC] += "_PythiaMonash";
	    else if (MonashTune==1) pathin[isMC] += "_PythiaRopes";
	    else if (MonashTune==3) pathin[isMC] += "_EPOSLHC";
	    else if (MonashTune==4 || MonashTune==5){
	      if (isMC==1) pathin[isMC] += "_PythiaMonash";
	      else if (isMC==2) pathin[isMC] += "_PythiaRopes";
	      else if (isMC==3) pathin[isMC] += "_EPOSLHC";
	    }
	    pathin[isMC]+= ".root";
	  }

	  cout << "\n\e[35mPathin: " << pathin[isMC] << "\e[39m"<< endl;
	  TFile *filein = new TFile(pathin[isMC], "");
	  if (!filein) {cout << "Input file not available " << endl; return;}

	  //input histo
	  if (isMC ==0) {
	    NameHisto13TeV = SPlotType[ParticleType] + "_" + RegionType[ireg] + "_13TeV_Stat";
	    NameHisto5TeV = SPlotType[ParticleType] + "_" + RegionType[ireg] + "_5TeV_Stat";
	    NameHistoSist13TeV = SPlotType[ParticleType] + "_" + RegionType[ireg] + "_13TeV_Sist";
	    NameHistoSist5TeV = SPlotType[ParticleType] + "_" + RegionType[ireg] + "_5TeV_Sist";
	    NameHistoSistMultUnCorr13TeV = SPlotType[ParticleType] + "_" + RegionType[ireg] + "_13TeV_SistMultUnCorr";
	    NameHistoSistMultUnCorr5TeV = SPlotType[ParticleType] + "_" + RegionType[ireg] + "_5TeV_SistMultUnCorr";
	  }
	  else {
	      NameHisto13TeV = "Yield_"+ SPlotTypeBis[ParticleType] + "_" + RegionTypeBis[ireg] + "_StatErr";
	      NameHistoSist13TeV = "Yield_"+ SPlotTypeBis[ParticleType] + "_" + RegionTypeBis[ireg] + "_SistErr";
	  }

	  //cout << "NameInputHisto " << NameHisto13TeV << endl;
	  //cout << "NameInputHistoSist " << NameHistoSist13TeV << endl;
	  if (Coll==0) histoYield[ireg][Coll][isMC][ParticleType]= (TH1F*) filein->Get(NameHisto13TeV);
	  else 	histoYield[ireg][Coll][isMC][ParticleType]= (TH1F*) filein->Get(NameHisto5TeV);
	  if (!histoYield[ireg][Coll][isMC][ParticleType]) {cout <<"no histo " << endl; return;}

	  if (Coll==0)	histoYieldSist[ireg][Coll][isMC][ParticleType]= (TH1F*) filein->Get(NameHistoSist13TeV);
	  else histoYieldSist[ireg][Coll][isMC][ParticleType]= (TH1F*) filein->Get(NameHistoSist5TeV);
	  if (!histoYieldSist[ireg][Coll][isMC][ParticleType]) {cout <<"no histo " << endl; return;}


	  if (isMC==0){
	    if (Coll==0)	histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]= (TH1F*) filein->Get(NameHistoSistMultUnCorr13TeV);
	    else histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]= (TH1F*) filein->Get(NameHistoSistMultUnCorr5TeV);
	    if (!histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]) {cout <<"no histo " << endl; return;}
	    histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType] = (TH1F*)histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->Clone(NameHistoSistMultUnCorr13TeV + "_Stat");
	    Float_t Err=0;
	    Float_t MinMult=0;
	    Float_t MaxMult=0;
	    for (Int_t b=1; b<= histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetNbinsX(); b++){
	      Err = histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetBinError(b);
	      histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->SetBinError(b, sqrt(pow(Err, 2) + pow(histoYield[ireg][Coll][isMC][ParticleType]->GetBinError(b), 2)));
	      MinMult = histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][0]);
	      MaxMult = histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][7]);
	      if (b==MinMult) {
		MinValueYield[ireg][Coll][isMC][ParticleType] = histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetBinContent(b);
		MinValueError[ireg][Coll][isMC][ParticleType] = histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetBinError(b);
	      }
	      else if (b==MaxMult) {
		MaxValueYield[ireg][Coll][isMC][ParticleType] = histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetBinContent(b);
		MaxValueError[ireg][Coll][isMC][ParticleType] = histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetBinError(b);
	      }
	    }
	    cout << "Ratio Max/Min: " << MaxValueYield[ireg][Coll][isMC][ParticleType]/MinValueYield[ireg][Coll][isMC][ParticleType]  << " +- " <<  sqrt(pow(MaxValueError[ireg][Coll][isMC][ParticleType]/MaxValueYield[ireg][Coll][isMC][ParticleType] , 2) + pow(MinValueError[ireg][Coll][isMC][ParticleType]/ MinValueYield[ireg][Coll][isMC][ParticleType], 2) ) *  MaxValueYield[ireg][Coll][isMC][ParticleType]/MinValueYield[ireg][Coll][isMC][ParticleType]<< endl;
	  }


	  cout << "Region: " << Region[ireg]<< endl;
	  if (Coll ==0 ) cout << "Got histos for  13 TeV (MB + HM) " << endl;
	  else  cout << "Got histos for 5 TeV " << endl;
	  for (Int_t b=1; b<=  histoYield[ireg][Coll][isMC][ParticleType]->GetNbinsX(); b++){
	    if (histoYield[ireg][Coll][isMC][ParticleType]->GetBinContent(b) != 0)	 {
	      cout << "dNdeta (central value of bin) " << histoYield[ireg][Coll][isMC][ParticleType]->GetBinCenter(b) << ": " << histoYield[ireg][Coll][isMC][ParticleType]->GetBinContent(b) << " +- " << histoYield[ireg][Coll][isMC][ParticleType]->GetBinError(b) << " (stat.) +-  " <<  histoYieldSist[ireg][Coll][isMC][ParticleType]->GetBinError(b) << " (syst.) " << endl;
	      if (isMC==0) cout << " fraction of mult. uncorr. uncertainty: " << histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetBinError(b)/histoYieldSist[ireg][Coll][isMC][ParticleType]->GetBinError(b) << endl;
	    }
	  }

	  /* TGRAPH creation */
	  if (isMC==0) nummoltDataorMC = 8;
	  else nummoltDataorMC = nummolt;
	  for (Int_t m=0; m< nummoltDataorMC; m++){
	    //cout << "m " << m << " " << dNdEta[isMC][Coll][m]<< endl;
	    Yields[ireg][Coll][isMC][ParticleType][m] = histoYield[ireg][Coll][isMC][ParticleType]->GetBinContent(histoYield[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m]));
	    YieldsErrors[ireg][Coll][isMC][ParticleType][m] = sqrt(pow(histoYield[ireg][Coll][isMC][ParticleType]->GetBinError(histoYield[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m])), 2) +  pow(histoYieldSist[ireg][Coll][isMC][ParticleType]->GetBinError(histoYield[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m])),2));
	    YieldsErrorsStat[ireg][Coll][isMC][ParticleType][m] = histoYield[ireg][Coll][isMC][ParticleType]->GetBinError(histoYield[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m]));
	    YieldsErrorsSist[ireg][Coll][isMC][ParticleType][m] = histoYieldSist[ireg][Coll][isMC][ParticleType]->GetBinError(histoYieldSist[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m]));
	    if (isMC==0) YieldsErrorsSistMultUnCorr[ireg][Coll][isMC][ParticleType][m] = histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->GetBinError(histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m]));
	    if (ParticleType==1) {
	      if (Coll==0) cout << "ratio K0s/ Xi, mult " << m << ": " << Yields[ireg][Coll][isMC][0][m]/Yields[ireg][Coll][isMC][1][m] << endl;
	      Yields[ireg][Coll][isMC][ParticleType][m] *= ScalingXi;
	      YieldsErrors[ireg][Coll][isMC][ParticleType][m] *= ScalingXi;
	      YieldsErrorsStat[ireg][Coll][isMC][ParticleType][m] *= ScalingXi;
	      YieldsErrorsSist[ireg][Coll][isMC][ParticleType][m] *= ScalingXi;
	      if (isMC==0) YieldsErrorsSistMultUnCorr[ireg][Coll][isMC][ParticleType][m] *= ScalingXi;
	    }
	  }
	  if (isMC==0) ghistoYield[ireg][Coll][isMC][ParticleType] = new TGraphAsymmErrors(nummoltDataorMC,dNdEta[isMC][Coll],Yields[ireg][Coll][isMC][ParticleType],0, 0,YieldsErrorsStat[ireg][Coll][isMC][ParticleType], YieldsErrorsStat[ireg][Coll][isMC][ParticleType]);
	  else ghistoYield[ireg][Coll][isMC][ParticleType] = new TGraphAsymmErrors(nummoltDataorMC,dNdEta[isMC][Coll],Yields[ireg][Coll][isMC][ParticleType],dNdEtaErrorL[isMC][Coll],dNdEtaErrorR[isMC][Coll],YieldsErrors[ireg][Coll][isMC][ParticleType], YieldsErrors[ireg][Coll][isMC][ParticleType]);
	  ghistoYield[ireg][Coll][isMC][ParticleType]->SetName(Form("ghistoYield_reg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	  ghistoYieldSist[ireg][Coll][isMC][ParticleType] = new TGraphAsymmErrors(nummoltDataorMC,dNdEta[isMC][Coll],Yields[ireg][Coll][isMC][ParticleType],dNdEtaErrorL[isMC][Coll],dNdEtaErrorR[isMC][Coll], YieldsErrorsSist[ireg][Coll][isMC][ParticleType], YieldsErrorsSist[ireg][Coll][isMC][ParticleType]);
	  ghistoYieldSist[ireg][Coll][isMC][ParticleType]->SetName(Form("ghistoYieldSist_reg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	  if (isMC==0){
	    ghistoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType] = new TGraphAsymmErrors(nummoltDataorMC,dNdEta[isMC][Coll],Yields[ireg][Coll][isMC][ParticleType],dNdEtaErrorL[isMC][Coll],dNdEtaErrorR[isMC][Coll], YieldsErrorsSistMultUnCorr[ireg][Coll][isMC][ParticleType], YieldsErrorsSistMultUnCorr[ireg][Coll][isMC][ParticleType]);
	    ghistoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->SetName(Form("ghistoYieldSistMultUnCorr_reg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	    StyleTGraphErrors(ghistoYield[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], MarkerType[ireg][Coll], MarkerSize[ireg], LineStyle[isMC]);
	    StyleTGraphErrors(ghistoYieldSist[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], MarkerType[ireg][Coll], MarkerSize[ireg], LineStyle[isMC]);
	    StyleTGraphErrors(ghistoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], MarkerType[ireg][Coll], MarkerSize[ireg], LineStyle[isMC]);
	  }
	  else {
	    if (ChosenRegion==0)  StyleTGraphErrors(ghistoYield[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], 1, 1, 1);
	    //	  if (ChosenRegion==0)  StyleTGraphErrors(ghistoYield[ireg][Coll][isMC][ParticleType], ColorDiff[ireg][Coll], 1, 1, LineStyle[isMC]);
	    else StyleTGraphErrors(ghistoYield[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], 1, 1, LineStyle[isMC]);
	  }

	  //Splines to TGraphs
	  splineY[ireg][Coll][isMC][ParticleType] = new TSpline3(Form("Spline_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), ghistoYield[ireg][Coll][isMC][ParticleType]);
	  sp3= (TSpline3*) splineY[ireg][Coll][isMC][ParticleType] ->Clone(Form("sp3Spline_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	  fsplineY[ireg][Coll][isMC][ParticleType]= new TF1(Form("fSpline_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), spline, 0, 40);
	  fsplineY[ireg][Coll][isMC][ParticleType]->SetLineColor(ColorDiff[ParticleType][Coll]);
	  fsplineY[ireg][Coll][isMC][ParticleType]->SetLineStyle(1);

	  //Fit to yields
	  fitY[ireg][Coll][isMC][ParticleType]= new TF1(Form("fitY_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), "pol0", 0, 40);
	  fitY[ireg][Coll][isMC][ParticleType]->SetLineColor(ColorDiff[ParticleType][Coll]);
	  fitY[ireg][Coll][isMC][ParticleType]->SetLineStyle(7);
	  cout << "Fitting the yields with statistical uncertainty" << endl;
	  histoYield[ireg][Coll][isMC][ParticleType]->Fit(fitY[ireg][Coll][isMC][ParticleType], "R0");

	  if (isMC==0){
	    cout << "Fitting the yields with statistical + syst. uncorrelated uncertainty" << endl;
	    histoYieldStatSistMultUnCorr[ireg][Coll][isMC][ParticleType]->Fit(fitY[ireg][Coll][isMC][ParticleType], "R0");
	  }
	
	  //create histo in order to perform division of TF1s
	  hMCDataRatio[ireg][Coll][isMC][ParticleType]= new TH1F (Form("hMCDataRatio_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), Form("fitY_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC), 1000, 0, 40);
	  //	if (ireg!=0)	hMCDataRatio[ireg][Coll][isMC][ParticleType]->Add(fitY[ireg][Coll][isMC][ParticleType]);
	  if (ireg!=0)	hMCDataRatio[ireg][Coll][isMC][ParticleType]->Add(fsplineY[ireg][Coll][isMC][ParticleType]);
	  else hMCDataRatio[ireg][Coll][isMC][ParticleType]->Add(fsplineY[ireg][Coll][isMC][ParticleType]);
      
	  canvas->cd();
	  gStyle->SetLegendBorderSize(0);
	  gStyle->SetLegendFillColor(0);
	  gStyle->SetLegendFont(42);
	  pad1->Draw();
	  pad1->cd();

	  ghistoYieldLine[ireg][Coll][isMC][ParticleType]=(TGraphAsymmErrors*) 	  ghistoYield[ireg][Coll][isMC][ParticleType]->Clone(Form("ghistoYieldLine_%i_%i_%i_%i", ireg, Coll, isMC, ParticleType));
	  ghistoYieldLine[ireg][Coll][isMC][ParticleType]->SetLineColor(ColorDiffLine[ParticleType][Coll]);
	  ghistoYieldLine[ireg][Coll][isMC][ParticleType]->SetLineStyle(LineStyle[isMC]);
	  if (isMC ==1 || isMC==2 || isMC==3){
	    ghistoYield[ireg][Coll][isMC][ParticleType]->SetFillColor(ColorDiff[ParticleType][Coll]);
	    if (isMC==1){
	      ghistoYield[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]);
	      ghistoYield[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(ColorDiff[ParticleType][Coll],AlphaColor[isMC]);
	    }
	    else if (isMC==2) {
	      ghistoYield[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]);
	      ghistoYield[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(ColorDiff[ParticleType][Coll],AlphaColor[isMC]);
	    }
	    else {
	      ghistoYield[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]);
	      ghistoYield[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(ColorDiff[ParticleType][Coll],AlphaColor[isMC]);
	    }
	  }

	  //	  StyleTGraphErrors(ghistoYield[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], 1, 1, LineStyle[isMC]);

	  if (isMC==1)	  ColorDiff[ParticleType][Coll] =  ColorDiffMC[ParticleType][Coll];
	  else if (isMC==2) 	  ColorDiff[ParticleType][Coll] =  ColorDiffMCRopes[ParticleType][Coll];

	  StyleHisto(histoYield[ireg][Coll][isMC][ParticleType], Low, Up, ColorDiff[ParticleType][Coll], MarkerType[ireg][Coll], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	  StyleHisto(histoYieldSist[ireg][Coll][isMC][ParticleType], Low, Up, ColorDiff[ParticleType][Coll], MarkerType[ireg][Coll], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	  if (isMC==0) 	StyleHisto(histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType], Low, Up, ColorDiff[ParticleType][Coll], MarkerType[ireg][Coll], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);

	  //legend
	  if (NLoopRegion==0 && Coll==0 && isMC!=0){
	    ghistoYieldGrey[ireg][Coll][isMC][ParticleType] = (TGraphAsymmErrors*)ghistoYield[ireg][Coll][isMC][ParticleType]->Clone(Form("ghistoYieldGrey_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	    ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetLineColor(kGray +3);
	    ghistoYieldRed[ireg][Coll][isMC][ParticleType] = (TGraphAsymmErrors*) ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->Clone(Form("ghistoYieldRed_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	    ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetLineColor(ColorDiff[ParticleType][Coll]);
	    if (isMC==1) {
	      //ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillStyle(3001); //TEST
	      ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]); 
	      ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(ColorDiff[ParticleType][Coll],AlphaColor[isMC]);
	    }
	    else if (isMC==2) {
	      //ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillStyle(3008); //TEST
	      ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]); 
	      ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(ColorDiff[ParticleType][Coll],AlphaColor[isMC]);
	    }
	    else {
	      //ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillStyle(3002); //TEST
	      ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]); 
	      ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(ColorDiff[ParticleType][Coll],AlphaColor[isMC]);
	    }
	    ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetLineStyle(1);
	    ghistoYieldRed[ireg][Coll][isMC][ParticleType]->SetLineWidth(0);

	    if (/*PlotType==0 ||*/ ChosenRegion==0){
	      //if (PlotType==0 || ChosenRegion==0){
	      if (isMC==1) {
		//ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillStyle(3001); //TEST
		ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]); 
		ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(kGray+3,AlphaColor[isMC]);
	      }
	      else if (isMC==2){
		//ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillStyle(3008); //TEST
		ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]); 
		ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(kGray+3,AlphaColor[isMC]);
	      }
	      else {
		//ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillStyle(3002); //TEST
		ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillStyle(FillStyle[isMC]); 
		ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(kGray+3,AlphaColor[isMC]);
	      }
	      ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->SetLineWidth(0);
	      if (ParticleType==0)  legendMCTypes->AddEntry(ghistoYieldGrey[ireg][Coll][isMC][ParticleType], MCType[isMC], "f");
	    }
	    else {
	      if (ParticleType==0) legendMCTypes->AddEntry(ghistoYieldGrey[ireg][Coll][isMC][ParticleType], MCType[isMC], "l");
	    }
	  }

	  if (NLoopRegion==0 && NTypeMC==0){
	    fHistYieldStatGrey[Coll]= (TH1F*) histoYield[ireg][Coll][isMC][ParticleType]->Clone(Form("fHistYieldStatBlack_%i", Coll));
	    fHistYieldStatGrey[Coll]->SetLineColor(ColorEnergy[Coll]);
	    fHistYieldStatGrey[Coll]->SetMarkerColor(ColorEnergy[Coll]);
	    fHistYieldSistGrey[Coll]= (TH1F*) histoYieldSist[ireg][Coll][isMC][ParticleType]->Clone(Form("fHistYieldSistBlack_%i", Coll));
	    fHistYieldSistGrey[Coll]->SetLineColor(ColorEnergy[Coll]);
	    fHistYieldSistGrey[Coll]->SetMarkerColor(ColorEnergy[Coll]);
	    fHistYieldSistMultUnCorrGrey[Coll]= (TH1F*) histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->Clone(Form("fHistYieldSistMultUnCorrBlack_%i", Coll));
	    fHistYieldSistMultUnCorrGrey[Coll]->SetLineColor(ColorEnergy[Coll]);
	    fHistYieldSistMultUnCorrGrey[Coll]->SetMarkerColor(ColorEnergy[Coll]);

	    if (Coll==0){
	      /*
		TLegendEntry* L1 = legendEnergyBox1->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 13 TeV", "pe");
		L1->SetTextAlign(32);
		L1->SetTextSize(0.042);
	      */
	      //fHistYieldStatGrey[Coll]->SetMarkerStyle(53);
	      fHistYieldStatGrey[Coll]->SetMarkerStyle(28);
	      fHistYieldStatGrey[Coll]->SetMarkerSize(2.5);
	      if (ParticleType==0){
		legendEnergyBox1->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 13 TeV", "p");
		legendEnergyBox1->SetTextSize(0.042);
		legendEnergyBox2->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 13 TeV", "p");
		legendEnergyBox2->SetTextSize(0.042);
	      }

	      //legendStatBox->AddEntry(fHistYieldStatGrey[Coll], "stat.", "le");
	      //legendStatBox->AddEntry(fHistYieldSistGrey[Coll], "syst.", "ef");
	      //legendStatBox->AddEntry(fHistYieldSistMultUnCorrGrey[Coll], "syst. uncorr.", "ef");

	    }
	    else {
	      /*
		TLegendEntry* L2 = legendEnergyBox2->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
		L2->SetTextAlign(32);
		L2->SetTextSize(0.042);
	      */
	      //fHistYieldStatGrey[Coll]->SetMarkerStyle(20);
	      fHistYieldStatGrey[Coll]->SetMarkerStyle(34);
	      fHistYieldStatGrey[Coll]->SetMarkerSize(2.5);
	      if (ParticleType==0){
		legendEnergyBox1->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 5.02 TeV", "p");
		legendEnergyBox1->SetTextSize(0.042);
		legendEnergyBox2->AddEntry(fHistYieldStatGrey[Coll], "pp, #sqrt{#it{s}} = 5.02 TeV", "p");
		legendEnergyBox2->SetTextSize(0.042);
	      }
	    }
	  }

	  if (NTypeMC == 0){
	    cout << "Coll " << Coll << endl;
	    cout << "MC: " << isMC << endl;
	    cout << "Region: " << ireg << endl;
	    legendEnergyBoxHelloLower->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "", "p");
	    if (Coll==0){
	      if (ParticleType==0) LegendParticle->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "h#minusK_{S}^{0} correlation,  |#it{#eta}^{K^{0}_{S}}| < 0.8", "p");
	      else LegendParticle->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "h#minus#Xi correlation #times 29, |#it{#eta}^{#Xi}| < 0.8", "p");
	    }
	    if (ParticleType == 0){
	      legendEnergyBoxColor1->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "", "p");
	      legendEnergyBoxColorBis1->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "", "p");
	    }
	    else {
	      if (Coll==0) {
		legendEnergyBoxColor3->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "pp, #sqrt{#it{s}} = 13 TeV", "p");
		legendEnergyBoxColorBis3->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "pp, #sqrt{#it{s}} = 13 TeV", "p");
	      }
	      else {
		legendEnergyBoxColor3->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "pp, #sqrt{#it{s}} = 5.02 TeV", "p");
		legendEnergyBoxColorBis3->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "pp, #sqrt{#it{s}} = 5.02 TeV", "p");
	      }
	    }
	  }

	  if (Coll==0 && NTypeMC ==0){

	    if (ireg==0) lReAll1Bis[ireg] = legendRegionJet->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], sRegion[ireg], "");
	    else if (ireg==1)lReAll1Bis[ireg] = legendRegionBulk->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], sRegion[ireg], "");
	    else lReAll1Bis[ireg] = legendRegionFull->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], sRegion[ireg], "");

	    lReAll1Bis[ireg]->SetTextSize(0.045);
	    //      lReAll1Bis[ireg]->SetTextAlign(32);

	    if (ireg==0) lReAll2Bis[ireg]=      legendRegionJet->AddEntry("", sRegion1[ireg], "");
	    else    if (ireg==1) lReAll2Bis[ireg]=      legendRegionBulk->AddEntry("", sRegion1[ireg], "");
	    else    lReAll2Bis[ireg]=      legendRegionFull->AddEntry("", sRegion1[ireg], "");
	    lReAll2Bis[ireg]->SetTextSize(0.035);

	    lReAll1[ireg] = legendRegionAll->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], sRegion[ireg], "p");
	    lReAll1[ireg]->SetTextSize(0.048);
	    lReAll2[ireg]= legendRegionAll->AddEntry("", sRegion1[ireg], "");
	    lReAll2[ireg]->SetTextSize(0.038);

	    if (NLoopRegion==0){
	      fHistYieldStatBlack= (TH1F*)      histoYield[ireg][Coll][isMC][ParticleType]->Clone("fHistYieldStatBlack");
	      fHistYieldSistBlack= (TH1F*)      histoYieldSist[ireg][Coll][isMC][ParticleType]->Clone("fHistYieldSistBlack");
	      fHistYieldSistMultUnCorrBlack= (TH1F*)      histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->Clone("fHistYieldSistMultUnCorrBlack");
	      fHistYieldStatBlack->SetLineColor(1);
	      fHistYieldStatBlack->SetMarkerColor(1);
	      fHistYieldStatBlack->SetMarkerStyle(20);
	      fHistYieldSistBlack->SetLineColor(1);
	      fHistYieldSistBlack->SetMarkerColor(1);
	      fHistYieldSistBlack->SetMarkerStyle(20);
	      fHistYieldSistMultUnCorrBlack->SetLineColor(1);
	      fHistYieldSistMultUnCorrBlack->SetMarkerColor(1);
	      fHistYieldSistMultUnCorrBlack->SetFillColor(1);
	      fHistYieldSistMultUnCorrBlack->SetFillStyle(3001);
	      fHistYieldSistMultUnCorrBlack->SetMarkerStyle(20);
	      if (ParticleType==0){
		legendStatBoxBis->AddEntry(fHistYieldSistMultUnCorrBlack, "syst. uncorr.", "ef");
		legendStatBoxBis->AddEntry(fHistYieldSistBlack, "syst.", "ef");
		legendStatBoxBis->AddEntry(fHistYieldStatBlack, "stat.", "le");
		legendStatBox->AddEntry(fHistYieldStatBlack, "stat.", "le");
		legendStatBox->AddEntry(fHistYieldSistBlack, "syst.", "ef");
		legendStatBox->AddEntry(fHistYieldSistMultUnCorrBlack, "syst. uncorr.", "ef");
	      }
	    }
	  }

	  if (Coll==0 && NTypeMC==0){
	    if (ParticleType==0){
	    legendStatBoxColor->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "stat.", "el");
	    legendStatBoxColor->AddEntry(histoYieldSist[ireg][Coll][isMC][ParticleType], "syst.", "ef");
	    histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->SetFillStyle(3001);
	    histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->SetFillColor(ColorDiff[ParticleType][Coll]);
	    legendStatBoxColor->AddEntry(histoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType], "syst. uncorr.", "ef");

	    TLegendEntry * lRe1 = legendOneRegion->AddEntry("", sRegion[ireg], "");
	    lRe1->SetTextSize(0.06);
	    lRe1->SetTextAlign(32);
	    TLegendEntry *leR2=      legendOneRegion->AddEntry("", sRegion1[ireg], "");
	    leR2->SetTextSize(0.04);
	    leR2->SetTextAlign(32);
	    }
	  }
	  //	if (Coll==0)	LegendColor->AddEntry("", NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");

	  if (NTypeMC==0 && ParticleType==0){
	    if (Coll==1)	legendEnergyBoxColor->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
	    else 	legendEnergyBoxColor->AddEntry(histoYield[ireg][Coll][isMC][ParticleType], "pp, #sqrt{#it{s}} = 13 TeV", "pe");
	    legendEnergyBoxColor->SetTextSize(0.042);
	  }

	  if (isMC ==1 || isMC==2 || isMC==3){
	    ghistoYield[ireg][Coll][isMC][ParticleType]->DrawClone("same e3");
	    ghistoYieldLine[ireg][Coll][isMC][ParticleType]->DrawClone("same lX");
	  }
	  else {
	    histoYieldDummy->DrawClone("same");
	    ghistoYieldSist[ireg][Coll][isMC][ParticleType]->SetFillStyle(0);
	    ghistoYieldSist[ireg][Coll][isMC][ParticleType]->SetFillColor(ColorDiff[ParticleType][Coll]);
	    ghistoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->SetFillStyle(3001);
	    ghistoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->SetFillColor(ColorDiff[ParticleType][Coll]);
	    if (ParticleType==1){
	      ghistoYieldSist[ireg][Coll][isMC][ParticleType]->DrawClone("same p2");
	      ghistoYieldSistMultUnCorr[ireg][Coll][isMC][ParticleType]->DrawClone("same p2");
	      ghistoYield[ireg][Coll][isMC][ParticleType]->DrawClone("same e");
	      ghistoYieldSist[ireg][Coll][isMC][0]->DrawClone("same p2");
	      ghistoYieldSistMultUnCorr[ireg][Coll][isMC][0]->DrawClone("same p2");
	      ghistoYield[ireg][Coll][isMC][0]->DrawClone("same e");
	    }
	  }

	  if (isMC == MCindexForPlotting){
	    LegendYields->Draw("");
	    LegendParticle->Draw("");
	    //LegendColor->Draw("");
	    //legendStatBoxColor->Draw("");
	    legendStatBox->Draw("");
	    legendOneRegion->Draw("");
	    legendEnergyBoxColor1->Draw("");
	    legendEnergyBoxColor3->Draw("");
	    legendMCTypes->Draw("");
	  }

	  //MC / DATA Ratios
	  if (isMC >0){
	    //Ratios between splines (jet) or fits (OOJ, full) -> to be used to fill a TGraph
	    hMCDataRatio[ireg][Coll][isMC][ParticleType]->Divide(hMCDataRatio[ireg][Coll][0][ParticleType]);
	  }
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

  TH1F* hdummy = new TH1F ("hdummy", "hdummy",  1000, 0, UpperValueX); 
  StyleHisto(hdummy, LowRatio, UpRatio, 1, 1, titleX, titleYMCToData, "" , 1, 0, UpperValueX, xOffset, yOffset, 1);
  SetFont(hdummy);
  SetHistoTextSize(hdummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hdummy, tickXL, tickYL);
  hdummy->GetYaxis()->SetDecimals(kTRUE);
  hdummy->GetYaxis()->SetNdivisions(505, "false");
  hdummy->GetYaxis()->CenterTitle(true);

  TF1 * lineAt1 = new TF1 ("lineAt1", "pol0", 0, UpperValueX);
  lineAt1->SetParameter(0,1);
  lineAt1->SetLineColor(1);
  lineAt1->SetLineStyle(2);

  cout << "\n\n//DEFINE RATIO data/data (=1) with error given by the sum in quadrature of stat and syst uncertainty " << endl;
  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    for (Int_t ParticleType=0; ParticleType < numParticles; ParticleType++){ //K0s and Xi
      cout << "Particle Type: " << ParticleType << endl;
      for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
	if (Coll==1) continue; //for the time being, only 13 TeV
	for (Int_t isMC = 0; isMC <= NumberOfMCs; isMC++){
	  if (isMC!=0) continue;
	  for (Int_t m=0; m< 8; m++){
	    YieldRatiosDATA[ireg][Coll][isMC][ParticleType][m] = 1;
	    YieldRatiosErrorsDATA[ireg][Coll][isMC][ParticleType][m] = sqrt(pow(histoYield[ireg][Coll][isMC][ParticleType]->GetBinError(histoYield[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m])), 2)+ pow(histoYieldSist[ireg][Coll][isMC][ParticleType]->GetBinError(histoYield[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m])), 2)) / histoYield[ireg][Coll][isMC][ParticleType]->GetBinContent(histoYield[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m]));
	    cout << YieldRatiosDATA[ireg][Coll][isMC][ParticleType][m] << " +- " <<  YieldRatiosErrorsDATA[ireg][Coll][isMC][ParticleType][m] << endl;
	  }
	  //cout << "Define TGraph " << endl;
	  ghistoYieldRatioDATA[ireg][Coll][isMC][ParticleType] = new TGraphAsymmErrors(8,dNdEta[isMC][Coll],YieldRatiosDATA[ireg][Coll][isMC][ParticleType],0,0,YieldRatiosErrorsDATA[ireg][Coll][isMC][ParticleType],YieldRatiosErrorsDATA[ireg][Coll][isMC][ParticleType]);
	  ghistoYieldRatioDATA[ireg][Coll][isMC][ParticleType]->SetName(Form("ghistoYieldRatioDATA_reg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	  StyleTGraphErrors(ghistoYieldRatioDATA[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], 1, 1, LineStyle[isMC]);
	}
      }
    }
  }

  NLoopRegion=-1;
  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    NLoopRegion++;
    for (Int_t ParticleType=0; ParticleType < numParticles; ParticleType++){ //K0s and Xi
      for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
	if (Coll==1) continue; //for the time being, only 13 TeV
	for (Int_t isMC = 0; isMC <= NumberOfMCs; isMC++){

	  if (MonashTune==1 && isMC!=2) continue; //only Ropes 
	  if (MonashTune==2 && isMC!=1) continue; //only Pythia default
	  if (MonashTune==3 && isMC!=3) continue; //only EPOSLHC
	  if (MonashTune==4 && isMC==3) continue; //only Pythia

	  StyleHisto(hMCDataRatio[ireg][Coll][isMC][ParticleType], LowRatio, UpRatio, ColorDiff[ParticleType][Coll], 1, titleX, titleYMCToData, "" , 1, 0, 40, xOffset, yOffset, 1);

	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetYaxis()->SetTitleSize(0.07);
	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetYaxis()->SetTitleOffset(0.6);

	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetYaxis()->SetLabelSize(0.10);
	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetYaxis()->SetLabelOffset(0.009);

	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetXaxis()->SetTitleSize(0.10);
	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetXaxis()->SetTitleOffset(1.3);

	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetXaxis()->SetLabelSize(0.10);
	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetXaxis()->SetLabelOffset(0.02);

	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetYaxis()->SetNdivisions(305);
	  hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetYaxis()->SetTickLength(0.02);

	  //hMCDataRatio[ireg][Coll][isMC][ParticleType]->Draw("same");

	  hdummy->DrawClone("same");

	  for (Int_t m=0; m< nummolt; m++){
	    YieldRatios[ireg][Coll][isMC][ParticleType][m] = hMCDataRatio[ireg][Coll][isMC][ParticleType]->GetBinContent(hMCDataRatio[ireg][Coll][isMC][ParticleType]->FindBin(dNdEta[isMC][Coll][m]));
	    YieldRatiosErrors[ireg][Coll][isMC][ParticleType][m] = 0; 
	    //	  cout << "m " << dNdEta[isMC][Coll][m] << " " << YieldRatios[ireg][Coll][isMC][ParticleType][m] << endl;
	  }
	  ghistoYieldRatio[ireg][Coll][isMC][ParticleType] = new TGraphAsymmErrors(nummolt,dNdEta[isMC][Coll],YieldRatios[ireg][Coll][isMC][ParticleType], dNdEtaErrorL[isMC][Coll],  dNdEtaErrorR[isMC][Coll], YieldRatiosErrors[ireg][Coll][isMC][ParticleType], YieldRatiosErrors[ireg][Coll][isMC][ParticleType]);
	  ghistoYieldRatio[ireg][Coll][isMC][ParticleType]->SetName(Form("ghistoYieldRatio_reg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	  StyleTGraphErrors(ghistoYieldRatio[ireg][Coll][isMC][ParticleType], ColorDiff[ParticleType][Coll], 1, 1, LineStyle[isMC]);

	  if (isMC!=0){
	    ghistoYieldRatio[ireg][Coll][isMC][ParticleType]->Draw("same");
	  }
	  else {
	    //ghistoYieldRatioDATA[ireg][Coll][isMC][ParticleType]->SetFillColorAlpha(ColorDiff[ParticleType][Coll],0.3);
	    ghistoYieldRatioDATA[ireg][Coll][isMC][ParticleType]->SetFillStyle(1001);
	    ghistoYieldRatioDATA[ireg][Coll][isMC][ParticleType]->SetFillColor(ColorDiffPale[ParticleType][Coll]);
	    if (ParticleType==1) {
	      ghistoYieldRatioDATA[ireg][Coll][isMC][1]->DrawClone("same 3");
	      ghistoYieldRatioDATA[ireg][Coll][isMC][0]->DrawClone("same 3"); //drawn afterwards since it is smaller
	    }
	  }
	  lineAt1->Draw("same");

	  if (Coll==0 && isMC!=0 && ParticleType==0){
	      ghistoYieldGreyBis[ireg][Coll][isMC][ParticleType] = (TGraphAsymmErrors*)ghistoYieldGrey[ireg][Coll][isMC][ParticleType]->Clone(Form("ghistoYieldGreyBis_ireg%i_Coll%i_isMC%i", ireg, Coll, isMC));
	      ghistoYieldGreyBis[ireg][Coll][isMC][ParticleType]->SetLineStyle(LineStyle[isMC]);
	      ghistoYieldGreyBis[ireg][Coll][isMC][ParticleType]->SetLineWidth(2);
	      ghistoYieldGreyBis[ireg][Coll][isMC][ParticleType]->SetLineColor(kGray + 3);
	      smalllegendMCTypes->AddEntry(ghistoYieldGreyBis[ireg][Coll][isMC][ParticleType], MCType[isMC], "l");
	  }
	  smalllegendMCTypes->Draw("");
	}
      }
    }
  }

  TString OutputDir = "FiguresForPaper/";
  canvas->SaveAs(OutputDir+"CompareDatavsMC_K0sXiJet"+"_" + RegionType[ChosenRegion]+ MCTypeBis[MonashTune-1]+".pdf");
  canvas->SaveAs(OutputDir+"CompareDatavsMC_K0sXiJet"+"_" + RegionType[ChosenRegion]+ MCTypeBis[MonashTune-1]+".eps");
  canvas->SaveAs(OutputDir+"CompareDatavsMC_K0sXiJet"+"_" + RegionType[ChosenRegion]+ MCTypeBis[MonashTune-1]+".png");
  
  TFile *fileout = new TFile("CompareDatavsMC.root", "RECREATE");
  fileout->WriteTObject(canvas);
  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    for (Int_t ParticleType=0; ParticleType < numParticles; ParticleType++){ //K0s and Xi
      for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
	if (Coll==1) continue; //for the time being, only 13 TeV
	for (Int_t isMC = 0; isMC <= NumberOfMCs; isMC++){
	  if (isMC==0) continue;
	  if (MonashTune==1 && isMC!=2) continue; //only Ropes 
	  if (MonashTune==2 && isMC!=1) continue; //only Pythia default
	  if (MonashTune==3 && isMC!=3) continue; //only EPOSLHC
	  ghistoYieldRatio[ireg][Coll][isMC][ParticleType]->Write();
	}
      }
    }
  }
  fileout->Close();

  cout <<"\nI created the file: CompareDatavsMC.root" << endl;
}
