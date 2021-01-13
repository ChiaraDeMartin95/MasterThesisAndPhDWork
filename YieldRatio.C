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

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(1.5);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->SetTitle(title);
}
void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); //1.2
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void ErrRatioCorr(TH1F* hNum, TH1F* hDenom, TH1F* hRatio, Bool_t FullCorr){
  //FullCorr == 1 means ro = 1;
  //FullCorr == 0 means ro = 
  Float_t Err1=0;
  Float_t Err2=0;
  Float_t ErrC=0;
  Float_t Err=0;
  for (Int_t b=1; b<=hNum->GetNbinsX();b++){
    if (hNum->GetBinContent(b)==0 ||hDenom->GetBinContent(b)==0){
      hRatio->SetBinError(b,0);
      continue;
    }
    Err1=pow(hNum->GetBinError(b)/hNum->GetBinContent(b),2);
    Err2=pow(hDenom->GetBinError(b)/hDenom->GetBinContent(b),2);
    if (FullCorr==0){
      ErrC=pow(hDenom->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    else {
      ErrC=hDenom->GetBinError(b) * hNum->GetBinError(b)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    Err=sqrt(Err1+Err2-ErrC);
    hRatio->SetBinError(b,Err*hRatio->GetBinContent(b));
  }
  //  return hRatio;                                                                                  
}

void YieldRatio(Bool_t isNormCorr=1, Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=0,   Int_t israp=0,TString yearK0s="1617_hK0s", TString yearK0sMC = "1617MC_hK0s", TString yearXi = "Run2DataRed_MECorr_hXi", TString yearXiMC = "AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"*//*"Run2DataRed_hXi"/"2016kehjl_hK0s"*/,  TString Path1 =""/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/,    TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=0, Bool_t ZeroYieldLowPt=0){

  TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
  TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
  TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
  TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
  TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
  TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 0.85", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 2.0", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};

  gStyle->SetOptStat(0);

  TString hhCorr[2]={"", "_hhCorr"};

  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9
  const Int_t numPtTrigger=1;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;

  const Int_t numsysDPhi =3;

  //latexnames
  //  TLatex LatexTitle="ALICE Preliminary";
  Int_t PtV0Min = 1; //0 el

  cout << "ok " << endl;
  TString DataOrMC[2]={"Data", "MC"};

  cout << "ok " << endl;
  const Int_t numtipo=10;
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Xi"};
  TString Stipo[numtipo]={"K^{0}_{S}", "#Xi"};
  //  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  cout << "ok " << endl;
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltBis[nummolt+1]={"0#minus5", "5#minus10", "10#minus30", "30#minus50", "50#minus100", "0#minus100"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};

  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};

  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};

  for (Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Nmolt[m] = Nmolt0[m];
      Smolt[m] = Smolt0[m];
      SmoltLegend[m] = SmoltLegend0[m];
    }
    else     if (MultBinning==1){
      Nmolt[m] = Nmolt1[m];
      Smolt[m] = Smolt1[m];
      SmoltLegend[m] = SmoltLegend1[m];
    }
    else     if (MultBinning==2){
      Nmolt[m] = Nmolt2[m];
      Smolt[m] = Smolt2[m];
      SmoltLegend[m] = SmoltLegend2[m];
    }

  }


  TString SPtV01[numPtV0]={"0-0.1","0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0, 0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0-0.1", "0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0,0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};
  Int_t numPtV0Max=numPtV0;
  TString SErrorSpectrum[3]={"stat.","syst. uncorr.","syst. corr."};
  TString SSystJet[2]={"BC [-1.0, 1.0]", "BC [-1.2, 1.2]"};
  TString SSystBulk[3]={"BC [1.0, 2.0]", "BC [2.0, 4.28]", "BC [1.0, 4.28]"};

  //  Int_t Marker[numSyst]={20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  //  Int_t MarkerBetter[numSyst]={1,22, 32, 1, 29, 1,   3,  34, 33, 1, 20, 21, 22};
  //  Int_t Color[numSyst]= {1,  2,  3,  4,  5,  6,  7,  7,  4, 10,  6,  1,  2};
  //  Int_t Color[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 628,868,1, 909, 1,801,418,860,  1, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};
  //  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Ymarker = 20;
  Int_t YmarkerRegion[3] = {20, 21, 33};
  Int_t YmarkerRegionBis[3] = {24, 25, 27};
  Float_t YoffsetSpectra[2] = {1.6, 1.65};
  Int_t ColorMult[nummolt+1] ={628,801,418,867,601,1};
  //Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Float_t FracTrig[nummolt+1] ={0.27, 0.18, 0.11, 0.05, 0.013,0.06};
  Float_t size[nummolt+1] ={2, 2, 2.5, 2.5, 2.5, 2};
  Float_t MarkerSize[nummolt+1] ={2, 2, 3};
  Int_t MarkerMult[nummolt+1] ={20,21,33,34,29,25};
  Int_t Color[3] ={628, 418, 600};
  Int_t ColorType[2] ={628,881};

  TString stringout;

  stringout = Dir+"/DATA"+year0+"/YieldRatio" +hhCorr[ishhCorr];
  stringout +=Srap[israp];
  stringout +=SSkipAssoc[SkipAssoc];
  stringout+=   Form("_PtMin%.1f", PtTrigMin);
  if (ZeroYieldLowPt)  stringout += "_ZeroYLowPtJet";
  if (isNormCorr) stringout+="_isNormCorr";
  stringout += "_.root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  //here I start the newly  created part!!

  TString  TitleYPtRatio = "N_{#Xi}/N_{K^{0}_{S}}";
  TString  TitleYPtDRatio = "N_{#Xi}/N_{K^{0}_{S}}_{mult} / N_{#Xi}/N_{K^{0}_{S}}_{0-100%}_{0-100%}";
  TString titleX=  "#it{p}_{T} (GeV/#it{c})";
  TString titleY=  "1/#Delta#it{#eta} #Delta#it{#varphi} 1/#it{N}_{trigg} d#it{N}/d#it{p}_{T} [(GeV/#it{c})^{-1}]";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 
  TString TitleYYieldRatio="#it{N}_{#Xi}/#it{N}_{K^{0}_{S}}";
  TString TitleYieldRatio="#Xi/K^{0}_{S} yield ratio vs multiplicity";
  TString titleYieldX="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
  TString titleYieldY="N/N_{trigg} 1/#Delta#eta #Delta#it{#varphi}";
  // TString titleYieldY="N/N_{Evt} 1/#Delta#eta #Delta#it{#varphi}";
  TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}"};
  //  TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{Evt} 1/#Delta#it{#eta} #Delta#it{#varphi}", "#it{N}_{#Xi}/#it{N}_{Evt} 1/#Delta#it{#eta} #Delta#it{#varphi}"};
  TString titleYield[2]={"K^{0}_{S}", "#Xi"};
  TString titleYieldYRelErr = {"Relative uncertainty"};
  TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

  TLegend *legendPhi = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legendErrorAll;
  TLegend *legendError = new TLegend(0.6, 0.7, 0.9, 0.9);
  TLegend *legendError2 = new TLegend(0.6, 0.7, 0.9, 0.9);

  TCanvas* canvasYield[2];
  TCanvas* canvasYieldFinal[2];
  TCanvas* canvasPtvsMult[2];
  TCanvas* canvasPtvsMultMethodComp[2];
  TCanvas*  canvasYieldSeparate[3][2];
  TCanvas*  canvasPtSpectra[3][2];
  TCanvas*  canvasPtSpectraMultRatio[3][2];
  TCanvas*  canvasYieldErrSeparate[3][2];
  TCanvas*  canvasPtvsMultErrSeparate[3][2];
  TCanvas* canvasYieldRatio = new TCanvas ("canvasYieldRatio", "canvasYieldRatio", 1300, 800);
  TCanvas* canvasYieldJet = new TCanvas ("canvasYieldJet", "canvasYieldJet", 1300, 800);

  TCanvas * 	canvasPtSpectraRatio= new TCanvas("canvasPtSpectraRatio", "canvasPtSpectraRatio", 1300, 800);
  canvasPtSpectraRatio->Divide(3,2);
TCanvas * 	canvasPtSpectraRatioAllMult= new TCanvas("canvasPtSpectraRatioAllMult", "canvasPtSpectraRatioAllMult", 1300, 400);
  canvasPtSpectraRatioAllMult->Divide(3,1);
    
TCanvas * 	canvasPtSpectraK0s= new TCanvas("canvasPtSpectraK0s", "canvasPtSpectraK0s", 1300, 800);
  canvasPtSpectraK0s->Divide(3,2);

  for (Int_t type=0; type<2; type++){
    for (Int_t iregion=0; iregion<3; iregion++){
      canvasYieldSeparate[iregion][type] = new TCanvas("canvasYield"+RegionType[iregion]+tipo[type], "canvasYield"+RegionType[iregion]+tipo[type], 1300, 800);
      canvasYieldErrSeparate[iregion][type] = new TCanvas("canvasRelErrorYield"+RegionType[iregion]+tipo[type], "canvasRelErrorYield"+RegionType[iregion]+tipo[type], 1300, 800);
      canvasPtvsMultErrSeparate[iregion][type] = new TCanvas("canvasRelErrorPtvsMult"+RegionType[iregion]+tipo[type], "canvasRelErrorPtvsMult"+RegionType[iregion]+tipo[type], 1300, 800);
      canvasPtSpectra[iregion][type] = new TCanvas("canvasPtSpectra"+RegionType[iregion]+tipo[type], "canvasPtSpectra"+RegionType[iregion]+tipo[type], 900, 1000);
      canvasPtSpectraMultRatio[iregion][type] = new TCanvas("canvasPtSpectraMultRatio"+RegionType[iregion]+tipo[type], "canvasPtSpectraMultRatio"+RegionType[iregion]+tipo[type], 1300, 800);
    }
  }
  Float_t LimSupYield[2] ={0.25, 0.02};
  Float_t  LimSupMultRatio[2]= {3.1,3.9};
  //  Float_t LimSupYieldF[2][3] ={{0.07, 0.25, 0.25},{0.003, 0.02, 0.02}};
  Float_t LimSupYieldF[2][3] ={{0.0699, 0.24999, 0.24999},{0.002999, 0.01999, 0.01999}};
  Float_t LimSupYieldFAll[2][3] ={{0.2699, 0.26999, 0.26999},{0.002999, 0.01999, 0.01999}};
  //  Float_t LimSupYieldFAll[2][3] ={{0.1, 0.1, 0.1},{0.006, 0.006, 0.006}}; This for Yield/event
  
  //Float_t LimSupSpectra[2][3] ={{0.02, 0.2, 0.2},{0.001, 0.01, 0.01}};
  //Float_t LimSupSpectra[2][3] ={{100, 100, 100},{10, 10, 10}};
  Float_t LimSupSpectra[2][3] ={{299.999, 299.999, 299.999},{29.999, 29.999, 29.999}};
  Float_t LimInfSpectra[2][3] ={{5*10e-7, 5*10e-7, 5*10e-7},{1.01*10e-7, 1.01*10e-7, 1.01*10e-7}};
  Float_t LimSupYieldErr[2] ={0.2, 0.3};
  Float_t LimInfYield[2]={10e-5, 10e-5};
  Float_t LimInfPt[2]={0.8, 1.2};
  Float_t LimSupPt[2]={2.8, 5};
  //  Float_t LimInfYield[2]={10e-6, 10e-7}; //this for yield/event
  //  Float_t LimInfRatioYield = 10e-5;
  Float_t LimInfRatioYield = 0.015;
  Float_t LimSupRatioYield = 0.13999; //0.12
  Float_t ScaleFactor[nummolt+1]={128,64, 32,16,8,1};
  TString sScaleFactor[nummolt+1]={" (x2^{7})"," (x2^{6})"," (x2^{5})"," (x2^{4})"," (x2^{3})",""};
  const Int_t numregions=3;

  TString PathIn[2];
  TString PathIn0[2];
  TString PathIn1[2];
  TString PathInYield[2];
  //  TString PathInDPhiProj[2];
  TFile * fileInYield[2];
  //  TFile * fileInDPhiProj[2];
  for (Int_t type=0; type<2; type++){
    if (type==0) PathIn0[type] =Form("_PtBinning%i",1);
    if(type>=0){
      if (!ishhCorr)      PathIn0[type] +="_"+tipo[type];
      PathIn0[type] +=Srap[israp];
      PathIn0[type] +=SSkipAssoc[SkipAssoc];
    }
    PathIn1[type]= PathIn0[type]+ Form("_SysPhi0_PtMin%.1f_", PtTrigMin);
    PathIn0[type]+= Form("_PtMin%.1f_", PtTrigMin);
  }

  TH1F* fHistYieldStat[2][numregions];
  TH1F* fHistYieldSist[2][numregions];
  TH1F* fHistPtvsMultStat[2][numregions];
  TH1F* fHistPtvsMultSist[2][numregions];
  TH1F* fHistPtvsMultStatFS[2][numregions];
  TH1F* fHistPtvsMultSistFS[2][numregions];
  TH1F* fHistPtvsMultStatMethodComp[2][numregions];
  TH1F* fHistPtvsMultSistMethodComp[2][numregions];
  TH1F* fHistYieldStatBlack[2][numregions];
  TH1F* fHistYieldSistBlack[2][numregions];
  TH1F* fHistYieldStatErr[2][numregions];
  TH1F* fHistYieldSistErr[2][numregions];
  TH1F* fHistPtvsMultStatErr[2][numregions];
  TH1F* fHistPtvsMultSistErr[2][numregions];
  TH1F* fHistYieldSistNoExtr[2][numregions];
  TH1F* fHistYieldStatJet[2][numregions];
  TH1F* fHistYieldSistJet[2][numregions];
  TH1F* fHistYieldSistNoExtrJet[2][numregions];
  TH1F* fHistYieldSistLowRelErr[2][numregions];
  TH1F* fHistYieldSistUpRelErr[2][numregions];
  TH1F* fHistYieldStatRatio[numregions];
  TH1F* fHistYieldSistRatio[numregions];
  TH1F* fHistYieldTotErrRatio[numregions];
  TH1F* fHistYieldSistNoExtrRatio[numregions];
  TH1F* fHistYieldStatRelErr[2][numregions];
  TH1F* fHistYieldSistRelErr[2][numregions];

  TF1* pol0Ratio[numregions];
  TF1* pol1Ratio[numregions];

  TH1F*	fHistYieldDatiPubblicati;
  TH1F*	fHistYieldDatiPubblicatiSistUncorr;
  TH1F*	fHistYieldDatiPubblicatiStat;
  TH1F*	fHistYieldDatiPubblicatiSistUncorrErr;
  TH1F*	fHistYieldDatiPubblicatiStatErr;
  TH1F*	fHistYieldDatiPubblicatiDenom;
  TH1F*	fHistYieldDatiPubblicatiDenomSistUncorr;
  TH1F*	fHistYieldDatiPubblicatiDenomStat;
  TH1F*	fHistYieldDatiPubblicatiDenomSistUncorrErr;
  TH1F*	fHistYieldDatiPubblicatiDenomStatErr;
  TH1F*	fHistYieldDatiPubblicatiRatio;
  TH1F*	fHistYieldDatiPubblicatiRatioSistUncorr;
  TH1F*	fHistYieldDatiPubblicatiRatioStat;

  TH1F*   fHistYieldErroriDatiPubblicatiStat; 
  TH1F*   fHistYieldErroriDatiPubblicatiSistUncorr; 
  TH1F*   fHistYieldErroriDatiPubblicati; 
  TH1F*   fHistYieldErroriDatiPubblicatiDenomStat; 
  TH1F*   fHistYieldErroriDatiPubblicatiDenomSistUncorr; 
  TH1F*   fHistYieldErroriDatiPubblicatiDenom; 

  TFile *filedatipubbl;
  TDirectoryFile *dir;
  filedatipubbl= new TFile("HEPData-1574358449-v1-Table_8c.root", "");
  if (!filedatipubbl) return;
  cout << " I got the files where published data are stored " << endl;
  dir  = (TDirectoryFile*)filedatipubbl->Get("Table 8c");
  if (!dir) return;

  TFile *filedatipubblNoExtr;
  TDirectoryFile *dirNE;
  filedatipubblNoExtr= new TFile("HEPData-ins1748157-v1-Table_9c.root", "");
  if (!filedatipubblNoExtr) return;
  cout << " I got the files where published data are stored " << endl;
  dirNE  = (TDirectoryFile*)filedatipubblNoExtr->Get("Table 9c");
  if (!dirNE) return;
  TH1F * fHistYieldPubNoExtrSyst[2]; 
  TH1F * fHistYieldPubNoExtr[2]; 
  TH1F * fHistYieldPubNoExtrRelSyst[2]; 
  TString nameYield[2] = {"Hist1D_y1", "Hist1D_y3"};
  TString nameYieldErr[2] = {"Hist1D_y1_e2", "Hist1D_y3_e2"};
  for (Int_t i =0; i<2 ; i++){
  fHistYieldPubNoExtr[i] = (TH1F*) dirNE->Get(nameYield[i]);
  fHistYieldPubNoExtrSyst[i] = (TH1F*) dirNE->Get(nameYieldErr[i]);
  fHistYieldPubNoExtrRelSyst[i] = (TH1F*)   fHistYieldPubNoExtrSyst[i]->Clone("fHistYieldPubNoExtrRelSystK0s");
  fHistYieldPubNoExtrRelSyst[i]->Divide(  fHistYieldPubNoExtr[i]);
  }

  Int_t bcorr=0;
  TH1F *fHistSpectrumStat[nummolt+1][2][numregions];
  TH1F *fHistSpectrumSist[nummolt+1][2][numregions];
  TH1F *fHistSpectrumStatScaled[nummolt+1][2][numregions];
  TH1F *fHistSpectrumSistScaled[nummolt+1][2][numregions];
  TH1F *fHistSpectrumStatScaledB[nummolt+1][2][numregions];
  TH1F *fHistSpectrumSistScaledB[nummolt+1][2][numregions];
  TH1F *fHistSpectrumStatRatio[nummolt+1][numregions];
  TH1F *fHistSpectrumSistRatio[nummolt+1][numregions];
  TH1F *fHistSpectrumStatDoubleRatio[nummolt+1][numregions];
  TH1F *fHistSpectrumSistDoubleRatio[nummolt+1][numregions];
  TH1F *fHistSpectrumStatMultRatio[nummolt+1][2][numregions];
  TH1F *fHistSpectrumSistMultRatio[nummolt+1][2][numregions];
  TSpline3*  splineK0s[nummolt+1][numregions];
  TF1 * lineat1 = new TF1 ("pol0", "pol0",0,8);
  lineat1->SetParameter(0,1);
  lineat1->SetLineColor(kBlack);
  Float_t AvgStat[nummolt+1][numregions]={0};
  Float_t AvgSist[nummolt+1][numregions]={0};


  gStyle->SetLegendBorderSize(1);
  TLegend *legendYield=new TLegend(0.6, 0.75, 0.74, 0.9);
  TLegend *legendRegion=new TLegend(0.75, 0.75, 0.9, 0.9);
  TLegend *legendJet=new TLegend(0.75, 0.75, 0.9, 0.9);

  TLegend *legendMult = new TLegend(0.6, 0.7, 0.9, 0.9);
  legendMult->SetHeader("Multiplicity classes");
  TLegend *legendAllMult = new TLegend(0.6, 0.6, 0.9, 0.9);
  legendAllMult->SetHeader("Multiplicity classes");
  TLegend *legendPY = new TLegend(0.5, 0.7, 0.9, 0.9);
  //ok!  TLegend *legendRegionAllF=new TLegend(0.06, 0.43, 0.5, 0.71);
  //  TLegend *legendRegionAllF=new TLegend(0.06, 0.43, 0.5, 0.71); //ok for margin 0.25 and no markers

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);
  TLegend *legendRegionAllF=new TLegend(0.15, 0.43, 0.59, 0.71);
  legendRegionAllF->SetFillStyle(0);
  legendRegionAllF->SetMargin(0.07);
  TLegendEntry * lReAll1[3];
  TLegendEntry *lReAll2[3];
  /*
  TLegend *legendRegionAllFJet=new TLegend(0.49, 0.82, 0.95, 0.905);
  legendRegionAllFJet->SetFillStyle(0);
  legendRegionAllFJet->SetMargin(1);
  TLegend *legendRegionAllFBulk=new TLegend(0.495, 0.72, 0.955, 0.805);
  legendRegionAllFBulk->SetFillStyle(0);
  legendRegionAllFBulk->SetMargin(1.12);
  TLegend *legendRegionAllFAll=new TLegend(0.4987, 0.62, 0.9587, 0.705);
  legendRegionAllFAll->SetFillStyle(0);
  legendRegionAllFAll->SetMargin(1.2);
  */
  TLegend *legendRegionAllFJet=new TLegend(0.15, 0.82, 0.61, 0.905);
  legendRegionAllFJet->SetFillStyle(0);
  legendRegionAllFJet->SetMargin(0.1);
  TLegend *legendRegionAllFBulk=new TLegend(0.15, 0.72, 0.61, 0.805);
  legendRegionAllFBulk->SetFillStyle(0);
  legendRegionAllFBulk->SetMargin(0.1);
  TLegend *legendRegionAllFAll=new TLegend(0.15, 0.62, 0.61, 0.705);
  legendRegionAllFAll->SetFillStyle(0);
  legendRegionAllFAll->SetMargin(0.1);

  TLegendEntry * lReAll1Bis[3];
  TLegendEntry *lReAll2Bis[3];

  TGraphAsymmErrors * gYield;
  Float_t   mult[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};

  //published yields ******************************************************

  TString yParticleNum = "y3"; //Xi
  fHistYieldDatiPubblicati=(TH1F*)dir->Get("Hist1D_"+yParticleNum);
  fHistYieldDatiPubblicatiStat=(TH1F*)  fHistYieldDatiPubblicati->Clone("Hist1D_Stat_"+yParticleNum);
  fHistYieldDatiPubblicatiSistUncorr=(TH1F*)  fHistYieldDatiPubblicati->Clone("Hist1D_SistUncorr_"+yParticleNum);
  fHistYieldErroriDatiPubblicatiStat=(TH1F*)dir->Get("Hist1D_"+yParticleNum+"_e1");  //this is the stat\istical uncertainty                                                                                   
  fHistYieldErroriDatiPubblicatiSistUncorr=(TH1F*)dir->Get("Hist1D_"+yParticleNum+"_e3"); //this is only the uncorr systematic error                                                                          
  fHistYieldErroriDatiPubblicati=(TH1F*)dir->Get("Hist1D_"+yParticleNum+"_e2"); //this is only the tota l systematic error (the statistic is 1/10 and the uncorrelated systematic is approx half)   

  TString yParticleDenom = "y1"; //K0s
  fHistYieldDatiPubblicatiDenom=(TH1F*)dir->Get("Hist1D_"+yParticleDenom);
  fHistYieldDatiPubblicatiDenomStat=(TH1F*)  fHistYieldDatiPubblicatiDenom->Clone("Hist1D_Stat_"+yParticleDenom);
  fHistYieldDatiPubblicatiDenomSistUncorr=(TH1F*)  fHistYieldDatiPubblicatiDenom->Clone("Hist1D_SistUncorr_"+yParticleDenom);
  fHistYieldErroriDatiPubblicatiDenomStat=(TH1F*)dir->Get("Hist1D_"+yParticleDenom+"_e1");  //this is the statistical uncertainty                                                                             
  fHistYieldErroriDatiPubblicatiDenomSistUncorr=(TH1F*)dir->Get("Hist1D_"+yParticleDenom+"_e3"); //this is only the uncorr systematic error                                                                  
  fHistYieldErroriDatiPubblicatiDenom=(TH1F*)dir->Get("Hist1D_"+yParticleDenom+"_e2"); //this is only the tota l systematic error (the statistic is 1/10 and the uncorrelated systematic is approx half)   

  for(Int_t k=1; k <= fHistYieldDatiPubblicati->GetNbinsX(); k++){
    fHistYieldDatiPubblicati->SetBinError(k,     fHistYieldErroriDatiPubblicati->GetBinContent(k));
    fHistYieldDatiPubblicatiDenom->SetBinError(k,     fHistYieldErroriDatiPubblicatiDenom->GetBinContent(k));
    fHistYieldDatiPubblicatiStat->SetBinError(k,     fHistYieldErroriDatiPubblicatiStat->GetBinContent(k));
    fHistYieldDatiPubblicatiDenomStat->SetBinError(k,     fHistYieldErroriDatiPubblicatiDenomStat->GetBinContent(k));
    fHistYieldDatiPubblicatiSistUncorr->SetBinError(k,     fHistYieldErroriDatiPubblicatiSistUncorr->GetBinContent(k));
    fHistYieldDatiPubblicatiDenomSistUncorr->SetBinError(k,     fHistYieldErroriDatiPubblicatiDenomSistUncorr->GetBinContent(k));

  }

  fHistYieldDatiPubblicati->Sumw2();
  fHistYieldDatiPubblicatiSistUncorr->Sumw2();
  fHistYieldDatiPubblicatiStat->Sumw2();
  fHistYieldDatiPubblicatiDenom->Sumw2();
  fHistYieldDatiPubblicatiDenomSistUncorr->Sumw2();
  fHistYieldDatiPubblicatiDenomStat->Sumw2();
  fHistYieldDatiPubblicatiRatio = (TH1F*)   fHistYieldDatiPubblicati->Clone("fHistYieldDatiPubblicatiRatio");
  fHistYieldDatiPubblicatiRatioStat = (TH1F*)   fHistYieldDatiPubblicatiStat->Clone("fHistYieldDatiPubblicatiRatioStat");
  fHistYieldDatiPubblicatiRatioSistUncorr = (TH1F*)   fHistYieldDatiPubblicatiSistUncorr->Clone("fHistYieldDatiPubblicatiRatioSistUncorr");
  fHistYieldDatiPubblicatiRatio->Divide(  fHistYieldDatiPubblicatiDenom);
  fHistYieldDatiPubblicatiRatioSistUncorr->Divide(  fHistYieldDatiPubblicatiDenomSistUncorr);
  fHistYieldDatiPubblicatiRatioStat->Divide(  fHistYieldDatiPubblicatiDenomStat);

  TH1F*   fHistYieldErroriRelDatiPubblicatiDenom = (TH1F*)fHistYieldErroriDatiPubblicatiDenom->Clone("fHistYieldErroriRelDatiPubblicatiDenom");
  TH1F*   fHistYieldErroriRelDatiPubblicati = (TH1F*)fHistYieldErroriDatiPubblicati->Clone("fHistYieldErroriRelDatiPubblicati");

  //************ end of published yields *********************************

  //************ Peng Yao yields *****************************************
  cout << "Peng Yao ratio: get the files... " << endl;
  TFile * filePY = new TFile("PengYaoXiK0sRatiopp.root", "");
  if (!filePY) {cout << "file PengYaoXiK0sRatiopp.root missing "<< endl; return;}
  TDirectoryFile *dirPY = (TDirectoryFile*)filePY->Get("Xi_toKRatio");
  if (!dirPY) {cout << " dir not there " << endl; return;}
  TH1F *fHistPengYaoRatio = (TH1F*) dirPY->FindObject("hJER");
  if (!fHistPengYaoRatio) {cout << " no histo " << endl; return;}
  TH1F *fHistPengYaoRatioOJ = (TH1F*) dirPY->FindObject("hUER");
  if (!fHistPengYaoRatioOJ) {cout << " no histo " << endl; return;}

  TCanvas * canvasPtSpectraRatioJetPY = new TCanvas("canvasPtSpectraRatioJetPY", "canvasPtSpectraRatioJetPY", 1300, 800);
  canvasPtSpectraRatioJetPY ->Divide(2,1);

  //*************end of PengYao yield ************************************

  cout << "\n**************************"<<endl;
  cout << "first part:I plot the yields vs mult for all the three regions " << endl;
  for (Int_t type=1; type>=0; type--){
    canvasYield[type]    = new TCanvas (Form("canvasYield_%i", type), Form("canvasYield_%i",type), 1300, 800);
    canvasYieldFinal[type]    = new TCanvas (Form("canvasYieldFinal_%i", type), Form("canvasYieldFinal_%i",type), 1300, 800);
    canvasPtvsMult[type]    = new TCanvas (Form("canvasPtMult_%i", type), Form("canvasPtMult_%i",type), 1300, 800);
    canvasPtvsMultMethodComp[type]    = new TCanvas (Form("canvasPtMultMethodComp_%i", type), Form("canvasPtMultMethodComp_%i",type), 1300, 800);

    for (Int_t iregion=0; iregion<numregions; iregion++){
      PathInYield[type] = Dir+"/DATA"+year0+"/PtSpectraBis" +hhCorr[ishhCorr] + PathIn0[type] +  RegionType[iregion];
      if (ZeroYieldLowPt && iregion==0 && type==1) PathInYield[type] = Dir+"/DATA"+year0+"/PtSpectraBis" +hhCorr[ishhCorr] + PathIn0[type] +  RegionType[iregion]+ "_ZeroYLowPtJet"; 
      //      PathInDPhiProj[type] = Dir+"/DATA"+year0+"/PtSpectra" +hhCorr[ishhCorr] + PathIn1[type] +  RegionType[iregion]+ ".root";
      if (isNormCorr)       PathInYield[type]+="_isNormCorr";
      PathInYield[type]+=".root";
      cout << PathInYield[type]<< endl;
      //      cout << PathInDPhiProj[type]<< endl;
      fileInYield[type] = new TFile(PathInYield[type], "");
      //      fileInDPhiProj[type] = new TFile(PathInDPhiProj[type], "");

      fHistYieldStat[type][iregion]=(TH1F*)      fileInYield[type]->Get("fHistYieldStat");
      fHistYieldSist[type][iregion]=(TH1F*)      fileInYield[type]->Get("fHistYieldSist");
      fHistYieldSistNoExtr[type][iregion]=(TH1F*)      fileInYield[type]->Get("fHistYieldSistNoExtr");
      if (type==1 && iregion==0){
      fHistYieldSistLowRelErr[type][iregion] = (TH1F*) fileInYield[type]->Get("fHistYieldSistLowRelErr");
      fHistYieldSistUpRelErr[type][iregion] = (TH1F*) fileInYield[type]->Get("fHistYieldSistUpRelErr");
      }
      fHistPtvsMultStat[type][iregion]=(TH1F*) fileInYield[type]->Get("fHistPtvsMultStat");
      fHistPtvsMultSist[type][iregion]=(TH1F*) fileInYield[type]->Get("fHistPtvsMultSist");
      fHistPtvsMultStatFS[type][iregion]=(TH1F*) fileInYield[type]->Get("fHistPtvsMultStatFromSpectrum");
      fHistPtvsMultSistFS[type][iregion]=(TH1F*) fileInYield[type]->Get("fHistPtvsMultSistFromSpectrum");
      fHistPtvsMultStatMethodComp[type][iregion]=(TH1F*)      fHistPtvsMultSistFS[type][iregion]->Clone("fHistPtvsMultStatMethodComp");
      fHistPtvsMultSistMethodComp[type][iregion]=(TH1F*)      fHistPtvsMultSistFS[type][iregion]->Clone("fHistPtvsMultSistMethodComp");


      if (!fHistYieldStat[type][iregion]) {cout << " missing histo 1 " << endl; return;}
      if (!fHistYieldSist[type][iregion]) {cout << " missing histo 1 " << endl; return;}
      if (!fHistYieldSistNoExtr[type][iregion]) {cout << " missing histo 1 " << endl; return;}

      fHistYieldStat[type][iregion]->SetName(Form("fHistYieldStat_type%i_%region%i", type, iregion));
      fHistYieldSist[type][iregion]->SetName(Form("fHistYieldSist_type%i_%region%i", type, iregion));
      fHistYieldSistNoExtr[type][iregion]->SetName(Form("fHistYieldSistNoExtr_type%i_%region%i", type, iregion));

      fHistYieldStatErr[type][iregion]=(TH1F*) fHistYieldStat[type][iregion]-> Clone("fHistYieldStatErr");
      fHistYieldSistErr[type][iregion]=(TH1F*) fHistYieldSist[type][iregion]-> Clone("fHistYieldSistErr");

      fHistPtvsMultStatErr[type][iregion]=(TH1F*) fHistPtvsMultStat[type][iregion]-> Clone("fHistPtvsMultStatErr");
      fHistPtvsMultSistErr[type][iregion]=(TH1F*) fHistPtvsMultSist[type][iregion]-> Clone("fHistPtvsMultSistErr");

      for (Int_t b=1; b<=fHistYieldStatErr[type][iregion]->GetNbinsX(); b++){
	if (fHistYieldStatErr[type][iregion]->GetBinContent(b)!=0){
	  fHistYieldStatErr[type][iregion]->SetBinContent(b, fHistYieldStatErr[type][iregion]->GetBinError(b)/fHistYieldStatErr[type][iregion]->GetBinContent(b));
	  fHistYieldSistErr[type][iregion]->SetBinContent(b, fHistYieldSistErr[type][iregion]->GetBinError(b)/fHistYieldSistErr[type][iregion]->GetBinContent(b));
	  fHistYieldStatErr[type][iregion]->SetBinError(b, 0);
	  fHistYieldSistErr[type][iregion]->SetBinError(b, 0);
	  fHistPtvsMultStatErr[type][iregion]->SetBinContent(b, fHistPtvsMultStatErr[type][iregion]->GetBinError(b)/fHistPtvsMultStatErr[type][iregion]->GetBinContent(b));
	  fHistPtvsMultSistErr[type][iregion]->SetBinContent(b, fHistPtvsMultSistErr[type][iregion]->GetBinError(b)/fHistPtvsMultSistErr[type][iregion]->GetBinContent(b));
	  fHistPtvsMultStatErr[type][iregion]->SetBinError(b, 0);
	  fHistPtvsMultSistErr[type][iregion]->SetBinError(b, 0);

	}
      }

      if (type==1  && iregion==0){
      fHistYieldSistLowRelErr[type][iregion]->SetName(Form("fHistYieldSistLowRelErr_type%i_%region%i", type, iregion));
      fHistYieldSistUpRelErr[type][iregion]->SetName(Form("fHistYieldSistUpRelErr_type%i_%region%i", type, iregion));
      }

      if (type==1 && iregion==0)   {
	gYield=(TGraphAsymmErrors*)fileInYield[type]->Get("Graph");
	if (!	gYield) return;
      }


      if (type==1) {
	fHistYieldStatRatio[iregion]=(TH1F*)      fHistYieldStat[type][iregion]->Clone(Form("fHistYieldStatRatio_iregion%i", iregion));
	fHistYieldSistRatio[iregion]=(TH1F*)      fHistYieldSist[type][iregion]->Clone(Form("fHistYieldSistRatio_iregion%i", iregion));
	fHistYieldSistNoExtrRatio[iregion]=(TH1F*)      fHistYieldSistNoExtr[type][iregion]->Clone(Form("fHistYieldSistNoExtrRatio_iregion%i", iregion));
      }

      else {
	fHistYieldStatRatio[iregion]->Divide(fHistYieldStat[type][iregion]);
	fHistYieldSistRatio[iregion]->Divide(fHistYieldSist[type][iregion]);
	fHistYieldSistNoExtrRatio[iregion]->Divide(fHistYieldSistNoExtr[type][iregion]);

	fHistYieldTotErrRatio[iregion]= (TH1F*)fHistYieldSistRatio[iregion]->Clone(Form("fHistYieldTotErrRatio_iregion%i", iregion));
	for (Int_t b=1; b <= fHistYieldTotErrRatio[iregion]->GetNbinsX(); b++){
	  fHistYieldTotErrRatio[iregion]->SetBinContent(b, fHistYieldStatRatio[iregion]->GetBinContent(b));
	  fHistYieldTotErrRatio[iregion]->SetBinError(b, sqrt(pow(fHistYieldStatRatio[iregion]->GetBinError(b),2) + pow(fHistYieldSistRatio[iregion]->GetBinError(b),2)));
	}
      }

      //graph for ratio
      Float_t Xl[nummolt+1]= {0};
      Float_t bin[nummolt+1]= {0};
      Float_t Xh[nummolt+1]= {0};
      Float_t multctrbin[nummolt+1] ={0} ;
      Float_t YieldRatio[nummolt+1]= {0};
      Float_t YieldRatioErrSistLow[nummolt+1]= {0};
      Float_t YieldRatioErrSistUp[nummolt+1]= {0};
      TGraphAsymmErrors * gYieldRatio;

      if (iregion==0 && type==0){
	for (Int_t m=0; m<nummolt+1; m++){
	  bin[m] =   fHistYieldStat[type][iregion]->FindBin(mult[m]);
	  multctrbin[m] =   fHistYieldStat[type][iregion]->GetXaxis()->GetBinCenter(bin[m]);
	  YieldRatio[m] =      fHistYieldStatRatio[iregion]->GetBinContent( bin[m]);
  
	  YieldRatioErrSistLow[m] = sqrt(pow(fHistYieldSistLowRelErr[1][iregion]->GetBinContent( bin[m]),2) + pow( fHistYieldSist[0][iregion]->GetBinError(bin[m])/ fHistYieldSist[0][iregion]->GetBinContent(bin[m]),2)) *  fHistYieldSistRatio[iregion]->GetBinContent(bin[m]);
	  YieldRatioErrSistUp[m] = sqrt(pow(fHistYieldSistUpRelErr[1][iregion]->GetBinContent( bin[m]),2) + pow( fHistYieldSist[0][iregion]->GetBinError(bin[m])/ fHistYieldSist[0][iregion]->GetBinContent(bin[m]),2)) *  fHistYieldSistRatio[iregion]->GetBinContent(bin[m]);

	}

	gYieldRatio= new TGraphAsymmErrors(nummolt,multctrbin,YieldRatio,Xl, Xh, YieldRatioErrSistLow, YieldRatioErrSistUp);

	gYieldRatio->SetMarkerColor(kBlue);
	gYieldRatio->SetLineColor(kBlue);
	gYieldRatio->SetMarkerStyle(1);
      }

      canvasYield[type]->cd();
//      fHistYieldStat[type][iregion]->GetYaxis()->SetRangeUser(LimInfYield[type], LimSupYield[type])  ;
//      fHistYieldSist[type][iregion]->GetYaxis()->SetRangeUser(LimInfYield[type], LimSupYield[type])  ;
//      fHistYieldSistNoExtr[type][iregion]->GetYaxis()->SetRangeUser(LimInfYield[type], LimSupYield[type])  ;

      StyleHisto(fHistYieldStat[type][iregion], LimInfYield[type], LimSupYield[type], Color[iregion], 1, titleYieldX, titleYieldYType[type], titleYield[type] + " yield vs multiplicity");
      StyleHisto(fHistYieldSist[type][iregion], LimInfYield[type], LimSupYield[type], Color[iregion], 1, titleYieldX, titleYieldYType[type], titleYield[type] + " yield vs multiplicity");
      StyleHisto(fHistYieldSistNoExtr[type][iregion], LimInfYield[type], LimSupYield[type], Color[iregion], 1, titleYieldX, titleYieldYType[type], titleYield[type] + " yield vs multiplicity");
  
      fHistYieldStat[type][iregion]->GetYaxis()->SetTitleOffset(1.2);
      fHistYieldSist[type][iregion]->GetYaxis()->SetTitleOffset(1.2);
      fHistYieldSistNoExtr[type][iregion]->GetYaxis()->SetTitleOffset(1.2);
      fHistYieldStat[type][iregion]->Draw("same e");  
      fHistYieldSist[type][iregion]->SetFillStyle(0);
      fHistYieldSist[type][iregion]->Draw("same e2");
      fHistYieldSistNoExtr[type][iregion]->SetFillStyle(3001);  
      fHistYieldSistNoExtr[type][iregion]->SetFillColorAlpha(Color[iregion], 1);
      fHistYieldSistNoExtr[type][iregion]->Draw("same e2");  

      if (type==1) {
	legendRegion->AddEntry(fHistYieldStat[type][iregion],SRegionType[iregion] , "ple");
      }
      if (iregion==0 && type==1) {
	legendYield->AddEntry(fHistYieldSist[type][iregion], "syst.", "ef");
	legendYield->AddEntry(fHistYieldSistNoExtr[type][iregion], "syst. (no extr)", "ef");
	legendYield->AddEntry(fHistYieldStat[type][iregion], "stat.", "pel");
      }


      //drawing histos separately
      canvasYieldSeparate[iregion][type]->cd();
      canvasYieldSeparate[iregion][type]->SetFillColor(0);
      canvasYieldSeparate[iregion][type]->SetTickx(1);
      canvasYieldSeparate[iregion][type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);
      //      gPad->SetTicks();
      //      TLatex name1 = TLatex(0.35, 0.8, "#scale[0.8]{#color[3]{ciao}}"); 
      TLegend *legendYieldF=new TLegend(0.18, 0.52, 0.35, 0.67);
      legendYieldF->AddEntry(fHistYieldSist[type][iregion], "syst. error", "ef");
      //	legendYieldF->AddEntry(fHistYieldSistNoExtr[type][iregion], "syst. (no extr)", "ef");
      legendYieldF->AddEntry(fHistYieldStat[type][iregion], "stat. error", "pe");

      fHistYieldStatBlack[type][iregion]= (TH1F*)       fHistYieldStat[type][iregion]->Clone("fHistYieldStatBlack"); 
      fHistYieldSistBlack[type][iregion]= (TH1F*)       fHistYieldSist[type][iregion]->Clone("fHistYieldSistBlack"); 
      fHistYieldStatBlack[type][iregion]->SetLineColor(1);
      fHistYieldStatBlack[type][iregion]->SetMarkerColor(1);
      fHistYieldStatBlack[type][iregion]->SetMarkerStyle(Ymarker);
      fHistYieldSistBlack[type][iregion]->SetLineColor(1);
      fHistYieldSistBlack[type][iregion]->SetMarkerColor(1);
      fHistYieldSistBlack[type][iregion]->SetMarkerStyle(Ymarker);

      TLegend *legendStatBox=new TLegend(0.7, 0.76, 0.9, 0.88);
      TLegendEntry*sb1=(TLegendEntry*)      legendStatBox->AddEntry(fHistYieldSistBlack[type][iregion], "syst. error", "ef");
      legendStatBox->AddEntry(fHistYieldStatBlack[type][iregion], "stat. error", "pe");

      //      TLegend *legendStatBoxBis=new TLegend(0.2, 0.53, 0.35, 0.63);
      TLegend *legendStatBoxBis=new TLegend(0.2, 0.48, 0.35, 0.58);
      legendStatBoxBis->AddEntry(fHistYieldSistBlack[type][iregion], "syst. error", "ef");
      legendStatBoxBis->AddEntry(fHistYieldStatBlack[type][iregion], "stat. error", "pe");

      TLegend *legendRegionF=new TLegend(0.52, 0.8, 0.9, 0.9);
      //legendRegionF->SetFillStyle(0);
      TLegendEntry * lRe1 = legendRegionF->AddEntry("", sRegion[iregion], "");
      lRe1->SetTextSize(0.06);
      lRe1->SetTextAlign(32);
      TLegendEntry *leR2=      legendRegionF->AddEntry("", sRegion1[iregion], "");
      leR2->SetTextSize(0.045);
      leR2->SetTextAlign(32);

      //      TLegend *Legend1=new TLegend(0.06,0.72,0.5,0.93); ok if margin is set to defaul
      TLegend *Legend1=new TLegend(0.16,0.72,0.5,0.93); //0.6 in place of 0.53
      //Legend1->SetTextAlign(13);
      Legend1->SetMargin(0);
      //      Legend1->SetFillStyle(0);
      Legend1->AddEntry("", "#bf{ALICE Preliminary}", "");
      Legend1->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
      Legend1->AddEntry(""/*(TObject*)0*/, NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");

      //TLegend *Legend2=new TLegend(0.06,0.75,0.5,0.93); // ok if margin is set to default
      TLegend *Legend2=new TLegend(0.16,0.75,0.5,0.93);
      //Legend1->SetTextAlign(13);
      //      Legend2->SetFillStyle(0);
      Legend2->SetMargin(0);
      Legend2->AddEntry("", "#bf{ALICE Preliminary}", "");
      Legend2->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
      Legend2->AddEntry(""/*(TObject*)0*/, NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");

      StyleHistoYield(fHistYieldStat[type][iregion], LimInfYield[type]+10e-7, LimSupYieldF[type][iregion], Color[iregion], YmarkerRegion[iregion], titleYieldX, titleYieldYType[type],"", MarkerSize[iregion],1.2,1.25 );
      StyleHistoYield(fHistYieldSist[type][iregion], LimInfYield[type]+10e-7, LimSupYieldF[type][iregion], Color[iregion], YmarkerRegion[iregion], titleYieldX, titleYieldYType[type], "",  MarkerSize[iregion], 1.2, 1.25);
      StyleHistoYield(fHistYieldSistNoExtr[type][iregion], LimInfYield[type]+10e-7, LimSupYieldF[type][iregion], Color[iregion], YmarkerRegion[iregion], titleYieldX, titleYieldYType[type], "",  MarkerSize[iregion], 1.2, 1.25);

      fHistYieldStat[type][iregion]->DrawClone("same e0x0");  
      fHistYieldSist[type][iregion]->SetFillStyle(0);
      fHistYieldSist[type][iregion]->DrawClone("same e2");
      fHistYieldSistNoExtr[type][iregion]->SetFillStyle(3001);  
      fHistYieldSistNoExtr[type][iregion]->SetFillColorAlpha(Color[iregion], 1);
      //      fHistYieldSistNoExtr[type][iregion]->Draw("same e2");  

      Legend1->Draw("");
      legendYieldF->Draw("");
      legendRegionF->Draw("");

      //Drawing Yields altogether (final)
      canvasYieldFinal[type]->cd();
      //      gPad->SetLogy();
      canvasYieldFinal[type]->SetFillColor(0);
      canvasYieldFinal[type]->SetTickx(1);
      canvasYieldFinal[type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);

      StyleHistoYield(fHistYieldStat[type][iregion], LimInfYield[type]+10e-7, LimSupYieldFAll[type][1], Color[iregion], YmarkerRegion[iregion], titleYieldX, titleYieldYType[type],"",  MarkerSize[iregion], 1.2 , 1.25);
      StyleHistoYield(fHistYieldSist[type][iregion], LimInfYield[type]+10e-7, LimSupYieldFAll[type][1], Color[iregion], YmarkerRegion[iregion], titleYieldX, titleYieldYType[type], "",  MarkerSize[iregion], 1.2, 1.25);
      StyleHistoYield(fHistYieldSistNoExtr[type][iregion], LimInfYield[type]+10e-7, LimSupYieldFAll[type][1], Color[iregion], YmarkerRegion[iregion], titleYieldX, titleYieldYType[type], "",  MarkerSize[iregion], 1.2,1.25);

      if(type==0){
	//      lReAll1[iregion] = legendRegionAllF->AddEntry("", sRegion[iregion], "");
	lReAll1[iregion] = legendRegionAllF->AddEntry(fHistYieldStat[type][iregion], sRegion[iregion], "p");
	//lReAll1[iregion] = legendRegionAllF->AddEntry(fHistYieldStat[type][iregion],"#color[0]{|}"+sRegion[iregion]+ "#color[0]{|}", "p");
	lReAll1[iregion]->SetTextSize(0.048);
	//      lReAll1[iregion]->SetTextAlign(32);
	lReAll2[iregion]=      legendRegionAllF->AddEntry("", sRegion1[iregion], "");
	lReAll2[iregion]->SetTextSize(0.038);
	//      lReAll2[iregion]->SetTextAlign(32);
      }
      
      /*
      Int_t m=nummolt;
      for (Int_t b=0; b<=fHistYieldStat[type][iregion]->GetNbinsX(); b++){
	if (	fHistYieldStat[type][iregion]->GetBinContent(b)!=0) m--;
	else continue;
	cout << " trigger fraction" << FracTrig[m] << endl;
	fHistYieldStat[type][iregion]->SetBinContent(b,       fHistYieldStat[type][iregion]->GetBinContent(b)*FracTrig[m]);
	fHistYieldStat[type][iregion]->SetBinError(b, fHistYieldStat[type][iregion]->GetBinError(b)*FracTrig[m]);
	fHistYieldSist[type][iregion]->SetBinContent(b,       fHistYieldSist[type][iregion]->GetBinContent(b)*FracTrig[m]);
	fHistYieldSist[type][iregion]->SetBinError(b, fHistYieldSist[type][iregion]->GetBinError(b)*FracTrig[m]);

      }
      */
      fHistYieldStat[type][iregion]->DrawClone("same e0x0");  
      fHistYieldSist[type][iregion]->SetFillStyle(0);
      fHistYieldSist[type][iregion]->DrawClone("same e2");
      fHistYieldSistNoExtr[type][iregion]->SetFillStyle(3001);  
      fHistYieldSistNoExtr[type][iregion]->SetFillColorAlpha(Color[iregion], 1);
      //      fHistYieldSistNoExtr[type][iregion]->Draw("same e2");  
      Legend2->Draw("");
      legendStatBox->Draw("");
      legendRegionAllF->Draw("");


      canvasYieldErrSeparate[iregion][type]->cd();
      for (Int_t b=1; b<=fHistYieldErroriDatiPubblicatiDenom->GetNbinsX(); b++){
	if (fHistYieldErroriDatiPubblicatiDenom->GetBinContent(b)!=0){
	  fHistYieldErroriRelDatiPubblicatiDenom->SetBinContent(b, fHistYieldErroriDatiPubblicatiDenom->GetBinContent(b)/fHistYieldDatiPubblicatiDenom->GetBinContent(b));
	  fHistYieldErroriRelDatiPubblicatiDenom->SetBinError(b,0);
	}
      }
      for (Int_t b=1; b<=fHistYieldErroriDatiPubblicati->GetNbinsX(); b++){
	if (fHistYieldErroriDatiPubblicati->GetBinContent(b)!=0){
	  fHistYieldErroriRelDatiPubblicati->SetBinContent(b, fHistYieldErroriDatiPubblicati->GetBinContent(b)/fHistYieldDatiPubblicati->GetBinContent(b));
	  fHistYieldErroriRelDatiPubblicati->SetBinError(b,0);
	}
      }

      StyleHisto(fHistYieldStatErr[type][iregion], LimInfYield[type], LimSupYieldErr[type], Color[iregion], 33, titleYieldX, titleYieldYRelErr, "");
      StyleHisto(fHistYieldSistErr[type][iregion], LimInfYield[type], LimSupYieldErr[type], Color[iregion], 24, titleYieldX, titleYieldYRelErr, "");
      StyleHisto(  fHistYieldErroriRelDatiPubblicatiDenom, LimInfYield[type], LimSupYieldErr[type], 922, 33, titleYieldX, titleYieldYRelErr, "");
      StyleHisto(  fHistYieldErroriRelDatiPubblicati, LimInfYield[type], LimSupYieldErr[type], 922, 33, titleYieldX, titleYieldYRelErr, "");
      for (Int_t i=0;i<2; i++){
      StyleHisto(fHistYieldPubNoExtrRelSyst[i], LimInfYield[type], LimSupYieldErr[type], 1, 33, titleYieldX, titleYieldYRelErr, "");
      }


      TLegend* legendYieldErr = new TLegend(0.6, 0.6, 0.9, 0.9);
      legendYieldErr->AddEntry(fHistYieldStatErr[type][iregion], "stat.", "pl");
      legendYieldErr->AddEntry(fHistYieldSistErr[type][iregion], "syst.", "pl");
      legendYieldErr->AddEntry(fHistYieldPubNoExtrRelSyst[0], "syst. pub. no extr.", "pl");
      legendYieldErr->AddEntry(fHistYieldErroriRelDatiPubblicatiDenom, "syst. pub.", "pl");
      fHistYieldStatErr[type][iregion]->Draw("same p"); 
      fHistYieldSistErr[type][iregion]->Draw("same p"); 
      if (type==0) {
	fHistYieldErroriRelDatiPubblicatiDenom->Draw("same p hist");
	fHistYieldPubNoExtrRelSyst[0]->Draw("same p hist");
      }
      else {
	fHistYieldErroriRelDatiPubblicati->Draw("same p hist");
	fHistYieldPubNoExtrRelSyst[0]->Draw("same p hist");
      }
      legendYieldErr->Draw("");
      
      canvasPtvsMultErrSeparate[iregion][type]->cd();
      StyleHisto(fHistPtvsMultStatErr[type][iregion], LimInfYield[type], LimSupYieldErr[type], Color[iregion], 33, titleYieldX, titleYieldYRelErr, "");
      StyleHisto(fHistPtvsMultSistErr[type][iregion], LimInfYield[type], LimSupYieldErr[type], Color[iregion], 24, titleYieldX, titleYieldYRelErr, "");

      TLegend* legendPtvsMultErr = new TLegend(0.6, 0.6, 0.9, 0.9);
      legendPtvsMultErr->AddEntry(fHistPtvsMultStatErr[type][iregion], "stat.", "pl");
      legendPtvsMultErr->AddEntry(fHistPtvsMultSistErr[type][iregion], "syst.", "pl");
      //      legendYieldErr->AddEntry(fHistYieldErroriRelDatiPubblicatiDenom, "syst. pub.", "pl");
      fHistPtvsMultStatErr[type][iregion]->Draw("same p"); 
      fHistPtvsMultSistErr[type][iregion]->Draw("same p"); 
      legendPtvsMultErr->Draw("");


      if (type==0){
	canvasYieldRatio->cd();
	canvasYieldRatio->SetFillColor(0);
	canvasYieldRatio->SetTickx(1);
	canvasYieldRatio->SetTicky(1);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendFont(42);

	//	TLegend *Legend1B=new TLegend(0.1,0.72,0.38,0.92);
	TLegend *Legend1B=new TLegend(0.62,0.72,0.9,0.92);
	Legend1B->SetFillStyle(0);
	TLegendEntry* E1Bis =      Legend1B->AddEntry("", "#bf{ALICE Preliminary}", "");
	E1Bis->SetTextAlign(32);
	TLegendEntry* E2Bis =        Legend1B->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
	E2Bis->SetTextAlign(32);
	TLegendEntry* E3Bis =   	Legend1B->AddEntry(""/*(TObject*)0*/, "#it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
	E3Bis->SetTextAlign(32);

	//	Legend1B->AddEntry(""/*(TObject*)0*/, NameP[1]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV", "");

	//	fHistYieldStatRatio[iregion]->SetTitle(TitleYieldRatio);
	//	fHistYieldSistRatio[iregion]->SetTitle(TitleYieldRatio);
	//	fHistYieldSistNoExtrRatio[iregion]->SetTitle(TitleYieldRatio);
	//	gYieldRatio->Draw("same AP");

	//	fHistYieldStatRatio[iregion]->GetXaxis()->SetRangeUser(0,27);
	//	fHistYieldSistRatio[iregion]->GetXaxis()->SetRangeUser(0,27);
	//	fHistYieldSistNoExtrRatio[iregion]->GetXaxis()->SetRangeUser(0,27);

	StyleHistoYield(	fHistYieldDatiPubblicatiRatio, LimInfRatioYield, LimSupRatioYield, 857, YmarkerRegion[iregion], titleYieldX, TitleYYieldRatio,  "",  MarkerSize[iregion],1.4, 1.25);
	StyleHistoYield(	fHistYieldDatiPubblicatiRatioStat, LimInfRatioYield, LimSupRatioYield, 857, YmarkerRegion[iregion], titleYieldX, TitleYYieldRatio,  "",  MarkerSize[iregion],1.4, 1.25);
	StyleHistoYield(	fHistYieldDatiPubblicatiRatioSistUncorr, LimInfRatioYield, LimSupRatioYield, 598, YmarkerRegion[iregion], titleYieldX, TitleYYieldRatio,  "",  MarkerSize[iregion], 1.4,1.25);
	StyleHistoYield(	fHistYieldStatRatio[iregion]	, LimInfRatioYield, LimSupRatioYield, Color[iregion], YmarkerRegion[iregion], titleYieldX, TitleYYieldRatio,  "",  MarkerSize[iregion],1.4, 0.9);
	StyleHistoYield( 	fHistYieldSistRatio[iregion], LimInfRatioYield, LimSupRatioYield, Color[iregion], YmarkerRegion[iregion], titleYieldX, TitleYYieldRatio,  "",  MarkerSize[iregion],1.4, 0.9);
	StyleHistoYield( 	fHistYieldSistNoExtrRatio[iregion], LimInfRatioYield, LimSupRatioYield, Color[iregion], YmarkerRegion[iregion], titleYieldX, TitleYYieldRatio,  "",  MarkerSize[iregion], 1.4,0.9);
	StyleHistoYield(	fHistYieldTotErrRatio[iregion]	, LimInfRatioYield, LimSupRatioYield, Color[iregion], YmarkerRegion[iregion], titleYieldX, TitleYYieldRatio,  "",  MarkerSize[iregion],1.4, 0.9);
	if (iregion==0)	lReAll1Bis[iregion] = legendRegionAllFJet->AddEntry( fHistYieldStatRatio[iregion], sRegion[iregion], "p");
	else if (iregion==1)lReAll1Bis[iregion] = legendRegionAllFBulk->AddEntry( fHistYieldStatRatio[iregion], sRegion[iregion], "p");
	else lReAll1Bis[iregion] = legendRegionAllFAll->AddEntry( fHistYieldStatRatio[iregion], sRegion[iregion], "p");

	lReAll1Bis[iregion]->SetTextSize(0.045);
	//	lReAll1Bis[iregion]->SetTextAlign(32);
	if (iregion==0)	lReAll2Bis[iregion]=      legendRegionAllFJet->AddEntry("", sRegion1[iregion], "");
	else	if (iregion==1)	lReAll2Bis[iregion]=      legendRegionAllFBulk->AddEntry("", sRegion1[iregion], "");
	else	lReAll2Bis[iregion]=      legendRegionAllFAll->AddEntry("", sRegion1[iregion], "");
	lReAll2Bis[iregion]->SetTextSize(0.035);
	//	lReAll2Bis[iregion]->SetTextAlign(32);

	fHistYieldStatRatio[iregion]->GetYaxis()->SetTitleSize(0.07);
	fHistYieldSistRatio[iregion]->GetYaxis()->SetTitleSize(0.07);

	pol0Ratio[iregion] = new TF1 (Form("pol0Ratio_%i", iregion), "pol0", 0, 30);
	pol1Ratio[iregion] = new TF1 (Form("pol1Ratio_%i", iregion), "pol1", 0, 30);
	pol0Ratio[iregion]->SetLineColor(Color[iregion]);
	pol0Ratio[iregion]->SetLineStyle(1);
	pol1Ratio[iregion]->SetLineColor(Color[iregion]);
	pol1Ratio[iregion]->SetLineStyle(2);


	fHistYieldStatRatio[iregion]->DrawClone("same e0x0");  
	fHistYieldSistRatio[iregion]->SetFillStyle(0);
	fHistYieldSistRatio[iregion]->DrawClone("same e2");

	fHistYieldTotErrRatio[iregion]->Fit(pol0Ratio[iregion], "R0");
	fHistYieldTotErrRatio[iregion]->Fit(pol1Ratio[iregion], "R0");
	pol0Ratio[iregion]->Draw("same");
	pol1Ratio[iregion]->Draw("same");
	fHistYieldTotErrRatio[iregion]->SetLineColorAlpha(Color[iregion], 0);

	if (iregion==2){
	fHistYieldStatRatio[iregion]->DrawClone("same e0x0");  
	fHistYieldSistRatio[iregion]->SetFillStyle(0);
	fHistYieldSistRatio[iregion]->DrawClone("same e2");
	//	cout << "\n\n\n hey thee " << 	fHistYieldSistRatio[1]->GetBinContent(3)<< endl;
	fHistYieldStatRatio[1]->DrawClone("same e0x0");  
	fHistYieldSistRatio[1]->DrawClone("same e2");

	}

	fHistYieldSistNoExtrRatio[iregion]->SetFillStyle(3001);  
	fHistYieldSistNoExtrRatio[iregion]->SetFillColorAlpha(Color[iregion], 1);
	//	fHistYieldSistNoExtrRatio[iregion]->DrawClone("same e2");  

	fHistYieldDatiPubblicatiRatio->SetFillStyle(0);
	fHistYieldDatiPubblicatiRatioSistUncorr->SetFillStyle(0);
	fHistYieldDatiPubblicatiRatio->GetXaxis()->SetRangeUser(0,27);
	fHistYieldDatiPubblicatiRatioStat->GetXaxis()->SetRangeUser(0,27);
	fHistYieldDatiPubblicatiRatioSistUncorr->GetXaxis()->SetRangeUser(0,27);
	//			fHistYieldDatiPubblicatiRatio->Draw("same e2");
	//			fHistYieldDatiPubblicatiRatioSistUncorr->Draw("same e2");
	//			fHistYieldDatiPubblicatiRatioStat->Draw("same");
	//fHistYieldStatRatio[iregion]->GetXaxis()->SetRangeUser(0,25);
	//fHistYieldSistRatio[iregion]->GetXaxis()->SetRangeUser(0,25);
	//fHistYieldSistNoExtrRatio[iregion]->GetXaxis()->SetRangeUser(0,25);

	//	gYieldRatio->Draw("same P");
	Legend1B->Draw("");
	legendRegionAllFJet->Draw("");
	legendRegionAllFBulk->Draw("");
	legendRegionAllFAll->Draw("");
	legendStatBoxBis->Draw("");

      }

      Int_t nscale =25;
      TGraphAsymmErrors*	gYieldJet;
      Float_t YieldJetErrSistUp[nummolt+1]={0};
      Float_t YieldJetErrSistLow[nummolt+1]={0};
      Float_t YieldScaled[nummolt+1]={0};
      if (iregion==0){
	canvasYieldJet->cd();
	fHistYieldStatJet[type][iregion]=(TH1F*)fHistYieldStat[type][iregion]->Clone(Form("fHistYieldStatJet_type%i", type));
	fHistYieldSistJet[type][iregion]=(TH1F*)fHistYieldSist[type][iregion]->Clone(Form("fHistYieldSistJet_type%i", type));
	fHistYieldSistNoExtrJet[type][iregion]=(TH1F*)fHistYieldSistNoExtr[type][iregion]->Clone(Form("fHistYieldSistNoExtrJet_type%i", type));
	if(type==1){
	  fHistYieldStatJet[type][iregion]->Scale(nscale);
	  fHistYieldSistJet[type][iregion]->Scale(nscale);
	  fHistYieldSistNoExtrJet[type][iregion]->Scale(nscale);
	}
	//	fHistYieldStat[type][iregion]->GetYaxis()->SetRangeUser(10e-6, 0.1);
	StyleHisto(fHistYieldStatJet[type][iregion], LimInfYield[type], 0.07, ColorType[type], 1, titleYieldX, titleYieldY, "In-jet yield vs multiplicity");
	StyleHisto(fHistYieldSistJet[type][iregion], LimInfYield[type], 0.07, ColorType[type], 1, titleYieldX, titleYieldY,"In-jet yield vs multiplicity");
	StyleHisto(fHistYieldSistNoExtrJet[type][iregion], LimInfYield[type], 0.07,  ColorType[type], 1, titleYieldX, titleYieldY,  "In-jet yield vs multiplicity");
	if(type==0)	legendJet->AddEntry(fHistYieldStatJet[type][iregion],Stipo[type] , "ple");
	else 	legendJet->AddEntry(fHistYieldStatJet[type][iregion],Stipo[type]+ Form(" (x%i)",nscale) , "ple");

	fHistYieldStatJet[type][iregion]->GetYaxis()->SetTitleOffset(1.2);
	fHistYieldSistJet[type][iregion]->GetYaxis()->SetTitleOffset(1.2);
	fHistYieldSistNoExtrJet[type][iregion]->GetYaxis()->SetTitleOffset(1.2);

	if (type==1){
	  for (int i=0;i<gYield->GetN();i++) {
	    YieldScaled[i]=	  gYield->GetY()[i] * nscale;
	    YieldJetErrSistUp[i] = gYield->GetErrorYhigh(i) * nscale;
	    YieldJetErrSistLow[i]  =  gYield->GetErrorYlow(i)*nscale;
	  }
	  gYieldJet= new TGraphAsymmErrors(nummolt,multctrbin,YieldScaled,Xl, Xh, YieldJetErrSistLow, YieldJetErrSistUp);	
	}

	fHistYieldStatJet[type][iregion]->DrawClone("same e");  
	fHistYieldSistJet[type][iregion]->SetFillStyle(0);
	fHistYieldSistJet[type][iregion]->DrawClone("same e2");
	fHistYieldSistNoExtrJet[type][iregion]->SetFillStyle(3001);  
	fHistYieldSistNoExtrJet[type][iregion]->SetFillColorAlpha(ColorType[type], 1);
	fHistYieldSistNoExtrJet[type][iregion]->DrawClone("same e2");
	if (type==0){
	  //	gYieldJet->Draw("same P");
	}
	if (type==0)	legendYield->Draw("");
	if (type==0)	legendJet->Draw("");

      }

      //      for (Int_t m=0; m<nummolt+1; m++){
      for (Int_t m=nummolt; m>=0; m--){
	cout << " ratio between spectra; m= "<< m  << endl;
	fHistSpectrumSist[m][type][iregion]=(TH1F*)fileInYield[type]->Get("fHistSpectrumSistAll_"+Smolt[m]);
	fHistSpectrumStat[m][type][iregion]=(TH1F*)fileInYield[type]->Get("fHistSpectrum_"+Smolt[m]);
	
	if (!fHistSpectrumStat[m][type][iregion]) { cout << " no hist spectrum stat" << endl;  return;}
	if (!fHistSpectrumSist[m][type][iregion]) { cout << " no hist spectrum sist" << endl;  return;}

	cout << "\n\n m " << m << " type " << type << " iregion " << iregion << endl;
	for (Int_t b=1; b<= fHistSpectrumStat[m][type][iregion]->GetNbinsX(); b++){
	  cout << b << " " << fHistSpectrumStat[m][type][iregion]->GetBinContent(b) << endl;
	  if (type==1 && m==4 && b==fHistSpectrumStat[m][type][iregion]->GetNbinsX() && iregion!=2) {
	    fHistSpectrumStat[m][type][iregion]->SetBinContent(b, 0);
	    fHistSpectrumStat[m][type][iregion]->SetBinError(b, 0);
	    fHistSpectrumSist[m][type][iregion]->SetBinContent(b, 0);
	    fHistSpectrumSist[m][type][iregion]->SetBinError(b, 0);
	  }
	}

	fHistSpectrumSist[m][type][iregion]->SetName("fHistSpectrumSist_"+Smolt[m]+Form("_type%i_region%i", type, iregion));
	fHistSpectrumStat[m][type][iregion]->SetName("fHistSpectrumStat_"+Smolt[m]+Form("_type%i_region%i", type, iregion));
	fHistSpectrumSistScaled[m][type][iregion]=(TH1F*)	fHistSpectrumSist[m][type][iregion]->Clone("fHistSpectrumSistScaled_"+Smolt[m]+Form("_type%i_region%i", type,iregion));
	fHistSpectrumStatScaled[m][type][iregion]=(TH1F*)	fHistSpectrumStat[m][type][iregion]->Clone("fHistSpectrumStatScaled_"+Smolt[m]+Form("_type%i_region%i", type,iregion));

	if (type==1){
	  fHistSpectrumSistRatio[m][iregion]=(TH1F*)	fHistSpectrumSist[m][type][iregion]->Clone("fHistSpectrumSistRatio_"+Smolt[m]+Form("_region%i", iregion));
	  fHistSpectrumStatRatio[m][iregion]=(TH1F*)	fHistSpectrumStat[m][type][iregion]->Clone("fHistSpectrumStatRatio_"+Smolt[m]+Form("_region%i", iregion));

	}
	else {
	  splineK0s[m][iregion] = new TSpline3(fHistSpectrumStat[m][type][iregion] , Form("splineK0s_m%i_region%i",m, iregion),0.1, 3 );
     
	  //I calculate the average stat and sist error in the pt region where I do the spline

	  Int_t Counter=0;
	  for (Int_t b=1; b<=  fHistSpectrumStat[m][type][iregion]->GetNbinsX(); b++){
	    if (fHistSpectrumStat[m][type][iregion]->GetBinContent(b)==0) continue;
	    if ( fHistSpectrumStat[m][type][iregion]->GetXaxis()->GetBinLowEdge(b) <2.){
	      Counter++;
	      AvgStat[m][iregion] += fHistSpectrumStat[m][0][iregion]->GetBinError(b)/ fHistSpectrumStat[m][0][iregion]->GetBinContent(b);
	      AvgSist[m][iregion] += fHistSpectrumSist[m][0][iregion]->GetBinError(b)/ fHistSpectrumStat[m][0][iregion]->GetBinContent(b);
	    }
	  }
	  AvgStat[m][iregion]=  AvgStat[m][iregion]/Counter;
	  AvgSist[m][iregion] =  AvgSist[m][iregion]/Counter;
	  cout << " m " << m << "iregion " << iregion << " avg stat rel error" << AvgStat[m][iregion]<< endl;

	  for (Int_t b=2; b<=  fHistSpectrumStatRatio[m][iregion]->GetNbinsX(); b++){
	    if ((iregion!=2 && m==4 && b==fHistSpectrumStatRatio[m][iregion]->GetNbinsX()) || (b==2 && iregion==0)){
	      fHistSpectrumStatRatio[m][iregion]-> SetBinContent(b,0);
	      fHistSpectrumSistRatio[m][iregion]-> SetBinContent(b,0);
	      continue;
	    }
	  
	    if ( fHistSpectrumStatRatio[m][iregion]->GetXaxis()->GetBinLowEdge(b) <2.){
	      fHistSpectrumStatRatio[m][iregion]-> SetBinContent(b,fHistSpectrumStat[m][1][iregion]->GetBinContent(b)/splineK0s[m][iregion]->Eval(fHistSpectrumStat[m][1][iregion]->GetBinCenter(b)));
	      fHistSpectrumStatRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumStat[m][1][iregion]->GetBinError(b)/fHistSpectrumStat[m][1][iregion]->GetBinContent(b),2) + pow(AvgStat[m][iregion],2)) *  fHistSpectrumStatRatio[m][iregion]->GetBinContent(b) ); 
	      fHistSpectrumSistRatio[m][iregion]-> SetBinContent(b,fHistSpectrumSist[m][1][iregion]->GetBinContent(b)/splineK0s[m][iregion]->Eval(fHistSpectrumSist[m][1][iregion]->GetBinCenter(b)));
	      fHistSpectrumSistRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumSist[m][1][iregion]->GetBinError(b)/fHistSpectrumSist[m][1][iregion]->GetBinContent(b),2) + pow(AvgSist[m][iregion],2)) *   fHistSpectrumSistRatio[m][iregion]->GetBinContent(b) ); 

	    }
	    else {
	      bcorr = fHistSpectrumStat[m][0][iregion]->FindBin(fHistSpectrumStat[m][1][iregion]->GetBinCenter(b));
	      fHistSpectrumStatRatio[m][iregion]-> SetBinContent(b,fHistSpectrumStat[m][1][iregion]->GetBinContent(b)/fHistSpectrumStat[m][0][iregion]->GetBinContent(bcorr));
	      fHistSpectrumStatRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumStat[m][1][iregion]->GetBinError(b)/fHistSpectrumStat[m][1][iregion]->GetBinContent(b),2) + pow(fHistSpectrumStat[m][0][iregion]->GetBinError(bcorr)/fHistSpectrumStat[m][0][iregion]->GetBinContent(bcorr),2))* fHistSpectrumStatRatio[m][iregion]->GetBinContent(b));
	      fHistSpectrumSistRatio[m][iregion]-> SetBinContent(b,fHistSpectrumSist[m][1][iregion]->GetBinContent(b)/fHistSpectrumSist[m][0][iregion]->GetBinContent(bcorr));
	      fHistSpectrumSistRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumSist[m][1][iregion]->GetBinError(b)/fHistSpectrumSist[m][1][iregion]->GetBinContent(b),2) + pow(fHistSpectrumSist[m][0][iregion]->GetBinError(bcorr)/fHistSpectrumSist[m][0][iregion]->GetBinContent(bcorr),2))* fHistSpectrumSistRatio[m][iregion]->GetBinContent(b));

	    }
	  } // end loop pt bin
	  canvasPtSpectraK0s->cd(m+1);
	  fHistSpectrumStat[m][0][iregion]->Draw("same ");
	  fHistSpectrumSist[m][0][iregion]->SetFillStyle(0);
	  fHistSpectrumSist[m][0][iregion]->Draw("same e2");
	  splineK0s[m][iregion]->SetLineColor(Color[iregion]);
	  splineK0s[m][iregion]->Draw("same");

	  canvasPtSpectraRatio->cd(iregion+1);
	  gPad->SetLeftMargin(0.15);
	  if (iregion==0){
	    fHistSpectrumStatRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(0.6),0);
	    fHistSpectrumStatRatio[m][iregion]->SetBinError(fHistSpectrumStatRatio[m][iregion]->FindBin(0.6),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(0.6),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinError(fHistSpectrumSistRatio[m][iregion]->FindBin(0.6),0);
	  }
	  if ((iregion==0 || iregion==1) && m==4){
	    fHistSpectrumStatRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(5),0);
	    fHistSpectrumStatRatio[m][iregion]->SetBinError(fHistSpectrumStatRatio[m][iregion]->FindBin(5),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(5),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinError(fHistSpectrumSistRatio[m][iregion]->FindBin(5),0);
	  }
	  StyleHisto(fHistSpectrumStatRatio[m][iregion], 0.00001, 0.4,  ColorMult[m], 33, titleX, TitleYPtRatio,  "");
	  //StyleHisto(fHistSpectrumStatRatio[m][iregion], 0.00001, 0.4,  ColorMult[m], 33, titleX, "",  "");
	  //	  StyleHisto(fHistSpectrumSistRatio[m][iregion], 0.00001, 0.4,  ColorMult[m], 33, titleX,  TitleYPtRatio,  TitleYYieldRatio + " spectra ratio " + RegionType[iregion]);
	  StyleHisto(fHistSpectrumSistRatio[m][iregion], 0.00001, 0.4,  ColorMult[m], 33, titleX,  TitleYPtRatio, "");
	  //	  StyleHisto(fHistSpectrumSistRatio[m][iregion], 0.00001, 0.4,  ColorMult[m], 33, titleX,  "", "");
	  fHistSpectrumStatRatio[m][iregion]->SetMarkerSize(0.8);
	  fHistSpectrumSistRatio[m][iregion]->SetMarkerSize(0.8);
	  if (m==0 || m==4 || m==nummolt){
	    if (iregion==0)	    legendMult->AddEntry(	  fHistSpectrumStatRatio[m][iregion], SmoltLegend[m], "pl");
	    fHistSpectrumStatRatio[m][iregion]->DrawClone("same e");
	  fHistSpectrumSistRatio[m][iregion]->SetFillStyle(0);
	  fHistSpectrumSistRatio[m][iregion]->DrawClone("same e2");
	  if (m==nummolt && iregion==0) legendMult->Draw("");
	  }

	  canvasPtSpectraRatioAllMult->cd(iregion+1);
	  gPad->SetLeftMargin(0.15);
	    if (iregion==0)	    legendAllMult->AddEntry(	  fHistSpectrumStatRatio[m][iregion], SmoltLegend[m], "pl");
	    fHistSpectrumStatRatio[m][iregion]->DrawClone("same e");
	  fHistSpectrumSistRatio[m][iregion]->SetFillStyle(0);
	  fHistSpectrumSistRatio[m][iregion]->DrawClone("same e2");
	  if (m==nummolt && iregion==0) legendAllMult->Draw("");

	  cout << "Peng Yao ratio " << endl;
	  if (iregion==0){
	  canvasPtSpectraRatioJetPY->cd(1);
	  gPad->SetLeftMargin(0.15);
	   
	    if (m==nummolt){
	      fHistPengYaoRatio->SetTitle(SRegionType[iregion]);
	      fHistPengYaoRatio->Scale(2);
	      StyleHisto(fHistPengYaoRatio, 0.00001, 0.15, 857, 33, titleX,  TitleYPtRatio,  "Near-side jet");
	      fHistPengYaoRatio->GetYaxis()->SetRangeUser(10e-5,0.15);
	      fHistPengYaoRatio->Draw("");
	      legendPY->AddEntry(fHistPengYaoRatio, "Jet finder", "pl");

	    }
	    if (m==0 || m==4 || m==nummolt){
	      fHistSpectrumStatRatio[m][iregion]-> SetTitle(SRegionType[iregion]);
	      fHistSpectrumSistRatio[m][iregion]-> SetTitle(SRegionType[iregion]);
	      fHistSpectrumStatRatio[m][iregion]->	   GetYaxis()->SetRangeUser(10-5,0.15);
	      fHistSpectrumSistRatio[m][iregion]->	   GetYaxis()->SetRangeUser(10-5,0.15);
	   	   
	      legendPY->AddEntry(	  fHistSpectrumStatRatio[m][iregion], "This analysis "+ SmoltLegend[m], "pl");
	      fHistSpectrumStatRatio[m][iregion]->DrawClone("same e");
	      fHistSpectrumSistRatio[m][iregion]->Draw("same e2");
	    }
	    if (m==nummolt && iregion==0) legendPY->Draw("");
	  }
	  if (iregion==1){
	  canvasPtSpectraRatioJetPY->cd(2);
	  gPad->SetLeftMargin(0.15);
	   
	    if (m==nummolt){
	      fHistPengYaoRatioOJ->Scale(2);
	      fHistPengYaoRatioOJ->SetTitle(SRegionType[iregion]);
	      StyleHisto(fHistPengYaoRatioOJ, 0.00001, 0.4, 857, 33, titleX,  TitleYPtRatio,  "Out-of-jet");
	      fHistPengYaoRatioOJ->GetYaxis()->SetRangeUser(10e-5,0.5);
	      fHistPengYaoRatioOJ->Draw("");

	    }
	    if ( m==4 || m==0 || m==nummolt){
	      fHistSpectrumStatRatio[m][iregion]->	   GetYaxis()->SetRangeUser(10-5,0.5);
	      fHistSpectrumSistRatio[m][iregion]->	   GetYaxis()->SetRangeUser(10-5,0.5);

	      fHistSpectrumStatRatio[m][iregion]->DrawClone("same e");
	      fHistSpectrumSistRatio[m][iregion]->Draw("same e2");
	    }
	    if (m==nummolt && iregion==0) legendPY->Draw("");
	  }


	  for (Int_t b=1; b<=  fHistSpectrumStatRatio[m][iregion]->GetNbinsX(); b++){
	    bcorr = fHistSpectrumStat[m][0][iregion]->FindBin(fHistSpectrumStat[m][1][iregion]->GetBinCenter(b));	
	    cout << "\n iregion " << iregion << " m: " << m << " low edge bin: " <<	    fHistSpectrumStatRatio[m][iregion]->GetXaxis()->GetBinLowEdge(b)<< endl;
	    cout << " num info " << fHistSpectrumStat[m][1][iregion]->GetBinContent(b) <<endl;
	    cout << " denum info " << fHistSpectrumStat[m][0][iregion]->GetBinContent(bcorr) << " or spline " << splineK0s[m][iregion]->Eval(fHistSpectrumSist[m][1][iregion]->GetBinCenter(b))<<endl;
	    cout << fHistSpectrumStatRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumStatRatio[m][iregion]->GetBinError(b)<< endl;
	    cout << fHistSpectrumSistRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumSistRatio[m][iregion]->GetBinError(b)<< endl;
	  }

	  fHistSpectrumStatDoubleRatio[m][iregion]=(TH1F*) fHistSpectrumStatRatio[m][iregion]->Clone("fHistSpectrumStatDoubleRatio_"+Smolt[m]+Form("_region%i", iregion));
	  fHistSpectrumSistDoubleRatio[m][iregion]=(TH1F*) fHistSpectrumSistRatio[m][iregion]->Clone("fHistSpectrumSistDoubleRatio_"+Smolt[m]+Form("_region%i", iregion));
	  
	  fHistSpectrumSistDoubleRatio[m][iregion]->Sumw2();
	  fHistSpectrumSistDoubleRatio[m][iregion]->Divide(	  fHistSpectrumSistRatio[5][iregion]);
	  fHistSpectrumStatDoubleRatio[m][iregion]->Sumw2();
	  fHistSpectrumStatDoubleRatio[m][iregion]->Divide(	  fHistSpectrumStatRatio[5][iregion]);

	  if (iregion==0){
	  for (Int_t b=1; b<=  fHistSpectrumStatDoubleRatio[m][iregion]->GetNbinsX(); b++){
	    cout << "\n iregion " << iregion << " m: " << m << " low edge bin: " <<	    fHistSpectrumStatDoubleRatio[m][iregion]->GetXaxis()->GetBinLowEdge(b)<< endl;
	    cout << fHistSpectrumStatDoubleRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumStatDoubleRatio[m][iregion]->GetBinError(b)<< endl;
	    cout << fHistSpectrumSistDoubleRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumSistDoubleRatio[m][iregion]->GetBinError(b)<< endl;
	  }
	  }
	  canvasPtSpectraRatio->cd(iregion+4);
	  gPad->SetLeftMargin(0.15);
	  //	  StyleHisto(fHistSpectrumStatDoubleRatio[m][iregion], 0.0001, 2,  ColorMult[m], 33, titleX,  TitleYPtDRatio,  TitleYYieldRatio + " spectra ratio " + RegionType[iregion]);
	  StyleHisto(fHistSpectrumStatDoubleRatio[m][iregion], 0.0001, 2,  ColorMult[m], 33, titleX, "", "");
	  //	  StyleHisto(fHistSpectrumSistDoubleRatio[m][iregion], 0.0001, 2,  ColorMult[m], 33, titleX,  TitleYPtDRatio,  TitleYYieldRatio + " spectra ratio " + RegionType[iregion]);
	  StyleHisto(fHistSpectrumSistDoubleRatio[m][iregion], 0.0001, 2,  ColorMult[m], 33, titleX, "", "");
	  fHistSpectrumStatDoubleRatio[m][iregion]->SetMarkerSize(1.2);
	  fHistSpectrumSistDoubleRatio[m][iregion]->SetMarkerSize(1.2);
	  //	  if (m!=nummolt){
	  if (m==0 || m==4){
	  fHistSpectrumStatDoubleRatio[m][iregion]->Draw("same p");
	  fHistSpectrumSistDoubleRatio[m][iregion]->SetFillStyle(0);
	  fHistSpectrumSistDoubleRatio[m][iregion]->Draw("same e2");
	  //	  fHistSpectrumSistDoubleRatio[m][iregion]->Draw("same p hist");
	  lineat1->Draw("same");
	  }
	}//end type 0 else
      }//end loop on mult


      //drawing spectra
      canvasPtSpectra[iregion][type]->cd();
      gPad->SetLogy();
      canvasPtSpectra[iregion][type]->SetFillColor(0);
      canvasPtSpectra[iregion][type]->SetTickx(1);
      canvasPtSpectra[iregion][type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.2);
      gPad->SetBottomMargin(0.13);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);

      TLegend *Legend1Bis=new TLegend(0.3,0.8,0.92,0.92);
      Legend1Bis->SetFillStyle(0);
      Legend1Bis->SetTextAlign(32);
      Legend1Bis->AddEntry("", "#bf{ALICE Preliminary}", "");
      Legend1Bis->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
      Legend1Bis->AddEntry("", NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
      //      Legend1Bis->AddEntry("", sRegionBlack[iregion], "");

      TLegend *Legend1NoPrel[2]; 
      for (Int_t type=0; type<=1; type++){
	Legend1NoPrel[type]=new TLegend(0.3,0.83,0.92,0.92);
	Legend1NoPrel[type]->SetFillStyle(0);
	Legend1NoPrel[type]->SetTextAlign(32);
	Legend1NoPrel[type]->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
	Legend1NoPrel[type]->AddEntry("", NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
	}

      TLegend *legendRegionBlack=new TLegend(0.65, 0.7, 0.92, 0.78);
      TLegendEntry* lReAll1b;
      TLegendEntry* lReAll2b;
      lReAll1b = legendRegionBlack->AddEntry("", sRegionBlack[iregion], "");
      lReAll1b->SetTextSize(0.036);
      lReAll1b->SetTextAlign(32);
      lReAll2b=      legendRegionBlack->AddEntry("", sRegion1[iregion], "");
      lReAll2b->SetTextSize(0.029);
      lReAll2b->SetTextAlign(32);
      

      fHistSpectrumSistScaledB[0][type][iregion]= (TH1F*) 	fHistSpectrumSistScaled[0][type][iregion]->Clone("sistblack");
      fHistSpectrumSistScaledB[0][type][iregion]->SetLineColor(1);
      fHistSpectrumSistScaledB[0][type][iregion]->SetMarkerColor(1);
      fHistSpectrumSistScaledB[0][type][iregion]->SetMarkerStyle(20);
      TLegend *legendStatBoxK0s=new TLegend(0.25, 0.33, 0.47, 0.42);
      legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[0][type][iregion], "stat. error", "pe");
      legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[0][type][iregion], "syst. error", "ef");

      TLegend *legendStatBoxXi=new TLegend(0.23, 0.83, 0.42, 0.92);
      legendStatBoxXi->AddEntry(fHistSpectrumSistScaledB[0][type][iregion], "stat. error", "pe");
      legendStatBoxXi->AddEntry(fHistSpectrumSistScaledB[0][type][iregion], "syst. error", "ef");

      TLegend *LegendMolt=new TLegend(0.25,0.16,0.9,0.31);
      LegendMolt->SetNColumns(3);
      LegendMolt->SetFillStyle(0);
      LegendMolt->SetHeader("V0M Multiplicity Percentile");
      TLegendEntry *lheader = (TLegendEntry*)LegendMolt->GetListOfPrimitives()->First();
      //      TLegendEntry * lheader = (TLegendEntry*)LegendMolt->SetHeader("V0M Multiplicity Percentile");
      lheader-> SetTextSize(0.033);
      //LegendMolt->SetTextAlign(32);

      for (Int_t m=0; m<nummolt+1; m++){

	fHistSpectrumStatScaled[m][type][iregion]->Scale(ScaleFactor[m]);
	fHistSpectrumSistScaled[m][type][iregion]->Scale(ScaleFactor[m]);
	StyleHistoYield(fHistSpectrumStatScaled[m][type][iregion], LimInfSpectra[type][iregion], LimSupSpectra[type][iregion], ColorMult[m], MarkerMult[m], titleX, titleY,"", size[m],1.15, YoffsetSpectra[type] );
	StyleHistoYield(fHistSpectrumSistScaled[m][type][iregion], LimInfSpectra[type][iregion], LimSupSpectra[type][iregion], ColorMult[m], MarkerMult[m], titleX, titleY,"", size[m], 1.15,  YoffsetSpectra[type]);

	LegendMolt->AddEntry(fHistSpectrumStatScaled[m][type][iregion],SmoltBis[m]+"%"+sScaleFactor[m]+" ","pef");
      fHistSpectrumStatScaled[m][type][iregion]->Draw("same e0x0");  
      fHistSpectrumSistScaled[m][type][iregion]->SetFillStyle(0);
      fHistSpectrumSistScaled[m][type][iregion]->Draw("same e2");
      }

      //LegendMolt->AddEntry("","Uncertainties","");
      LegendMolt->Draw("");
      Legend1Bis->Draw("");
      legendRegionBlack->Draw("");
      if(type==0)legendStatBoxK0s->Draw("");
      else      if(type==1)legendStatBoxXi->Draw("");

      canvasPtSpectraMultRatio[iregion][type]->cd();
      canvasPtSpectraMultRatio[iregion][type]->SetFillColor(0);
      canvasPtSpectraMultRatio[iregion][type]->SetTickx(1);
      canvasPtSpectraMultRatio[iregion][type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.05);
      gPad->SetBottomMargin(0.13);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);

      TLine * lineat1Mult = new TLine(0, 1,8,1);
      lineat1Mult->SetLineColor(1);
      lineat1Mult->SetLineStyle(2);

      fHistSpectrumStatMultRatio[nummolt][type][iregion]= (TH1F*) fHistSpectrumStat[nummolt][type][iregion]->Clone("fHistSpectrumStatMultRatio_"+Smolt[nummolt]+Form("_type%i_region%i", type,iregion));
      fHistSpectrumSistMultRatio[nummolt][type][iregion]= (TH1F*) fHistSpectrumSist[nummolt][type][iregion]->Clone("fHistSpectrumSistMultRatio_"+Smolt[nummolt]+Form("_type%i_region%i", type,iregion));
      for (Int_t m=0; m<nummolt; m++){
	fHistSpectrumStatMultRatio[m][type][iregion]= (TH1F*) fHistSpectrumStat[m][type][iregion]->Clone("fHistSpectrumStatMultRatio_"+Smolt[m]+Form("_type%i_region%i", type,iregion));
	fHistSpectrumSistMultRatio[m][type][iregion]= (TH1F*) fHistSpectrumSist[m][type][iregion]->Clone("fHistSpectrumSistMultRatio_"+Smolt[m]+Form("_type%i_region%i", type,iregion));
	fHistSpectrumStatMultRatio[m][type][iregion]->Divide(	fHistSpectrumStatMultRatio[nummolt][type][iregion]);
	fHistSpectrumSistMultRatio[m][type][iregion]->Divide(	fHistSpectrumSistMultRatio[nummolt][type][iregion]);
	StyleHistoYield(fHistSpectrumStatMultRatio[m][type][iregion], 10e-5, LimSupMultRatio[type], ColorMult[m], MarkerMult[m], titleX, "Ratio to 0-100% class","", size[m],1.15, YoffsetSpectra[type] );
	StyleHistoYield(fHistSpectrumSistMultRatio[m][type][iregion], 10e-5, LimSupMultRatio[type], ColorMult[m], MarkerMult[m], titleX,  "Ratio to 0-100% class","", size[m], 1.15,  YoffsetSpectra[type]);

      fHistSpectrumStatMultRatio[m][type][iregion]->Draw("same e0x0");  
      fHistSpectrumSistMultRatio[m][type][iregion]->SetFillStyle(0);
      fHistSpectrumSistMultRatio[m][type][iregion]->Draw("same e2");
      lineat1Mult->Draw("same");
      }
      //      legendMult->Draw("");
      Legend1NoPrel[type]->Draw("");
      legendRegionBlack->Draw("");
      //if(type==0)legendStatBoxK0s->Draw("");
      //      else      if(type==1)legendStatBoxXi->Draw("");

      
      //Drawing pt vs mult altogether
      canvasPtvsMult[type]->cd();
      //      gPad->SetLogy();
      canvasPtvsMult[type]->SetFillColor(0);
      canvasPtvsMult[type]->SetTickx(1);
      canvasPtvsMult[type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);

      StyleHistoYield(fHistPtvsMultStat[type][iregion], LimInfPt[type]+10e-7, LimSupPt[type], Color[iregion], YmarkerRegion[iregion], titleYieldX, titlePtvsMultYType[type],"",  MarkerSize[iregion], 1.2 , 1.25);
      StyleHistoYield(fHistPtvsMultSist[type][iregion], LimInfPt[type]+10e-7, LimSupPt[type], Color[iregion], YmarkerRegion[iregion], titleYieldX, titlePtvsMultYType[type], "",  MarkerSize[iregion], 1.2, 1.25);
      StyleHistoYield(fHistPtvsMultStatFS[type][iregion], LimInfPt[type]+10e-7, LimSupPt[type], Color[iregion], YmarkerRegionBis[iregion], titleYieldX, titlePtvsMultYType[type],"",  MarkerSize[iregion], 1.2 , 1.25);
      StyleHistoYield(fHistPtvsMultSistFS[type][iregion], LimInfPt[type]+10e-7, LimSupPt[type], Color[iregion], YmarkerRegionBis[iregion], titleYieldX, titlePtvsMultYType[type], "",  MarkerSize[iregion], 1.2, 1.25);

      fHistPtvsMultStat[type][iregion]->DrawClone("same e0x0");  
      fHistPtvsMultSist[type][iregion]->SetFillStyle(0);
      fHistPtvsMultSist[type][iregion]->DrawClone("same e2");
      fHistPtvsMultStatFS[type][iregion]->DrawClone("same e0x0");  
      fHistPtvsMultSistFS[type][iregion]->SetFillStyle(0);
      fHistPtvsMultSistFS[type][iregion]->DrawClone("same e2");


      canvasPtvsMultMethodComp[type]->cd();
      canvasPtvsMultMethodComp[type]->SetFillColor(0);
      canvasPtvsMultMethodComp[type]->SetTickx(1);
      canvasPtvsMultMethodComp[type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);

      StyleHistoYield(fHistPtvsMultStatMethodComp[type][iregion], 0.7+10e-7, 1.4, Color[iregion], YmarkerRegion[iregion], titleYieldX, "<p_{T}> ratio (Fit/Spectrum)","",  MarkerSize[iregion], 1.2 , 1.25);
      StyleHistoYield(fHistPtvsMultSistMethodComp[type][iregion], 0.7+10e-7, 1.4, Color[iregion], YmarkerRegion[iregion], titleYieldX, "<p_{T}> ratio (Fit/Spectrum)", "",  MarkerSize[iregion], 1.2, 1.25);

      fHistPtvsMultStatMethodComp[type][iregion]->Divide(fHistPtvsMultStat[type][iregion]);
      fHistPtvsMultSistMethodComp[type][iregion]->Divide(fHistPtvsMultSist[type][iregion]);
      ErrRatioCorr(fHistPtvsMultStatFS[type][iregion],fHistPtvsMultStat[type][iregion],  fHistPtvsMultStatMethodComp[type][iregion], 1);
      ErrRatioCorr(fHistPtvsMultSistFS[type][iregion],fHistPtvsMultSist[type][iregion],  fHistPtvsMultSistMethodComp[type][iregion], 1);
      fHistPtvsMultStatMethodComp[type][iregion]->DrawClone("same e0x0");  
      fHistPtvsMultSistMethodComp[type][iregion]->SetFillStyle(0);
      fHistPtvsMultSistMethodComp[type][iregion]->DrawClone("same e2");

      //      Legend2->Draw("");
      //      legendStatBox->Draw("");
      //      legendRegionAllF->Draw("");

    }//end loop region

    canvasYield[type]->cd();
    legendRegion->Draw("");
    legendYield->Draw("");

    if(type==0){
    canvasYieldRatio->cd();
    //    legendRegion->Draw("");
    //    legendYield->Draw("");
    }
    //    fileout->WriteTObject(canvasYield[type]);
    //    fileout->WriteTObject(canvasYieldFinal[type]);
  }//end loop type

  fileout->WriteTObject(canvasYieldJet);
  fileout->WriteTObject(canvasYieldRatio);
  fileout->WriteTObject(canvasPtSpectraRatioJetPY);
  fileout->WriteTObject(canvasPtSpectraRatio);
  fileout->WriteTObject(canvasPtSpectraRatioAllMult);
  fileout->WriteTObject(canvasPtSpectraK0s);

  TString DirPicture = "PictureForNote/";
  if (!isNormCorr) DirPicture = "PictureForNote/IsNotNormCorr_";

  canvasPtSpectraRatioJetPY->SaveAs(DirPicture+"PtSpectraRatioK0sXi_PengYaoComp.pdf");  	
  canvasPtSpectraRatio->SaveAs(DirPicture+"PtSpectraRatioK0sXi.pdf");  	
  canvasPtSpectraRatioAllMult->SaveAs(DirPicture+"PtSpectraRatioK0sXiAllMult.pdf");  	
  canvasYieldJet->SaveAs(DirPicture+"InJetYieldvsMultK0sXi.pdf");  	
  if (ZeroYieldLowPt){
  canvasYieldRatio->SaveAs(DirPicture+"YieldRatiovsMult_ZeroYieldLowPtJet.pdf");  	
  }
  else {
  canvasYieldRatio->SaveAs(DirPicture+"YieldRatiovsMult.eps");  	
  canvasYieldRatio->SaveAs(DirPicture+"YieldRatiovsMult.pdf");  	
  }
  //  canvasYield[0]->SaveAs(DirPicture+"AllYieldvsMultK0s.pdf");  	
  // canvasYield[1]->SaveAs(DirPicture+"AllYieldvsMultXi.pdf");  	
  canvasYieldFinal[0]->SaveAs(DirPicture+"AllYieldvsMultK0s.eps");  	
  canvasYieldFinal[0]->SaveAs(DirPicture+"AllYieldvsMultK0s.pdf");  	
  canvasYieldFinal[1]->SaveAs(DirPicture+"AllYieldvsMultXi.eps");  	
  canvasYieldFinal[1]->SaveAs(DirPicture+"AllYieldvsMultXi.pdf");  	

  fileout->WriteTObject(  canvasYieldFinal[0]);
  fileout->WriteTObject(  canvasYieldFinal[1]);

  for (Int_t type=0; type<2; type++){
    for (Int_t iregion=0; iregion<3; iregion++){
      canvasYieldSeparate[iregion][type]->SaveAs(DirPicture+"FinalYieldvsMult"+tipo[type]+RegionType[iregion]+".eps");
      canvasPtSpectra[iregion][type]->SaveAs(DirPicture+"FinalPtSpectra"+tipo[type]+RegionType[iregion]+".eps");
      canvasYieldErrSeparate[iregion][type]->SaveAs(DirPicture+"FinalYieldvsMultErr"+tipo[type]+RegionType[iregion]+".eps");

      canvasYieldSeparate[iregion][type]->SaveAs(DirPicture+"FinalYieldvsMult"+tipo[type]+RegionType[iregion]+".pdf");
      canvasPtSpectra[iregion][type]->SaveAs(DirPicture+"FinalPtSpectra"+tipo[type]+RegionType[iregion]+".pdf");
      canvasPtSpectraMultRatio[iregion][type]->SaveAs(DirPicture+"FinalPtSpectraMultRatio"+tipo[type]+RegionType[iregion]+".pdf");
      canvasYieldErrSeparate[iregion][type]->SaveAs(DirPicture+"FinalYieldvsMultErr"+tipo[type]+RegionType[iregion]+".pdf");
      canvasPtvsMultErrSeparate[iregion][type]->SaveAs(DirPicture+"FinalPtvsMultvsMultErr"+tipo[type]+RegionType[iregion]+".pdf");

      fileout->WriteTObject(canvasYieldSeparate[iregion][type]);
      fileout->WriteTObject(canvasPtSpectra[iregion][type]);
      fileout->WriteTObject(canvasPtSpectraMultRatio[iregion][type]);
      fileout->WriteTObject(canvasYieldErrSeparate[iregion][type]);
      fileout->WriteTObject(canvasPtvsMultErrSeparate[iregion][type]);
    }
    canvasPtvsMult[type]->SaveAs(DirPicture+"FinalPtvsMult"+tipo[type]+".pdf");
    canvasPtvsMultMethodComp[type]->SaveAs(DirPicture+"FinalPtvsMultMethodComp"+tipo[type]+".pdf");
    fileout->WriteTObject(canvasPtvsMult[type]);
  }

  fileout->Close();

  cout << "\n\n**fit to Xi/K0s with pol0 and pol1" << endl;
  cout << "\nChi/NDF pol0 fit :" << endl;
  for (Int_t iregion=0; iregion<numregions; iregion++){
    cout << RegionType[iregion] << " " <<  pol0Ratio[iregion]->GetChisquare() << "/"<< pol0Ratio[iregion]->GetNDF() << endl;
  } 
  cout << "\nChi/NDF pol1 fit :" << endl;
  for (Int_t iregion=0; iregion<numregions; iregion++){
    cout << RegionType[iregion] << " " <<  pol1Ratio[iregion]->GetChisquare() << "/"<< pol1Ratio[iregion]->GetNDF() << endl;
  } 

  for (Int_t type=0; type<2; type++){
    cout << "\n\nstarting from the file(s) "  <<  PathInYield[type] << endl; //" and " <<  PathInDPhiProj[type] << endl;
  }
  cout << " I have created the file " << stringout << endl;
}

