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
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include <Macros/constants.h>

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
//TString titledNdetaTrigg="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}^{#it{p}_{T, trigg}>3 GeV/#it{c}}";
TString titledNdetaTrigg="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5, #it{p}_{T,trigg}>3 GeV/#it{c}}";
TString titleYield[2]={"K^{0}_{S}", "#Xi"};
TString TitleYYieldRatio="#it{N}_{#Xi}/#it{N}_{K^{0}_{S}}";

//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}"};
TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})", "#it{N}_{#Xi}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})"};
//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)"}; //to much space around the parenthesis
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
//TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
//TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};
//TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
TString sRegion1[3]={"|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 1.1", "0.86 < |#Delta#it{#eta}| < 1.2, 1.1 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};

const Int_t numMolt = 8; //(5 MB points + 3 HM points)

void RatiosXiK0s3Systems( Int_t PlotType =0, Int_t ChosenRegion = -1,  Float_t ScalingFactorXiK0s = 1/*0.8458/1.08747*/, Int_t isPreliminary =0, Bool_t isChangesIncluded=1, Bool_t isFit=0, Bool_t isdNdEtaTriggered=1, Bool_t MaterialBudgetCorr =1){

  //NB: "ChangesIncluded" labels the new file produced by choosing Sidebands for K0s and Xi in HM, !SkipAssoc, and larger dEta choice for K0s

  Bool_t isMultCorrEval = 1; //yield with sist uncertainties uncorrelated with multiplicity are taken in input
  if (PlotType>=3) isMultCorrEval = 0;
  else isMultCorrEval = 1;

  Bool_t isWingsCorrectionApplied = 0;
  if (PlotType==0 || PlotType==1 || PlotType==3){
    cout <<"Do you want to analyse the files with the wings correction applied? Type 1 if you DO want" << endl;
    cin >> isWingsCorrectionApplied;
  }

  //isPreliminary = 1: preliminary plots for RATIO (not corrected by norm factor)
  //isPreliminary = 2: preliminary plots for YIELDS (corrected by norm factor)

  if (PlotType==0 && isPreliminary==2) return;
  if (PlotType>0 && isPreliminary==1) return;

  //PlotType = 0: Xi/K0s ratio
  //PlotType = 1: K0s yield vs mult 
  //PlotType = 2: Xi yield vs mult 
  //PlotType = 3: K0s pt vs mult 
  //PlotType = 4: Xi pt vs mult 

  if (PlotType==1 || PlotType ==3) {
    //    sRegion1[0] = "|#Delta#it{#eta}| < 0.85, |#Delta#it{#varphi}| < 1.09";
    //    sRegion1[1] = "0.85 < |#Delta#it{#eta}| < 1.2, 1.09 < #Delta#it{#varphi} <1.8";
  }

  Float_t UpperValueX= 45;

  //Set titles
  TString titleY;
  if (PlotType==0) titleY = TitleYYieldRatio;
  else if (PlotType==1) titleY = titleYieldYType[0];
  else if (PlotType==2) titleY = titleYieldYType[1];
  else if (PlotType==3) titleY = titlePtvsMultYType[0];
  else if (PlotType==4) titleY = titlePtvsMultYType[1];

  TString  titleX = titledNdetaTrigg;
  Float_t xOffset =1.2;
  Float_t yOffset =1.25;
  //  Float_t MarkerSize[3] ={2, 2, 3};
  Float_t MarkerSize[3] ={1.5, 2, 3};
  Int_t MarkerType[3] = {20, 21, 33};

  //  TString titleYToOOJ = "#Xi/K^{0}_{S} (JET) / #Xi/K^{0}_{S} (OOJ)";
  //  TString titleYToOOJ = "Near-side jet / out-of-jet";
  TString titleYToOOJ = "Toward / Transverse";

  //ChosenRegion == -1 //all regions (JET, OOJ, FULL) are plotted
  if (ChosenRegion>2) {cout << "Choose a valid region! " << endl; return;}
  const Int_t numColls =2; //pp13 TeV (MB + HM), pp5TeV
  const Int_t numRegions =3;
  const Int_t numTypes =2;
  TString CollisionsComp[2] = {"_HMMultBinning1_vsHM", "_vs5TeV"};
  TString SPlotType[5] = {"Ratio", "K0s", "Xi", "K0spt", "Xipt"};
  TString SisPreliminary[3] = {"", "_isPreliminaryForRatio", "_isPreliminaryForYields"};
  TString SisFit[2] = {"", "_isFit"};

  gStyle->SetOptStat(0);
  TString Region[numRegions] = {"Jet", "Bulk", "All"};
  TString RegionBis[numRegions] = {"Jet", "Bulk", "Inclusive"};
  Float_t Up= 0.14-10e-6;
  Float_t Low = 0.02+10e-6;
  if (PlotType==1) {
    Low = 10e-6; Up = 0.50-10e-6; //0.45
    if (ChosenRegion==0) {Low = 0.015+10e-6; Up = 0.06-10e-6;} //0.035 
  }
  if (PlotType==2) {
    Low = 10e-6; Up = 0.04-10e-6;
    if (ChosenRegion==0) {Up = 0.0045-10e-6;}
  }
  if (PlotType==3) {
    Low = 0.8+10e-4; Up = 2.4999;
    if (ChosenRegion==0) {Low = 1.5; Up = 2.5;}
    else if (ChosenRegion>0) {Low = 0.8; Up = 1.3;}
  }
  if (PlotType==4) {
    Low = 1+10e-4; Up = 4.999;
    if (ChosenRegion==0) {Low = 1.5; Up = 5;}
    else if (ChosenRegion>0) {Low = 1.0; Up = 2.2;}
  }
  Int_t Color[numRegions] = {628,418,601};
  Int_t ColorEnergy[numColls] = {922,920};
  Int_t ColorDiff[numRegions][numColls] = {{634,628}, {418,829} , {601, 867}};
  Int_t Marker[numTypes] = {33, 21};
  Float_t Size[numTypes] = {2, 1.4};
  TString ParticleType[numTypes]={"K0s", "Xi"};

  TF1* fitPol1[3];
  TF1* fitPol0[3];

  TF1 * fitToRatioToOOJ = new TF1("fitToRatioToOOJ", "pol0", 0, UpperValueX);
  fitToRatioToOOJ->SetLineColor(1);
  fitToRatioToOOJ->SetLineStyle(2);
  fitToRatioToOOJ->SetLineWidth(2);

  TH1F *histoYield[numTypes][numRegions][numColls];
  TH1F *histoYieldAllErr[numTypes][numRegions][numColls];
  TH1F *histoYieldRatio[numRegions][numColls];
  TH1F *histoYieldRatioAllErr[numRegions][numColls]; //stat + syst in quadrature
  TH1F *histoYieldRatioToOOJ[numColls];
  TH1F *histoYieldSist[numTypes][numRegions][numColls];
  TH1F *histoYieldSistMultUnCorr[numTypes][numRegions][numColls];
  TH1F *histoYieldRatioSist[numRegions][numColls];
  TH1F *histoYieldRatioSistMultUnCorr[numRegions][numColls];
  TH1F *histoYieldRatioSistToOOJ[numColls];
  TH1F *histoYieldRatioSistMultUnCorrToOOJ[numColls];
  TH1F*  fHistYieldStatBlack;
  TH1F*  fHistYieldSistBlack;
  TH1F*  fHistYieldSistMultUnCorrBlack;
  TH1F*  fHistYieldStatGrey[numColls];

  TGraphAsymmErrors *gYieldRatio[numRegions][numColls];
  TGraphAsymmErrors *gYieldRatioToOOJ[numColls];
  TGraphAsymmErrors *gYieldRatioSist[numRegions][numColls];
  TGraphAsymmErrors *gYieldRatioToOOJSist[numColls];
  TGraphAsymmErrors *gYieldRatioSistMultUnCorr[numRegions][numColls];
  TGraphAsymmErrors *gYieldRatioToOOJSistMultUnCorr[numColls];
  Float_t  gY[numRegions][numColls][numMolt] = {0};
  Float_t  gYError[numRegions][numColls][numMolt] = {0};
  Float_t  gYErrorSist[numRegions][numColls][numMolt] = {0};
  Float_t  gYErrorSistMultUnCorr[numRegions][numColls][numMolt] = {0};
  Float_t  gYRatioToOOJ[numColls][numMolt] = {0};
  Float_t  gYRatioToOOJError[numColls][numMolt] = {0};
  Float_t  gYRatioToOOJErrorSist[numColls][numMolt] = {0};
  Float_t  gYRatioToOOJErrorSistMultUnCorr[numColls][numMolt] = {0};

  //DUMMY HISTO FOR AXES
  TH1F*histoYieldDummy= new TH1F("histoYieldDummy", "histoYieldDummy", 100, 0, UpperValueX);
  StyleHisto(histoYieldDummy, Low, Up, 1, 1, titleX, titleY, "" , 1,0, UpperValueX, xOffset, yOffset, 1);
  if (PlotType==0) {
    histoYieldDummy->GetYaxis()->SetTitleSize(0.08);//0.07
    histoYieldDummy->GetYaxis()->SetTitleOffset(0.7);//0.7
    histoYieldDummy->GetYaxis()->SetTickLength(0.02);
  }

  TH1F* hdummyRatioToOOJ = new TH1F ("hdummyRatioToOOJ", "hdummyRatioToOOJ",  1000, 0, UpperValueX);
  StyleHisto(hdummyRatioToOOJ, 0.42, 0.92, 1, 1, titleX, titleYToOOJ, "" , 1, 0, UpperValueX, xOffset, yOffset, 1);
  hdummyRatioToOOJ->GetYaxis()->SetTitleSize(0.07);
  hdummyRatioToOOJ->GetYaxis()->SetTitleOffset(0.7);
  hdummyRatioToOOJ->GetYaxis()->SetLabelSize(0.10);
  hdummyRatioToOOJ->GetYaxis()->SetLabelOffset(0.009);
  hdummyRatioToOOJ->GetXaxis()->SetTitleSize(0.10);
  hdummyRatioToOOJ->GetXaxis()->SetTitleOffset(1.3);
  hdummyRatioToOOJ->GetXaxis()->SetLabelSize(0.10);
  hdummyRatioToOOJ->GetXaxis()->SetLabelOffset(0.02);
  hdummyRatioToOOJ->GetYaxis()->SetNdivisions(305);
  hdummyRatioToOOJ->GetYaxis()->SetTickLength(0.02);

  Float_t dNdEta[numColls][numMolt] = {0};
  Float_t dNdEtaErrorL[numColls][numMolt] = {0};
  Float_t dNdEtaErrorR[numColls][numMolt] = {0};
  for (Int_t m=0; m<numMolt; m++){
    dNdEta[0][m] =  dNdEtaFinal13TeV[m];
    dNdEtaErrorL[0][m] =  dNdEtaFinal13TeV_ErrorL[m];
    dNdEtaErrorR[0][m] =  dNdEtaFinal13TeV_ErrorR[m];
    dNdEta[1][m] =  dNdEtaFinal5TeV[m];
    dNdEtaErrorL[1][m] =  dNdEtaFinal5TeV_ErrorL[m];
    dNdEtaErrorR[1][m] =  dNdEtaFinal5TeV_ErrorR[m];
  }

  TString NameHisto13TeV="histoYieldComparison";
  TString NameHistoSist13TeV="histoYieldSistComparison";
  TString NameHistoSistMultUnCorr13TeV="histoYieldSistMultUnCorrComparison";
  TString NameHisto5TeV="fHistYield_pp5TeV";
  TString NameHistoSist5TeV="fHistYield_pp5TeV_Sist";
  TString NameHistoSistMultUnCorr5TeV="fHistYield_pp5TeV_SistMultUnCorr";
  TString NameHistoFinal[numTypes][numRegions][numColls];
  TString pathin ="";
  TString YieldOrAvgPt[5] = {"Yield", "Yield", "Yield", "AvgPt","AvgPt"};
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  TCanvas * canvasRatioToOOJ = new TCanvas("canvasRatioToOOJ", "canvasRatioToOOJ", 1300, 800);
  TCanvas * canvasSplit = new TCanvas("canvasSplit", "canvasSplit", 1300, 1100);
  TPad* pad1 = new TPad( "pad1" ,"pad1" ,0 ,0.36 ,1 ,1);
  TPad* pad2 = new TPad( "pad2" ,"pad2" ,0 ,0.01 ,1 ,0.35);
  canvasSplit->SetFillColor(0);

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
  TLegendEntry* E1Ratio =      LegendRatio->AddEntry("", "#bf{This work}", "");
  TLegendEntry* E3Ratio =      LegendRatio->AddEntry(""/*(TObject*)0*/, "#it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  E1Ratio->SetTextSize(0.05);
  E3Ratio->SetTextSize(0.05);
  E1Ratio->SetTextAlign(32);
  E3Ratio->SetTextAlign(32);

  TLegend *LegendColor=new TLegend(0.16,0.80,0.5,0.92); 
  LegendColor->SetMargin(0);
  //  LegendColor->AddEntry("", "#bf{ALICE Preliminary}", "");
  LegendColor->AddEntry("", "#bf{This work}", "");
  //LegendColor->AddEntry("", "", "");

  TLegend *LegendYields=new TLegend(0.16,0.75,0.5,0.93);
  LegendYields->SetMargin(0);
  //  LegendYields->AddEntry("", "#bf{ALICE Preliminary}", "");
  LegendYields->AddEntry("", "#bf{This work}", "");
  //LegendYields->AddEntry("", "", "");
  //  LegendYields->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");

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

  TLegend *legendEnergyBoxBis1=new TLegend(0.66, 0.69, 0.89, 0.76);
  TLegend *legendEnergyBoxBis2=new TLegend(0.66, 0.62, 0.89, 0.69);

  TLegend *legendStatBox=new TLegend(0.73, 0.79, 0.93, 0.93);
  TLegend *legendStatBoxBis=new TLegend(0.2, 0.42, 0.35, 0.58);
  TLegend *legendStatBoxColor=new TLegend(0.18, 0.56, 0.35, 0.71);

  Int_t NLoopType =-1;
  Int_t NLoopRegion =-1;
  for (Int_t type=1; type>=0; type--){
    cout << "Type: " << type << endl;
    if ((PlotType == 1 || PlotType ==3) && type==1) continue;
    else  if ((PlotType == 2 || PlotType ==4) && type==0) continue;
    NLoopType++;
    LegendYields->AddEntry(""/*(TObject*)0*/, NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
    for (Int_t ireg=0; ireg<numRegions; ireg++){
      if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
      NLoopRegion++;
      cout << "\n\n*****" << ParticleType[type] << " " << Region[ireg]<< endl;
      for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
	pathin = "Compare" +YieldOrAvgPt[PlotType] + "DifferentCollisions";
	pathin +=  CollisionsComp[Coll];
	pathin +="_"+ ParticleType[type] + Region[ireg];
	if (isPreliminary ==1) pathin +=  "_isPreliminaryForRatio";
	else if (isPreliminary==2) pathin +=  "_isPreliminaryForYields";
	if (isChangesIncluded)  pathin +=  "_ChangesIncluded";
	if (isdNdEtaTriggered) pathin +="_isdNdEtaTriggered";
	if (type==0 && isPreliminary==0) pathin += "_EffCorr";
	if (isWingsCorrectionApplied) {
	  if (type==0)  pathin += "_WingsCorrApplied";
	}
	if (MaterialBudgetCorr){
	  pathin += "_MatBudgetCorrFAST";
	}
	if (isMultCorrEval) pathin += "_MultCorrSistEval";
	pathin +=  ".root";
	cout << "\n\e[35mPathin: " << pathin << "\e[39m"<< endl;
	TFile *filein = new TFile(pathin, "");
	if (!filein) {cout << "Input file not available " << endl; return;}

	NameHistoFinal[type][ireg][Coll]= Form("histoYield_Reg%i_Type%i_Coll%i", ireg, type, Coll);
	if (Coll==0)	histoYield[type][ireg][Coll]= (TH1F*) filein->Get(NameHisto13TeV);
	else 	histoYield[type][ireg][Coll]= (TH1F*) filein->Get(NameHisto5TeV);
	if (!histoYield[type][ireg][Coll]) {cout <<"no histo " << endl; return;}
	cout <<  histoYield[type][ireg][Coll]->GetXaxis()->GetXmax() << endl;
	histoYield[type][ireg][Coll]->SetName(NameHistoFinal[type][ireg][Coll]);

	histoYieldAllErr[type][ireg][Coll] = (TH1F*) histoYield[type][ireg][Coll]->Clone(NameHistoFinal[type][ireg][Coll]+ "_AllErr");

	if (Coll==0)	histoYieldSist[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoSist13TeV);
	else histoYieldSist[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoSist5TeV);
	if (!histoYieldSist[type][ireg][Coll]) {cout <<"no histo " << endl; return;}
	histoYieldSist[type][ireg][Coll]->SetName(NameHistoFinal[type][ireg][Coll]+"Sist");

	if (isMultCorrEval){
	  if (Coll==0)	histoYieldSistMultUnCorr[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoSistMultUnCorr13TeV);
	  else histoYieldSistMultUnCorr[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoSistMultUnCorr5TeV);
	  if (!histoYieldSistMultUnCorr[type][ireg][Coll]) {cout <<"no histo " << endl; return;}
	  histoYieldSistMultUnCorr[type][ireg][Coll]->SetName(NameHistoFinal[type][ireg][Coll]+"SistMultUnCorr");
	}

	cout << "Type: " << ParticleType[type] <<  " region: " << Region[ireg]<< endl;
	if (Coll ==0 ) cout << "Got histos for  13 TeV (MB + HM) " << endl;
	else  cout << "Got histos for 5 TeV " << endl;
	for (Int_t b=1; b<=  histoYield[type][ireg][Coll]->GetNbinsX(); b++){
	  if (histoYield[type][ireg][Coll]->GetBinContent(b) != 0)	 {
	    cout << histoYield[type][ireg][Coll]->GetBinContent(b) << " +- " << histoYield[type][ireg][Coll]->GetBinError(b) << " (stat.) +-  " << histoYieldSist[type][ireg][Coll]->GetBinError(b) << " (syst.) " << endl;
	  }
	}

	if ((type==1 && (PlotType == 0 || PlotType == 2 || PlotType == 4))  || (type==0 && (PlotType == 1 || PlotType == 3))) {
	  histoYieldRatio[ireg][Coll]= (TH1F*)      histoYield[type][ireg][Coll]->Clone(Form("Ratio_reg%i_Coll%i_Ratio",ireg, Coll));
	  histoYieldRatioSist[ireg][Coll]= (TH1F*)      histoYieldSist[type][ireg][Coll]->Clone(Form("Ratio_reg%i_Coll%i_RatioSist",ireg, Coll));
	  if (isMultCorrEval)  histoYieldRatioSistMultUnCorr[ireg][Coll]= (TH1F*)      histoYieldSistMultUnCorr[type][ireg][Coll]->Clone(Form("Ratio_reg%i_Coll%i_RatioSistMultUnCorr",ireg, Coll));
	  histoYieldRatioAllErr[ireg][Coll]= (TH1F*)      histoYieldSist[type][ireg][Coll]->Clone(Form("Ratio_reg%i_Coll%i_RatioAllErr",ireg, Coll));
	  histoYieldRatio[ireg][Coll]->Sumw2();
	  histoYieldRatioSist[ireg][Coll]->Sumw2();
	  if (isMultCorrEval) histoYieldRatioSistMultUnCorr[ireg][Coll]->Sumw2();
	  histoYieldRatioAllErr[ireg][Coll]->Sumw2();
	  for (Int_t i=1; i<= histoYieldRatioAllErr[ireg][Coll]->GetNbinsX(); i++){
	    if (histoYieldRatioAllErr[ireg][Coll]->GetBinContent(i) ==0) continue;
	    else histoYieldRatioAllErr[ireg][Coll]->SetBinError(i, sqrt(pow(histoYieldRatio[ireg][Coll]->GetBinError(i), 2) + pow(histoYieldRatioSist[ireg][Coll]->GetBinError(i),2) ));
	  }
	}
	else {
	  if (PlotType==0){
	    histoYieldRatio[ireg][Coll]->Divide(histoYield[type][ireg][Coll]);
	    histoYieldRatioSist[ireg][Coll]->Divide(histoYieldSist[type][ireg][Coll]);
	    if (isMultCorrEval)   histoYieldRatioSistMultUnCorr[ireg][Coll]->Divide(histoYieldSistMultUnCorr[type][ireg][Coll]);
	    histoYieldRatioAllErr[ireg][Coll]->Divide(histoYield[type][ireg][Coll]);
	    for (Int_t i=1; i<= histoYieldRatioAllErr[ireg][Coll]->GetNbinsX(); i++){
	      if (histoYieldRatioAllErr[ireg][Coll]->GetBinContent(i) ==0) continue;
	      else histoYieldRatioAllErr[ireg][Coll]->SetBinError(i, sqrt(pow(histoYieldRatio[ireg][Coll]->GetBinError(i), 2) + pow(histoYieldRatioSist[ireg][Coll]->GetBinError(i),2) ));
	    }
	    if (ireg ==0) {
	      histoYieldRatio[ireg][Coll]->Scale(ScalingFactorXiK0s);
	      histoYieldRatioSist[ireg][Coll]->Scale(ScalingFactorXiK0s);
	      if (isMultCorrEval)   histoYieldRatioSistMultUnCorr[ireg][Coll]->Scale(ScalingFactorXiK0s);
	      histoYieldRatioToOOJ[Coll]= (TH1F*)      histoYieldRatio[ireg][Coll]->Clone(Form("RatioToOOJ_Coll%i",Coll));
	      if (!histoYieldRatioToOOJ[Coll]) return;
	      histoYieldRatioSistToOOJ[Coll]= (TH1F*)      histoYieldRatioSist[ireg][Coll]->Clone(Form("RatioToOOJ_Coll%i_Sist", Coll));
	      if (!histoYieldRatioSistToOOJ[Coll]) return;
	      if (isMultCorrEval){
		histoYieldRatioSistMultUnCorrToOOJ[Coll]= (TH1F*)      histoYieldRatioSistMultUnCorr[ireg][Coll]->Clone(Form("RatioToOOJ_Coll%i_SistMultUnCorr", Coll));
		if (!histoYieldRatioSistMultUnCorrToOOJ[Coll]) return;
	      }
	    }
	    else if (ireg==1){
	      histoYieldRatioToOOJ[Coll]->Divide(histoYieldRatio[ireg][Coll]);
	      histoYieldRatioSistToOOJ[Coll]->Divide(histoYieldRatioSist[ireg][Coll]);
	      if (isMultCorrEval) histoYieldRatioSistMultUnCorrToOOJ[Coll]->Divide(histoYieldRatioSistMultUnCorr[ireg][Coll]);
	    }
	  }
	}

	if (!(type==1 && PlotType == 0)){
	  //SetValues on TGRAPH
	  for (Int_t m=0; m< numMolt; m++){
	    gY[ireg][Coll][m] = histoYieldRatio[ireg][Coll]->GetBinContent(histoYieldRatio[ireg][Coll]->FindBin(dNdEta[Coll][m]));
	    gYError[ireg][Coll][m] = histoYieldRatio[ireg][Coll]->GetBinError(histoYieldRatio[ireg][Coll]->FindBin(dNdEta[Coll][m]));
	    gYErrorSist[ireg][Coll][m] = histoYieldRatioSist[ireg][Coll]->GetBinError(histoYieldRatio[ireg][Coll]->FindBin(dNdEta[Coll][m]));
	    if (isMultCorrEval)	    gYErrorSistMultUnCorr[ireg][Coll][m] = histoYieldRatioSistMultUnCorr[ireg][Coll]->GetBinError(histoYieldRatio[ireg][Coll]->FindBin(dNdEta[Coll][m]));
	  }
	  gYieldRatio[ireg][Coll] = new TGraphAsymmErrors(numMolt,dNdEta[Coll],gY[ireg][Coll] , 0,0, gYError[ireg][Coll], gYError[ireg][Coll]);
	  gYieldRatioSist[ireg][Coll] = new TGraphAsymmErrors(numMolt,dNdEta[Coll],gY[ireg][Coll] ,dNdEtaErrorL[Coll],dNdEtaErrorR[Coll], gYErrorSist[ireg][Coll], gYErrorSist[ireg][Coll]);
	  if (isMultCorrEval) gYieldRatioSistMultUnCorr[ireg][Coll] = new TGraphAsymmErrors(numMolt,dNdEta[Coll],gY[ireg][Coll] ,dNdEtaErrorL[Coll],dNdEtaErrorR[Coll], gYErrorSistMultUnCorr[ireg][Coll], gYErrorSistMultUnCorr[ireg][Coll]);

	  StyleTGraphErrors(gYieldRatio[ireg][Coll], ColorDiff[ireg][Coll], MarkerType[ireg], MarkerSize[ireg], 1);
	  StyleTGraphErrors(gYieldRatioSist[ireg][Coll], ColorDiff[ireg][Coll], MarkerType[ireg], MarkerSize[ireg], 1);
	  if (isMultCorrEval) StyleTGraphErrors(gYieldRatioSistMultUnCorr[ireg][Coll], ColorDiff[ireg][Coll], MarkerType[ireg], MarkerSize[ireg], 1);

	  for (Int_t i =0; i<gYieldRatioSist[ireg][Coll]->GetN(); i++){
	    double x;
	    double y;
	    gYieldRatioSist[ireg][Coll]->GetPoint(i, x, y);
	    cout << " dNdeta " << x << " value " << y << endl;
	  }

	  if (PlotType==0 && ireg==1){
	    for (Int_t m=0; m< numMolt; m++){
	      gYRatioToOOJ[Coll][m] = histoYieldRatioToOOJ[Coll]->GetBinContent(histoYieldRatioToOOJ[Coll]->FindBin(dNdEta[Coll][m]));
	      gYRatioToOOJError[Coll][m] = histoYieldRatioToOOJ[Coll]->GetBinError(histoYieldRatioToOOJ[Coll]->FindBin(dNdEta[Coll][m]));
	      gYRatioToOOJErrorSist[Coll][m] = histoYieldRatioSistToOOJ[Coll]->GetBinError(histoYieldRatioToOOJ[Coll]->FindBin(dNdEta[Coll][m]));
	      if (isMultCorrEval) gYRatioToOOJErrorSistMultUnCorr[Coll][m] = histoYieldRatioSistMultUnCorrToOOJ[Coll]->GetBinError(histoYieldRatioToOOJ[Coll]->FindBin(dNdEta[Coll][m]));
	    }
	    gYieldRatioToOOJ[Coll] = new TGraphAsymmErrors(numMolt,dNdEta[Coll],gYRatioToOOJ[Coll] , 0, 0, gYRatioToOOJError[Coll], gYRatioToOOJError[Coll]);
	    gYieldRatioToOOJSist[Coll] = new TGraphAsymmErrors(numMolt,dNdEta[Coll],gYRatioToOOJ[Coll] ,dNdEtaErrorL[Coll],dNdEtaErrorR[Coll], gYRatioToOOJErrorSist[Coll], gYRatioToOOJErrorSist[Coll]);
	    if (isMultCorrEval)    gYieldRatioToOOJSistMultUnCorr[Coll] = new TGraphAsymmErrors(numMolt,dNdEta[Coll],gYRatioToOOJ[Coll] ,dNdEtaErrorL[Coll],dNdEtaErrorR[Coll], gYRatioToOOJErrorSistMultUnCorr[Coll], gYRatioToOOJErrorSistMultUnCorr[Coll]);

	    StyleTGraphErrors(gYieldRatioToOOJ[Coll], ColorDiff[0][Coll], MarkerType[0], MarkerSize[0], 1);
	    StyleTGraphErrors(gYieldRatioToOOJSist[Coll], ColorDiff[0][Coll], MarkerType[0], MarkerSize[0], 1);
	    if (isMultCorrEval) StyleTGraphErrors(gYieldRatioToOOJSistMultUnCorr[Coll], ColorDiff[0][Coll], MarkerType[0], MarkerSize[0], 1);
	  }
	}
	//End of SetValues on TGRAPH

	if (Coll==0 && isFit){
	  fitPol1[ireg] = new TF1(Form("pol1_%i", ireg), "pol1", 0, UpperValueX);
	  fitPol0[ireg] = new TF1(Form("pol0_%i", ireg), "pol0", 0, UpperValueX);
	  fitPol1[ireg]->SetLineColor(Color[ireg]);
	  fitPol0[ireg]->SetLineColor(Color[ireg]);
	  fitPol0[ireg]->SetLineStyle(10);
	  histoYieldRatioAllErr[ireg][Coll]->Fit(fitPol1[ireg], "R+");
	  histoYieldRatioAllErr[ireg][Coll]->Fit(fitPol0[ireg], "R+");
	  legendFit->AddEntry(fitPol0[ireg], Form("pol0 Chi2/NDF %.3f/%i", fitPol0[ireg]->GetChisquare(),fitPol0[ireg]->GetNDF()), "l" );
	  legendFit->AddEntry(fitPol1[ireg], Form("pol1 Chi2/NDF %.3f/%i", fitPol1[ireg]->GetChisquare(),fitPol1[ireg]->GetNDF()), "l" );
	}
	if (type==0){
	  cout << "\nXi/K0s ratio: " << endl;
	  for (Int_t b=1; b<=  histoYieldRatio[ireg][Coll]->GetNbinsX(); b++){
	    if (histoYieldRatio[ireg][Coll]->GetBinContent(b) != 0)	 {
	      cout << histoYieldRatio[ireg][Coll]->GetBinContent(b) << " +- " << histoYieldRatio[ireg][Coll]->GetBinError(b) << " (stat.) +-  " << histoYieldRatioSist[ireg][Coll]->GetBinError(b) << " (syst.) " << endl;
	    }
	  }
	}

	if (type==0 && ireg == 1 && PlotType==0){
	  canvasRatioToOOJ->cd();
	  canvasRatioToOOJ->SetFillColor(0);
	  canvasRatioToOOJ->SetTickx(1);
	  canvasRatioToOOJ->SetTicky(1);
	  gPad->SetTopMargin(0.05);
	  gPad->SetLeftMargin(0.13);
	  gPad->SetBottomMargin(0.15);
	  gPad->SetRightMargin(0.05);
	  gStyle->SetLegendBorderSize(0);
	  gStyle->SetLegendFillColor(0);
	  gStyle->SetLegendFont(42);
	  StyleHisto(histoYieldRatioToOOJ[Coll], 0.42, 0.92, ColorDiff[0][Coll], MarkerType[0], titleX, titleYToOOJ, "" , 1,0, 40, xOffset, yOffset, MarkerSize[0]);
	  StyleHisto(histoYieldRatioSistToOOJ[Coll], 0.42, 0.92, ColorDiff[0][Coll], MarkerType[0], titleX, titleYToOOJ, "" , 1,0, 40, xOffset, yOffset, MarkerSize[0]);
	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetTitleOffset(0.7);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetTitleOffset(0.7);

	  if (Coll==0)	  histoYieldRatioToOOJ[Coll]->Fit(fitToRatioToOOJ, "R+");
	  if (Coll==0)	{
	    legendFitRatio1 = legendFitRatioToOOJ->AddEntry(fitToRatioToOOJ, "pol0 fit", "l");
	    legendFitRatio1->SetTextAlign(22);
	    legendFitRatio1->SetTextSize(0.08);
	  }
	  histoYieldRatioToOOJ[Coll]->DrawClone("same e");
	  histoYieldRatioSistToOOJ[Coll]->SetFillStyle(0);
	  histoYieldRatioSistToOOJ[Coll]->DrawClone("same e2");

	  canvasSplit->cd();
	  gStyle->SetLegendBorderSize(0);
	  gStyle->SetLegendFillColor(0);
	  gStyle->SetLegendFont(42);
	  pad2->Draw();
	  pad2->cd();

	  if (Coll==0){
	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetTitleOffset(0.7);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetTitleOffset(0.7);

	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetLabelSize(0.10);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetLabelSize(0.10);
	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetLabelOffset(0.009);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetLabelOffset(0.009);

	  histoYieldRatioToOOJ[Coll]->GetXaxis()->SetTitleSize(0.10);
	  histoYieldRatioSistToOOJ[Coll]->GetXaxis()->SetTitleSize(0.10);
	  histoYieldRatioToOOJ[Coll]->GetXaxis()->SetTitleOffset(1.3);
	  histoYieldRatioSistToOOJ[Coll]->GetXaxis()->SetTitleOffset(1.3);

	  histoYieldRatioToOOJ[Coll]->GetXaxis()->SetLabelSize(0.10);
	  histoYieldRatioSistToOOJ[Coll]->GetXaxis()->SetLabelSize(0.10);
	  histoYieldRatioToOOJ[Coll]->GetXaxis()->SetLabelOffset(0.02);
	  histoYieldRatioSistToOOJ[Coll]->GetXaxis()->SetLabelOffset(0.02);

	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetNdivisions(305);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetNdivisions(305);
	  histoYieldRatioToOOJ[Coll]->GetYaxis()->SetTickLength(0.02);
	  histoYieldRatioSistToOOJ[Coll]->GetYaxis()->SetTickLength(0.02);
	  }

	  hdummyRatioToOOJ->Draw("same");
	  gYieldRatioToOOJ[Coll]->Draw("same p e");
	  gYieldRatioToOOJSist[Coll]->SetFillStyle(0);
	  gYieldRatioToOOJSist[Coll]->SetFillColor(ColorDiff[0][Coll]);
	  gYieldRatioToOOJSist[Coll]->Draw("same p2");
	  if (isMultCorrEval){
	    gYieldRatioToOOJSistMultUnCorr[Coll]->SetFillStyle(3001);
	    //	    gYieldRatioToOOJSistMultUnCorr[Coll]->SetFillColorAlpha(ColorDiff[0][Coll], 0.03);	 
	    gYieldRatioToOOJSistMultUnCorr[Coll]->SetFillColor(ColorDiff[0][Coll]);
	    gYieldRatioToOOJSistMultUnCorr[Coll]->Draw("same p2");
	  }
	  fitToRatioToOOJ->Draw("same");
	  
	  //histoYieldRatioToOOJ[Coll]->DrawClone("same e0x0");
	  histoYieldRatioSistToOOJ[Coll]->SetFillStyle(0);
	  //histoYieldRatioSistToOOJ[Coll]->DrawClone("same e2");
	  
	  legendFitRatioToOOJ->Draw("");
	}

	canvas->cd();
	canvas->SetFillColor(0);
	canvas->SetTickx(1);
	canvas->SetTicky(1);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendFont(42);
	StyleHisto(histoYieldRatio[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	StyleHisto(histoYieldRatioSist[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	if (isMultCorrEval) StyleHisto(histoYieldRatioSistMultUnCorr[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);

	if (PlotType==0) {
	  histoYieldRatio[ireg][Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatioSist[ireg][Coll]->GetYaxis()->SetTitleSize(0.07);
	  if (isMultCorrEval)  histoYieldRatioSistMultUnCorr[ireg][Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatio[ireg][Coll]->GetYaxis()->SetTitleOffset(0.7);
	  histoYieldRatioSist[ireg][Coll]->GetYaxis()->SetTitleOffset(0.7);
	  if (isMultCorrEval)  histoYieldRatioSistMultUnCorr[ireg][Coll]->GetYaxis()->SetTitleOffset(0.7);
	  histoYieldRatio[ireg][Coll]->GetYaxis()->SetTickLength(0.02);
	  histoYieldRatioSist[ireg][Coll]->GetYaxis()->SetTickLength(0.02);
	  if (isMultCorrEval)  histoYieldRatioSistMultUnCorr[ireg][Coll]->GetYaxis()->SetTickLength(0.02);
	}
	//legend
	if (NLoopType ==0 && NLoopRegion==0){
	  fHistYieldStatGrey[Coll]= (TH1F*)      histoYieldRatio[ireg][Coll]->Clone("fHistYieldStatBlack");
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
	if (NLoopType==0 && Coll==0){
	  if (ireg==0) lReAll1Bis[ireg] = legendRegionJet->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  else if (ireg==1)lReAll1Bis[ireg] = legendRegionBulk->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  else lReAll1Bis[ireg] = legendRegionFull->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  lReAll1Bis[ireg]->SetTextSize(0.045);
	  //      lReAll1Bis[ireg]->SetTextAlign(32);

	  if (ireg==0) lReAll2Bis[ireg]=      legendRegionJet->AddEntry("", sRegion1[ireg], "");
	  else    if (ireg==1) lReAll2Bis[ireg]=      legendRegionBulk->AddEntry("", sRegion1[ireg], "");
	  else    lReAll2Bis[ireg]=      legendRegionFull->AddEntry("", sRegion1[ireg], "");
	  lReAll2Bis[ireg]->SetTextSize(0.035);

	  lReAll1[ireg] = legendRegionAll->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  lReAll1[ireg]->SetTextSize(0.048);
	  lReAll2[ireg]= legendRegionAll->AddEntry("", sRegion1[ireg], "");
	  lReAll2[ireg]->SetTextSize(0.038);

	  if (NLoopRegion==0){
	    fHistYieldStatBlack= (TH1F*)      histoYieldRatio[ireg][Coll]->Clone("fHistYieldStatBlack");
	    fHistYieldSistBlack= (TH1F*)      histoYieldRatioSist[ireg][Coll]->Clone("fHistYieldSistBlack");
	    if (isMultCorrEval) fHistYieldSistMultUnCorrBlack= (TH1F*)      histoYieldRatioSistMultUnCorr[ireg][Coll]->Clone("fHistYieldSistMultUnCorrBlack");
	    fHistYieldStatBlack->SetLineColor(1);
	    fHistYieldStatBlack->SetMarkerColor(1);
	    fHistYieldStatBlack->SetMarkerStyle(20);
	    fHistYieldSistBlack->SetLineColor(1);
	    fHistYieldSistBlack->SetMarkerColor(1);
	    fHistYieldSistBlack->SetMarkerStyle(20);
	    if (isMultCorrEval) {
	    fHistYieldSistMultUnCorrBlack->SetLineColor(1);
	    fHistYieldSistMultUnCorrBlack->SetMarkerColor(1);
	    fHistYieldSistMultUnCorrBlack->SetMarkerStyle(20);
	    fHistYieldSistMultUnCorrBlack->SetFillStyle(3001);
	    //	    fHistYieldSistMultUnCorrBlack->SetFillColorAlpha(1, 0.03);
	    fHistYieldSistMultUnCorrBlack->SetFillColor(1);
	    }

	    legendStatBoxBis->AddEntry(fHistYieldStatBlack, "stat.", "pe");
	    legendStatBoxBis->AddEntry(fHistYieldSistBlack, "syst.", "ef");
	    if (isMultCorrEval)  legendStatBoxBis->AddEntry(fHistYieldSistMultUnCorrBlack, "syst. uncorr.", "ef");
	    legendStatBox->AddEntry(fHistYieldStatBlack, "stat.", "pe");
	    legendStatBox->AddEntry(fHistYieldSistBlack, "syst.", "ef");
	    if (isMultCorrEval)  legendStatBox->AddEntry(fHistYieldSistMultUnCorrBlack, "syst. uncorr.", "ef");

	  }
	}
	if (Coll==1){
	  legendStatBoxColor->AddEntry(histoYieldRatio[ireg][Coll], "stat.", "pe");
	  legendStatBoxColor->AddEntry(histoYieldRatioSist[ireg][Coll], "syst.", "ef");
	  if (isMultCorrEval) {
	    histoYieldRatioSistMultUnCorr[ireg][Coll]->SetFillStyle(3001);
	    //	    histoYieldRatioSistMultUnCorr[ireg][Coll]->SetFillColorAlpha(ColorDiff[ireg][Coll], 0.03);
	    histoYieldRatioSistMultUnCorr[ireg][Coll]->SetFillColor(ColorDiff[ireg][Coll]);
	    legendStatBoxColor->AddEntry(histoYieldRatioSistMultUnCorr[ireg][Coll], "syst. uncorr.", "ef");
	  }

	  TLegendEntry * lRe1 = legendOneRegion->AddEntry("", sRegion[ireg], "");
	  lRe1->SetTextSize(0.06);
	  lRe1->SetTextAlign(32);
	  TLegendEntry *leR2=      legendOneRegion->AddEntry("", sRegion1[ireg], "");
	  leR2->SetTextSize(0.04);
	  leR2->SetTextAlign(32);
	}
	if (Coll==0)	LegendColor->AddEntry(""/*(TObject*)0*/, NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");

	if (Coll==1)	legendEnergyBoxColor->AddEntry(histoYieldRatio[ireg][Coll], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
	else 	legendEnergyBoxColor->AddEntry(histoYieldRatio[ireg][Coll], "pp, #sqrt{#it{s}} = 13 TeV", "pe");

	histoYieldDummy->Draw("same");
	if (!(type==1 && PlotType == 0)){
	  gYieldRatio[ireg][Coll]->Draw("same pe");
	  gYieldRatioSist[ireg][Coll]->SetFillStyle(0);
	  gYieldRatioSist[ireg][Coll]->SetFillColor(ColorDiff[ireg][Coll]);
	  gYieldRatioSist[ireg][Coll]->Draw("same p2");
	  if (isMultCorrEval){
	    gYieldRatioSistMultUnCorr[ireg][Coll]->SetFillStyle(3001);
	    //	    gYieldRatioSistMultUnCorr[ireg][Coll]->SetFillColorAlpha(ColorDiff[ireg][Coll], 0.03);
	    gYieldRatioSistMultUnCorr[ireg][Coll]->SetFillColor(ColorDiff[ireg][Coll]);
	    gYieldRatioSistMultUnCorr[ireg][Coll]->Draw("same p2");
	  }
	  for (Int_t i =0; i<gYieldRatioSist[ireg][Coll]->GetN(); i++){
	    double x;
	    double y;
	    gYieldRatioSist[ireg][Coll]->GetPoint(i, x, y);
	    cout << " dNdeta " << x << " value " << y << endl;
	  }
	}
	//histoYieldRatio[ireg][Coll]->Draw("same e0x0");
	histoYieldRatioSist[ireg][Coll]->SetFillStyle(0);
	//	if (isMultCorrEval) histoYieldRatioSistMultUnCorr[ireg][Coll]->SetFillStyle(1001);
	//histoYieldRatioSist[ireg][Coll]->Draw("same e2");
	if (isFit){
	  fitPol0[ireg]->Draw("same");
	  fitPol1[ireg]->Draw("same");
	}

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

	canvasSplit->cd();
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendFont(42);
	pad1->Draw();
	pad1->cd();
	if (PlotType==0 && type==0){
	  histoYieldDummy->Draw("same");
	  //	  gYieldRatio[ireg][Coll]->Draw("same p e0x0");
	  gYieldRatio[ireg][Coll]->Draw("same pe");
	  gYieldRatioSist[ireg][Coll]->SetFillStyle(0);
	  gYieldRatioSist[ireg][Coll]->SetFillColor(ColorDiff[ireg][Coll]);
	  gYieldRatioSist[ireg][Coll]->Draw("same p2");
	  if (isMultCorrEval){
	    gYieldRatioSistMultUnCorr[ireg][Coll]->SetFillStyle(3001);
	    gYieldRatioSistMultUnCorr[ireg][Coll]->SetFillColorAlpha(ColorDiff[ireg][Coll], 0.03);
	    gYieldRatioSistMultUnCorr[ireg][Coll]->SetFillColor(ColorDiff[ireg][Coll]);
	    gYieldRatioSistMultUnCorr[ireg][Coll]->Draw("same p2");
	  }
	  //histoYieldRatio[ireg][Coll]->Draw("same e0x0");
	  histoYieldRatioSist[ireg][Coll]->SetFillStyle(0);
	  //histoYieldRatioSist[ireg][Coll]->Draw("same e2");
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
      }
    }
  }

  if (isFit){
    cout << "\nPol0/pol1 fit results: " << endl;
    for (Int_t ireg=0; ireg<numRegions; ireg++){
      if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
      cout << Region[ireg] << ":\npol1: " << fitPol1[ireg]->GetParameter(1) << " +- " << fitPol1[ireg]->GetParError(1) << " Chi2/NDF " << fitPol1[ireg]->GetChisquare()<<"/"<<fitPol1[ireg]->GetNDF() << " = " << fitPol1[ireg]->GetChisquare()/fitPol1[ireg]->GetNDF()<< "\npol0: " << fitPol0[ireg]->GetParameter(0) << " +- " << fitPol0[ireg]->GetParError(0) << " Chi2/NDF " <<  fitPol0[ireg]->GetChisquare()<< "/" << fitPol0[ireg]->GetNDF()<< " = " <<  fitPol0[ireg]->GetChisquare()/fitPol0[ireg]->GetNDF()<< endl;
    }
  }
  if (ChosenRegion < 0 && PlotType==0){
    cout << "\n\nXi/K0s (JET) / Xi/0s (OOJ) " << endl;
    cout << "Chi2/NDF " << fitToRatioToOOJ->GetChisquare()<< "/" <<fitToRatioToOOJ->GetNDF()<< " =" << fitToRatioToOOJ->GetChisquare()/fitToRatioToOOJ->GetNDF()<< endl;
  }
  if (ChosenRegion<0){
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+SisFit[isFit] +".pdf");
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+SisFit[isFit] +".eps");
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+SisFit[isFit] +".png");
    //canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+SisFit[isFit] +".jpeg"); //non bene
    if (PlotType==0) {
      canvasRatioToOOJ->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+"_RatioToOOJ.pdf");
      canvasSplit->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+"_WithRatioToOOJ.pdf");
      canvasRatioToOOJ->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+"_RatioToOOJ.eps");
      canvasSplit->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary]+"_WithRatioToOOJ.eps");
    }
  }
  else {
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+Region[ChosenRegion]+SisPreliminary[isPreliminary]+".pdf");
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+Region[ChosenRegion]+SisPreliminary[isPreliminary]+".eps");
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+Region[ChosenRegion]+SisPreliminary[isPreliminary]+".png");
    //    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+Region[ChosenRegion]+SisPreliminary[isPreliminary]+".jpeg"); //non bene
  }

  TString fileoutHistosName = "RatiosXiK0s3Systems_" + SPlotType[PlotType];
  if (ChosenRegion>=0) fileoutHistosName += "_" + RegionType[ChosenRegion];
  if (isdNdEtaTriggered) fileoutHistosName +="_isdNdEtaTriggered";
  if (isWingsCorrectionApplied) fileoutHistosName += "_WingsCorrApplied";
  if (MaterialBudgetCorr) fileoutHistosName += "_MatBudgetCorrFAST";
  if (isMultCorrEval)   fileoutHistosName += "_MultCorrSistEval";
  fileoutHistosName += ".root";
  TFile * fileoutHistos = new TFile(fileoutHistosName, "RECREATE");
  canvas->Write();
  canvasSplit->Write();
  canvasRatioToOOJ->Write();
  cout <<" Saving histos in file " << endl;
  for (Int_t ireg=0; ireg<numRegions; ireg++){
    if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
    NLoopRegion++;
    for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
      if (Coll==0){
	histoYieldRatio[ireg][Coll]->SetName(SPlotType[PlotType] + "_" + RegionType[ireg] + "_13TeV_Stat");
	histoYieldRatioSist[ireg][Coll]->SetName(SPlotType[PlotType] + "_" + RegionType[ireg] + "_13TeV_Sist");
	if (isMultCorrEval) histoYieldRatioSistMultUnCorr[ireg][Coll]->SetName(SPlotType[PlotType] + "_" + RegionType[ireg] + "_13TeV_SistMultUnCorr");
      }
      else {
	histoYieldRatio[ireg][Coll]->SetName(SPlotType[PlotType] + "_" + RegionType[ireg] + "_5TeV_Stat");
	histoYieldRatioSist[ireg][Coll]->SetName(SPlotType[PlotType] + "_" + RegionType[ireg] + "_5TeV_Sist");
	if (isMultCorrEval) histoYieldRatioSistMultUnCorr[ireg][Coll]->SetName(SPlotType[PlotType] + "_" + RegionType[ireg] + "_5TeV_SistMultUnCorr");
      }
      histoYieldRatio[ireg][Coll]->Write();
      histoYieldRatioSist[ireg][Coll]->Write();
      if (isMultCorrEval) histoYieldRatioSistMultUnCorr[ireg][Coll]->Write();
    }
  }
  cout << "I have created the file " << fileoutHistosName << endl;
  fileoutHistos->Close();

  TFile* fileout = new TFile("FinalPlot.root", "RECREATE");
  canvas->Write();
  fileout->Close();
}
