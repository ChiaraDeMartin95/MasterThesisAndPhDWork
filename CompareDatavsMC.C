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
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>

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

TString TitleYieldRatio="#Xi/K^{0}_{S} yield ratio vs multiplicity";
TString titleRelUnc = "Relative uncertainty";
TString titleMult = "Multiplicity class ";
TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
TString titleYield[2]={"K^{0}_{S}", "#Xi"};
TString TitleYYieldRatio="#it{N}_{#Xi}/#it{N}_{K^{0}_{S}}";

//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}"};
TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})", "#it{N}_{#Xi}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})"};
//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)"}; //to much space around the parenthesis
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString RegionTypeBis[3] = {"Jet", "OOj", "Full"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
//TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
//TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};
//TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 0.85", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 2.0", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};


void CompareDatavsMC( Int_t PlotType =0, Int_t ChosenRegion = -1,  Float_t ScalingFactorXiK0s = 1/*0.8458/1.08747*/, Int_t isPreliminary =0, Bool_t isChangesIncluded=1, Bool_t isFit=0){

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

  //Set titles
  TString titleY;
  if (PlotType==0) titleY = TitleYYieldRatio;
  else if (PlotType==1) titleY = titleYieldYType[0];
  else if (PlotType==2) titleY = titleYieldYType[1];
  else if (PlotType==3) titleY = titlePtvsMultYType[0];
  else if (PlotType==4) titleY = titlePtvsMultYType[1];

  TString  titleX = titledNdeta;
  Float_t xOffset =1.2;
  Float_t yOffset =1.25;
  Float_t MarkerSize[3] ={2, 2, 3};
  Int_t MarkerType[3] = {20, 21, 33};

  //  TString titleYToOOJ = "#Xi/K^{0}_{S} (JET) / #Xi/K^{0}_{S} (OOJ)";
  //  TString titleYToOOJ = "Near-side jet / out-of-jet";
  TString titleYToOOJ = "Toward / Transverse";

  //ChosenRegion == -1 //all regions (JET, OOJ, FULL) are plotted
  if (ChosenRegion>2) {cout << "Choose a valid region! " << endl; return;}
  const Int_t numColls =2; //pp13 TeV (MB + HM), pp5TeV
  const Int_t numRegions =3;
  const Int_t numDataorMC =2;
  TString CollisionsComp[2] = {"_HMMultBinning1_vsHM", "_vs5TeV"};
  TString SPlotType[5] = {"Ratio", "K0s", "Xi", "K0spt", "Xipt"};
  TString SPlotTypeBis[5] = {"", "K0s", "Xi", "K0s", "Xi"};
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
  Int_t ColorDiffMC[numRegions][numColls] = {{907, 907}, {812,812} , {870, 870}};
  Int_t Marker[2] = {33, 21};
  Float_t Size[2] = {2, 1.4};

  TF1* fitPol1[3];
  TF1* fitPol0[3];

  TF1 * fitToRatioToOOJ = new TF1("fitToRatioToOOJ", "pol0", 0, 45);
  fitToRatioToOOJ->SetLineColor(1);
  fitToRatioToOOJ->SetLineStyle(2);
  fitToRatioToOOJ->SetLineWidth(2);

  TH1F *histoYield[numRegions][numColls][numDataorMC];
  TH1F *histoYieldSist[numRegions][numColls][numDataorMC];
  TH1F*  fHistYieldStatBlack;
  TH1F*  fHistYieldSistBlack;
  TH1F*  fHistYieldStatGrey[numColls];
  TString NameHisto13TeV="histoYieldComparison";
  TString NameHistoSist13TeV="histoYieldSistComparison";
  TString NameHisto5TeV="fHistYield_pp5TeV";
  TString NameHistoSist5TeV="fHistYield_pp5TeV_Sist";
  TString NameHistoFinal[numRegions][numColls][numDataorMC];
  TString pathin[2] ={""};
  TString YieldOrAvgPt[5] = {"Yield", "Yield", "Yield", "AvgPt","AvgPt"};
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  TPad* pad1 = new TPad( "pad1" ,"pad1" ,0 ,0.36 ,1 ,1);
  TPad* pad2 = new TPad( "pad2" ,"pad2" ,0 ,0.01 ,1 ,0.35);

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
  TLegendEntry* E1Ratio =      LegendRatio->AddEntry("", "#bf{ALICE Preliminary}", "");
  TLegendEntry* E3Ratio =           LegendRatio->AddEntry(""/*(TObject*)0*/, "#it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  E3Ratio->SetTextAlign(32);

  TLegend *LegendColor=new TLegend(0.16,0.80,0.5,0.92); 
  LegendColor->SetMargin(0);
  LegendColor->AddEntry("", "#bf{ALICE Preliminary}", "");
  //LegendColor->AddEntry("", "", "");

  TLegend *LegendYields=new TLegend(0.16,0.75,0.5,0.93);
  LegendYields->SetMargin(0);
  LegendYields->AddEntry("", "#bf{ALICE Preliminary}", "");
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

  TLegend *legendEnergyBoxBis1=new TLegend(0.68, 0.69, 0.89, 0.76);
  TLegend *legendEnergyBoxBis2=new TLegend(0.68, 0.62, 0.89, 0.69);

  TLegend *legendStatBox=new TLegend(0.73, 0.80, 0.93, 0.92);
  TLegend *legendStatBoxBis=new TLegend(0.2, 0.48, 0.35, 0.58);
  TLegend *legendStatBoxColor=new TLegend(0.18, 0.56, 0.35, 0.71);

  Int_t NLoopType =-1;
  Int_t NLoopRegion =-1;
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
      for (Int_t isMC = 0; isMC < 2; isMC++){
	cout << "\nisMC? " << isMC << endl;
	//Input file 
	if (isMC==0) {
	  pathin[isMC] = "RatiosXiK0s3Systems_" + SPlotType[PlotType];
	  if (ChosenRegion>=0) pathin[isMC] += "_" + RegionType[ChosenRegion];
	  pathin[isMC] += ".root";
	}
	else {
	  pathin[isMC] = "Results_pp13TeVMB_FastMCPrediction.root"; //temporary solution
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
	    else    NameHisto13TeV = "PtvsMult_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_StatErr";
	    if (PlotType<=2)   NameHistoSist13TeV = "Yield_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_SistErr";
	    else    NameHistoSist13TeV = "PtvsMult_"+ SPlotTypeBis[PlotType] + "_" + RegionTypeBis[ireg] + "_SistErr";	
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

	if (isMC==0){
	  StyleHisto(histoYield[ireg][Coll][isMC], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	  StyleHisto(histoYieldSist[ireg][Coll][isMC], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	}
	else {
	  StyleHisto(histoYield[ireg][Coll][isMC], Low, Up, ColorDiffMC[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	  StyleHisto(histoYieldSist[ireg][Coll][isMC], Low, Up, ColorDiffMC[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	}
	if (PlotType==0) {
	  histoYield[ireg][Coll][isMC]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldSist[ireg][Coll][isMC]->GetYaxis()->SetTitleSize(0.07);
	  histoYield[ireg][Coll][isMC]->GetYaxis()->SetTitleOffset(0.7);
	  histoYieldSist[ireg][Coll][isMC]->GetYaxis()->SetTitleOffset(0.7);
	  histoYield[ireg][Coll][isMC]->GetYaxis()->SetTickLength(0.02);
	  histoYieldSist[ireg][Coll][isMC]->GetYaxis()->SetTickLength(0.02);
	}

	/*
	//legend
	if (NLoopType ==0 && NLoopRegion==0){
	  fHistYieldStatGrey[Coll]= (TH1F*)      histoYield[ireg][Coll][isMC]->Clone("fHistYieldStatBlack");
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

	if (Coll==1){
	  legendStatBoxColor->AddEntry(histoYield[ireg][Coll][isMC], "stat. error", "pe");
	  legendStatBoxColor->AddEntry(histoYieldSist[ireg][Coll][isMC], "syst. error", "ef");

	  TLegendEntry * lRe1 = legendOneRegion->AddEntry("", sRegion[ireg], "");
	  lRe1->SetTextSize(0.06);
	  lRe1->SetTextAlign(32);
	  TLegendEntry *leR2=      legendOneRegion->AddEntry("", sRegion1[ireg], "");
	  leR2->SetTextSize(0.04);
	  leR2->SetTextAlign(32);
	}
	if (Coll==0)	LegendColor->AddEntry("", NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");

	if (Coll==1)	legendEnergyBoxColor->AddEntry(histoYield[ireg][Coll][isMC], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
	else 	legendEnergyBoxColor->AddEntry(histoYield[ireg][Coll][isMC], "pp, #sqrt{#it{s}} = 13 TeV", "pe");

	*/
	histoYield[ireg][Coll][isMC]->Draw("same e0x0");
	histoYieldSist[ireg][Coll][isMC]->SetFillStyle(0);
	histoYieldSist[ireg][Coll][isMC]->Draw("same e2");

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
      }
    }
  }

  if (ChosenRegion<0){
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType] +".pdf");
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType] +".eps");
  }
  else {
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType]+"_" + RegionType[ChosenRegion]+".pdf");
    canvas->SaveAs("CompareDatavsMC_"+SPlotType[PlotType]+"_" + RegionType[ChosenRegion]+".eps");
  }
}
