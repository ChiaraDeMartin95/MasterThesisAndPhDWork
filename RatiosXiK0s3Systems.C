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

TString titlePtYield=  "1/#Delta#it{#eta} #Delta#it{#varphi} 1/#it{N}_{trigg} d#it{N}/d#it{p}_{T} [(GeV/#it{c})^{-1}]";

TString titlePt=  "#it{p}_{T} (GeV/#it{c})";
TString titleYieldY="#it{N}/#it{N}_{trigg} 1/#Delta#eta #Delta#it{#varphi}";
TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}"};
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 1.8", "|#De\
lta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};


void RatiosXiK0s3Systems( Int_t PlotType =0, Int_t ChosenRegion = -1,  Float_t ScalingFactorXiK0s = 1/*0.8458/1.08747*/, Bool_t isPreliminary =0){

  //PlotType = 0: Xi/K0s ratio
  //PlotType = 1: K0s yield vs mult 
  //PlotType = 2: Xi yield vs mult 
  //PlotType = 3: K0s pt vs mult 
  //PlotType = 4: Xi pt vs mult 

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

  //ChosenRegion == -1 //all regions (JET, OOJ, FULL) are plotted
  if (ChosenRegion>2) {cout << "Choose a valid region! " << endl; return;}
  const Int_t numColls =2; //pp13 TeV (MB + HM), pp5TeV
  const Int_t numRegions =3;
  const Int_t numTypes =2;
  TString CollisionsComp[2] = {"_HMMultBinning1_vsHM", "_vs5TeV"};
  TString SPlotType[5] = {"Ratio", "K0s", "Xi", "K0spt", "Xipt"};
  TString SisPreliminary[2] = {"", "_isPreliminary"};

  gStyle->SetOptStat(0);
  TString Region[numRegions] = {"Jet", "Bulk", "All"};
  TString RegionBis[numRegions] = {"Jet", "Bulk", "Inclusive"};
  Float_t Up= 0.14-10e-6;
  Float_t Low = 0.02+10e-6;
  if (PlotType==1) {
    Low = 10e-8; Up = 0.45-10e-6;
    if (ChosenRegion==0) {Low = 0.015; Up = 0.05;} //0.035 
  }
  if (PlotType==2) {
    Low = 10e-8; Up = 0.038-10e-6;
    if (ChosenRegion==0) {Up = 0.003;}
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
  Int_t ColorDiff[numRegions][numColls] = {{634,628}, {418,829} , {601, 867}};
  Int_t Marker[numTypes] = {33, 21};
  Float_t Size[numTypes] = {2, 1.4};
  TString ParticleType[numTypes]={"K0s", "Xi"};
  TH1F *histoYield[numTypes][numRegions][numColls];
  TH1F *histoYieldRatio[numRegions][numColls];
  TH1F *histoYieldSist[numTypes][numRegions][numColls];
  TH1F *histoYieldRatioSist[numRegions][numColls];
  TH1F*  fHistYieldStatBlack;
  TH1F*  fHistYieldSistBlack;
  TString NameHisto13TeV="histoYieldComparison";
  TString NameHistoSist13TeV="histoYieldSistComparison";
  TString NameHisto5TeV="fHistYield_pp5TeV";
  TString NameHistoSist5TeV="fHistYield_pp5TeV_Sist";
  TString NameHistoFinal[numTypes][numRegions][numColls];
  TString pathin ="";
  TString YieldOrAvgPt[5] = {"Yield", "Yield", "Yield", "AvgPt","AvgPt"};
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  //legend
  //      TLegend *Legend1B=new TLegend(0.1,0.72,0.38,0.92);
  TLegend *Legend1B=new TLegend(0.62,0.72,0.9,0.92);
  Legend1B->SetFillStyle(0);
  TLegendEntry* E1Bis =      Legend1B->AddEntry("", "#bf{ALICE Preliminary}", "");
  //TLegendEntry* E1Bis =      Legend1B->AddEntry("", "", "");
  E1Bis->SetTextAlign(32);
  TLegendEntry* E2Bis;
  //  if (ispp5TeV)   E2Bis =        Legend1B->AddEntry("", "pp, #sqrt{#it{s}} = 5 TeV", "");
  E2Bis =       Legend1B->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  E2Bis->SetTextAlign(32);
  TLegendEntry* E3Bis =           Legend1B->AddEntry(""/*(TObject*)0*/, "#it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  E3Bis->SetTextAlign(32);

  TLegend *Legend2=new TLegend(0.16,0.75,0.5,0.93);
  //      Legend2->SetFillStyle(0);
  Legend2->SetMargin(0);
  //      Legend2->AddEntry("", "#bf{ALICE Preliminary}", "");
  Legend2->AddEntry("", "", "");
  Legend2->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");

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

  TLegend *legendRegionAllF=new TLegend(0.15, 0.43, 0.59, 0.71);
  legendRegionAllF->SetFillStyle(0);
  legendRegionAllF->SetMargin(0.07);
  TLegendEntry * lReAll1[3];
  TLegendEntry *lReAll2[3];

  TLegend *legendStatBox=new TLegend(0.7, 0.76, 0.9, 0.88);
  TLegend *legendStatBoxBis=new TLegend(0.2, 0.48, 0.35, 0.58);

  Int_t NLoopType =-1;
  Int_t NLoopRegion =-1;
  for (Int_t type=1; type>=0; type--){
    cout << "Type: " << type << endl;
    if ((PlotType == 1 || PlotType ==3) && type==1) continue;
    else  if ((PlotType == 2 || PlotType ==4) && type==0) continue;
    NLoopType++;
    Legend2->AddEntry(""/*(TObject*)0*/, NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
    for (Int_t ireg=0; ireg<numRegions; ireg++){
      if (ChosenRegion>=0 && ireg != ChosenRegion) continue;
      NLoopRegion++;
      cout << "\n\n*****" << ParticleType[type] << " " << Region[ireg]<< endl;
      for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
	pathin = "Compare" +YieldOrAvgPt[PlotType] + "DifferentCollisions";
	pathin +=  CollisionsComp[Coll];
	pathin +="_"+ ParticleType[type] + Region[ireg];
	if (isPreliminary) pathin +=  "_isPreliminary";
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

	if (Coll==0)	histoYieldSist[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoSist13TeV);
	else histoYieldSist[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoSist5TeV);
	if (!histoYieldSist[type][ireg][Coll]) {cout <<"no histo " << endl; return;}
	histoYieldSist[type][ireg][Coll]->SetName(NameHistoFinal[type][ireg][Coll]+"Sist");

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
	  histoYieldRatio[ireg][Coll]->Sumw2();
	  histoYieldRatioSist[ireg][Coll]->Sumw2();
	}
	else {
	  if (PlotType==0){
	    histoYieldRatio[ireg][Coll]->Divide(histoYield[type][ireg][Coll]);
	    histoYieldRatioSist[ireg][Coll]->Divide(histoYieldSist[type][ireg][Coll]);
	    if (ireg ==0) {
	      histoYieldRatio[ireg][Coll]->Scale(ScalingFactorXiK0s);
	      histoYieldRatioSist[ireg][Coll]->Scale(ScalingFactorXiK0s);
	    }
	  }
	}
	if (type==0){
	  cout << "\nXi/K0s ratio: " << endl;
	  for (Int_t b=1; b<=  histoYieldRatio[ireg][Coll]->GetNbinsX(); b++){
	    if (histoYieldRatio[ireg][Coll]->GetBinContent(b) != 0)	 {
	      cout << histoYieldRatio[ireg][Coll]->GetBinContent(b) << " +- " << histoYieldRatio[ireg][Coll]->GetBinError(b) << " (stat.) +-  " << histoYieldRatioSist[ireg][Coll]->GetBinError(b) << " (syst.) " << endl;
	    }
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

	StyleHisto(histoYieldRatio[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	StyleHisto(histoYieldRatioSist[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], MarkerType[ireg], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, MarkerSize[ireg]);
	if (PlotType==0) {
	  histoYieldRatio[ireg][Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatioSist[ireg][Coll]->GetYaxis()->SetTitleSize(0.07);
	  histoYieldRatio[ireg][Coll]->GetYaxis()->SetTitleOffset(0.7);
	  histoYieldRatioSist[ireg][Coll]->GetYaxis()->SetTitleOffset(0.7);
	}

	//legend
	if (NLoopType==0 && Coll==0){
	  if (ireg==0) lReAll1Bis[ireg] = legendRegionAllFJet->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  else if (ireg==1)lReAll1Bis[ireg] = legendRegionAllFBulk->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  else lReAll1Bis[ireg] = legendRegionAllFAll->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  lReAll1Bis[ireg]->SetTextSize(0.045);
	  //      lReAll1Bis[ireg]->SetTextAlign(32);

	  if (ireg==0) lReAll2Bis[ireg]=      legendRegionAllFJet->AddEntry("", sRegion1[ireg], "");
	  else    if (ireg==1) lReAll2Bis[ireg]=      legendRegionAllFBulk->AddEntry("", sRegion1[ireg], "");
	  else    lReAll2Bis[ireg]=      legendRegionAllFAll->AddEntry("", sRegion1[ireg], "");
	  lReAll2Bis[ireg]->SetTextSize(0.035);

	  lReAll1[ireg] = legendRegionAllF->AddEntry(histoYieldRatio[ireg][Coll], sRegion[ireg], "p");
	  lReAll1[ireg]->SetTextSize(0.048);
	  lReAll2[ireg]= legendRegionAllF->AddEntry("", sRegion1[ireg], "");
	  lReAll2[ireg]->SetTextSize(0.038);

	  if (NLoopRegion==0){
	    fHistYieldStatBlack= (TH1F*)      histoYieldRatio[ireg][Coll]->Clone("fHistYieldStatBlack");
	    fHistYieldSistBlack= (TH1F*)      histoYieldRatioSist[ireg][Coll]->Clone("fHistYieldSistBlack");
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

	histoYieldRatio[ireg][Coll]->Draw("same e0x0");
	histoYieldRatioSist[ireg][Coll]->SetFillStyle(0);
	histoYieldRatioSist[ireg][Coll]->Draw("same e2");

	if (PlotType==0){
	Legend1B->Draw("");
        legendRegionAllFJet->Draw("");
        legendRegionAllFBulk->Draw("");
        legendRegionAllFAll->Draw("");
        legendStatBoxBis->Draw("");
	}
	else if (PlotType==1 || PlotType==2){
	  Legend2->Draw("");
	  legendStatBox->Draw("");
	  legendRegionAllF->Draw("");
	}
      }
    }
  }

  if (ChosenRegion<0){
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+SisPreliminary[isPreliminary] +".pdf");
  }
  else     canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+Region[ChosenRegion]+SisPreliminary[isPreliminary]+".pdf");
}
