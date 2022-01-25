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
#include <TFile.h>
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(2);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.8);
  histo->SetTitle(title);
}

void RatiosXiK0s3Systems( Int_t PlotType =0, Int_t ChosenRegion = -1,  Float_t ScalingFactorXiK0s = 1/*0.8458/1.08747*/){

  //PlotType = 0: Xi/K0s ratio
  //PlotType = 1: K0s yield vs mult 
  //PlotType = 2: Xi yield vs mult 
  //PlotType = 3: K0s pt vs mult 
  //PlotType = 4: Xi pt vs mult 

  //ChosenRegion == -1 //all regions (JET, OOJ, FULL) are plotted
  if (ChosenRegion>2) {cout << "Choose a valid region! " << endl; return;}
  const Int_t numColls =2; //pp13 TeV (MB + HM), pp5TeV
  const Int_t numRegions =3;
  const Int_t numTypes =2;
  TString CollisionsComp[2] = {"_HMMultBinning1_vsHM", "_vs5TeV"};
  TString SPlotType[5] = {"Ratio", "K0s", "Xi", "K0spt", "Xipt"};

  gStyle->SetOptStat(0);
  TString Region[numRegions] = {"Jet", "Bulk", "All"};
  TString RegionBis[numRegions] = {"Jet", "Bulk", "Inclusive"};
  Float_t Up= 0.1;
  Float_t Low = 10e-4;
  if (PlotType==1) {
    Low = 10e-8; Up = 0.3;
    if (ChosenRegion==0) {Low = 0.015; Up = 0.035;}
  }
  if (PlotType==2) {
    Low = 10e-8; Up = 0.025;
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
  TString NameHisto13TeV="histoYieldComparison";
  TString NameHistoSist13TeV="histoYieldSistComparison";
  TString NameHisto5TeV="fHistYield_pp5TeV";
  TString NameHistoSist5TeV="fHistYield_pp5TeV_Sist";
  TString NameHistoFinal[numTypes][numRegions][numColls];
  TString pathin ="";
  TString YieldOrAvgPt[5] = {"Yield", "Yield", "Yield", "AvgPt","AvgPt"};
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);

  for (Int_t type=1; type>=0; type--){
    if ((PlotType == 1 || PlotType ==3) && type==1) continue;
    else  if ((PlotType == 2 || PlotType ==4) && type==0) continue;
    for (Int_t ireg=0; ireg<numRegions; ireg++){
      if (ireg != ChosenRegion) continue;
      cout << "\n\n*****" << ParticleType[type] << " " << Region[ireg]<< endl;
      for (Int_t Coll=0; Coll<2; Coll++){ //loop on: ppMB + ppHM, pp5TeV
	pathin = "Compare" +YieldOrAvgPt[PlotType] + "DifferentCollisions";
	pathin +=  CollisionsComp[Coll];
	pathin +="_"+ ParticleType[type] + Region[ireg]+ ".root";
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
	StyleHisto(histoYieldRatio[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], 33, "", "N_{#Xi}/N_{K_0^S}", "" , 0,0, 45);
	StyleHisto(histoYieldRatioSist[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], 33, "", "N_{#Xi}/N_{K_0^S}", "" , 0,0, 45);
	histoYieldRatio[ireg][Coll]->Draw("same");
	histoYieldRatioSist[ireg][Coll]->SetFillStyle(0);
	histoYieldRatioSist[ireg][Coll]->Draw("same e2");
      }
    }
  }

  if (ChosenRegion<0){
    canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+".pdf");
  }
  else     canvas->SaveAs("XiK0s3Systems"+SPlotType[PlotType]+Region[ChosenRegion]+".pdf");
}
