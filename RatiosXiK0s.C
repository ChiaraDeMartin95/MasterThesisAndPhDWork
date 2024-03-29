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

void RatiosXiK0s( Int_t isComp=0, Float_t ScalingFactorXiK0s = 1/*0.8458/1.08747*/, Int_t PlotType =0){

  //PlotType = 0: Xi/K0s ratio
  //PlotType = 1: K0s yield vs mult 
  //PlotType = 2: Xi yield vs mult 
  //PlotType = 3: K0s pt vs mult 
  //PlotType = 4: Xi pt vs mult 

  //isComp =0 : MB compared to HM
  //isComp =1 : MB compared to pp5TeV (same mult classes)
  //isComp =2 : MB compared to pp5TeV (only two mult classes for 5TeV)

  TString CollisionsComp[3] = {"_HMMultBinning1_vsHM", "_vs5TeV5Mult", "_vs5TeV"};
  TString SPlotType[5] = {"", "K0s", "Xi", "K0spt", "Xipt"};

  gStyle->SetOptStat(0);
  TString Region[3] = {"Jet", "Bulk", "All"};
  TString RegionBis[3] = {"Jet", "Bulk", "Inclusive"};
  Float_t Up= 0.1;
  Float_t Low = 10e-4;
  if (PlotType==1) {Low = 10e-8; Up = 0.3;}
  if (PlotType==2) {Low = 10e-8; Up = 0.025;}
  if (PlotType==3) {Low = 0.8+10e-4; Up = 2.4999;}
  if (PlotType==4) {Low = 1+10e-4; Up = 4.999;}
  Int_t Color[3] = {628,418,601};
  Int_t ColorDiff[3][2] = {{628,628}, {418,418} , {601, 601}};
  if (isComp!=0){
    ColorDiff[0][1] = 628;
    if (isComp==2)     ColorDiff[0][0] = 634;
    ColorDiff[1][1] = 829;
    ColorDiff[2][1] = 867;
  }
  Int_t Marker[2] = {33, 21};
  Float_t Size[2] = {2, 1.4};
  TString ParticleType[2]={"K0s", "Xi"};
  TH1F *histoYield[2][3][2];
  TH1F *histoYieldRatio[3][2];
  TH1F *histoYieldSist[2][3][2];
  TH1F *histoYieldRatioSist[3][2];
  TString NameHisto="histoYieldComparison";
  TString NameHistoSist="histoYieldSistComparison";
  TString NameHistoInput[2] = {"fHistYield_ppMB", "fHistYield_ppHM"};
  TString NameHistoInputSist[2] = {"fHistYield_ppMB_Sist", "fHistYield_ppHM_Sist"};
  if (isComp!=0) {
    NameHistoInput[1] = "fHistYield_pp5TeV";
    NameHistoInputSist[1] = "fHistYield_pp5TeV_Sist";
  }
  TString NameHistoFinal[2][3][2];
  TString pathin ="";
  TString YieldOrAvgPt[5] = {"Yield", "Yield", "Yield", "AvgPt","AvgPt"};
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);

  for (Int_t type=1; type>=0; type--){
    if ((PlotType == 1 || PlotType ==3) && type==1) continue;
    else  if ((PlotType == 2 || PlotType ==4) && type==0) continue;
    for (Int_t ireg=0; ireg<3; ireg++){
      cout << "\n\n*****" << ParticleType[type] << " " << Region[ireg]<< endl;
      pathin = "Compare" +YieldOrAvgPt[PlotType] + "DifferentCollisions";
      pathin +=  CollisionsComp[isComp];
      pathin +="_"+ ParticleType[type] + Region[ireg]+ ".root";
      cout << "Pathin: " << pathin << endl;
      TFile *filein = new TFile(pathin, "");
      if (!filein) {cout << "Input file not available " << endl; return;}
      for (Int_t Coll=0; Coll<2; Coll++){
	if (isComp==0 && Coll==1) continue;
	NameHistoFinal[type][ireg][Coll]= Form("histoYield_Reg%i_Type%i_Coll%i", ireg, type, Coll);
	if (isComp==0)	histoYield[type][ireg][Coll]= (TH1F*) filein->Get(NameHisto);
	else 	histoYield[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoInput[Coll]);
	if (!histoYield[type][ireg][Coll]) {cout <<"no histo " << endl; return;}
	cout <<  histoYield[type][ireg][Coll]->GetXaxis()->GetXmax() << endl;
	histoYield[type][ireg][Coll]->SetName(NameHistoFinal[type][ireg][Coll]);

	if (isComp==0)	histoYieldSist[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoSist);
	else histoYieldSist[type][ireg][Coll]= (TH1F*) filein->Get(NameHistoInputSist[Coll]);
	if (!histoYieldSist[type][ireg][Coll]) {cout <<"no histo " << endl; return;}
	histoYieldSist[type][ireg][Coll]->SetName(NameHistoFinal[type][ireg][Coll]+"Sist");

	cout << "Type: " << ParticleType[type] <<  " region: " << Region[ireg] << " System: " << Coll << endl;
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
	canvas->cd();
	StyleHisto(histoYieldRatio[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], 33, "", "N_{#Xi}/N_{K_0^S}", "" , 0,0, 45);
	StyleHisto(histoYieldRatioSist[ireg][Coll], Low, Up, ColorDiff[ireg][Coll], 33, "", "N_{#Xi}/N_{K_0^S}", "" , 0,0, 45);
	histoYieldRatio[ireg][Coll]->Draw("same");
	histoYieldRatioSist[ireg][Coll]->SetFillStyle(0);
	histoYieldRatioSist[ireg][Coll]->Draw("same e2");
      }
    }
  }

  canvas->SaveAs("XiK0sRatio"+ CollisionsComp[isComp]+SPlotType[PlotType]+".pdf");
}
