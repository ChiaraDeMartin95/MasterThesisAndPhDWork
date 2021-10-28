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

void RatiosXiK0s(Float_t ScalingFactorXiK0s = 0.8458/1.08747){

  gStyle->SetOptStat(0);
  TString Region[3] = {"Jet", "Bulk", "All"};
  TString RegionBis[3] = {"Jet", "Bulk", "Inclusive"};
  Float_t Up= 0.1;
  Float_t Low = 10e-4;
  Int_t Color[3] = {628,418,601};
  TString ParticleType[2]={"K0s", "Xi"};
  TH1F *histoYield[2][3];
  TH1F *histoYieldRatio[3];
  TH1F *histoYieldSist[2][3];
  TH1F *histoYieldRatioSist[3];
  TString NameHisto="histoYieldComparison";
  TString NameHistoSist="histoYieldSistComparison";
  TString NameHistoFinal[2][3];
  TString pathin ="";

  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);

  for (Int_t type=1; type>=0; type--){
    for (Int_t ireg=0; ireg<3; ireg++){
      cout << "\n\n*****" << ParticleType[type] << " " << Region[ireg]<< endl;
      pathin = "CompareYieldDifferentCollisions_HMMultBinning1";
      pathin += ParticleType[type];
      pathin += Region[ireg];
      pathin += ".root";
      cout << "Pathin: " << pathin << endl;
      TFile *filein = new TFile(pathin, "");
      if (!filein) {cout << "Input file not available " << endl; return;}
      NameHistoFinal[type][ireg]= Form("histoYield_Reg%i_Type%i", ireg, type);
      histoYield[type][ireg]= (TH1F*) filein->Get(NameHisto);
      if (!histoYield[type][ireg]) {cout <<"no histo " << endl; return;}
      histoYield[type][ireg]->SetName(NameHistoFinal[type][ireg]);
      histoYieldSist[type][ireg]= (TH1F*) filein->Get(NameHistoSist);
      if (!histoYieldSist[type][ireg]) {cout <<"no histo sist" << endl; return;}
      histoYieldSist[type][ireg]->SetName(NameHistoFinal[type][ireg]+"Sist");

      if (type==1){
	histoYieldRatio[ireg]= (TH1F*)      histoYield[type][ireg]->Clone(NameHisto[ireg]+"Ratio");
	histoYieldRatioSist[ireg]= (TH1F*)      histoYieldSist[type][ireg]->Clone(NameHisto[ireg]+"RatioSist");
	histoYieldRatio[ireg]->Sumw2();
	histoYieldRatioSist[ireg]->Sumw2();
      }
      else {
	histoYieldRatio[ireg]->Divide(histoYield[type][ireg]);
	histoYieldRatioSist[ireg]->Divide(histoYieldSist[type][ireg]);
	if (ireg ==0) {
	  histoYieldRatio[ireg]->Scale(ScalingFactorXiK0s);
	  histoYieldRatioSist[ireg]->Scale(ScalingFactorXiK0s);
	}
	canvas->cd();
	StyleHisto(histoYieldRatio[ireg], Low, Up, Color[ireg], 33, "", "N_{#Xi}/N_{K_0^S}", "" , 0,0, 45);
	StyleHisto(histoYieldRatioSist[ireg], Low, Up, Color[ireg], 33, "", "N_{#Xi}/N_{K_0^S}", "" , 0,0, 45);
	histoYieldRatio[ireg]->Draw("same");
	histoYieldRatioSist[ireg]->SetFillStyle(0);
	histoYieldRatioSist[ireg]->Draw("same e2");
      }
    }
  }
}
