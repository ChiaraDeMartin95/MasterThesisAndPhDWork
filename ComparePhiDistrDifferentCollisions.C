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
#include <TLine.h>
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>


void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.8);
  histo->SetTitle(title);
}

void ComparePhiDistrDifferentCollisions(Int_t TypeAnalysis=0, Int_t NSystems =2, Int_t type=0){

  const Int_t nummolt=5;
  const Int_t numPtV0 = 9;

  TString Region[3] = {"Jet", "Bulk", "All"};
  TString RegionName[3] = {"BulkSub", "Bulk", ""};
  Float_t Up[3] = {0.008, 0.014, 0.014};
  Float_t UpSpectrum[3] = {0.02, 0.014, 0.014};
  Float_t Low[3] = {-0.0004, -0.0001, -0.0001};
  Int_t Color[nummolt+1] = {882, 909, 634, 810, 797, 1};
  TFile *file[NSystems];
  TFile *fileSpectra[NSystems];
  TString pathin[NSystems] = {"", ""}; // pp MB - pp HM 
  pathin[0] = "FinalOutput/DATA2016/histo/AngularCorrelation1617_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppMB
  pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppHM

  TString pathinSpectra[NSystems] = {"", ""}; // pp MB - pp HM 
  pathinSpectra[0] = "FinalOutput/DATA2016/SystematicAnalysis1617_hK0s_PtBinning1_K0s_Eta0.8_";
  pathinSpectra[0] += Region[TypeAnalysis];
  pathinSpectra[0] += "Data_PtMin3.0_IsEtaEff.root"; //pp MB                                                      
  pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysis2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_";
  pathinSpectra[1] += Region[TypeAnalysis];
  pathinSpectra[1] += "Data_PtMin3.0_IsEtaEff.root";

  TString NameHisto[NSystems][nummolt+1][numPtV0];
  TString NameHistoSpectrum[NSystems][nummolt+1];
  TString NameHistoSpectrumFinal[NSystems][nummolt+1];

  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  canvas->Divide(5,2);
  gStyle->SetOptStat(0);

  TCanvas * canvasSpectrum = new TCanvas("canvasSpectrum", "canvasSpectrum", 1300, 800);
  gStyle->SetOptStat(0);
  canvasSpectrum->Divide(2,1);

  TLegend * legendMult = new TLegend(0.3, 0.3, 0.9, 0.9);

  TH1F * histoPhiDistr[NSystems][nummolt+1][numPtV0];
  TH1F * histoSpectrum[NSystems][nummolt+1];
  TH1F * histoSpectrumMB;
  TH1F * histoSpectrumFinal[NSystems][nummolt+1];
  TH1F * histoSpectrumMBFinal;
  TH1F * histoSpectrumRatio[NSystems][nummolt+1];

  //NPtV0
  Double_t NPtV0[numPtV0+1]={0};
  TString  SPtV0[numPtV0+1];
  Double_t NPtV0_V0[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString  SPtV0_V0[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0_Casc[numPtV0+1]={0,0.5,1,1.5,2,2.5,3,4,8,8};
  TString  SPtV0_Casc[numPtV0]={"", "0.5-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  for(Int_t v=0; v<=numPtV0; v++){
    if (type==0){
      NPtV0[v] = NPtV0_V0[v];
      if (v<numPtV0) SPtV0[v] = SPtV0_V0[v];
    }
    else{
      NPtV0[v] = NPtV0_Casc[v];
      if (v<numPtV0) SPtV0[v] = SPtV0_Casc[v];
    }
  }
  Int_t PtV0Min = 0;
  if (type!=0) PtV0Min = 1;
  Int_t numPtV0Max = numPtV0;
  if (type!=0) numPtV0Max = numPtV0-1;

  //Molt
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100};

  TLine *lineAt0 = new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2, 0);
  if (TypeAnalysis==0) lineAt0 = new TLine(-1.1, 0, 1.1, 0);
  lineAt0->SetLineColor(1);

  for (Int_t i=0; i< NSystems; i++){
    file[i] = new TFile(pathin[i], "");
    if (i==1) {
      Nmolt[1] = 0.001;
      Nmolt[2] = 0.005;
      Nmolt[3] = 0.01;
      Nmolt[4] = 0.05;
      Nmolt[5] = 0.1;
      Smolt[0] = "0-0.001";
      Smolt[1] = "0.001-0.005";
      Smolt[2] = "0.005-0.01";
      Smolt[3] = "0.01-0.05";
      Smolt[4] = "0.05-0.1";
      Color[0] = 416;
      Color[1] = 844;
      Color[2] = 861;
      Color[3] = 857;
      Color[4] = 601;
    }
    for (Int_t m=0; m< nummolt; m++){
      if (i==1 && m<3 && TypeAnalysis==0) continue;
      for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	NameHisto[i][m][v] = "ME_m"+ Smolt[m]+ "_v"+ SPtV0[v] + "_AC_phi_V0Sub_"+RegionName[TypeAnalysis] +"_EffCorr_TrCorr";
	histoPhiDistr[i][m][v] = (TH1F*) file[i]->Get(NameHisto[i][m][v]);
	if (!histoPhiDistr[i][m][v] ) return;

	Float_t DeltaPt = NPtV0[v+1] - NPtV0[v];
	histoPhiDistr[i][m][v] ->Scale(1./DeltaPt);
	if (TypeAnalysis==0)	histoPhiDistr[i][m][v]->GetXaxis()->SetRangeUser(-1.1, 1.1);
	StyleHisto(histoPhiDistr[i][m][v], Low[TypeAnalysis] , Up[TypeAnalysis], Color[m], 33, "#Delta#varphi", "dN/dp_{T} 1/N_{Trigg} 1/#Delta#eta", Form("%.1f < p_{T} < %.1f GeV/c", NPtV0[v], NPtV0[v+1]), 0,0,0);

	canvas->cd(v+1);
	gPad->SetLeftMargin(0.15);
	histoPhiDistr[i][m][v]->Draw("same");
	if (i==0 && m==0)	lineAt0->Draw("same");
      }
    }
  }

  //legend
  
  for (Int_t i=1; i>=0; i--){
    for (Int_t m=0; m< nummolt; m++){
      if (i==0){
	Smolt[0] = "0-5";
	Smolt[1] = "5-10";
	Smolt[2] = "10-30";
	Smolt[3] = "30-50";
	Smolt[4] = "50-100";
      }
      if (i==1 && m<3 && TypeAnalysis==0) continue;
      cout << " m " << Smolt[m] << endl;
      legendMult->AddEntry(histoPhiDistr[i][m][PtV0Min], Smolt[m] , "pl");
    }
  }
  canvas->cd(10);
  legendMult->Draw();
  
  //comparison of spectra
  for (Int_t i=0; i< NSystems; i++){
    fileSpectra[i] = new TFile(pathinSpectra[i], "");
    if (i==0){
      Color[0] = 882;
      Color[1] = 909;
      Color[2] = 634;
      Color[3] = 810;
      Color[4] = 797;
    }
    else if (i==1) {
      Nmolt[1] = 0.001;
      Nmolt[2] = 0.005;
      Nmolt[3] = 0.01;
      Nmolt[4] = 0.05;
      Nmolt[5] = 0.1;
      Smolt[0] = "0-0.001";
      Smolt[1] = "0.001-0.005";
      Smolt[2] = "0.005-0.01";
      Smolt[3] = "0.01-0.05";
      Smolt[4] = "0.05-0.1";
      Color[0] = 416;
      Color[1] = 844;
      Color[2] = 861;
      Color[3] = 857;
      Color[4] = 601;
    }
    for (Int_t m=0; m<=nummolt; m++){
      if (i==1 && m<3 && TypeAnalysis==0) continue;
      if (i==1 && m==nummolt) continue; //I am interested in 0-100% only
      cout << Smolt[m] << endl;
      NameHistoSpectrum[i][m] = Form("fHistSpectrumPart_m%i_syst0", m);
      NameHistoSpectrumFinal[i][m] = "fHistSpectrumPart_m_" + Smolt[m] + "syst0";
      histoSpectrum[i][m] = (TH1F*) fileSpectra[i]->Get(NameHistoSpectrum[i][m]);
      if (!histoSpectrum[i][m]) {cout << NameHistoSpectrum[i][m] << endl; return;}
      histoSpectrumFinal[i][m] = (TH1F*)       histoSpectrum[i][m]->Clone(NameHistoSpectrumFinal[i][m]);
      histoSpectrumFinal[i][m] ->Sumw2();

      if (m==0 && i==0){
	NameHistoSpectrum[i][nummolt] = Form("fHistSpectrumPart_m%i_syst0", nummolt);
	histoSpectrumMB = (TH1F*) fileSpectra[i]->Get(NameHistoSpectrum[i][nummolt]);
	if (!histoSpectrumMB) {cout << NameHistoSpectrum[i][nummolt] << endl; return;}
	histoSpectrumMBFinal = (TH1F*)       histoSpectrumMB->Clone("fhistoSpectrumMB");
	histoSpectrumMBFinal->Sumw2();
      }

      histoSpectrumRatio[i][m] = (TH1F*)       histoSpectrumFinal[i][m]->Clone(NameHistoSpectrumFinal[i][m]+ "_Ratio");
      histoSpectrumRatio[i][m] ->Divide(histoSpectrumMBFinal);

      StyleHisto(histoSpectrumFinal[i][m],3*10e-4, UpSpectrum[TypeAnalysis], Color[m], 33, "p_{T}", "dN/dp_{T} 1/N_{Trigg} 1/#Delta#eta #Delta#varphi", "", 0,0,0);

      canvasSpectrum->cd(1);
      gPad->SetLeftMargin(0.15);
      histoSpectrumFinal[i][m]->Draw("same");

      StyleHisto(histoSpectrumRatio[i][m],0.5, 1.5, Color[m], 33, "p_{T}", "Ratio to MB", "", 0,0,0);

      canvasSpectrum->cd(2);
      gPad->SetLeftMargin(0.15);
      if (!(m==nummolt && i==0))      histoSpectrumRatio[i][m]->Draw("same");

    }
  }

  TString outname = "ComparePhiDistrDifferentCollisions" + Region[TypeAnalysis] +".root";
  TFile * fileout = new TFile (outname, "RECREATE");
  canvas->Write();
  canvasSpectrum->Write();
  fileout->Close();
  canvas->SaveAs("ComparePhiDistrDifferentCollisions"+ Region[TypeAnalysis]+".pdf");
  cout << " I have produced the file ComparePhiDistrDifferentCollisions"<< Region[TypeAnalysis] << ".root" << endl;
}
