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

void ComparePhiDistrDifferentCollisions(Int_t TypeAnalysis=0,   const Int_t NSystems =2,  Int_t type=0, Int_t mRatio=0){

  const Int_t nummolt=5;
  const Int_t numPtV0 = 9;
  if (mRatio > nummolt) {cout << " mRatio should lie between 0 and 5 " << endl; return;}

  TString Region[3] = {"Jet", "Bulk", "All"};
  TString RegionName[3] = {"BulkSub_", "Bulk_", ""};
  Float_t Up[3] = {0.008, 0.03, 0.03};
  Float_t UpSpectrum[3] = {0.02, 0.25, 0.25};
  Float_t UpSpectrumRatio[3] = {1.5, 5, 5};
  Float_t LowSpectrumRatio[3] = {0.5, 0, 0};
  Float_t UpEff[3] = {0.4, 0.4, 0.4};
  Float_t Low[3] = {-0.0004, -0.0001, -0.0001};
  Int_t Color[nummolt+1] = {882, 909, 634, 810, 797, 1};
  TFile *file[NSystems];
  TFile *fileSpectra[NSystems];
  TFile *fileEff[NSystems];
  TString pathin[3] = {"", "", ""}; // pp MB - pp HM 
  pathin[0] = "FinalOutput/DATA2016/histo/AngularCorrelation1617_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppMB
  //  pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppHM
  pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelationAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppHM

  TString pathinSpectra[3] = {"", "", ""}; // pp MB - pp HM 
  pathinSpectra[0] = "FinalOutput/DATA2016/SystematicAnalysis1617_hK0s_PtBinning1_K0s_Eta0.8_";
  pathinSpectra[0] += Region[TypeAnalysis];
  pathinSpectra[0] += "Data_PtMin3.0_IsEtaEff.root"; //pp MB                                                      
  //  pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysis2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_";
  pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysisAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_";
  pathinSpectra[1] += Region[TypeAnalysis];
  pathinSpectra[1] += "Data_PtMin3.0_IsEtaEff.root";

  TString pathinEff[3] = {"", "", ""}; // pp MB - pp HM 
  pathinEff[0] = "FinalOutput/DATA2016/Efficiency/Efficiency1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root"; //ppMB
  //pathinEff[0] = "FinalOutput/DATA2016/Efficiency/Efficiency2018g4_extra_EtaEff_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
  pathinEff[1] = "FinalOutput/DATA2016/Efficiency/Efficiency2019h11c_extra_HM_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root"; //ppHM
  pathinEff[2] = "FinalOutput/DATA2016/Efficiency/Efficiency1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0_HM.root"; //ppHM 0-1% of the MB

  cout << " starting from files for MB: " << endl;
  cout << pathin[0] << "\n" << pathinSpectra[0] << "\n" << pathinEff[0] << endl;
  cout << " and from files for HM: " << endl;
  cout << pathin[1] << "\n" << pathinSpectra[1] << "\n" << pathinEff[1] << endl;

  TString NameHisto[NSystems][nummolt+1][numPtV0];
  TString NameHistoSpectrum[NSystems][nummolt+1];
  TString NameHistoSpectrumFinal[NSystems][nummolt+1];
  TString NameHistoEff[NSystems][nummolt+1];
  TString NameHistoEffFinal[NSystems][nummolt+1];
  TString NameHistoEffEta[NSystems][nummolt+1];
  TString NameHistoEffEtaFinal[NSystems][nummolt+1];
  TString NameHistoEff2D[NSystems][nummolt+1];
  TString NameHistoEff2DFinal[NSystems][nummolt+1];

  TCanvas * canvas= new TCanvas("canvas", "canvas", 1300, 800);
  canvas->Divide(5,2);
  gStyle->SetOptStat(0);

  TCanvas * canvasRatio= new TCanvas("canvasRatio", "canvasRatio", 1300, 800);
  canvasRatio->Divide(5,2);
  gStyle->SetOptStat(0);

  TCanvas * canvasSpectrum = new TCanvas("canvasSpectrum", "canvasSpectrum", 1300, 800);
  gStyle->SetOptStat(0);
  canvasSpectrum->Divide(2,1);

  TCanvas * canvasEff= new TCanvas("canvasEff", "canvasEff", 1300, 800);
  canvasEff->Divide(2,1);
  gStyle->SetOptStat(0);

  TCanvas * canvasEffEta= new TCanvas("canvasEffEta", "canvasEffEta", 1300, 800);
  canvasEffEta->Divide(2,1);
  gStyle->SetOptStat(0);

  TCanvas * canvasEff2D= new TCanvas("canvasEff2D", "canvasEff2D", 1300, 800);
  canvasEff2D -> Divide(5,2);
  TCanvas * canvasEff2DRatio= new TCanvas("canvasEff2DRatio", "canvasEff2DRatio", 1300, 800);
  canvasEff2DRatio -> Divide(5,2);

  TLegend * legendMult = new TLegend(0.3, 0.3, 0.9, 0.9);

  TH1F * histoPhiDistr[NSystems][nummolt+1][numPtV0];
  TH1F * histoPhiDistrRatio[NSystems][nummolt+1][numPtV0];
  TH1F * histoPhiDistrDenom[NSystems];
  TH1F * histoSpectrum[NSystems][nummolt+1];
  TH1F * histoSpectrumMB;
  TH1F * histoSpectrumFinal[NSystems][nummolt+1];
  TH1F * histoSpectrumMBFinal;
  TH1F * histoSpectrumRatio[NSystems][nummolt+1];
  TH1F * histoEff[NSystems][nummolt+1];
  TH1F * histoEffMB;
  TH1F * histoEffFinal[NSystems][nummolt+1];
  TH1F * histoEffMBFinal;
  TH1F * histoEffRatio[NSystems][nummolt+1];

  TH1F * histoEffEta[NSystems][nummolt+1];
  TH1F * histoEffEtaMB;
  TH1F * histoEffEtaFinal[NSystems][nummolt+1];
  TH1F * histoEffEtaMBFinal;
  TH1F * histoEffEtaRatio[NSystems][nummolt+1];

  TH2F * histoEff2D[NSystems][nummolt+1];
  TH2F * histoEff2DMB;
  TH2F * histoEff2DFinal[NSystems][nummolt+1];
  TH2F * histoEff2DMBFinal;
  TH2F * histoEff2DRatio[NSystems][nummolt+1];

  //NPtV0
  Double_t NPtV0[numPtV0+1]={0};
  TString  SPtV0[numPtV0];
  Double_t NPtV0_V0[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString  SPtV0_V0[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0_Casc[numPtV0+1]={0,0.5,1,1.5,2,2.5,3,4,8,8};
  TString  SPtV0_Casc[numPtV0]={"", "0.5-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  for(Int_t v=0; v<=numPtV0; v++){
    if (type==0){
      NPtV0[v] = NPtV0_V0[v];
      if (v<numPtV0) {
	SPtV0[v] = SPtV0_V0[v];
      }
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

  TLine *lineAt1 = new TLine(-TMath::Pi()/2, 1, 3*TMath::Pi()/2, 1);
  if (TypeAnalysis==0) lineAt1 = new TLine(-1.1, 1, 1.1, 1);
  lineAt1->SetLineColor(1);

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
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      //      if (i==1 && m<3 && TypeAnalysis==1) continue;
      for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	NameHisto[i][m][v] = "ME_m"+ Smolt[m]+ "_v"+ SPtV0[v] + "_AC_phi_V0Sub_"+RegionName[TypeAnalysis] +"EffCorr_TrCorr";
	if (i==0) NameHisto[i][mRatio][v] = "ME_m"+ Smolt[mRatio]+ "_v"+ SPtV0[v] + "_AC_phi_V0Sub_"+RegionName[TypeAnalysis] +"EffCorr_TrCorr";
	histoPhiDistr[i][m][v] = (TH1F*) file[i]->Get(NameHisto[i][m][v]);
	if (!histoPhiDistr[i][m][v] ) {cout << " phi distr file not found: " << NameHisto[i][m][v] << endl; return;}
	//histoPhiDistr[i][m][v]->Sumw2();

	Float_t DeltaPt = NPtV0[v+1] - NPtV0[v];
	histoPhiDistr[i][m][v] ->Scale(1./DeltaPt);
	if (TypeAnalysis==0)	histoPhiDistr[i][m][v]->GetXaxis()->SetRangeUser(-1.1, 1.1);
	StyleHisto(histoPhiDistr[i][m][v], Low[TypeAnalysis] , Up[TypeAnalysis], Color[m], 33, "#Delta#varphi", "dN/dp_{T} 1/N_{Trigg} 1/#Delta#eta", Form("%.1f < p_{T} < %.1f GeV/c", NPtV0[v], NPtV0[v+1]), 0,0,0);

	canvas->cd(v+1);
	gPad->SetLeftMargin(0.15);
	Float_t NormFactor = histoPhiDistr[i][m][v]->Integral();
	if (TypeAnalysis==1) {
	  histoPhiDistr[i][m][v]->Rebin(2);
	  NormFactor = histoPhiDistr[i][m][v]->Integral(1,2);
	  histoPhiDistr[i][m][v]->Scale(1./NormFactor);
	  histoPhiDistr[i][m][v]->GetYaxis()->SetRangeUser(0, 3); //0, 0.06 if full integral is employed
	}
	histoPhiDistr[i][m][v]->DrawClone("same");
	if (i==0 && m==0)	lineAt0->Draw("same");
	
	if (i==0 && m ==0){
	  histoPhiDistrDenom[v] = (TH1F*) file[i]->Get(NameHisto[i][mRatio][v]);
	  histoPhiDistrDenom[v]->SetName("PhiDistrRatio"+ SPtV0[v]);
	  histoPhiDistrDenom[v]->Rebin(2);
	  Float_t NormFactorDenom = histoPhiDistrDenom[v]->Integral();
	  if (TypeAnalysis==1) NormFactorDenom = histoPhiDistrDenom[v]->Integral(1,2);
	  histoPhiDistrDenom[v]->Scale(1./NormFactorDenom);
	}

	if(TypeAnalysis!=1)	{
	  histoPhiDistr[i][m][v]->Rebin(2);
	  histoPhiDistr[i][m][v]->Scale(1./NormFactor);
	}

	histoPhiDistrRatio[i][m][v] = (TH1F*) histoPhiDistr[i][m][v] -> Clone(NameHisto[i][m][v] + "Ratio");
	histoPhiDistrRatio[i][m][v]->Divide(histoPhiDistrDenom[v]);

	canvasRatio->cd(v+1);
	gPad->SetLeftMargin(0.15);
	histoPhiDistrRatio[i][m][v]->GetYaxis()->SetRangeUser(0.5, 1.5);
	histoPhiDistrRatio[i][m][v]->DrawClone("same");
	if (i==0 && m==0)	lineAt1->Draw("same");

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
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      //      if (i==1 && m<3 && TypeAnalysis==1) continue;
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
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      //      if (i==1 && m<3 && TypeAnalysis==1) continue;
      if (i==1 && m==nummolt) continue; //I am interested in 0-100% only
      NameHistoSpectrum[i][m] = Form("fHistSpectrumPart_m%i_syst0", m);
      NameHistoSpectrumFinal[i][m] = "fHistSpectrumPart_m_" + Smolt[m] + "syst0";
      histoSpectrum[i][m] = (TH1F*) fileSpectra[i]->Get(NameHistoSpectrum[i][m]);
      if (!histoSpectrum[i][m]) {cout << NameHistoSpectrum[i][m] << endl; return;}
      histoSpectrumFinal[i][m] = (TH1F*)       histoSpectrum[i][m]->Clone(NameHistoSpectrumFinal[i][m]);
      //      histoSpectrumFinal[i][m] ->Sumw2();

      if (m==0 && i==0){
	NameHistoSpectrum[i][mRatio] = Form("fHistSpectrumPart_m%i_syst0", mRatio);
	histoSpectrumMB = (TH1F*) fileSpectra[i]->Get(NameHistoSpectrum[i][mRatio]);
	if (!histoSpectrumMB) {cout << NameHistoSpectrum[i][mRatio] << endl; return;}
	histoSpectrumMBFinal = (TH1F*)       histoSpectrumMB->Clone("fhistoSpectrumMB");
	//	histoSpectrumMBFinal->Sumw2();
      }

      histoSpectrumRatio[i][m] = (TH1F*)       histoSpectrumFinal[i][m]->Clone(NameHistoSpectrumFinal[i][m]+ "_Ratio");
      histoSpectrumRatio[i][m] ->Divide(histoSpectrumMBFinal);

      StyleHisto(histoSpectrumFinal[i][m],3*10e-4, UpSpectrum[TypeAnalysis], Color[m], 33, "p_{T}", "dN/dp_{T} 1/N_{Trigg} 1/#Delta#eta #Delta#varphi", "", 0,0,0);

      canvasSpectrum->cd(1);
      gPad->SetLeftMargin(0.15);
      histoSpectrumFinal[i][m]->Draw("same");

      StyleHisto(histoSpectrumRatio[i][m],LowSpectrumRatio[TypeAnalysis], UpSpectrumRatio[TypeAnalysis], Color[m], 33, "p_{T}", "Ratio to MB", "", 0,0,0);

      canvasSpectrum->cd(2);
      gPad->SetLeftMargin(0.15);
      if (!(m==mRatio && i==0))      histoSpectrumRatio[i][m]->Draw("same");

    }
  }

  //comparison of efficiencies
  for (Int_t i=0; i< NSystems; i++){
    fileEff[i] = new TFile(pathinEff[i], "");
    if (i==0){
      Color[0] = 882;
      Color[1] = 909;
      Color[2] = 634;
      Color[3] = 810;
      Color[4] = 797;
      Smolt[0] = "0-5";
      Smolt[1] = "5-10";
      Smolt[2] = "10-30";
      Smolt[3] = "30-50";
      Smolt[4] = "50-100";

    }
    else if (i==1 || i==2) {
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
    if (i==2) Color[0] = 922;
    for (Int_t m=0; m<=nummolt; m++){
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      if (i==1 && m==nummolt) continue; //I am interested in 0-100% only
      if (i==2 && m!=0) continue; //I am interested in 0-1% only 

      NameHistoEff[i][m] = "fHistV0EfficiencyPtBins_"+Smolt[m];
      NameHistoEffFinal[i][m] = "fHistEff_m_" + Smolt[m];
      histoEff[i][m] = (TH1F*) fileEff[i]->Get(NameHistoEff[i][m]);
      if (!histoEff[i][m]) {cout << "1D eff " << NameHistoEff[i][m] << endl; return;}
      histoEffFinal[i][m] = (TH1F*)       histoEff[i][m]->Clone(NameHistoEffFinal[i][m]);

      NameHistoEffEta[i][m] = "fHistV0EfficiencyEta_"+Smolt[m];
      NameHistoEffEtaFinal[i][m] = "fHistEffEta_m_" + Smolt[m];
      histoEffEta[i][m] = (TH1F*) fileEff[i]->Get(NameHistoEffEta[i][m]);
      if (!histoEffEta[i][m]) {cout << "1D eff " << NameHistoEffEta[i][m] << endl; return;}
      histoEffEtaFinal[i][m] = (TH1F*)       histoEffEta[i][m]->Clone(NameHistoEffEtaFinal[i][m]);

      if (!(pathinEff[0].Index("2018g4_extra") !=1)) {
	NameHistoEff2D[i][m] = "fHistV0EfficiencyPtV0EtaV0PtBins_"+Smolt[m];
	NameHistoEff2DFinal[i][m] = "fHistEff2D_m_" + Smolt[m];
	histoEff2D[i][m] = (TH2F*) fileEff[i]->Get(NameHistoEff2D[i][m]);
	if (!histoEff2D[i][m]) {cout << NameHistoEff2D[i][m] << endl; return;}
	histoEff2DFinal[i][m] = (TH2F*)       histoEff2D[i][m]->Clone(NameHistoEff2DFinal[i][m]);
	histoEff2DFinal[i][m] -> GetZaxis()->SetRangeUser(0,0.5);
	cout << " i " << i << " m" << m << endl;
	cout << " num bins x axis " << histoEff2DFinal[i][m]->GetNbinsX() << endl;
	cout << " num bins y axis " << histoEff2DFinal[i][m]->GetNbinsY() << endl;
	//      histoEffFinal[i][m] ->Sumw2();
      }

      if (m==0 && i==0){
	NameHistoEff[i][mRatio] = "fHistV0EfficiencyPtBins_"+Smolt[mRatio];
	histoEffMB = (TH1F*) fileEff[i]->Get(NameHistoEff[i][mRatio]);
	if (!histoEffMB) {cout << NameHistoEff[i][mRatio] << endl; return;}
	histoEffMBFinal = (TH1F*)       histoEffMB->Clone("fhistoEffMB");
	//	histoEffMBFinal->Sumw2();
	NameHistoEffEta[i][mRatio] = "fHistV0EfficiencyEta_"+Smolt[mRatio];
	histoEffEtaMB = (TH1F*) fileEff[i]->Get(NameHistoEffEta[i][mRatio]);
	if (!histoEffEtaMB) {cout << NameHistoEffEta[i][mRatio] << endl; return;}
	histoEffEtaMBFinal = (TH1F*)       histoEffEtaMB->Clone("fhistoEffMB");

	if (!(pathinEff[0].Index("2018g4_extra") !=1)) {
	  NameHistoEff2D[i][mRatio] = "fHistV0EfficiencyPtV0EtaV0PtBins_"+Smolt[mRatio];
	  histoEff2DMB = (TH2F*) fileEff[i]->Get(NameHistoEff2D[i][mRatio]);
	  if (!histoEff2DMB) {cout << NameHistoEff2D[i][mRatio] << endl; return;}
	  histoEff2DMBFinal = (TH2F*)       histoEff2DMB->Clone("fhistoEff2DMB");
	}
      }

      histoEffRatio[i][m] = (TH1F*)       histoEffFinal[i][m]->Clone(NameHistoEffFinal[i][m]+ "_Ratio");
      histoEffRatio[i][m] ->Divide(histoEffMBFinal);

      histoEffEtaRatio[i][m] = (TH1F*)       histoEffEtaFinal[i][m]->Clone(NameHistoEffEtaFinal[i][m]+ "_Ratio");
      histoEffEtaRatio[i][m] ->Divide(histoEffEtaMBFinal);

      if (!(pathinEff[0].Index("2018g4_extra") !=1)) {
	histoEff2DRatio[i][m] = (TH2F*)       histoEff2DFinal[i][m]->Clone(NameHistoEff2DFinal[i][m]+ "_Ratio");
	histoEff2DRatio[i][m] ->Divide(histoEff2DMBFinal);
	histoEff2DRatio[i][m] -> GetZaxis()->SetRangeUser(0.8, 1.2);
      }

      StyleHisto(histoEffFinal[i][m],3*10e-4, UpEff[TypeAnalysis], Color[m], 33, "p_{T}", "Efficiency", "", 0,0,0);
      StyleHisto(histoEffEtaFinal[i][m],3*10e-4, UpEff[TypeAnalysis], Color[m], 33, "#eta", "Efficiency", "", 0,0,0);

      canvasEff->cd(1);
      gPad->SetLeftMargin(0.15);
      histoEffFinal[i][m]->Draw("same");
      legendMult->Draw();

      canvasEffEta->cd(1);
      gPad->SetLeftMargin(0.15);
      histoEffEtaFinal[i][m]->Draw("same p");
      legendMult->Draw();

      StyleHisto(histoEffRatio[i][m],0.5, 1.5, Color[m], 33, "p_{T}", "Ratio to MB", "", 0,0,0);
      StyleHisto(histoEffEtaRatio[i][m],0.5, 1.5, Color[m], 33, "#eta", "Ratio to MB", "", 0,0,0);

      canvasEff->cd(2);
      gPad->SetLeftMargin(0.15);
      if (!(m==mRatio && i==0))      histoEffRatio[i][m]->Draw("same");

      canvasEffEta->cd(2);
      gPad->SetLeftMargin(0.15);
      if (!(m==mRatio && i==0))      histoEffEtaRatio[i][m]->Draw("same p");

      if (!(pathinEff[0].Index("2018g4_extra") !=1)) {
	if (i==0)   canvasEff2D->cd(nummolt-m);
	else   canvasEff2D->cd(nummolt+m+1);
	gPad->SetLeftMargin(0.15);
	histoEff2DFinal[i][m]->GetYaxis()->SetRangeUser(-0.8, 0.8);
	histoEff2DFinal[i][m]->Draw("same colz");

	if (i==0)   canvasEff2DRatio->cd(nummolt-m);
	else   canvasEff2DRatio->cd(nummolt+m+1);
	gPad->SetLeftMargin(0.15);
	histoEff2DRatio[i][m]->GetYaxis()->SetRangeUser(-0.8, 0.8);
	histoEff2DRatio[i][m]->Draw("same colz");
      }
    }
  }

  TString NameFileout ="ComparePhiDistrDifferentCollisions";
  if (pathin[1].Index("2016k_HM_hK0s")!=-1) NameFileout += "_HM16k";
  if (pathinEff[0].Index("18g4_extra")!=-1) NameFileout += "_MBEff18g4extra";
  NameFileout += Region[TypeAnalysis]+ ".root";

  TFile * fileout = new TFile (NameFileout, "RECREATE");
  canvas->Write();
  canvasRatio->Write();
  canvasSpectrum->Write();
  canvasEff->Write();
  canvasEffEta->Write();
  canvasEff2D->Write();
  canvasEff2DRatio->Write();
  fileout->Close();
  canvas->SaveAs("ComparePhiDistrDifferentCollisions"+ Region[TypeAnalysis]+".pdf");

  cout << " I have produced the file " << NameFileout << endl;
}
