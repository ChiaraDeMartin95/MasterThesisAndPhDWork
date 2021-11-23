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

void ComparePhiDistrDifferentCollisions(Int_t TypeAnalysis=0,   const Int_t NSystems =2,  Int_t type=1, Int_t mRatio=0, Int_t isComp=2, Int_t MultBinning=0, Bool_t isNormFactorCorr=1){

  //isComp =0 : MB compared to HM
  //isComp =1 : MB compared to pp5TeV (same mult classes)
  //isComp =2 : MB compared to pp5TeV (only two mult classes for 5TeV)

  const Int_t nummolt=5;
  const Int_t numPtV0 = 9;
  if (mRatio > nummolt) {cout << " mRatio should lie between 0 and 5 " << endl; return;}

  TString Region[3] = {"Jet", "Bulk", "All"};
  TString RegionBis[3] = {"Jet", "Bulk", "Inclusive"};
  TString RegionName[3] = {"BulkSub_", "Bulk_", ""};
  TString SSystem[2] = {"MB", "HM"};
  if (isComp==1 || isComp==2) {
    SSystem[0] = "13 TeV";
    SSystem[1] = "5 TeV";
  }
  TString tipo[2] = {"K0s", "Xi"};

  //renages depend on choice of mRatio
  Float_t Up[3] = {0.008, 1.5, 4};
  Float_t UpRatio[3] = {1.5, 5, 5};
  Float_t LowRatio[3] = {0.5, 0.5, 0.5};
  Float_t LowSpectrum[3] = {10e-8, 10e-8, 10e-8};
  Float_t UpSpectrum[3] = {0.02, 0.25, 0.25};
  Float_t UpSpectrumRatio[3] = {1.5, 5, 5};
  if (mRatio==0 && (isComp==1 || isComp==2)) UpSpectrumRatio[1] = UpSpectrumRatio[2] = 1;
  if (type==1){
    UpSpectrum[1] = UpSpectrum[2] = 0.012;
    if (isComp!=0)     UpSpectrum[0] = 0.001;
  }
  Float_t LowSpectrumRatio[3] = {0.5, 0, 0};
  if (type==1 && isComp!=0) {
    LowSpectrumRatio[0]=0;
    UpSpectrumRatio[0]=2;
  }
  Float_t UpEff[3] = {0.4, 0.4, 0.4};
  Float_t Low[3] = {-0.0004, 0, 0.4};
  //----------------------------------
  Int_t Color0[nummolt+1] = {882, 909, 634, 810, 797, 1};
  Int_t Color1[nummolt+1] = {416, 844, 861, 857, 601, 1};
  Int_t Color[nummolt+1];
  TFile *file[NSystems];
  TFile *fileSpectra[NSystems];
  TFile *fileEff[NSystems];
  TString pathin[3] = {"", "", ""}; // pp MB - pp HM 
  //  pathin[0] = "FinalOutput/DATA2016/histo/AngularCorrelation1617_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppMB, old AODs
  pathin[0] = "FinalOutput/DATA2016/histo/AngularCorrelation1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppMB, new AODs
  //  pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppHM
  //  pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelationAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppHM
  //  pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelationAllhK0sHM_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root"; //ppHM
  if (isComp==0) {
    pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelationAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff_MultBinning1.root"; //ppHM new eff + only 3 mult classes
  }
  else {
    pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelation17pq_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff";
    if (isComp==2) pathin[1] += "_MultBinning3";
    pathin[1]+=".root";
  }

  if (type==1){
    //PRELIMINARY    pathin[0] = "FinalOutput/DATA2016/histo/AngularCorrelationRun2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output.root";//ppMB, 
    pathin[0] = "FinalOutput/DATA2016/histo/AngularCorrelation161718Full_AOD234_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff.root";//ppMB 
    if (TypeAnalysis==0)  {
      //PRELIMINARY      pathin[0] = "OOJComparisonRun2DataRed_MECorr_hXi_Xi_Eta0.8_sys0_PtTrigMin3.0_PtTrigMin0.2_Output.root";
      pathin[0] = "OOJComparison161718Full_AOD234_hXi_Xi_Eta0.8_sys0_PtTrigMin3.0_PtTrigMin0.2_Output_IsEtaEff.root";
      //      pathin[1] = "OOJComparison17pq_pp5TeV_hXi_pttrig0.15_Xi_Eta0.8_sys0_PtTrigMin3.0_PtTrigMin0.2_Output_IsEtaEff_MultBinning3.root"; //pp5TeV
    }

    if (isComp==0) {
      pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelation161718_HM_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff_MultBinning1.root"; //pp HM
    }
    else {
      pathin[1] = "FinalOutput/DATA2016/histo/AngularCorrelation17pq_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin3.0_Output_IsEtaEff_MultBinning3.root"; //pp5TeV
      if (TypeAnalysis==0)  pathin[1] = "OOJComparison17pq_pp5TeV_hXi_pttrig0.15_Xi_Eta0.8_sys0_PtTrigMin3.0_PtTrigMin0.2_Output_IsEtaEff_MultBinning3.root"; //pp5TeV
    }
  }

  TString pathinSpectra[3] = {"", "", ""}; // pp MB - pp HM 
  //ppMB old AODs  pathinSpectra[0] = "FinalOutput/DATA2016/SystematicAnalysis1617_hK0s_PtBinning1_K0s_Eta0.8_";
  //  pathinSpectra[0] = "FinalOutput/DATA2016/SystematicAnalysis1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_";
  //  pathinSpectra[0] += Region[TypeAnalysis];
  //  pathinSpectra[0] += "Data_PtMin3.0_IsEtaEff.root"; //pp MB 

  pathinSpectra[0] = "FinalOutput/DATA2016/PtSpectraBis_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_PtMin3.0_";
  pathinSpectra[0] += RegionBis[TypeAnalysis];
  pathinSpectra[0] += "_isNormCorrFullyComputed.root"; //pp MB 
  
  if (isComp==0){
    //  pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysis2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_";
    pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysisAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_";
    pathinSpectra[1] += Region[TypeAnalysis];
    pathinSpectra[1] += "Data_PtMin3.0_IsEtaEff_MultBinning1.root";
  }
  else{
    pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysis17pq_hK0s_PtBinning1_K0s_Eta0.8_";
    pathinSpectra[1] += Region[TypeAnalysis];
    pathinSpectra[1] += "Data_PtMin3.0_IsEtaEff";
    if (isComp==2)     pathinSpectra[1]+="_MultBinning3";
    if (isNormFactorCorr){
      pathinSpectra[1]+="_isNormCorr";
      if (TypeAnalysis==0 || TypeAnalysis==1)       pathinSpectra[1]+="_isNFFrom13TeV";
    }
    pathinSpectra[1]+= ".root";
  }
  if (type==1){
    /* OLD -- PRELIMINARY
    pathinSpectra[0] = "FinalOutput/DATA2016/SystematicAnalysisRun2DataRed_MECorr_hXi_Jet0.75_";
    if (TypeAnalysis==0)   pathinSpectra[0] += "OOJNoTriggerSmoothed_";
    pathinSpectra[0] += "Xi_Eta0.8_";
    pathinSpectra[0] += Region[TypeAnalysis];
    pathinSpectra[0] += "Data_PtMin3.0.root"; //pp MB 
    */
    /*
    */
    pathinSpectra[0] = "FinalOutput/DATA2016/SystematicAnalysis161718Full_AOD234_hXi_";
    if (TypeAnalysis==0)     pathinSpectra[0]+= "OOJNoTriggerSmoothed_";
    pathinSpectra[0] += "Xi_Eta0.8_";
    pathinSpectra[0] += Region[TypeAnalysis];
    pathinSpectra[0] += "Data_PtMin3.0_IsEtaEff";
    if (isNormFactorCorr)  pathinSpectra[0] += "_isNormCorr"; //pp MB                                                                       
    pathinSpectra[0] += ".root"; //pp MB                                                                       
    if (isComp==0){
      pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysis161718_HM_hXi";
      if (TypeAnalysis==0)  pathinSpectra[1] += "_OOJAllMult";
      pathinSpectra[1] += "_Xi_Eta0.8_";
      pathinSpectra[1] += Region[TypeAnalysis];
      pathinSpectra[1] += "Data_PtMin3.0_IsEtaEff_MultBinning1.root";
    }
    else {
    pathinSpectra[1] = "FinalOutput/DATA2016/SystematicAnalysis17pq_hXi_";
    if (TypeAnalysis==0)   pathinSpectra[1] += "OOJNoTriggerSmoothed_";
    pathinSpectra[1] += "Xi_Eta0.8_";
    pathinSpectra[1] += Region[TypeAnalysis];
    pathinSpectra[1] += "Data_PtMin3.0_IsEtaEff_MultBinning3";
    if (isNormFactorCorr)  pathinSpectra[1]+="_isNormCorr_isNFFrom13TeV";
    pathinSpectra[1]+=".root";
    }
  }

  TString pathinEff[3] = {"", "", ""}; // pp MB - pp HM 
  //  pathinEff[0] = "FinalOutput/DATA2016/Efficiency/Efficiency1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root"; //ppMB old AODs
  pathinEff[0] = "FinalOutput/DATA2016/Efficiency/Efficiency1617_GP_AOD235_With18c12b_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root"; //ppMB new AODs
  //pathinEff[0] = "FinalOutput/DATA2016/Efficiency/Efficiency2018g4_extra_EtaEff_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
  if (isComp==0){
  //  pathinEff[1] = "FinalOutput/DATA2016/Efficiency/Efficiency2019h11abc_extra_HM_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root"; //ppHM
    pathinEff[1] = "FinalOutput/DATA2016/Efficiency/Efficiency2019h11_HM_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0_MultBinning1.root"; //ppHM
  }
  else {
    pathinEff[1] = "FinalOutput/DATA2016/Efficiency/Efficiency17pq_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0";
    if (isComp==2) pathinEff[1] += "_MultBinning3";
    pathinEff[1] += ".root";
  }
  pathinEff[2] = "FinalOutput/DATA2016/Efficiency/Efficiency1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0_HM.root"; //ppHM 0-1% of the MB

  if (type==1){
    //old    pathinEff[0] = "FinalOutput/DATA2016/Efficiency/Efficiency161718_MD_EtaEff_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root"; //ppMB
    pathinEff[0] = "FinalOutput/DATA2016/Efficiency/Efficiency161718Full_AOD235_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root"; //ppMB
    if (isComp==0){
      pathinEff[1] = "FinalOutput/DATA2016/Efficiency/Efficiency161718_HM_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0_MultBinning1.root"; //pp HM
    }
    else if (isComp==2){
      pathinEff[1] = "FinalOutput/DATA2016/Efficiency/Efficiency17pq_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0_MultBinning3.root"; //pp5 TeV
    }
  }

  cout << "\nstarting from files for MB: " << endl;
  cout << pathin[0] << "\n" << pathinSpectra[0] << "\n" << pathinEff[0] << endl;
  cout << "\nand from files for HM (or pp 5 TeV): " << endl;
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

  gStyle->SetOptStat(0);
  //  gStyle->SetOptFit(1111);

  TCanvas * canvas= new TCanvas("canvas", "canvas", 1300, 800);
  canvas->Divide(5,2);

  Int_t EasyCompCounter[numPtV0][3]= {0};
  Int_t EasyCompCounterSpectrum[3]= {0};
  TCanvas * canvasEasyComp[3];
  TCanvas * canvasEasyCompRatio[3];
  TCanvas * canvasEasyCompSpectrum[3];
  for (Int_t t=0; t<3; t++){
    canvasEasyComp[t]= new TCanvas(Form("canvasEasyComp%i", t), Form("canvasEasyComp%i", t), 1300, 800);
    canvasEasyComp[t]->Divide(5,2);
    canvasEasyCompRatio[t]= new TCanvas(Form("canvasEasyCompRatio%i", t), Form("canvasEasyCompRatio%i", t), 1300, 800);
    canvasEasyCompRatio[t]->Divide(5,2);
    canvasEasyCompSpectrum[t]= new TCanvas(Form("canvasEasyCompSpectrum%i", t), Form("canvasEasyCompSpectrum%i", t), 1300, 800);
    canvasEasyCompSpectrum[t]->Divide(2,1);
  }

  TCanvas * canvasRatio[3];
  for (Int_t t=0; t<3; t++){
    canvasRatio[t] = new TCanvas(Form("canvasRatio%i",t), Form("canvasRatio%i",t), 1300, 800);
    canvasRatio[t]->Divide(5,2);
  }

  TCanvas * canvasSpectrum = new TCanvas("canvasSpectrum", "canvasSpectrum", 1300, 800);
  canvasSpectrum->Divide(2,1);

  TCanvas * canvasEff= new TCanvas("canvasEff", "canvasEff", 1300, 800);
  canvasEff->Divide(2,1);

  TCanvas * canvasEffEta= new TCanvas("canvasEffEta", "canvasEffEta", 1300, 800);
  canvasEffEta->Divide(2,1);

  TCanvas * canvasEff2D= new TCanvas("canvasEff2D", "canvasEff2D", 1300, 800);
  canvasEff2D -> Divide(5,2);
  TCanvas * canvasEff2DRatio= new TCanvas("canvasEff2DRatio", "canvasEff2DRatio", 1300, 800);
  canvasEff2DRatio -> Divide(5,2);

  TLegend * legendMult = new TLegend(0.35, 0.5, 0.9, 0.9);

  TH1F * histoPhiDistr[NSystems][nummolt+1][numPtV0];
  TH1F * histoPhiDistrClone[NSystems][nummolt+1][numPtV0];
  TH1F * histoPhiDistrCloneRatio[numPtV0][3];
  TH1F * histoPhiDistrRatio[NSystems][nummolt+1][numPtV0];
  TH1F * histoPhiDistrDenom[NSystems];
  TH1F * histoSpectrum[NSystems][nummolt+1];
  TH1F * histoSpectrumMB;
  TH1F * histoSpectrumFinal[NSystems][nummolt+1];
  TH1F * histoSpectrumMBFinal;
  TH1F * histoSpectrumRatio[NSystems][nummolt+1];
  TH1F * histoSpectrumRatioBis[3];
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
  TString  SmoltMB[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString  SmoltHM[nummolt+1]={"0-0.001", "0.001-0.005", "0.005-0.01", "0.01-0.05", "0.05-0.1", "_all"};
  TString  SmoltHMMultBinning1[nummolt+1]={"0-0a", "0-0b", "0-0.01", "0.01-0.05", "0.05-0.1", "_all"};
  TString  Smolt5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};

  TString  Smolt[nummolt+1];
  Double_t dNdetaMB13TeV[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};
  Double_t dNdetaHM13TeV[nummolt+1]={39.40, 36.89, 35.16, 32.57, 30.43, 36.29};
  Double_t dNdetaMB5TeV[nummolt+1] ={15.27, 11.91, 8.73, 5.78, 3.03, 5.49};
  Double_t dNdetaMB5TeVBis[nummolt+1] ={13.595, 4.91,0, 0, 0, 5.49};
  Double_t dNdeta[nummolt+1];
  
  Double_t NmoltMB[nummolt+1]={0,5,10,30,50,100};
  Double_t NmoltHM[nummolt+1]={0,0.001,0.005,0.01,0.05,0.1};
  Double_t Nmolt5TeV[nummolt+1]={0,10, 100, 100, 100, 100};
  Double_t Nmolt[nummolt+1];

  TLine *lineAt0 = new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2, 0);
  if (TypeAnalysis==0) lineAt0 = new TLine(-1.1, 0, 1.1, 0);
  lineAt0->SetLineColor(1);

  TLine *lineAt1 = new TLine(-TMath::Pi()/2, 1, 3*TMath::Pi()/2, 1);
  if (TypeAnalysis==0) lineAt1 = new TLine(-1.1, 1, 1.1, 1);
  lineAt1->SetLineColor(1);

  TLine *lineAt1Spectrum = new TLine(0, 1, 8, 1);
  lineAt1Spectrum->SetLineColor(1);

  for (Int_t i=0; i< NSystems; i++){
    file[i] = new TFile(pathin[i], "");
    if (!file[i]) {cout << pathin[i] << " does not exist " << endl; return;}
    for (Int_t m =0; m<nummolt+1; m++){
      if (isComp==0 && m<=1 && i==1) continue;
      if (i==0) {
	Nmolt[m] = NmoltMB[m];
	Smolt[m] = SmoltMB[m];
	dNdeta[m] = dNdetaMB13TeV[m];
	Color[m] = Color0[m];
      }
      else if (i==1) {
	Color[m] = Color1[m];
	if (isComp==0) {
	  Nmolt[m] = NmoltHM[m];
	  Smolt[m] = SmoltHM[m];
	  if (MultBinning==1) 	  Smolt[m] = SmoltHMMultBinning1[m];
	  dNdeta[m] = dNdetaHM13TeV[m];
	}
	else if (isComp==1) {
	  Nmolt[m] = NmoltMB[m];
	  Smolt[m] = SmoltMB[m];
	  dNdeta[m] = dNdetaMB5TeV[m];
	}
	else if (isComp==2) {
	  Nmolt[m] = Nmolt5TeV[m];
	  Smolt[m] = Smolt5TeV[m];
	  dNdeta[m] = dNdetaMB5TeVBis[m];
	}
      }
    }

    cout << "\nPhi distributions comparison for i=" << i<< endl;
    for (Int_t m=0; m< nummolt; m++){
      if (isComp==0 && m<=1 && i==1) continue;
      if (i==1 && (m==2 || m==3 || m==4) && isComp==2) continue; 
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      //      if (i==1 && m<3 && TypeAnalysis==1) continue;
      cout << "mult " << m << endl;
      for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	NameHisto[i][m][v] = "ME_m"+ Smolt[m]+ "_v"+ SPtV0[v] + "_AC_phi_V0Sub_"+RegionName[TypeAnalysis] +"EffCorr_TrCorr";
	if (i==0) NameHisto[i][mRatio][v] = "ME_m"+ Smolt[mRatio]+ "_v"+ SPtV0[v] + "_AC_phi_V0Sub_"+RegionName[TypeAnalysis] +"EffCorr_TrCorr";
	if (TypeAnalysis==0 && type==1 && isComp!=0){
	  if ( (NPtV0[v] >=0.1 && NPtV0[v] <=2.0)) {
	    NameHisto[i][m][v] = "ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmooth";
	    if (i==0) 	    NameHisto[i][mRatio][v] = "ME_m"+Smolt[mRatio]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmooth";
	  }
	  else {
	    NameHisto[i][m][v]= "ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_DefaultMethod";
	    if (i==0)  NameHisto[i][mRatio][v]= "ME_m"+Smolt[mRatio]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_DefaultMethod";
	  }
	}
	histoPhiDistr[i][m][v] = (TH1F*) file[i]->Get(NameHisto[i][m][v]);
	if (!histoPhiDistr[i][m][v] ) {cout << "Phi distr file not found: " << NameHisto[i][m][v] << endl; return;}
	//histoPhiDistr[i][m][v]->Sumw2();

	Float_t DeltaPt = NPtV0[v+1] - NPtV0[v];
	histoPhiDistr[i][m][v] ->Scale(1./DeltaPt);

	if (TypeAnalysis==0)	histoPhiDistr[i][m][v]->GetXaxis()->SetRangeUser(-1.1, 1.1);
	StyleHisto(histoPhiDistr[i][m][v], Low[TypeAnalysis] , Up[TypeAnalysis], Color[m], 33, "#Delta#varphi", "dN/dp_{T} 1/N_{Trigg} 1/#Delta#eta", Form("%.1f < p_{T} < %.1f GeV/c", NPtV0[v], NPtV0[v+1]), 0,0,0);

	histoPhiDistrClone[i][m][v] = (TH1F*) histoPhiDistr[i][m][v]->Clone(NameHisto[i][m][v]+"_Clone");
	if (TypeAnalysis!=0)	histoPhiDistrClone[i][m][v]->Rebin(2);
	StyleHisto(histoPhiDistrClone[i][m][v], 0, 1.5*(histoPhiDistrClone[i][m][v]->GetBinContent(histoPhiDistrClone[i][m][v]->GetMaximumBin())), Color[m], 33, "#Delta#varphi", "dN/dp_{T} 1/N_{Trigg} 1/#Delta#eta", Form("%.1f < p_{T} < %.1f GeV/c", NPtV0[v], NPtV0[v+1]), 0,0,0);
	if (TypeAnalysis==0) {
	  histoPhiDistrClone[i][m][v]->GetXaxis()->SetRangeUser(-1.1, 1.1);
	  histoPhiDistrClone[i][m][v]->GetYaxis()->SetRangeUser(Low[TypeAnalysis], 1.5*(histoPhiDistrClone[i][m][v]->GetBinContent(histoPhiDistrClone[i][m][v]->GetMaximumBin())));
	}

	Bool_t isPassed =0;
	for (Int_t t =0; t<3; t++){
	  canvasEasyComp[t]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  if (t==0) isPassed = ( (i==0 && (m==4)) || (i==1 && (m==4)) );
	  else if (t==1) isPassed =  ( (i==0 && (m==1)) || (i==1 && (m==0)) );
	  else if (t==2) isPassed =  ( (i==0 && (m==2)) || (i==1 && (m==1)) );
	  if (isComp==2 && t==0) isPassed = ( (i==0 && (m==3)) || (i==1 && (m==1)) );
	  if (isPassed){
	    //	    cout << " i " << i << " m " << m << endl;
	    histoPhiDistrClone[i][m][v]->DrawClone("same");
	    if (TypeAnalysis==0) lineAt0->Draw("same");

	    EasyCompCounter[v][t]++;
	    if (EasyCompCounter[v][t]==1)  histoPhiDistrCloneRatio[v][t] = (TH1F*) histoPhiDistrClone[i][m][v]->Clone(NameHisto[i][m][v]+Form("_CloneRatio%i", t));
	    else  histoPhiDistrCloneRatio[v][t] ->Divide( histoPhiDistrClone[i][m][v]);
	  }

	  canvasEasyCompRatio[t]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  if (EasyCompCounter[v][t]==2)	{
	    histoPhiDistrCloneRatio[v][t]->GetYaxis()->SetRangeUser(0.5, 1.5);
	    histoPhiDistrCloneRatio[v][t]->DrawClone("same");
	    lineAt1->Draw("same");
	  }
	}

	Float_t NormFactor = histoPhiDistr[i][m][v]->Integral();
	if (TypeAnalysis!=0)  NormFactor = histoPhiDistr[i][m][v]->Integral(1,2);

	if (TypeAnalysis!=0) 	histoPhiDistr[i][m][v]->Rebin(2);
	NormFactor = histoPhiDistr[i][m][v]->Integral(1,2);
	if (TypeAnalysis!=0)	histoPhiDistr[i][m][v]->Scale(1./NormFactor);

	if (i==0 && m==0){
	  histoPhiDistr[i][m][v]->SetName("dummy");
	  histoPhiDistrDenom[v] = (TH1F*) file[i]->Get(NameHisto[i][mRatio][v]);
	  histoPhiDistrDenom[v]->SetName("PhiDistrDenom"+ SPtV0[v]);
	  histoPhiDistr[i][m][v]->SetName(NameHisto[i][m][v]);
	  histoPhiDistrDenom[v] ->Scale(1./DeltaPt);
	  if (TypeAnalysis!=0)	  histoPhiDistrDenom[v]->Rebin(2);
	  Float_t NormFactorDenom = histoPhiDistrDenom[v]->Integral();
	  if (TypeAnalysis!=0) NormFactorDenom = histoPhiDistrDenom[v]->Integral(1,2);
	  if (TypeAnalysis!=0) histoPhiDistrDenom[v]->Scale(1./NormFactorDenom);
	}

	canvas->cd(v+1);
	gPad->SetLeftMargin(0.15);
	histoPhiDistr[i][m][v]->DrawClone("same");
	if (i==0 && m==0 && TypeAnalysis==0)	lineAt0->Draw("same");

	histoPhiDistrRatio[i][m][v] = (TH1F*) histoPhiDistr[i][m][v] -> Clone(NameHisto[i][m][v] + "Ratio");
	//	if (v==PtV0Min)	cout << "i " << i << " m " << m << " binning ratio: " << histoPhiDistrRatio[i][m][v]->GetNbinsX() << " denom: " << histoPhiDistrDenom[v]->GetNbinsX() << endl;
	histoPhiDistrRatio[i][m][v]->Divide(histoPhiDistrDenom[v]);
	//if (v==PtV0Min)	cout << "i " << i << " m " << m << " binning ratio: " << histoPhiDistrRatio[i][m][v]->GetNbinsX() << " denom: " << histoPhiDistrDenom[v]->GetNbinsX() << endl;

	canvasRatio[0]->cd(v+1);
	gPad->SetLeftMargin(0.15);
	histoPhiDistrRatio[i][m][v]->GetYaxis()->SetRangeUser(LowRatio[TypeAnalysis], UpRatio[TypeAnalysis]);
	if (!(i==0 && m==mRatio))	histoPhiDistrRatio[i][m][v]->DrawClone("same");
	if (i==1 && m==0)	lineAt1->Draw("same");

	if (TypeAnalysis==1){
	  canvasRatio[1]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  histoPhiDistrRatio[i][m][v]->GetXaxis()->SetRangeUser(-1.1, 1.1);
	  if (!(i==0 && m==mRatio))	histoPhiDistrRatio[i][m][v]->DrawClone("same");
	  if (i==1 && m==0)	lineAt1->Draw("same");

	  canvasRatio[2]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  histoPhiDistrRatio[i][m][v]->GetXaxis()->SetRangeUser(TMath::Pi()/2, 3./2*TMath::Pi());
	  if (!(i==0 && m==mRatio))	histoPhiDistrRatio[i][m][v]->DrawClone("same");
	  if (i==1 && m==0)	lineAt1->Draw("same");

	}
      }
    }
  }

  //legend
  
  for (Int_t i=1; i>=0; i--){
    for (Int_t m=0; m< nummolt; m++){
      if (isComp==0 && m<=1 && i==1) continue;
      if (i==1 && (m==2 || m==3 || m==4) && isComp==2) continue; 
      if (i==0) {
	Nmolt[m] = NmoltMB[m];
	Smolt[m] = SmoltMB[m];
	dNdeta[m] = dNdetaMB13TeV[m];
	Color[m] = Color0[m];
      }
      else if (i==1) {
	Color[m] = Color1[m];
	if (isComp==0) {
	  Nmolt[m] = NmoltHM[m];
	  Smolt[m] = SmoltHM[m];
	  if (MultBinning==1) 	  Smolt[m] = SmoltHMMultBinning1[m];
	  dNdeta[m] = dNdetaHM13TeV[m];
	}
	else if (isComp==1) {
	  Nmolt[m] = NmoltMB[m];
	  Smolt[m] = SmoltMB[m];
	  dNdeta[m] = dNdetaMB5TeV[m];
	}
	else if (isComp==2) {
	  Nmolt[m] = Nmolt5TeV[m];
	  Smolt[m] = Smolt5TeV[m];
	  dNdeta[m] = dNdetaMB5TeVBis[m];
	}
      }
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      //      if (i==1 && m<3 && TypeAnalysis==1) continue;
      if (isComp==0) legendMult->AddEntry(histoPhiDistr[i][m][PtV0Min], Smolt[m] , "pl");
      else if (isComp==1) legendMult->AddEntry(histoPhiDistr[i][m][PtV0Min], SSystem[i] + " "+ Smolt[m]+ "%"+Form(" dN/d#eta=%.2f", dNdeta[m]), "pl");
      else if (isComp==2) legendMult->AddEntry(histoPhiDistr[i][m][PtV0Min], SSystem[i] + " "+ Smolt[m]+ "%"+Form(" dN/d#eta=%.2f", dNdeta[m]), "pl");
    }
  }
  canvas->cd(10);
  legendMult->Draw();

  for (Int_t t = 0; t < 3; t++){
    canvasEasyComp[t]->cd(10);
    legendMult->Draw();
  }
  
  //comparison of spectra
  cout << "\nSpectra distributions comparison "<<endl;
  for (Int_t i=0; i< NSystems; i++){
    fileSpectra[i] = new TFile(pathinSpectra[i], "");
    if (!fileSpectra[i]) {cout << pathinSpectra[i] << " does not exist " << endl; return; }
    for (Int_t m=0; m<=nummolt; m++){
      if (isComp==0 && m<=1 && i==1) continue;
      if (i==1 && (m==2 || m==3 || m==4) && isComp==2) continue; 
      //      cout << "m " << m << endl;
      if (i==0) {
	Nmolt[m] = NmoltMB[m];
	Smolt[m] = SmoltMB[m];
	dNdeta[m] = dNdetaMB13TeV[m];
	Color[m] = Color0[m];
      }
      else if (i==1) {
	Color[m] = Color1[m];
	if (isComp==0) {
	  Nmolt[m] = NmoltHM[m];
	  Smolt[m] = SmoltHM[m];
	  if (MultBinning==1) 	  Smolt[m] = SmoltHMMultBinning1[m];
	  dNdeta[m] = dNdetaHM13TeV[m];
	}
	else if (isComp==1) {
	  Nmolt[m] = NmoltMB[m];
	  Smolt[m] = SmoltMB[m];
	  dNdeta[m] = dNdetaMB5TeV[m];
	}
	else if (isComp==2) {
	  Nmolt[m] = Nmolt5TeV[m];
	  Smolt[m] = Smolt5TeV[m];
	  dNdeta[m] = dNdetaMB5TeVBis[m];
	}
      }
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      //      if (i==1 && m<3 && TypeAnalysis==1) continue;
      if (i==1 && m==nummolt) continue; //I am interested in 0-100% only
      NameHistoSpectrum[i][m] = Form("fHistSpectrumPart_m%i_syst0", m);
      //      cout << pathinSpectra[i]<< endl;
      if (pathinSpectra[i].Index("PtSpectraBis")!=-1)  { NameHistoSpectrum[i][m] = "fHistSpectrum_" + Smolt[m];}
      NameHistoSpectrumFinal[i][m] = "fHistSpectrumPart_m_" + Smolt[m] + "syst0";
      //      cout << NameHistoSpectrum[i][m] << endl;
      histoSpectrum[i][m] = (TH1F*) fileSpectra[i]->Get(NameHistoSpectrum[i][m]);
      if (!histoSpectrum[i][m]) {cout << NameHistoSpectrum[i][m] <<" not present " <<  endl; return;}
      histoSpectrumFinal[i][m] = (TH1F*)       histoSpectrum[i][m]->Clone(NameHistoSpectrumFinal[i][m]);
      //      histoSpectrumFinal[i][m] ->Sumw2();

      if (m==0 && i==0){
	NameHistoSpectrum[i][mRatio] = Form("fHistSpectrumPart_m%i_syst0", mRatio);
	if (pathinSpectra[i].Index("PtSpectraBis")!=-1)  { NameHistoSpectrum[i][mRatio] = "fHistSpectrum_" + Smolt[mRatio];}
	histoSpectrumMB = (TH1F*) fileSpectra[i]->Get(NameHistoSpectrum[i][mRatio]);
	if (!histoSpectrumMB) {cout << NameHistoSpectrum[i][mRatio] << endl; return;}
	histoSpectrumMBFinal = (TH1F*)       histoSpectrumMB->Clone("fhistoSpectrumMB");
	//	histoSpectrumMBFinal->Sumw2();
      }

      histoSpectrumRatio[i][m] = (TH1F*)       histoSpectrumFinal[i][m]->Clone(NameHistoSpectrumFinal[i][m]+ "_Ratio");
      histoSpectrumRatio[i][m] ->Divide(histoSpectrumMBFinal);

      StyleHisto(histoSpectrumFinal[i][m],LowSpectrum[TypeAnalysis], UpSpectrum[TypeAnalysis], Color[m], 33, "p_{T}", "dN/dp_{T} 1/N_{Trigg} 1/#Delta#eta #Delta#varphi", "", 0,0,0);

      canvasSpectrum->cd(1);
      gPad->SetLeftMargin(0.15);
      histoSpectrumFinal[i][m]->SetStats(0);
      histoSpectrumFinal[i][m]->Draw("same");
      legendMult->Draw();

      StyleHisto(histoSpectrumRatio[i][m],LowSpectrumRatio[TypeAnalysis], UpSpectrumRatio[TypeAnalysis], Color[m], 33, "p_{T}", "Ratio to MB " + Smolt[mRatio] + "%", "", 0,0,0);

      canvasSpectrum->cd(2);
      gPad->SetLeftMargin(0.15);
      if (!(m==mRatio && i==0)) {
	if (m!=nummolt || TypeAnalysis!=0) histoSpectrumRatio[i][m]->Draw("same");
	histoSpectrumRatio[i][m]->SetStats(0);
	lineAt1Spectrum->Draw("same");
      }

      //---------comparison between specific mult classes
      cout << "\nComparison in specific multiplicity classes " << endl;
      Bool_t isPassed =0;
      for (Int_t t =0; t<3; t++){
	canvasEasyCompSpectrum[t]->cd(1);
	gPad->SetLeftMargin(0.15);
	if (t==0) isPassed = ( (i==0 && (m==4)) || (i==1 && (m==4)) );
	else if (t==1) isPassed =  ( (i==0 && (m==1)) || (i==1 && (m==0)) );
	else if (t==2) isPassed =  ( (i==0 && (m==2)) || (i==1 && (m==1)) );
	if (isComp==2){
	  if (t==0) isPassed = ( (i==0 && (m==3)) || (i==1 && (m==1)) );
	  else if (t==1) isPassed = ( (i==0 && (m==1)) || (i==1 && (m==0)) );
	  else if (t==2) isPassed = ( (i==0 && (m==2)) || (i==1 && (m==0)) );
	}
	if (isPassed){
	  //	  cout <<"->Comparison between (i="<< i<<" && m=" << m << ")" << endl;
	  histoSpectrumFinal[i][m]->DrawClone("same");
	  EasyCompCounterSpectrum[t]++;
	  if (EasyCompCounterSpectrum[t]==1)  histoSpectrumRatioBis[t] = (TH1F*) histoSpectrumFinal[i][m]->Clone(Form("fHistSpectrumPart_t%i_RatioBis", t));
	  else  histoSpectrumRatioBis[t] ->Divide( histoSpectrumFinal[i][m]);
	}

	canvasEasyCompSpectrum[t]->cd(2);
	gPad->SetLeftMargin(0.15);
	if (EasyCompCounterSpectrum[t]==2)	{
	  histoSpectrumRatioBis[t]->GetYaxis()->SetRangeUser(0.5, 1.5);
	  histoSpectrumRatioBis[t]->DrawClone("same");
	  lineAt1Spectrum->Draw("same");
	}
      }
      //---------comparison between specific mult classes
    }
  }

  //comparison of efficiencies
  cout << "\nEfficiency comparison" << endl;
  for (Int_t i=0; i< NSystems; i++){
    fileEff[i] = new TFile(pathinEff[i], "");
    if (!fileEff[i]) {cout << pathinEff[i]<< " does not exist " << endl;}
    for (Int_t m=0; m<=nummolt; m++){
      if (isComp==0 && m<=1 && i==1) continue;
      if (i==1 && (m==2 || m==3 || m==4) && isComp==2) continue; 
      if (i==0) {
	Nmolt[m] = NmoltMB[m];
	Smolt[m] = SmoltMB[m];
	dNdeta[m] = dNdetaMB13TeV[m];
	Color[m] = Color0[m];
      }
      else if (i==1) {
	Color[m] = Color1[m];
	if (isComp==0) {
	  Nmolt[m] = NmoltHM[m];
	  Smolt[m] = SmoltHM[m];
	  if (MultBinning==1) 	  Smolt[m] = SmoltHMMultBinning1[m];
	  dNdeta[m] = dNdetaHM13TeV[m];
	}
	else if (isComp==1) {
	  Nmolt[m] = NmoltMB[m];
	  Smolt[m] = SmoltMB[m];
	  dNdeta[m] = dNdetaMB5TeV[m];
	}
	else if (isComp==2) {
	  Nmolt[m] = Nmolt5TeV[m];
	  Smolt[m] = Smolt5TeV[m];
	  dNdeta[m] = dNdetaMB5TeVBis[m];
	}
      }
      else if (i==2) Color[0] = 922;
      //      if (i==1 && m<3 && TypeAnalysis==0) continue;
      if (i==1 && m==nummolt) continue; //I am interested in 0-100% only
      if (i==2 && m!=0) continue; //I am interested in 0-1% only 

      NameHistoEff[i][m] = "fHistV0EfficiencyPtBins_"+Smolt[m];
      //      if (EffEta)      NameHistoEff[i][m] = "fHistV0EfficiencyEta_"+Smolt[m];
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

      StyleHisto(histoEffRatio[i][m],0.5, 1.5, Color[m], 33, "p_{T}", "Ratio to MB "+ Smolt[mRatio]+"%", "", 0,0,0);
      StyleHisto(histoEffEtaRatio[i][m],0.5, 1.5, Color[m], 33, "#eta", "Ratio to MB " + Smolt[mRatio] +"%", "", 0,0,0);

      canvasEff->cd(2);
      gPad->SetLeftMargin(0.15);
      if (!(m==mRatio && i==0))  {
	histoEffRatio[i][m]->Draw("same");
	lineAt1Spectrum->Draw("same");
      }

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
  if (pathin[0].Index("AOD234")!=-1) NameFileout +="_13TeVNewAODs";
  if (pathin[1].Index("2016k_HM_hK0s")!=-1) NameFileout += "_HM16k";
  if (pathin[1].Index("17pq")!=-1) NameFileout += "_pp5TeV";
  if (pathinEff[0].Index("18g4_extra")!=-1) NameFileout += "_MBEff18g4extra";
  if (pathinEff[1].Index("2019_HM")!=-1) NameFileout += "_EffNewAOD";
  NameFileout += "_"+tipo[type];
  NameFileout += "_" +Region[TypeAnalysis];
  NameFileout += Form("_RatioToMBMult%i",mRatio);

  TFile * fileout = new TFile (NameFileout+ ".root", "RECREATE");
  canvas->Write();
  for (Int_t t=0; t<3; t++){
    canvasEasyComp[t]->Write();
    canvasEasyCompRatio[t]->Write();
    canvasEasyCompSpectrum[t]->Write();
  }
  canvasRatio[0]->Write();
  if (TypeAnalysis==1){
  canvasRatio[1]->Write();
  canvasRatio[2]->Write();
  }
  canvasSpectrum->Write();
  canvasEff->Write();
  canvasEffEta->Write();
  canvasEff2D->Write();
  canvasEff2DRatio->Write();
  fileout->Close();
  canvas->SaveAs(NameFileout +".pdf");
  canvasSpectrum->SaveAs(NameFileout +"_Spactra.pdf");

  cout <<"\n\e[35mInput files: " << endl;
  cout << "Efficiencies:\e[39m\n" << pathinEff[0] <<"\n" << pathinEff[1] << endl;
  cout << "\e[35mPhi distributions:\e[39m\n" << pathin[0] <<"\n" << pathin[1] << endl;
  cout << "\e[35mSpectra:\e[39m\n" << pathinSpectra[0] <<"\n" << pathinSpectra[1] << endl;
  cout << "\nRatios performed with respect to the mult class " << mRatio << " in the "<< SSystem[0] << endl;
  cout << "\nI have produced the file " << NameFileout << ".root and .pdf " << endl;
}
