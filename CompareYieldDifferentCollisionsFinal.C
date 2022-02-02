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

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, Int_t msize, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(msize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.2);
  histo->SetTitle(title);
}

void CompareYieldDifferentCollisionsFinal(Int_t TypeAnalysis=0,  Int_t type=0, Int_t isComp=0,  Bool_t isAvgPtvsMult=0,  Bool_t isPreliminary=0, Bool_t FitYields=0){

  if (isPreliminary && isAvgPtvsMult) return;
  const   Int_t NSystems =2;

  if (isAvgPtvsMult) FitYields = 0;
  //isComp =0 : MB compared to HM
  //isComp =1 : MB compared to pp5TeV (same mult classes)
  //isComp =2 : MB compared to pp5TeV (only two mult classes for 5TeV)

  Float_t ScalingFactor =1;
  //  if (isPreliminary && TypeAnalysis==0) ScalingFactor = 0.8458/1.08747; //=1./1.286
  if (TypeAnalysis==0 && isAvgPtvsMult==0)  ScalingFactor = 0.8458/1.08747; //=1./1.286
  //ScalingFactor =1./( 0.8458/1.08747); to adjust MB preliminary bins to the others
  TString SisPreliminary[2] = {"", "_isPreliminary"};
  TString YieldOrAvgPt[2] = {"Yield", "AvgPt"};
  TString SFitYields[2] = {"", "_FitToYields"};
  TString CollisionsComp[3] = {"_vsHM", "_vs5TeV5Mult", "_vs5TeV"};
  gStyle->SetOptStat(0);
  TString titleYieldX="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
  TString titleYieldY="N/N_{trigg} 1/#Delta#eta #Delta#it{#varphi}";
  if (isAvgPtvsMult) titleYieldY="#LT#it{p}_{T}#GT";
  TString tipo[2] = {"K0s", "Xi"};
  TString Region[3] = {"Jet", "Bulk", "All"};
  TString RegionBis[3] = {"Jet", "Bulk", "Inclusive"};
  Int_t UpRangeMult = 45;
  Int_t NbinsMult = 225; //450
  if (isComp!=0)  {
    UpRangeMult = 25;
    //NbinsMult = 125;
    NbinsMult = 250;
  }
  Float_t Up[3] = {0.05, 0.4, 0.4}; //0.035
  Float_t Low[3] = {0.015, 10e-4, 10e-4};
  Float_t LowPt[3] = {1.5,0.8,0.8};
  Float_t UpPt[3] = {2.5, 1.3, 1.3};
  if (type==1) {
    Low[0] = 10e-6; Up[0] = 0.003;
    Up[1] = 0.03;
    Up[2] = 0.03;
    LowPt[0] = 1.5;     LowPt[1] = 1;     LowPt[2] = 1;
    UpPt[0] = 5;     UpPt[1] = 2.2;     UpPt[2] = 2.2;
  }
  if (isComp!=0){
    if (type==0) {
      Up[1] = 0.25;
      Up[2] = 0.25;
    }
    if (type==1) {
      Low[0] = 10e-6; Up[0] = 0.003;
      Up[1] = 0.015;
      Up[2] = 0.015;
    }
  }
  if (isAvgPtvsMult){
    for (Int_t i =0; i<3; i++){
      Low[i] = LowPt[i];
      Up[i] = UpPt[i];
    }
  }
  Int_t Color[3] = {628,418,601};
  Int_t ColorDiff[3][2] = {{634,628}, {418,418} , {601, 601}};
  /*
  if (isComp!=0){
    ColorDiff[0][1] = 628;
    ColorDiff[1][1] = 419;
    ColorDiff[2][1] = 603;
  }
  */
  //  if (isComp!=0){
    ColorDiff[0][1] = 628;
    ColorDiff[1][1] = 829;
    ColorDiff[2][1] = 867;
    //  }
  Int_t Marker[2] = {33, 21};
  Float_t Size[2] = {2, 1.4};
  TFile *file[NSystems];
  TString pathin[NSystems] = {"", ""}; // pp MB - pp HM 
  TString ParticleType[2]={"K0s", "Xi"};
  if (type==0){
    pathin[0] = "FinalOutput/DATA2016/PtSpectraBisNew_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_PtMin3.0_";
    pathin[0] += RegionBis[TypeAnalysis];
    pathin[0] +="_isNormCorrFullyComputed_YieldMeanMacro_isErrorAssumedPtCorr.root";
    if (isPreliminary){
      pathin[0] = "FinalOutput/DATA2016/PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_";
      pathin[0] += RegionBis[TypeAnalysis];
      pathin[0] += "_isNormCorr_Preliminaries.root";
    }
    if (isComp==0){
      pathin[1] = "FinalOutput/DATA2016/PtSpectraBisNew_pp13TeVHM_PtBinning1_AllhK0sHM_RedNo16k_K0s_Eta0.8_PtMin3.0_";
      pathin[1] += RegionBis[TypeAnalysis];
      pathin[1] +="_isNormCorrFullyComputed_YieldMeanMacro_isErrorAssumedPtCorr_MultBinning1.root";
    }
    else {
      pathin[1] = "FinalOutput/DATA2016/PtSpectraBisNew_pp5TeV_PtBinning1_17pq_hK0s_K0s_Eta0.8_PtMin3.0_";
      pathin[1] += RegionBis[TypeAnalysis];
      pathin[1] +="_isNormCorrFullyComputed_YieldMeanMacro_isErrorAssumedPtCorr";
      if (isComp==2)     pathin[1]+="_MultBinning3";
      pathin[1]+= ".root";
    }
  }
  else if (type==1){
    pathin[0] = "FinalOutput/DATA2016/PtSpectraBisNew_161718Full_AOD234_hXi_Xi_Eta0.8_PtMin3.0_";
    pathin[0] += RegionBis[TypeAnalysis];
    pathin[0] += "_isNormCorrFullyComputed_YieldMeanMacro_isErrorAssumedPtCorr.root";
    if (isPreliminary){
      pathin[0] = "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_";
      pathin[0] += RegionBis[TypeAnalysis];
      pathin[0] += "_isNormCorr_Preliminaries.root";
    }
    if (isComp==0){
      pathin[1] = "FinalOutput/DATA2016/PtSpectraBisNew_pp13TeVHM_161718_HM_hXi_WithFlat16k_No18p_Xi_Eta0.8_PtMin3.0_";
      pathin[1] += RegionBis[TypeAnalysis];
      pathin[1] += "_isNormCorrFullyComputed_YieldMeanMacro_isErrorAssumedPtCorr_MultBinning1.root";
    }
    else {
      pathin[1] = "FinalOutput/DATA2016/PtSpectraBisNew_pp5TeV_17pq_hXi_Xi_Eta0.8_PtMin3.0_";
      pathin[1] += RegionBis[TypeAnalysis];
      pathin[1] += "_isNormCorrFullyComputed_YieldMeanMacro_isErrorAssumedPtCorr_MultBinning3.root";
    }
  }

  TF1 * pol1Common = new TF1("pol1Common", "pol1", 0, UpRangeMult);
  pol1Common->SetLineColor(Color[TypeAnalysis]);
  pol1Common->SetLineWidth(0.3);
  TF1 *pol1[NSystems];
  TH1F *histoYield[NSystems];
  TH1F *histoYieldSist[NSystems];
  TH1F *histoYieldFinal[NSystems];
  TH1F *histoYieldSistFinal[NSystems];
  TString NameHisto[NSystems] = {"fHistYieldStat", "fHistYieldStat"};
  TString NameHistoSist[NSystems] = {"fHistYieldSist", "fHistYieldSist"};
  TString NameHistoFinal[NSystems] = {"fHistYield_ppMB", "fHistYield_ppHM"};
  TString NameHistoFinalSist[NSystems] = {"fHistYield_ppMB_Sist", "fHistYield_ppHM_Sist"};
  if (isComp!=0) {
    NameHistoFinal[1] = "fHistYield_pp5TeV";
    NameHistoFinalSist[1] = "fHistYield_pp5TeV_Sist";
  }
  if (isAvgPtvsMult) {
    NameHisto[0] = "fHistPtvsMultStat";     NameHisto[1] = "fHistPtvsMultStat";
    NameHistoSist[0] = "fHistPtvsMultSist";     NameHistoSist[1] = "fHistPtvsMultSist";
  }
  TH1F * histoYieldComparison = new TH1F ("histoYieldComparison", "histoYieldComparison",NbinsMult, 0, UpRangeMult );
  TH1F * histoYieldSistComparison = new TH1F ("histoYieldSistComparison", "histoYieldSistComparison",NbinsMult, 0, UpRangeMult );
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  TCanvas * canvasBoth = new TCanvas("canvasBoth", "canvasBoth", 1300, 800);

  Int_t bin=0;
  for (Int_t b=1; b<= histoYieldComparison->GetNbinsX(); b++){
    histoYieldComparison -> SetBinContent(b, 0);
    histoYieldComparison -> SetBinError(b, 0);
    histoYieldSistComparison -> SetBinContent(b, 0);
    histoYieldSistComparison -> SetBinError(b, 0);
  }


  for (Int_t i=0; i< NSystems; i++){
    file[i] = new TFile(pathin[i], "");
    histoYield[i] = (TH1F*) file[i]->Get(NameHisto[i]);
    if (!histoYield[i]) { cout << "Histo not available " << endl; return; }
    histoYieldFinal[i]= (TH1F*)    histoYield[i] ->Clone(NameHistoFinal[i]);
    histoYieldSist[i] = (TH1F*) file[i]->Get(NameHistoSist[i]); 
    if (!histoYieldSist[i]) { cout << "Histo sist not available " << endl; return; }
    histoYieldSistFinal[i]= (TH1F*)    histoYieldSist[i] ->Clone(NameHistoFinalSist[i]);

    if (TypeAnalysis==0){
      if (isPreliminary){
	if (i==1){ //5 TeV and 13 TeV HM were done with dPhi = 1.08, whereas preliminary results were done with dphi = 0.85. In this way I adjust the dphi interval of the HM and 5 TeV results TO the 13 TeV preliminarie s(which cannot be changed)
	  histoYieldFinal[i]->Scale(1./ScalingFactor);
	  histoYieldSistFinal[i]->Scale(1./ScalingFactor);
	}
      }
      else {
	histoYieldFinal[i]->Scale(1./ScalingFactor);
	histoYieldSistFinal[i]->Scale(1./ScalingFactor);
      }
    }
    cout << "Bin width of final histos " << histoYieldFinal[i]->GetBinWidth(1) << endl;
  
    canvasBoth->cd();
    StyleHisto(histoYieldFinal[i], Low[TypeAnalysis] , Up[TypeAnalysis], ColorDiff[TypeAnalysis][i], Marker[i], Size[i], titleYieldX, titleYieldY, "" , 0,0, UpRangeMult);
    StyleHisto(histoYieldSistFinal[i], Low[TypeAnalysis] , Up[TypeAnalysis], ColorDiff[TypeAnalysis][i], Marker[i], Size[i], titleYieldX, titleYieldY, "" , 0,0, UpRangeMult);
    pol1[i] =     new TF1(Form("pol1_%i",i), "pol1", 2, UpRangeMult);
    if (i==0)     pol1[i]->SetLineStyle(9);
    else      pol1[i]->SetLineStyle(10);
    pol1[i]->SetLineColor(ColorDiff[TypeAnalysis][i]);
    histoYieldFinal[i]->Draw("same");
    cout << "Fit of yield " << i << endl;
    if (FitYields)    histoYieldFinal[i]->Fit(pol1[i], "R+");
    histoYieldSistFinal[i]->SetFillStyle(0);
    histoYieldSistFinal[i]->Draw("same e2");

    cout << "Bin width of YieldComparison " << histoYieldComparison->GetBinWidth(1) << endl;
    for (Int_t b=1; b<= histoYieldFinal[i]->GetNbinsX(); b++){
      if (histoYieldFinal[i]->GetBinContent(b)==0)	continue;
      bin = histoYieldComparison ->FindBin(histoYieldFinal[i]-> GetBinCenter(b));
      histoYieldComparison -> SetBinContent( bin , histoYieldFinal[i]->GetBinContent(b));
      histoYieldComparison -> SetBinError( bin , histoYieldFinal[i]->GetBinError(b));
      histoYieldSistComparison -> SetBinContent( bin , histoYieldSistFinal[i]->GetBinContent(b));
      histoYieldSistComparison -> SetBinError( bin , histoYieldSistFinal[i]->GetBinError(b));
    }
    for (Int_t b=1; b<= histoYieldComparison->GetNbinsX(); b++){
      if ( histoYieldComparison->GetBinContent(b) !=0){
	cout << "\nMultiplicity: "<<  histoYieldComparison->GetBinCenter(b) << " yield:" <<   histoYieldComparison->GetBinContent(b)<< " +- " << histoYieldComparison->GetBinError(b) << " (rel stat:" << histoYieldComparison->GetBinError(b)/histoYieldComparison->GetBinContent(b)<< " ) "  << endl;
	cout <<      histoYieldSistComparison->GetBinContent(b)<< " +- " << histoYieldSistComparison->GetBinError(b) << " (rel sist:" << histoYieldSistComparison->GetBinError(b)/histoYieldSistComparison->GetBinContent(b)<< " ) "  << endl;
      }
    }

    cout << "\nFit of total yield vs mult " << endl; 
    if (TypeAnalysis!=0 && FitYields && i==1)  {
      histoYieldComparison->Fit(pol1Common, "R+");
    }
  }

  TString NameFileout ="Compare" + YieldOrAvgPt[isAvgPtvsMult]+"DifferentCollisions";
  if (pathin[1].Index("2016k_HM_hK0s")!=-1) NameFileout += "_HM16k";
  if (pathin[1].Index("MultBinning1")!=-1) NameFileout += "_HMMultBinning1";
  if (pathin[0].Index("18g4extra")!=-1) NameFileout += "_MBEff18g4extra_NoEtaEff";
  NameFileout +=  CollisionsComp[isComp];
  NameFileout +="_"+ ParticleType[type] + Region[TypeAnalysis];
  if (isPreliminary)  NameFileout += "_isPreliminary";
  NameFileout += ".root";

  TString HistoTitle = "Yield vs multiplicity";
  if (isAvgPtvsMult) HistoTitle = "Avg pt vs multiplicity";
  StyleHisto(histoYieldComparison, Low[TypeAnalysis] , Up[TypeAnalysis], Color[TypeAnalysis], 33, 2, titleYieldX, titleYieldY, HistoTitle, 0,0, UpRangeMult);
  StyleHisto(histoYieldSistComparison, Low[TypeAnalysis] , Up[TypeAnalysis], Color[TypeAnalysis], 33, 2, titleYieldX, titleYieldY, HistoTitle, 0,0, UpRangeMult);

  canvas->cd();
  histoYieldComparison->Draw("");
  histoYieldSistComparison->SetFillStyle(0);
  histoYieldSistComparison->Draw("same e2");
  if (FitYields){
    pol1[0]->Draw("same");
    pol1[1]->Draw("same");
  }

  TFile * fileout = new TFile (NameFileout, "RECREATE");
  histoYieldComparison->Write("");
  histoYieldSistComparison->Write("");
  for (Int_t i=0; i<2; i++){
    histoYieldFinal[i]->Write();
    histoYieldSistFinal[i]->Write();
  }
  canvas->Write("");
  canvasBoth->Write("");
  fileout->Close();

  cout << "\n" << endl;
  canvasBoth->SaveAs("Compare" + YieldOrAvgPt[isAvgPtvsMult]+ "DifferentCollisions_"+ tipo[type]+Region[TypeAnalysis]+"_Superimp" + CollisionsComp[isComp]+SFitYields[FitYields]+SisPreliminary[isPreliminary]+".pdf");
  canvas->SaveAs("Compare" + YieldOrAvgPt[isAvgPtvsMult]+ "DifferentCollisions_"+ tipo[type]+Region[TypeAnalysis]+CollisionsComp[isComp] + SFitYields[FitYields]+SisPreliminary[isPreliminary]+".pdf");

  cout << "\n\e[35mStarting from the files: \e[39m" << endl;
  cout << pathin[0] << endl;
  cout << pathin[1] << endl;
  cout << "\n\e[35mI have produced the file\e[39m:\n" << NameFileout << "\n" << endl;
}
