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
#include "TLatex.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TFile.h>
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include "Macros/constants.h"
#include <Macros/ErrRatioCorr.C>
#include <Macros/BarlowVariable.C>

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title){
  histo->GetYaxis()->SetRangeUser(Low, Up);
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
  //
}

TString titleX=  "p_{T} (GeV/#it{c})";
TString titleY=  "#it{R}";
TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5, #it{p}_{T,trigg}>3 GeV/#it{c}}";

void MacroRRatio(Int_t type=8, Int_t TypeAnalysis=2, Bool_t ispp5TeV=0, Int_t TypeSel=0, Float_t NSigmaB=2, Bool_t isAlsoppHM=1, Bool_t isOnlyppHM=1){

  TString STypeSel[2] = {"_TopoSel", "_dPhiSel"};
  //TypeSel==0 --> topological selections (tight and loose)
  //TypeSel==1 --> dPhi selections 
  if (TypeAnalysis==2 && TypeSel==1) {
    cout << "Ehm... for full analysis you cannot change dPhi selection) " << endl;
    return;
  }
  if ((ispp5TeV || isAlsoppHM) && TypeSel==0) {
    cout << "Not implemented " << endl;
    return;
  } 
  const Int_t nummoltMax = 8;
  Int_t nummoltMaxEff=5;
  if (isAlsoppHM) nummoltMaxEff=8;
  else if (ispp5TeV) nummoltMaxEff=2;
  Int_t nummoltMinEff=0;
  if (isOnlyppHM) nummoltMinEff = 5;
  Int_t numVar = 2;

  /*
    if (TypeAnalysis==0) numVar += 1;
    else if (TypeAnalysis==1) numVar += 3;
  */
  if (TypeSel==1){
    if (TypeAnalysis==0) numVar = 1;
    else if (TypeAnalysis==1) numVar = 3;
  }
  TString tipo[]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
  TString Smolt[]={"50-100", "30-50", "10-30", "5-10", "0-5", "0.05-0.1", "0.01-0.05", "0-0.01"};
  if (ispp5TeV) {
    Smolt[0] = "10-100";
    Smolt[1] = "0-10";
  }
  Int_t Color[] = {401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  Int_t Marker[] = {20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  TString TitleVar[10] = {"Tight topological sel.", "Loose topological sel."};
  if (TypeSel==1){
    if (TypeAnalysis==0){
      TitleVar[0] = "|#Delta#varphi|<1.33";
    }
    else if (TypeAnalysis==1){
      TitleVar[0] = "0.96 < |#Delta#varphi| < 1.8";
      TitleVar[1] = "1.1 < |#Delta#varphi| < 2.05";
      TitleVar[2] = "0.96 < |#Delta#varphi| < 2.05";
    }
  }
  Float_t LimInfRatio = 0.95;
  Float_t LimSupRatio = 1.15;
  Float_t LimInfR = 0.97;
  Float_t LimSupR = 1.07;
  Float_t LimInfRelErr = 0.;
  Float_t LimSupRelErr = 0.03;

  TH1F * hSpectrumDef0100;
  TH1F * hYieldDef;
  TH1F * hSpectrumDef[nummoltMax];
  TH1F * hSpectrumVar0100[numVar];
  TH1F * hYieldVar[numVar];
  TH1F * hYieldRatio[numVar];
  TH1F * hYieldR[numVar];
  TH1F * hSpectrumRatio0100[numVar];
  TH1F * hSpectrumVar[numVar][nummoltMax];
  TH1F * hSpectrumRatio[numVar][nummoltMax];
  TH1F * hSpectrumR[numVar][nummoltMax];
  TH1F * hSpectrumRMaximum[nummoltMax];
  TH1F * hSpectrumRelErr[nummoltMax];
  TH1F * hSpectrumFracUncorr[nummoltMax];

  TString sFileDef = "";
  TString sFileVar = "";
  TString sFileDefYieldDPhi = "";
  TString sFileYieldDef = "";
  TString sFileYieldVar = "";
  TString sFileVarBase = "";
  TString sFileYieldVarBase = "";
  TFile* fileVar[numVar];
  TFile* fileYieldVar[numVar];

  TCanvas *canvasRatio[nummoltMax];
  TCanvas *canvasRatio0100 = new TCanvas("canvasRatio0100", "canvasRatio0100", 1000, 750);
  TCanvas *canvasYieldRatio= new TCanvas("canvasYieldRatio", "canvasYieldRatio", 1000, 750);
  TCanvas *canvasR[nummoltMax];
  TCanvas *canvasUnCorrRelUncertainty[nummoltMax];
  TCanvas *canvasFractionUnCorr[nummoltMax];
  TCanvas *canvasFractionUnCorrAllMult = new TCanvas("canvasFractionUnCorrAllMult", "canvasFractionUnCorrAllMult", 1000, 750);

  TLegend* legendALICE;
  TLegend* legendVar;
  TLegend* legendChi2;
  TLegend * legendMult = new TLegend(0.6, 0.1, 0.9, 0.3);

  TF1* linePtAtOne = new TF1("linePtAtOne", "pol0", 0, 8);
  linePtAtOne->SetLineStyle(9);
  linePtAtOne->SetLineWidth(1);
  linePtAtOne->SetLineColor(1);
  linePtAtOne->FixParameter(0, 1);
  TF1* lineMultAtOne = new TF1("lineMultAtOne", "pol0", 0, 45);
  lineMultAtOne->SetLineStyle(9);
  lineMultAtOne->SetLineWidth(1);
  lineMultAtOne->SetLineColor(1);
  lineMultAtOne->FixParameter(0, 1);

  TF1 * fitR[numVar][nummoltMax];

  TH1F* hSpectrumRelErrorSE[nummoltMax];
  TH1F* hSpectrumRelErrordPhi[nummoltMax];
  TFile * fileYieldMBDef;

  for (Int_t m=nummoltMinEff; m<nummoltMaxEff; m++){
    cout << "\n\e[36mMultiplicity class: " << Smolt[m] << "\e[37m"<<endl;
    sFileYieldDef = "FinalOutput/DATA2016/PtSpectraBisNew";
    if (type==8)  sFileYieldDef += "_161718Full_AOD234_hXi_Xi";
    else if (type==0) {
      if (ispp5TeV){
      }
      else {
	if (m<5) sFileYieldDef += "_PtBinning1_1617_AOD234_hK0s_K0s";
	else  sFileYieldDef += "_pp13TeVHM_PtBinning1_AllhK0sHM_RedNo16k_K0s";
      }
    }
    sFileYieldDef += "_Eta0.8_AllAssoc_PtMin3.0_";
    sFileYieldDef += RegionType[TypeAnalysis];
    sFileYieldDef += "_isNormCorrFullyComputed_isErrorAssumedPtCorr_ChangesIncluded_isdNdEtaTriggered";
    if (!ispp5TeV && m>=5)  sFileYieldDef += "_MultBinning1";
    sFileYieldDef += "_EffCorr_isWingsCorrectionAppliedNew";
    sFileYieldDef += "_MatBudgetCorr";
    if (type==8 || m>=5 || ispp5TeV) sFileYieldDef += "FAST";
    sFileYieldVarBase =  sFileYieldDef;
    //    sFileYieldDef += "_PForT.root"; //This file contains the histograms of the relative errors vs pt for the different sources
    sFileYieldDef += "_MultCorr.root"; //This file contains the histograms of the relative errors vs pt for the different sources
    TFile * fileYieldDef = new TFile(sFileYieldDef, "");
    cout << "Input file for relative errors vs pt: " << sFileYieldDef << endl;
    if (!fileYieldDef) {cout << "file " << sFileYieldDef << " not available " << endl; return;}
    if (m==nummoltMinEff) fileYieldMBDef = new TFile(sFileYieldDef, "");

    if (TypeSel==0) {
      hSpectrumRelErrorSE[m] = (TH1F*) fileYieldDef->Get("fHistSpectrumSistRelErrorSE_"+Smolt[m]);
      if (!hSpectrumRelErrorSE[m]) {cout << "Rel. error of topological selections not found" << endl; return; }
    }
    else if (TypeSel==1) {
      hSpectrumRelErrordPhi[m] = (TH1F*) fileYieldDef->Get("fHistSpectrumSistRelErrorDPhi_"+Smolt[m]);
      if (!hSpectrumRelErrordPhi[m]) {cout << "Rel. error of dPhi not found" << endl; return;} 
    }

    sFileDef = "FinalOutput/DATA2016/PtSpectraNew";
    if (type==8) {
      if (ispp5TeV){

      }
      else {
	if (m<5) sFileDef += "_161718Full_AOD234_hXi_Xi";
	else sFileDef += "";
      }
    }
    else if (type==0){
      if (ispp5TeV) sFileDef += "_17pq_hK0s_PtBinning1_K0s";
      else {
	if (m<5) sFileDef += "_1617_AOD234_hK0s_PtBinning1_K0s";
	else sFileDef += "_AllhK0sHM_RedNo16k_PtBinning1_K0s";
      }
    }
    sFileDef += "_Eta0.8_AllAssoc_SysPhi";
    sFileVarBase = sFileDef;
    sFileDef += "0_PtMin3.0_";
    sFileDef += RegionType[TypeAnalysis];
    if (!ispp5TeV && m>=5)  sFileDef += "_MultBinning1";
    if (type==0)     sFileDef += "_EffCorr_isWingsCorrectionAppliedNew";
    if (type==0 && m<5) sFileDef +="_MatBudgetCorr";
    sFileDefYieldDPhi = sFileDef + "_wYields.root";
    sFileDef += ".root";
    TFile * fileDef = new TFile(sFileDef, "");
    cout << "Input file with default spectra + yields (not integrated with pt): "<< sFileDef << endl;
    if (!fileDef) {cout << "file " << sFileDef << " not available " << endl; return;}

    if (m==nummoltMinEff) {
      if (isOnlyppHM)      hSpectrumDef0100 = (TH1F*)fileDef->Get("fHistSpectrum_0-0.1");
      else       hSpectrumDef0100 = (TH1F*)fileDef->Get("fHistSpectrum__all");
    }

    hSpectrumDef[m] =  (TH1F*)fileDef->Get("fHistSpectrum_"+ Smolt[m]);
    if (!hSpectrumDef0100) {cout << "Missing histo: hSpectrumDef0100" << endl; return; }
    if (!hSpectrumDef[m]) {cout << "Missing histo: hSpectrumDef[m]" << endl; return; }

    hSpectrumDef[m]->SetName("fHistSpectrum_"+ Smolt[m]+ "_Def"); 
    hSpectrumDef0100->SetName("fHistSpectrum__all_Def");

    canvasRatio[m]= new TCanvas(Form("canvasRatio_m%i", m), Form("canvasRatio_m%i", m), 1000, 750);
    canvasR[m]= new TCanvas(Form("canvasR_m%i", m), Form("canvasR_m%i", m), 1000, 750);
    canvasUnCorrRelUncertainty[m]= new TCanvas(Form("canvasUnCorrRelUncertainty_m%i", m), Form("canvasUnCorrRelUncertainty_m%i", m), 1000, 750);
    canvasFractionUnCorr[m]= new TCanvas(Form("canvasFractionUnCorr_m%i", m), Form("canvasFractionUnCorr_m%i", m), 1000, 750);

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);

    legendALICE= new TLegend(0.05, 0.7, 0.4, 0.9);
    legendALICE->SetFillStyle(0);
    legendALICE->SetTextAlign(12);
    legendALICE->SetTextSize(0.03);
    legendALICE->AddEntry("", "#bf{This work}", "");
    if (ispp5TeV) legendALICE->AddEntry("", "pp, #sqrt{#it{s}} = 5.02 TeV", "");
    else  legendALICE->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");

    TLegend *legendALICEBis = (TLegend*) legendALICE->Clone("legendALICEBis");
    legendALICE->AddEntry("", "V0M Multiplicity Percentile "+Smolt[m]+ " %", "");

    legendVar= new TLegend(0.6, 0.7, 0.9, 0.9);
    legendVar->SetFillStyle(0);
    legendVar->SetTextAlign(12);
    legendVar->SetTextSize(0.03);

    legendChi2= new TLegend(0.6, 0.6, 0.9, 0.7);
    legendChi2->SetFillStyle(0);
    legendChi2->SetTextAlign(12);
    legendChi2->SetTextSize(0.03);

    Int_t phiSel=0;
    for (Int_t var=0; var < numVar; var++){
      //      if (var>=2) phiSel = var-1;
      if (TypeSel==0) phiSel=0;
      else {
	phiSel=var+1;
	if (TypeAnalysis==0) phiSel+=1;
      }
      sFileVar = sFileVarBase +Form("%i_PtMin3.0_", phiSel);
      sFileVar += RegionType[TypeAnalysis];
      if (!ispp5TeV && m>=5)  sFileVar += "_MultBinning1";
      if (type==0) sFileVar += "_EffCorr_isWingsCorrectionAppliedNew";
      if (type==0 && m<5) sFileVar += "_MatBudgetCorr";
      if (TypeSel==0){
	if (var==0) sFileVar += "_SysV0Tightest";
	else if (var==1)  sFileVar += "_SysV0Loosest";
      }
      sFileVar += ".root";
      cout << "\n\e[35mVariation n: " << var << " file name: \e[39m\n" << sFileVar << endl;

      fileVar[var] = new TFile(sFileVar, "");
      if (!fileVar[var]) {cout << "file " << sFileVar << " not available " << endl; return;}

      if (m==nummoltMinEff){
	if (isOnlyppHM)	hSpectrumVar0100[var] = (TH1F*)fileVar[var]->Get("fHistSpectrum_0-0.1");
	else 	hSpectrumVar0100[var] = (TH1F*)fileVar[var]->Get("fHistSpectrum__all");
      }

      hSpectrumVar[var][m] =  (TH1F*)fileVar[var]->Get("fHistSpectrum_"+ Smolt[m]);
      if (!hSpectrumVar0100[var]) {cout << "Missing histo: hSpectrumVar0100" << endl; return; }
      if (!hSpectrumVar[var][m]) {cout << "Missing histo: hSpectrumVar[m]" << endl; return; }
      hSpectrumVar[var][m]->SetName("fHistSpectrum_"+ Smolt[m]+ Form("_Var%i", var)); 
      hSpectrumVar0100[var]->SetName(Form("fHistSpectrum__all_Var%i", var));

      hSpectrumRatio[var][m] = (TH1F*)hSpectrumVar[var][m]->Clone("fHistSpectrumRatio_"+ Smolt[m]+ Form("_Var%i", var)); 
      hSpectrumRatio0100[var] = (TH1F*)hSpectrumVar0100[var]->Clone(Form("fHistSpectrumRatio__all_Var%i", var)); 

      hSpectrumRatio[var][m]->Divide(hSpectrumDef[m]);
      hSpectrumRatio0100[var]->Divide(hSpectrumDef0100);
      StyleHisto(hSpectrumRatio[var][m], LimInfRatio, LimSupRatio, Color[var], Marker[var], titleX, "#it{N}_{var}/#it{N}_{def}", "");
      ErrRatioCorr(hSpectrumVar[var][m], hSpectrumDef[m], hSpectrumRatio[var][m], 0);
      ErrRatioCorr(hSpectrumVar0100[var], hSpectrumDef0100, hSpectrumRatio0100[var], 0);

      hSpectrumR[var][m] = (TH1F*) hSpectrumRatio[var][m]->Clone(Form("RRatio_Var%i_", var)+ Smolt[m]);
      hSpectrumR[var][m]->Divide(hSpectrumRatio0100[var]);

      StyleHisto(hSpectrumR[var][m], LimInfR, LimSupR, Color[var], Marker[var], titleX, titleY, "");
      ErrRatioCorr(hSpectrumRatio[var][m], hSpectrumRatio0100[var], hSpectrumR[var][m], 0);

      cout << "Variation " << TitleVar[var] << " content: " << hSpectrumR[var][m]->GetBinContent(2) << endl;
      fitR[var][m] = new TF1(Form("fitR_var%i_m%i", var, m), "pol0", 0, 8);
      fitR[var][m]->SetLineColor(Color[var]);
      hSpectrumR[var][m]->Fit(fitR[var][m], "RQ0");
      cout << "\n------------Fit results-------------" << endl;
      cout << "Var " << TitleVar[var] << ", mult. " << Smolt[m]<< "%"  << endl;
      cout << fitR[var][m]->GetParameter(0) <<"+-" << fitR[var][m]->GetParError(0) << " (" << TMath::Abs(1-fitR[var][m]->GetParameter(0))/fitR[var][m]->GetParError(0)<<" sigma from 1)" <<endl;
      cout << fitR[var][m]->GetChisquare() << "/" <<  fitR[var][m]->GetNDF() << " (" << fitR[var][m]->GetChisquare()/fitR[var][m]->GetNDF() << ")" << endl;
      cout << "------------------------------------" << endl;
      legendChi2->AddEntry(fitR[var][m], Form("Chi2/NDF = %.2f/%i", fitR[var][m]->GetChisquare(), fitR[var][m]->GetNDF()), "l");

      //      cout << "RMaximum defintion " << endl;
      if (Smolt[m] == "0-5" && TypeAnalysis==1) NSigmaB = 1;
      else NSigmaB = 2;
      if (var==0){
	//cout << "Var 0 " << endl;
	hSpectrumRMaximum[m] = (TH1F*)  hSpectrumR[var][m]->Clone(Form("RRatio_Maximum_m%i", m));
	for (Int_t b=1; b<= hSpectrumRMaximum[m]->GetNbinsX(); b++){
	  if ( hSpectrumRMaximum[m]->GetBinContent(b)!=0){
	    //cout << "pt " <<  hSpectrumRMaximum[m]->GetBinCenter(b) << " R: " <<  hSpectrumRMaximum[m]->GetBinContent(b) << endl;
	    //cout << TMath::Abs(1-hSpectrumRMaximum[m]->GetBinContent(b)) << " <? " << NSigmaB*hSpectrumRMaximum[m]->GetBinError(b) << endl;
	    if (TMath::Abs(1-hSpectrumRMaximum[m]->GetBinContent(b)) <  NSigmaB*hSpectrumRMaximum[m]->GetBinError(b)) {//Barlow check not passed
	      //cout << "pt " << hSpectrumRMaximum[m]->GetBinCenter(b) << " Not Barlow passed " << endl;
	      hSpectrumRMaximum[m]->SetBinContent(b,1);
	      hSpectrumRMaximum[m]->SetBinError(b,0);
	    }
	  }
	}
      }
      if (var>0){
	//cout << "Var 1 " << endl;
	for (Int_t b=1; b<= hSpectrumRMaximum[m]->GetNbinsX(); b++){
	  if ( hSpectrumR[var][m]->GetBinContent(b)!=0){
	    if (TMath::Abs(hSpectrumR[var][m]->GetBinContent(b)-1) >  TMath::Abs(hSpectrumRMaximum[m]->GetBinContent(b)-1)){
	      if (TMath::Abs(1-hSpectrumR[var][m]->GetBinContent(b)) >  NSigmaB*hSpectrumR[var][m]->GetBinError(b)) {//Barlow check passed
		//cout << "R: " << hSpectrumR[var][m]->GetBinContent(b) << "+-" << hSpectrumR[var][m]->GetBinError(b)<<  endl;
		//cout << "pt " << hSpectrumRMaximum[m]->GetBinCenter(b) << " Barlow passed " << endl;
		hSpectrumRMaximum[m]->SetBinContent(b,hSpectrumR[var][m]->GetBinContent(b));
		hSpectrumRMaximum[m]->SetBinError(b,hSpectrumR[var][m]->GetBinError(b));
	      }
	    }
	  }
	}
      }

      legendVar->AddEntry(hSpectrumRatio[var][m], TitleVar[var], "pl");

      //      cout << "Saving Ratio " << endl;
      canvasRatio[m]->cd();
      hSpectrumRatio[var][m]->Draw("same");
      linePtAtOne->Draw("same");
      if (var==numVar-1)  {
	legendALICE->Draw("same");
	legendVar->Draw("same");
      }

      //cout << "Saving RMaximum " << endl;
      StyleHisto(hSpectrumRMaximum[m], LimInfR, LimSupR, 1, 33, titleX, titleY, "");
      canvasR[m]->cd();
      hSpectrumR[var][m]->DrawClone("same");
      hSpectrumRMaximum[m]->Draw("same");
      fitR[var][m]->Draw("same");
      linePtAtOne->Draw("same");
      legendChi2->Draw("same");
      if (var==numVar-1)  {
	legendALICE->Draw("same");
	legendVar->Draw("same");
      }

      if (m==nummoltMinEff){
	canvasRatio0100->cd();
	StyleHisto(hSpectrumRatio0100[var], LimInfRatio, LimSupRatio, Color[var], 33, titleX, titleY, "");
	hSpectrumRatio0100[var]->Draw("same");
	if (var==numVar-1)  {
	  legendALICE->Draw("same");
	  legendVar->Draw("same");
	}
      }
    }

    canvasRatio[m]->SaveAs("canvasRatio_"+ tipo[type]+"_"+Smolt[m]+"_" + RegionType[TypeAnalysis] + STypeSel[TypeSel] +".pdf");
    canvasR[m]->SaveAs("canvasR_"+ tipo[type]+"_"+Smolt[m]+"_" + RegionType[TypeAnalysis] + STypeSel[TypeSel] +".pdf");

    //cout << "Saving Rel Errors " << endl;
    hSpectrumRelErr[m]= (TH1F*)    hSpectrumRMaximum[m]->Clone(Form("hSpectrumRelErr_m%i", m));
    for (Int_t b=1; b<= hSpectrumRMaximum[m]->GetNbinsX(); b++){
      if ( hSpectrumRMaximum[m]->GetBinContent(b)!=0){
	if (TypeSel==0)	hSpectrumRelErr[m]->SetBinContent(b, TMath::Abs(hSpectrumRMaximum[m]->GetBinContent(b)-1)/sqrt(12));
	else 	hSpectrumRelErr[m]->SetBinContent(b, TMath::Abs(hSpectrumRMaximum[m]->GetBinContent(b)-1)/2);
      }
    }
    hSpectrumRelErr[m]->SetMarkerStyle(33);
    hSpectrumRelErr[m]->SetLineColor(1);
    hSpectrumRelErr[m]->Smooth(1,"R");
    StyleHisto(hSpectrumRelErr[m], LimInfRelErr, LimSupRelErr, 1, 33, titleX, "Relative uncertainty", "");

    canvasUnCorrRelUncertainty[m]->cd();
    hSpectrumRelErr[m]->Draw("same p");
    if (TypeSel==0)    hSpectrumRelErrorSE[m]->Draw("same");
    else if (TypeSel==1) hSpectrumRelErrordPhi[m]->Draw("same");
    canvasUnCorrRelUncertainty[m]->SaveAs("canvasRelErr_"+ tipo[type]+"_"+Smolt[m]+"_" + RegionType[TypeAnalysis] + STypeSel[TypeSel] +".pdf");

    hSpectrumFracUncorr[m] = (TH1F*) hSpectrumRelErr[m]->Clone("hSpectrumFracUncorr_" + Smolt[m]);
    if (TypeSel==0)    hSpectrumFracUncorr[m]->Divide(hSpectrumRelErrorSE[m]);
    else     hSpectrumFracUncorr[m]->Divide(hSpectrumRelErrordPhi[m]);
    StyleHisto(hSpectrumFracUncorr[m], 0, 1, 1, 33, titleX, "Fraction of uncertainty uncorr with mult", "");
    for (Int_t b=1; b<=hSpectrumFracUncorr[m]->GetNbinsX(); b++){
      if (hSpectrumFracUncorr[m]->GetBinContent(b) > 1)      hSpectrumFracUncorr[m]->SetBinContent(b,1);
    }
    hSpectrumFracUncorr[m]->Smooth(1, "R");
    canvasFractionUnCorr[m]->cd();
    hSpectrumFracUncorr[m]->DrawClone("same hist");
    canvasFractionUnCorr[m]->SaveAs("canvasFracUncorr_"+ tipo[type]+"_"+Smolt[m]+"_" + RegionType[TypeAnalysis] + STypeSel[TypeSel] +".pdf");

    canvasFractionUnCorrAllMult->cd();
    StyleHisto(hSpectrumFracUncorr[m], 0, 1, Color[m], 33, titleX, "Fraction of uncertainty uncorr with mult", "");
    legendMult->AddEntry(hSpectrumFracUncorr[m], Smolt[m], "l");
    hSpectrumFracUncorr[m]->SetLineWidth(2);
    hSpectrumFracUncorr[m]->Draw("same hist");
  }

  canvasFractionUnCorrAllMult->cd();
  legendMult->Draw("same");
  canvasFractionUnCorrAllMult->SaveAs("canvasFracUncorr_"+ tipo[type]+"_AllMult_" + RegionType[TypeAnalysis] + STypeSel[TypeSel] +".pdf");

  //Yield variation vs topo sel
  if (TypeSel==0){
    cout << "\n\nYield vs multiplicity for different topo selections " << endl;
    hYieldDef = (TH1F*)fileYieldMBDef->Get("fHistYieldStat");
    if (!hYieldDef) {cout << "Missing histo: hYieldDef" << endl; return; }
    hYieldDef->SetName("hYieldDef");

    for (Int_t var=0; var < numVar; var++){
      if (var==0) sFileYieldVar = sFileYieldVarBase + "_SysV0Tightest";
      else if (var>=1)  sFileYieldVar = sFileYieldVarBase + "_SysV0Loosest"; 
      sFileYieldVar += "_PForT.root";  //This file contains the histograms of the relative errors vs pt for the different sources
      cout << "Variation n: " << var << " file name: " << sFileYieldVar<< endl;

      fileYieldVar[var] = new TFile(sFileYieldVar, "");
      if (!fileYieldVar[var]) {cout << "file " << sFileYieldVar << " not available " << endl; return;}
      hYieldVar[var] = (TH1F*)fileYieldVar[var]->Get("fHistYieldStat");
      if (!hYieldVar[var]) {cout << "Missing histo: hYieldVar" << endl; return; }
      hYieldVar[var]->SetName(Form("fHistYieldStat_Var%i", var));
      hYieldRatio[var] = (TH1F*)hYieldVar[var]->Clone(Form("hYieldRatio_Var%i", var)); 
      hYieldRatio[var]->Divide(hYieldDef);
      ErrRatioCorr(hYieldVar[var], hYieldDef, hYieldRatio[var], 0);
      /*
      hYieldR[var] = (TH1F*)hYieldVar[var]->Clone(Form("hYieldRatio_Var%i", var)); 
      hYieldR[var]->Divide(hYieldDef);
      ErrRatioCorr(hYieldVar[var], hYieldDef, hYieldRatio[var], 0);
      */
      StyleHisto(hYieldRatio[var], 0.95, 1.05, Color[var], Marker[var], titledNdeta, "#it{N}_{var}/#it{N}_{def}", "");
      canvasYieldRatio->cd();
      hYieldRatio[var]->Draw("same");
      hYieldDef->Draw("same");
      lineMultAtOne->Draw("same");
      if (var==numVar-1)  {
	//	legendALICEBis->Draw("same");
	legendVar->Draw("same");
      }
    }
  }

  //Yield variation vs phi (yield w/o extrapolation)
  Int_t phiSel=0;
  if (TypeSel==1){
    cout << "\n\n\e[35mYield vs multiplicity for different dPhi selections (yield not extrapolated at low pt)\e[35m" << endl;
    TFile * fileDefYieldDPhi = new TFile(sFileDefYieldDPhi, "");
    if (!fileDefYieldDPhi) {cout << "file " << sFileDefYieldDPhi << " not available " << endl; return;}
    cout << "\n\e[35mDefault file name: \e[39m\n" << sFileDefYieldDPhi << endl;
    hYieldDef = (TH1F*)fileDefYieldDPhi->Get("fHistYieldStatD");
    if (!hYieldDef) {cout << "Missing histo: hYieldDef" << endl; return; }
    hYieldDef->SetName("hYieldDef");

    for (Int_t var=0; var < numVar; var++){
      phiSel = var+1;
      if (TypeAnalysis==0) phiSel+=1;
      sFileYieldVar = sFileVarBase +Form("%i_PtMin3.0_", phiSel);
      sFileYieldVar += RegionType[TypeAnalysis];
      if (!ispp5TeV && isAlsoppHM)  sFileYieldVar += "_MultBinning1";
      if (type==0) sFileYieldVar += "_EffCorr_isWingsCorrectionAppliedNew";
      if (!ispp5TeV && !isAlsoppHM) sFileYieldVar +="_MatBudgetCorr";
      sFileYieldVar += "_wYields.root";
      cout << "\n\e[35mVariation n: " << var << " file name: \e[39m\n" << sFileYieldVar << endl;
      fileYieldVar[var] = new TFile(sFileYieldVar, "");
      if (!fileYieldVar[var]) {cout << "file " << sFileYieldVar << " not available " << endl; return;}
      hYieldVar[var] = (TH1F*)fileYieldVar[var]->Get("fHistYieldStatD");
      if (!hYieldVar[var]) {cout << "Missing histo: hYieldVar" << endl; return; }
      hYieldVar[var]->SetName(Form("fHistYieldStat_Var%i", var));
      hYieldRatio[var] = (TH1F*)hYieldVar[var]->Clone(Form("hYieldRatio_Var%i", var)); 
      hYieldRatio[var]->Divide(hYieldDef);
      ErrRatioCorr(hYieldVar[var], hYieldDef, hYieldRatio[var], 0);
      StyleHisto(hYieldRatio[var], 0.95, 1.05, Color[var], Marker[var], titledNdeta, "#it{N}_{var}/#it{N}_{def}", "");
      canvasYieldRatio->cd();
      hYieldRatio[var]->Draw("same");
      hYieldDef->Draw("same");
      lineMultAtOne->Draw("same");
      if (var==numVar-1)  {
	legendVar->Draw("same");
      }
    }
  }

  canvasYieldRatio->SaveAs("canvasYieldRatio_"+ tipo[type]+"_"+ RegionType[TypeAnalysis] + STypeSel[TypeSel] + ".pdf");
  canvasRatio0100->SaveAs("canvasRatio0100_"+ tipo[type]+"_"+ RegionType[TypeAnalysis] + STypeSel[TypeSel] + ".pdf");
  cout << "\n\n\e[35mDefault path in input:\e[39m " << endl;
  cout << sFileDef << endl;
  if (TypeSel==0)  cout << sFileYieldDef << endl;
  cout << "\n\n" << endl;

  TString sFileOut = "UncertaintiesUncorrMult_" + tipo[type]+"_"+ RegionType[TypeAnalysis] + STypeSel[TypeSel];
  if (isOnlyppHM) sFileOut += "-isOnlyppHM";
  sFileOut += ".root";
  TFile * fileout = new TFile(sFileOut, "RECREATE");
  for (Int_t m=nummoltMinEff; m<nummoltMaxEff; m++){
    hSpectrumFracUncorr[m]->Write();
  }
  fileout->Close();
  cout << "Ho anche creato il file: " << sFileOut << endl;

}
