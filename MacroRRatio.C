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

void MacroRRatio(Int_t type=8, Int_t TypeAnalysis=2, Bool_t ispp5TeV=0, Int_t TypeSel=1){

  //TypeSel==0 --> topological selections (tight and loose)
  //TypeSel==1 --> dPhi selections 
  if (TypeAnalysis==2 && TypeSel==1) return;
 
  const Int_t nummoltMax = 5;
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
  TString SmoltLegend[]={"50-100", "30-50", "10-30", "5-10", "0-5", "0.05-0.1", "0.01-0.05", "0-0.01"};
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
  Float_t LimInfR = 0.95;
  Float_t LimSupR = 1.15;

  TH1F * hSpectrumDef0100;
  TH1F * hYieldDef;
  TH1F * hSpectrumDef[nummoltMax];
  TH1F * hSpectrumVar0100[numVar];
  TH1F * hYieldVar[numVar];
  TH1F * hYieldRatio[numVar];
  TH1F * hSpectrumRatio0100[numVar];
  TH1F * hSpectrumVar[numVar][nummoltMax];
  TH1F * hSpectrumRatio[numVar][nummoltMax];
  TH1F * hSpectrumR[numVar][nummoltMax];
  TH1F * hSpectrumRMaximum[nummoltMax];
  TH1F * hSpectrumRelErr[nummoltMax];

  TString sFileDef = "";
  TString sFileVar = "";
  TString sFileYieldDef = "";
  TString sFileYieldVar = "";
  TString sFileVarBase = "";
  TString sFileYieldVarBase = "";
  TFile* fileVar[numVar];
  TFile* fileYieldVar[numVar];

  TCanvas *canvasRatio[nummoltMax];
  TCanvas *canvasYieldRatio= new TCanvas("canvasYieldRatio", "canvasYieldRatio", 1000, 750);
  TCanvas *canvasR[nummoltMax];
  TCanvas *canvasUnCorrRelUncertainty[nummoltMax];

  TLegend* legendALICE;
  TLegend* legendVar;
  TLegend* legendChi2;

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

  for (Int_t m=0; m<nummoltMax; m++){
    cout << "\n\e[36mMultiplicity class: " << SmoltLegend[m] << "\e[37m"<<endl;
    sFileDef = "FinalOutput/DATA2016/PtSpectraNew";
    if (type==8) {
      if (m<5) sFileDef += "_161718Full_AOD234_hXi_Xi";
      else sFileDef += "";
    }
    sFileDef += "_Eta0.8_AllAssoc_SysPhi";
    sFileVarBase = sFileDef;
    sFileDef += "0_PtMin3.0_";
    sFileDef += RegionType[TypeAnalysis];
    sFileDef += ".root";
    TFile * fileDef = new TFile(sFileDef, "");
    if (!fileDef) {cout << "file " << sFileDef << " not available " << endl; return;}

    if (m==0) {
      hSpectrumDef0100 = (TH1F*)fileDef->Get("fHistSpectrum__all");
    }

    hSpectrumDef[m] =  (TH1F*)fileDef->Get("fHistSpectrum_"+ Smolt[m]);
    if (!hSpectrumDef0100) {cout << "Missing histo: hSpectrumDef0100" << endl; return; }
    if (!hSpectrumDef[m]) {cout << "Missing histo: hSpectrumDef[m]" << endl; return; }

    hSpectrumDef[m]->SetName("fHistSpectrum_"+ Smolt[m]+ "_Def"); 
    hSpectrumDef0100->SetName("fHistSpectrum__all_Def");

    canvasRatio[m]= new TCanvas(Form("canvasRatio_m%i", m), Form("canvasRatio_m%i", m), 1000, 750);
    canvasR[m]= new TCanvas(Form("canvasR_m%i", m), Form("canvasR_m%i", m), 1000, 750);
    canvasUnCorrRelUncertainty[m]= new TCanvas(Form("canvasUnCorrRelUncertainty_m%i", m), Form("canvasUnCorrRelUncertainty_m%i", m), 1000, 750);

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
    legendALICE->AddEntry("", "V0M Multiplicity Percentile "+SmoltLegend[m]+ " %", "");

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
      if (TypeSel==0)phiSel=0;
      else phiSel=var;
      sFileVar = sFileVarBase +Form("%i_PtMin3.0_", phiSel);
      sFileVar += RegionType[TypeAnalysis];
      if (var==0) sFileVar += "_SysV0Tightest";
      else if (var==1)  sFileVar += "_SysV0Loosest";
      sFileVar += ".root";
      cout << "\n\e[35mVariation n: " << var << " file name: \e[39m\n" << sFileVar << endl;

      fileVar[var] = new TFile(sFileVar, "");
      if (!fileVar[var]) {cout << "file " << sFileVar << " not available " << endl; return;}

      if (m==0){
	hSpectrumVar0100[var] = (TH1F*)fileVar[var]->Get("fHistSpectrum__all");
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
      StyleHisto(hSpectrumRatio[var][m], LimInfR, LimSupR, Color[var], Marker[var], titleX, "#it{N}_{var}/#it{N}_{def}", "");
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

      if (var==0)   {
	hSpectrumRMaximum[m] = (TH1F*)  hSpectrumR[var][m]->Clone(Form("RRatio_Maximum_m%i", m));
	for (Int_t b=1; b<= hSpectrumRMaximum[m]->GetNbinsX(); b++){
	  if ( hSpectrumRMaximum[m]->GetBinContent(b)!=0){
	    if (TMath::Abs(1-hSpectrumRMaximum[m]->GetBinContent(b)) < hSpectrumRMaximum[m]->GetBinError(b)) {//Barlow check not passed
	      hSpectrumRMaximum[m]->SetBinContent(b,1);
	    }
	  }
	}
      }
      if (var>0){
	for (Int_t b=1; b<= hSpectrumRMaximum[m]->GetNbinsX(); b++){
	  if ( hSpectrumR[var][m]->GetBinContent(b)!=0){
	    if (TMath::Abs(hSpectrumR[var][m]->GetBinContent(b)-1) >  TMath::Abs(hSpectrumRMaximum[m]->GetBinContent(b)-1)){
	      if (TMath::Abs(1-hSpectrumR[var][m]->GetBinContent(b)) > hSpectrumR[var][m]->GetBinError(b)) {//Barlow check passed
		hSpectrumRMaximum[m]->SetBinContent(b,hSpectrumR[var][m]->GetBinContent(b));
	      }
	    }
	  }
	}
      }

      legendVar->AddEntry(hSpectrumRatio[var][m], TitleVar[var], "pl");

      canvasRatio[m]->cd();
      // hSpectrumRatio[var][m]->Draw("same");
      linePtAtOne->Draw("same");
      if (var==numVar-1)  {
	legendALICE->Draw("same");
	legendVar->Draw("same");
      }
      canvasRatio[m]->Close();

      StyleHisto(hSpectrumRMaximum[m], 0, 0.2, 1, 33, titleX, titleY, "");
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
    }

    canvasRatio[m]->SaveAs("canvasRatio_"+ tipo[type]+"_"+Smolt[m]+"_" + RegionType[TypeAnalysis] + ".pdf");
    canvasR[m]->SaveAs("canvasR_"+ tipo[type]+"_"+Smolt[m]+"_" + RegionType[TypeAnalysis] + ".pdf");

    hSpectrumRelErr[m]= (TH1F*)    hSpectrumRMaximum[m]->Clone(Form("hSpectrumRelErr_m%i", m));
    for (Int_t b=1; b<= hSpectrumRMaximum[m]->GetNbinsX(); b++){
      if ( hSpectrumRMaximum[m]->GetBinContent(b)!=0){
	hSpectrumRelErr[m]->SetBinContent(b, TMath::Abs(hSpectrumRMaximum[m]->GetBinContent(b)-1));
      }
    }
    canvasUnCorrRelUncertainty[m]->cd();
    hSpectrumRelErr[m]->Draw("same");
  }


  //Yield variation vs topo sel
  if (TypeSel==0){
    sFileYieldDef = "FinalOutput/DATA2016/PtSpectraBisNew";
    if (type==8){
      sFileYieldDef += "_161718Full_AOD234_hXi_Xi";
    }
    sFileYieldDef += "_Eta0.8_AllAssoc_PtMin3.0_";
    sFileYieldDef += RegionType[TypeAnalysis];
    sFileYieldDef += "_isNormCorrFullyComputed_isErrorAssumedPtCorr_ChangesIncluded_isdNdEtaTriggered_MatBudgetCorrFAST";
    sFileYieldVarBase =  sFileYieldDef;
    sFileYieldDef += ".root";
    TFile * fileYieldDef = new TFile(sFileYieldDef, "");
    if (!fileYieldDef) {cout << "file " << sFileYieldDef << " not available " << endl; return;}
    hYieldDef = (TH1F*)fileYieldDef->Get("fHistYieldStat");
    if (!hYieldDef) {cout << "Missing histo: hYieldDef" << endl; return; }
    hYieldDef->SetName("hYieldDef");

    for (Int_t var=0; var < numVar; var++){
      if (var==0) sFileYieldVar = sFileYieldVarBase + "_SysV0Tightest";
      else if (var>=1)  sFileYieldVar = sFileYieldVarBase + "_SysV0Loosest"; //ATTENTION: I am not comparing yields with different dPhi
      sFileYieldVar += ".root";
      cout << "Variation n: " << var << " file name: " << sFileYieldVar<< endl;

      fileYieldVar[var] = new TFile(sFileYieldVar, "");
      if (!fileYieldVar[var]) {cout << "file " << sFileYieldVar << " not available " << endl; return;}
      hYieldVar[var] = (TH1F*)fileYieldVar[var]->Get("fHistYieldStat");
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
	//	legendALICEBis->Draw("same");
	legendVar->Draw("same");
      }
    }
  }

  canvasYieldRatio->SaveAs("canvasYieldRatio_"+ tipo[type]+"_"+ RegionType[TypeAnalysis] + ".pdf");
  cout << "\n\n\e[35mDefault path in input:\e[39m " << endl;
  cout << sFileDef << endl;
  if (TypeSel==0)  cout << sFileYieldDef << endl;
  cout << "\n\n" << endl;
}
