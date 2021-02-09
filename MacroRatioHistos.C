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
#include <TFile.h>
#include <TLegend.h>
Double_t SetEfficiencyError(Int_t k, Int_t n){
  return sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
}


void ErrRatioCorr(TH1F* hNum, TH1F* hDenom, TH1F* hRatio, Bool_t FullCorr){
  //FullCorr == 1 means ro = 1;                                                                                                           
  //FullCorr == 0 means ro =                                                                                                              
  Float_t Err1=0;
  Float_t Err2=0;
  Float_t ErrC=0;
  Float_t Err=0;
  for (Int_t b=1; b<=hNum->GetNbinsX();b++){
    if (hNum->GetBinContent(b)==0 ||hDenom->GetBinContent(b)==0){
      hRatio->SetBinError(b,0);
      continue;
    }
    Err1=pow(hNum->GetBinError(b)/hNum->GetBinContent(b),2);
    Err2=pow(hDenom->GetBinError(b)/hDenom->GetBinContent(b),2);
    if (FullCorr==0){
      ErrC=pow(hDenom->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    else {
      ErrC=hDenom->GetBinError(b) * hNum->GetBinError(b)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    Err=sqrt(Err1+Err2-ErrC);
    hRatio->SetBinError(b,Err*hRatio->GetBinContent(b));
  }
  //  return hRatio;                                                                                                                      
}

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); //1.2                                                                                       
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void MacroRatioHistos(Int_t RunVar=4, TString  CommonFileName = ""/*name common to all files */,  TString OutputName =""){

  //RunVar should be increased when you want to do a new comparison; also, OutputName And CommonFileName should be updated below!!
  Int_t numFiles=0;
  if (RunVar==0){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelationRun2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3.root";
    numFiles=8;
  }
  else if (RunVar==1){
    CommonFileName= "FinalOutput/DATA2016/histo/AngularCorrelation2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName="hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOB.root";
    numFiles=4;
  }
  else if (RunVar==2){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOB_Comparison.root";
    numFiles=8;
  }

  else if (RunVar==3){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius_Comparison.root";
    numFiles=12;
  }

  else if (RunVar==4){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius.root";
    numFiles=6;
  }

  else if (RunVar==5){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius_SkipAllAssocComp.root";
    numFiles=8;
  }

  //this macros superimpose in a same canvas pad same histos found in different files. A loop over the files is done. The ratio of the histos to the histo found in the first file is performed. You can choose if histos have to be considered fully correlated, uncorrelated, or if the histos are obtaiend from a subsample of the data used to obtain the histo found in the first file (i.e. partial correlation)

  //  TString VarName[numFiles] = {"2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin4.0", "Run2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin4.0"}; //might be defined within the loop
  TString VarName[numFiles] = {""}; //might be defined within the loop
  TString histoName= "fHistQA6";
  Float_t numDef=3;
  Float_t num=0;
  Int_t numEff=0;

  //histo style selections
  Float_t Low=10e-8;
  Float_t Up=0.003;
  Float_t LowRatio=10e-8;
  Float_t UpRatio=0.6;
  Int_t color[10]={1,402 , 628, /*905, 881,*/601, 867, 418, 905, 881 };
  Int_t style =33;
  TString titleX = "Multiplicity class";
  TString titleY="#hXi events / #INT7 events";
  TString title="Fraction of INT7 events containing trigger particle and Xi";
  TString titleRatio="Ratio to #it{p}_{T}^{trigg} > 3 GeV/#it{c}";
 
  if (RunVar==2){
    titleRatio = "Ratio of only TOF to (SPD+TOF) pileup rej";
    LowRatio = 0.5;
    UpRatio = 1;
  }

  else  if (RunVar==3){
    titleRatio = "Ratio of Run3 to Run 2 selections";
    LowRatio = 0.5;
    UpRatio = 1;
  }

  else  if (RunVar==5){
    titleRatio = "Ratio of SkipAssoc to AllAssoc with Run3Sel";
    LowRatio = 0.8;
    UpRatio = 1.3;
  }

  TLegend * legend = new TLegend (0.6, 0.7, 0.9, 0.9);
  legend->SetHeader("#it{p}_{T}^{trigg} > ");
  TString LegendName[numFiles]={""};

  TH1F * histo[numFiles];
  TH1F * histoRatio[numFiles];

  TString InputName="";
  TFile * InputFile;
  TFile * OutputFile= new TFile (OutputName, "RECREATE");

  gStyle->SetOptStat(0);
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  canvas->Divide(2,1);

  for (Int_t i=0; i<numFiles; i++){
    num = numDef+i;
    numEff= num;
    if (RunVar==2 || RunVar==3 || RunVar==5) {
      if (i<numFiles/2)      numEff = num;
      else if (i>=numFiles/2) numEff=num-numFiles/2;
    }
    cout << numEff << endl;
    if (RunVar==0 || RunVar==1 || RunVar==4)    VarName[i] = Form("PtMin%.1f", num);
    else if (RunVar==2){
      if (i>=0 && i<numFiles/2) VarName[i] = Form("Run2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
      else if (i>=numFiles/2) VarName[i] =  Form("2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==3){
      if (i>=0 && i<numFiles/2) VarName[i] = Form("Run2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
      else if (i>=numFiles/2) VarName[i] =  Form("2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==5){
      if (i<=numFiles/2) VarName[i]+= "_AllAssoc";
      VarName[i] += Form("_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    InputName = CommonFileName + VarName[i];
    //    if (RunVar==4) InputName += "_IsOnlypiKpemu";
    InputName+="_IsEstimateRun3.root";
    cout << " loop n. " << i << " file name: " << InputName << endl;
    InputFile = new TFile (InputName, "");
    if (!InputFile) return;
    histo[i] = (TH1F*)InputFile->Get(histoName);
    if (!histo[i]) {cout << "histogram is not there " << endl; return;}
    histo[i]->Sumw2();

    canvas->cd(1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    if ((RunVar==2 || RunVar==3 || RunVar==5) && i>=numFiles/2) style = 27;
    StyleHistoYield(histo[i], Low, Up, color[numEff-3], style, titleX, titleY, title, 1, 1.2, 1.4);
    LegendName[i] = Form("%.1f GeV/#it{c}", num);
    if (RunVar==2){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}, only TOF PU rej", (float)numEff);
    }
    if (RunVar==3){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}, Run3 sel", (float)numEff);
    }
    if (RunVar==5){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}, AllAsso", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
    }

    legend->AddEntry(histo[i], LegendName[i], "pl");
    histo[i]->Draw("same");
    if (i==numFiles-1) legend->Draw("");

    canvas->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    histoRatio[i] = (TH1F*) histo[i]->Clone(histoName + "_Ratio");

    if (RunVar==2 || RunVar==3 || RunVar==5){
      if (i>=numFiles/2){
	histoRatio[i]->Divide(histo[i-numFiles/2]);
	ErrRatioCorr(histo[i], histo[i-numFiles/2], histoRatio[i], 0);
      }
    }
    else {
      if (i!=0){
	histoRatio[i]->Divide(histo[0]);
	ErrRatioCorr(histo[i], histo[0], histoRatio[i], 0);
      }
    }

    /*
    for (Int_t b=1; b<= histoRatio[i]->GetNbinsX(); b++){
      cout << "error before b:"<< b << " "  << histoRatio[i]->GetBinError(b) << endl;
    }
    
    if (RunVar==2 || RunVar==3){
      if (i>=numFiles/2){
	ErrRatioCorr(histo[i], histo[i-numFiles/2], histoRatio[i], 0);
      }
    }
    else {
      if (i!=0){
	ErrRatioCorr(histo[i], histo[0], histoRatio[i], 0);
      }
    }

    for (Int_t b=1; b<= histoRatio[i]->GetNbinsX(); b++){
      cout << "error after b:"<< b << " "  << histoRatio[i]->GetBinError(b) << endl;
    }
    */

    StyleHistoYield(histoRatio[i], LowRatio, UpRatio, color[numEff-3], style, titleX, "Ratio", titleRatio, 1, 1.2, 1.4);
    if (RunVar==2 || RunVar==3 || RunVar==5) {
      if (i>=numFiles/2)      histoRatio[i]->Draw("same");
    }
    else {
    if (i!=0)  histoRatio[i]->Draw("same");
    if (i==numFiles-1) legend->Draw("");
    }

  }

  OutputFile->WriteTObject(canvas);
  OutputFile->Close();
  cout << "I produced the output file " << OutputName << endl;

}
