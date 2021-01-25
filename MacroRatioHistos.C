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

void MacroRatioHistos(const Int_t numFiles = 8, TString  CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelationRun2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_", /*name common to all files */  TString OutputName ="hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3.root"){

  //this macros superimpose in a same canvas pad same histos found in different files. A loop over the files is done. The ratio of the histos to the histo found in the first file is performed. You can choose if histos have to be considered fully correlated, uncorrelated, or if the histos are obtaiend from a subsample of the data used to obtain the histo found in the first file (i.e. partial correlation)

  TString VarName[numFiles] = {""}; //might be defined within the loop
  TString histoName= "fHistQA6";
  Float_t numDef=3;
  Float_t num=0;

  //histo style selections
  Float_t Low=10e-8;
  Float_t Up=0.003;
  Int_t color[numFiles]={1,402 , 628, 905, 881,601, 867, 418 };
  Int_t style =33;
  TString titleX = "Multiplicity class";
  TString titleY="#hXi events / #INT7 events";
  TString title="Fraction of INT7 events containing trigger particle and Xi";
  TString titleRatio="Ratio to #it{p}_{T}^{trigg} > 3 GeV/#it{c}";
 
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
    VarName[i] = Form("PtMin%.1f", num);
    InputName = CommonFileName + VarName[i]+"_IsEstimateRun3.root";
    InputFile = new TFile (InputName, "");
    if (!InputFile) return;
    histo[i] = (TH1F*)InputFile->Get(histoName);
    if (!histo[i]) {cout << "histogram is not there " << endl; return;}
    histo[i]->Sumw2();

    canvas->cd(1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    StyleHistoYield(histo[i], Low, Up, color[i], style, titleX, titleY, title, 1, 1.2, 1.4);
    LegendName[i] = Form("%.1f GeV/#it{c}", num);
    legend->AddEntry(histo[i], LegendName[i], "pl");
    histo[i]->Draw("same");
    if (i==numFiles-1) legend->Draw("");

    canvas->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    histoRatio[i] = (TH1F*) histo[i]->Clone(histoName + "_Ratio");
    if (i!=0)    histoRatio[i]->Divide(histo[i-1]);
    StyleHistoYield(histoRatio[i], Low, 0.6, color[i], style, titleX, "Ratio", titleRatio, 1, 1.2, 1.4);
    if (i!=0)  histoRatio[i]->Draw("same");
    if (i==numFiles-1) legend->Draw("");

  }

  OutputFile->WriteTObject(canvas);
  cout << "I produced the output file " << OutputName << endl;

}
