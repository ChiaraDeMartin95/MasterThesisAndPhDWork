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
#include <TLegendEntry.h>
#include <TFile.h>
#include "TFitResult.h"
#include "Macros/SetEfficiencyError.C"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);

}

TString TitleYieldRatio="#Xi/K^{0}_{S} yield ratio vs multiplicity";
TString titleRelUnc = "Relative uncertainty";
TString titleMult = "Multiplicity class ";
TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
TString titleYield[2]={"K^{0}_{S}", "#Xi"};
TString TitleYYieldRatio="#it{N}_{#Xi}/#it{N}_{K^{0}_{S}}";

//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}"};
TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})", "#it{N}_{#Xi}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})"};
//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)"}; //to much space around the parenthesis
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
//TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
//TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};
//TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 0.85", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 2.0", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};


void FractionEventswTrigger(){

  TString SInputFile[3] = {"FinalOutput/AnalysisResults1617_AOD234_hK0s.root", "FinalOutput/AnalysisResults17pq_hK0s.root", "FinalOutput/AnalysisResultsAllhK0sHM_RedNo16k.root"};
  TString DirName[3] = {"MyTask_PtTrigMin3.0_PtTrigMax15.0", "MyTask_PtTrigMin3.0_PtTrigMax15.0", "MyTask_PtTrigMin3.0_PtTrigMax15.0"};
  TString ListName[3] = {"MyOutputContainer_hK0s_Task_", "MyOutputContainer_hK0s_Task_", "MyOutputContainer_hK0s_Task_"};

  TH1F * hEventsvsMult[3];
  TH1F * hEventsvsMult_EvwTrigger[3];

  Int_t nummolt = 5;
  Float_t Nmolt[3][nummolt+1]={{0,5,10,30,50,100}, {0, 10, 100, 100, 100, 100}, {0, 0.01, 0.05, 0.1, 0.1, 0.1}};
  Float_t dNdEta[3][nummolt+1]={{21.2, 16.17, 11.4625, 7.135, 3.33, 6.94}, {13.89, 6.95, 0, 0, 0, 0}, {36.29, 32.57, 30.43, 0, 0, 0}};
  TH1F * hRatio = new TH1F ("hRatio", "hRatio",450, 0, 45);
  TH1F * hRatio5TeV = new TH1F ("hRatio5TeV", "hRatio5TeV",450, 0, 45);
  Float_t LimInf =0; 
  Float_t LimSup =0; 
  Float_t FracEvWTrigger =0; 
  TString System[3] = {"ppMB13TeV",  "pp5TeV", "ppHM13TeV"};
  TString NameEventsvsMult[3] = {"AllEvents_ppMB13TeV",  "AllEvents_pp5TeV", "AllEvents_ppHM13TeV"};
  TString NameEventsvsMult_EvwTrigger[3] = {"EventsWTrigger_ppMB13TeV", "EventsWTrigger_pp5TeV",  "EventsWTrigger_ppHM13TeV"};

  TFile * InputFile;
  for (Int_t i =0; i<3; i++){
    cout <<"Analysing : " << System[i] << endl;
    InputFile = new TFile(SInputFile[i], "");
    if (!InputFile) return;
    TDirectoryFile * dir = (TDirectoryFile*)InputFile->Get(DirName[i]); 
    if (!dir) return;
    TList * list = (TList*) dir->Get(ListName[i]);
    if (!list) return;
    hEventsvsMult[i] = (TH1F*)list->FindObject("fHist_multiplicityAllSelEvents");
    if (!hEventsvsMult[i]) return;
    hEventsvsMult[i]->SetName(NameEventsvsMult[i]);
    hEventsvsMult_EvwTrigger[i] = (TH1F*)list->FindObject("fHist_multiplicity_EvwTrigger");
    if (!hEventsvsMult_EvwTrigger[i]) return;
    hEventsvsMult_EvwTrigger[i]->SetName(NameEventsvsMult_EvwTrigger[i]);
    for (Int_t m = 0; m< nummolt; m++){
      if (i==1 && m>1) continue;
      if (i==2 && m>2) continue;
      LimInf = hEventsvsMult_EvwTrigger[i]->FindBin(Nmolt[i][m]);
      LimSup = hEventsvsMult_EvwTrigger[i]->FindBin(Nmolt[i][m+1]);
      FracEvWTrigger = hEventsvsMult_EvwTrigger[i]->Integral(LimInf, LimSup)/hEventsvsMult[i]->Integral(LimInf, LimSup);
      if (i==0 || i ==2){
	hRatio->SetBinContent(hRatio->FindBin(dNdEta[i][m]), FracEvWTrigger);
	hRatio->SetBinError(hRatio->FindBin(dNdEta[i][m]), SetEfficiencyError(hEventsvsMult_EvwTrigger[i]->Integral(LimInf, LimSup), hEventsvsMult[i]->Integral(LimInf, LimSup)));
      }
      else {
	hRatio5TeV->SetBinContent(hRatio5TeV->FindBin(dNdEta[i][m]), FracEvWTrigger);
	hRatio5TeV->SetBinError(hRatio5TeV->FindBin(dNdEta[i][m]), SetEfficiencyError(hEventsvsMult_EvwTrigger[i]->Integral(LimInf, LimSup), hEventsvsMult[i]->Integral(LimInf, LimSup)));
      }
      cout << "m " << Nmolt[i][m] << " (dNdEta = " << dNdEta[i][m]<< ")  Fraction of events w trigger particle: " << FracEvWTrigger << endl;
      //      cout << hRatio->GetBinContent(hRatio->FindBin(dNdEta[i][m])) << " +- " << hRatio->GetBinError(hRatio->FindBin(dNdEta[i][m])) << endl;
    } 
  }


  //Set titles
  TString titleY = "Fraction of events w trigger particle";
  TString  titleX = titledNdeta;
  Float_t xOffset =1.2;
  Float_t yOffset =1.25;
  Float_t MarkerSize[3] ={2, 2, 3};
  Int_t MarkerType[3] = {20, 21, 33};

  gStyle->SetOptStat(0);
  Int_t Color = 628; 
  //  Int_t ColorDiff[numRegions][numColls] = {{634,628}, {418,829} , {601, 867}};
  //  Int_t Marker[numTypes] = {33, 21};
  //  Float_t Size[numTypes] = {2, 1.4};

  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  canvas->cd();
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
 
  StyleHisto(hRatio, 0, 0.6, 601, 33, titleX, titleY, "" , 1,0, 40, xOffset, yOffset, 2);
  StyleHisto(hRatio5TeV, 0, 0.6, 867, 33, titleX, titleY, "" , 1,0, 40, xOffset, yOffset, 2);
  hRatio->Draw("");
  hRatio5TeV->Draw("same");

  TLegend *legendEnergyBoxColor=new TLegend(0.16, 0.62, 0.39, 0.74);
  legendEnergyBoxColor->AddEntry(hRatio5TeV, "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
  legendEnergyBoxColor->AddEntry(hRatio, "pp, #sqrt{#it{s}} = 13 TeV", "pe");
  legendEnergyBoxColor->Draw("");

  cout << "Saving canvas" << endl;
  canvas->SaveAs("FractionOfEventsWithTrigger.pdf");
  
  TString fileName = "FractionOfEventsWithTrigger.root";
  cout << "Writing canvas on file " << fileName <<  endl;
  TFile* fileout = new TFile(fileName, "RECREATE");
  canvas->Write();
  fileout->Close();
  
}
