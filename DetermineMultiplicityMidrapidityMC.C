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
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>

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

void DetermineMultiplicityMidrapidity(TString PathIn = "FinalOutput/AnalysisResultsEPOSLHC_7BEvForhXi_MCTruth.root"/*"FinalOutput/AnalysisResultsEPOSLHC_3BEvForhK0s_MCTruth.root"/*"FinalOutput/AnalysisResultsPythiaMonash_MCTruth.root"/*"FinalOutput/AnalysisResultsMonashTuneForMult.root"/*"FinalOutput/AnalysisResultsPythiaRopes_MCTruth.root"*/, Bool_t isMultDistributionComp=0){

  cout <<"Input file: " << PathIn << endl;
  TFile * filein = new TFile (PathIn, "");
  TDirectoryFile * dir = (TDirectoryFile*)filein->Get("MyTask_PtTrigMin3.0_PtTrigMax15.0");
  TList *list = (TList*)dir->Get("MyOutputContainer_hK0s_Task_K0s");
  TH2F *hMultiplicityV0vsMidrapidity = (TH2F*)list->FindObject("fHistMultForwardvsMidRap");

  TString titledNdetaTrigg="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5, #it{p}_{T,trigg}>3 GeV/#it{c}}";

  TCanvas * canvasMult2D = new TCanvas("canvasMult2D", "canvasMult2D", 1300, 1000);
  TLegend *Legend=new TLegend(0.13,0.78,0.32,0.93);
  Legend->SetMargin(0);
  Legend->SetTextSize(0.03);
  Legend->AddEntry("", "#bf{This work}", "");
  Legend->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  Legend->AddEntry("", "EPOS LHC simulation", "");

  canvasMult2D->SetFillColor(0);
  canvasMult2D->SetTickx(1);
  canvasMult2D->SetTicky(1);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetOptStat(0);
  gPad->SetLogz();

  hMultiplicityV0vsMidrapidity->Rebin2D(1);
  hMultiplicityV0vsMidrapidity->GetXaxis()->SetRangeUser(0, 100);
  hMultiplicityV0vsMidrapidity->GetXaxis()->SetTitleOffset(1.5);
  hMultiplicityV0vsMidrapidity->GetYaxis()->SetTitleOffset(1.2);
  hMultiplicityV0vsMidrapidity->GetXaxis()->SetTitle(titledNdetaTrigg);
  hMultiplicityV0vsMidrapidity->GetYaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{V0A+V0C acceptance, #it{p}_{T,trigg}>3 GeV/#it{c}}");
  hMultiplicityV0vsMidrapidity->SetTitle("");
  hMultiplicityV0vsMidrapidity->Draw("colz");
  Legend->Draw("");

  canvasMult2D->SaveAs("MultV0vsMultMidrapidity.pdf");

  /*
  const Int_t nummolt = 12;
  TString Smolt[nummolt+1]={"0-15", "15-30", "30-45", "45-54", "54-63", "63-72", "72-81", "81-90", "90-105", "105-120", "120-160", "160-300", "0-300"};
  Double_t Nmolt[nummolt+1]={0, 15, 30, 45, 54, 63, 72, 81, 90, 105, 120, 160, 300};

  const Int_t nummolt = 9;
  TString Smolt[nummolt+1]={"0-13", "13-21", "21-28", "28-37", "37-47", "47-55", "55-64", "64-81", "81-300", "0-300"};
  Double_t Nmolt[nummolt+1]={0, 13, 21, 28, 37, 47, 55, 64, 81,300};
  */
  
  const Int_t nummolt = 10;
  TString Smolt[nummolt+1]={"0-30", "30-45", "45-54", "54-63", "63-72", "72-81", "81-90", "90-105", "105-120", "120-300", "0-300"};
  Double_t Nmolt[nummolt+1]={0,30, 45, 54, 63, 72, 81, 90, 105, 120, 300};
  
  /*
  const Int_t nummolt = 6;
  TString Smolt[nummolt+1]={"120-140", "140-160", "160-180", "180-200", "200-240","240-300", "0-300"};
  Double_t Nmolt[nummolt+1]={120, 140, 160, 180, 200, 240, 300};
  
  const Int_t nummolt = 5;
  TString Smolt[nummolt+1]={"120-140", "140-160", "160-180", "180-300", "0-300"};
  Double_t Nmolt[nummolt+1]={120, 140, 160, 180, 300};
*/
  TH1F*  hMultiplicityMidrapidity[nummolt+1];
  Float_t FractionOfTheTotal[nummolt+1]={0};
  Float_t TotalEntries = hMultiplicityV0vsMidrapidity->GetEntries();
  Float_t Percentile = 100;
  for (Int_t m=0; m<= nummolt; m++){
    if (m<nummolt)    hMultiplicityMidrapidity[m] = (TH1F*) hMultiplicityV0vsMidrapidity->ProjectionX("hMultiplicityMidrapidity" + Smolt[m], hMultiplicityV0vsMidrapidity->GetYaxis()->FindBin(Nmolt[m]+0.001), hMultiplicityV0vsMidrapidity->GetYaxis()->FindBin(Nmolt[m+1]-0.001));
    else   hMultiplicityMidrapidity[m] = (TH1F*) hMultiplicityV0vsMidrapidity->ProjectionX("hMultiplicityMidrapidity" + Smolt[m], 0, Nmolt[m]);
    FractionOfTheTotal[m] = hMultiplicityMidrapidity[m]->GetEntries()/TotalEntries;
    cout << "\e[39mNumber of particles in V0 acceptance:\e[35m " << Smolt[m] << "\e[39m ---> dN/deta (|eta|<0.5): \e[35m" << hMultiplicityMidrapidity[m]->GetMean() <<" +- " <<hMultiplicityMidrapidity[m]->GetMeanError() << endl;
    cout << "Percentile: " << Percentile << "-";
    Percentile -= 100*FractionOfTheTotal[m];
    cout << Percentile << endl;
  }

  //Ratio between Mult distribution in events NT>0 and all INEL>0 events (with pt,trig > 0)
  if (isMultDistributionComp){
    TCanvas * canvasRatio = new TCanvas("canvasRatio", "canvasRatio", 800, 500);
    canvasRatio->Divide(2,1);
    TList *listINEL = (TList*)dir->Get("MyOutputContainer_hK0s_Task_K0s_Inclusive");
    if (!listINEL) return;
    TH2F *hMultiplicityV0vsMidrapidityINEL = (TH2F*)listINEL->FindObject("fHistMultForwardvsMidRap");
    if (!hMultiplicityV0vsMidrapidityINEL) return;
    TH1F *hMultiplicityV0INEL = (TH1F*)hMultiplicityV0vsMidrapidityINEL->ProjectionY("hMultiplicityV0INEL");
    TH1F *hMultiplicityV0 = (TH1F*)hMultiplicityV0vsMidrapidity->ProjectionY("hMultiplicityV0");
    hMultiplicityV0INEL->Scale(1./hMultiplicityV0INEL->GetEntries());
    hMultiplicityV0INEL->SetLineColor(kBlack);
    hMultiplicityV0->Scale(1./hMultiplicityV0->GetEntries());
    hMultiplicityV0->SetLineColor(kGray);
    TH1F *hMultiplicityV0Ratio = (TH1F*)hMultiplicityV0->Clone("hMultiplicityV0Ratio");
    hMultiplicityV0Ratio->Divide(hMultiplicityV0INEL);

    canvasRatio->cd(1);
    hMultiplicityV0INEL->Draw("");
    hMultiplicityV0->Draw("same");

    canvasRatio->cd(2);
    hMultiplicityV0Ratio->Draw("");

  }


}
