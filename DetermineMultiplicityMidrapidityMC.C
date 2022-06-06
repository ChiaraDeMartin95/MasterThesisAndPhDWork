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

void DetermineMultiplicityMidrapidity(TString PathIn = "FinalOutput/AnalysisResultsPythiaRopes_MCTruth.root"){

  cout <<"Input file: " << PathIn << endl;
  TFile * filein = new TFile (PathIn, "");
  TDirectoryFile * dir = (TDirectoryFile*)filein->Get("MyTask_PtTrigMin3.0_PtTrigMax15.0");
  TList *list = (TList*)dir->Get("MyOutputContainer_hK0s_Task_K0s");
  TH2F *hMultiplicityV0vsMidrapidity = (TH2F*)list->FindObject("fHistMultForwardvsMidRap");

  const Int_t nummolt = 12;
  TString Smolt[nummolt+1]={"0-15", "15-30", "30-45", "45-54", "54-63", "63-72", "72-81", "81-90", "90-105", "105-120", "120-160", "160-300", "0-300"};
  Double_t Nmolt[nummolt+1]={0, 15, 30, 45, 54, 63, 72, 81, 90, 105, 120, 160, 300};

  TH1F*  hMultiplicityMidrapidity[nummolt+1];
  for (Int_t m=0; m<= nummolt; m++){
    if (m<nummolt)    hMultiplicityMidrapidity[m] = (TH1F*) hMultiplicityV0vsMidrapidity->ProjectionX("hMultiplicityMidrapidity" + Smolt[m], hMultiplicityV0vsMidrapidity->GetYaxis()->FindBin(Nmolt[m]+0.001), hMultiplicityV0vsMidrapidity->GetYaxis()->FindBin(Nmolt[m+1]-0.001));
    else   hMultiplicityMidrapidity[m] = (TH1F*) hMultiplicityV0vsMidrapidity->ProjectionX("hMultiplicityMidrapidity" + Smolt[m], 0, Nmolt[m]);
    cout << "\e[39mNumber of particles in V0 acceptance:\e[35m " << Smolt[m] << "\e[39m ---> dN/deta (|eta|<0.5): \e[35m" << hMultiplicityMidrapidity[m]->GetMean() <<" +- " <<hMultiplicityMidrapidity[m]->GetMeanError() << endl;
  }
}
