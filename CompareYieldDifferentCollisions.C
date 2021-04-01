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

void CompareYieldDifferentCollisions(Int_t TypeAnalysis=0, Int_t NSystems =2, Int_t type=0){

  TString Region[3] = {"Jet", "Bulk", "All"};
  Float_t Up[3] = {0.04, 0.4, 0.4};
  Float_t Low[3] = {0.02, 0, 0};
  Int_t Color[3] = {628,418,601};
  TFile *file[NSystems];
  TString pathin[NSystems] = {"", ""}; // pp MB - pp HM 
  pathin[0] = "FinalOutput/DATA2016/SystematicAnalysis1617_hK0s_PtBinning1_K0s_Eta0.8_";
  pathin[0] += Region[TypeAnalysis];
  pathin[0] += "Data_PtMin3.0_IsEtaEff.root"; //pp MB
  //pathin[0] += "Data_PtMin3.0_Eff18g4extra.root"; //pp MB
  //pathin[1] = "FinalOutput/DATA2016/SystematicAnalysis2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_";
  pathin[1] = "FinalOutput/DATA2016/SystematicAnalysisAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_";
  pathin[1] += Region[TypeAnalysis];
  pathin[1] += "Data_PtMin3.0_IsEtaEff.root"; 

  TF1 * pol1 = new TF1("pol1", "pol1", 0, 45);
  TH1F *histoYield[NSystems];
  TH1F *histoYieldFinal[NSystems];
  TString NameHisto[NSystems] = {"fHistYieldvsErrSoloStat", "fHistYieldvsErrSoloStat"};
  TString NameHistoFinal[NSystems] = {"fHistYield_ppMB", "fHistYield_ppHM"};
  TH1F * histoYieldComparison = new TH1F ("histoYieldComparison", "histoYieldComparison",225, 0, 45 );
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);

  Int_t bin=0;
  for (Int_t b=1; b<= histoYieldComparison->GetNbinsX(); b++){
    histoYieldComparison -> SetBinContent(b, 0);
    histoYieldComparison -> SetBinError(b, 0);
  }


  for (Int_t i=0; i< NSystems; i++){
    file[i] = new TFile(pathin[i], "");
    histoYield[i] = (TH1F*) file[i]->Get(NameHisto[i]);
    histoYieldFinal[i]= (TH1F*)    histoYield[i] ->Clone(NameHistoFinal[i]);
    for (Int_t b=1; b<= histoYieldFinal[i]->GetNbinsX(); b++){
      if (histoYieldFinal[i]->GetBinContent(b) == 0) continue;
      cout << histoYieldFinal[i]->GetBinContent(b) << endl;
    }

    for (Int_t b=1; b<= histoYieldFinal[i]->GetNbinsX(); b++){
      if (histoYieldFinal[i]->GetBinContent(b)==0)	continue;
      bin = histoYieldComparison ->FindBin(histoYieldFinal[i]-> GetBinCenter(b));
      //      if (histoYieldFinal[i]-> GetBinCenter(b) > 34 && TypeAnalysis==0) continue;
      histoYieldComparison -> SetBinContent( bin , histoYieldFinal[i]->GetBinContent(b));
      histoYieldComparison -> SetBinError( bin , histoYieldFinal[i]->GetBinError(b));
    }
    for (Int_t b=1; b<= histoYieldComparison->GetNbinsX(); b++){
      if ( histoYieldComparison->GetBinContent(b) !=0){
      cout <<      histoYieldComparison->GetBinContent(b) << endl;
      }
    }
    //    if (TypeAnalysis!=0)    histoYieldComparison->Fit(pol1);
  }

  TString NameFileout ="CompareYieldDifferentCollisions";
  if (pathin[1].Index("2016k_HM_hK0s")!=-1) NameFileout += "_HM16k";
  if (pathin[0].Index("18g4extra")!=-1) NameFileout += "_MBEff18g4extra_NoEtaEff";
  NameFileout += Region[TypeAnalysis]+ ".root";

  TFile * fileout = new TFile (NameFileout, "RECREATE");
  StyleHisto(histoYieldComparison, Low[TypeAnalysis] , Up[TypeAnalysis], Color[TypeAnalysis], 33, "", "N/N_{Trigg} 1/#Delta#eta #Delta#phi", "Yield vs multiplicity" , 0,0, 45);
  histoYieldComparison->Draw("");
  canvas->Draw("");
  histoYieldComparison->Write("");
  canvas->Write("");
  fileout->Close();
  canvas->SaveAs("CompareYieldDifferentCollisions"+ Region[TypeAnalysis]+".pdf");
  cout << " I have produced the file " << NameFileout << endl;
}
