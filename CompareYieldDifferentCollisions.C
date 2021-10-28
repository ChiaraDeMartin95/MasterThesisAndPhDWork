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

  gStyle->SetOptStat(0);
  TString Region[3] = {"Jet", "Bulk", "All"};
  TString RegionBis[3] = {"Jet", "Bulk", "Inclusive"};
  Float_t Up[3] = {0.035, 0.4, 0.4};
  Float_t Low[3] = {0.01, 10e-4, 10e-4};
  if (type==1) {
    Low[0] = 10e-6; Up[0] = 0.003;
    Up[1] = 0.03;
    Up[2] = 0.03;
  }
  Int_t Color[3] = {628,418,601};
  TFile *file[NSystems];
  TString pathin[NSystems] = {"", ""}; // pp MB - pp HM 
  TString ParticleType[2]={"K0s", "Xi"};
  if (type==0){
    /*
    pathin[0] = "FinalOutput/DATA2016/SystematicAnalysis1617_hK0s_PtBinning1_K0s_Eta0.8_";
    pathin[0] += Region[TypeAnalysis];
    pathin[0] += "Data_PtMin3.0_IsEtaEff.root"; //pp MB
    */
    pathin[0] = "FinalOutput/DATA2016/PtSpectraBis_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_PtMin3.0_";
    pathin[0] += RegionBis[TypeAnalysis];
    pathin[0] += "_isNormCorrFullyComputed.root"; //pp MB 

    //pathin[0] += "Data_PtMin3.0_Eff18g4extra.root"; //pp MB
    //pathin[1] = "FinalOutput/DATA2016/SystematicAnalysis2016k_HM_hK0s_PtBinning1_K0s_Eta0.8_";
    pathin[1] = "FinalOutput/DATA2016/SystematicAnalysisAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_";
    pathin[1] += Region[TypeAnalysis];
    pathin[1] += "Data_PtMin3.0_IsEtaEff_MultBinning1.root"; 

  }
  else if (type==1){
    pathin[0] = "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_";
    pathin[0] += RegionBis[TypeAnalysis];
    if (TypeAnalysis==1) pathin[0] +="Blue";
    pathin[0] += "_isNormCorrFullyComputed_isEtaEff_isErrorAssumedPtCorr.root";
    /*
    pathin[0] = "FinalOutput/DATA2016/SystematicAnalysis161718Full_AOD234_hXi_";
    if (TypeAnalysis==0)     pathin[0]+= "OOJNoTriggerSmoothed_";
    pathin[0] += "Xi_Eta0.8_";
    pathin[0] += Region[TypeAnalysis];
    //    if (TypeAnalysis==1)    pathin[0] += "Blue";
    pathin[0] += "Data_PtMin3.0_IsEtaEff_isNormCorr.root"; //pp MB
    */
    pathin[1] = "FinalOutput/DATA2016/SystematicAnalysis161718_HM_hXi";
    if (TypeAnalysis==0)    pathin[1] += "_OOJAllMult";
    pathin[1] += "_Xi_Eta0.8_";
    pathin[1] += Region[TypeAnalysis];
    pathin[1] += "Data_PtMin3.0_IsEtaEff_MultBinning1.root"; 
  }

  TF1 * pol1 = new TF1("pol1", "pol1", 0, 45);
  pol1->SetLineColor(Color[TypeAnalysis]);
  pol1->SetLineWidth(0.3);
  TH1F *histoYield[NSystems];
  TH1F *histoYieldSist[NSystems];
  TH1F *histoYieldFinal[NSystems];
  TH1F *histoYieldSistFinal[NSystems];
  TString NameHisto[NSystems] = {"fHistYieldvsErrSoloStat", "fHistYieldvsErrSoloStat"};
  TString NameHistoSist[NSystems] = {"fHistYieldSist", "fHistYieldSist"};
  if (pathin[0].Index("PtSpectraBis")!=-1) NameHisto[0] = "fHistYieldStat";
  TString NameHistoFinal[NSystems] = {"fHistYield_ppMB", "fHistYield_ppHM"};
  TString NameHistoFinalSist[NSystems] = {"fHistYield_ppMB_Sist", "fHistYield_ppHM_Sist"};
  TH1F * histoYieldComparison = new TH1F ("histoYieldComparison", "histoYieldComparison",225, 0, 45 );
  TH1F * histoYieldSistComparison = new TH1F ("histoYieldSistComparison", "histoYieldSistComparison",225, 0, 45 );
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);

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
    if (!histoYield[i]) return;
    histoYieldFinal[i]= (TH1F*)    histoYield[i] ->Clone(NameHistoFinal[i]);
    if (i==0)    histoYieldSist[i] = (TH1F*) file[i]->Get(NameHistoSist[i]);
    else   {
      histoYieldSist[i] = (TH1F*) file[i]->Get(NameHisto[i]);
      for (Int_t b=1; b<= histoYield[i]->GetNbinsX(); b++){
	if (	histoYieldSist[i]->GetBinContent(b) ==0) continue;
	histoYieldSist[i]->SetBinError(b, histoYieldSist[0]->GetBinError( histoYieldSist[0]->FindBin(21.3))/ histoYieldSist[0]->GetBinContent( histoYieldSist[0]->FindBin(21.3))*histoYieldSist[i]->GetBinContent(b));
      }
    }
    //    if (!histoYieldSist[i]) return;
    histoYieldSistFinal[i]= (TH1F*)    histoYieldSist[i] ->Clone(NameHistoFinalSist[i]);
    for (Int_t b=1; b<= histoYieldFinal[i]->GetNbinsX(); b++){
      if (histoYieldFinal[i]->GetBinContent(b) == 0) continue;
      //      cout << histoYieldFinal[i]->GetBinContent(b) << endl;
    }

    for (Int_t b=1; b<= histoYieldFinal[i]->GetNbinsX(); b++){
      if (histoYieldFinal[i]->GetBinContent(b)==0)	continue;
      bin = histoYieldComparison ->FindBin(histoYieldFinal[i]-> GetBinCenter(b));
      //      if (histoYieldFinal[i]-> GetBinCenter(b) > 34 && TypeAnalysis==0) continue;
      histoYieldComparison -> SetBinContent( bin , histoYieldFinal[i]->GetBinContent(b));
      histoYieldComparison -> SetBinError( bin , histoYieldFinal[i]->GetBinError(b));
      histoYieldSistComparison -> SetBinContent( bin , histoYieldSistFinal[i]->GetBinContent(b));
      histoYieldSistComparison -> SetBinError( bin , histoYieldSistFinal[i]->GetBinError(b));
    }
    for (Int_t b=1; b<= histoYieldComparison->GetNbinsX(); b++){
      if ( histoYieldComparison->GetBinContent(b) !=0){
	cout << "\n"<<     histoYieldComparison->GetBinContent(b)<< " +- " << histoYieldComparison->GetBinError(b) << " (rel stat:" << histoYieldComparison->GetBinError(b)/histoYieldComparison->GetBinContent(b)<< " ) "  << endl;
	cout <<      histoYieldSistComparison->GetBinContent(b)<< " +- " << histoYieldSistComparison->GetBinError(b) << " (rel sist:" << histoYieldSistComparison->GetBinError(b)/histoYieldSistComparison->GetBinContent(b)<< " ) "  << endl;
      }
    }
    //    if (TypeAnalysis!=0)    histoYieldComparison->Fit(pol1);
  }

  TString NameFileout ="CompareYieldDifferentCollisions";
  if (pathin[1].Index("2016k_HM_hK0s")!=-1) NameFileout += "_HM16k";
  if (pathin[1].Index("MultBinning1")!=-1) NameFileout += "_HMMultBinning1";
  if (pathin[0].Index("18g4extra")!=-1) NameFileout += "_MBEff18g4extra_NoEtaEff";
  NameFileout += ParticleType[type] + Region[TypeAnalysis]+ ".root";

  TFile * fileout = new TFile (NameFileout, "RECREATE");
  StyleHisto(histoYieldComparison, Low[TypeAnalysis] , Up[TypeAnalysis], Color[TypeAnalysis], 33, "", "N/N_{Trigg} 1/#Delta#eta #Delta#phi", "Yield vs multiplicity" , 0,0, 45);
  StyleHisto(histoYieldSistComparison, Low[TypeAnalysis] , Up[TypeAnalysis], Color[TypeAnalysis], 33, "", "N/N_{Trigg} 1/#Delta#eta #Delta#phi", "Yield vs multiplicity" , 0,0, 45);
  histoYieldComparison->Draw("");
  histoYieldSistComparison->SetFillStyle(0);
  histoYieldSistComparison->Draw("same e2");
  canvas->Draw("");
  histoYieldComparison->Write("");
  histoYieldSistComparison->Write("");
  canvas->Write("");
  fileout->Close();
  canvas->SaveAs("CompareYieldDifferentCollisions"+ Region[TypeAnalysis]+".pdf");
  cout << "Starting from the files: " << endl;
  cout << pathin[0] << endl;
  cout << pathin[1] << endl;
  cout << " I have produced the file " << NameFileout << endl;
}
