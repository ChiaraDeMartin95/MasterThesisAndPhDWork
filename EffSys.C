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

void EffSys(){
  TString file = "AngularCorrelation2018d8_MCEff_15runs_6thtry";
  TString PathInBis;
  TFile *fileinbis;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=5;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;

  TString tipo[numtipo]={"kK0s", "bo"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString SysT[numSysTrigger]={"DCAz < 1","DCAz < 2","DCAz < 0.5"};
  TString SysV0[numSysV0]={"default", "cosTP> 0.997", "ctau <3 ", "YK0s < 0.5", "Lrejection 10 MeV","Lrejection 5 MeV", "V0dca< 0.3"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};

  Int_t Marker[numSysV0]={7,4,20,22,29, 35};
  Int_t Color[numSysV0]={2,3,4,6,7,9, 27};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};

  TCanvas *canvasSysT=new TCanvas ("canvasSysT", "canvasSysT", 1500, 800);
  TCanvas *canvasSysV0=new TCanvas ("canvasSysV0", "canvasSysV0", 1500, 800);
  TCanvas *canvasSysV0Bis=new TCanvas ("canvasSysV0Bis", "canvasSysV0Bis", 1500, 800);
  canvasSysT->Divide(2,1);
  canvasSysV0->Divide(3,2);
  canvasSysV0Bis->Divide(3,2);
  TCanvas *canvasSysPercT=new TCanvas ("canvasSysPercT", "canvasSysPercT", 1500, 800);
  TCanvas *canvasSysPercV0=new TCanvas ("canvasSysPercV0", "canvasSysPercV0", 1500, 800);
  TCanvas *canvasSysPercV0Bis=new TCanvas ("canvasSysPercV0Bis", "canvasSysPercV0Bis", 1500, 800);
  TCanvas *canvasSysPercV0Zoom=new TCanvas ("canvasSysPercV0Zoom", "canvasSysPercV0Zoom", 1500, 800);
  canvasSysPercT->Divide(2,1);
  canvasSysPercV0->Divide(3,2);
  canvasSysPercV0Zoom->Divide(3,2);
  canvasSysPercV0Bis->Divide(3,2);

  TH1D *histoEff;
  TH1D* HistContMolt;
  TH1D *histoEffMaster;
  TH1D* HistContMoltMaster;
  TH1D *histoEffPerc;
  TH1D* HistContMoltPerc;
  
  auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->SetHeader("Selezioni applicate");     
  auto legendV0 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legendV0->SetHeader("Selezioni applicate");     

  for(Int_t i=0; i < numSysTrigger; i++){
    PathInBis="histo/" + file + Form("_MC_Efficiency_SysT%i_SysV0%i.root",i, 0);
    fileinbis = new TFile(PathInBis);

    if(i==0) {
      histoEffMaster=(TH1D*)fileinbis->Get("HistoTriggerEfficiency");
      HistContMoltMaster=(TH1D*)fileinbis->Get("HistContTriggerMolt");
    }

    histoEff=(TH1D*)fileinbis->Get("HistoTriggerEfficiency");
    histoEffPerc=(TH1D*)histoEff->Clone("HistoTriggerEfficiencyPerc");
    canvasSysT->cd(1);
    histoEff->GetYaxis()->SetRangeUser(0.1,0.12);
    histoEff->SetMarkerStyle(Marker[i]);
    histoEff->SetLineColor(Color[i]);
    histoEff->SetMarkerColor(Color[i]);
    legend->AddEntry(histoEff,SysT[i],"pel");   
    histoEff->Draw("samee");
    if(i==numSysTrigger-1) legend->Draw();

    for(Int_t j=1; j<= histoEff->GetNbinsX(); j++){
      histoEffPerc->SetBinContent(j,(histoEff->GetBinContent(j)-histoEffMaster->GetBinContent(j))/histoEffMaster->GetBinContent(j));
    }
    canvasSysPercT->cd(1);
    histoEffPerc->GetYaxis()->SetRangeUser(-0.002,0.002);
    histoEffPerc->SetMarkerStyle(Marker[i]);
    histoEffPerc->SetLineColor(Color[i]);
    histoEffPerc->SetMarkerColor(Color[i]);
    histoEffPerc->Draw("samee");
    if(i==numSysTrigger-1) legend->Draw();

    HistContMolt=(TH1D*)fileinbis->Get("HistContTriggerMolt");
    HistContMoltPerc=(TH1D*)    HistContMolt->Clone("HistContTriggerMoltPerc");
    canvasSysT->cd(2);
    HistContMolt->GetYaxis()->SetRangeUser(0,0.03);
    HistContMolt->SetMarkerStyle(Marker[i]);
    HistContMolt->SetLineColor(Color[i]);
    HistContMolt->SetMarkerColor(Color[i]);
    HistContMolt->Draw("same");
    if(i==numSysTrigger-1) legend->Draw();

    for(Int_t j=1; j<= HistContMolt->GetNbinsX(); j++){
      HistContMoltPerc->SetBinContent(j,(HistContMolt->GetBinContent(j)-HistContMoltMaster->GetBinContent(j))/HistContMoltMaster->GetBinContent(j));
    }
    canvasSysPercT->cd(2);
    HistContMoltPerc->GetYaxis()->SetRangeUser(-0.2,0.2);
    HistContMoltPerc->SetMarkerStyle(Marker[i]);
    HistContMoltPerc->SetLineColor(Color[i]);
    HistContMoltPerc->SetMarkerColor(Color[i]);
    HistContMoltPerc->Draw("same");
    if(i==numSysTrigger-1) legend->Draw();

  }


  for(Int_t molt=0; molt< nummolt+1; molt++){
    for(Int_t i=0; i < numSysV0     ; i++){
      PathInBis="histo/" + file + Form("_MC_Efficiency_SysT%i_SysV0%i.root",0, i);
      fileinbis = new TFile(PathInBis);
      if(i==0) {
	histoEffMaster=(TH1D*)fileinbis->Get("fHistV0EfficiencyPtBins_" + Smolt[molt]);
      }
      histoEff=(TH1D*)fileinbis->Get("fHistV0EfficiencyPtBins_" + Smolt[molt]);
      histoEffPerc=(TH1D*)histoEff->Clone("HistoV0EfficiencyPerc");
      
      canvasSysV0->cd(molt+1);    
      histoEff->GetYaxis()->SetRangeUser(0,0.5);
      histoEff->SetMarkerStyle(Marker[i]);
      histoEff->SetLineColor(Color[i]);
      histoEff->SetMarkerColor(Color[i]);
      if (molt==0) legendV0->AddEntry(histoEff,SysV0[i],"pel");   
      histoEff->Draw("samee");
      if(i==numSysV0-1) legendV0->Draw();

      for(Int_t j=1; j<= histoEff->GetNbinsX(); j++){
	if(histoEffMaster->GetBinContent(j)!=0){
	histoEffPerc->SetBinContent(j,(histoEff->GetBinContent(j)-histoEffMaster->GetBinContent(j))/histoEffMaster->GetBinContent(j));
	}
	else 	histoEffPerc->SetBinContent(j,0);
      }
      canvasSysPercV0->cd(molt+1);  
      histoEffPerc->GetYaxis()->SetRangeUser(-0.6,0.8);
      histoEffPerc->SetMarkerStyle(Marker[i]);
      histoEffPerc->SetLineColor(Color[i]);
      histoEffPerc->SetMarkerColor(Color[i]);
      histoEffPerc->Draw("samee");
      if(i==numSysV0-1) legendV0->Draw();

      // if (i!=3){
      // canvasSysPercV0Zoom->cd(molt+1);  
      // histoEffPerc->GetYaxis()->SetRangeUser(-0.2,0.2);
      // histoEffPerc->SetMarkerStyle(Marker[i]);
      // histoEffPerc->SetLineColor(Color[i]);
      // histoEffPerc->SetMarkerColor(Color[i]);
      // histoEffPerc->Draw("samee");
      // if(i==numSysV0-1) legendV0->Draw();
      // }
    }
  }

  for(Int_t molt=0; molt< nummolt+1; molt++){
    for(Int_t i=0; i < 1     ; i++){
      PathInBis="histo/" + file + Form("_MC_Efficiency_SysT%i_SysV0%i.root",0, i);
      fileinbis = new TFile(PathInBis);
      if(i==0) {
	HistContMoltMaster=(TH1D*)fileinbis->Get("HistContV0PtBins_"+Smolt[molt]);
      }
      HistContMolt    =(TH1D*)fileinbis->Get("HistContV0PtBins_"+Smolt[molt]);
      HistContMoltPerc=(TH1D*)    HistContMolt->Clone("HistContV0MoltPerc");
      
      canvasSysV0Bis->cd(molt+1);
      HistContMolt->GetYaxis()->SetRangeUser(0,0.03);
      HistContMolt->SetMarkerStyle(Marker[i]);
      HistContMolt->SetLineColor(Color[i]);
      HistContMolt->SetMarkerColor(Color[i]);
      HistContMolt->Draw("same");
      //  if(i==numSysV0-1) legendV0->Draw();
      legendV0->Draw();

      for(Int_t j=1; j<= HistContMolt->GetNbinsX(); j++){
	if(HistContMoltMaster->GetBinContent(j)!=0){
	HistContMoltPerc->SetBinContent(j,(HistContMolt->GetBinContent(j)-HistContMoltMaster->GetBinContent(j))/HistContMoltMaster->GetBinContent(j));
	}
	else 	HistContMoltPerc->SetBinContent(j,0);
      }
      canvasSysPercV0Bis->cd(molt+1);     
      HistContMoltPerc->GetYaxis()->SetRangeUser(-0.2,0.2);
      HistContMoltPerc->SetMarkerStyle(Marker[i]);
      HistContMoltPerc->SetLineColor(Color[i]);
      HistContMoltPerc->SetMarkerColor(Color[i]);
      HistContMoltPerc->Draw("same");
      //if(i==numSysV0-1) legend->Draw();
      legendV0->Draw();
    }
  }
}

