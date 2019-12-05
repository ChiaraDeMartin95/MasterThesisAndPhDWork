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

void EffCfrPtMin(Bool_t ishhCorr =1,Int_t type=0,  Int_t sysTrigger=0, Int_t sysV0=0, TString data="2018f1_extra", TString year0="2016"){
  //onlyTriggerWithHighestPt
  const Int_t numPtTMin=10;
  const Int_t nummolt=5;
  Float_t PtTrigMinArray[numPtTMin]={3., 4., 5, 6., 7., 8., 9., 10., 11. , 12.};
  Int_t Color[numPtTMin]={1,401, 801, 628, 909, 881, 860, 868, 841, 418};
  Int_t Marker[numPtTMin]={7,20,20,22,29,25, 7, 20, 22, 29};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString nameRes_1D[1][6];
  TString PathIn;  
  TFile *filein;
  TH1F *HistoTriggerEfficiency[numPtTMin];
  TH1F *HistContTrigger[numPtTMin];
  TH1F *HistResolution[numPtTMin];
  TH1F* HistV0EfficiencyPtBins[numPtTMin][nummolt+1];
  TH1F* fHistTriggerEfficiencyEta[numPtTMin][nummolt+1];
  TH1F* HistContV0PtBins[numPtTMin][nummolt+1];
  TString TorV[2]={"Trigger", "V0"};
  TString Var[3]={"Pt", "Phi", "Eta"};

  TCanvas * canvasTrigger=new TCanvas ("canvasTrigger", "canvasTrigger",1500,800);
  //TCanvas * canvasTriggerMolt=new TCanvas ("canvasTriggerMolt", "canvasTriggerMolt",1500,800);
  TCanvas * canvasV0Eff=new TCanvas ("canvasV0Eff", "canvasV0Eff",1500,800);
  TCanvas * canvasV0Cont=new TCanvas ("canvasV0Cont", "canvasV0Cont",1500,800);
  TCanvas * canvasResol=new TCanvas ("canvasResol", "canvasResol",1500,800);
  canvasTrigger->Divide(2,1);
  //  canvasTriggerMolt->Divide(3,2);
  canvasV0Eff->Divide(3,2);
  canvasV0Cont->Divide(3,2);
  canvasResol->Divide(3,2);
  auto legend = new TLegend(0.6, 0.1, 0.9, 0.4);
  legend->SetHeader("Minimum Pt Trigger");     


  for(Int_t pT=0;  pT<numPtTMin; pT++){
    PathIn="FinalOutput/DATA" + year0 + "/Efficiency/Efficiency"+data+ Form("_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMinArray[pT]);
    if (ishhCorr)     PathIn="FinalOutput/DATA" + year0 + "/Efficiency/Efficiency"+data+ Form("_hhCorr_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMinArray[pT]);
    filein = new TFile(PathIn, "");
    HistoTriggerEfficiency[pT]=(TH1F*)filein->Get("HistoTriggerEfficiency");
    HistContTrigger[pT]=(TH1F*)filein->Get("HistContTriggerMolt");
 cout << "pT trigger minimo " << PtTrigMinArray[pT] << "  " << HistoTriggerEfficiency[pT]->GetBinError(1) << " bin content " << HistoTriggerEfficiency[pT]->GetBinContent(1)<< " relativo: " << HistoTriggerEfficiency[pT]->GetBinError(1)/HistoTriggerEfficiency[pT]->GetBinContent(1)<< endl;
  canvasTrigger->cd(1);
  gPad->SetLeftMargin(0.15);
  HistoTriggerEfficiency[pT]->GetYaxis()->SetRangeUser(0,1);
    HistoTriggerEfficiency[pT]->SetMarkerStyle(Marker[pT]);
  HistoTriggerEfficiency[pT]->GetYaxis()->SetTitle("#epsilon_{Trigg}");
  HistoTriggerEfficiency[pT]->GetYaxis()->SetTitleOffset(1);
  HistoTriggerEfficiency[pT]->GetYaxis()->SetTitleSize(0.05);
  HistoTriggerEfficiency[pT]->SetStats(0);
  HistoTriggerEfficiency[pT]->SetLineColor(Color[pT]);
  HistoTriggerEfficiency[pT]->SetMarkerColor(Color[pT]);
  HistoTriggerEfficiency[pT]->Draw("samep");
  legend->AddEntry(  HistoTriggerEfficiency[pT],Form("p_{T, Min}>%.1f",PtTrigMinArray[pT]), "pel");
  if(pT ==numPtTMin-1)     legend->Draw();

  canvasTrigger->cd(2);
    gPad->SetLeftMargin(0.15);
HistContTrigger[pT]->GetYaxis()->SetRangeUser(0,0.02);
HistContTrigger[pT]->GetYaxis()->SetTitle("C_{Trigg}");
HistContTrigger[pT]->GetYaxis()->SetTitleOffset(1.5);
HistContTrigger[pT]->GetYaxis()->SetTitleSize(0.05);
HistContTrigger[pT]->SetStats(0);
HistContTrigger[pT]->SetLineColor(Color[pT]);
HistContTrigger[pT]->SetMarkerColor(Color[pT]);
HistContTrigger[pT]->Draw("same");

  if(pT ==numPtTMin-1)     legend->Draw();

  for(Int_t mult=0; mult< nummolt+1; mult++){
    HistV0EfficiencyPtBins[pT][mult]=(TH1F*)filein->Get("fHistV0EfficiencyPtBins_"+Smolt[mult]);
    canvasV0Eff->cd(mult+1);
    HistV0EfficiencyPtBins[pT][mult]->GetYaxis()->SetRangeUser(0,0.3);
    if (ishhCorr)     HistV0EfficiencyPtBins[pT][mult]->GetYaxis()->SetRangeUser(0,1);
    HistV0EfficiencyPtBins[pT][mult]->GetYaxis()->SetTitle("#epsilon_{Assoc}");
    HistV0EfficiencyPtBins[pT][mult]->GetYaxis()->SetTitleOffset(1);
    HistV0EfficiencyPtBins[pT][mult]->GetYaxis()->SetTitleSize(0.05);
    HistV0EfficiencyPtBins[pT][mult]->SetStats(0);
    HistV0EfficiencyPtBins[pT][mult]->SetLineColor(Color[pT]);
    HistV0EfficiencyPtBins[pT][mult]->SetMarkerColor(Color[pT]);
    HistV0EfficiencyPtBins[pT][mult]->SetMarkerStyle(Marker[pT]);
    HistV0EfficiencyPtBins[pT][mult]->Draw("same");
    if(pT ==numPtTMin-1)     legend->Draw();

    HistContV0PtBins[pT][mult]=(TH1F*)filein->Get("HistContV0PtBins_"+Smolt[mult]);
    canvasV0Cont->cd(mult+1);
    gPad->SetLeftMargin(0.15);
    HistContV0PtBins[pT][mult]->GetYaxis()->SetRangeUser(0,0.01);
    if (ishhCorr)     HistContV0PtBins[pT][mult]->GetYaxis()->SetRangeUser(0,0.08);
    HistContV0PtBins[pT][mult]->GetYaxis()->SetTitle("C_{Assoc}");
    HistContV0PtBins[pT][mult]->GetYaxis()->SetTitleOffset(1.5);
    HistContV0PtBins[pT][mult]->GetYaxis()->SetTitleSize(0.05);
    HistContV0PtBins[pT][mult]->SetStats(0);
    HistContV0PtBins[pT][mult]->SetLineColor(Color[pT]);
    HistContV0PtBins[pT][mult]->SetMarkerColor(Color[pT]);
    HistContV0PtBins[pT][mult]->SetMarkerStyle(Marker[pT]);
    HistContV0PtBins[pT][mult]->Draw("same");
    if(pT ==numPtTMin-1)     legend->Draw();

    /*
    fHistTriggerEfficiencyEta[pT][mult]=(TH1F*)filein->Get("fHistTriggerEfficiencyEta_"+Smolt[mult]);
    canvasTriggerMolt->cd(mult+1);
    fHistTriggerEfficiencyEta[pT][mult]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyEta[pT][mult]->GetYaxis()->SetTitle("#C_{Assoc}");
    fHistTriggerEfficiencyEta[pT][mult]->GetYaxis()->SetTitleOffset(1);
    fHistTriggerEfficiencyEta[pT][mult]->GetYaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyEta[pT][mult]->SetStats(0);
    fHistTriggerEfficiencyEta[pT][mult]->SetLineColor(Color[pT]);
    fHistTriggerEfficiencyEta[pT][mult]->SetMarkerColor(Color[pT]);
    fHistTriggerEfficiencyEta[pT][mult]->SetMarkerStyle(Marker[pT]);
    fHistTriggerEfficiencyEta[pT][mult]->Draw("same");
    if(pT ==numPtTMin-1)     legend->Draw();
    */
  }//endl loop on mult

  Int_t i=0;
    for(Int_t m=0; m< 3; m++){
    for(Int_t t=0; t< 2; t++){
      if (t==0) i=m;
      if (t==1) i=m+3;
      nameRes_1D[0][i]="fHistResolution" + TorV[t] + Var[m] + Smolt[0];
    }
    }
    //    cout << "Minimum pT of trigger particle is " <<PtTrigMinArray[pT]<< endl;
  for(Int_t i=0; i<6; i++){
    HistResolution[pT]=(TH1F*)filein->Get(nameRes_1D[0][i]);
    // if (i==0) cout << "Resolution of trigger particle (pT in GeV, Eta, Phi)" << endl;
    // if (i==3) cout << "Resolution of associated particle (pT in GeV, Eta, Phi)" << endl;
    // cout << " Molt class: " << Smolt[0] << ", Mean: "<<      HistResolution[pT]->GetMean()<< " RMS: "<< HistResolution[pT]->GetRMS()<< endl ;
  canvasResol->cd(i+1);
    gPad->SetLeftMargin(0.15);
HistResolution[pT]->GetYaxis()->SetRangeUser(0,0.3);
 if (i==0){
   HistResolution[pT]->GetYaxis()->SetRangeUser(0,0.3);
   HistResolution[pT]->Rebin(10);
 }
HistResolution[pT]->GetYaxis()->SetTitle("Resolution");
HistResolution[pT]->GetYaxis()->SetTitleOffset(1.5);
HistResolution[pT]->GetYaxis()->SetTitleSize(0.05);
HistResolution[pT]->SetStats(0);
HistResolution[pT]->SetLineColor(Color[pT]);
HistResolution[pT]->SetMarkerColor(Color[pT]);
HistResolution[pT]->Draw("same");
  if(pT ==numPtTMin-1)     legend->Draw();
  }//end loop on type of resolution histo

  }//end loop on minimum pT of trigger particles



   
  for(Int_t i=0; i<6; i++){
      if (i==0) cout << "Resolution of trigger particle pT in GeV" << endl;
      if (i==1) cout << "Resolution of trigger particle Eta" << endl;
      if (i==2) cout << "Resolution of trigger particle Phi" << endl;
      if (i==3) cout << "Resolution of associated particle (pT in GeV)" << endl;
      if (i==4) cout << "Resolution of associated particle (Eta)" << endl;
      if (i==5) cout << "Resolution of associated particle (Phi)" << endl;
    for (Int_t pT=0; pT< numPtTMin; pT++){      
      PathIn="FinalOutput/DATA" + year0 + "/Efficiency/Efficiency"+data+ Form("_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMinArray[pT]);
      if (ishhCorr)     PathIn="FinalOutput/DATA" + year0 + "/Efficiency/Efficiency"+data+ Form("_hhCorr_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMinArray[pT]);
      filein = new TFile(PathIn, "");
      HistResolution[pT]=(TH1F*)filein->Get(nameRes_1D[0][i]);
      cout << " Molt class: " << Smolt[0] << "for minimum pT trigger part " << PtTrigMinArray[pT] << ", Mean: "<<      HistResolution[pT]->GetMean()<< " RMS: "<< HistResolution[pT]->GetRMS()<< endl;
    }
  }
}
