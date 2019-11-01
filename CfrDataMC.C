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

void CfrDataMC(TString year0="2016"){

  TString PathIn;
 TString PathOut;
  TFile *filein;
 

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;

  Int_t Color[6]={2,2,418,418, 4, 4};
  Int_t ColorBis[6]={2,418,4,2, 418, 4};

  TString Sys[6]={"Jet DATA", "Jet MC", "OJ DATA", "OJ MC", "J+OJ DATA", "J+OJ MC"};
  TString SysErr[2]={"syst. DATA", "syst. MC"};
  TString StatErr[2]={"stat. DATA", "stat. MC"};
  TString StatErrBis[3]={"stat. in-jet", "stat. out-of-jet", "stat. inclusive"};  
  TString SysErrBis[3]={"sist. in-jet", "sist. out-of-jet", "sist. inclusive"};
  TString SysErrMB[2]={"syst. DATA 0-100 %", "syst. MC 0-100 %"};
  TString StatErrMB[2]={"stat. DATA 0-100 %", "stat. MC 0-100 %"};
    TString Region[3]={"Jet", "Bulk", "All"};
  TString DATAorMC[2]={"DATA", "MC"};

  TString tipo[numtipo]={"kK0s", "bo"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  TString SmoltBis[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-1.5", "1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};
  TString SPeriod[numPeriod]={"2016", "2017"};//, "2018"};
  TString SRun[numPeriod]={"2016k", "2017h"};//, "2018m"};
   // TString SPeriod[numPeriod]={"2017", "2018", "2016"};
   // TString SRun[numPeriod]={"2017h", "2018m", "2016k"};

  //20,21,33
  Int_t Marker[6]={33,4,33,4, 33, 4};
  Int_t MarkerOrigin[6]={20,21,33,20,21,33};
  Int_t MarkerBis[2]={23,24};
  Int_t MarkerTris[2]={31, 33};
  Int_t MarkerMBStat[6]={20,4,20,4, 20, 4};
  Int_t MarkerMBSist[6]={21,25,21,25, 21, 25};
  // Int_t Color[numSysV0]={2,3,4,6,7,9, 27};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};
  TF1* pol0[nummolt+1][numPtV0][numPeriod];
  TF1* pol0bis[nummolt+1][numPtV0][numPeriod];

  TH1D*  fHistYieldvsErrSoloStat[6];
  TH1D*  fHistYieldvsErrSoloSist[6];
  TH1D*  fHistYieldvsErrSoloStatMB[6];
  TH1D*  fHistYieldvsErrSoloSistMB[6];
  TH1D*  fHistYieldvsErrTotal[6];
  TH1D*  fHistYieldvsErrSoloStatRatio[6];
  TH1D*  fHistYieldvsErrSoloSistRatio[6];
  TH1D*  fHistYieldvsErrSoloStatRatioMB[6];
  TH1D*  fHistYieldvsErrSoloSistRatioMB[6];
  TH1D*  fHistYieldvsErrTotalRatio[6];
  TH1D*  fHistYieldvsErrErrSoloStat[6];
  TH1D*  fHistYieldvsErrErrSoloSist[6];
  TH1D*  fHistYieldvsErrErrSoloStatMB[6];
  TH1D*  fHistYieldvsErrErrSoloSistMB[6];
  TH1D*  fHistYieldvsErrErrTotal[6];

  TLegend* legend[3];
  TLegend* legendtype[3];
  TLegend* legendbis[3];
  Int_t msize=2;

  auto legenderr=new TLegend(0.6,0.6, 0.9, 0.9);
  auto legenderrbis=new TLegend(0.6,0.6, 0.9, 0.9);

 
   // TF1 *power[4]
   //   power[j] = new TF1("power"+Sys[j], "[1]*x**[0]",1,30);
   // power[j]->SetLineColor(Color[j]);
   // fHistYieldvsErr[j]->Fit(power[j],"R+");

  TCanvas *canvas[3];
  Int_t col=0;
  for(Int_t l=0; l<3; l++){
  Int_t count=0;
  canvas[l] = new TCanvas(Form("canvas%i", l), Region[l], 1300, 1000);
  legend[l]=new TLegend(0.7,0.1, 0.9, 0.3);
  legendtype[l]=new TLegend(0.7,0.35, 0.9, 0.45);
  legendbis[l]=new TLegend(0.7,0.1, 0.9, 0.3);
  cout << l << endl;
  for(Int_t j=0; j<2; j++){
    if(l==0){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJet.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJetMC.root";
    col=j;
    }
    if(l==1){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulk.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulkMC.root";
    col=j+2;
    }
    if (l==2){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAll.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAllMC.root";
    col=j+4;
    }
    filein=new TFile(PathIn, "");
    fHistYieldvsErrSoloSistMB[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSistMB");
    fHistYieldvsErrSoloStatMB[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStatMB");
    fHistYieldvsErrSoloStatRatio[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStatRatio");
    fHistYieldvsErrSoloStatRatio[j]->SetTitle("");
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerColor(Color[col]);
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerStyle(Marker[j]);
    fHistYieldvsErrSoloStatRatio[j]->SetLineColor(Color[col]);
    fHistYieldvsErrSoloStatRatio[j]->GetYaxis()->SetRangeUser(0.1,2);
    fHistYieldvsErrSoloStatRatio[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerSize(msize);
    if (j==0)    fHistYieldvsErrSoloStatRatio[j]->SetMarkerSize(3);
    fHistYieldvsErrSoloStatRatio[j]->Draw("samee");
       // fHistYieldvsErrTotal[j]=(TH1D*)filein->Get("fHistYieldvsErrTot");
       // fHistYieldvsErrTotal[j]->SetMarkerColor(Color[j]);
       // fHistYieldvsErrTotal[j]->SetMarkerStyle(Marker[j]);
       // fHistYieldvsErrTotal[j]->SetLineColor(Color[j]);
       // fHistYieldvsErrTotal[j]->Draw("same");
    // fHistYieldvsErrSoloStatRatioMB[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStatRatioMB");
    // fHistYieldvsErrSoloStatRatioMB[j]->SetMarkerColor(1);
    // fHistYieldvsErrSoloStatRatioMB[j]->SetMarkerStyle(MarkerMBSist[j]);
    // fHistYieldvsErrSoloStatRatioMB[j]->SetLineColor(1);
    // fHistYieldvsErrSoloStatRatioMB[j]->Draw("samee");
      
    fHistYieldvsErrSoloSistRatio[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSistRatio");
    fHistYieldvsErrSoloSistRatio[j]->SetTitle("");
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerColor(Color[col]);
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerStyle(Marker[j]);
    fHistYieldvsErrSoloSistRatio[j]->SetLineColor(Color[col]);
    fHistYieldvsErrSoloSistRatio[j]->SetFillStyle(0);
    fHistYieldvsErrSoloSistRatio[j]->GetYaxis()->SetRangeUser(0.1,2);
    fHistYieldvsErrSoloSistRatio[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerSize(msize);
    if (j==0)    fHistYieldvsErrSoloSistRatio[j]->SetMarkerSize(3);
     fHistYieldvsErrSoloSistRatio[j]->Draw("samee2");
      
    // fHistYieldvsErrSoloSistRatioMB[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSistRatioMB");
    // fHistYieldvsErrSoloSistRatioMB[j]->SetMarkerColor(1);
    // fHistYieldvsErrSoloSistRatioMB[j]->SetMarkerStyle(MarkerMBSist[j]);
    // fHistYieldvsErrSoloSistRatioMB[j]->SetLineColor(1);
    // fHistYieldvsErrSoloSistRatioMB[j]->SetFillStyle(0);
    // fHistYieldvsErrSoloSistRatioMB[j]->Draw("samee2");
        
    legend[l]->AddEntry(    fHistYieldvsErrSoloStatRatio[j], Sys[col], "pl");
    // if (j==1) legend[l]->AddEntry(    fHistYieldvsErrSoloStatRatio[0], Sys[j+2*l] + " 0-100 %", "pl");
    // if (j==1) legend[l]->AddEntry(    fHistYieldvsErrSoloStatRatio[j], Sys[col] + " 0-100 %", "pl");
     if(j==1)    legendtype[l]->AddEntry(fHistYieldvsErrSoloStatMB[j], "stat. ", "el");
    if(j==1)    legendtype[l]->AddEntry(fHistYieldvsErrSoloStatMB[j], "syst. ", "f");
    if(j==1)legend[l]->Draw();
    if(j==1)legendtype[l]->Draw();
    // fHistYieldvsErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStat");
    // fHistYieldvsErrSoloSist[j]->SetMarkerColor(Color[j]);
    // fHistYieldvsErrSoloSist[j]->Draw("samee");
  }
  canvas[l]->SaveAs("FinalOutput/DATA"+year0+"/Yield"+Region[l]+".pdf");
  }

  /*
 TCanvas *canvasratio[3];
  for(Int_t l=0; l<3; l++){
  Int_t count=0;
  canvasratio[l] = new TCanvas(Form("canvasratio%i", l), "ratio_"+Region[l], 1300, 1000);
  cout << l << endl;
  for(Int_t j=0; j<2; j++){
    if(l==0){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJet.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJetMC.root";
    col=j;
    }
    if(l==1){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulk.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulkMC.root";
    col=j+2;
    }
    if (l==2){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAll.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAllMC.root";
    col=j+4;
    }
    filein=new TFile(PathIn, "");

    fHistYieldvsErrSoloStatRatio[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStatRatio");
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerColor(Color[col]);
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerStyle(Marker[j]);
    fHistYieldvsErrSoloStatRatio[j]->SetLineColor(Color[col]);
    fHistYieldvsErrSoloStatRatio[j]->Draw("samee");
       // fHistYieldvsErrTotal[j]=(TH1D*)filein->Get("fHistYieldvsErrTot");
       // fHistYieldvsErrTotal[j]->SetMarkerColor(Color[j]);
       // fHistYieldvsErrTotal[j]->SetMarkerStyle(Marker[j]);
       // fHistYieldvsErrTotal[j]->SetLineColor(Color[j]);
       // fHistYieldvsErrTotal[j]->Draw("same");
        
    fHistYieldvsErrSoloSistRatio[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSistRatio");
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerColor(Color[col]);
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerStyle(Marker[j]);
    fHistYieldvsErrSoloSistRatio[j]->SetLineColor(Color[col]);
    fHistYieldvsErrSoloSistRatio[j]->SetFillStyle(0);
    fHistYieldvsErrSoloSistRatio[j]->Draw("samee2");
        

    if(j==0)legend[l]->Draw();
    // fHistYieldvsErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStat");
    // fHistYieldvsErrSoloSist[j]->SetMarkerColor(Color[j]);
    // fHistYieldvsErrSoloSist[j]->Draw("samee");
  }
  }
  */

 TCanvas *canvasratiobis[3];
  for(Int_t l=0; l<2; l++){
  Int_t count=0;
  canvasratiobis[l] = new TCanvas(Form("canvasratiobis%i", l), "ratio_"+DATAorMC[l], 1300, 1000);
  cout << l << endl;
  for(Int_t j=0; j<3; j++){
    if(l==0){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJet.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulk.root";
    if (j==2)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAll.root";
    col=j;
    }
    if(l==1){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJetMC.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulkMC.root";
    if (j==2)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAllMC.root";
    col=j+3;
    }
    filein=new TFile(PathIn, "");

    fHistYieldvsErrSoloStatRatio[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStatRatio");
    fHistYieldvsErrSoloStatRatio[j]->SetTitle("");
    fHistYieldvsErrSoloStatRatio[j]->GetYaxis()->SetRangeUser(0.03, 2.2);
    fHistYieldvsErrSoloStatRatio[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerColor(Color[2*j]);
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerStyle(MarkerOrigin[col]);
    fHistYieldvsErrSoloStatRatio[j]->SetMarkerSize(msize);
    if (j==2)    fHistYieldvsErrSoloStatRatio[j]->SetMarkerSize(3);
    fHistYieldvsErrSoloStatRatio[j]->SetLineColor(Color[2*j]);
    fHistYieldvsErrSoloStatRatio[j]->Draw("samee");
       // fHistYieldvsErrTotal[j]=(TH1D*)filein->Get("fHistYieldvsErrTot");
       // fHistYieldvsErrTotal[j]->SetMarkerColor(Color[j]);
       // fHistYieldvsErrTotal[j]->SetMarkerStyle(Marker[j]);
       // fHistYieldvsErrTotal[j]->SetLineColor(Color[j]);
       // fHistYieldvsErrTotal[j]->Draw("same");
        
    fHistYieldvsErrSoloSistRatio[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSistRatio");
    fHistYieldvsErrSoloSistRatio[j]->SetTitle("");
    fHistYieldvsErrSoloSistRatio[j]->GetYaxis()->SetRangeUser(0.03, 2.2);
    fHistYieldvsErrSoloSistRatio[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerColor(Color[2*j]);
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerStyle(MarkerOrigin[col]);
    fHistYieldvsErrSoloSistRatio[j]->SetLineColor(Color[2*j]);
    fHistYieldvsErrSoloSistRatio[j]->SetFillStyle(0);
    fHistYieldvsErrSoloSistRatio[j]->SetMarkerSize(msize);
    if (j==2)    fHistYieldvsErrSoloSistRatio[j]->SetMarkerSize(3);
    fHistYieldvsErrSoloSistRatio[j]->Draw("samee2");
        
    if(l==0)    legendbis[l]->AddEntry(    fHistYieldvsErrSoloStatRatio[j], Sys[2*j], "pl");
    if(l==1)    legendbis[l]->AddEntry(    fHistYieldvsErrSoloStatRatio[j], Sys[2*j+1], "pl");
    //    legend[l]->AddEntry(    fHistYieldvsErrSoloStatRatio[j], Sys[j], "pel");
    if(j==0)legendbis[l]->Draw();
    if(j==0)legendtype[l]->Draw();
    // fHistYieldvsErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStat");
    // fHistYieldvsErrSoloSist[j]->SetMarkerColor(Color[j]);
    // fHistYieldvsErrSoloSist[j]->Draw("samee");
  }
    canvasratiobis[l]->SaveAs("FinalOutput/DATA2016/Yield"+DATAorMC[l]+".pdf");
  }

  /*
 TCanvas *canvasbis[2];
  for(Int_t l=0; l<2; l++){
  Int_t count=0;
  canvasbis[l] = new TCanvas(Form("canvasbis%i", l), "NotRatio", 1300, 1000);
  cout << l << endl;
  for(Int_t j=0; j<3; j++){
    if(l==0){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJet.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulk.root";
    if (j==2)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAll.root";
    col=0;
    }
    if(l==1){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJetMC.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulkMC.root";
    if (j==2)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAllMC.root";
    col=1;
    }
    filein=new TFile(PathIn, "");
    fHistYieldvsErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSist");
    fHistYieldvsErrSoloStat[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStat");
    fHistYieldvsErrSoloStat[j]->SetMarkerColor(Color[2*j]);
    fHistYieldvsErrSoloStat[j]->SetMarkerStyle(Marker[col]);
    fHistYieldvsErrSoloStat[j]->SetLineColor(Color[2*j]);
    fHistYieldvsErrSoloStat[j]->Draw("samee");
    fHistYieldvsErrSoloSist[j]->SetMarkerColor(Color[2*j]);
    fHistYieldvsErrSoloSist[j]->SetMarkerStyle(Marker[col]);
    fHistYieldvsErrSoloSist[j]->SetLineColor(Color[2*j]);
    fHistYieldvsErrSoloSist[j]->SetFillStyle(0);
    fHistYieldvsErrSoloSist[j]->Draw("samee2");
        
    if(j==0)legendbis[l]->Draw();
  }
  }
  */

  cout << "\n\n\ndisegno errori relativi in tre canvas (1. JET, 2. OJ 3.J+OJ)" << endl;
  TCanvas *canvaserr[3];


  for(Int_t l=0; l<3; l++){
  Int_t count=0;
    canvaserr[l] = new TCanvas(Form("canvaserr%i", l), Region[l], 1300, 1000);
    cout << l << endl;
   
  for(Int_t j=0; j<2; j++){
    //    cout << l << "  " << j << endl;
    count++;
    if (l==0){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJet.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJetMC.root";
    }
    if (l==1){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulk.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulkMC.root";
    }
    if (l==2){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAll.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAllMC.root";
    }
    filein=new TFile(PathIn, "");
    
    fHistYieldvsErrErrSoloStat[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloStat");
    fHistYieldvsErrErrSoloStat[j]->SetTitle("");
    fHistYieldvsErrErrSoloStat[j]->SetMarkerColor(2);
    fHistYieldvsErrErrSoloStat[j]->SetMarkerStyle(Marker[j]);
    fHistYieldvsErrErrSoloStat[j]->SetLineColor(2);
    fHistYieldvsErrErrSoloStat[j]->Draw("sameep");

    fHistYieldvsErrErrSoloStat[j]->GetXaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloStat[j]->GetXaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloStat[j]->GetYaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloStat[j]->GetYaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloStat[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrErrSoloStat[j]->SetMarkerSize(msize);
      

    // fHistYieldvsErrTotal[j]=(TH1D*)filein->Get("fHistYieldvsErrTot");
       // fHistYieldvsErrTotal[j]->SetMarkerColor(Color[j]);
       // fHistYieldvsErrTotal[j]->SetMarkerStyle(Marker[j]);
       // fHistYieldvsErrTotal[j]->SetLineColor(Color[j]);
       // fHistYieldvsErrTotal[j]->Draw("same");

    fHistYieldvsErrErrSoloStatMB[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloStatMB");
    fHistYieldvsErrErrSoloStatMB[j]->SetTitle("");
    fHistYieldvsErrErrSoloStatMB[j]->SetMarkerColor(1);
    fHistYieldvsErrErrSoloStatMB[j]->SetMarkerStyle(MarkerMBStat[j]);
    fHistYieldvsErrErrSoloStatMB[j]->SetLineColor(1);
    fHistYieldvsErrErrSoloStatMB[j]->SetMarkerSize(msize);
    fHistYieldvsErrErrSoloStatMB[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrErrSoloStatMB[j]->Draw("sameep");

        
    fHistYieldvsErrErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloSist");
    fHistYieldvsErrErrSoloSist[j]->SetTitle("");
    fHistYieldvsErrErrSoloSist[j]->SetMarkerColor(4);
    fHistYieldvsErrErrSoloSist[j]->SetMarkerStyle(Marker[j]);
    fHistYieldvsErrErrSoloSist[j]->SetLineColor(4);
    fHistYieldvsErrErrSoloSist[j]->SetFillStyle(0);
    fHistYieldvsErrErrSoloSist[j]->SetMarkerSize(msize);
    fHistYieldvsErrErrSoloSist[j]->Draw("sameep");

    fHistYieldvsErrErrSoloSist[j]->GetXaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloSist[j]->GetXaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloSist[j]->GetYaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloSist[j]->GetYaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloSist[j]->GetYaxis()->SetRangeUser(0,0.14);
    fHistYieldvsErrErrSoloSist[j]->GetXaxis()->SetRangeUser(0,25);

    fHistYieldvsErrErrSoloSistMB[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloSistMB");
    fHistYieldvsErrErrSoloSistMB[j]->SetTitle("");
    fHistYieldvsErrErrSoloSistMB[j]->SetMarkerColor(1);
    fHistYieldvsErrErrSoloSistMB[j]->SetMarkerStyle(MarkerMBSist[j]);
    fHistYieldvsErrErrSoloSistMB[j]->SetLineColor(1);
    fHistYieldvsErrErrSoloSistMB[j]->SetMarkerSize(msize);
    fHistYieldvsErrErrSoloSistMB[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrErrSoloSistMB[j]->Draw("sameep");

        
    if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloStat[j], StatErr[j], "pl");
    if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloSist[j], SysErr[j], "pl");
    if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloStatMB[j], StatErrMB[j], "pl");
    if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloSistMB[j], SysErrMB[j], "pl");
    if(count==2)legenderr->Draw();
    
  }
  canvaserr[l] ->SaveAs("FinalOutput/DATA2016/RelErr"+Region[l]+".pdf"); 
  }


  cout << "\n\n\ndisegno errori relativi in due canvas (1. DATA, 2. MC)" << endl;
 TCanvas *canvaserrall[2];

 for(Int_t l=0; l<2; l++){
  Int_t count=0;
    canvaserrall[l] = new TCanvas(Form("canvaserrall%i", l), Form("canvaserrall%i", l), 1300, 1000);
    cout << l << endl;
   
  for(Int_t j=0; j<3; j++){
    //    cout << l << "  " << j << endl;
    count++;
    if (l==0){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJet.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulk.root";
    if (j==2)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAll.root";
    }
    if (l==1){
    if (j==0)   PathIn="FinalOutput/DATA2016/SystematicAnalysisJetMC.root";
    if (j==1)   PathIn="FinalOutput/DATA2016/SystematicAnalysisBulkMC.root";
    if (j==2)   PathIn="FinalOutput/DATA2016/SystematicAnalysisAllMC.root";
    }
    filein=new TFile(PathIn, "");
    
    fHistYieldvsErrErrSoloStat[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloStat");
    fHistYieldvsErrErrSoloStat[j]->SetTitle("");
    fHistYieldvsErrErrSoloStat[j]->SetMarkerColor(ColorBis[j]);
    fHistYieldvsErrErrSoloStat[j]->SetMarkerStyle(MarkerTris[0]);
    fHistYieldvsErrErrSoloStat[j]->SetLineColor(ColorBis[j]);
    fHistYieldvsErrErrSoloStat[j]->Draw("sameep");

    fHistYieldvsErrErrSoloStat[j]->GetXaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloStat[j]->GetXaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloStat[j]->GetYaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloStat[j]->GetYaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloStat[j]->GetXaxis()->SetRangeUser(0,25);
    fHistYieldvsErrErrSoloStat[j]->SetMarkerSize(msize);
      

    // fHistYieldvsErrTotal[j]=(TH1D*)filein->Get("fHistYieldvsErrTot");
       // fHistYieldvsErrTotal[j]->SetMarkerColor(Color[j]);
       // fHistYieldvsErrTotal[j]->SetMarkerStyle(Marker[j]);
       // fHistYieldvsErrTotal[j]->SetLineColor(Color[j]);
       // fHistYieldvsErrTotal[j]->Draw("same");

    fHistYieldvsErrErrSoloStatMB[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloStatMB");
    fHistYieldvsErrErrSoloStatMB[j]->SetTitle("");
    fHistYieldvsErrErrSoloStatMB[j]->SetMarkerColor(1);
    fHistYieldvsErrErrSoloStatMB[j]->SetMarkerStyle(MarkerTris[0]);
    fHistYieldvsErrErrSoloStatMB[j]->SetLineColor(1);
    fHistYieldvsErrErrSoloStatMB[j]->SetMarkerSize(msize);
    fHistYieldvsErrErrSoloStatMB[j]->GetXaxis()->SetRangeUser(0,25);
    //fHistYieldvsErrErrSoloStatMB[j]->Draw("sameep");

        
    fHistYieldvsErrErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloSist");
    fHistYieldvsErrErrSoloSist[j]->SetTitle("");
    fHistYieldvsErrErrSoloSist[j]->SetMarkerColor(ColorBis[j]);
    fHistYieldvsErrErrSoloSist[j]->SetMarkerStyle(MarkerTris[1]);
    fHistYieldvsErrErrSoloSist[j]->SetLineColor(ColorBis[j]);
    fHistYieldvsErrErrSoloSist[j]->SetFillStyle(0);
    fHistYieldvsErrErrSoloSist[j]->SetMarkerSize(msize);
    fHistYieldvsErrErrSoloSist[j]->Draw("sameep");

    fHistYieldvsErrErrSoloSist[j]->GetXaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloSist[j]->GetXaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloSist[j]->GetYaxis()->SetTitleOffset(1.1);
    fHistYieldvsErrErrSoloSist[j]->GetYaxis()->SetTitleSize(0.042);
    fHistYieldvsErrErrSoloSist[j]->GetYaxis()->SetRangeUser(0,0.14);
    fHistYieldvsErrErrSoloSist[j]->GetXaxis()->SetRangeUser(0,25);

    fHistYieldvsErrErrSoloSistMB[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloSistMB");
    fHistYieldvsErrErrSoloSistMB[j]->SetTitle("");
    fHistYieldvsErrErrSoloSistMB[j]->SetMarkerColor(1);
    fHistYieldvsErrErrSoloSistMB[j]->SetMarkerStyle(MarkerTris[1]);
    fHistYieldvsErrErrSoloSistMB[j]->SetLineColor(1);
    fHistYieldvsErrErrSoloSistMB[j]->SetMarkerSize(msize);
    fHistYieldvsErrErrSoloSistMB[j]->GetXaxis()->SetRangeUser(0,25);
    //fHistYieldvsErrErrSoloSistMB[j]->Draw("sameep");

        
    if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloStat[j], StatErrBis[j], "pl");
    if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloSist[j], SysErrBis[j], "pl");
    // if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloStatMB[j], StatErrMB[j], "pl");
    // if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloSistMB[j], SysErrMB[j], "pl");
    if(count==3)legenderrbis->Draw();
    
  }


 }
  canvaserrall[0] ->SaveAs("FinalOutput/DATA2016/RelErrAllDATA.pdf"); 
  canvaserrall[1] ->SaveAs("FinalOutput/DATA2016/RelErrAllMC.pdf");   
}

