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


void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title){
  histo->GetYaxis()->SetRangeUser(Low, Up);
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

void PtSpectraBis(Int_t type=0,  Int_t TypeAnalysis=1,Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=0,   Int_t israp=0,TString year="1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"*//*"Run2DataRed_hXi"/"2016kehjl_hK0s"*/,  TString Path1 ="_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/,  Bool_t isEfficiency=1,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=0, Int_t PtBinning=1, TString FitFixed="Boltzmann"/*"Fermi-Dirac"*/,   Int_t sys=0){

  if (TypeAnalysis>2) {cout << "sys errors not yet implemented for these regions " << endl; return;}

  if (type==0){
    PtBinning=1;
    if (!isMC) year= "1617_hK0s";
    else  year= "1617MC_hK0s";
  }
  else if (type==8){
    PtBinning=0;
    if (!isMC) year="Run2DataRed_MECorr_hXi";
    else  year=  "AllMC_hXi";
  }


  TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};

  if (ishhCorr && type!=0){
    type=0;
  }// {cout << " for hh correlation type should be set to zero in order to get the correct files " << endl; return;}

  const   Int_t NumberTypeAnalysis=11;
  cout << "here's the meaning of different values of TypeAnalysis:" << endl;
  cout << 0 << "in-jet production " << endl;
  cout << 1 << "out-of-jet production  (Delta Phi between 1 and 2)" << endl;
  cout << 2 << "inclusive production (from JetBulkEffCorr)" << endl;
  cout << 3 << "inclusive production not scaled by DeltaEta and DeltaPhi (for comparison with published data) (from JetBulkEffCorr)" << endl;

  if (TypeAnalysis>10) return;

  if (year != "2018f1_extra" && year != "2016k" && year != "2018f1_extra_onlyTriggerWithHighestPt" && year != "2016k_onlyTriggerWithHighestPt") {
    //    cout << "output file should be changed: it must include the year name to avoid overwriting output files " << endl;
    //    return;
  } 

  gStyle->SetOptStat(0);

  TString hhCorr[2]={"", "_hhCorr"};
  TFile *fileinbis;
  TFile *filein;

  TString PathIn1;
  TString file = year;
  if (PtBinning!=0)  file+=Form("_PtBinning%i", PtBinning);
  file+=Path1;

  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9
  const Int_t numPtTrigger=1;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;

  const Int_t numsysDPhi =3;

  Int_t PtV0Min = 0; //0 el
  if (type>0 || ishhCorr )   PtV0Min =1;
  if (type==0 && PtBinning!=0) PtV0Min=1;

  cout << "ok " << endl;
  TString JetOrBulk[NumberTypeAnalysis]={"Jet", "Bulk", "All", "AllNS", "AwaySide", "AllButJet", "JetBFit", "JetZYAM", "ASBFit", "ASZYAM", "JetNS"};
  TString DataOrMC[2]={"Data", "MC"};
  Int_t jet=TypeAnalysis;

  cout << "ok " << endl;
  const Int_t numtipo=10;
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  cout << "ok " << endl;
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  Float_t SpectrumSup[numtipo][NumberTypeAnalysis]={{0.015,0.2, 0.2,0.8,0.2,0.2,0.015, 0.015, 0.03, 0.03, 0.03},{0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001, 0.002, 0.002, 0.002},{0.015,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001}};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};

  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};

  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};

  for (Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Nmolt[m] = Nmolt0[m];
      Smolt[m] = Smolt0[m];
      SmoltLegend[m] = SmoltLegend0[m];
    }
    else     if (MultBinning==1){
      Nmolt[m] = Nmolt1[m];
      Smolt[m] = Smolt1[m];
      SmoltLegend[m] = SmoltLegend1[m];
    }
    else     if (MultBinning==2){
      Nmolt[m] = Nmolt2[m];
      Smolt[m] = Smolt2[m];
      SmoltLegend[m] = SmoltLegend2[m];
    }

  }



  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", "", ""};  
  //TString SPtV0[numPtV0]={"", "0.5-1", "0.5-1",  "1-1.5","1.5-2", "2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0-0.5"};
    SPtV0[1]={"0.5-1"};
  }

  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100, 100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", "", ""};
  SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }


  TString SPtV01[numPtV0]={"0-0.1","0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0, 0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0-0.1", "0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0,0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};
  Int_t numPtV0Max=numPtV0;
  if (PtBinning==1 || PtBinning==2) numPtV0Max = numPtV0;
  else numPtV0Max = numPtV0-2;

  if (PtBinning==1){
    for(Int_t v=0; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=0; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }
  cout << "ok " << endl;
  TString SErrorSpectrum[3]={"stat.","syst. uncorr.","syst. corr."};
  TString SSystJet[2]={"BC [-1.0, 1.0]", "BC [-1.2, 1.2]"};
  TString SSystBulk[3]={"BC [1.0, 2.0]", "BC [2.0, 4.28]", "BC [1.0, 4.28]"};

  //  Int_t Marker[numSyst]={20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  //  Int_t MarkerBetter[numSyst]={1,22, 32, 1, 29, 1,   3,  34, 33, 1, 20, 21, 22};
  //  Int_t Color[numSyst]= {1,  2,  3,  4,  5,  6,  7,  7,  4, 10,  6,  1,  2};
  //  Int_t Color[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 628,868,1, 909, 1,801,418,860,  1, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};
  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Color[3] ={628, 418, 600};

  TString stringout;
  TString PathIn0;
  TString PathIn;
  PathIn0 = Dir+"/DATA"+year0+"/PtSpectra" +hhCorr[ishhCorr];
  if (PtBinning>0) PathIn0 +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      PathIn0 +="_"+tipo[type];
    PathIn0 +=Srap[israp];
    PathIn0 +=SSkipAssoc[SkipAssoc];
  }

  stringout = Dir+"/DATA"+year0+"/PtSpectraBis" +hhCorr[ishhCorr];
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+=   Form("_PtMin%.1f_", PtTrigMin);
  stringout+= RegionType[TypeAnalysis];
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  cout << " definition integral regions " << endl;
  //********************************************************************* 
  //definition of integral regions in deltaPhi distribution
  //********************************************************************* 
  Float_t ALowBin[3]={-1}; 
  Float_t AUpBin[3]={1};
  Float_t ALowBinFit[3]={-1};
  Float_t AUpBinFit[3]={1};
  Float_t DeltaPhiWidth[3]={0};
  Float_t DeltaPhiWidthApprox[3]={0};

  if (TypeAnalysis==0 || TypeAnalysis==10){
    if (!ishhCorr && type==0){// || type==8)){
      ALowBinFit[0]=	ALowBin[0]={-0.8};
      AUpBinFit[0]=	AUpBin[0]={0.8};
      
      ALowBinFit[1]=	ALowBin[1]={-1.2};
      AUpBinFit[1]=	AUpBin[1]={1.2};
    }
    else if (!ishhCorr && (type==4 || type==5 || type==8)){
      ALowBinFit[0]=	ALowBin[0]={-0.7};
      AUpBinFit[0]=	AUpBin[0]={0.7};
      
      ALowBinFit[1]=	ALowBin[1]={-1.1};
      AUpBinFit[1]=	AUpBin[1]={1.1};
    }
    if (ishhCorr){
      ALowBinFit[0]=	ALowBin[0]={-1.2};
      AUpBinFit[0]=	AUpBin[0]={1.2};
      
      ALowBinFit[1]=	ALowBin[1]={-1.4};
      AUpBinFit[1]=	AUpBin[1]={1.4};
    }
    
  }
  else if (TypeAnalysis==1){
    ALowBinFit[0]=	ALowBin[0]={1};
    AUpBinFit[0]=	AUpBin[0]={2};

    ALowBinFit[1]=	ALowBin[1]={2};
    AUpBinFit[1]=	AUpBin[1]={4.28};

    ALowBinFit[2]=	ALowBin[2]={1};
    AUpBinFit[2]=	AUpBin[2]={4.28};

  }
  else if (TypeAnalysis==2 || TypeAnalysis==3){
    ALowBinFit[0]=	ALowBin[0]=-0.5*TMath::Pi()+0.01;
    AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
  }
  else if (TypeAnalysis==4){
    ALowBinFit[0]=	ALowBin[0]=TMath::Pi()-1;
    AUpBinFit[0]=	AUpBin[0]=TMath::Pi()+1;
  }
  else if (TypeAnalysis==5){
    ALowBinFit[0]=	ALowBin[0]=1;
    AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
  }
  else if (TypeAnalysis==6 || TypeAnalysis==7){
    if (type==0){
      ALowBinFit[0]=	ALowBin[0]=-0.8;
      AUpBinFit[0]=	AUpBin[0]=0.8;
    }
    else if (type==8){
      ALowBinFit[0]=	ALowBin[0]=-0.7;
      AUpBinFit[0]=	AUpBin[0]=0.7;
    }
  }
  else if (TypeAnalysis==8 || TypeAnalysis==9){
    if (type==0){
      ALowBinFit[0]=	ALowBin[0]=1.6;
      AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
    }
    else if (type==8){
      ALowBinFit[0]=	ALowBin[0]=1.6;
      AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
    }
  }

  //low and upper ranges for the fit****************************
  //************************************************************

  Double_t LowRange[nummolt+1]= {2, 2, 1, 1.5, 1.5, 1.5}; 
  Double_t UpRange[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeSpectrumPart[nummolt+1]= {0};
  Double_t UpRangeSpectrumPart[nummolt+1]= {0};

  Double_t LowRangeJet[nummolt+1]= {2, 1.5, 1, 1.5, 1.5, 1}; 
  Double_t UpRangeJet[nummolt+1]= {8,8,8,8,8,8}; 

  //  Double_t LowRangeBulk[nummolt+1]= {1, 1, 1, 1, 1, 1}; 
  Double_t LowRangeBulk[nummolt+1]= {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; 
  Double_t UpRangeBulk[nummolt+1]= {4,4,4,4,4,4};

  //  Double_t LowRangeAll[nummolt+1]=  {1, 1, 1, 1, 1, 1}; 
  Double_t LowRangeAll[nummolt+1]=  {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; 
  Double_t UpRangeAll[nummolt+1]= {4,4,4,4,4,4};

  Double_t LowRangeAS[nummolt+1]= {1, 1, 1, 1, 1, 1}; 
  Double_t UpRangeAS[nummolt+1]=  {4,4,4,4,4,4};

  Double_t LowRangeAllButJet[nummolt+1]=  {1, 1, 1, 1, 1, 1}; 
  Double_t UpRangeAllButJet[nummolt+1]={4,4,4,4,4,4};

  Double_t LowRangeJetBFit[nummolt+1]= {1.5, 1.5, 1,1, 1.5,1}; 
  Double_t UpRangeJetBFit[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeJetZYAM[nummolt+1]= {1.5 , 1.5, 0.5,1, 1.5,1}; 
  Double_t UpRangeJetZYAM[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeASBFit[nummolt+1]= {1.5 , 1.5, 0.5,1, 1.5,1}; 
  Double_t UpRangeASBFit[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeASZYAM[nummolt+1]= {1.5 , 1.5, 0.5,1, 1.5,1}; 
  Double_t UpRangeASZYAM[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeJet1[nummolt+1]= {2, 2, 1.5, 1, 1, 1}; 
  Double_t UpRangeJet1[nummolt+1]= {8,8,8,8,8,8}; 
  Double_t LowRangeJet1MC[nummolt+1]= {2, 1.5, 1.5, 1.5, 1.5, 1.5}; 
  Double_t UpRangeJet1MC[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeJet2[nummolt+1]= {2, 1.5, 1, 1, 1, 1}; 
  Double_t UpRangeJet2[nummolt+1]= {8,8,8,8,8,8}; 
  Double_t LowRangeJet2MC[nummolt+1]= {2, 1.5, 1.5, 1.5, 1.5, 1.5}; 
  Double_t UpRangeJet2MC[nummolt+1]= {8,8,8,8,8,8}; 

  if (!ishhCorr && type==0){
    for (Int_t m =0; m<nummolt+1; m++){
      if (m==4 && !isMC)     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1.0;
      else   LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.5;
      LowRangeBulk[m]= 0; //0
      LowRangeAll[m]= 0; //0 
      UpRangeAll[m]= 4; 
    }
    if (PtBinning==1) {
      for (Int_t m =0; m<nummolt+1; m++){
	if (!isMC){
	  if (m==4 || m==3)     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.5;
	  else if (m==5 )     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0;
	  else  if (m==2)  {LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0;}
	  else   {LowRangeJet[m] = 0.8;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.5;}
	  if (m==4) UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=3;
	  else	if (m==3) UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
	  else UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
	  LowRangeBulk[m]= 0.1; 
	  LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 0.1; 
	  UpRangeAll[m]= 2; 
	  UpRangeBulk[m]= 2; 
	  UpRangeASBFit[m]= 	  UpRangeASZYAM[m]= 2; 
	}
      }
    }
    if (PtBinning==1 && year=="1617_hK0s") {

      for (Int_t m =0; m<nummolt+1; m++){
	if (PtTrigMin==3){
	  if (!isMC){
	    if (m==5|| m==2 || m==4 )     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	    else  {LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;}
	    UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;	  
	    if (m==4) UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
	    LowRangeBulk[m]= 0.1; 
	    LowRangeAll[m]= 0.1; 
	    UpRangeAll[m]= 2; 
	    UpRangeBulk[m]= 2; 
	    LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 0.1; 
	    UpRangeASBFit[m]= 	  UpRangeASZYAM[m]= 2; 
	  }
	}
     
	else if (PtTrigMin==4){
	  if (!isMC){
	    if (m==5)     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	    else  {LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;}
	    UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;	  
	    LowRangeBulk[m]= 0.1; 
	    LowRangeAll[m]= 0.1; 
	    UpRangeAll[m]= 2; 
	    UpRangeBulk[m]= 2; 
	    LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 0.1; 
	    UpRangeASBFit[m]= 	  UpRangeASZYAM[m]= 2; 
	  }

	}
      }
    }


  }

  if (MultBinning==1 && type==8){
    for (Int_t m=0; m<nummolt+1; m++){
      if (!isMC)      LowRangeJet[m] = LowRangeJet1[m];
      else if (isMC)      LowRangeJet[m] = LowRangeJet1MC[m];
    }
  }

  if (MultBinning==2 && type==8){
    for (Int_t m=0; m<nummolt+1; m++){
      if (!isMC)      LowRangeJet[m] = LowRangeJet2[m];
      else if (isMC)      LowRangeJet[m] = LowRangeJet2MC[m];
    }
  }

  if (MultBinning==0 && type==8){
    for (Int_t m=0; m<nummolt+1; m++){
      if (isMC)  {
	if (m==3)     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1.5;
	else  if (m==0 || m==4 || m==5)  {LowRangeJet[m] = 1.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1;}
	else  if (m==1) {LowRangeJet[m] = 2.0;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1;}
	else  if (m==2) {LowRangeJet[m] = LowRangeJetBFit[m]=1.5;       LowRangeJetZYAM[m]= 1;}
	UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;

	LowRangeBulk[m]= 1; //not enought statistics below 
	LowRangeAll[m]= 1; 
	UpRangeAll[m]= 4; 
	UpRangeBulk[m]= 4; 
	if (m==5) {
	  LowRangeBulk[m]= 0.5;
	  LowRangeAll[m]= 0.5; 
	  UpRangeAll[m]= 8; 
	  UpRangeBulk[m]= 8; 
	}
      }
      else {
	LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 1.5; 
	if (m==5 || m==3) 	LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 1; 
	UpRangeASBFit[m]= 	  UpRangeASZYAM[m]= 4; 
      }
    }
  }


  //end of low and upper ranges for the fit****************************
  //************************************************************

  Int_t PtBinMin[nummolt+1]={0};
  for (Int_t m=0; m<nummolt+1; m++){
    if (TypeAnalysis==0 || TypeAnalysis==10){
      LowRange[m] =  LowRangeJet[m];
      UpRange[m] =  UpRangeJet[m];
    }
    if (TypeAnalysis==1){
      LowRange[m] =  LowRangeBulk[m];
      UpRange[m] =  UpRangeBulk[m];
    }
    if (TypeAnalysis==2 || TypeAnalysis==3){
      LowRange[m] =  LowRangeAll[m];
      UpRange[m] =  UpRangeAll[m];
    }
    if (TypeAnalysis==4){
      LowRange[m] =  LowRangeAS[m];
      UpRange[m] =  UpRangeAS[m];
    }
    if (TypeAnalysis==5){
      LowRange[m] =  LowRangeAllButJet[m];
      UpRange[m] =  UpRangeAllButJet[m];
    }
    if (TypeAnalysis==6){
      LowRange[m] =  LowRangeJetBFit[m];
      UpRange[m] =  UpRangeJetBFit[m];
    }
    if (TypeAnalysis==7){
      LowRange[m] =  LowRangeJetZYAM[m];
      UpRange[m] =  UpRangeJetZYAM[m];
    }
    if (TypeAnalysis==8){
      LowRange[m] =  LowRangeASBFit[m];
      UpRange[m] =  UpRangeASBFit[m];
    }
    if (TypeAnalysis==9){
      LowRange[m] =  LowRangeASZYAM[m];
      UpRange[m] =  UpRangeASZYAM[m];
    }

    //bin counting is done at pt> LowRangeSpectrumPart, integral of fit function at pt <  LowRangeSpectrumPart
    if  (TypeAnalysis==1 || TypeAnalysis==2){
      if (type>0 && !isMC)    LowRangeSpectrumPart[m] = 0.5;
      else if (type>0 && isMC && m!=5)    LowRangeSpectrumPart[m] = 1.;
      else     LowRangeSpectrumPart[m] = 0.1;
    }
    else  LowRangeSpectrumPart[m] = LowRange[m];
    //    else  LowRangeSpectrumPart[m] = 8;
    cout << " ok " << endl;
    for(Int_t v=0; v < numPtV0Max; v++){
      if (LowRangeSpectrumPart[m]<=NPtV0[v]) break;
      PtBinMin[m]++;
      cout << "LowRange " << LowRangeSpectrumPart[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
    }
  }

  //syst1_limit indica gli effetti sistematici sullo spettro in pT
  Int_t syst1_limit=0;
  if(!TypeAnalysis==1) syst1_limit= 1; //2 for BC systematic
  if (TypeAnalysis==1) syst1_limit= 1;
  if (TypeAnalysis==2) syst1_limit= 1; //non considero sistematico associato a bin counting per total
  if (TypeAnalysis>2) syst1_limit=1;


  Int_t sysV0=0;    
  if(isMC && isEfficiency) file = year + "_MCEff" +Path1;
  //  file+=Form("_PtBinning%i", PtBinning);
  //  file+= Path1;

  //here I start the newly  created part!!

  TString titleX=  "p_{T} (GeV/c)";
  TString titleY=  "1/#Delta #eta #Delta #phi 1/N_{trigg} dN/dp_{T}";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 

  TH1F* fHistSpectrumStat[nummolt+1];
  TH1F* fHistSpectrumSist[nummolt+1];
  TH1F* fHistSpectrum_max[nummolt+1];
  TH1F* fHistSpectrum_min[nummolt+1];

  TH1F* fHistSpectrumStatDPhi[nummolt+1][numsysDPhi];
  TH1F* fHistSpectrumSistDPhi[nummolt+1];
  TH1F* fHistSpectrumSistOOJ[nummolt+1];
  TH1F* fHistSpectrumSistAll[nummolt+1];

  //histos for relative uncertainty
  TH1F* fHistSpectrumStatRelError[nummolt+1]; //stat
  TH1F* fHistSpectrumSistRelError[nummolt+1]; //sist on the DeltaPhi projections
  TH1F* fHistSpectrumSistRelErrorDPhi[nummolt+1]; //sist assoc to choice of DeltaPhi interval
  TH1F* fHistSpectrumSistRelErrorAll[nummolt+1]; //total sist

  TCanvas* canvasPtSpectra = new TCanvas ("canvasPtSpectra", "canvasPtSpectra", 1300, 800);
  canvasPtSpectra->Divide(3,2);

  TCanvas* canvasPtSpectraFit = new TCanvas ("canvasPtSpectraFit", "canvasPtSpectraFit", 1300, 800);
  canvasPtSpectraFit->Divide(3,2);

  TCanvas* canvasPtSpectraAll = new TCanvas ("canvasPtSpectraAll", "canvasPtSpectraAll", 1300, 800);
  canvasPtSpectraAll->Divide(3,2);
  TCanvas* canvasBarlow = new TCanvas ("canvasBarlow", "canvasBarlow", 1300, 800);
  canvasBarlow->Divide(3,2);
  TCanvas* canvasPtSpectraRelError = new TCanvas ("canvasPtSpectraRelError", "canvasPtSpectraRelError", 1300, 800);
  canvasPtSpectraRelError->Divide(3,2);
  TCanvas* canvasPtSpectraRelErrorAll = new TCanvas ("canvasPtSpectraRelErrorAll", "canvasPtSpectraRelErrorAll", 1300, 800);
  canvasPtSpectraRelErrorAll->Divide(3,2);

  Bool_t BarlowPassed[nummolt+1][numsysDPhi]={0};
  Int_t BarlowSign[nummolt+1][numsysDPhi]={0};
  TH1F* hBarlowVar[nummolt+1][numsysDPhi];
  Float_t BarlowVar[nummolt+1][numPtV0][numsysDPhi]={0};
  Float_t ErrDPhi[nummolt+1][numPtV0]={0};

  Int_t ColorsysDPhi[12] = { 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};

  Float_t LimSupYield=0.1;
  if (type==8){
    LimSupYield=0.01;
  }

  Float_t LimSup=0.2;
  if (type==0){
    if (TypeAnalysis==0) LimSup =0.02;
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSup =0.001;
    else   if (TypeAnalysis==1) LimSup =0.01;
    else   if (TypeAnalysis==2) LimSup =0.01;
  }
  Float_t LimSupError=0.01;
  if (type==0) {
    if (TypeAnalysis==0) LimSupError =0.15;
    else   if (TypeAnalysis==1) LimSupError =0.06;
    else   if (TypeAnalysis==2) LimSupError =0.02;
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSupError =0.4;
    else   if (TypeAnalysis==1) LimSupError =0.01;
    else   if (TypeAnalysis==2) LimSupError =0.01;
  }

  //first part: I get default spectra with stat uncertainty
  cout << "\n**************************"<<endl;
  cout << "first part: I get default spectra with stat uncertainty\n" << endl;
  TString  PathInDef=PathIn0;
  PathInDef+=   Form("_SysPhi%i_PtMin%.1f_", 0, PtTrigMin);
  PathInDef+= RegionType[TypeAnalysis];
  PathInDef += ".root";
  cout << "\n\n" << PathInDef << endl;
  TFile *  fileinDef = new TFile(PathInDef, "");
  if (!fileinDef) {cout << PathInDef << " does not exist" << endl; return;}
  for(Int_t m=0; m<nummolt+1; m++){
    //    cout << " m " << m << endl;
    fHistSpectrumStat[m]=(TH1F*)fileinDef->Get("fHistSpectrum_"+Smolt[m]);
    fHistSpectrumSist[m]=(TH1F*)fileinDef->Get("fHistSpectrumSist_"+Smolt[m]);
    if (!fHistSpectrumStat[m]) {cout << "fHistSpectrumStat does not exist" << endl; return;}
    if (!fHistSpectrumSist[m]) {cout << "fHistSpectrumSist does not exist" << endl; return;}

    fHistSpectrumSistRelErrorAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorAll_"+Smolt[m]);
    fHistSpectrumStatRelError[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumStatRelError_"+Smolt[m]);
    fHistSpectrumSistRelErrorDPhi[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorDPhi_"+Smolt[m]);

    fHistSpectrum_max[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMax_"+Smolt[m]);
    fHistSpectrum_min[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMin_"+Smolt[m]);
    fHistSpectrumSistDPhi[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDPhi_"+Smolt[m]);
    fHistSpectrumSistOOJ[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistOOJ_"+Smolt[m]);
    fHistSpectrumSistAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistAll_"+Smolt[m]);

    for(Int_t v=PtV0Min; v <    fHistSpectrumSistAll[m]->GetNbinsX() ; v++){
      //      cout << "v " << v << endl;
      if (fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSist[m]->GetBinError(v+1)/ fHistSpectrumSist[m]->GetBinContent(v+1)); //this works for inclusive, for jet and OOJ I will change it
	fHistSpectrumStatRelError[m]->SetBinContent(v+1, fHistSpectrumStat[m]->GetBinError(v+1)/ fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
      //      cout << "sist rel error " <<    fHistSpectrumSistRelErrorAll[m]->GetBinContent(v+1) << endl;
      fHistSpectrumStatRelError[m]->SetBinError(v+1,0);
      //      cout << "stat rel error " <<    fHistSpectrumStatRelError[m]->GetBinContent(v+1) << endl;

    }
  } //end loop m

  //second part: I evaluate syst uncertainty associated to choice of DeltaPhi (for jet and out of jet)
  Int_t NSign=2;
  if (TypeAnalysis!=2){
    cout << "\n**************************"<<endl;
    cout << " second part: I evaluate syst uncertainty associated to choice of DeltaPhi (for jet and out of jet)\n " << endl;

    cout << "2a. Evaluating Barlow significance of DeltaPhi change" << endl;
    for(Int_t m=0; m<nummolt+1; m++){
      cout << " m " << m << endl;
      for (Int_t sysDPhi=1; sysDPhi<numsysDPhi; sysDPhi++){
	PathIn=PathIn0;
	PathIn+=   Form("_SysPhi%i_PtMin%.1f_", sysDPhi, PtTrigMin);
	PathIn+= RegionType[TypeAnalysis];
	PathIn+= ".root";
	cout << "" << PathIn << endl;
	filein = new TFile(PathIn, "");
	fHistSpectrumStatDPhi[m][sysDPhi]=(TH1F*)filein->Get("fHistSpectrum_"+Smolt[m]);
	fHistSpectrumStatDPhi[m][sysDPhi]->SetName(Form("fHistSpectrumStat_sysDPhi%i", sysDPhi)+Smolt[m]);

	StyleHisto(fHistSpectrumStatDPhi[m][sysDPhi], 0, LimSup, ColorsysDPhi[sysDPhi], 6, titleX, titleY,  title+SmoltLegend[m]);
	StyleHisto(fHistSpectrumStat[m], 0, LimSup, Color[TypeAnalysis], 6, titleX, titleY,  title+SmoltLegend[m]);

	canvasPtSpectraAll->cd(m+1);
	//      cout << "I'm drawing on the canvas " << endl;
	gPad->SetLeftMargin(0.15);
	if (sysDPhi==1) fHistSpectrumStat[m]->Draw("");
	fHistSpectrumStatDPhi[m][sysDPhi]->Draw("same");

	hBarlowVar[m][sysDPhi]=(TH1F*)    fHistSpectrumStatDPhi[m][sysDPhi]->Clone("fHistBarlowVar_"+Smolt[m]);
	for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	  //	cout << "\nv: " << v << endl;
	  //	cout << "...histo...dphi " << sysDPhi << " " << fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1)<< endl;
	  BarlowVar[m][v][sysDPhi] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) + pow(fHistSpectrumStat[m]->GetBinError(v+1),2)));
	  if (TMath::Abs(BarlowVar[m][v][sysDPhi])>2) BarlowSign[m][sysDPhi]++;
	  hBarlowVar[m][sysDPhi] ->SetBinContent(v+1, BarlowVar[m][v][sysDPhi]) ;
	  hBarlowVar[m][sysDPhi] ->SetBinError(v+1, 0) ;
	}//end loop v
	if (BarlowSign[m][sysDPhi]>= NSign) BarlowPassed[m][sysDPhi]=1;
	StyleHisto(hBarlowVar[m][sysDPhi], -5, 5, ColorsysDPhi[sysDPhi], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m]);
	if ( BarlowPassed[m][sysDPhi])       StyleHisto(hBarlowVar[m][sysDPhi], -5, 5, ColorsysDPhi[sysDPhi], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m]);

	cout << "Barlow passed " << BarlowPassed[m][sysDPhi] << endl;
	canvasBarlow->cd(m+1);
	hBarlowVar[m][sysDPhi] ->Draw("same p");

	if (!BarlowPassed[m][sysDPhi]) continue;
	//      cout << "defining max and min for " << endl;
	for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	  //	cout << " v: " << v << endl;
	  //	cout <<"before: " <<  fHistSpectrum_min[m]->GetBinContent(v+1)<< " << " << fHistSpectrumStat[m]->GetBinContent(v+1)<< " << " <<fHistSpectrum_max[m]->GetBinContent(v+1)<<endl;
	  //	cout << "dphi " << sysDPhi << " " << fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1)<< endl;
	  if(fHistSpectrum_max[m]->GetBinContent(v+1) < fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1)){
	    fHistSpectrum_max[m]->SetBinContent(v+1, fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1));
	  }
	  if(fHistSpectrum_min[m]->GetBinContent(v+1) > fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1)){
	    fHistSpectrum_min[m]->SetBinContent(v+1, fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1));
	  }
	  //cout <<"after: " <<  fHistSpectrum_min[m]->GetBinContent(v+1)<< " << " << fHistSpectrumStat[m]->GetBinContent(v+1)<< " << " <<fHistSpectrum_max[m]->GetBinContent(v+1)<<endl;


	} //end loop v


      } //end loop sysDPhi

      cout << "2b. calculating syst error associated to DeltaPhi choice" << endl;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	//      cout << " v: " << v << endl;
	ErrDPhi[m][v]= TMath::Abs(	  fHistSpectrum_min[m]->GetBinContent(v+1) - 	  fHistSpectrum_max[m]->GetBinContent(v+1))/2;
	//      cout << fHistSpectrum_min[m]->GetBinContent(v+1)<< " << " << fHistSpectrumStat[m]->GetBinContent(v+1)<< " << " <<fHistSpectrum_max[m]->GetBinContent(v+1)<<endl;
	//      cout << "err dphi " << ErrDPhi[m][v] << endl;
	fHistSpectrumSistDPhi[m]->SetBinError(v+1,ErrDPhi[m][v]);

	if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1,ErrDPhi[m][v]/ fHistSpectrumStat[m]->GetBinContent(v+1));
	}
	else       fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1,0);
	fHistSpectrumSistRelErrorDPhi[m]->SetBinError(v+1,0);

	fHistSpectrumSistAll[m]->SetBinError(v+1,sqrt(pow(fHistSpectrumSistDPhi[m]->GetBinError(v+1),2) + pow(fHistSpectrumSist[m]->GetBinError(v+1),2)));

	if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSistAll[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	}
	else       fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1,0);
	fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
	cout << " stat rel error " <<       fHistSpectrumStatRelError[m]->GetBinContent(v+1)<< endl;
	cout << " sist rel error assoc to DeltaPhi choice " <<       fHistSpectrumSistRelErrorDPhi[m]->GetBinContent(v+1)<< endl;
      }

      //    cout << "drawing sist rel error " << endl;
      StyleHisto(fHistSpectrumSistRelErrorDPhi[m], 0, LimSupError,807 , 6, titleX, titleYRel, title+SmoltLegend[m]);
      canvasPtSpectraRelErrorAll->cd(m+1);
      fHistSpectrumSistRelErrorDPhi[m]->Draw("same");
    }//end loop on m
  }

  //third part: for the jet, I evaluate sys associated to type of OOJ subtraction (for the jet)
  if (TypeAnalysis==0 && type==0){
    cout << "\n*******************************************"<<endl;
    cout << "3. sysy associated to pol0 OOJ subtraction for hK0s"<< endl;
  }


  //I plot spectra with stat + total syst uncertainty
  // I plot stat + syst (total, DeltaPhi assoc., OOJ sub assoc.) relative uncertainty

  for(Int_t m=0; m<nummolt+1; m++){

    canvasPtSpectra->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistAll[m], 0, LimSup, Color[TypeAnalysis], 6, titleX, titleY, title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumStat[m], 0, LimSup, Color[TypeAnalysis], 6, titleX, titleY,  title+SmoltLegend[m]);
    fHistSpectrumStat[m]->Draw("same");
    fHistSpectrumSistAll[m]->SetFillStyle(0);
    fHistSpectrumSistAll[m]->Draw("same e2");

    canvasPtSpectraRelErrorAll->cd(m+1);
    StyleHisto(fHistSpectrumSistRelErrorAll[m], 0, LimSupError,Color[TypeAnalysis] , 27, titleX, titleYRel, title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumStatRelError[m], 0, LimSupError,Color[TypeAnalysis] , 33, titleX, titleYRel, title+SmoltLegend[m]);
    fHistSpectrumSistRelErrorAll[m]->Draw("same");
    fHistSpectrumStatRelError[m]->Draw("same");

    canvasPtSpectraRelError->cd(m+1);
    fHistSpectrumSistRelErrorAll[m]->Draw("same");
    fHistSpectrumStatRelError[m]->Draw("same");
  } //end loop m


  //fourth part: fit to obtain pt-integrated yield vs mult
  AliPWGFunc pwgfunc;
  const Int_t numfittipo=4;
  TLegend *legendfit=new TLegend(0.6, 0.6, 0.9, 0.9);
  TString   nameMTscaling[nummolt+1][numfittipo];
  Int_t ColorFit[numfittipo+1]={860, 881, 868, 628, 419};
  TFitResultPtr fFitResultPtr0[nummolt+1][numfittipo];
  TF1* fit_MTscaling[nummolt+1][numfittipo];
  TString       nameFit[numfittipo]={"mT-scaling", "Boltzmann", "Fermi-Dirac", "Levi"};
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  Int_t factor=1;

  Float_t    YieldExtrHighPt[nummolt+1][numfittipo]={0};
  Float_t    YieldExtrLowPt[nummolt+1][numfittipo]={0};

  Float_t    YieldExtrHighPtAvg[nummolt+1]={0};
  Float_t    YieldExtrLowPtAvg[nummolt+1]={0};
  Float_t    YieldErrStatHighPtAvg[nummolt+1]={0};
  Float_t    YieldErrStatLowPtAvg[nummolt+1]={0};
  Float_t    YieldExtr[nummolt+1]={0};
  Float_t    YieldExtrMax[nummolt+1]={0};
  Float_t    YieldExtrMin[nummolt+1]={0};
  Float_t    YieldExtrErrStat[nummolt+1]={0};
  Float_t    YieldExtrErrSist[nummolt+1]={0};
  Float_t    YieldExtrErrSistUp[nummolt+1]={0};
  Float_t    YieldExtrErrSistLow[nummolt+1]={0};
  Float_t    YieldSpectrum[nummolt+1]={0};
  Float_t    YieldSpectrumErrStat[nummolt+1]={0};
  Float_t    YieldSpectrumErrSist[nummolt+1]={0};
  Float_t    Yield[nummolt+1]={0};
  Float_t    YieldErrStat[nummolt+1]={0};
  Float_t    YieldErrSist[nummolt+1]={0};
  Float_t    YieldErrSistUp[nummolt+1]={0};
  Float_t    YieldErrSistLow[nummolt+1]={0};

  for(Int_t m=0; m<nummolt+1; m++){
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (type==8 && TypeAnalysis==10) factor=2; //otherwise fit is not properly done                 
      pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
      nameMTscaling[m][typefit] = Form("fitMTscaling_m%i_fit%i",m, typefit);
      cout << nameFit[typefit]<< endl;
      if (typefit==0)      fit_MTscaling[m][typefit]=    pwgfunc.GetMTExp(massParticle[type], 0.1, 0.04*factor, nameMTscaling[m][typefit]); //mass, T, norm, name                                
      if (typefit==1)      fit_MTscaling[m][typefit]=    pwgfunc.GetBoltzmann(massParticle[type],0.1, 0.04*factor, nameMTscaling[m][typefit]);
      if (typefit==2)      fit_MTscaling[m][typefit]=    pwgfunc.GetFermiDirac(massParticle[type],0.1, 0.04*factor, nameMTscaling[m][typefit]);
      if (typefit==3)  {
	fit_MTscaling[m][typefit]=    pwgfunc.GetLevi(massParticle[type],0.1, 0.03, 0.04*factor, nameMTscaling[m][typefit]);  //norm, n, T, mass (but the function must be called with these parameters in inverse order) 
	fit_MTscaling[m][typefit]->SetParLimits(0, 0,fHistSpectrumStat[m]->GetBinContent(fHistSpectrumStat[m]->GetMaximumBin())*0.5*10);
	fit_MTscaling[m][typefit]->SetParLimits(2, 0.1, 10);
	fit_MTscaling[m][typefit]->SetParLimits(1, 2, 30);
	if (type==8)         fit_MTscaling[m][typefit]->SetParLimits(1, 15,30);
	fit_MTscaling[m][typefit]->SetParameter(2, 0.7);
      }
      if (m==0 )    legendfit->AddEntry( fit_MTscaling[m][typefit], nameFit[typefit],"l");

      fit_MTscaling[m][typefit]->SetLineColor(ColorFit[typefit]);
      fit_MTscaling[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtr0[m][typefit]=       fHistSpectrumStat[m]->Fit(    fit_MTscaling[m][typefit],"SR0");
      fit_MTscaling[m][typefit]->SetRange(0,8);

      canvasPtSpectraFit->cd(m+1);
      fHistSpectrumStat[m]->Draw("same");
      fit_MTscaling[m][typefit]->Draw("same");
      legendfit->Draw("");

      //calculatin yields
      YieldExtrLowPt[m][typefit]= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPt[m][typefit] = fit_MTscaling[m][typefit]->Integral(8,50);

      if (typefit==0){
	YieldExtrMax[m] =  YieldExtrLowPt[m][typefit]+YieldExtrHighPt[m][typefit];
	YieldExtrMin[m] =  YieldExtrLowPt[m][typefit]+YieldExtrHighPt[m][typefit];
      }
      if ((YieldExtrLowPt[m][typefit]+YieldExtrHighPt[m][typefit]) > YieldExtrMax[m]) {
	YieldExtrMax[m] = YieldExtrLowPt[m][typefit]+YieldExtrHighPt[m][typefit];
      }
      if ((YieldExtrLowPt[m][typefit]+YieldExtrHighPt[m][typefit]) < YieldExtrMin[m]) {
	YieldExtrMin[m] = YieldExtrLowPt[m][typefit]+YieldExtrHighPt[m][typefit];
      }

      YieldExtrLowPtAvg[m]+= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvg[m]+= fit_MTscaling[m][typefit]->Integral(8,50);     

      YieldErrStatLowPtAvg[m]+= pow(    fit_MTscaling[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
      YieldErrStatHighPtAvg[m]+= pow(    fit_MTscaling[m][typefit]->IntegralError(8,50,fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);

    } //end loop typefit


    YieldExtrHighPtAvg[m]=YieldExtrHighPtAvg[m]/numfittipo;
    YieldExtrLowPtAvg[m]=YieldExtrLowPtAvg[m]/numfittipo;
    YieldErrStatHighPtAvg[m]=sqrt(YieldErrStatHighPtAvg[m]/numfittipo);
    YieldErrStatLowPtAvg[m]=sqrt(YieldErrStatLowPtAvg[m]/numfittipo);

    YieldExtr[m] =     YieldExtrHighPtAvg[m]+    YieldExtrLowPtAvg[m];
    YieldExtrErrStat[m] = sqrt(    pow(YieldErrStatHighPtAvg[m],2)+    pow(YieldErrStatLowPtAvg[m],2));
    YieldExtrErrSist[m] = (YieldExtrMax[m]-YieldExtrMin[m])/2;
    YieldExtrErrSistUp[m] = YieldExtrMax[m]-YieldExtr[m];
    YieldExtrErrSistLow[m] = YieldExtr[m];

    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      YieldSpectrum[m] +=  ( fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1));
      YieldSpectrumErrStat[m] += pow ( fHistSpectrumStat[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1),2);
      YieldSpectrumErrSist[m] += pow ( fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1),2);
    }
    YieldSpectrumErrStat[m]=sqrt( YieldSpectrumErrStat[m]);
    YieldSpectrumErrSist[m]=sqrt( YieldSpectrumErrSist[m]);

    Yield[m] =  YieldSpectrum[m]+YieldExtr[m];
    YieldErrStat[m] = sqrt(    pow(YieldSpectrumErrStat[m],2) + pow(YieldExtrErrStat[m],2));
    YieldErrSist[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSist[m],2));
    YieldErrSistUp[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistUp[m],2));
    YieldErrSistLow[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistLow[m],2));

    cout << " m " << m << endl;
    cout << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) +- "<< YieldErrSist[m] << " (sist) "<< endl;
    cout << " (stat rel: " <<  YieldErrStat[m]/ Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldErrSist[m]/ Yield[m] << ")= " << endl;
    cout << "from spectrum: " <<   YieldSpectrum[m] << "+- " << YieldSpectrumErrStat[m]<<" (stat) +- " <<YieldSpectrumErrSist[m]<< " (sist)"  << endl;
    cout << "from extr: " <<     YieldExtr[m] << "+- " <<     YieldExtrErrStat[m]<< " (stat)  +- " << YieldExtrErrSist[m]<< " (sist) " << endl;
    if (TypeAnalysis==0 && type==8)     cout << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) + "<< YieldErrSistUp[m]<< " - " << YieldErrSistLow[m] << " (sist) "<< endl;

  }//end loop m

  TH1D* fHistYieldStat=new TH1D ("fHistYieldStat","fHistYieldStat",300,0,30);
  TH1D* fHistYieldSist=new TH1D ("fHistYieldSist","fHistYieldSist",300,0,30);
  Double_t   mult[nummolt+1]={26.02, 20.02, 14.97, 10.69, 4.4};

  for(Int_t m=0; m<nummolt+1; m++){
  fHistYieldStat->SetBinContent(fHistYieldStat->FindBin(mult[m]),Yield[m]);
  fHistYieldStat->SetBinError(fHistYieldStat->FindBin(mult[m]),YieldErrStat[m]);
  fHistYieldSist->SetBinContent(fHistYieldSist->FindBin(mult[m]),Yield[m]);
  fHistYieldSist->SetBinError(fHistYieldSist->FindBin(mult[m]),YieldErrSist[m]);
  }

  TString titleYieldX="dN_{ch}/d#eta";
  canvasYield->cd();
  StyleHisto(fHistYield[m], 0, LimSupYield, Color[TypeAnalysis], 33, titleYieldX, titleYieldY, titleYield );
  fHistYieldStat->Draw("e");
  fHistYieldSist->SetFillStyle(0);
  fHistYieldSist->Draw("same e2");
  //fifth part: syst associated to bkg used for OOJ subtraction  (only for Xi and jet)

  cout << "\n\n going to write on file " << endl;   
  fileout->WriteTObject(canvasYield);
  fileout->WriteTObject(fHistYieldvsMultStat);
  fileout->WriteTObject(fHistYieldvsMultSist);
  fileout->WriteTObject(canvasPtSpectraAll);
  fileout->WriteTObject(canvasPtSpectra);
  fileout->WriteTObject(canvasPtSpectraFit);
  fileout->WriteTObject(canvasBarlow);
  fileout->WriteTObject(canvasPtSpectraRelError);
  fileout->WriteTObject(canvasPtSpectraRelErrorAll);
  for(Int_t m=0; m<nummolt+1; m++){
    fileout->WriteTObject(      fHistSpectrumSistAll[m]);
    fileout->WriteTObject(      fHistSpectrumStat[m]);
  }

  fileout->Close();
  cout << "starting from the file(s) "  << PathInDef<< endl;
  cout << " I have created the file " << stringout << endl;
}

