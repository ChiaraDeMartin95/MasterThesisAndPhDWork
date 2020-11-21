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
#include "TGraphAsymmErrors.h"
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

void PtSpectraBis(Int_t type=0,  Int_t TypeAnalysis=0,Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=0,   Int_t israp=0,TString year="1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"/"Run2DataRed_hXi"/*"2016kehjl_hK0s"*/,  TString Path1 ="_Jet0.75"/*"_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/,  Bool_t isEfficiency=0,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=0, Int_t PtBinning=1, TString FitFixed="Boltzmann"/*"Fermi-Dirac"*/,   Int_t sys=0, Bool_t ZeroYieldLowPt=0, Bool_t TwoFitFunctions=0){

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
  TString RegionTypeOld[3] = {"Jet", "Bulk", "All"};

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

  TF1 * lineat2= new TF1 ("lineat2", "pol0", 0,8);
  lineat2->FixParameter(0,2);
  lineat2->SetLineColor(kBlack);
  lineat2->SetLineWidth(1);
  TF1 * lineatm2= new TF1 ("lineatm2", "pol0", 0,8);
  lineatm2->FixParameter(0,-2);
  lineatm2->SetLineColor(1);
  lineatm2->SetLineWidth(1);

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
  if (ZeroYieldLowPt && TypeAnalysis==0) stringout +="_ZeroYLowPtJet";
  //  stringout += "_LowPtExtr0.5";
  if (TwoFitFunctions) stringout += "_TwoFitFunctions";
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  cout << " definition integral regions " << endl;
  //low and upper ranges for the fit****************************
  //************************************************************

  Double_t LowRange[nummolt+1]= {2, 2, 1, 1.5, 1.5, 1.5}; 
  Double_t UpRange[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeSpectrumPart[nummolt+1]= {0};
  Double_t UpRangeSpectrumPart[nummolt+1]= {0};

  //Double_t LowRangeJet[nummolt+1]= {2, 1.5, 1, 1.5, 1.5, 1}; 
  Double_t LowRangeJet[nummolt+1]= {1,1,1,1,1,1};
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
	    //LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1; for LowPtExtr0.5
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
	if (m==0 || m==1)	UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
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
    else  {
      LowRangeSpectrumPart[m] = LowRange[m];
      if (type==8)       LowRangeSpectrumPart[m] = LowRange[m];
    }
    //    else  LowRangeSpectrumPart[m] = 8;
    cout << " ok " << endl;
    for(Int_t v=0; v < numPtV0Max; v++){
      if (LowRangeSpectrumPart[m]<=NPtV0[v]) break;
      PtBinMin[m]++;
      cout << "LowRange " << LowRangeSpectrumPart[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
      cout << "PtBinMin " << PtBinMin[m] << endl;
    }

    UpRangeSpectrumPart[m] =8;
    if (type==8 && m==4 && (TypeAnalysis==0 || TypeAnalysis==1))    UpRangeSpectrumPart[m] =4;
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
  TString titleY=  "1/#Delta#eta #Delta#phi 1/N_{trigg} dN/dp_{T}";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 

  TLegend *legendPhi = new TLegend(0.6, 0.6, 0.9, 0.9);

  /*
  if (type==0){
    legendErrorAll = new TLegend(0.6, 0.6, 0.9, 0.9);
  }
  else if (type==8){
    legendErrorAll = new TLegend(0.6, 0.1, 0.9, 0.4);
  }
  */
  TLegend *legendError = new TLegend(0.6, 0.7, 0.9, 0.9);
  TLegend *legendError2 = new TLegend(0.6, 0.7, 0.9, 0.9);

  TH1F* fHistSpectrumStat[nummolt+1];
  TH1F* fHistSpectrumStatpol0[nummolt+1];
  TH1F* fHistSpectrumStatOOJSubDef[nummolt+1];
  TH1F* fHistSpectrumStatLeadTrackDef[nummolt+1];
  TH1F* fHistSpectrumStatMCChoiceDef[nummolt+1];
  TH1F* fHistSpectrumSist[nummolt+1];
  TH1F* fHistSpectrum_max[nummolt+1];
  TH1F* fHistSpectrum_min[nummolt+1];

  TH1F* fHistSpectrumStatDPhi[nummolt+1][numsysDPhi];
  TH1F* fHistSpectrumSistDPhi[nummolt+1];
  TH1F* fHistSpectrumSistOOJ[nummolt+1];
  TH1F* fHistSpectrumSistAll[nummolt+1];
  TH1F* fHistSpectrumSistLeadTrack[nummolt+1];
  TH1F* fHistSpectrumSistMCChoice[nummolt+1];

  //histos for relative uncertainty
  TH1F* fHistSpectrumStatRelError[nummolt+1]; //stat
  TH1F* fHistSpectrumSistRelError[nummolt+1]; //sist on the DeltaPhi projections
  TH1F* fHistSpectrumSistRelErrorDPhi[nummolt+1]; //sist assoc to choice of DeltaPhi interval
  TH1F* fHistSpectrumSistRelErrorAll[nummolt+1]; //total sist
  TH1F* fHistSpectrumSistRelErrorDCAz[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorSE[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorDeltaEta[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMB[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorOOBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorIBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorLeadTrack[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMCChoice[nummolt+1]; 

  TCanvas* canvasYield = new TCanvas ("canvasYield", "canvasYield", 1300, 800);
  TCanvas* canvasYieldErr = new TCanvas ("canvasYieldErr", "canvasYieldErr", 1300, 800);

  TCanvas* canvasPtSpectra = new TCanvas ("canvasPtSpectra", "canvasPtSpectra", 1300, 800);
  canvasPtSpectra->Divide(3,2);
  TCanvas* canvaspol0 = new TCanvas ("canvaspol0", "canvaspol0", 1300, 800);
  canvaspol0->Divide(3,2);
  TCanvas* canvasOOJSubDef = new TCanvas ("canvasOOJSubDef", "canvasOOJSubDef", 1300, 800);
  canvasOOJSubDef->Divide(3,2);
  TCanvas* canvasLeadTrackDef = new TCanvas ("canvasLeadTrackDef", "canvasLeadTrackDef", 1300, 800);
  canvasLeadTrackDef->Divide(3,2);
  TCanvas* canvasMCChoiceDef = new TCanvas ("canvasMCChoiceDef", "canvasMCChoiceDef", 1300, 800);
  canvasMCChoiceDef->Divide(3,2);

  TCanvas* canvasPtSpectraFit = new TCanvas ("canvasPtSpectraFit", "canvasPtSpectraFit", 1300, 800);
  canvasPtSpectraFit->Divide(3,2);

  TCanvas* canvasPtSpectraFitBis = new TCanvas ("canvasPtSpectraFitBis", "canvasPtSpectraFitBis", 1300, 800);
  canvasPtSpectraFitBis->Divide(3,2);

  TCanvas* canvasPtSpectraAll = new TCanvas ("canvasPtSpectraAll", "canvasPtSpectraAll", 1300, 800);
  canvasPtSpectraAll->Divide(3,2);
  TCanvas* canvasBarlow = new TCanvas ("canvasBarlow", "canvasBarlow", 1300, 800);
  canvasBarlow->Divide(3,2);
  TCanvas* canvasBarlowpol0 = new TCanvas ("canvasBarlowpol0", "canvasBarlowpol0", 1300, 800);
  canvasBarlowpol0->Divide(3,2);
  TCanvas* canvasBarlowOOJSubDef = new TCanvas ("canvasBarlowOOJSubDef", "canvasBarlowOOJSubDef", 1300, 800);
  canvasBarlowOOJSubDef->Divide(3,2);
  TCanvas* canvasBarlowLeadTrackDef = new TCanvas ("canvasBarlowLeadTrackDef", "canvasBarlowLeadTrackDef", 1300, 800);
  canvasBarlowLeadTrackDef->Divide(3,2);
  TCanvas* canvasBarlowMCChoiceDef = new TCanvas ("canvasBarlowMCChoiceDef", "canvasBarlowMCChoiceDef", 1300, 800);
  canvasBarlowMCChoiceDef->Divide(3,2);


  TCanvas* canvasPtSpectraRelError = new TCanvas ("canvasPtSpectraRelError", "canvasPtSpectraRelError", 1300, 800);
  canvasPtSpectraRelError->Divide(3,2);
  TCanvas* canvasPtSpectraRelErrorAll = new TCanvas ("canvasPtSpectraRelErrorAll", "canvasPtSpectraRelErrorAll", 1300, 800);
  canvasPtSpectraRelErrorAll->Divide(3,2);

  Bool_t BarlowPassed[nummolt+1][numsysDPhi]={0};
  Int_t BarlowSign[nummolt+1][numsysDPhi]={0};
  TH1F* hBarlowVar[nummolt+1][numsysDPhi];
  Float_t BarlowVar[nummolt+1][numPtV0][numsysDPhi]={0};
  Bool_t BarlowPassedpol0[nummolt+1]={0};
  Int_t BarlowSignpol0[nummolt+1]={0};
  TH1F* hBarlowVarpol0[nummolt+1];
  Float_t BarlowVarpol0[nummolt+1][numPtV0]={0};
  Bool_t BarlowPassedOOJSubDef[nummolt+1]={0};
  Int_t BarlowSignOOJSubDef[nummolt+1]={0};
  TH1F* hBarlowVarOOJSubDef[nummolt+1];
  Float_t BarlowVarOOJSubDef[nummolt+1][numPtV0]={0};
  Bool_t BarlowPassedLeadTrackDef[nummolt+1]={0};
  Int_t BarlowSignLeadTrackDef[nummolt+1]={0};
  TH1F* hBarlowVarLeadTrackDef[nummolt+1];
  Float_t BarlowVarLeadTrackDef[nummolt+1][numPtV0]={0};
  Bool_t BarlowPassedMCChoiceDef[nummolt+1]={0};
  Int_t BarlowSignMCChoiceDef[nummolt+1]={0};
  TH1F* hBarlowVarMCChoiceDef[nummolt+1];
  Float_t BarlowVarMCChoiceDef[nummolt+1][numPtV0]={0};

  Float_t ErrDPhi[nummolt+1][numPtV0]={0};

  Int_t ColorsysDPhi[12] = {807,  881, 867,909, 881,860,868,841,  418, 881, 7, 1};

  Float_t LimSupYield=0.1;
  if (type==8){
    if (TypeAnalysis==0)    LimSupYield=0.003;
    else   LimSupYield=0.02;
  }
  if (type==0){
    if (TypeAnalysis==0)    LimSupYield=0.045;
    else    LimSupYield=0.3;
  }

  Float_t LimInfYield=0.1;
  if (type==8){
    if (TypeAnalysis==0)    LimInfYield=10e-8;
    else   LimInfYield=10e-5;
  }
  if (type==0){
    if (TypeAnalysis==0)    LimInfYield=0.02;
    else    LimInfYield=10e-5;
  }

  Float_t LimSupYieldErr=0.04;
  if (type==0){
    if (TypeAnalysis==1) LimSupYieldErr=0.02;
    if (TypeAnalysis==2) LimSupYieldErr=0.02;
  }
  if (type==8){
    if (TypeAnalysis==0) LimSupYieldErr=0.45;
    if (TypeAnalysis==1) LimSupYieldErr=0.1;
    if (TypeAnalysis==2) LimSupYieldErr=0.1;
  }

  Float_t LimSup=0.2;
  Float_t LimInf=0;
  if (type==0){
    if (TypeAnalysis==0) LimSup =0.02;
    LimInf =10e-5;
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSup =0.001;
    else   if (TypeAnalysis==1) LimSup =0.01;
    else   if (TypeAnalysis==2) LimSup =0.01;
    LimInf =10e-10;
  }
  Float_t LimSupError=0.01;
  Float_t LimInfError=0;
  Float_t LimSupErrorLog=0.01;
  if (type==0) {
    if (TypeAnalysis==0){ LimSupError =0.3; LimSupErrorLog =0.5;}
    else   if (TypeAnalysis==1) LimSupError =0.15;
    else   if (TypeAnalysis==2) LimSupError =0.05;
    LimInfError=10e-5;
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSupError =1;
    else   if (TypeAnalysis==1) LimSupError =0.3;
    else   if (TypeAnalysis==2) LimSupError =0.1;
    LimInfError=10e-5;
  }

  TH1F*  hDeltaPhiLimit;
  Float_t  LowBinDPhi[numsysDPhi] ={0};
  Float_t   UpBinDPhi[numsysDPhi]={0}  ;

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
    fHistSpectrumSistRelErrorDCAz[m]    =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorDCAz_"+Smolt[m]); 
    fHistSpectrumSistRelErrorSE[m]      =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorSE_"+Smolt[m]);
    fHistSpectrumSistRelErrorDeltaEta[m]=(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorDeltaEta_"+Smolt[m]);
    if (!fHistSpectrumSistRelErrorDCAz[m]    ) return;
    if (!fHistSpectrumSistRelErrorSE[m]      ) return;
    if (!fHistSpectrumSistRelErrorDeltaEta[m]) return;

    fHistSpectrumStat[m]=(TH1F*)fileinDef->Get("fHistSpectrum_"+Smolt[m]);
    fHistSpectrumSist[m]=(TH1F*)fileinDef->Get("fHistSpectrumSist_"+Smolt[m]);
    hDeltaPhiLimit= (TH1F*) fileinDef-> Get("DeltaPhiLimit");
    if (!hDeltaPhiLimit) {cout << "DetaPhiLimit not there " << endl;return;}
    LowBinDPhi[0] =	hDeltaPhiLimit->GetBinContent(1);
    UpBinDPhi[0] =	hDeltaPhiLimit->GetBinContent(2);

    if (!fHistSpectrumStat[m]) {cout << "fHistSpectrumStat does not exist" << endl; return;}
    if (!fHistSpectrumSist[m]) {cout << "fHistSpectrumSist does not exist" << endl; return;}

    fHistSpectrumSistRelErrorAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorAll_"+Smolt[m]);
    fHistSpectrumStatRelError[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumStatRelError_"+Smolt[m]);
    fHistSpectrumSistRelErrorDPhi[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorDPhi_"+Smolt[m]);
    fHistSpectrumSistRelErrorMB[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMB_"+Smolt[m]);
    fHistSpectrumSistRelErrorOOBPU[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorOOBPU_"+Smolt[m]);
    fHistSpectrumSistRelErrorIBPU[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorIBPU_"+Smolt[m]);
    fHistSpectrumSistRelErrorLeadTrack[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorLeadTrack_"+Smolt[m]);
    fHistSpectrumSistRelErrorMCChoice[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMCChoice_"+Smolt[m]);

    fHistSpectrum_max[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMax_"+Smolt[m]);
    fHistSpectrum_min[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMin_"+Smolt[m]);
    fHistSpectrumSistDPhi[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDPhi_"+Smolt[m]);
    fHistSpectrumSistOOJ[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistOOJ_"+Smolt[m]);
    fHistSpectrumSistAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistAll_"+Smolt[m]);
    fHistSpectrumSistLeadTrack[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistLeadTrack_"+Smolt[m]);
    fHistSpectrumSistMCChoice[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistMCChoice_"+Smolt[m]);

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
  Int_t NSign=3;
  Float_t  Sigma2OOBPU[nummolt+1] ={0};
  Float_t  Sigma2IBPU[nummolt+1] ={0};
  Float_t  Sigma2MB[nummolt+1]={0};   
  Float_t OOBPU = 0.012; //relative syst. uncertainty on pT spectrum associated to OOB PileUp taken from Fiorella's paper and "averaged" over pT
  Float_t IBPU=0.02; //same but for Inbunch pileup
  Float_t MB = 0.01; //same but for material budget
  if (type==8){
    OOBPU = IBPU = MB = 0.02;
  }

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
	hDeltaPhiLimit= (TH1F*) filein-> Get("DeltaPhiLimit");
	if (!hDeltaPhiLimit) {cout << "DetaPhiLimit not there " << endl;return;}
	LowBinDPhi[sysDPhi] =	hDeltaPhiLimit->GetBinContent(1);
	UpBinDPhi[sysDPhi] =	hDeltaPhiLimit->GetBinContent(2);
	fHistSpectrumStatDPhi[m][sysDPhi]=(TH1F*)filein->Get("fHistSpectrum_"+Smolt[m]);
	fHistSpectrumStatDPhi[m][sysDPhi]->SetName(Form("fHistSpectrumStat_sysDPhi%i", sysDPhi)+Smolt[m]);

	StyleHisto(fHistSpectrumStatDPhi[m][sysDPhi], 0, LimSup, ColorsysDPhi[sysDPhi], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
	StyleHisto(fHistSpectrumStat[m], 0, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m],  0, 0, 0);

	canvasPtSpectraAll->cd(m+1);
	//      cout << "I'm drawing on the canvas " << endl;
	gPad->SetLeftMargin(0.15);
	if(m==0 && sysDPhi==1) legendPhi->AddEntry(fHistSpectrumStat[m],Form("%.2f-%.2f", LowBinDPhi[0], UpBinDPhi[0]) ,"l");
	if(m==0) legendPhi->AddEntry(fHistSpectrumStatDPhi[m][sysDPhi],Form("%.2f-%.2f", LowBinDPhi[sysDPhi], UpBinDPhi[sysDPhi]) ,"l");
	if (sysDPhi==1) fHistSpectrumStat[m]->Draw("");
	fHistSpectrumStatDPhi[m][sysDPhi]->Draw("same");
	if (sysDPhi==numsysDPhi-1) legendPhi->Draw("same");

	hBarlowVar[m][sysDPhi]=(TH1F*)    fHistSpectrumStatDPhi[m][sysDPhi]->Clone("fHistBarlowVar_"+Smolt[m]);
	for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	  if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
	  if (v==1 && TypeAnalysis==0 && type==8) continue;
	  if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
	  //	cout << "\nv: " << v << endl;
	  //	cout << "...histo...dphi " << sysDPhi << " " << fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1)<< endl;
	  BarlowVar[m][v][sysDPhi] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) - pow(fHistSpectrumStatDPhi[m][sysDPhi]->GetBinError(v+1),2)));
	  if (TMath::Abs(BarlowVar[m][v][sysDPhi])>2) BarlowSign[m][sysDPhi]++;
	  hBarlowVar[m][sysDPhi] ->SetBinContent(v+1, BarlowVar[m][v][sysDPhi]) ;
	  hBarlowVar[m][sysDPhi] ->SetBinError(v+1, 0) ;
	  cout << "barlow var v" << v << " "  << 	  hBarlowVar[m][sysDPhi] ->GetBinContent(v+1)<<" " <<	  hBarlowVar[m][sysDPhi] ->GetBinError(v+1) << endl;

	}//end loop v
	if (BarlowSign[m][sysDPhi]>= NSign) BarlowPassed[m][sysDPhi]=1;
	StyleHisto(hBarlowVar[m][sysDPhi], -5, 5, ColorsysDPhi[sysDPhi], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m],  0, 0, 0);
	if ( BarlowPassed[m][sysDPhi])       StyleHisto(hBarlowVar[m][sysDPhi], -5, 5, ColorsysDPhi[sysDPhi], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m],  0, 0, 0);

	cout << "Barlow passed " << BarlowPassed[m][sysDPhi] << endl;
	canvasBarlow->cd(m+1);
	gPad->SetLeftMargin(0.15);
	hBarlowVar[m][sysDPhi] ->Draw("same p");

	if (!BarlowPassed[m][sysDPhi]) continue;
	//      cout << "defining max and min for " << endl;
	for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	  if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
	  if (v==1 && TypeAnalysis==0 && type==8) continue;
	  if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
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
	if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
	if (v==1 && TypeAnalysis==0 && type==8) continue;
	if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
	//      cout << " v: " << v << endl;
	ErrDPhi[m][v]= TMath::Abs(	  fHistSpectrum_min[m]->GetBinContent(v+1) - 	  fHistSpectrum_max[m]->GetBinContent(v+1))/2;
	cout << fHistSpectrum_min[m]->GetBinContent(v+1)<< " << " << fHistSpectrumStat[m]->GetBinContent(v+1)<< " << " <<fHistSpectrum_max[m]->GetBinContent(v+1)<<endl;
	cout << "err dphi " << ErrDPhi[m][v] << endl;
	fHistSpectrumSistDPhi[m]->SetBinError(v+1,ErrDPhi[m][v]);

	if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1,ErrDPhi[m][v]/ fHistSpectrumStat[m]->GetBinContent(v+1));
	}
	//	else       fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1,0);
	fHistSpectrumSistRelErrorDPhi[m]->SetBinError(v+1,0);
      }
      //end new
    }//end loop on m
  } //end of type analysis

  //here syst related to pt < pt,trigg 
  //    cout << "\n*******************************************"<<endl;
  //    cout<< "hXi: systematic effect associated to selection pt,assoc < pt,Trigg  " << endl;
  Float_t   YieldSpectrumErrLeadTrack[nummolt+1]={0};
  TString    PathInLeadTrackDef;
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistLeadTrack[m]->SetBinError(v+1,0);
    }
  }
  if (type==8){
    cout << "\n*******************************************"<<endl;
    cout<< "hXi: systematic effect associated to selection pt,assoc < pt,Trigg  " << endl;
    PathInLeadTrackDef = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    //    PathInLeadTrackDef +=Path1; //it was without year                                                        
    if (type==8 && TypeAnalysis==0) PathInLeadTrackDef += "_OOJSmoothedBis";
    if(type>=0){
      PathInLeadTrackDef +="_"+tipo[type];
      PathInLeadTrackDef +=Srap[israp];
      PathInLeadTrackDef +="_AllAssoc";
    }
    //    if (IsHighPtExtr) PathInLeadTrackDef+="_HighPtExtr";
    PathInLeadTrackDef+= hhCorr[ishhCorr]+"_"+RegionTypeOld[TypeAnalysis] +DataOrMC[isMC] + Form("_PtMin%.1f.root", PtTrigMin);
    cout << "\n\n" << PathInLeadTrackDef << endl;
    TFile *  fileinLeadTrackDef = new TFile(PathInLeadTrackDef, "");
    if (!fileinLeadTrackDef) {cout << PathInLeadTrackDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummolt+1; m++){
      //    cout << " m " << m << endl;
      fHistSpectrumStatLeadTrackDef[m]    =(TH1F*)fileinLeadTrackDef->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      if (!fHistSpectrumStatLeadTrackDef[m]) {cout << " I was looking for histo in " << PathInLeadTrackDef  << endl; return;}

      hBarlowVarLeadTrackDef[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistBarlowVarLeadTrackDef_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (NPtV0[v] <3) continue;
	if (m== 4 && v==numPtV0Max-1 && TypeAnalysis==0 && type==8) continue;
	if (fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1) ==0) continue;
	BarlowVarLeadTrackDef[m][v] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) -pow(fHistSpectrumStatLeadTrackDef[m]->GetBinError(v+1),2)));
	if (TMath::Abs(BarlowVarLeadTrackDef[m][v])>2) BarlowSignLeadTrackDef[m]++;
	hBarlowVarLeadTrackDef[m] ->SetBinContent(v+1, BarlowVarLeadTrackDef[m][v]) ;
	hBarlowVarLeadTrackDef[m] ->SetBinError(v+1, 0) ;
	cout << "barlow var v" << v << " "  << 	  hBarlowVarLeadTrackDef[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarLeadTrackDef[m] ->GetBinError(v+1) << endl;

      }//end loop v

      if (BarlowSignLeadTrackDef[m]>= 2) BarlowPassedLeadTrackDef[m]=1;
      StyleHisto(hBarlowVarLeadTrackDef[m], -5, 5, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);
      if (BarlowPassedLeadTrackDef[m]) StyleHisto(hBarlowVarLeadTrackDef[m], -5, 5, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);

      cout << "Barlow passed " << BarlowPassedLeadTrackDef[m] << endl;
      canvasBarlowLeadTrackDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      hBarlowVarLeadTrackDef[m] ->Draw("same p");

      canvasLeadTrackDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStat[m]->Draw("ep");
      StyleHisto(fHistSpectrumStatLeadTrackDef[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStatLeadTrackDef[m]->Draw("same ep");

      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (NPtV0[v] <3 || TypeAnalysis==2){
	  fHistSpectrumSistLeadTrack[m]->SetBinError(v+1,0);
	  continue;
	}
	if (v==numPtV0Max-1 && m==4 && TypeAnalysis==0 && type==8) continue;
	if (fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1) ==0) continue;
	//	if(TMath::Abs(BarlowVarLeadTrackDef[m][v])>2){
	if (BarlowPassedLeadTrackDef[m]){
	  fHistSpectrumSistLeadTrack[m]->SetBinError(v+1,	  (fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1))/2);
	  //YieldSpectrumErrLeadTrack[m]+= TMath::Abs(   fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1));
	}
      }
    }//end loop m
  }//end loop syst associated to choice pt <pt,trigg

  //    cout << "\n*******************************************"<<endl;
  //    cout<< "systematic effect associated to MC used to calculate efficiency  " << endl;
  Float_t   YieldSpectrumErrMCChoice[nummolt+1]={0};
  TString    PathInMCChoiceDef;
  TLegend * legendMCChoice= new TLegend(0.6, 0.6, 0.9, 0.9);
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistMCChoice[m]->SetBinError(v+1,0);
    }
  }
  if (type==0){
    cout << "\n*******************************************"<<endl;
    cout<< "hXi: systematic effect associated to choice of MC used to calculat K0s/Xi efficiency  " << endl;
    PathInMCChoiceDef = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    //    PathInMCChoiceDef +=Path1; //it was without year                                                    
    if (PtBinning>0) PathInMCChoiceDef +=Form("_PtBinning%i",PtBinning);
    if (type==8 && TypeAnalysis==0) PathInMCChoiceDef += "_OOJSmoothedBis";
    if(type>=0){
      PathInMCChoiceDef +="_"+tipo[type];
      PathInMCChoiceDef +=Srap[israp];
      PathInMCChoiceDef +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInMCChoiceDef+="_HighPtExtr";
    PathInMCChoiceDef+= hhCorr[ishhCorr]+"_"+RegionTypeOld[TypeAnalysis] +DataOrMC[isMC] + Form("_PtMin%.1f_EPOS.root", PtTrigMin);
    cout << "\n\n" << PathInMCChoiceDef << endl;
    TFile *  fileinMCChoiceDef = new TFile(PathInMCChoiceDef, "");
    if (!fileinMCChoiceDef) {cout << PathInMCChoiceDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummolt+1; m++){
      //    cout << " m " << m << endl;
      fHistSpectrumStatMCChoiceDef[m]    =(TH1F*)fileinMCChoiceDef->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      if (!fHistSpectrumStatMCChoiceDef[m]) {cout << " I was looking for histo in " << PathInMCChoiceDef  << endl; return;}

      hBarlowVarMCChoiceDef[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistBarlowVarMCChoiceDef_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	//	if (NPtV0[v] <3) continue;
	if (m== 4 && v==numPtV0Max-1 && TypeAnalysis==0 && type==8) continue;
	if (fHistSpectrumStatMCChoiceDef[m]->GetBinContent(v+1) ==0) continue;
	BarlowVarMCChoiceDef[m][v] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatMCChoiceDef[m]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) -pow(fHistSpectrumStatMCChoiceDef[m]->GetBinError(v+1),2)));
	if (TMath::Abs(BarlowVarMCChoiceDef[m][v])>2) BarlowSignMCChoiceDef[m]++;
	hBarlowVarMCChoiceDef[m] ->SetBinContent(v+1, BarlowVarMCChoiceDef[m][v]) ;
	hBarlowVarMCChoiceDef[m] ->SetBinError(v+1, 0) ;
	cout << "barlow var v" << v << " "  << 	  hBarlowVarMCChoiceDef[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarMCChoiceDef[m] ->GetBinError(v+1) << endl;

      }//end loop v

      if (BarlowSignMCChoiceDef[m]>= 3) BarlowPassedMCChoiceDef[m]=1;
      StyleHisto(hBarlowVarMCChoiceDef[m], -5, 5, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);
      if (BarlowPassedMCChoiceDef[m]) StyleHisto(hBarlowVarMCChoiceDef[m], -5, 5, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);

      cout << "Barlow passed " << BarlowPassedMCChoiceDef[m] << endl;
      canvasBarlowMCChoiceDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      hBarlowVarMCChoiceDef[m] ->Draw("same p");
      lineat2->Draw("same");
      lineatm2->Draw("same");

      canvasMCChoiceDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendMCChoice->AddEntry(fHistSpectrumStat[m], "default", "ple");
      fHistSpectrumStat[m]->Draw("ep");
      StyleHisto(fHistSpectrumStatMCChoiceDef[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendMCChoice->AddEntry(fHistSpectrumStatMCChoiceDef[m], "EPOS", "ple");
      fHistSpectrumStatMCChoiceDef[m]->Draw("same ep");
      legendMCChoice->Draw("");

      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistMCChoice[m]->SetBinError(v+1,0);
	if (v==numPtV0Max-1 && m==4 && TypeAnalysis==0 && type==8) continue;
	if (fHistSpectrumStatMCChoiceDef[m]->GetBinContent(v+1) ==0) continue;
	//	if(TMath::Abs(BarlowVarMCChoiceDef[m][v])>2){
	if (BarlowPassedMCChoiceDef[m]){
	  fHistSpectrumSistMCChoice[m]->SetBinError(v+1,	  (fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatMCChoiceDef[m]->GetBinContent(v+1))/2);
	  //YieldSpectrumErrMCChoice[m]+= TMath::Abs(   fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatMCChoiceDef[m]->GetBinContent(v+1));
	}
      }
    }//end loop m
  }//end loop syst associated to MC used to calculate efficiency


  //end pt < pt,trig

  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      Sigma2OOBPU[m] =pow(OOBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
      Sigma2IBPU[m] =pow(IBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
      Sigma2MB[m]= pow(MB *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
      fHistSpectrumSistAll[m]->SetBinError(v+1,sqrt(pow(fHistSpectrumSistDPhi[m]->GetBinError(v+1),2) + pow(fHistSpectrumSist[m]->GetBinError(v+1),2) + Sigma2OOBPU[m]+ Sigma2IBPU[m] + Sigma2MB[m] + pow(fHistSpectrumSistLeadTrack[m]->GetBinError(v+1),2)+  pow(fHistSpectrumSistMCChoice[m]->GetBinError(v+1),2)));

      if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSistAll[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumSistRelErrorMB[m]->SetBinContent(v+1,MB);
	fHistSpectrumSistRelErrorIBPU[m]->SetBinContent(v+1,IBPU);
	fHistSpectrumSistRelErrorOOBPU[m]->SetBinContent(v+1,OOBPU);
	fHistSpectrumSistRelErrorLeadTrack[m]->SetBinContent(v+1,fHistSpectrumSistLeadTrack[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumSistRelErrorMCChoice[m]->SetBinContent(v+1,fHistSpectrumSistMCChoice[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      //	else       fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1,0);
      fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMB[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorIBPU[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorOOBPU[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorLeadTrack[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMCChoice[m]->SetBinError(v+1,0);
      //      cout << " stat rel error " <<       fHistSpectrumStatRelError[m]->GetBinContent(v+1)<< endl;
      //      cout << " sist rel error assoc to DeltaPhi choice " <<       fHistSpectrumSistRelErrorDPhi[m]->GetBinContent(v+1)<< endl;
    }

    //    cout << "drawing sist rel error " << endl;
    canvasPtSpectraRelErrorAll->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumStatRelError[m], LimInfError, LimSupError, Color[TypeAnalysis], 33, titleX, titleYRel ,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorAll[m], LimInfError, LimSupError, Color[TypeAnalysis], 27, titleX, titleYRel ,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorSE[m], LimInfError, LimSupError, 807, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorDCAz[m], LimInfError, LimSupError, 867, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorDeltaEta[m], LimInfError, LimSupError, 881, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorDPhi[m], LimInfError, LimSupError, 922, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorMB[m], LimInfError, LimSupError, 401, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorOOBPU[m], LimInfError, LimSupError, 825, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorIBPU[m], LimInfError, LimSupError, 631, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorLeadTrack[m], LimInfError, LimSupError, 630, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorMCChoice[m], LimInfError, LimSupError, 907, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    gPad->SetLogy();

    TLegend *legendErrorAll= new TLegend(0.6, 0.1, 0.9, 0.4);
    legendErrorAll->AddEntry(fHistSpectrumStatRelError[m], "stat.", "pl");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorSE[m], "syst. signal extr.", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDCAz[m], "syst. DCAzTrigger", "l");
    if(TypeAnalysis!=2)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDeltaEta[m], "syst. #Delta#eta region", "l");
    if(TypeAnalysis!=2)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDPhi[m], "syst. #Delta#phi region", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMB[m], "Material Budget", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorOOBPU[m], "Out-of-bunch PU", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorIBPU[m], "In-bunch PU", "l");
    if(BarlowPassedLeadTrackDef[m])   legendErrorAll->AddEntry(fHistSpectrumSistRelErrorLeadTrack[m], "p_{T} < p_{T}^{Trigg}", "l");
    if(BarlowPassedMCChoiceDef[m])   legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMCChoice[m], "MCChoice", "l");

    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumSistRelErrorMB[m]->Draw("same");
    fHistSpectrumSistRelErrorIBPU[m]->Draw("same");
    fHistSpectrumSistRelErrorOOBPU[m]->Draw("same");
    fHistSpectrumSistRelErrorSE[m]->Draw("same");
    fHistSpectrumSistRelErrorDCAz[m]->Draw("same");
    fHistSpectrumSistRelErrorLeadTrack[m]->Draw("same");
    fHistSpectrumSistRelErrorMCChoice[m]->Draw("same");
    if (TypeAnalysis!=2) {
      fHistSpectrumSistRelErrorDeltaEta[m]->Draw("same");
      fHistSpectrumSistRelErrorDPhi[m]->Draw("same");
    }
    fHistSpectrumStatRelError[m]->Draw("same p");
    legendErrorAll->Draw("");
    /*
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      cout <<"m " << m << "v " << v << " " << NPtV0[v] <<   " DCAz " << fHistSpectrumSistRelErrorDCAz[m]->GetBinContent(v+1) << "signal extr " <<  fHistSpectrumSistRelErrorSE[m]->GetBinContent(v+1) << " DeltaEta " <<  fHistSpectrumSistRelErrorDeltaEta[m]  ->GetBinContent(v+1) <<  " DPhi " << fHistSpectrumSistRelErrorDPhi[m]->GetBinContent(v+1)<<endl;
      cout << " spectrum sist " << fHistSpectrumSist[m]  ->GetBinError(v+1)/fHistSpectrumSist[m]  ->GetBinContent(v+1)<< endl;
      cout << " sum in quadrature of above values + OOB/IB/MB" << sqrt(pow(fHistSpectrumSistRelErrorDCAz[m]->GetBinContent(v+1),2) + pow(fHistSpectrumSistRelErrorSE[m]->GetBinContent(v+1),2) +pow(fHistSpectrumSistRelErrorDeltaEta[m]->GetBinContent(v+1),2) + pow(fHistSpectrumSistRelErrorDPhi[m]->GetBinContent(v+1),2) + pow(OOBPU,2) + pow(IBPU,2) +pow(MB,2))<< endl;
      cout << " total " << fHistSpectrumSistRelErrorAll[m]->GetBinContent(v+1) << endl;         
    }
    */
  }
	
  //third part: for the jet, I evaluate sys associated to type of OOJ subtraction (for the jet)
  TString    PathInpol0;
  if (TypeAnalysis==0 && type==0){
    cout << "\n*******************************************"<<endl;
    cout << "3. sysy associated to pol0 OOJ subtraction for hK0s"<< endl;

    PathInpol0 = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    PathInpol0 +=Path1; //it was without year                                                        
    if(type>=0){
      PathInpol0 +="_"+tipo[type];
      PathInpol0 +=Srap[israp];
      PathInpol0 +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInpol0+="_HighPtExtr";
    PathInpol0+= hhCorr[ishhCorr]+"_JetBFit" +DataOrMC[isMC] + Form("_PtMin%.1f_Try1.root", PtTrigMin);
    cout << "\n\n" << PathInpol0 << endl;
    TFile *  fileinpol0 = new TFile(PathInpol0, "");
    if (!fileinpol0) {cout << PathInpol0 << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummolt+1; m++){
      //    cout << " m " << m << endl;
      fHistSpectrumStatpol0[m]    =(TH1F*)fileinpol0->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
    if (!fHistSpectrumStatpol0[m]) {cout << " I was looking for histo in " << PathInpol0  << endl; return;}

    hBarlowVarpol0[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistBarlowVarpol0_"+Smolt[m]);
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
      if (v==1 && TypeAnalysis==0 && type==8) continue;
      if (NPtV0[v] = 4 && TypeAnalysis==0 && type==8 && m==4) continue;
      BarlowVarpol0[m][v] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatpol0[m]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) -pow(fHistSpectrumStatpol0[m]->GetBinError(v+1),2)));
      if (TMath::Abs(BarlowVarpol0[m][v])>2) BarlowSignpol0[m]++;
      hBarlowVarpol0[m] ->SetBinContent(v+1, BarlowVarpol0[m][v]) ;
      hBarlowVarpol0[m] ->SetBinError(v+1, 0) ;
      cout << "barlow var v" << v << " "  << 	  hBarlowVarpol0[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarpol0[m] ->GetBinError(v+1) << endl;

    }//end loop v

    if (BarlowSignpol0[m]>= NSign) BarlowPassedpol0[m]=1;
    StyleHisto(hBarlowVarpol0[m], -5, 5, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0,0,0);
    if (BarlowPassedpol0[m]) StyleHisto(hBarlowVarpol0[m], -5, 5, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0,0,0);

    cout << "Barlow passed " << BarlowPassedpol0[m] << endl;
    canvasBarlowpol0->cd(m+1);
    gPad->SetLeftMargin(0.15);
    hBarlowVarpol0[m] ->Draw("same p");

    canvaspol0->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m],0,0,0);
    fHistSpectrumStat[m]->Draw("ep");
    StyleHisto(fHistSpectrumStatpol0[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0,0,0);
    fHistSpectrumStatpol0[m]->Draw("same ep");

    }
  }

  //I plot spectra with stat + total syst uncertainty
  // I plot stat + syst (total, DeltaPhi assoc., OOJ sub assoc.) relative uncertainty

  for(Int_t m=0; m<nummolt+1; m++){

    canvasPtSpectra->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistAll[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
    fHistSpectrumStat[m]->Draw("same");
    fHistSpectrumSistAll[m]->SetFillStyle(0);
    fHistSpectrumSistAll[m]->Draw("same e2");
    if (m==0){
    legendError2->AddEntry(fHistSpectrumStat[m], "stat.", "ple");
    legendError2->AddEntry(fHistSpectrumSist[m], "syst.", "fe");
    }
    legendError2->Draw("");

    canvasPtSpectraRelErrorAll->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistRelErrorAll[m], LimInfError, LimSupError,Color[TypeAnalysis] , 27, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumStatRelError[m], LimInfError, LimSupError,Color[TypeAnalysis] , 33, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);
    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");

    canvasPtSpectraRelError->cd(m+1);
    gPad->SetLeftMargin(0.15);
    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");
    if (m==0){
      legendError->AddEntry(    fHistSpectrumStatRelError[m], "stat.", "pl");
      legendError->AddEntry(    fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
    }
    legendError->Draw("");
  } //end loop m

  //*************************************************************************
  //for Xi in-jet production: systematic effect associated to ooj subtraction
  //*************************************************************************
  Float_t   YieldSpectrumErrOOJSub[nummolt+1]={0};
  TString    PathInOOJSubDef;
  if (TypeAnalysis==0 && type==8){
    cout << "\n*******************************************"<<endl;
    cout<< "hXi: systematic effect associated to OOJ subtraction " << endl;
    PathInOOJSubDef = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    PathInOOJSubDef +=Path1; //it was without year                                                        
    //PathInOOJSubDef +="_Jet0.75"; //it was without year                                                        
    if(type>=0){
      PathInOOJSubDef +="_"+tipo[type];
      PathInOOJSubDef +=Srap[israp];
      PathInOOJSubDef +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInOOJSubDef+="_HighPtExtr";
    PathInOOJSubDef+= hhCorr[ishhCorr]+"_Jet" +DataOrMC[isMC] + Form("_PtMin%.1f.root", PtTrigMin);
    cout << "\n\n" << PathInOOJSubDef << endl;
    TFile *  fileinOOJSubDef = new TFile(PathInOOJSubDef, "");
    if (!fileinOOJSubDef) {cout << PathInOOJSubDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummolt+1; m++){
      //    cout << " m " << m << endl;
      fHistSpectrumStatOOJSubDef[m]    =(TH1F*)fileinOOJSubDef->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      if (!fHistSpectrumStatOOJSubDef[m]) {cout << " I was looking for histo in " << PathInOOJSubDef  << endl; return;}

      hBarlowVarOOJSubDef[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistBarlowVarOOJSubDef_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (NPtV0[v] >2.5) continue;
	if (v==1 && TypeAnalysis==0 && type==8) continue;
	if (fHistSpectrumStatOOJSubDef[m]->GetBinContent(v+1) ==0) continue;
	BarlowVarOOJSubDef[m][v] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatOOJSubDef[m]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) -pow(fHistSpectrumStatOOJSubDef[m]->GetBinError(v+1),2)));
	if (TMath::Abs(BarlowVarOOJSubDef[m][v])>2) BarlowSignOOJSubDef[m]++;
	hBarlowVarOOJSubDef[m] ->SetBinContent(v+1, BarlowVarOOJSubDef[m][v]) ;
	hBarlowVarOOJSubDef[m] ->SetBinError(v+1, 0) ;
	cout << "barlow var v" << v << " "  << 	  hBarlowVarOOJSubDef[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarOOJSubDef[m] ->GetBinError(v+1) << endl;

      }//end loop v

      if (BarlowSignOOJSubDef[m]>= 1) BarlowPassedOOJSubDef[m]=1;
      StyleHisto(hBarlowVarOOJSubDef[m], -5, 5, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);
      if (BarlowPassedOOJSubDef[m]) StyleHisto(hBarlowVarOOJSubDef[m], -5, 5, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);

      cout << "Barlow passed " << BarlowPassedOOJSubDef[m] << endl;
      canvasBarlowOOJSubDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      hBarlowVarOOJSubDef[m] ->Draw("same p");

      canvasOOJSubDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStat[m]->Draw("ep");
      StyleHisto(fHistSpectrumStatOOJSubDef[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStatOOJSubDef[m]->Draw("same ep");

      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (NPtV0[v] >2.5) continue;
	if (v==1 && TypeAnalysis==0 && type==8) continue;
	if (fHistSpectrumStatOOJSubDef[m]->GetBinContent(v+1) ==0) continue;
	if(TMath::Abs(BarlowVarOOJSubDef[m][v])>2){
	  YieldSpectrumErrOOJSub[m]+= TMath::Abs(   fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatOOJSubDef[m]->GetBinContent(v+1));
	}
      }
    }//end loop m
  }//end loop ooj sub for hXi


  //fourth part: fit to obtain pt-integrated yield vs mult
  AliPWGFunc pwgfunc;
  const Int_t numfittipo=4;
  TLegend *legendfit=new TLegend(0.6, 0.6, 0.9, 0.9);
  TString   nameMTscaling[nummolt+1][numfittipo];
  Int_t ColorFit[numfittipo+1]={860, 881, 868, 628, 419};
  TFitResultPtr fFitResultPtr0[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtr1[nummolt+1][numfittipo];
  TF1* fit_MTscaling[nummolt+1][numfittipo];
  TF1* fit_MTscalingBis[nummolt+1][numfittipo];
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
  Float_t    YieldExtrMaxLowPt[nummolt+1]={0};
  Float_t    YieldExtrMinLowPt[nummolt+1]={0};
  Float_t    YieldExtrMaxHighPt[nummolt+1]={0};
  Float_t    YieldExtrMinHighPt[nummolt+1]={0};
  Float_t    YieldExtrErrStat[nummolt+1]={0};
  Float_t    YieldExtrErrSist[nummolt+1]={0};
  Float_t    YieldExtrErrSistLowPt[nummolt+1]={0};
  Float_t    YieldExtrErrSistHighPt[nummolt+1]={0};
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

      fit_MTscalingBis[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Bis");

      fit_MTscaling[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fit_MTscalingBis[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      if (TwoFitFunctions) {
	fit_MTscaling[m][typefit]->SetRange(LowRange[m], 1.5);
	fit_MTscalingBis[m][typefit]->SetRange(2, UpRange[m]);
      }
      fFitResultPtr0[m][typefit]=       fHistSpectrumStat[m]->Fit(    fit_MTscaling[m][typefit],"SR0");
      fit_MTscaling[m][typefit]->SetRange(0,50);
      canvasPtSpectraFit->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStat[m]->Draw("same");
      fit_MTscaling[m][typefit]->Draw("same");
      legendfit->Draw("");

      fFitResultPtr1[m][typefit]=       fHistSpectrumStat[m]->Fit(    fit_MTscalingBis[m][typefit],"SR0");
      fit_MTscalingBis[m][typefit]->SetRange(0,50);
      canvasPtSpectraFitBis->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStat[m]->Draw("same");
      fit_MTscalingBis[m][typefit]->Draw("same");
      legendfit->Draw("");

      //calculatin yields
      YieldExtrLowPt[m][typefit]= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPt[m][typefit] = fit_MTscalingBis[m][typefit]->Integral(UpRangeSpectrumPart[m],50);

      if (typefit==0){
	YieldExtrMaxLowPt[m] =  YieldExtrLowPt[m][typefit];
	YieldExtrMaxHighPt[m] =  YieldExtrHighPt[m][typefit];
	YieldExtrMinLowPt[m] =  YieldExtrLowPt[m][typefit];
	YieldExtrMinHighPt[m] =  YieldExtrHighPt[m][typefit];

      }
      if ((YieldExtrLowPt[m][typefit]) > YieldExtrMaxLowPt[m]) {
	YieldExtrMaxLowPt[m] = YieldExtrLowPt[m][typefit];
      }
      if ((YieldExtrLowPt[m][typefit]) < YieldExtrMinLowPt[m]) {
	YieldExtrMinLowPt[m] = YieldExtrLowPt[m][typefit];
      }

      if ((YieldExtrHighPt[m][typefit]) > YieldExtrMaxHighPt[m]) {
	YieldExtrMaxHighPt[m] = YieldExtrHighPt[m][typefit];
      }
      if ((YieldExtrHighPt[m][typefit]) < YieldExtrMinHighPt[m]) {
	YieldExtrMinHighPt[m] = YieldExtrHighPt[m][typefit];
      }

      YieldExtrLowPtAvg[m]+= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvg[m]+= fit_MTscalingBis[m][typefit]->Integral(UpRangeSpectrumPart[m],50);     

      YieldErrStatLowPtAvg[m]+= pow(    fit_MTscaling[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
      YieldErrStatHighPtAvg[m]+= pow(    fit_MTscalingBis[m][typefit]->IntegralError(UpRangeSpectrumPart[m],50,fFitResultPtr1[m][typefit] ->GetParams(),(fFitResultPtr1[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);

    } //end loop typefit


    YieldExtrHighPtAvg[m]=YieldExtrHighPtAvg[m]/numfittipo;
    YieldExtrLowPtAvg[m]=YieldExtrLowPtAvg[m]/numfittipo;
    YieldErrStatHighPtAvg[m]=sqrt(YieldErrStatHighPtAvg[m])/numfittipo;
    YieldErrStatLowPtAvg[m]=sqrt(YieldErrStatLowPtAvg[m])/numfittipo;

    YieldExtr[m] =     YieldExtrHighPtAvg[m]+    YieldExtrLowPtAvg[m];
    YieldExtrErrStat[m] = sqrt(    pow(YieldErrStatHighPtAvg[m],2)+    pow(YieldErrStatLowPtAvg[m],2));
    YieldExtrErrSistHighPt[m] = (YieldExtrMaxHighPt[m]-YieldExtrMinHighPt[m])/2;
    YieldExtrErrSistLowPt[m] = (YieldExtrMaxLowPt[m]-YieldExtrMinLowPt[m])/2;
    YieldExtrErrSistUp[m] = YieldExtrMaxLowPt[m]-YieldExtrLowPtAvg[m];
    YieldExtrErrSistLow[m] = YieldExtrLowPtAvg[m];
    //    YieldExtrErrSist[m]=sqrt(pow(    YieldExtrErrSistUp[m],2) + pow(    YieldExtrErrSistLow[m],2)) ;
    YieldExtrErrSist[m]=sqrt(pow(    YieldExtrErrSistHighPt[m],2) + pow(    YieldExtrErrSistLowPt[m],2)) ;

    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      YieldSpectrum[m] +=  ( fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1));
      YieldSpectrumErrStat[m] += pow ( fHistSpectrumStat[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1),2);
      YieldSpectrumErrSist[m] += pow ( fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1),2);
     
    }
    YieldSpectrumErrStat[m]=sqrt( YieldSpectrumErrStat[m]);
    YieldSpectrumErrSist[m]=sqrt( YieldSpectrumErrSist[m]);

    Yield[m] =  YieldSpectrum[m]+YieldExtr[m];
    YieldErrStat[m] = sqrt(    pow(YieldSpectrumErrStat[m],2) + pow(YieldExtrErrStat[m],2));
    YieldErrSist[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistHighPt[m],2) +  pow(YieldExtrErrSistLowPt[m],2) + pow(YieldSpectrumErrOOJSub[m],2));
    YieldErrSistUp[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistUp[m],2) +  pow(YieldExtrErrSistHighPt[m],2)  + pow(YieldSpectrumErrOOJSub[m],2));
    YieldErrSistLow[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistLow[m],2) +  pow(YieldExtrErrSistHighPt[m],2)  + pow(YieldSpectrumErrOOJSub[m],2));

    if (ZeroYieldLowPt){
      Yield[m] =  YieldSpectrum[m]+ YieldExtrHighPtAvg[m];
      YieldErrStat[m] = sqrt(    pow(YieldSpectrumErrStat[m],2) + pow(YieldErrStatHighPtAvg[m],2));
      YieldErrSist[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistHighPt[m],2) + pow(YieldSpectrumErrOOJSub[m],2));
    }

    cout << " m " << m << endl;
    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      cout << " v " << v << " " << NPtV0[v] << " yield: " << fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1) << endl;
    }
    cout << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) +- "<< YieldErrSist[m] << " (sist) "<< endl;
    cout << " (stat rel: " <<  YieldErrStat[m]/ Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldErrSist[m]/ Yield[m] << ")= " << endl;
    cout << "from spectrum: " <<   YieldSpectrum[m] << "+- " << YieldSpectrumErrStat[m]<<" (stat) +- " <<YieldSpectrumErrSist[m]<< " (sist)"  << endl;
    if (type==8 && TypeAnalysis==0)    cout << "from OOJ Sub " <<      YieldSpectrumErrOOJSub[m] << endl;
    cout << "from extr: " <<     YieldExtr[m] << "+- " <<     YieldExtrErrStat[m]<< " (stat)  +- " << YieldExtrErrSist[m]<< " (sist) " << endl;
    cout << " extr at low pt " << YieldExtrLowPtAvg[m]<<" " << "err extrap high pt " << YieldExtrHighPtAvg[m]<< endl; 
    cout << " extr at low pt " << YieldExtrLowPtAvg[m]<<" " << "err extrap low pt " << YieldExtrErrSistLow[m]<< endl; 
    if (TypeAnalysis==0 && type==8)     cout << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) + "<< YieldErrSistUp[m]<< " - " << YieldErrSistLow[m] << " (sist) "<< endl;

    cout << " fraction of yield from low pt extrapolation " <<  YieldExtrLowPtAvg[m]/Yield[m] << endl;
    cout << " fraction of yield from high pt extrapolation " <<  YieldExtrHighPtAvg[m]/Yield[m] << endl;
    cout << "fit range " << LowRange[m] << "- " << UpRange[m] << endl;
  }//end loop m


  TLegend *legendYield=new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend *legendYieldErr=new TLegend(0.6, 0.6, 0.9, 0.9);

  TF1 * pol0 = new TF1 ("pol0", "pol0", 0, 30);
  TString titleYieldX="dN_{ch}/d#eta";
  TString titleYieldYType[2]={"N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi", "N_{#Xi}/N_{Trigg} 1/#Delta#eta #Delta#phi"};
  TString titleYieldY;
  if(type==0) titleYieldY=titleYieldYType[0];
  else if(type==8) titleYieldY=titleYieldYType[1];
  TString titleYield[3]={"In-jet", "Out-of-jet", "Inclusive"};
\
  TH1F* fHistYieldStat=new TH1F ("fHistYieldStat","fHistYieldStat",150,0,30);
  TH1F* fHistYieldSist=new TH1F ("fHistYieldSist","fHistYieldSist",150,0,30);
  TH1F* fHistYieldSistNoExtr=new TH1F ("fHistYieldSistNoExtr","fHistYieldSistNoExtr",150,0,30);
  TH1F* fHistYieldStatRelErr=new TH1F ("fHistYieldStatRelErr","fHistYieldStatRelErr",150,0,30);
  TH1F* fHistYieldSistRelErr=new TH1F ("fHistYieldSistRelErr","fHistYieldSistRelErr",150,0,30);
  TH1F* fHistYieldSistNoExtrRelErr=new TH1F ("fHistYieldSistNoExtrRelErr","fHistYieldSistNoExtrRelErr",150,0,30);
  TH1F* fHistYieldSistLowRelErr=new TH1F ("fHistYieldSistLowRelErr","fHistYieldSistLowRelErr",150,0,30);
  TH1F* fHistYieldSistUpRelErr=new TH1F ("fHistYieldSistUpRelErr","fHistYieldSistUpRelErr",150,0,30);
  Float_t   mult[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};

  for(Int_t m=0; m<nummolt; m++){
  fHistYieldStat->SetBinContent(fHistYieldStat->FindBin(mult[m]),Yield[m]);
  fHistYieldStat->SetBinError(fHistYieldStat->FindBin(mult[m]),YieldErrStat[m]);
  fHistYieldSist->SetBinContent(fHistYieldSist->FindBin(mult[m]),Yield[m]);
  fHistYieldSist->SetBinError(fHistYieldSist->FindBin(mult[m]),YieldErrSist[m]);
  fHistYieldSistNoExtr->SetBinContent(fHistYieldSist->FindBin(mult[m]),Yield[m]);
  fHistYieldSistNoExtr->SetBinError(fHistYieldSist->FindBin(mult[m]), YieldSpectrumErrSist[m]);

  fHistYieldStatRelErr->SetBinContent(fHistYieldStat->FindBin(mult[m]),YieldErrStat[m]/Yield[m]);
  fHistYieldStatRelErr->SetBinError(fHistYieldStat->FindBin(mult[m]),0);
  fHistYieldSistRelErr->SetBinContent(fHistYieldSist->FindBin(mult[m]),YieldErrSist[m]/Yield[m]);
  fHistYieldSistRelErr->SetBinError(fHistYieldSist->FindBin(mult[m]),0);
  fHistYieldSistNoExtrRelErr->SetBinContent(fHistYieldSist->FindBin(mult[m]),YieldSpectrumErrSist[m]/Yield[m]);
  fHistYieldSistNoExtrRelErr->SetBinError(fHistYieldSist->FindBin(mult[m]),0);
  fHistYieldSistLowRelErr->SetBinContent(fHistYieldSist->FindBin(mult[m]),YieldErrSistLow[m]/Yield[m]);
  fHistYieldSistLowRelErr->SetBinError(fHistYieldSist->FindBin(mult[m]),0);
  fHistYieldSistUpRelErr->SetBinContent(fHistYieldSist->FindBin(mult[m]),YieldErrSistUp[m]/Yield[m]);
  fHistYieldSistUpRelErr->SetBinError(fHistYieldSist->FindBin(mult[m]),0);

  fHistYieldStat->SetMarkerSize(1.5);
  fHistYieldSist->SetMarkerSize(1.5);
  fHistYieldSistNoExtr->SetMarkerSize(1.5);
  fHistYieldStatRelErr->SetMarkerSize(2);
  fHistYieldSistRelErr->SetMarkerSize(2);
  fHistYieldSistNoExtrRelErr->SetMarkerSize(2);
  fHistYieldSistLowRelErr->SetMarkerSize(2);
  fHistYieldSistUpRelErr->SetMarkerSize(2);

  }

  //  TGraphAsymmErrors* gYield = new TGraphAsymmErrors(nummolt+1,mult,Yield,0.,0.,YieldErrSistLow,YieldErrSistUp);
  Float_t Xl[nummolt+1]= {0};
  Float_t Xh[nummolt+1]= {0};
  Float_t multctrbin[nummolt+1] ={0} ;
  for (Int_t m=0; m<nummolt+1; m++){
    multctrbin[m] =   fHistYieldStat->GetXaxis()->GetBinCenter(  fHistYieldStat->FindBin(mult[m]));
  }

  TGraphAsymmErrors* gYield = new TGraphAsymmErrors(nummolt,multctrbin,Yield,Xl, Xh, YieldErrSistLow, YieldErrSistUp);
  gYield->SetTitle(titleYield[TypeAnalysis]+" yield vs multiplicity");
  //gYield->SetMarkerColor(Color[TypeAnalysis]);
  gYield->SetMarkerColor(kBlue);
  //  gYield->SetLineColor(Color[TypeAnalysis]);
  gYield->SetLineColor(kBlue);
  gYield->SetMarkerStyle(1);

  canvasYield->cd();
  StyleHisto(fHistYieldStat, LimInfYield, LimSupYield, Color[TypeAnalysis], 1, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity",1, 0, 25);
  StyleHisto(fHistYieldSist, LimInfYield, LimSupYield, Color[TypeAnalysis], 1, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity", 1, 0, 25);
  StyleHisto(fHistYieldSistNoExtr, LimInfYield, LimSupYield, Color[TypeAnalysis], 1, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity", 1, 0, 25);
  fHistYieldStat->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldSist->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldSistNoExtr->GetYaxis()->SetTitleOffset(1.2);
  if (TypeAnalysis==0){
    fHistYieldStat->Fit(pol0, "R0");
  }

  fHistYieldStat->Draw("e");
  //  if (TypeAnalysis==0 && type==8)  gYield->Draw("same p");
  //  fHistYieldSistNoExtr->SetFillStyle(9);
  fHistYieldSistNoExtr->SetFillStyle(3001);
  fHistYieldSistNoExtr->SetFillColorAlpha(Color[TypeAnalysis], 1);
  fHistYieldSistNoExtr->Draw("same e2");
  fHistYieldSist->SetFillStyle(0);
  //  if (!(TypeAnalysis==0 && type==8))  fHistYieldSist->Draw("same e2");
  fHistYieldSist->Draw("same e2");
  gYield->SetFillStyle(0);
  legendYield->AddEntry(fHistYieldSist, "syst.", "ef");
  legendYield->AddEntry(fHistYieldSistNoExtr, "syst. (no extr)", "ef");
  legendYield->AddEntry(fHistYieldStat, "stat.", "pel");
  legendYield->Draw("");
  canvasYield->SaveAs("provay.pdf");

  canvasYieldErr->cd();
  StyleHisto(fHistYieldStatRelErr, 10e-5, LimSupYieldErr, Color[TypeAnalysis], 33, titleYieldX, titleYRel,titleYRel+" of "+ titleYield[TypeAnalysis] + " yield vs multiplicity", 1, 0, 25);
  fHistYieldStatRelErr->GetYaxis()->SetTitleOffset(1.2);
  StyleHisto(fHistYieldSistRelErr, 10e-5, LimSupYieldErr, Color[TypeAnalysis], 27, titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, 25);
  StyleHisto(fHistYieldSistUpRelErr, 10e-5, LimSupYieldErr, Color[TypeAnalysis], 24, titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, 25);
  StyleHisto(fHistYieldSistLowRelErr, 10e-5, LimSupYieldErr, Color[TypeAnalysis], 25, titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, 25);
  StyleHisto(fHistYieldSistNoExtrRelErr, 10e-5, LimSupYieldErr, 881, 27,  titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, 25);
  legendYieldErr->AddEntry(fHistYieldStatRelErr, "stat.", "pl");
  //  if (!(TypeAnalysis==0 && type==8))   legendYieldErr->AddEntry(fHistYieldSistRelErr, "syst.", "pl");
  legendYieldErr->AddEntry(fHistYieldSistRelErr, "syst.", "pl");
  legendYieldErr->AddEntry(fHistYieldSistNoExtrRelErr, "syst. (no extr)", "pl");
  fHistYieldSistRelErr->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldStatRelErr->Draw("e p");
  fHistYieldSistNoExtrRelErr->Draw("same p");
  /*
  if (!(TypeAnalysis==0 && type==8))  fHistYieldSistRelErr->Draw("same p");
  else {
    fHistYieldSistUpRelErr->Draw("same p");
    fHistYieldSistLowRelErr->Draw("same p");
    legendYieldErr->AddEntry(fHistYieldSistUpRelErr, "syst. up", "pl");
    legendYieldErr->AddEntry(fHistYieldSistLowRelErr, "syst. low", "pl");
  }
  */
  fHistYieldSistRelErr->Draw("same p");
  legendYieldErr->Draw("");

  //fifth part: syst associated to bkg used for OOJ subtraction  (only for Xi and jet)

  cout << "\n\n going to write on file " << endl;   
  fileout->WriteTObject(canvasYield);
  fileout->WriteTObject(canvasYieldErr);
  fileout->WriteTObject(fHistYieldStat);
  fileout->WriteTObject(fHistYieldSist);
  fileout->WriteTObject(fHistYieldSistNoExtr);
  fileout->WriteTObject(fHistYieldStatRelErr);
  fileout->WriteTObject(fHistYieldSistRelErr);
  fileout->WriteTObject(fHistYieldSistNoExtrRelErr);
  if (TypeAnalysis==0 && type==8) {
    fileout->  WriteTObject(  fHistYieldSistUpRelErr);
    fileout->  WriteTObject(  fHistYieldSistLowRelErr);
    fileout->WriteTObject(gYield);
    fileout->WriteTObject(canvasOOJSubDef);
    fileout->WriteTObject(canvasBarlowOOJSubDef);
  }
  if (type==8){
  fileout->WriteTObject(canvasLeadTrackDef);
  fileout->WriteTObject(canvasBarlowLeadTrackDef);
  }
  if (type==0){
  fileout->WriteTObject(canvasMCChoiceDef);
  fileout->WriteTObject(canvasBarlowMCChoiceDef);
  }
  fileout->WriteTObject(canvasPtSpectraAll);
  fileout->WriteTObject(canvaspol0);
  fileout->WriteTObject(canvasPtSpectra);
  fileout->WriteTObject(canvasPtSpectraFit);
  fileout->WriteTObject(canvasPtSpectraFitBis);
  fileout->WriteTObject(canvasBarlow);
  fileout->WriteTObject(canvasBarlowpol0);
  fileout->WriteTObject(canvasPtSpectraRelError);
  fileout->WriteTObject(canvasPtSpectraRelErrorAll);
  for(Int_t m=0; m<nummolt+1; m++){
    fileout->WriteTObject(      fHistSpectrumSistAll[m]);
    fileout->WriteTObject(      fHistSpectrumStat[m]);
  }
  if (ZeroYieldLowPt && TypeAnalysis==0){
  canvasYield->SaveAs("PictureForNote/YieldvsMult"+tipo[type]+RegionType[TypeAnalysis]+"_ZeroYieldLowPt.pdf");
  canvasYieldErr->SaveAs("PictureForNote/YieldvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+"_ZeroYieldLowPt.pdf");
  }
  else {
  canvasYield->SaveAs("PictureForNote/YieldvsMult"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasYieldErr->SaveAs("PictureForNote/YieldvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  }
  canvasPtSpectra->SaveAs("PictureForNote/PtSpectra"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraFit->SaveAs("PictureForNote/PtSpectraFit"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraFit->SaveAs("PictureForNote/PtSpectraFitBis"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraRelError->SaveAs("PictureForNote/PtSpectraRelErr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraRelErrorAll->SaveAs("PictureForNote/PtSpectraRelErrAll"+tipo[type]+RegionType[TypeAnalysis]+".pdf");


  fileout->Close();

  if (TypeAnalysis==0)  cout <<" q " <<  pol0->GetParameter(0) << " red chisq " << pol0->GetChisquare()/pol0->GetNDF() << endl;
  cout << "starting from the file(s) "  << PathInDef<< endl;
  if (type==0 && TypeAnalysis==0)  cout << "starting from the file(s) "  << PathInpol0<< endl;
  cout << " I have created the file " << stringout << endl;
}

