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
#include <TRandom.h>
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

void PtSpectraBis(Int_t type=0,  Int_t TypeAnalysis=0,Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=0,   Int_t israp=0, TString year="1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"/"Run2DataRed_hXi"/*"2016kehjl_hK0s"*/, TString yearEtaEff = "1617_AOD234_hK0s",  TString Path1 ="_Jet0.75"/*"_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/,  Bool_t isEfficiency=0,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=0, Int_t PtBinning=0, TString FitFixed="Boltzmann"/*"Fermi-Dirac"*/,   Int_t sys=0, Bool_t ZeroYieldLowPt=0, Bool_t TwoFitFunctions=0, Bool_t isNormCorr=0, Bool_t isReanalysisWithEtaEff=0, Bool_t isReanalysisEtaEffComp=0, Bool_t isNormCorrFullyComputed=0, Int_t YieldComp=0, Bool_t isPreliminary=0){

  if (!isPreliminary){
    isReanalysisWithEtaEff=0;
    YieldComp=0;
    isReanalysisEtaEffComp=0;
    isNormCorrFullyComputed=1;
  }

  //isReanalysisWithEtaEff = 1 : spectra and yields are obtained when using eta-dependent efficiency
  //isReanalysisEtaEffComp = 1 : a comparison between spectra and yields obtained when using eta-dependent efficiency is done 
  if (isReanalysisEtaEffComp && isNormCorr) {cout << "when comparing with eta-eff corrected results, I should not multiply by normalisation factor " << endl; return;}
  if (isNormCorrFullyComputed==1) isNormCorr=0; //the two ways of normalization are different
  if (YieldComp==1){ //I compare the results corrected by the complete norm factor to the preliminary ones
    isNormCorr=1;
    isNormCorrFullyComputed=0;
    isReanalysisWithEtaEff=0;
  }

  if (YieldComp==2) {// I compare results corrected by complete norm factor and obtained using eta-dependence of efficiency to the preliminary ones
    isNormCorr=1;
    isNormCorrFullyComputed=0;
    isReanalysisWithEtaEff=0;
  }
  if (YieldComp>2) return; 

  if (TypeAnalysis>2) {cout << "sys errors not yet implemented for these regions " << endl; return;}

  if (isPreliminary){
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
  }
  else {
    if (type==0){
      PtBinning=1;
      if (!isMC) year= "1617_AOD234_hK0s";
      else  year= "";
    }
    else if (type==8){
      PtBinning=0;
      if (!isMC) year="";
      else  year=  "";
    }
  }

  TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
  TString RegionTypeOld[3] = {"Jet", "Bulk", "All"};
  TString RegionTypeNew[3] = {"Jet", "Bulk", "Full"};

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

  TString JetOrBulk[NumberTypeAnalysis]={"Jet", "Bulk", "All", "AllNS", "AwaySide", "AllButJet", "JetBFit", "JetZYAM", "ASBFit", "ASZYAM", "JetNS"};
  TString DataOrMC[2]={"Data", "MC"};
  Int_t jet=TypeAnalysis;

  const Int_t numtipo=10;
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
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
  if (!isPreliminary) PathIn0 +="_" + year;
  if (PtBinning>0) PathIn0 +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      PathIn0 +="_"+tipo[type];
    PathIn0 +=Srap[israp];
    PathIn0 +=SSkipAssoc[SkipAssoc];
  }

  stringout = Dir+"/DATA"+year0+"/";
  stringout += "PtSpectraBis" +hhCorr[ishhCorr];
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  //  if (year!="1617_hK0s" && year!="Run2DataRed_MECorr_hXi") stringout+= "_"+year;
  if (yearEtaEff=="1617_AOD234_hK0s" && isReanalysisWithEtaEff) stringout+= "_"+yearEtaEff;
  if (!isPreliminary) stringout+= "_"+year;
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
  if (isNormCorr) stringout+="_isNormCorr";
  if (isNormCorrFullyComputed) stringout+="_isNormCorrFullyComputed";
  if (isReanalysisWithEtaEff) stringout+="_isEtaEff";
  if (YieldComp ==1) stringout+="_isNormCorrFullyComputedComp";
  if (YieldComp ==2) stringout+="_isNormCorrFullyComputedAndEtaEffComp";
  //  stringout += "_Fixed";
  TString PathOutPictures = stringout;
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  cout << "\nDefinition integral regions... " << endl;
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
    if (PtBinning==1 && year=="1617_hK0s") {
      for (Int_t m =0; m<nummolt+1; m++){
	if (PtTrigMin==3){
	  if (!isMC){
	    //	    if (m==5|| m==2 || m==4 )     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	    //	    else  {LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;}
	    LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	    UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;	  
	    if (m==4) UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
	    //	    LowRangeBulk[m]= 0.1; 
	    LowRangeBulk[m]= 0.5; 
	    LowRangeAll[m]= 0.5; 
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
    
    if (PtBinning==1 && yearEtaEff=="1617_AOD234_hK0s" && isReanalysisWithEtaEff) {
      for (Int_t m=0; m<nummolt+1; m++){
	if (PtTrigMin==3){
	  LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	  UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=3;
	  LowRangeBulk[m]= 0.5; 
	  LowRangeAll[m]= 0.5; 
	  UpRangeAll[m]= 2; 
	  UpRangeBulk[m]= 2; 
	}
      }
    }
    else if (PtBinning==1 && yearEtaEff=="1617_AOD234_hK0s" && !isPreliminary) {
      for (Int_t m=0; m<nummolt+1; m++){
	LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=3;
	LowRangeBulk[m]= 0.5; 
	LowRangeAll[m]= 0.5; 
	UpRangeAll[m]= 2.5; 
	UpRangeBulk[m]= 2.5; 
      }
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
	if (m==0 || m==1 || m==4)	UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
	if (isNormCorrFullyComputed) LowRangeJet[m] = 1.5;
	//if (m==3) LowRangeBulk[m] = 1;
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
      //      else     LowRangeSpectrumPart[m] = 0.5;
    }
    else  {
      LowRangeSpectrumPart[m] = LowRange[m];
      if (type==8)       LowRangeSpectrumPart[m] = LowRange[m];
    }
    //    else  LowRangeSpectrumPart[m] = 8;
    for(Int_t v=0; v < numPtV0Max; v++){
      if (LowRangeSpectrumPart[m]<=NPtV0[v]) break;
      PtBinMin[m]++;
      cout << "LowRange " << LowRangeSpectrumPart[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
      cout << "PtBinMin " << PtBinMin[m] << endl;
    }

    UpRangeSpectrumPart[m] =8;
    if (type==8 && m==4 && (TypeAnalysis==0 || TypeAnalysis==1))    UpRangeSpectrumPart[m] =4;
    if (type==0 && m==4 && yearEtaEff=="1617_AOD234_hK0s" &&  isReanalysisWithEtaEff)    UpRangeSpectrumPart[m] =4;
  }

  cout << "...Done! " << endl;

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
  TH1F* fHistSpectrumStatUp[nummolt+1];
  TH1F* fHistSpectrumStatDown[nummolt+1];
  TH1F* fHistSpectrumTemp[nummolt+1];
  TH1F* fHistSpectrumStatpol0[nummolt+1];
  TH1F* fHistSpectrumStatOOJSubDef[nummolt+1];
  TH1F* fHistSpectrumStatLeadTrackDef[nummolt+1];
  TH1F* fHistSpectrumStatMCChoiceDef[nummolt+1];
  TH1F* fHistSpectrumStatFakeSBDef[nummolt+1];
  TH1F* fHistSpectrumStatEtaEff[nummolt+1];
  TH1F* fHistSpectrumStatEtaEffRatio[nummolt+1];
  TH1F* fHistSpectrumStatEtaEffRatioRef[nummolt+1];
  TH1F* fHistSpectrumStatNormCorr[nummolt+1];
  TH1F* fHistSpectrumStatNormCorrRatio[nummolt+1];
  TH1F* fHistSpectrumStatNormCorrRatioRef[nummolt+1];
  TH1F* fHistSpectrumSist[nummolt+1];
  TH1F* fHistSpectrum_max[nummolt+1];
  TH1F* fHistSpectrum_min[nummolt+1];

  TH1F* fHistSpectrumStatDPhi[nummolt+1][numsysDPhi];
  TH1F* fHistSpectrumSistDPhi[nummolt+1];
  TH1F* fHistSpectrumSistOOJ[nummolt+1];
  TH1F* fHistSpectrumSistAll[nummolt+1];
  TH1F* fHistSpectrumSistLeadTrack[nummolt+1];
  TH1F* fHistSpectrumSistMCChoice[nummolt+1];
  TH1F* fHistSpectrumSistFakeSB[nummolt+1];

  //histos for relative uncertainty
  TH1F* fHistSpectrumStatRelError[nummolt+1]; //stat
  TH1F* fHistSpectrumSistRelError[nummolt+1]; //sist on the DeltaPhi projections
  TH1F* fHistSpectrumSistRelErrorDPhi[nummolt+1]; //sist assoc to choice of DeltaPhi interval
  TH1F* fHistSpectrumSistRelErrorAll[nummolt+1]; //total sist
  TH1F* fHistSpectrumSistRelErrorDCAz[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorPurity[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorSE[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorDeltaEta[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMB[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorOOBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorIBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorLeadTrack[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMCChoice[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorFakeSB[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMCclosure[nummolt+1]; 

  TCanvas* canvasYield = new TCanvas ("canvasYield", "canvasYield", 1300, 800);
  TCanvas* canvasYieldEtaEff = new TCanvas ("canvasYieldEtaEff", "canvasYieldEtaEff", 1300, 800);
  TCanvas* canvasYieldEtaEffRatio = new TCanvas ("canvasYieldEtaEffRatio", "canvasYieldEtaEffRatio", 1300, 800);
  TCanvas* canvasYieldNormCorr = new TCanvas ("canvasYieldNormCorr", "canvasYieldNormCorr", 1300, 800);
  TCanvas* canvasYieldNormCorrRatio = new TCanvas ("canvasYieldNormCorrRatio", "canvasYieldNormCorrRatio", 1300, 800);
  TCanvas* canvasPtvsMult = new TCanvas ("canvasPtvsMult", "canvasPtvsMult", 1300, 800);
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
  TCanvas* canvasFakeSBDef = new TCanvas ("canvasFakeSBDef", "canvasFakeSBDef", 1300, 800);
  canvasFakeSBDef->Divide(3,2);

  TCanvas* canvasEtaEff = new TCanvas ("canvasEtaEff", "canvasEtaEff", 1300, 800);
  canvasEtaEff->Divide(3,2);
  TCanvas* canvasEtaEffRatio = new TCanvas ("canvasEtaEffRatio", "canvasEtaEffRatio", 1300, 800);
  canvasEtaEffRatio->Divide(3,2);

  TCanvas* canvasNormCorr = new TCanvas ("canvasNormCorr", "canvasNormCorr", 1300, 800);
  canvasNormCorr->Divide(3,2);
  TCanvas* canvasNormCorrRatio = new TCanvas ("canvasNormCorrRatio", "canvasNormCorrRatio", 1300, 800);
  canvasNormCorrRatio->Divide(3,2);

  TCanvas* canvasPtSpectraFit = new TCanvas ("canvasPtSpectraFit", "canvasPtSpectraFit", 1300, 800);
  canvasPtSpectraFit->Divide(3,2);
  TCanvas* canvasPtSpectraFitUp = new TCanvas ("canvasPtSpectraFitUp", "canvasPtSpectraFitUp", 1300, 800);
  canvasPtSpectraFitUp->Divide(3,2);
  TCanvas* canvasPtSpectraFitDown = new TCanvas ("canvasPtSpectraFitDown", "canvasPtSpectraFitDown", 1300, 800);
  canvasPtSpectraFitDown->Divide(3,2);

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
  TCanvas* canvasBarlowFakeSBDef = new TCanvas ("canvasBarlowFakeSBDef", "canvasBarlowFakeSBDef", 1300, 800);
  canvasBarlowFakeSBDef->Divide(3,2);


  TCanvas* canvasPtSpectraRelError = new TCanvas ("canvasPtSpectraRelError", "canvasPtSpectraRelError", 1300, 800);
  canvasPtSpectraRelError->Divide(3,2);
  TCanvas* canvasPtSpectraRelErrorAll = new TCanvas ("canvasPtSpectraRelErrorAll", "canvasPtSpectraRelErrorAll", 1300, 800);
  canvasPtSpectraRelErrorAll->Divide(3,2);

  // cout << "\nprendo histo per confronto con dati pubblicati " << endl;                                        
  TString PathDatiPubblicati ="";                                                                                 
  if (type==0) PathDatiPubblicati = "HEPData-ins1748157-v1-Table_1.root";                                         
  else if (type==8) PathDatiPubblicati = "HEPData-1583750454-v1-Table_3.root";                                    
  TFile *filedatipubblicati = new TFile(PathDatiPubblicati, "");                                                  
  if (!filedatipubblicati) {cout << "file dati pubblicati not there " << endl; return;}                           
  TString STable = "Table 3";                                                                                     
  if (type==0) STable = "Table 1";                                                                                
  TDirectoryFile *dirspectra = (TDirectoryFile*)filedatipubblicati->Get(STable);                                  
  if (!dirspectra)  {cout << "directory dati pubblicati not there " << endl; return;}              

  TH1F *   hspectrum[11];
  TH1F *   hspectrum1[11];
  TH1F *   hspectrum2[11];
  TH1F *   hspectrum3[11];

  for (Int_t i=0; i<11; i++){
    hspectrum[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i", i+1));                                            
    hspectrum1[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e1", i+1)); //stat
    hspectrum2[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e2", i+1)); //syst total                                
    hspectrum3[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e3", i+1)); //sist Uncorr                           

    if (!hspectrum[i] ||     !hspectrum1[i] || !hspectrum2[i]|| !hspectrum3[i] ) { cout << "histo is missing " << endl; return;}                                                                                              
  }


  TH1F*  fHistSpectrumSistRelErrorPublished = (TH1F*) hspectrum2[10]->Clone("fHistSpectrumSistRelErrorPublished");
  fHistSpectrumSistRelErrorPublished->Divide(hspectrum[10]);
  for (Int_t b=0; b<= fHistSpectrumSistRelErrorPublished->GetNbinsX(); b++){
    fHistSpectrumSistRelErrorPublished->SetBinError(b,0);
  }

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
  Bool_t BarlowPassedFakeSBDef[nummolt+1]={0};
  Int_t BarlowSignFakeSBDef[nummolt+1]={0};
  TH1F* hBarlowVarFakeSBDef[nummolt+1];
  Float_t BarlowVarFakeSBDef[nummolt+1][numPtV0]={0};

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
    if (TypeAnalysis==0)    LimInfYield=0.015;
    else    LimInfYield=10e-5;
  }

  Float_t LimInfPtvsMult=10e-5;
  Float_t LimSupPtvsMult=3;
  if (type==0){
    if (TypeAnalysis==0) LimSupPtvsMult=4; //2.3
    if (TypeAnalysis==1) LimSupPtvsMult=2; //1.2
    if (TypeAnalysis==2) LimSupPtvsMult=2; //1.2
  }
  if (type==8){
    if (TypeAnalysis==0) LimSupPtvsMult=8; //5
    if (TypeAnalysis==1) LimSupPtvsMult=4; //2.4
    if (TypeAnalysis==2) LimSupPtvsMult=4;
  }
  if (type==0){
    if (TypeAnalysis==0) LimInfPtvsMult=1.7;
    if (TypeAnalysis==1) LimInfPtvsMult=0.8;
    if (TypeAnalysis==2) LimInfPtvsMult=0.8;
  }
  if (type==8){
    if (TypeAnalysis==0) LimInfPtvsMult=2;
    if (TypeAnalysis==1) LimInfPtvsMult=1.2;
    if (TypeAnalysis==2) LimInfPtvsMult=1.2;
  }

  Float_t LimSupYieldErr=0.04;
  if (type==0){
    if (TypeAnalysis==1) LimSupYieldErr=0.02;
    if (TypeAnalysis==2) LimSupYieldErr=0.02;
    LimSupYieldErr=0.2;
  }
  if (type==8){
    if (TypeAnalysis==0) LimSupYieldErr=0.45;
    if (TypeAnalysis==1) LimSupYieldErr=0.1;
    if (TypeAnalysis==2) LimSupYieldErr=0.1;
    LimSupYieldErr=0.4;
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
    //    LimInf =10e-10;
    LimInf =10e-6;
  }

  Float_t LimSupNormComp=1.2;
  Float_t LimInfNormComp=0.8;
  if (type==8) {
    LimSupNormComp=1.4;
    LimInfNormComp=0.6;
  }

  Float_t LimSupError=0.01;
  Float_t LimInfError=0;
  Float_t LimSupErrorLog=0.01;
  if (type==0) {
    if (TypeAnalysis==0){ LimSupError =0.6; LimSupErrorLog =0.5;}
    else   if (TypeAnalysis==1) LimSupError =0.15;
    else   if (TypeAnalysis==2) LimSupError =0.05;
    LimSupError=0.3;
    LimInfError=10e-5;
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSupError =1;
    else   if (TypeAnalysis==1) LimSupError =0.3;
    else   if (TypeAnalysis==2) LimSupError =0.1;
    LimSupError=0.5;
    LimInfError=10e-5;
  }

  TH1F*  hDeltaPhiLimit;
  Float_t  LowBinDPhi[numsysDPhi] ={0};
  Float_t   UpBinDPhi[numsysDPhi]={0}  ;

  //first part: I get default spectra with stat uncertainty
  cout << "\n**************************"<<endl;
  cout << "First part: I get default spectra with stat uncertainty";
  TString  PathInDef=PathIn0;
  PathInDef+=   Form("_SysPhi%i_PtMin%.1f_", 0, PtTrigMin);
  PathInDef+= RegionType[TypeAnalysis];
  PathInDef += ".root";
  cout << " from the file: " << PathInDef << endl;
  TFile *  fileinDef = new TFile(PathInDef, "");
  if (!fileinDef) {cout << PathInDef << " does not exist" << endl; return;}
  for(Int_t m=0; m<nummolt+1; m++){
    //    cout << " m " << m << endl;
    fHistSpectrumSistRelErrorDCAz[m]    =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorDCAz_"+Smolt[m]); 
    fHistSpectrumSistRelErrorPurity[m]    =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorPurity_"+Smolt[m]); 
    fHistSpectrumSistRelErrorSE[m]      =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorSE_"+Smolt[m]);
    fHistSpectrumSistRelErrorDeltaEta[m]=(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorDeltaEta_"+Smolt[m]);
    if (!fHistSpectrumSistRelErrorDCAz[m]    ) return;
    if (!fHistSpectrumSistRelErrorPurity[m]    ) return;
    if (!fHistSpectrumSistRelErrorSE[m]      ) return;
    if (!fHistSpectrumSistRelErrorDeltaEta[m]) return;

    fHistSpectrumStat[m]=(TH1F*)fileinDef->Get("fHistSpectrum_"+Smolt[m]);
    fHistSpectrumStatUp[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumUp_"+Smolt[m]);
    fHistSpectrumStatDown[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumDown_"+Smolt[m]);
    fHistSpectrumTemp[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumTemp_"+Smolt[m]);
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
    fHistSpectrumSistRelErrorFakeSB[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorFakeSB_"+Smolt[m]);


    fHistSpectrum_max[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMax_"+Smolt[m]);
    fHistSpectrum_min[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMin_"+Smolt[m]);
    fHistSpectrumSistDPhi[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDPhi_"+Smolt[m]);
    fHistSpectrumSistOOJ[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistOOJ_"+Smolt[m]);
    fHistSpectrumSistAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistAll_"+Smolt[m]);
    fHistSpectrumSistLeadTrack[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistLeadTrack_"+Smolt[m]);
    fHistSpectrumSistMCChoice[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistMCChoice_"+Smolt[m]);
    fHistSpectrumSistFakeSB[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistFakeSB_"+Smolt[m]);

    for(Int_t v=PtV0Min; v <    fHistSpectrumSistAll[m]->GetNbinsX() ; v++){
      //      cout << "v " << v << endl;
      if (fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSist[m]->GetBinError(v+1)/ fHistSpectrumSist[m]->GetBinContent(v+1)); //this works for inclusive, for jet and OOJ I will change it
	fHistSpectrumStatRelError[m]->SetBinContent(v+1, fHistSpectrumStat[m]->GetBinError(v+1)/ fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      else       if (type==0 && NPtV0[v] == 0.1){
	fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
	//      cout << "sist rel error " <<    fHistSpectrumSistRelErrorAll[m]->GetBinContent(v+1) << endl;
	fHistSpectrumStatRelError[m]->SetBinError(v+1,0);
      }
      else{
	fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
	fHistSpectrumStatRelError[m]->SetBinError(v+1,0);
      }
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
    cout << "Second part: I evaluate syst uncertainty associated to choice of DeltaPhi (for jet and out of jet)\n " << endl;

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
	if (sysDPhi==1) fHistSpectrumStat[m]->DrawClone("");
	fHistSpectrumStatDPhi[m][sysDPhi]->Draw("same");
	if (sysDPhi==numsysDPhi-1) legendPhi->Draw("same");

	hBarlowVar[m][sysDPhi]=(TH1F*)    fHistSpectrumStatDPhi[m][sysDPhi]->Clone("fHistBarlowVar_"+Smolt[m]);
	for ( Int_t b=1; b<=hBarlowVar[m][sysDPhi]->GetNbinsX(); b++){
	  hBarlowVar[m][sysDPhi]->SetBinContent(b,0);
	  hBarlowVar[m][sysDPhi]->SetBinError(b,0);
	}
	for(Int_t v=PtBinMin[m]; v < numPtV0Max; v++){
	  //if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
	  if (v==1 && TypeAnalysis==0 && type==8) continue;
	  if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
	  if (NPtV0[v] == 4 && type==0 && m==4) continue;
	  //	cout << "\nv: " << v << endl;
	  //	cout << "...histo...dphi " << sysDPhi << " " << fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1)<< endl;
	  BarlowVar[m][v][sysDPhi] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) - pow(fHistSpectrumStatDPhi[m][sysDPhi]->GetBinError(v+1),2)));
	  if (TMath::Abs(BarlowVar[m][v][sysDPhi])>2) BarlowSign[m][sysDPhi]++;
	  hBarlowVar[m][sysDPhi] ->SetBinContent(v+1, BarlowVar[m][v][sysDPhi]) ;
	  hBarlowVar[m][sysDPhi] ->SetBinError(v+1, 0) ;
	  cout << "barlow var " << NPtV0[v] << " "  << 	  hBarlowVar[m][sysDPhi] ->GetBinContent(v+1)<<" " <<	  hBarlowVar[m][sysDPhi] ->GetBinError(v+1) << endl;

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
	if (isPreliminary){
	  if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
	  if (v==1 && TypeAnalysis==0 && type==8) continue;
	  if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
	}
	else {
	  if (v==1 && TypeAnalysis==0 && type==0) continue;
	  if (NPtV0[v]==4 && type==0 && m==4) continue;
	}
	//      cout << " v: " << v << endl;
	ErrDPhi[m][v]= TMath::Abs(	  fHistSpectrum_min[m]->GetBinContent(v+1) - 	  fHistSpectrum_max[m]->GetBinContent(v+1))/2;
	//	cout << fHistSpectrum_min[m]->GetBinContent(v+1)<< " << " << fHistSpectrumStat[m]->GetBinContent(v+1)<< " << " <<fHistSpectrum_max[m]->GetBinContent(v+1)<<endl;
	//	cout << "err dphi " << ErrDPhi[m][v] << endl;
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
  //cout << "\n*******************************************"<<endl;
  //cout<< "hXi: systematic effect associated to selection pt,assoc < pt,Trigg  " << endl;
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
	//cout << "barlow var v" << v << " "  << 	  hBarlowVarLeadTrackDef[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarLeadTrackDef[m] ->GetBinError(v+1) << endl;

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
      fHistSpectrumStat[m]->DrawClone("ep");
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

  cout << "\n*******************************************"<<endl;
  cout<< "Systematic effect associated to MC used to calculate efficiency  " << endl;
  Float_t   YieldSpectrumErrMCChoice[nummolt+1]={0};
  TString    PathInMCChoiceDef;
  TLegend * legendMCChoice= new TLegend(0.6, 0.6, 0.9, 0.9);
  Int_t IsSignMCChoice=0;
  Int_t Varm = 0;
  TString PathMCChoiceSyst="";

  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistMCChoice[m]->SetBinError(v+1,0);
    }
  }

  if (isPreliminary){
    if (type==0 || type==8){
      // cout << "\n*******************************************"<<endl;
      // cout<< "Systematic effect associated to choice of MC used to calculat K0s/Xi efficiency  " << endl;
      PathInMCChoiceDef = Dir+"/DATA"+year0+"/SystematicAnalysis" + year;
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
	gPad->SetLogy();
	gPad->SetLeftMargin(0.15);
	StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
	if (m==0)      legendMCChoice->AddEntry(fHistSpectrumStat[m], "PYTHIA8", "ple");
	fHistSpectrumStat[m]->DrawClone("ep");
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
	    IsSignMCChoice++;
	    Varm = m;
	    fHistSpectrumSistMCChoice[m]->SetBinError(v+1,	  (fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatMCChoiceDef[m]->GetBinContent(v+1))/2);
	    //YieldSpectrumErrMCChoice[m]+= TMath::Abs(   fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatMCChoiceDef[m]->GetBinContent(v+1));
	  }
	}
      }//end loop m
      for(Int_t m=0; m<nummolt+1; m++){
	for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	  if (IsSignMCChoice>0){
	    fHistSpectrumSistMCChoice[m]->SetBinError(v+1,	  (fHistSpectrumStat[Varm]->GetBinContent(v+1) -    fHistSpectrumStatMCChoiceDef[Varm]->GetBinContent(v+1))/2/fHistSpectrumSistMCChoice[Varm]->GetBinContent(v+1)*fHistSpectrumSistMCChoice[m]->GetBinContent(v+1));
	  }
	}
      }
    }//end loop syst associated to MC used to calculate efficiency
  }
  else {
    PathMCChoiceSyst = "FinalOutput/DATA2016/PtSpectraBis";
    if (PtBinning>0) PathMCChoiceSyst += "_PtBinning1";
    PathMCChoiceSyst += "_" + tipo[type];
    PathMCChoiceSyst +=  Srap[israp];
    PathMCChoiceSyst += "_" + SSkipAssoc[SkipAssoc];
    PathMCChoiceSyst +=Form("PtMin%.1f_",PtTrigMin ) + RegionType[TypeAnalysis]+ ".root";
    TFile * fileMCChoiceSyst = new TFile (PathMCChoiceSyst, "");
    if (!fileMCChoiceSyst) return;
    for(Int_t m=0; m<nummolt+1; m++){
      fHistSpectrumSistRelErrorMCChoice[m]= (TH1F*)  fileMCChoiceSyst->Get("fHistSpectrumSistRelErrorMCChoice_"+Smolt[m]);
      if (!fHistSpectrumSistRelErrorMCChoice[m]) {cout << "Spectrum with syst. uncertainty related to MC choice not there " << endl; return; }
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistMCChoice[m]->SetBinError(v+1, fHistSpectrumSistRelErrorMCChoice[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinContent(v+1));
	if ( fHistSpectrumStat[m] ->GetBinContent(v+1) ==0) fHistSpectrumSistRelErrorMCChoice[m]->SetBinContent(v+1, 0);
      }
    }
  }

  //diagnostic cout:
  for(Int_t m=0; m<nummolt+1; m++){
    if (m!=nummolt) continue; //choose the multiplicity for diagnosis
    cout << "\nSyst. uncert. associated to spectrum (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " +- " << fHistSpectrumSistMCChoice[m]->GetBinError(b)<< " (" << fHistSpectrumSistMCChoice[m]->GetBinError(b) / fHistSpectrumStat[m]->GetBinContent(b)  << ") " << endl;
    }
  }

  //    cout << "\n*******************************************"<<endl;
  //    cout<< "systematic effect associated to how fake K0s/Xi are sub  " << endl;
  Float_t   YieldSpectrumErrFakeSB[nummolt+1]={0};
  TString    PathInFakeSBDef;
  TString PathFakeSBSyst="";
  TLegend * legendFakeSB= new TLegend(0.6, 0.6, 0.9, 0.9);
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistFakeSB[m]->SetBinError(v+1,0);
    }
  }
  if (type==0){
    cout << "\n*******************************************"<<endl;
    cout<< "Systematic effect associated to how fake K0s/Xi are sub  " << endl;
    PathInFakeSBDef = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    if (PtBinning>0) PathInFakeSBDef +=Form("_PtBinning%i",PtBinning);
    if(type>=0){
      PathInFakeSBDef +="_"+tipo[type];
      PathInFakeSBDef +=Srap[israp];
      PathInFakeSBDef +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInFakeSBDef+="_HighPtExtr";
    PathInFakeSBDef+= hhCorr[ishhCorr]+"_"+RegionTypeOld[TypeAnalysis] +DataOrMC[isMC] + Form("_PtMin%.1f", PtTrigMin);
    PathInFakeSBDef += "_Sidebands";
    if (!isPreliminary) PathInFakeSBDef+="_IsEtaEff";
    PathInFakeSBDef += ".root";
    //    cout << "\nFrom file: " << PathInFakeSBDef << endl;
    TFile *  fileinFakeSBDef = new TFile(PathInFakeSBDef, "");
    if (!fileinFakeSBDef) {cout << PathInFakeSBDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummolt+1; m++){
      //    cout << " m " << m << endl;
      fHistSpectrumStatFakeSBDef[m]    =(TH1F*)fileinFakeSBDef->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      if (!fHistSpectrumStatFakeSBDef[m]) {cout << " I was looking for histo in " << PathInFakeSBDef  << endl; return;}

      hBarlowVarFakeSBDef[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistBarlowVarFakeSBDef_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	//	if (NPtV0[v] <3) continue;
	if (fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1) ==0) continue;
	BarlowVarFakeSBDef[m][v] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) -pow(fHistSpectrumStatFakeSBDef[m]->GetBinError(v+1),2)));
	if (TMath::Abs(BarlowVarFakeSBDef[m][v])>2) BarlowSignFakeSBDef[m]++;
	hBarlowVarFakeSBDef[m] ->SetBinContent(v+1, BarlowVarFakeSBDef[m][v]) ;
	hBarlowVarFakeSBDef[m] ->SetBinError(v+1, 0) ;
	//	cout << "barlow var v" << v << " "  << 	  hBarlowVarFakeSBDef[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarFakeSBDef[m] ->GetBinError(v+1) << endl;

      }//end loop v

      if (BarlowSignFakeSBDef[m]>= 3) BarlowPassedFakeSBDef[m]=1;
      StyleHisto(hBarlowVarFakeSBDef[m], -5, 5, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);
      if (BarlowPassedFakeSBDef[m]) StyleHisto(hBarlowVarFakeSBDef[m], -5, 5, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);

      //      cout << "Barlow passed " << BarlowPassedFakeSBDef[m] << endl;
      canvasBarlowFakeSBDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      hBarlowVarFakeSBDef[m] ->Draw("same p");
      lineat2->Draw("same");
      lineatm2->Draw("same");

      canvasFakeSBDef->cd(m+1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendFakeSB->AddEntry(fHistSpectrumStat[m], "Default", "ple");
      fHistSpectrumStat[m]->DrawClone("ep");
      StyleHisto(fHistSpectrumStatFakeSBDef[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendFakeSB->AddEntry(fHistSpectrumStatFakeSBDef[m], "Sidebands", "ple");
      fHistSpectrumStatFakeSBDef[m]->Draw("same ep");
      legendFakeSB->Draw("");

      //      cout << " going to set errors " << endl;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistFakeSB[m]->SetBinError(v+1,0);
	if (fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1) ==0) continue;
	//	if(TMath::Abs(BarlowVarFakeSBDef[m][v])>2){
	//	if (BarlowPassedFakeSBDef[m]){
	if (m==nummolt)	  fHistSpectrumSistFakeSB[m]->SetBinError(v+1,	  (fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1))/2);
	//	}
      }
    }//end loop m
    for(Int_t m=0; m<nummolt; m++){
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistFakeSB[m]->SetBinError(v+1,	  (fHistSpectrumStat[5]->GetBinContent(v+1) -    fHistSpectrumStatFakeSBDef[5]->GetBinContent(v+1))/2/ fHistSpectrumSistFakeSB[5]->GetBinContent(v+1)* fHistSpectrumSistFakeSB[m]->GetBinContent(v+1));
      }
    }
  }//end loop syst associated to how fake K0s/Xi are subtracted

  if (type==8){
    cout << " getting rel uncert from K0s" << endl;
    TString SfileinSB = "FinalOutput/DATA2016/PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_";
    SfileinSB+= RegionType[TypeAnalysis];
    SfileinSB+="_isNormCorr.root";
    cout << SfileinSB << endl;
    TFile * fileinSB = new TFile(SfileinSB, "");
    TH1F *fHistSpectrumSistFakeSBK0s =  (TH1F*) fileinSB->Get("fHistSpectrumSistFakeSB_"+Smolt[5]);
    if (!fHistSpectrumSistFakeSBK0s) return;
    Float_t     ErrorSistRelFakeSBK0s[numPtV0]={0};
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      ErrorSistRelFakeSBK0s[v]=  fHistSpectrumSistFakeSBK0s->GetBinError( fHistSpectrumSistFakeSBK0s->FindBin(NPtV0[v]+0.0001)) /fHistSpectrumSistFakeSBK0s->GetBinContent( fHistSpectrumSistFakeSBK0s->FindBin(NPtV0[v]+0.0001));
      for(Int_t m=0; m<nummolt+1; m++){
	fHistSpectrumSistFakeSB[m]->SetBinError(v+1,      fHistSpectrumSistFakeSB[m]->GetBinContent(v+1)* ErrorSistRelFakeSBK0s[v]);
      }
    }
  }

  if (!isPreliminary) {
    PathFakeSBSyst = "FinalOutput/DATA2016/PtSpectraBis";
    if (PtBinning>0) PathFakeSBSyst += "_PtBinning1";
    PathFakeSBSyst += "_" + tipo[type];
    PathFakeSBSyst +=  Srap[israp];
    PathFakeSBSyst += "_" + SSkipAssoc[SkipAssoc];
    PathFakeSBSyst +=Form("PtMin%.1f_",PtTrigMin ) + RegionType[TypeAnalysis]+ ".root";
    cout << "\nFrom file: " << PathFakeSBSyst << endl;
    TFile * fileFakeSBSyst = new TFile (PathFakeSBSyst, "");
    if (!fileFakeSBSyst) return;
    for(Int_t m=0; m<nummolt+1; m++){
      fHistSpectrumSistRelErrorFakeSB[m]= (TH1F*)  fileFakeSBSyst->Get("fHistSpectrumSistRelErrorFakeSB_"+Smolt[m]);
      if (!fHistSpectrumSistRelErrorFakeSB[m]) {cout << "Spectrum with syst. uncertainty related to MC choice not there " << endl; return; }
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistFakeSB[m]->SetBinError(v+1, fHistSpectrumSistRelErrorFakeSB[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinContent(v+1));
	if ( fHistSpectrumStat[m] ->GetBinContent(v+1) ==0) fHistSpectrumSistRelErrorFakeSB[m]->SetBinContent(v+1, 0);
      }
    }
  }
  
  //diagnostic cout:
  for(Int_t m=0; m<nummolt+1; m++){
    if (m!=nummolt) continue; //choose the multiplicity for diagnosis
    cout << "\nSyst. uncert. associated to spectrum (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " +- " << fHistSpectrumSistFakeSB[m]->GetBinError(b) << " (" << fHistSpectrumSistFakeSB[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ") " << endl;
    }
  }

  TH1F* fHistRelErrorClosure[nummolt+1];
  TH1F* fHistRelErrorClosureK0s[nummolt+1];
  TString SfileMCClosure;
  if (isPreliminary){
    cout << "\n\n error on closure ****************"<< endl;
    //error associated to closure
    if (type==0) SfileMCClosure= "FinalOutput/DATA2016/MCClosureCheck_RecoToHybrid2018f1_extra_Reco_hK0s_vs_2018f1_extra_Hybrid_hK0s_K0s_Eta0.8_PtMin3.0_IsParticleTrue_isEffMassSel.root";
    else  SfileMCClosure= "FinalOutput/DATA2016/MCClosureCheck_RecoToHybridAllMC_hXi_vs_2018f1g4_extra_hXi_Hybrid_Xi_Eta0.8_PtMin3.0_IsParticleTrue.root";
    TFile* fileMCClosure= new TFile(SfileMCClosure, "");
    if (!fileMCClosure) {cout << " no MC closure file " << endl; return;}
    for(Int_t m=0; m<nummolt+1; m++){
      if (type==8) {
	fHistRelErrorClosure[m] = (TH1F*) fileMCClosure->Get("SpectrumRelErrorClosure"+RegionTypeOld[TypeAnalysis]+Form("_m%i",5));
	fHistRelErrorClosure[m]->SetName("fSpectrumRelErrorClosure"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));
      }
      if (type==0){
	fHistRelErrorClosure[m] = (TH1F*) fHistSpectrumStat[m]->Clone("fSpectrumRelErrorClosure"+RegionTypeOld[TypeAnalysis]+Form("_m%i",5));
	fHistRelErrorClosure[m] ->SetName("fSpectrumRelErrorClosure"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));
	fHistRelErrorClosureK0s[m] = (TH1F*) fileMCClosure->Get("SpectrumRelErrorClosure"+RegionTypeOld[TypeAnalysis]+Form("_m%i",5));
	fHistRelErrorClosureK0s[m] ->SetName("SpectrumRelErrorClosure"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));
	for (Int_t b=1; b<= fHistRelErrorClosure[m]->GetNbinsX(); b++){
	  //	cout << "before " << fHistRelErrorClosure[m]->GetBinContent(b) << endl;
	  //cout << "I should put " << 	fHistRelErrorClosureK0s[m]->GetBinContent( fHistRelErrorClosureK0s[m]->FindBin( fHistRelErrorClosure[m]->GetBinCenter(b))) << endl; 
	  fHistRelErrorClosure[m]->SetBinContent(b, fHistRelErrorClosureK0s[m]->GetBinContent( fHistRelErrorClosureK0s[m]->FindBin( fHistRelErrorClosure[m]->GetBinCenter(b))));
	  fHistRelErrorClosure[m]->SetBinError(b,0);
	  //cout << "after " <<fHistRelErrorClosure[m]->GetBinContent(b) << endl;
	}
      }
      if (!fHistRelErrorClosure[m]) return;
      fHistSpectrumSistRelErrorMCclosure[m]=(TH1F*)fHistRelErrorClosure[m]->Clone("fHistSpectrumSistRelErrorMCclosure_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  fHistSpectrumSistRelErrorMCclosure[m]->SetBinContent(v+1, TMath::Abs(fHistSpectrumSistRelErrorMCclosure[m]->GetBinContent(v+1))/2);
	}
	else 	fHistSpectrumSistRelErrorMCclosure[m]->SetBinContent(v+1,0);
      }
      fHistSpectrumSistRelErrorMCclosure[m]->Smooth();
      fHistSpectrumSistRelErrorMCclosure[m]->SetBinContent(1,0);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	//new entry! 
	if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  fHistSpectrumSistRelErrorMCclosure[m]->SetBinContent(v+1,0.05);
	}
	if (TypeAnalysis!=2 && type==8 && m==4 && NPtV0[v]==4){
	  fHistSpectrumSistRelErrorMCclosure[m]->SetBinContent(v+1,0);
	  cout <<	fHistSpectrumSistRelErrorMCclosure[m]->GetBinContent(v+1)<<endl;
	}
      }
    }
  }
  else {
    for(Int_t m=0; m<nummolt+1; m++){
      fHistSpectrumSistRelErrorMCclosure[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMCclosure_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistRelErrorMCclosure[m]->SetBinContent(v+1, 0);
	fHistSpectrumSistRelErrorMCclosure[m]->SetBinError(v+1, 0);
      }
    }
  }

  if (isReanalysisWithEtaEff){
    for(Int_t m=0; m<nummolt+1; m++){
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistRelErrorMCclosure[m]->SetBinContent(v+1,0);
      }
    }
  }

  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      Sigma2OOBPU[m] =pow(OOBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
      Sigma2IBPU[m] =pow(IBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
      Sigma2MB[m]= pow(MB *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
      fHistSpectrumSistAll[m]->SetBinError(v+1,sqrt(pow(fHistSpectrumSistDPhi[m]->GetBinError(v+1),2) + pow(fHistSpectrumSist[m]->GetBinError(v+1),2) + Sigma2OOBPU[m]+ Sigma2IBPU[m] + Sigma2MB[m] + pow(fHistSpectrumSistLeadTrack[m]->GetBinError(v+1),2)+  pow(fHistSpectrumSistMCChoice[m]->GetBinError(v+1),2)+ pow(fHistSpectrumSistFakeSB[m]->GetBinError(v+1),2) +  pow(TMath::Abs(fHistSpectrumSistRelErrorMCclosure[m]->GetBinContent(v+1))*fHistSpectrumStat[m]->GetBinContent(v+1),2)));
      //      fHistSpectrumSistAll[m]->SetBinError(v+1, fHistSpectrumSistDPhi[m]->GetBinError(v+1) + fHistSpectrumSist[m]->GetBinError(v+1) +sqrt( Sigma2OOBPU[m])+sqrt( Sigma2IBPU[m]) +sqrt( Sigma2MB[m]) + fHistSpectrumSistLeadTrack[m]->GetBinError(v+1)+  fHistSpectrumSistMCChoice[m]->GetBinError(v+1)+ fHistSpectrumSistFakeSB[m]->GetBinError(v+1));

      if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSistAll[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumSistRelErrorMB[m]->SetBinContent(v+1,MB);
	fHistSpectrumSistRelErrorIBPU[m]->SetBinContent(v+1,IBPU);
	fHistSpectrumSistRelErrorOOBPU[m]->SetBinContent(v+1,OOBPU);
	fHistSpectrumSistRelErrorLeadTrack[m]->SetBinContent(v+1,fHistSpectrumSistLeadTrack[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
       	fHistSpectrumSistRelErrorMCChoice[m]->SetBinContent(v+1,fHistSpectrumSistMCChoice[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumSistRelErrorFakeSB[m]->SetBinContent(v+1,fHistSpectrumSistFakeSB[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	
      }

      fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMB[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorIBPU[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorOOBPU[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorLeadTrack[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMCChoice[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorFakeSB[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMCclosure[m]->SetBinError(v+1,0);
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
    StyleHisto(fHistSpectrumSistRelErrorPurity[m], LimInfError, LimSupError, 834, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorDeltaEta[m], LimInfError, LimSupError, 881, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorDPhi[m], LimInfError, LimSupError, 922, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorMB[m], LimInfError, LimSupError, 401, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorOOBPU[m], LimInfError, LimSupError, 825, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorIBPU[m], LimInfError, LimSupError, 631, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorLeadTrack[m], LimInfError, LimSupError, 630, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorMCChoice[m], LimInfError, LimSupError, 907, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorFakeSB[m], LimInfError, LimSupError, 600, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    StyleHisto(fHistSpectrumSistRelErrorMCclosure[m], LimInfError, LimSupError, 855, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
    gPad->SetLogy();

    TLegend *legendErrorAll= new TLegend(0.6, 0.1, 0.9, 0.4);
    legendErrorAll->AddEntry(fHistSpectrumStatRelError[m], "stat.", "pl");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorSE[m], "syst. signal extr.", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDCAz[m], "syst. DCAzTrigger", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorPurity[m], "syst. purity", "l");
    if(TypeAnalysis!=2)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDeltaEta[m], "syst. #Delta#eta region", "l");
    if(TypeAnalysis!=2)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDPhi[m], "syst. #Delta#phi region", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMB[m], "Material Budget", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorOOBPU[m], "Out-of-bunch PU", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorIBPU[m], "In-bunch PU", "l");
    if (isPreliminary)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMCclosure[m], "MC closure", "l");
    if(BarlowPassedLeadTrackDef[m])   legendErrorAll->AddEntry(fHistSpectrumSistRelErrorLeadTrack[m], "p_{T} < p_{T}^{Trigg}", "l");
    //    if(BarlowPassedMCChoiceDef[m])  
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMCChoice[m], "MCChoice", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorFakeSB[m], "Fake "+tipo[type], "l");

    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumSistRelErrorMB[m]->Draw("same");
    fHistSpectrumSistRelErrorIBPU[m]->Draw("same");
    fHistSpectrumSistRelErrorOOBPU[m]->Draw("same");
    fHistSpectrumSistRelErrorSE[m]->Draw("same");
    fHistSpectrumSistRelErrorDCAz[m]->Draw("same");
    fHistSpectrumSistRelErrorPurity[m]->Draw("same");
    fHistSpectrumSistRelErrorLeadTrack[m]->Draw("same");
    fHistSpectrumSistRelErrorMCChoice[m]->Draw("same");
    fHistSpectrumSistRelErrorFakeSB[m]->Draw("same");
    fHistSpectrumSistRelErrorMCclosure[m]->Draw("same");
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
    cout << "3. Syst. associated to pol0 OOJ subtraction for hK0s"<< endl;

    PathInpol0 = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    if (!isPreliminary && PtBinning>0) PathInpol0 += "_PtBinning1";
    if (isPreliminary)    PathInpol0 +=Path1; //it was without year                                                        
    if(type>=0){
      PathInpol0 +="_"+tipo[type];
      PathInpol0 +=Srap[israp];
      PathInpol0 +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInpol0+="_HighPtExtr";
    PathInpol0+= hhCorr[ishhCorr]+"_JetBFit" +DataOrMC[isMC] + Form("_PtMin%.1f", PtTrigMin);
    if (isPreliminary) PathInpol0 += "_Try1";
    if (!isPreliminary) PathInpol0 += "_IsEtaEff";
    PathInpol0 +=".root";
    cout << "\nFrom file: " << PathInpol0 << endl;
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
	if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
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
      fHistSpectrumStat[m]->DrawClone("ep");
      StyleHisto(fHistSpectrumStatpol0[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0,0,0);
      fHistSpectrumStatpol0[m]->Draw("same ep");

    }
  }


  //******ReAnalysis with eta-dependent efficiency**********++

  //I plot spectra with stat + total syst uncertainty
  // I plot stat + syst (total, DeltaPhi assoc., OOJ sub assoc.) relative uncertainty

  //***********average eff used **********
  if (isPreliminary){
    for(Int_t m=0; m<nummolt+1; m++){
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if ( fHistSpectrumStat[m]->GetBinContent(v+1)==0)  fHistSpectrumStatMCChoiceDef[m]->SetBinContent(v+1,0);
      }
      fHistSpectrumStatMCChoiceDef[m] ->Scale(0.5);
      fHistSpectrumStat[m] ->Scale(0.5);
      fHistSpectrumStat[m]->Add(    fHistSpectrumStatMCChoiceDef[m]);
      fHistSpectrumSistAll[m] ->Scale(0.5);
      fHistSpectrumSistAll[m]->Add(    fHistSpectrumStatMCChoiceDef[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistAll[m]->SetBinError(v+1, fHistSpectrumSistRelErrorAll[m]->GetBinContent(v+1)*fHistSpectrumSistAll[m]->GetBinContent(v+1));
      }
    }
  }

  TString  PathInEtaEff;
  TFile *  fileinEtaEff;
  PathInEtaEff = Dir+"/DATA"+year0+"/SystematicAnalysis" +yearEtaEff;
  if (PtBinning>0) PathInEtaEff +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    PathInEtaEff +="_"+tipo[type];
    PathInEtaEff +=Srap[israp];
    PathInEtaEff +=SSkipAssoc[SkipAssoc];
  }
  //    if (IsHighPtExtr) PathInEtaEff+="_HighPtExtr";
  PathInEtaEff+= hhCorr[ishhCorr]+"_"+RegionTypeOld[TypeAnalysis] +DataOrMC[isMC] + Form("_PtMin%.1f_IsEtaEff.root", PtTrigMin);
  fileinEtaEff = new TFile(PathInEtaEff, "");
    
  if (isReanalysisEtaEffComp){
    cout << "\n*******************************************"<<endl;
    cout<< "Reanalysis with eta-dependent efficiency " << endl;
    cout << "From file: " << PathInEtaEff << endl;
    TLegend * legendEtaEff= new TLegend(0.6, 0.6, 0.9, 0.9);
    if (!fileinEtaEff) {cout << PathInEtaEff << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummolt+1; m++){
      //    cout << " m " << m << endl;
      fHistSpectrumStatEtaEff[m]    =(TH1F*)fileinEtaEff->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      if (!fHistSpectrumStatEtaEff[m]) {cout << " I was looking for histo in " << PathInEtaEff  << endl; return;}

      cout << " first canvas " << endl;
      canvasEtaEff->cd(m+1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendEtaEff->AddEntry(fHistSpectrumStat[m], "Default", "ple");
      fHistSpectrumStat[m]->DrawClone("ep");
      fHistSpectrumSistAll[m]->SetFillStyle(0);
      fHistSpectrumSistAll[m]->DrawClone("same e2");
      StyleHisto(fHistSpectrumStatEtaEff[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendEtaEff->AddEntry(fHistSpectrumStatEtaEff[m], "Reanalysis", "ple");
      fHistSpectrumStatEtaEff[m]->Draw("same ep");
      legendEtaEff->Draw("");

      cout << " second canvas " << endl;
      canvasEtaEffRatio->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatEtaEffRatio[m] =(TH1F*) fHistSpectrumStatEtaEff[m] -> Clone(Form("fHistSpectrumEtaEffRatio_m%i", m));
      fHistSpectrumStatEtaEffRatioRef[m] =(TH1F*) fHistSpectrumStat[m] -> Clone(Form("fHistSpectrumEtaEffRatioRef_m%i", m));
      fHistSpectrumStatEtaEffRatio[m]->Divide(fHistSpectrumStat[m]);
      fHistSpectrumStatEtaEffRatioRef[m]->Divide(fHistSpectrumStat[m]);

      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
	if (v==1 && TypeAnalysis==0 && type==8) continue;
	if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
	fHistSpectrumStatEtaEffRatio[m]->SetBinError(v+1,0);
	fHistSpectrumStatEtaEffRatioRef[m]->SetBinError(v+1, sqrt(pow(fHistSpectrumStat[m]->GetBinError(v+1),2) + pow(fHistSpectrumSistAll[m]->GetBinError(v+1),2)) / fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      StyleHisto(fHistSpectrumStatEtaEffRatio[m],0.85, 1.15, Color[TypeAnalysis],33, titleX, "Ratio to preliminary",  title+SmoltLegend[m], 0, 0, 0);
      StyleHisto(fHistSpectrumStatEtaEffRatioRef[m],0.85, 1.15, 1, 1, titleX, "Ratio to preliminary",  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStatEtaEffRatioRef[m]->SetFillStyle(3001);
      fHistSpectrumStatEtaEffRatioRef[m]->SetFillColorAlpha(kGray+1, 1);
      fHistSpectrumStatEtaEffRatioRef[m]->Draw("p e2");
      fHistSpectrumStatEtaEffRatio[m]->Draw("same p");
    }//end loop m
  }//end loop to compare results obtained when considering eta dependence of efficiency

  //*************** from now on spectra obtained using eta-dependent efficiency are taken! *******
  if (isReanalysisWithEtaEff){
    if (!fileinEtaEff) {cout << PathInEtaEff << " does not exist" << endl; return;}
    for(Int_t m=0; m<nummolt+1; m++){
      fHistSpectrumStat[m]    =(TH1F*)fileinEtaEff->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      fHistSpectrumStat[m] ->SetName("fHistSpectrum_"+Smolt[m]);
      if (!fHistSpectrumStat[m]) {cout << " I was looking for histo in " << PathInEtaEff  << endl; return;}
      for (Int_t b=1; b<=       fHistSpectrumSistAll[m]->GetNbinsX(); b++){
	fHistSpectrumSistAll[m]->SetBinError(b,       fHistSpectrumSistAll[m]->GetBinError(b)/ fHistSpectrumSistAll[m]->GetBinContent(b)* fHistSpectrumStat[m]->GetBinContent(b));
	fHistSpectrumSistAll[m]->SetBinContent(b,       fHistSpectrumStat[m]->GetBinContent(b));
      }
    }
  }

  cout << "\n //*************spectra normalization ******************" << endl;
  /*
    TString SfileNormCorr ="FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue2018f1_extra_Hybrid_hK0s_vs_2018f1_extra_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0_IsParticleTrue_isEffMassSel.root";
    if(type==8)SfileNormCorr="";
  */
  TString SfileNormCorr ="FinalOutput/DATA2016/SystematicAnalysis2018f1_extra_Reco_hK0s_K0s_Eta0.8_JetMC_PtMin3.0_IsParticleTrue_IsEfficiencyMassSel.root";
  TFile* fileNormCorr;
  //  TH1F * fHistNormCorr[nummolt+1];
  TH1F * fHistNormCorr;
  TH1F * fHistNormCorrAllMult;
  TLegend * legendNormCorr= new TLegend(0.6, 0.6, 0.9, 0.9);
  //  if (isNormCorr && type==0){
  if (isNormCorr){
    cout << " Normalization factor for comparison with event generators" << endl;
    cout << " from file " << SfileNormCorr << endl;
    fileNormCorr=new TFile(SfileNormCorr,"");
    if (!fileNormCorr) {cout << "file norm not there " << endl; return;}
    /*
      for(Int_t m=0; m<nummolt+1; m++){
      fHistNormCorr[m] = (TH1F*)fileNormCorr->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+ Form("_m%i",m));
      if (!fHistNormCorr[m]) return;
      fHistNormCorr[m] ->Sumw2();
      fHistSpectrumStat[m]->Divide(fHistNormCorr[m]);
      fHistSpectrumSistAll[m]->Divide(fHistNormCorr[m]);
      cout << " spectrum " <<  fHistSpectrumStat[m]->GetNbinsX() << endl;
      cout << " norm hist "<<     fHistNormCorr[m]->GetNbinsX()<< endl;
      }
    */
    //new
    fHistNormCorr = (TH1F*)fileNormCorr->Get("fHistEventLoss");
    if (!fHistNormCorr) {cout << "histo not there " << endl; return;}
    fHistNormCorrAllMult = (TH1F*)fileNormCorr->Get("fHistEventLossMultAll");
    if (!fHistNormCorrAllMult) {cout <<"histo not there " << endl; return;}
    Float_t Factorfn=1.071;
    for(Int_t m=0; m<nummolt; m++){
      fHistSpectrumStat[m]->Scale(fHistNormCorr->GetBinContent(m+1)*Factorfn);
      fHistSpectrumSistAll[m]->Scale(fHistNormCorr->GetBinContent(m+1)*Factorfn);
    }
    fHistSpectrumStat[nummolt]->Scale(fHistNormCorrAllMult->GetBinContent(1)*Factorfn);
    fHistSpectrumSistAll[nummolt]->Scale(fHistNormCorrAllMult->GetBinContent(1)*Factorfn);
  }

  TString SfileNormCorrFC ="FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0.root";
  if (type==8)  SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue161718_hXi_Hybrid_vs_161718_hXi_Xi_Eta0.8_PtMin3.0.root";
  TFile* fileNormCorrFC;
  TH1F * fHistNormCorrFC[nummolt+1];
  TH1F * fHistNormCorrAllMultFC[nummolt+1];
  TString SfileNormCorrFCComp ="FinalOutput/DATA2016/PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_Jet_isNormCorrFullyComputed_Fixed.root";
  if (type==8) SfileNormCorrFCComp= "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_Jet_isNormCorrFullyComputed.root";
  TFile* fileNormCorrFCComp;

  if (isNormCorrFullyComputed==1){ // In this case I perform a correction with the *comnplete* normalisation factor. When this is set to 0, I take the output produced here and compare it to the preliminary results
    cout << "Normalization factor for comparison with event generators ---- fully computed!" << endl;
    cout << "From file: " << SfileNormCorrFC << endl;
    fileNormCorrFC=new TFile(SfileNormCorrFC,"");
    if (!fileNormCorrFC) {cout << "file norm not there " << endl; return;}
    for(Int_t m=0; m<nummolt+1; m++){
      fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));
      fHistSpectrumStat[m]->Divide(fHistNormCorrFC[m]);
      for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ //I only have to *scale* sist error
	//	cout <<       fHistNormCorrFC[m]->GetBinContent(b) << endl;
	fHistSpectrumSistAll[m]->SetBinContent(b,fHistSpectrumSistAll[m]->GetBinContent(b)/fHistNormCorrFC[m]->GetBinContent(b));
	fHistSpectrumSistAll[m]->SetBinError(b,fHistSpectrumSistAll[m]->GetBinError(b)/fHistNormCorrFC[m]->GetBinContent(b));
      }
    }
  }
  else if (isNormCorr && YieldComp!=0) { //comparison between preliminary results and results obtained by correcting for the full norm factor
    if (YieldComp==1) {
      cout << " comparison between preliminary results and results obtained by correcting for the full norm factor " << endl;
      if (type==0)  SfileNormCorrFCComp ="FinalOutput/DATA2016/PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_"+ RegionType[TypeAnalysis] +"_isNormCorrFullyComputed.root";
      else if (type==8) SfileNormCorrFCComp= "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_" + RegionType[TypeAnalysis]+ "_isNormCorrFullyComputed.root";
    }
    else if (YieldComp==2) {
      cout << " comparison between preliminary results and results obtained by correcting for the full norm factor and by eta eff" << endl;
      //      if (type==0) SfileNormCorrFCComp ="FinalOutput/DATA2016/PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_"+ RegionType[TypeAnalysis]+"_isNormCorrFullyComputed_isEtaEff.root";
      if (type==0) SfileNormCorrFCComp ="FinalOutput/DATA2016/PtSpectraBis_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_PtMin3.0_"+ RegionType[TypeAnalysis]+"_isNormCorrFullyComputed_isEtaEff.root";
      else if (type==8) SfileNormCorrFCComp= "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_" + RegionType[TypeAnalysis]+ "_isNormCorrFullyComputed_isEtaEff.root";
    }
    fileNormCorrFCComp=new TFile(SfileNormCorrFCComp,"");
    if (!fileNormCorrFCComp) {cout << "file norm not there " << endl; return;}
    for(Int_t m=0; m<nummolt+1; m++){
      fHistSpectrumStatNormCorr[m]    =(TH1F*)fileNormCorrFCComp->Get("fHistSpectrum_"+Smolt[m]); 
      if (!fHistSpectrumStatNormCorr[m]) {cout << " I was looking for histo fHistSpectrum_"<<Smolt[m]<< " in " << SfileNormCorrFCComp  << endl; return;}

      cout << " first canvas " << endl;
      canvasNormCorr->cd(m+1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendNormCorr->AddEntry(fHistSpectrumStat[m], "Default", "ple");
      fHistSpectrumStat[m]->DrawClone("ep");
      fHistSpectrumSistAll[m]->SetFillStyle(0);
      fHistSpectrumSistAll[m]->DrawClone("same e2");
      StyleHisto(fHistSpectrumStatNormCorr[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==0)      legendNormCorr->AddEntry(fHistSpectrumStatNormCorr[m], "Reanalysis", "ple");
      fHistSpectrumStatNormCorr[m]->Draw("same ep");
      legendNormCorr->Draw("");

      cout << " second canvas " << endl;
      canvasNormCorrRatio->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatNormCorrRatio[m] =(TH1F*) fHistSpectrumStatNormCorr[m] -> Clone(Form("fHistSpectrumNormCorrRatio_m%i", m));
      fHistSpectrumStatNormCorrRatioRef[m] =(TH1F*) fHistSpectrumStat[m] -> Clone(Form("fHistSpectrumNormCorrRatioRef_m%i", m));
      fHistSpectrumStatNormCorrRatio[m]->Divide(fHistSpectrumStat[m]);
      fHistSpectrumStatNormCorrRatioRef[m]->Divide(fHistSpectrumStat[m]);

      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (v==1 && TypeAnalysis==0 && type==0 && (m==0 || m==1 || m==3)) continue;
	if (v==1 && TypeAnalysis==0 && type==8) continue;
	if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
	fHistSpectrumStatNormCorrRatio[m]->SetBinError(v+1,0);
	fHistSpectrumStatNormCorrRatioRef[m]->SetBinError(v+1, sqrt(pow(fHistSpectrumStat[m]->GetBinError(v+1),2) + pow(fHistSpectrumSistAll[m]->GetBinError(v+1),2)) / fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      StyleHisto(fHistSpectrumStatNormCorrRatio[m],0.7, 1.3, Color[TypeAnalysis],33, titleX, "Ratio to preliminary",  title+SmoltLegend[m], 0, 0, 0);
      StyleHisto(fHistSpectrumStatNormCorrRatioRef[m],0.7, 1.3, 1, 1, titleX, "Ratio to preliminary",  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStatNormCorrRatioRef[m]->SetFillStyle(3001);
      fHistSpectrumStatNormCorrRatioRef[m]->SetFillColorAlpha(kGray+1, 1);
      fHistSpectrumStatNormCorrRatioRef[m]->Draw("p e2");
      fHistSpectrumStatNormCorrRatio[m]->Draw("same p");
    }//end loop m
  }

  cout << "...end of spectra normalization *********************** \n" << endl;

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
    cout << "mult. : " << m << endl;
    for (Int_t b= PtBinMin[m] ; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << NPtV0[b-1] << " " << fHistSpectrumStat[m]->GetBinContent(b) << ", rel stat: " << fHistSpectrumStat[m]->GetBinError(b)/ fHistSpectrumStat[m]->GetBinContent(b) << endl;
      cout << NPtV0[b-1] << " " << fHistSpectrumSistAll[m]->GetBinContent(b) << ", rel stat: " << fHistSpectrumSistAll[m]->GetBinError(b)/ fHistSpectrumSistAll[m]->GetBinContent(b) << endl;
    }
    canvasPtSpectraRelErrorAll->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistRelErrorAll[m], LimInfError, LimSupError,Color[TypeAnalysis] , 27, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumStatRelError[m], LimInfError, LimSupError,Color[TypeAnalysis] , 33, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumSistRelErrorPublished, LimInfError, LimSupError, 922, 33, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);

    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");

    canvasPtSpectraRelError->cd(m+1);
    gPad->SetLeftMargin(0.15);
    if (m==nummolt) fHistSpectrumSistRelErrorPublished->Draw("same");
    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");
    if (m==0){
      legendError->AddEntry(    fHistSpectrumStatRelError[m], "stat.", "pl");
      legendError->AddEntry(    fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
    }
    if (m==nummolt)      legendError->AddEntry(    fHistSpectrumSistRelErrorPublished, "syst. published", "l");
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
  TFitResultPtr fFitResultPtrUp[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrDown[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtr1[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrTemp[nummolt+1][numfittipo];
  TF1* fit_MTscaling[nummolt+1][numfittipo];
  TF1* fit_MTscalingUp[nummolt+1][numfittipo];
  TF1* fit_MTscalingDown[nummolt+1][numfittipo];
  TF1* fit_MTscalingBis[nummolt+1][numfittipo];
  TF1* fit_MTscalingTemp[nummolt+1][numfittipo];
  TString       nameFit[numfittipo]={"mT-scaling", "Boltzmann", "Fermi-Dirac", "Levi"};
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  Int_t factor=1;

  Float_t   AvgPt[nummolt+1]={0};
  Float_t   AvgPtFS[nummolt+1]={0};
  Float_t   AvgPtFSNum[nummolt+1]={0};
  Float_t   AvgPtFSDenom[nummolt+1]={0};
  Float_t   AvgPtUp[nummolt+1]={0};
  Float_t   AvgPtDown[nummolt+1]={0};
  Float_t   AvgPtMax[nummolt+1]={0};
  Float_t   AvgPtMin[nummolt+1]={0};
  Float_t   AvgTemp[3] = {0};
  Float_t   AvgPtFit[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitMinFit[nummolt+1]={0};
  Float_t   AvgPtFitMaxFit[nummolt+1]={0};
  Float_t   AvgPtFitUp[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitDown[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitTemp[nummolt+1][numfittipo]={0};
  Float_t   AvgPtSistErr[nummolt+1]={0};
  Float_t   AvgPtStatErr[nummolt+1]={0};
  Float_t   AvgPtSistErrFS[nummolt+1]={0};
  Float_t   AvgPtStatErrFS[nummolt+1]={0};
  Float_t   AvgPtSistErrFit[nummolt+1]={0};
  Float_t   AvgPtStatErrFit[nummolt+1][numfittipo]={0};
  Float_t   YieldExtrHighPt[nummolt+1][numfittipo]={0};
  Float_t   YieldExtrLowPt[nummolt+1][numfittipo]={0};

  Float_t    YieldExtrHighPtAvg[nummolt+1]={0};
  Float_t    YieldExtrHighPtAvgUp[nummolt+1]={0};
  Float_t    YieldExtrHighPtAvgDown[nummolt+1]={0};
  Float_t    YieldExtrLowPtAvg[nummolt+1]={0};
  Float_t    YieldExtrLowPtAvgUp[nummolt+1]={0};
  Float_t    YieldExtrLowPtAvgDown[nummolt+1]={0};
  Float_t    YieldErrStatHighPtAvg[nummolt+1]={0};
  Float_t    YieldErrStatLowPtAvg[nummolt+1]={0};
  Float_t    YieldExtr[nummolt+1]={0};
  Float_t    YieldExtrMaxLowPt[nummolt+1]={0};
  Float_t    YieldExtrMinLowPt[nummolt+1]={0};
  Float_t    YieldExtrMaxHighPt[nummolt+1]={0};
  Float_t    YieldExtrMinHighPt[nummolt+1]={0};
  Float_t    YieldExtrErrStat[nummolt+1]={0};
  Float_t    YieldExtrErrSist[nummolt+1]={0};
  Float_t    YieldExtrErrSistFourFit[nummolt+1]={0};
  Float_t    YieldExtrErrSistLowPt[nummolt+1]={0};
  Float_t    YieldExtrErrSistHighPt[nummolt+1]={0};
  Float_t    YieldExtrErrSistUp[nummolt+1]={0};
  Float_t    YieldExtrErrSistLow[nummolt+1]={0};
  Float_t    YieldErrSystExtrLowPt[nummolt+1]={0};
  Float_t    YieldErrSystExtrHighPt[nummolt+1]={0};
  Float_t    YieldSpectrum[nummolt+1]={0};
  Float_t    YieldSpectrumErrStat[nummolt+1]={0};
  Float_t    YieldSpectrumErrSist[nummolt+1]={0};
  Float_t    Yield[nummolt+1]={0};
  Float_t    YieldErrStat[nummolt+1]={0};
  Float_t    YieldErrSist[nummolt+1]={0};
  Float_t    YieldErrSistUp[nummolt+1]={0};
  Float_t    YieldErrSistLow[nummolt+1]={0};
  Float_t    YieldErrNormFactor[nummolt+1]={0};
  Float_t    YieldErrRelNormFactor[nummolt+1]={0, 0, 0, 0.015, 0.138, 0.08};

  TCanvas* canvasFitResult = new TCanvas ("canvasFitResult", "canvasFitResult", 1300, 800);
  TH1F * hFitResult[numfittipo];
  TH1F*  fHistAvgPtDistr[nummolt+1][numfittipo];

  for (Int_t typefit =0; typefit<numfittipo; typefit++){
    hFitResult[typefit] = new TH1F ("hFitResult"+nameFit[typefit], "hFitResult",nummolt, Nmolt);
  }
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
      if (m==0)    legendfit->AddEntry( fit_MTscaling[m][typefit], nameFit[typefit],"l");

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

      hFitResult[typefit] ->SetBinContent(m+1,fit_MTscaling[m][typefit]->GetChisquare()/fit_MTscaling[m][typefit]->GetNDF());

      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumStatUp[m]->SetBinContent(v+1,fHistSpectrumSistAll[m]->GetBinError(v+1)+fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumStatDown[m]->SetBinContent(v+1,-fHistSpectrumSistAll[m]->GetBinError(v+1)+fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      //fit +1sigma sistematica
      fit_MTscalingUp[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Up");      
      fit_MTscalingUp[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrUp[m][typefit]=       fHistSpectrumStatUp[m]->Fit(    fit_MTscalingUp[m][typefit],"SR0");
      fit_MTscalingUp[m][typefit]->SetRange(0,50);
      canvasPtSpectraFitUp->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatUp[m]->Draw("same");
      fit_MTscalingUp[m][typefit]->Draw("same");
      legendfit->Draw("");

      //fit -1sigma sistematica
      fit_MTscalingDown[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Down");      
      fit_MTscalingDown[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrDown[m][typefit]=       fHistSpectrumStatDown[m]->Fit(    fit_MTscalingDown[m][typefit],"SR0");
      fit_MTscalingDown[m][typefit]->SetRange(0,50);
      canvasPtSpectraFitDown->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatDown[m]->Draw("same");
      fit_MTscalingDown[m][typefit]->Draw("same");
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
      YieldExtrLowPtAvgUp[m]+= fit_MTscalingUp[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvgUp[m]+= fit_MTscalingUp[m][typefit]->Integral(UpRangeSpectrumPart[m],50);     
      YieldExtrLowPtAvgDown[m]+= fit_MTscalingDown[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvgDown[m]+= fit_MTscalingDown[m][typefit]->Integral(UpRangeSpectrumPart[m],50);     

      YieldErrStatLowPtAvg[m]+= pow(    fit_MTscaling[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
      YieldErrStatHighPtAvg[m]+= pow(    fit_MTscalingBis[m][typefit]->IntegralError(UpRangeSpectrumPart[m],50,fFitResultPtr1[m][typefit] ->GetParams(),(fFitResultPtr1[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
      //      YieldErrStatLowPtAvgUp[m]+= pow(    fit_MTscalingUp[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtrUp[m][typefit] ->GetParams(),(fFitResultPtrUp[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
      //      YieldErrStatHighPtAvgUp[m]+= pow(    fit_MTscalingUp[m][typefit]->IntegralError(UpRangeSpectrumPart[m],50,fFitResultPtrUp[m][typefit] ->GetParams(),(fFitResultPtrUp[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
      //      YieldErrStatLowPtAvgDown[m]+= pow(    fit_MTscalingDown[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtrDown[m][typefit] ->GetParams(),(fFitResultPtrDown[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
      //      YieldErrStatHighPtAvgDown[m]+= pow(    fit_MTscalingDown[m][typefit]->IntegralError(UpRangeSpectrumPart[m],50,fFitResultPtrDown[m][typefit] ->GetParams(),(fFitResultPtrDown[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2);

      //calculate average pt vs mult from fit
      Int_t numInt = 200;
      Float_t Pti=0;
      Float_t DeltaPt=	20./numInt;

      //...and systematic uncertainty of pt
      for (Int_t i=0; i<= numInt; i++){
	Pti = 20./numInt*i;
	AvgPtFit[m][typefit] += fit_MTscalingBis[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingBis[m][typefit]->Integral(0,20);
	//	cout << " AVERAGE PT at iteration " << i << " of "<< numInt << " (pt = "  << Pti << "): " << 	AvgPtFit[m][typefit] << endl;
	AvgPtFitUp[m][typefit] += fit_MTscalingUp[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingUp[m][typefit]->Integral(0,20);
	//	cout << " AVERAGE PT at iteration  (spectrum+1SigmaSist)" << i << " of "<< numInt << " (pt = "  << Pti << "): " << 	AvgPtFitUp[m][typefit] << endl;
        AvgPtFitDown[m][typefit] += fit_MTscalingDown[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingDown[m][typefit]->Integral(0,20);
	//	cout << " AVERAGE PT (spectrum-1SigmaSist) at iteration " << i << " of "<< numInt << " (pt = "  << Pti << "): " << 	AvgPtFitDown[m][typefit] << endl;
      } 

      AvgPt[m] += AvgPtFit[m][typefit];
      AvgPtUp[m] += AvgPtFitUp[m][typefit];
      AvgPtDown[m] += AvgPtFitDown[m][typefit];

      //...and statistical uncertainty of pt
      fit_MTscalingTemp[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Temp");
      Int_t IterNum =100;
      fHistAvgPtDistr[m][typefit] = new TH1F (Form("fHistAvgPtDistr_m%i_typefit%i", m, typefit), Form("fHistAvgPtDistr_m%i_typefit%i", m, typefit), 1000, LimInfPtvsMult, LimSupPtvsMult);

      for (Int_t i= 0; i<IterNum;i++ ){
	Float_t Temp=0;
	for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
	  //extract a random number
	  gRandom->SetSeed(i*10+b);
	  Temp=	  gRandom->Gaus(0, fHistSpectrumStat[m]->GetBinError(b));
	  fHistSpectrumTemp[m]->SetBinContent(b,fHistSpectrumStat[m]->GetBinContent(b) + Temp);
	  fHistSpectrumTemp[m]->SetBinError(b,fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b)*fHistSpectrumTemp[m]->GetBinContent(b));
	  //	  cout <<" stat " <<  fHistSpectrumStat[m]->GetBinContent(b) << " temp " <<  fHistSpectrumTemp[m]->GetBinContent(b) << endl;
	}
	fit_MTscalingTemp[m][typefit]->SetRange(LowRange[m], UpRange[m]);
	fFitResultPtrTemp[m][typefit]=       fHistSpectrumTemp[m]->Fit(fit_MTscalingTemp[m][typefit],"SR0Q");
	fit_MTscalingTemp[m][typefit]->SetRange(0,50);
	AvgPtFitTemp[m][typefit]=0;
	for (Int_t l=0; l<= numInt; l++){
	  Pti = 20./numInt*l;
	  AvgPtFitTemp[m][typefit] += fit_MTscalingTemp[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingTemp[m][typefit]->Integral(0,20);
	  //	cout << " AVERAGE PT at iteration " << i << " of "<< numInt << " (pt = "  << Pti << "): " << 	AvgPtFit[m][typefit] << endl;
	}
	fHistAvgPtDistr[m][typefit] ->Fill(AvgPtFitTemp[m][typefit]);
      } //end number of iterations

      AvgPtStatErrFit[m][typefit] = 	fHistAvgPtDistr[m][typefit]->GetRMS()/	fHistAvgPtDistr[m][typefit]->GetMean()* AvgPtFit[m][typefit];
      AvgPtStatErr[m]+=pow( AvgPtStatErrFit[m][typefit],2);
    } //end loop typefit

      //calculate average pt vs mult from spectrum
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      AvgPtFSNum[m] +=  fHistSpectrumStat[m]->GetBinContent(b) *  fHistSpectrumStat[m]->GetBinCenter(b)* fHistSpectrumStat[m]->GetBinWidth(b);
      AvgPtFSDenom[m] +=  fHistSpectrumStat[m]->GetBinContent(b) *  fHistSpectrumStat[m]->GetBinWidth(b);
    }

    Float_t DeltaP =0;
    Float_t AvP =0;
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      DeltaP = fHistSpectrumStat[m]->GetBinWidth(b);
      AvP = fHistSpectrumStat[m]->GetBinCenter(b);
      AvgPtStatErrFS[m] += pow((DeltaP*AvP*AvgPtFSDenom[m]-pow(DeltaP,2)*AvP*fHistSpectrumStat[m]->GetBinContent(b))/pow(AvgPtFSDenom[m],2) *fHistSpectrumStat[m]->GetBinError(b) ,2);
      AvgPtSistErrFS[m] += pow((DeltaP*AvP*AvgPtFSDenom[m]-pow(DeltaP,2)*AvP*fHistSpectrumSist[m]->GetBinContent(b))/pow(AvgPtFSDenom[m],2) *fHistSpectrumSist[m]->GetBinError(b) ,2);
    }

    AvgPtFS[m]=    AvgPtFSNum[m]/AvgPtFSDenom[m];
    AvgPtSistErrFS[m] = sqrt(      AvgPtSistErrFS[m]);
    AvgPtStatErrFS[m] = sqrt(      AvgPtStatErrFS[m]);

    //syst uncertainty associated to the choice of the fit function
    AvgPtFitMaxFit[m] = 0;
    AvgPtFitMinFit[m] = 1000;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if(AvgPtFitMaxFit[m] < AvgPtFit[m][typefit]){
	AvgPtFitMaxFit[m] = AvgPtFit[m][typefit];
      }
      if(AvgPtFitMinFit[m] > AvgPtFit[m][typefit]){
	AvgPtFitMinFit[m] = AvgPtFit[m][typefit];
      }
    }

    AvgPt[m] = AvgPt[m]/numfittipo;
    AvgPtUp[m] = AvgPtUp[m]/numfittipo;
    AvgPtDown[m] = AvgPtDown[m]/numfittipo;
    AvgPtStatErr[m] =      sqrt(AvgPtStatErr[m])/numfittipo;
    AvgPtSistErrFit[m] =  (AvgPtFitMaxFit[m] - AvgPtFitMinFit[m])/2;

    AvgTemp[0] = AvgPt[m];
    AvgTemp[1] = AvgPtUp[m];
    AvgTemp[2] = AvgPtDown[m];
    AvgPtMax[m]= 0;
    AvgPtMin[m]= 1000;

    for (Int_t i=0; i<3; i++){
      if (AvgTemp[i] >= AvgPtMax[m]) {
	AvgPtMax[m] = AvgTemp[i];
      } 
      if (AvgTemp[i] <= AvgPtMin[m]) {
	AvgPtMin[m] = AvgTemp[i];
      } 
    }
    AvgPtSistErr[m] = TMath::Abs((AvgPtMax[m]- AvgPtMin[m]))/2;

    cout << "\n\n" << endl;
    cout << "******** avg pt obtained with different fit functions " << endl;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      cout << "\n **" << endl;
      cout << " avg pt with fit " << nameFit[typefit]<< ": " << AvgPtFit[m][typefit]<< endl;
      cout << "avg pt temp " << AvgPtFitTemp[m][typefit]<< endl;
      cout << "m: " << m << ", typefit: " << nameFit[typefit] << ", RMS: " << fHistAvgPtDistr[m][typefit]->GetRMS() << ", MEAN: " << fHistAvgPtDistr[m][typefit]->GetMean()<< endl;
      cout << " bin width (check it is much smaller than stat error " << fHistAvgPtDistr[m][typefit]->GetXaxis()->GetBinWidth(1)<< endl;
    }
    cout << "*********" << endl;
    cout << " avg pt: " <<  AvgPt[m] << " avg pt up: " <<  AvgPtUp[m] << " avg pt down: " << AvgPtDown[m]<< endl;
    cout << " avg pt Max: " <<  AvgPtMax[m] << " avg pt Min: " << AvgPtMin[m]<< endl;
    cout << " AVERAGE PT for the mult  " << m << ": "	<<AvgPt[m]<< " +- " << AvgPtStatErr[m] << " (stat.) +- " << AvgPtSistErr[m] << " (syst. no function choice) +- "<<AvgPtSistErrFit[m] <<" (syst. function choice)"<<  endl;
    cout << " Rel errors on avg pt for the mult  " << m << ": " << AvgPtStatErr[m]/AvgPt[m] << " (stat.) +- " << AvgPtSistErr[m]/AvgPt[m] << " (syst. no function choice) +- "<< AvgPtSistErrFit[m]/AvgPt[m]<< " (syst. function choice)"<<  endl;
    cout << " av pt from spectrum " << endl;
    cout << " AVERAGE PT for the mult  " << m << ": "	<<AvgPtFS[m]<< " +- " << AvgPtStatErrFS[m] << " (stat.) +- " << AvgPtSistErrFS[m] << " (syst. no function choice) +- "<< endl; 
    // end of calculation of pt vs mult

    YieldExtrHighPtAvg[m]=YieldExtrHighPtAvg[m]/numfittipo;
    YieldExtrLowPtAvg[m]=YieldExtrLowPtAvg[m]/numfittipo;
    YieldExtrHighPtAvgUp[m]=YieldExtrHighPtAvgUp[m]/numfittipo;
    YieldExtrLowPtAvgUp[m]=YieldExtrLowPtAvgUp[m]/numfittipo;
    YieldExtrHighPtAvgDown[m]=YieldExtrHighPtAvgDown[m]/numfittipo;
    YieldExtrLowPtAvgDown[m]=YieldExtrLowPtAvgDown[m]/numfittipo;
    YieldErrStatHighPtAvg[m]=sqrt(YieldErrStatHighPtAvg[m])/numfittipo;
    YieldErrStatLowPtAvg[m]=sqrt(YieldErrStatLowPtAvg[m])/numfittipo;


    YieldExtr[m] =     YieldExtrHighPtAvg[m]+    YieldExtrLowPtAvg[m];
    YieldExtrErrStat[m] = sqrt(    pow(YieldErrStatHighPtAvg[m],2)+    pow(YieldErrStatLowPtAvg[m],2));

    YieldExtrErrSistHighPt[m] = (YieldExtrMaxHighPt[m]-YieldExtrMinHighPt[m])/2;
    YieldExtrErrSistLowPt[m] = (YieldExtrMaxLowPt[m]-YieldExtrMinLowPt[m])/2;
    YieldErrSystExtrLowPt[m] = (YieldExtrLowPtAvgUp[m]- YieldExtrLowPtAvgDown[m])/2;
    YieldErrSystExtrHighPt[m] = (YieldExtrHighPtAvgUp[m]- YieldExtrHighPtAvgDown[m])/2;

    YieldExtrErrSistUp[m] = YieldExtrMaxLowPt[m]-YieldExtrLowPtAvg[m];
    YieldExtrErrSistLow[m] = YieldExtrLowPtAvg[m];
    //    YieldExtrErrSist[m]=sqrt(pow(    YieldExtrErrSistUp[m],2) + pow(    YieldExtrErrSistLow[m],2)) ;
    YieldExtrErrSistFourFit[m]=sqrt(pow(    YieldExtrErrSistHighPt[m],2) + pow(    YieldExtrErrSistLowPt[m],2));
    YieldExtrErrSist[m]=sqrt(pow(  YieldErrSystExtrLowPt[m],2) + pow(  YieldErrSystExtrHighPt[m] ,2)) ;

    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      YieldSpectrum[m] +=  ( fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1));
      YieldSpectrumErrStat[m] += pow ( fHistSpectrumStat[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1),2);
      //errlin
      YieldSpectrumErrSist[m] += pow ( fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1),2);
      //errlin      YieldSpectrumErrSist[m] += fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1);
     
    }
    YieldSpectrumErrStat[m]=sqrt( YieldSpectrumErrStat[m]);
    //errlin 
    YieldSpectrumErrSist[m]=sqrt( YieldSpectrumErrSist[m]);

    Yield[m] =  YieldSpectrum[m]+YieldExtr[m];
    YieldErrNormFactor[m] = YieldErrRelNormFactor[m]*Yield[m];
    YieldErrStat[m] = sqrt(    pow(YieldSpectrumErrStat[m],2) + pow(YieldExtrErrStat[m],2));
    //errlin    YieldErrSist[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistHighPt[m],2) +  pow(YieldExtrErrSistLowPt[m],2) + pow(YieldSpectrumErrOOJSub[m],2) +  pow(YieldErrSystExtrLowPt[m],2) + pow( YieldErrSystExtrHighPt[m],2));
    //way1      YieldErrSist[m] =   YieldSpectrumErrSist[m] + YieldExtrErrSistHighPt[m] +  YieldExtrErrSistLowPt[m] + YieldSpectrumErrOOJSub[m] +  YieldErrSystExtrLowPt[m] +  YieldErrSystExtrHighPt[m];
    
    //prelimnary error YieldErrSist[m] =   YieldSpectrumErrSist[m] + YieldSpectrumErrOOJSub[m] +  YieldErrSystExtrLowPt[m] +  YieldErrSystExtrHighPt[m];
    YieldErrSist[m] = sqrt(pow(YieldSpectrumErrSist[m],2) + pow(YieldSpectrumErrOOJSub[m],2) +  pow(YieldErrSystExtrLowPt[m],2) + pow( YieldErrSystExtrHighPt[m],2));
    if (isNormCorr)    YieldErrSist[m]= sqrt(pow(      YieldErrSist[m],2) + pow(YieldExtrErrSistLowPt[m],2) + pow(YieldExtrErrSistHighPt[m],2) + pow(YieldErrNormFactor[m],2));
    else     YieldErrSist[m]= sqrt(pow(      YieldErrSist[m],2) + pow(YieldExtrErrSistLowPt[m],2) + pow(YieldExtrErrSistHighPt[m],2) );

    YieldErrSistUp[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistUp[m],2) +  pow(YieldExtrErrSistHighPt[m],2)  + pow(YieldSpectrumErrOOJSub[m],2) + pow( YieldErrSystExtrHighPt[m],2));
    YieldErrSistLow[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistLow[m],2) +  pow(YieldExtrErrSistHighPt[m],2)  + pow(YieldSpectrumErrOOJSub[m],2) +  pow(YieldErrSystExtrLowPt[m],2));

    if (ZeroYieldLowPt){
      Yield[m] =  YieldSpectrum[m]+ YieldExtrHighPtAvg[m];
      YieldErrStat[m] = sqrt(    pow(YieldSpectrumErrStat[m],2) + pow(YieldErrStatHighPtAvg[m],2));
      YieldErrSist[m] = sqrt(    pow(YieldSpectrumErrSist[m],2) + pow(YieldExtrErrSistHighPt[m],2) + pow(YieldSpectrumErrOOJSub[m],2) + pow( YieldErrSystExtrHighPt[m],2));
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
    cout << "from extr: " <<     YieldExtr[m] << "+- " <<     YieldExtrErrStat[m]<< " (stat)  +- " << YieldExtrErrSistFourFit[m]<< " (sist 4 fit) " << YieldExtrErrSist[m] << " (sist extr) "<< endl;
    cout << "rel uncertainties " << "syst on low pt extr " <<  YieldExtrErrSist[m]/Yield[m] << " syst on choice fit function " <<  YieldExtrErrSistFourFit[m]/Yield[m] <<endl;
    cout << " extr at low pt " << YieldExtrLowPtAvg[m]<<" " << "err extrap high pt " << YieldExtrHighPtAvg[m]<< endl; 
    cout << " extr at low pt " << YieldExtrLowPtAvg[m]<<" " << "err extrap low pt " << YieldExtrErrSistLow[m]<< endl; 
    if (TypeAnalysis==0 && type==8)     cout << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) + "<< YieldErrSistUp[m]<< " - " << YieldErrSistLow[m] << " (sist) "<< endl;

    cout << " fraction of yield from low pt extrapolation " <<  YieldExtrLowPtAvg[m]/Yield[m] << endl;
    cout << " fraction of yield from high pt extrapolation " <<  YieldExtrHighPtAvg[m]/Yield[m] << endl;
    cout << "fit range " << LowRange[m] << "- " << UpRange[m] << endl;
    cout << "!extrapolation range " << LowRangeSpectrumPart[m] << "- " << UpRangeSpectrumPart[m] << endl;
  }//end loop m


  TLegend *legendYield=new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend *legendPtvsMult=new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend *legendYieldErr=new TLegend(0.6, 0.6, 0.9, 0.9);

  TF1 * pol0 = new TF1 ("pol0", "pol0", 0, 30);
  TString titleYieldX="dN_{ch}/d#eta";
  TString titleYieldYType[2]={"N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi", "N_{#Xi}/N_{Trigg} 1/#Delta#eta #Delta#phi"};
  TString titlePtvsMultYType[2]={"<p^{K^{0}_{S}}_{T}> (GeV/c)", "<p^{#Xi}_{T} (GeV/c)>"};
  TString titleYieldY;
  if(type==0) titleYieldY=titleYieldYType[0];
  else if(type==8) titleYieldY=titleYieldYType[1];
  TString titlePtvsMultY;
  if(type==0) titlePtvsMultY=titlePtvsMultYType[0];
  else if(type==8) titlePtvsMultY=titlePtvsMultYType[1];
  TString titleYield[3]={"In-jet", "Out-of-jet", "Inclusive"};
  TH1F* fHistYieldStat=new TH1F ("fHistYieldStatD","fHistYieldStatD",150,0,30);
  TH1F* fHistYieldSist=new TH1F ("fHistYieldSist","fHistYieldSist",150,0,30);
  TH1F* fHistPtvsMultStat = new TH1F ("fHistPtvsMultStat","fHistPtvsMultStat",150,0,30);
  TH1F* fHistPtvsMultSist = new TH1F ("fHistPtvsMultSist","fHistPtvsMultSist",150,0,30);
  TH1F* fHistPtvsMultStatFromSpectrum = new TH1F ("fHistPtvsMultStatFromSpectrum","fHistPtvsMultStatFromSpectrum",150,0,30);
  TH1F* fHistPtvsMultSistFromSpectrum = new TH1F ("fHistPtvsMultSistFromSpectrum","fHistPtvsMultSistFromSpectrum",150,0,30);
  TH1F* fHistYieldSistNoExtr=new TH1F ("fHistYieldSistNoExtr","fHistYieldSistNoExtr",150,0,30);
  TH1F* fHistYieldStatRelErr=new TH1F ("fHistYieldStatRelErr","fHistYieldStatRelErr",150,0,30);
  TH1F* fHistYieldSistRelErr=new TH1F ("fHistYieldSistRelErr","fHistYieldSistRelErr",150,0,30);
  TH1F* fHistYieldSistNoExtrRelErr=new TH1F ("fHistYieldSistNoExtrRelErr","fHistYieldSistNoExtrRelErr",150,0,30);
  TH1F* fHistYieldSistLowRelErr=new TH1F ("fHistYieldSistLowRelErr","fHistYieldSistLowRelErr",150,0,30);
  TH1F* fHistYieldSistUpRelErr=new TH1F ("fHistYieldSistUpRelErr","fHistYieldSistUpRelErr",150,0,30);

  TH1F* fHistYieldStatEtaEff=new TH1F ("fHistYieldStatEtaEff","fHistYieldStatEtaEff",150,0,30);
  TH1F* fHistYieldStatEtaEffRatio;
  TH1F* fHistYieldStatEtaEffRatioRef;
  TH1F* fHistYieldStatNormCorr=new TH1F ("fHistYieldStatNormCorr","fHistYieldStatNormCorr",150,0,30);
  TH1F* fHistYieldStatNormCorrRatio;
  TH1F* fHistYieldStatNormCorrRatioRef;

  Float_t   mult[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};

  Int_t LimSupChi=120;
  if (type==1) LimSupChi=100;
  canvasFitResult->cd();
  for (Int_t typefit=0; typefit<numfittipo; typefit++){
    StyleHisto( hFitResult[typefit], 0,LimSupChi, ColorFit[typefit], 33, "Multiplicity class", "#chi^{2}/NDF", "#chi^{2} vs multiplicity", 0, 0, 0);
    hFitResult[typefit]->GetYaxis()->SetTitleOffset(1.2);
    hFitResult[typefit]->Draw("same");
    if (typefit==0)    legendfit->Draw("");
  }

  for(Int_t m=0; m<nummolt; m++){
    fHistPtvsMultStat->SetBinContent(fHistPtvsMultStat->FindBin(mult[m]),AvgPt[m]);
    fHistPtvsMultSist->SetBinContent(fHistPtvsMultStat->FindBin(mult[m]),AvgPt[m]);
    fHistPtvsMultStat->SetBinError(fHistPtvsMultStat->FindBin(mult[m]),AvgPtStatErr[m]);
    fHistPtvsMultSist->SetBinError(fHistPtvsMultStat->FindBin(mult[m]),sqrt(pow(AvgPtSistErr[m],2) + pow(AvgPtSistErrFit[m],2)));
    fHistPtvsMultStatFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(mult[m]),AvgPtFS[m]);
    fHistPtvsMultSistFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(mult[m]),AvgPtFS[m]);
    fHistPtvsMultStatFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(mult[m]),AvgPtStatErrFS[m]);
    fHistPtvsMultSistFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(mult[m]),AvgPtSistErrFS[m]);

  }

  for(Int_t m=0; m<nummolt; m++){
    fHistYieldStat->SetBinContent(fHistYieldStat->FindBin(mult[m]),Yield[m]);
    fHistYieldStat->SetBinError(fHistYieldStat->FindBin(mult[m]),YieldErrStat[m]);
    fHistYieldSist->SetBinContent(fHistYieldSist->FindBin(mult[m]),Yield[m]);
    fHistYieldSist->SetBinError(fHistYieldSist->FindBin(mult[m]),YieldErrSist[m]);
    fHistYieldSistNoExtr->SetBinContent(fHistYieldSist->FindBin(mult[m]),Yield[m]);
    if (isNormCorr)  fHistYieldSistNoExtr->SetBinError(fHistYieldSist->FindBin(mult[m]), sqrt(pow(YieldSpectrumErrSist[m],2) +pow(YieldErrNormFactor[m],2)));
    else   fHistYieldSistNoExtr->SetBinError(fHistYieldSist->FindBin(mult[m]), YieldSpectrumErrSist[m]);

    fHistYieldStatRelErr->SetBinContent(fHistYieldStat->FindBin(mult[m]),YieldErrStat[m]/Yield[m]);
    fHistYieldStatRelErr->SetBinError(fHistYieldStat->FindBin(mult[m]),0);
    fHistYieldSistRelErr->SetBinContent(fHistYieldSist->FindBin(mult[m]),YieldErrSist[m]/Yield[m]);
    fHistYieldSistRelErr->SetBinError(fHistYieldSist->FindBin(mult[m]),0);
    if (isNormCorr)  fHistYieldSistNoExtrRelErr->SetBinContent(fHistYieldSist->FindBin(mult[m]), sqrt(pow(YieldSpectrumErrSist[m],2) + pow(YieldErrNormFactor[m],2))/Yield[m]);
    else   fHistYieldSistNoExtrRelErr->SetBinContent(fHistYieldSist->FindBin(mult[m]), YieldSpectrumErrSist[m]/Yield[m]);
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

  fHistYieldStat->DrawClone("e");
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

  if (isReanalysisEtaEffComp && !isNormCorr){ //the yield obtained from reanalysis has not been multiplied by norm factor
    fHistYieldStatEtaEff    =(TH1F*)fileinEtaEff->Get("fHistYieldvsErrSoloStat"); 
    fHistYieldStatEtaEff->SetName("fHistYieldvsErrSoloStatEtaEff");
    if (!fHistYieldStatEtaEff) {cout << " I was looking for histo in " << PathInEtaEff  << endl; return;}
    //    fHistYieldStatEtaEff->Rebin(2); //binning is different

    canvasYieldEtaEff->cd();
    StyleHisto(fHistYieldStat, LimInfYield, LimSupYield, Color[TypeAnalysis], 33, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity",1, 0, 25);
    StyleHisto(fHistYieldStatEtaEff, LimInfYield, LimSupYield, 1, 33, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity",1, 0, 25);
    fHistYieldStat->DrawClone("e");
    fHistYieldStatEtaEff  ->DrawClone("e same");  

    canvasYieldEtaEffRatio->cd();
    fHistYieldStatEtaEffRatio = (TH1F*) fHistYieldStatEtaEff->Clone("fHistYieldvsErrSoloStatEtaEffRatio");
    fHistYieldStatEtaEffRatioRef = (TH1F*) fHistYieldStat->Clone("fHistYieldvsErrSoloStatEtaEffRatioRef");
    fHistYieldStatEtaEffRatioRef->Divide(fHistYieldStat);
    fHistYieldStatEtaEffRatio->Divide(fHistYieldStat); 
    cout << "\n\n*********************\n\n\n" << endl;
    for (Int_t b=1; b<= fHistYieldStatEtaEff->GetNbinsX(); b++){
      fHistYieldStatEtaEffRatio->SetBinError(b, 0);
      if (fHistYieldStat->GetBinContent(b) !=0){
	fHistYieldStatEtaEffRatioRef->SetBinError(b, sqrt(pow(fHistYieldStat->GetBinError(b),2) + pow(fHistYieldSist->GetBinError(b),2) ) / fHistYieldStat->GetBinContent(b));
      }
    }

    StyleHisto(fHistYieldStatEtaEffRatio, 0.8, 1.2, Color[TypeAnalysis], 33, titleYieldX,"Ratio", titleYield[TypeAnalysis] + " yield ratio to preliminary vs multiplicity",1, 0, 25);
    StyleHisto(fHistYieldStatEtaEffRatioRef, 0.8, 1.2, 1, 1, titleYieldX,"Ratio", titleYield[TypeAnalysis] + " yield ratio to preliminary vs multiplicity",1, 0, 25);
    fHistYieldStatEtaEffRatioRef->SetFillStyle(3001);
    fHistYieldStatEtaEffRatioRef->SetFillColorAlpha(kGray+1, 1);
    fHistYieldStatEtaEffRatioRef->DrawClone("p e2");
    fHistYieldStatEtaEffRatio   ->DrawClone("p same");  

  }

  if (isNormCorr && YieldComp!=0){ //the yield obtained from reanalysis has not been multiplied by norm factor
    for (Int_t b=1; b<= fHistYieldStat->GetNbinsX(); b++){
      if (fHistYieldStat->GetBinContent(b)!=0){
	cout <<     fHistYieldStat->GetBinContent(b) << endl;
      }
    }
    fHistYieldStatNormCorr=(TH1F*)fileNormCorrFCComp->Get("fHistYieldStat"); 
    if (!fHistYieldStatNormCorr) {cout << " I was looking for histo in " << SfileNormCorrFCComp  << endl; return;}
    fHistYieldStatNormCorr->SetName("fHistYieldStatNormCorrFullyComputed");
    for (Int_t b=1; b<= fHistYieldStat->GetNbinsX(); b++){
      if (fHistYieldStat->GetBinContent(b)!=0){
	cout <<     fHistYieldStatNormCorr->GetBinContent(b) << endl;
      }
    }

    canvasYieldNormCorr->cd();
    StyleHisto(fHistYieldStat, LimInfYield, LimSupYield, Color[TypeAnalysis], 33, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity",1, 0, 25);
    StyleHisto(fHistYieldStatNormCorr, LimInfYield, LimSupYield, 1, 33, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity",1, 0, 25);
    fHistYieldStat->DrawClone("e");
    fHistYieldStatNormCorr  ->DrawClone("e same");  

    canvasYieldNormCorrRatio->cd();
    fHistYieldStatNormCorrRatio = (TH1F*) fHistYieldStatNormCorr->Clone("fHistYieldvsErrSoloStatNormCorrRatio");
    fHistYieldStatNormCorrRatioRef = (TH1F*) fHistYieldStat->Clone("fHistYieldvsErrSoloStatNormCorrRatioRef");
    fHistYieldStatNormCorrRatioRef->Divide(fHistYieldStat);
    fHistYieldStatNormCorrRatio->Divide(fHistYieldStat); 
    cout << "\n\n*********************\n\n\n" << endl;
    for (Int_t b=1; b<= fHistYieldStatNormCorr->GetNbinsX(); b++){
      fHistYieldStatNormCorrRatio->SetBinError(b, 0);
      if (fHistYieldStat->GetBinContent(b) !=0){
	fHistYieldStatNormCorrRatioRef->SetBinError(b, sqrt(pow(fHistYieldStat->GetBinError(b),2) + pow(fHistYieldSist->GetBinError(b),2) ) / fHistYieldStat->GetBinContent(b));
      }
    }

    StyleHisto(fHistYieldStatNormCorrRatio, 0.8, 1.2, Color[TypeAnalysis], 33, titleYieldX,"Ratio", titleYield[TypeAnalysis] + " yield ratio to preliminary vs multiplicity",1, 0, 25);
    StyleHisto(fHistYieldStatNormCorrRatioRef, 0.8, 1.2, 1, 1, titleYieldX,"Ratio", titleYield[TypeAnalysis] + " yield ratio to preliminary vs multiplicity",1, 0, 25);
    fHistYieldStatNormCorrRatioRef->SetFillStyle(3001);
    fHistYieldStatNormCorrRatioRef->SetFillColorAlpha(kGray+1, 1);
    fHistYieldStatNormCorrRatioRef->DrawClone("p e2");
    fHistYieldStatNormCorrRatio   ->DrawClone("p same");  

  }
  fHistYieldStat->SetName("fHistYieldStat");
  fHistYieldStat->SetTitle("fHistYieldStat");

  canvasPtvsMult->cd();
  StyleHisto(fHistPtvsMultStat, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 1, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity "+ RegionTypeNew[TypeAnalysis],1, 0, 25);
  StyleHisto(fHistPtvsMultSist, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 1, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" + RegionTypeNew[TypeAnalysis], 1, 0, 25);
  StyleHisto(fHistPtvsMultStatFromSpectrum, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 4, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity "+ RegionTypeNew[TypeAnalysis],1, 0, 25);
  StyleHisto(fHistPtvsMultSistFromSpectrum, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 4, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" + RegionTypeNew[TypeAnalysis], 1, 0, 25);
     
  fHistPtvsMultStat->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSist->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultStatFromSpectrum->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSistFromSpectrum->GetYaxis()->SetTitleOffset(1.2);

  if (TypeAnalysis==0){
    fHistPtvsMultStat->Fit(pol0, "R0");
  }

  fHistPtvsMultStat->Draw("e");
  fHistPtvsMultSist->SetFillStyle(0);
  fHistPtvsMultSist->Draw("same e2");
  fHistPtvsMultStatFromSpectrum->Draw("same e");
  fHistPtvsMultSistFromSpectrum->SetFillStyle(0);
  fHistPtvsMultSistFromSpectrum->Draw("same e2");
  gYield->SetFillStyle(0);
  legendPtvsMult->AddEntry(fHistPtvsMultSist, "syst.", "ef");
  legendPtvsMult->AddEntry(fHistPtvsMultStat, "stat.", "pel");
  legendPtvsMult->Draw("");

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

  fileout->WriteTObject(canvasEtaEff);
  fileout->WriteTObject(canvasEtaEffRatio);
  if (YieldComp!=0){
    fileout->WriteTObject(canvasNormCorr);
    fileout->WriteTObject(canvasNormCorrRatio);
    fileout->WriteTObject(canvasYieldNormCorr);
    fileout->WriteTObject(canvasYieldNormCorrRatio);
  }
  fileout->WriteTObject(canvasYield);
  fileout->WriteTObject(canvasYieldEtaEff);
  fileout->WriteTObject(canvasYieldEtaEffRatio);
  fileout->WriteTObject(canvasYieldErr);
  fileout->WriteTObject(canvasPtvsMult);
  fileout->WriteTObject(fHistYieldStat);
  fileout->WriteTObject(fHistYieldSist);
  fileout->WriteTObject(fHistYieldSistNoExtr);
  fileout->WriteTObject(fHistPtvsMultStat);
  fileout->WriteTObject(fHistPtvsMultSist);
  fileout->WriteTObject(fHistPtvsMultStatFromSpectrum);
  fileout->WriteTObject(fHistPtvsMultSistFromSpectrum);
  fileout->WriteTObject(fHistYieldStatRelErr);
  fileout->WriteTObject(fHistYieldSistRelErr);
  fileout->WriteTObject(fHistYieldSistNoExtrRelErr);
  for(Int_t m=0; m<nummolt+1; m++){
    fileout->WriteTObject( fHistAvgPtDistr[m][0]);
  }
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
  fileout->WriteTObject(canvasMCChoiceDef);
  fileout->WriteTObject(canvasBarlowMCChoiceDef);
  fileout->WriteTObject(canvasFakeSBDef);
  fileout->WriteTObject(canvasBarlowFakeSBDef);
  fileout->WriteTObject(canvasPtSpectraAll);
  fileout->WriteTObject(canvaspol0);
  fileout->WriteTObject(canvasPtSpectra);
  fileout->WriteTObject(canvasPtSpectraFit);
  fileout->WriteTObject(canvasPtSpectraFitUp);
  fileout->WriteTObject(canvasPtSpectraFitDown);
  fileout->WriteTObject(canvasPtSpectraFitBis);
  fileout->WriteTObject(canvasBarlow);
  fileout->WriteTObject(canvasBarlowpol0);
  fileout->WriteTObject(canvasPtSpectraRelError);
  fileout->WriteTObject(canvasPtSpectraRelErrorAll);
  fileout->WriteTObject(canvasFitResult);
  for(Int_t m=0; m<nummolt+1; m++){
    fileout->WriteTObject(     fHistSpectrumSistRelErrorMCclosure[m]);
    fileout->WriteTObject(    fHistSpectrumSistFakeSB[m]);
    fileout->WriteTObject(      fHistSpectrumSistAll[m]);
    fileout->WriteTObject(      fHistSpectrumStat[m]);
    fileout->WriteTObject(    fHistSpectrumSistRelErrorMCChoice[m]);
    fileout->WriteTObject(    fHistSpectrumSistRelErrorFakeSB[m]);
    if (isNormCorr)    fileout->WriteTObject(      fHistNormCorr);
  }

  //  TString DirPicture = "PictureForNote/";
  //if (!isNormCorr) DirPicture = "PictureForNote/IsNotNormCorr_";
  TString DirPicture = PathOutPictures;
  canvasEtaEff->SaveAs(DirPicture+"PtSpectraEtaEffComp"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasEtaEffRatio->SaveAs(DirPicture+"PtSpectraEtaEffCompRatio"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  if (YieldComp!=0){
    canvasNormCorr->SaveAs(DirPicture+"PtSpectraNormCorrComp.pdf");
    canvasNormCorrRatio->SaveAs(DirPicture+"PtSpectraNormCorrCompRatio.png");
    canvasYieldNormCorr->SaveAs(DirPicture+"YieldvsMultNormCorrComp.png");
    canvasYieldNormCorrRatio->SaveAs(DirPicture+"YieldvsMultNormCorrCompRatio.png");
  }
  if (ZeroYieldLowPt && TypeAnalysis==0){
    canvasYield->SaveAs(DirPicture+"YieldvsMult"+tipo[type]+RegionType[TypeAnalysis]+"_ZeroYieldLowPt.pdf");
    canvasYieldErr->SaveAs(DirPicture+"YieldvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+"_ZeroYieldLowPt.pdf");
  }
  else {
    canvasYield->SaveAs(DirPicture+"YieldvsMult"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
    canvasYieldErr->SaveAs(DirPicture+"YieldvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  }
  canvasYieldEtaEff->SaveAs(DirPicture+"YieldvsMultEtaEffComp"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasYieldEtaEffRatio->SaveAs(DirPicture+"YieldvsMultEtaEffCompRatio"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtvsMult->SaveAs(DirPicture+"AvgPtvsMult"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectra->SaveAs(DirPicture+"PtSpectra"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraFit->SaveAs(DirPicture+"PtSpectraFit"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraFit->SaveAs(DirPicture+"PtSpectraFitBis"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraRelError->SaveAs(DirPicture+"PtSpectraRelErr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraRelErrorAll->SaveAs(DirPicture+"PtSpectraRelErrAll"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasFitResult->SaveAs(DirPicture+"ChiSquare"+tipo[type]+RegionType[TypeAnalysis]+".pdf");

  fileout->Close();

  if (TypeAnalysis==0)  cout <<" q " <<  pol0->GetParameter(0) << " red chisq " << pol0->GetChisquare()/pol0->GetNDF() << endl;

  if (isReanalysisWithEtaEff) cout << "\n\n**The syst uncertainty associated to the non-closure has been removed! " << endl;

  cout << "\nStarting from the file(s):"<< endl;
  cout  << "->for default spectra: ";
  if (isReanalysisWithEtaEff) cout  << PathInEtaEff<< endl;
  else cout  << PathInDef<< endl;
  if (type==0 && TypeAnalysis==0)  cout << "->to get uncertainty related to OOJ subtraction: " << PathInpol0<< endl;
  if (isReanalysisWithEtaEff || isReanalysisEtaEffComp) cout << "->to get the pt spectra corrected by eta-dependent efficiency " << PathInEtaEff << endl;
  cout << "To get uncertainty related to sidebands: "<< PathInFakeSBDef<< endl;
  if (isPreliminary) cout << "To get uncertainty related to choice of MC " <<   PathInMCChoiceDef<< endl;

  if (YieldComp!=0)  cout << "Path where spectra to be compared with default ones are found: " << SfileNormCorrFCComp << endl;
  else if (isReanalysisEtaEffComp)  cout << "Path where spectra to be compared with default ones are found: " <<PathInEtaEff << endl;

  if (isNormCorrFullyComputed) cout << "File from where normalisation factor is taken: " << SfileNormCorrFC << endl;

  if (!isPreliminary){
    cout << "\nWith respect to preliminaries: " << endl;
    cout << "1) The final spectra are obtained using PYTHIA8 efficiency, and are not an average of the spectra obtaiend with EPOS and PYTHIS8, becasue I couldn't compute 2D efficiency using EPOS. The relative uncertainty associated to the choice of MC is taken from Preliminary results, namely from the file:\n" << PathMCChoiceSyst<< endl;
    cout << "Syst. associated to fake K0s/Xi subtraction taken from Preliminaries, from the file: " << PathFakeSBSyst << endl; 
  }

  cout << "\nI have created the file:\n" << stringout << "\n" <<endl;
}

