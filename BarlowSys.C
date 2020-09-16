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

void BarlowSys(Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Int_t TypeAnalysis=3, Bool_t isMC=0,   Int_t israp=0,TString year=/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"*/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"*//*"Run2DataRed_hXi"/"2016kehjl_hK0s"*/,  TString Path1 ="_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, Int_t avoidthissyst=3,Int_t avoidthissystbis=20,Int_t avoidthissysttris=21, Int_t type=8,  Bool_t isEfficiency=1,   TString Dir="FinalOutput",TString year0="2016", Int_t sysTrigger=0, Float_t numSigmaCorr=2, Bool_t MasterThesisAnalysis=0, Bool_t isEnlargedDeltaEtaPhi=0,Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1,   Bool_t SkipAssoc=1, Int_t MultBinning=0, Int_t PtBinning=0, TString FitFixed="Fermi-Dirac"){

  if (ishhCorr && type!=0){
    type=0;
  }// {cout << " for hh correlation type should be set to zero in order to get the correct files " << endl; return;}

  const   Int_t NumberTypeAnalysis=8;
  cout << "here's the meaning of different values of TypeAnalysis:" << endl;
  cout << 0 << "in-jet production " << endl;
  cout << 1 << "out-of-jet production  (Delta Phi between 1 and 2)" << endl;
  cout << 2 << "inclusive production (from JetBulkEffCorr)" << endl;
  cout << 3 << "inclusive production not scaled by DeltaEta and DeltaPhi (for comparison with published data) (from JetBulkEffCorr)" << endl;
  cout << 4 << "away side (from JetBulkEffCorr)" << endl;
  cout << 5 << "out-of-jet production (including away side, from JetBulkEffCorr) " << endl;

  if (TypeAnalysis>7) return;
  Bool_t isBulk=0; Bool_t isTotal=0;
  if (TypeAnalysis==0) {isBulk=0; isTotal=0;}
  if (TypeAnalysis==1) {isBulk=1; isTotal=0;}
  if (TypeAnalysis==2) {isBulk=0; isTotal=1;}

  if (year != "2018f1_extra" && year != "2016k" && year != "2018f1_extra_onlyTriggerWithHighestPt" && year != "2016k_onlyTriggerWithHighestPt") {
    //    cout << "output file should be changed: it must include the year name to avoid overwriting output files " << endl;
    //    return;
  } 
  if (isBulk && isTotal) return;
  if (ishhCorr) avoidthissyst=25; //è un valore che non esiste
  gStyle->SetOptStat(0);

  if (sysTrigger!=0){
    cout << "sysTrigger must be zero" << endl;
    return;
  }

  if (isTotal) {
    avoidthissystbis=10;
    avoidthissysttris=11;
  }

  TString hhCorr[2]={"", "_hhCorr"};
  TFile *fileinbis;
  TFile *filein;
  TString PathIn1;
  TString file = year;
  if (PtBinning!=0)  file+=Form("_PtBinning%i", PtBinning);
  file+=Path1;
  TString PathInBis =  "FinalOutput/AnalysisResults" + year/* + Path1*/  + ".root";
  if (isMC && isEfficiency)  PathInBis =  "FinalOutput/AnalysisResults" + year  + "_MCEff" /*+ Path1*/ +".root";
  if (ishhCorr && !isMC)  PathInBis =  "FinalOutput/AnalysisResults" + year  +"_hhCorr" +Path1 + ".root";
  if (ishhCorr && isMC)  PathInBis =  "FinalOutput/AnalysisResults" + year  + "_hhCorr_MCEff" + Path1 + ".root";
  cout << "path in (task output) "<< PathInBis << endl;
  fileinbis=new TFile(PathInBis,"");
  if (!fileinbis){cout << PathInBis << " does not exist " << endl; return;}

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9
  const Int_t numPtTrigger=1;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;

  const Int_t numSysV0 = 7;//num eff sistematici associati alla selezione delle particelle associate
  const Int_t numSysV0hh = 3; //num eff sistematici associati alla selezione delle particelle associate (nel caso ishhCorr)
  const Int_t numSyst  = 13; //all systematics together(numsysv0 + numsystphi + systematic related to bin counting region (1))
  const Int_t numSysthh  = 9;
  const Int_t numsystPhi  = 5; //includo tutti i sistematici anche se i primi 3 non vengono usati in analisi hh (si tratta sia dei sistematici legati alla scelta delle regioni di massa invariante ('centrale' e sidebands) (non per hh) che dei sistematici legati alla scelta delle regioni in DeltaEta di jet e out of jet) 
  cout << "ok " << endl;
  Int_t PtV0Min = 0; //0 el
  if (type>0 || ishhCorr )   PtV0Min =1;
  if (type==0 && PtBinning!=0) PtV0Min=1;
  Int_t numSysV0Global = 0; 
  Int_t numSystGlobal  = 0;
  Int_t numsystPhiGlobal  = 0;
  if (!ishhCorr){
    /*    if (type==0){
	  numSysV0Global = numSysV0; 
	  numSystGlobal  = numSyst;
	  numsystPhiGlobal  = numsystPhi;
	  }*/
    //    if (type!=0){    
    numSysV0Global = 1; 
    numSystGlobal  = 2;
    numsystPhiGlobal  = 0;
    //    }
  }
  if (ishhCorr){
    numSysV0Global = numSysV0hh;
    numSystGlobal  = numSyst;
    numsystPhiGlobal  = numsystPhi;
  }

  cout << "ok " << endl;
  TString JetOrBulk[NumberTypeAnalysis]={"Jet", "Bulk", "All", "AllNS", "AwaySide", "AllButJet", "JetBFit", "JetZYAM"};
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

  Float_t SpectrumSup[numtipo][NumberTypeAnalysis]={{0.015,0.2, 0.2,0.8,0.2,0.2,0.015, 0.015},{0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001}};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};

  TString SSysT[numSysTrigger]={"DCAz < 1","DCAz < 2","DCAz < 0.5"};
  TString SSysV0[numSysV0]={"default", "cosTP> 0.997", "ctau <3 ", "YK0s < 0.5", "Lrejection 10 MeV","no Lrejection 5 MeV", "V0dca< 0.3"}; //all except no Lrejecxtion 5 MeV are done with the default cut on (i.e. Lrejection 5 MeV)
  TString SSysV0hh[numSysV0hh]={"DCAz < 1","DCAz < 2","DCAz < 0.5"}; 
  TString SSyst[numSyst]={"default", "cosTP> 0.997", "ctau <3 ", "YK0s < 0.5", "Lrejection 10 MeV","no Lrejection 5 MeV", "V0dca< 0.3", "Sideband 5sigma", "central 4sigma", "", "jet0.5", "bulk",  "BC 1.2"};
  TString SSysthh[numSysthh]={"DCAz < 1","DCAz < 2","DCAz < 0.5", " ", " ", " ",  "jet0.5", "bulk",  "BC 1.2"};
  TString SSystBetter[numSyst]={"default","cos(#theta_{P})", "c#tau", "YK0s", "#Lambda-rejection","", "V0 DCA_{PV}", "Sidebands region", "Peak region", "", "Jet region", "Out-of-jet region",  "Bin counting range"};
  TString SSystBetterhh[numSysthh]={"default","DCAz < 2","DCAz < 0.5", " ", " ", " ", "Jet region", "Out-of-jet region",  "Bin counting range"};

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
  TString SFit[3]={"exponential", "power law", "retta"};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};
  Int_t inutile[numSyst]=     {0,1,    1, 0,  1,  0,  1,  1,  1,  0,  1,  1,  0};
  Int_t Marker[numSyst]={20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  Int_t MarkerBetter[numSyst]={1,22, 32, 1, 29, 1,   3,  34, 33, 1, 20, 21, 22};
  //  Int_t Color[numSyst]= {1,  2,  3,  4,  5,  6,  7,  7,  4, 10,  6,  1,  2};
  Int_t Color[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 628,868,1, 909, 1,801,418,860,  1, 881, 7, 1};
  Int_t ColorBetter[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};
  cout << "ok " << endl;
  TCanvas *canvasPhiSys[4*(nummolt+1)];
  TCanvas *canvasSpectrumSys[2*(nummolt+1)];
  TCanvas *canvasSpectrumSysBis[2*(nummolt+1)];
  TLegend *legend_Bpassed[nummolt+1][numPtV0];
  TLegend *legend_Bpassed_Spectrum[nummolt+1];      
  TLegend *legend_corr[nummolt+1][numPtV0];      
  TLegend *legend_Corr_Spectrum[nummolt+1];      
  TLegend *legend_UnCorr[nummolt+1][numPtV0];      
  TLegend *legendCorrBis[nummolt+1];      
  TLegend *legendErrorSpectrum;
  TLegend *legendBulk;
  auto legend3= new TLegend(0.55, 0.55, 0.9, 0.9);
  legend3->SetHeader("Selezioni applicate");     
  // auto legend4= new TLegend(0.55, 0.55, 0.9, 0.9); //ter
  auto legend4= new TLegend(0.55, 0.55, 0.9, 0.9); //noter
  //  legend4->SetHeader("Selezioni applicate");     
  
  auto legendYield= new TLegend(0.7, 0.1, 0.9, 0.3);
  legendYield->SetHeader("Fit functions");     
  auto legendYieldError= new TLegend(0.7, 0.1, 0.9, 0.3);
  legendYieldError->SetHeader("Type of error");     

  Int_t sys=0;
  TH1D* fHistSpectrum_master[nummolt+1];
  TH1D* fHistSpectrum_masterOnlyStat[nummolt+1];
  TH1D* fHistSpectrum_masterOnlyStatRel[nummolt+1];
  TH1D* fHistSpectrum_masterSystCorr[nummolt+1];
  TH1D* fHistSpectrum_masterSystUnCorr[nummolt+1];
  TH1D* fHistSpectrum_masterTotalUnCorr[nummolt+1];
  TH1D* fHistSpectrum_masterTotal[nummolt+1];
  TH1D* fHistSpectrum[nummolt+1][numSyst];  
  TH1D* fHistSpectrumPart[nummolt+1][numSyst];  
  TH1D* fHistSpectrum_Corr[nummolt+1][numSyst];
  TH1D* fHistSpectrum_ratio[nummolt+1][numSyst];
  TH1D* fHistSpectrum_Barlow[nummolt+1][numSyst];
  TH1D* fHistSigmaBarlowSpectrum[nummolt+1][numSyst];
  TH1D* fHistSpectrumErrorSist[nummolt+1][numSyst];
  TH1D* fHistSpectrumErrorSistUnCorr[nummolt+1][numSyst][10];
  TH1D* fHistSpectrumErrorSistCorr[nummolt+1][numSyst][10];
  TH1D* fHistSpectrumErrorSistCorrSI[nummolt+1][10];
  TH1D* fHistSpectrumErrorSistUnCorrSI[nummolt+1][50];
  TH1D* fHistSpectrumErrorSistPhi[nummolt+1];
  TH1D* fHistSpectrumErrorSistBC[nummolt+1];

  TH1D* fHistPhiDistr[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_solostat[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_master[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_ratio[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_Barlow[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_Corr[nummolt+1][numPtV0][numsystPhi + numSysV0];  

  TF1*  gauss[nummolt+1][numPtV0][numSyst];
  TF1*  gaussint[nummolt+1][numPtV0][numSyst];

  Double_t YieldPerErrore[nummolt+1][numSyst]={0};
  Double_t YieldDefault[nummolt+1][numSyst]= {0};

  Double_t Sigma[nummolt+1][numPtV0][numSyst]={0};
  Double_t SigmaSpectrum[nummolt+1][numSyst]={0};
  Double_t SigmaSystSpectrumUnCorr[nummolt+1][numPtV0]={0};
  Double_t SigmaSystSpectrumUnCorrBC[nummolt+1][numPtV0]={0};
  Double_t SigmaSystSpectrumCorrBC[nummolt+1][numSyst]={0};
  Double_t SigmaSystSpectrumCorr[nummolt+1][numSyst]={0};
  Double_t SigmaSystSpectrumCorrOK[nummolt+1]={0}; //effetto correlato in pt è propagato correttamente su yield vs mult
  Double_t SigmaTotalSpectrumUnCorr[nummolt+1][numSyst]={0};
  Double_t SigmaTotalSpectrum[nummolt+1][numSyst]={0};
  Double_t Mean[nummolt+1][numPtV0][numSyst]={0}; 
  Double_t MeanSpectrum[nummolt+1][numSyst]={0};
  Bool_t   Correlated[nummolt+1][numPtV0][numSyst]={0};  
  Bool_t   CorrelatedSpectrum[nummolt+1][numSyst]={0};  
  Bool_t   CorrelatedBis[nummolt+1][numSyst]={0};
  Int_t    CorrelatedBisSys[nummolt+1]={0};
  Bool_t   BarlowPassed[nummolt+1][numPtV0][numSyst]={0};
  Bool_t   BarlowPassedSpectrum[nummolt+1][numSyst]={0};
  Int_t    NumberBarlowPassed[nummolt+1][numPtV0]={0};
  Int_t    NumberBarlowPassedSpectrum[nummolt+1]={0};
  Int_t    NumberCorr[nummolt+1][numPtV0]={0};  Int_t    NumberCorrSpectrum[nummolt+1]={0};
  Double_t SigmaBarlow[nummolt+1][numPtV0][numSyst][50]={0};
  Double_t SigmaBarlowMax[numSyst]={0};
  Double_t SigmaBarlowMin[numSyst]={1000000};
  Double_t SigmaBarlowAverage[numSyst]={0};

  Int_t   NTrigger[nummolt+1]={0}; //total number of trigger particles 
  Int_t   NTriggerV0[nummolt+1]={0}; //total number of trigger particles only in events with V > 0
  Float_t NSpectrum[nummolt+1][numPtV0][numSyst]={0}; 

  Float_t NSpectrumFinal[nummolt+1][numPtV0][100]={0}; 
  Float_t NSpectrumErrorFinal[nummolt+1][numPtV0][100]={0}; 
  Float_t NSpectrumErrorSistUnCorrPiuStat[nummolt+1][numPtV0][numSyst]={0}; 

  Float_t NSpectrumErrorSistUnCorr[nummolt+1][numPtV0][numSyst][50]={0}; //errore sistematico scorrelato associato a m, pt e eff.sist in delta phi
  Float_t NSpectrumErrorSistUnCorrSI[nummolt+1][numPtV0][50]={0}; //errore sistematico scorrelato associato a m, pt 
  Float_t NSpectrumErrorSistSI[nummolt+1][numPtV0][50]={0}; 
  Float_t NSpectrumErrorSistCorr[nummolt+1][numPtV0][numSyst][50]={0}; //errore sistematico correlato associato a m, pt e eff.sist in delta phi
  Float_t NSpectrumErrorSistCorrSI[nummolt+1][numPtV0][50]={0}; //errore sistematico correlato associato a m, pt 
  Float_t NSpectrumErrorSoloStat[nummolt+1][numPtV0][numSyst]={0}; 
  Float_t NSpectrumErrorSoloStatFinal[nummolt+1][numPtV0][100]={0};

  TH1D* fHistSpectrumBulk_master[nummolt+1];

  TH1D* fHistSigmaSyst[nummolt+1][numPtV0][numSyst];
  TH1D* fHistSigmaSystNSmoothed[nummolt+1][numPtV0][numSyst];
  Double_t SigmaSyst[nummolt+1][numPtV0][50]={0};
  
  TH1F*   fHistV0EfficiencyPtBins[nummolt];
  TH1F*   HistContV0PtBins[nummolt];

  TH1F * fHistBinCountingRegionAndPtFit = new TH1F ("fHistBinCountingRegionAndPtFit", "fHistBinCountingRegionAndPtFit", 15, 0, 15);
  fHistBinCountingRegionAndPtFit->GetXaxis()->SetBinLabel(1,"BCLow");
  fHistBinCountingRegionAndPtFit->GetXaxis()->SetBinLabel(2,"BCUp");
  for(Int_t m=0; m<nummolt+1; m++){
    fHistBinCountingRegionAndPtFit->GetXaxis()->SetBinLabel(3+2*m,Form("FitLow_%i", m));
    fHistBinCountingRegionAndPtFit->GetXaxis()->SetBinLabel(2+2*(m+1),Form("FitUp_%i", m));
  }

  cout << "ok " << endl;
  TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Lambda", "Xi","Xi", "Omega", "Omega", "Xi", "Omega"};
  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask"+dirinputtype[type]);
  if (!dir) {cout << " input dir not found " << endl; return;}
  TList *list = (TList*)dir->Get("MyOutputContainer");
  if (!list) {cout << " input list not found " << endl; return;}
  cout << "ok 1" << endl;
  TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtMaxvsMultBefAll");
  if (MasterThesisAnalysis)  fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtvsMultBefAll");
  cout << "ok 2" << endl;
  TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMin+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMax-0.00001) );
  cout << "ok 3" << endl;
  TString stringout;
  cout << "ok 8" << endl;

  stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
  //  stringout +=Form("_PtBinning%i", PtBinning);
  stringout +=Path1; //it was without year
  if(type>=0){
    stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+= hhCorr[ishhCorr]+"_" + JetOrBulk[TypeAnalysis]+DataOrMC[isMC] + Form("_PtMin%.1f.root", PtTrigMin);

  if (isEnlargedDeltaEtaPhi)  stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +Path1+ hhCorr[ishhCorr]+"_" + JetOrBulk[TypeAnalysis]+DataOrMC[isMC] + Form("_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);

  if(MasterThesisAnalysis){
    stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +hhCorr[ishhCorr]+Form("_Jet.root", PtTrigMin);
    if (isMC && isEfficiency) stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +hhCorr[ishhCorr]+Form("_JetMC.root", PtTrigMin);
    if (isBulk) stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +hhCorr[ishhCorr]+Form("_Bulk.root", PtTrigMin);
    if (isBulk && isMC &&isEfficiency) stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +hhCorr[ishhCorr]+Form("_BulkMC.root", PtTrigMin);
    if (isTotal) stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +hhCorr[ishhCorr]+Form("_All.root", PtTrigMin);
    if (isTotal && isMC &&isEfficiency) stringout = Dir+"/DATA"+year0+"/SystematicAnalysis" +hhCorr[ishhCorr]+Form("_AllMC.root", PtTrigMin);
  }  

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

  if (TypeAnalysis==0){
    if (!ishhCorr && type==0){// || type==8)){
      ALowBinFit[0]=	ALowBin[0]={-1};
      AUpBinFit[0]=	AUpBin[0]={1};
      
      ALowBinFit[1]=	ALowBin[1]={-1.2};
      AUpBinFit[1]=	AUpBin[1]={1.2};
    }
    else if (!ishhCorr && (type==4 || type==5 || type==8)){
      ALowBinFit[0]=	ALowBin[0]={-0.7};
      AUpBinFit[0]=	AUpBin[0]={0.7};
      
      ALowBinFit[1]=	ALowBin[1]={-1.1};
      AUpBinFit[1]=	AUpBin[1]={1.1};
    }
    if (!ishhCorr && isEnlargedDeltaEtaPhi){
      ALowBinFit[0]=	ALowBin[0]={-1.2};
      AUpBinFit[0]=	AUpBin[0]={1.2};
      
      ALowBinFit[1]=	ALowBin[1]={-1.4};
      AUpBinFit[1]=	AUpBin[1]={1.4};
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
  else if (TypeAnalysis==6){
    ALowBinFit[0]=	ALowBin[0]=-1;
    AUpBinFit[0]=	AUpBin[0]=1;
  }
  else if (TypeAnalysis==7){
    ALowBinFit[0]=	ALowBin[0]=-1;
    AUpBinFit[0]=	AUpBin[0]=1;
  }


  //low and upper ranges for the fit****************************
  //************************************************************

  Double_t LowRange[nummolt+1]= {2, 2, 1, 1.5, 1.5, 1.5}; 
  Double_t UpRange[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeSpectrumPart[nummolt+1]= {0};

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
	  LowRangeBulk[m]= 0; 
	  LowRangeAll[m]= 0; 
	  UpRangeAll[m]= 2; 
	  UpRangeBulk[m]= 2; 
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

  //end of low and upper ranges for the fit****************************
  //************************************************************

  Int_t PtBinMin[nummolt+1]={0};
  for (Int_t m=0; m<nummolt+1; m++){
    if (TypeAnalysis==0){
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

    //bin counting is done at pt> LowRangeSpectrumPart, integral of fit function at pt <  LowRangeSpectrumPart
    if  (TypeAnalysis==1 || TypeAnalysis==2){
      if (type>0)    LowRangeSpectrumPart[m] = 0.5;
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


  fHistBinCountingRegionAndPtFit->SetBinContent(1,ALowBinFit[0]);
  fHistBinCountingRegionAndPtFit->SetBinContent(2,AUpBinFit[0]);
  for (Int_t m=0; m<nummolt+1; m++){
    fHistBinCountingRegionAndPtFit->SetBinContent(3+m*2,LowRange[m]);
    fHistBinCountingRegionAndPtFit->SetBinContent(2+(m+1)*2,UpRange[m]);
  }

  //syst1_limit indica gli effetti sistematici sullo spettro in pT
  Int_t syst1_limit=0;
  if(!isBulk) syst1_limit= 1; //2 for BC systematic
  if (isBulk) syst1_limit= 1;
  if (isTotal) syst1_limit= 1; //non considero sistematico associato a bin counting per total
  if (TypeAnalysis>2) syst1_limit=1;

  //********************************************************************* 
  //**************calcolo numero particelle di trigger*******************
  //********************************************************************* 
  for(Int_t m=0; m<nummolt+1; m++){

    for(Int_t syst=0; syst<numSystGlobal; syst++){
      fHistSpectrum[m][syst]=new TH1D ("fHistSpectrum_"+Smolt[m]+Form("_Sys%i",syst),"fHistSpectrum_"+Smolt[m]+Form("_Sys%i",syst), numPtV0Max, NPtV0) ;     

    }

    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.001); j<=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<=fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    cout << "n trigger in mult range (all triggers)    " << m << "  " <<  NTrigger[m] <<   endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////

  Int_t sysV0=0;    
  if(isMC && isEfficiency) file = year + "_MCEff";
  //  file+=Form("_PtBinning%i", PtBinning);
  //file+= Path1;
 
  for(Int_t m=0; m<nummolt+1; m++){
    //    if (m==0) continue;
    //   if (m>=3) continue;     //added  
    //cout << "\n\n****************************************************\n" <<m << endl;
    //    cout << "\n****************************************************\n" << endl;
   
    for(Int_t syst=0; syst<numSysV0Global+numsystPhiGlobal; syst++){
      cout << "********************************* syst = " << syst << endl;
      if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
      if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
      if (syst< numSysV0Global){
	sysV0=syst;
	sys=0;
      }
      else {
	sys=syst-numSysV0Global+1;
	sysV0=0;
      }
      cout << "syst " << syst << " sysV0 " << sysV0 << " sys " << sys << endl;    
      if (ishhCorr && sys!=0 && (sys<4 || sys>5)) {
	cout << "questo valore di sys non è previsto per l'analisi hhCorr" << endl;
	continue;
      }

      if (isMC && !isEfficiency) PathIn1=Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file+  Form("_Sys%i", sys)+"_Output.root";
      else {
	PathIn1 = Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file;
	if(type>=0){
	  PathIn1 +="_"+tipo[type];
	  PathIn1 +=Srap[israp];
	  PathIn1 +=SSkipAssoc[SkipAssoc];
	}
	PathIn1+= hhCorr[ishhCorr]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sys, PtTrigMin)+"_Output.root";
	if (isEnlargedDeltaEtaPhi) PathIn1 = Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file  + hhCorr[ishhCorr]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sys, PtTrigMin)+"_DeltaEtaPhiEnlarged_Output.root";
	if (MasterThesisAnalysis)	PathIn1 = Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file  + hhCorr[ishhCorr]+Form("_SysT%i_SysV0%i_Sys%i", sysTrigger, sysV0, sys)+"_Output.root";
      }
      cout << "\n\n" << PathIn1 << endl;
      filein = new TFile(PathIn1, "");
      if (!filein) {cout << filein << "does not exist " << endl ; return; }
     
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	cout << "getting histograms....m: " << m << "  v: " << v << " syst: " << syst << endl;
	if(syst==0){
	  if (TypeAnalysis==6) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_EffCorr_BulkSubBulkFit");
	  else 	  if (TypeAnalysis==7) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_EffCorr_BulkSubZYAM");
	  else 	  if (!isBulk) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	  if (isBulk) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");
	  if (isTotal || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr");
	  if (TypeAnalysis==3) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled");
	  fHistPhiDistr_master[m][v]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]);
	}

	if (TypeAnalysis==6) fHistPhiDistr[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_EffCorr_BulkSubBulkFit");
	else 	  if (TypeAnalysis==7) fHistPhiDistr[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_EffCorr_BulkSubZYAM");

	else	if (!isBulk)	fHistPhiDistr[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	if (isBulk)	fHistPhiDistr[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");
	if (isTotal || TypeAnalysis==4 || TypeAnalysis==5)	fHistPhiDistr[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr");
	if (TypeAnalysis==3) fHistPhiDistr[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled");

	if (TypeAnalysis==6) fHistPhiDistr_solostat[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_EffCorr_BulkSubBulkFit");
	else 	  if (TypeAnalysis==7) fHistPhiDistr_solostat[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_EffCorr_BulkSubZYAM");

	else	if (!isBulk) fHistPhiDistr_solostat[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	if (isBulk) fHistPhiDistr_solostat[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");
	if (isTotal || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_solostat[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr");
	if (TypeAnalysis==3) fHistPhiDistr_solostat[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled");

	if(!fHistPhiDistr[m][v][syst]         ) return;
	if(!fHistPhiDistr_solostat[m][v][syst]) return;
	if(!fHistPhiDistr_master[m][v])         return;

	fHistPhiDistr[m][v][syst]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v] +Form("_syst%i", syst));
	fHistPhiDistr_solostat[m][v][syst]->SetName("PhiDistr_solostat_m"+Smolt[m]+"_v"+SPtV0[v] +Form("_syst%i", syst));
	fHistPhiDistr_ratio[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistPhiDistr_ratio_m%i_v%i_syst%i",m,v,syst));
	fHistPhiDistr_Barlow[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistPhiDistr_Barlow_m%i_v%i_syst%i",m,v,syst));

	if(syst==0){	
	  fHistPhiDistr_ratio[m][v][syst]={0};
	  fHistPhiDistr_Barlow[m][v][syst]={0};
	}

	if(syst!=0){	
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX(); j++){
	    // cout << "num bin" << fHistPhiDistr[m][v][syst]->GetNbinsX()-5 << endl;
	    if(fHistPhiDistr_master[m][v]->GetBinContent(j)!=0){
	      fHistPhiDistr_ratio[m][v][syst]->SetBinContent(j,fHistPhiDistr[m][v][syst]->GetBinContent(j)/fHistPhiDistr_master[m][v]->GetBinContent(j));
	      fHistPhiDistr_ratio[m][v][syst]->SetBinError(j,0);
	    }
	    if(sqrt(TMath::Abs(pow(fHistPhiDistr[m][v][syst]->GetBinError(j),2)-pow(fHistPhiDistr_master[m][v]->GetBinError(j),2)))!=0){
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinContent(j,(fHistPhiDistr[m][v][syst]->GetBinContent(j)-fHistPhiDistr_master[m][v]->GetBinContent(j))/sqrt(TMath::Abs(pow(fHistPhiDistr[m][v][syst]->GetBinError(j),2)-pow(fHistPhiDistr_master[m][v]->GetBinError(j),2))));
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinError(j,0);
	    }
	    else{
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinContent(j,0);
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinError(j,0);
	      cout << "mult " << m << " v " << v << " sys" << syst<< "j " << j << " denominatore Barlow = 0" << endl;
	    }
	  }
	}
	
	Int_t count=0;
	Correlated[m][v][0]=kTRUE; 
	BarlowPassed[m][v][0]=kTRUE;
	fHistPhiDistr_Corr[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistPhiDistr_Corr_m%i_v%i_syst%i",m,v,syst));
	if (syst!=0) {
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){ //escludo bin ai bordi
	    if (TMath::Abs(fHistPhiDistr_Barlow[m][v][syst]->GetBinContent(j)) >= 2) count++;
	  }
	  //	  if (count > 0.05* (fHistPhiDistr_Barlow[m][v][syst]->GetNbinsX()-5)) BarlowPassed[m][v][syst]=kTRUE;
	  if (count >= 6) BarlowPassed[m][v][syst]=kTRUE;
	  if ( BarlowPassed[m][v][syst]==kTRUE) {
	    //	    cout << "systematic effect n. " << SSyst[syst]<< " for m " << Smolt[m] << " and v " << SPtV0[v] << " has failed the Barlow check within 2 sigmas and is not significant" << endl;
	    NumberBarlowPassed[m][v]+=1;
	    for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){
	      Mean[m][v][syst] += fHistPhiDistr_ratio[m][v][syst]->GetBinContent(j);
	    }
	    Mean[m][v][syst]=Mean[m][v][syst]/fHistPhiDistr[m][v][syst]->GetNbinsX();
	    for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){
	      Sigma[m][v][syst] += pow(fHistPhiDistr_ratio[m][v][syst]->GetBinContent(j)-Mean[m][v][syst],2);
	    }
	    Sigma[m][v][syst]=sqrt(Sigma[m][v][syst]/(fHistPhiDistr[m][v][syst]->GetNbinsX()-5)/(fHistPhiDistr[m][v][syst]->GetNbinsX()-5+1));
	    if(TMath::Abs((Mean[m][v][syst]-1))>numSigmaCorr*Sigma[m][v][syst]) Correlated[m][v][syst]=kTRUE;
	    if (Correlated[m][v][syst]==kTRUE){
	      NumberCorr[m][v]+=1;
	    }
	  
	    Int_t numPhiBins = fHistPhiDistr[m][v][syst]->GetNbinsX();

	    for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){
	      SigmaBarlow[m][v][syst][j]=(fHistPhiDistr[m][v][syst]->GetBinContent(j)-fHistPhiDistr_master[m][v]->GetBinContent(j))/sqrt(12.);
	      //	if (TMath::Abs(SigmaBarlow[m][v][syst][j]) > 1000) cout << m << " " << v << " " << syst <<  " " << j<< " " <<fHistPhiDistr[m][v][syst]->GetBinContent(j)<< "  " << fHistPhiDistr_master[m][v]->GetBinContent(j) << "  " << SigmaBarlow[m][v][syst][j] << endl;
	    }
	    // i sistematici 4 e 5 sono tra loro correlati, prendo la massima variazione come errore
	    if (Correlated[m][v][syst]==kFALSE){
	      for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX(); j++){//errore pari a zero in tutti i bin
		fHistPhiDistr_Corr[m][v][syst]->SetBinContent(j,0);
	      }
	    }	      
	  }
	}
      } //end loop on pt v0
    }//end loop on syst
   
    //this part is used for those variables to which more than one selection has been applied (es. lambda rejection in hV0 analysis, DCAz in hh analysis)
    if (!ishhCorr){
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (BarlowPassed[m][v][4]==kTRUE && BarlowPassed[m][v][5]==kTRUE){
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][0]->GetNbinsX()-5; j++){
	    if (SigmaBarlow[m][v][4][j]>=SigmaBarlow[m][v][5][j]) SigmaBarlow[m][v][5][j]=0;
	    else {SigmaBarlow[m][v][4][j]=SigmaBarlow[m][v][5][j]; SigmaBarlow[m][v][5][j]=0;}
	  }
	}
	if (BarlowPassed[m][v][4]==kFALSE && BarlowPassed[m][v][5]==kTRUE){
	  BarlowPassed[m][v][4]=kTRUE;
	  if (Correlated[m][v][5]) Correlated[m][v][4]=kTRUE;
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][0]->GetNbinsX()-5; j++){
	    SigmaBarlow[m][v][4][j]=SigmaBarlow[m][v][5][j]; 
	    SigmaBarlow[m][v][5][j]=0;
	  }
	}
      }
    }

    if (ishhCorr){
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (BarlowPassed[m][v][1]==kTRUE && BarlowPassed[m][v][2]==kTRUE){
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][0]->GetNbinsX()-2; j++){
	    if (SigmaBarlow[m][v][1][j]>=SigmaBarlow[m][v][2][j]) SigmaBarlow[m][v][2][j]=0;
	    else {SigmaBarlow[m][v][1][j]=SigmaBarlow[m][v][2][j]; SigmaBarlow[m][v][2][j]=0;}
	  }
	}
	if (BarlowPassed[m][v][1]==kFALSE && BarlowPassed[m][v][2]==kTRUE){
	  BarlowPassed[m][v][1]=kTRUE;
	  if (Correlated[m][v][2]) Correlated[m][v][1]=kTRUE;
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][0]->GetNbinsX()-2; j++){
	    SigmaBarlow[m][v][1][j]=SigmaBarlow[m][v][2][j]; 
	    SigmaBarlow[m][v][2][j]=0;
	  }
	}
      }
    }
   
    Int_t numPhiBins = fHistPhiDistr[0][1][0]->GetNbinsX();

    //**************************************************************************************************************** 
    cout << "\n\ncalcolo l'errore sistematico (non correlato in deltaPhi) e lo sommo all'errore statistico in quadratura..." << endl;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      for(Int_t syst=0; syst<numSysV0Global+numsystPhiGlobal; syst++){
	cout << " open loop on syst " << endl;

	if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	cout << "  getting histograms....m: " << m << "  v: " << v << " syst: " << syst << endl;

	if (!ishhCorr) cout << " type of systematic error considered " << SSyst[syst] << endl; 
	if (ishhCorr) cout << " type of systematic error considered " << SSysthh[syst] << endl; 
	//	cout << "ishhCorr " <<  ishhCorr << endl;

	//cout << "ok" << endl;

     	fHistSigmaSyst[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistSigmaSyst_m%i_v%i_syst%i",m,v,syst));
	fHistSigmaSystNSmoothed[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistSigmaSyst_notsmoothed_m%i_v%i_syst%i",m,v,syst));

	for(Int_t j=1; j <= numPhiBins; j++){
	  fHistSigmaSyst[m][v][syst]->SetBinContent(j, 0);
	  fHistSigmaSyst[m][v][syst]->SetBinError(j, 0);
	  fHistSigmaSystNSmoothed[m][v][syst]->SetBinContent(j, 0);
	  fHistSigmaSystNSmoothed[m][v][syst]->SetBinError(j, 0);
	}
	
	for(Int_t j=1; j <= numPhiBins; j++){
	  //   cout << " m" << m << " v " << v << " syst " << syst << "j " << j << " stat error (relative, not abs. value): " << fHistPhiDistr[m][v][syst]->GetBinError(j)/fHistPhiDistr[m][v][syst]->GetBinContent(j) << " sist error (relative, abs. value): " << TMath::Abs(SigmaBarlow[m][v][syst][j])/fHistPhiDistr[m][v][syst]->GetBinContent(j) << endl; 
	  fHistSigmaSyst[m][v][syst]->SetBinContent(j, 0);
	  fHistSigmaSyst[m][v][syst]->SetBinError(j, 0);
	  //cout << " j" << j << " v " << v << " syst " << syst << endl; 
	  if (j < numPhiBins-5){
	    fHistSigmaSystNSmoothed[m][v][syst]->SetBinContent(j, TMath::Abs(SigmaBarlow[m][v][syst][j]/fHistPhiDistr_master[m][v]->GetBinContent(j)) );
	    //	    fHistSigmaSyst[m][v][syst]->SetBinContent(j, TMath::Abs(SigmaBarlow[m][v][syst][j]));
	    fHistSigmaSyst[m][v][syst]->SetBinContent(j, SigmaBarlow[m][v][syst][j]);
	    fHistSigmaSyst[m][v][syst]->SetBinError(j, 0);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetBinError(j, 0);
	  }
	}
	fHistSigmaSyst[m][v][syst]->Smooth();
	cout << " end loop on syst considered " << endl;
      }//end loop on syst
      cout << " end loop on all syst " << endl;

      for(Int_t j=1; j < numPhiBins-5; j++){
	for(Int_t syst=0; syst<numSysV0Global+numsystPhiGlobal; syst++){
	  if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	  if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
     	  if (Correlated[m][v][syst]==kTRUE)continue;
	  SigmaSyst[m][v][j]+= pow(fHistSigmaSyst[m][v][syst]->GetBinContent(j),2);
	}
      
	SigmaSyst[m][v][j]= sqrt(SigmaSyst[m][v][j]);
	fHistPhiDistr[m][v][0]->SetBinError(j,sqrt(pow(fHistPhiDistr_master[m][v]->GetBinError(j),2) +pow(SigmaSyst[m][v][j],2)));  
	fHistSigmaSystNSmoothed[m][v][0]->SetBinContent(j,sqrt(pow(fHistPhiDistr_master[m][v]->GetBinError(j),2) +pow(SigmaSyst[m][v][j],2))/TMath::Abs(fHistPhiDistr_master[m][v]->GetBinContent(j) ));
	//	fHistSigmaSystNSmoothed[m][v][0]->SetBinContent(j,sqrt(pow(fHistPhiDistr_master[m][v]->GetBinError(j),2) )/TMath::Abs(fHistPhiDistr_master[m][v]->GetBinContent(j) ));
	fHistSigmaSystNSmoothed[m][v][0]->SetBinError(j,0);
      }
    }//end loop on v

    for(Int_t syst=0; syst< numSysV0Global + numsystPhiGlobal; syst++){    
      Int_t count_corr=0;
      for(Int_t v=PtV0Min; v< numPtV0Max; v++){    
	if(Correlated[m][v][syst]==kFALSE) continue; 
	count_corr++;
	//	cout << count_corr << endl;
      }
      if (count_corr==numPtV0Max) {
	CorrelatedBis[m][syst]=kTRUE;
	CorrelatedBisSys[m]++;
      }
    }
    
    fHistSpectrum_masterOnlyStat[m]=(TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_"+Smolt[m]+"_eta<0.5");
    cout << "\n\ncalcolo l'errore sistematico (correlato in deltaPhi)" << endl;
    //    cout << numSystGlobal- numsystPhiGlobal-numSysV0Global+1 << endl;
    for (Int_t syst1=0; syst1<syst1_limit; syst1++){ //to consider systematics which are only related to bin counting
      //      cout << " syst1 " << syst1 << endl;
      for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){       
	fHistSpectrumErrorSistUnCorr[m][syst][syst1]=(TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrumSistUnCorr_"+Smolt[m]+Form("_syst%i_syst1%i", syst, syst1));
	fHistSpectrumErrorSistCorr[m][syst][syst1]=(TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrumSistCorr_"+Smolt[m]+Form("_syst%i_syst1%i", syst, syst1));
      }
      fHistSpectrumErrorSistUnCorrSI[m][syst1]=(TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrumSistUnCorr_"+Smolt[m]+Form("_syst1%i", syst1));
      fHistSpectrumErrorSistCorrSI[m][syst1]=(TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrumSistCorr_"+Smolt[m]+Form("_syst1%i", syst1));

      //      cout << "ho def istogrammi" << endl;
      for(Int_t v=PtV0Min; v< numPtV0Max; v++){ 
	fHistSpectrumErrorSistUnCorrSI[m][syst1]->SetBinContent(v+1,0);
	fHistSpectrumErrorSistCorrSI[m][syst1]->SetBinContent(v+1,0);
	fHistSpectrumErrorSistUnCorrSI[m][syst1]->SetBinError(v+1,0);
	fHistSpectrumErrorSistCorrSI[m][syst1]->SetBinError(v+1,0);
	//	cout << "syst " << syst1 << " v " <<  v<<  endl;	
	for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){       
	  fHistSpectrumErrorSistUnCorr[m][syst][syst1]->SetBinContent(v+1,0);
	  fHistSpectrumErrorSistCorr[m][syst][syst1]->SetBinContent(v+1,0);
	  fHistSpectrumErrorSistUnCorr[m][syst][syst1]->SetBinError(v+1,0);
	  fHistSpectrumErrorSistCorr[m][syst][syst1]->SetBinError(v+1,0);
	  NSpectrum[m][v][syst]=0;      
	  NSpectrumErrorSoloStat[m][v][syst]=0;
	  //  cout << "syst " << syst << "  " << 	NSpectrum[m][v][0]<<endl;//"  " <<	NSpectrumErrorSoloStat[m][v][syst]<<  endl;	
	}
	
	//	 cout <<  	NSpectrum[m][v][0]<<"  " <<	NSpectrumErrorSoloStat[m][v][0]<<  endl;	
	NSpectrumErrorSistCorrSI[m][v][syst1]=0;
	
	for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	  // cout << "gli indici vanno da 0 a " <<numsystPhiGlobal + numSysV0Global-1 << " compresi" << endl;
	  // cout << "\n\n\n ************************************syst "<<syst << endl;
	  if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	  if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  //inserisco errori sist correlati in delta phi nello spettro relativo a syst =0	  
	  if ((Correlated[m][v][syst]==kTRUE && syst< 12) ){
	    //	    cout << "m " << m << " v" << v << " sono un sist correlato in delta phi" << endl;
	    for(Int_t j=fHistPhiDistr_Corr[m][v][syst]->FindBin(ALowBinFit[syst1]); j<= fHistPhiDistr_Corr[m][v][syst]->FindBin(AUpBinFit[syst1]); j++){
	      //	      if (j== fHistPhiDistr_Corr[m][v][syst]->FindBin(ALowBinFit[syst1])) cout << "\n\n\nbin content " << fHistPhiDistr_Corr[m][v][0]->GetBinContent(j)<< endl;
	      //	      if (j== fHistPhiDistr_Corr[m][v][syst]->FindBin(AUpBinFit[syst1])) cout << "\n\n\nbin content " << fHistPhiDistr_Corr[m][v][0]->GetBinContent(j)<< endl;
	      if (syst ==0) {
		NSpectrumErrorSoloStat[m][v][syst]+=pow(fHistPhiDistr_master[m][v]->GetBinError(j),2);
	      }
	      NSpectrum[m][v][syst]+=(fHistPhiDistr_Corr[m][v][0]->GetBinContent(j)+ 	    fHistSigmaSyst[m][v][syst]->GetBinContent(j));
	      //	   if (syst==0)	cout << "  v " << v << "syst " << syst << "  " << 	NSpectrum[m][v][0]<< endl;

	    }
	    //	    cout << " nspectrumerrossistuncorr per syst " << 	syst << "  " <<      NSpectrumErrorSistUnCorr[m][v][syst][syst1]<< endl;
       
	    NSpectrumErrorSistCorr[m][v][syst][syst1]=TMath::Abs(NSpectrum[m][v][0] -NSpectrum[m][v][syst]);
	    fHistSpectrumErrorSistCorr[m][syst][syst1]->SetBinContent(v+1,TMath::Abs(NSpectrum[m][v][0] -NSpectrum[m][v][syst]));
	    //	    cout << "get binc " << fHistSpectrumErrorSistCorr[m][syst][syst1]->GetBinContent(v+1);

	  }

	  if ((Correlated[m][v][syst]==kFALSE &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){
	    for(Int_t j=fHistPhiDistr_Corr[m][v][0]->FindBin(ALowBinFit[syst1]); j<= fHistPhiDistr_Corr[m][v][0]->FindBin(AUpBinFit[syst1]); j++){
	      NSpectrumErrorSistUnCorr[m][v][syst][syst1]+=pow(fHistSigmaSyst[m][v][syst]->GetBinContent(j), 2);
	    }
	    NSpectrumErrorSistUnCorr[m][v][syst][syst1]=sqrt(NSpectrumErrorSistUnCorr[m][v][syst][syst1]);
	    fHistSpectrumErrorSistUnCorr[m][syst][syst1]->SetBinContent(v+1,NSpectrumErrorSistUnCorr[m][v][syst][syst1]);
	  }
	} //end loop syst

	/*
	  cout << "\n\ntry -1 per pT V0 = " << v << " e molt " << m << endl;
	  for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	  if (syst==9 || syst == avoidthissyst || syst == avoidthissystbis || syst==avoidthissysttris) continue;
	  cout << "syst n. " << syst << endl;	
	  if ((Correlated[m][v][syst]==kFALSE &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){	  
	  cout << "errore sist scorrelato in delta phi " << NSpectrumErrorSistUnCorr[m][v][syst][syst1]<< endl;
	  }
	  if ((Correlated[m][v][syst] &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){
	  cout << "errore sist correlato in delta phi " << NSpectrumErrorSistCorr[m][v][syst][syst1]<< endl;
	  }
	  }
	*/

      }//end loop v

      // cout << "ora faccio smooth errori sistematici dello spettro in pt " << endl;
      //***effettuo smooth errori sistematici correlati e non correlati*******+
      for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){       
	fHistSpectrumErrorSistCorr[m][syst][syst1]->Smooth();
	fHistSpectrumErrorSistUnCorr[m][syst][syst1]->Smooth();
      }
      //      cout << "ho fatto smooth " << endl;
      for(Int_t v=PtV0Min; v< numPtV0Max; v++){ 
	for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	  if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	  if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
    
	  if ((Correlated[m][v][syst]==kTRUE && syst< 12) ){
	    NSpectrumErrorSistCorr[m][v][syst][syst1]= fHistSpectrumErrorSistCorr[m][syst][syst1]->GetBinContent(v+1);
	    NSpectrumErrorSistCorrSI[m][v][syst1]+= pow(fHistSpectrumErrorSistCorr[m][syst][syst1]->GetBinContent(v+1),2); 
	  }
	  if ((Correlated[m][v][syst]==kFALSE &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){
	    NSpectrumErrorSistUnCorr[m][v][syst][syst1]= fHistSpectrumErrorSistUnCorr[m][syst][syst1]->GetBinContent(v+1);
	    NSpectrumErrorSistUnCorrSI[m][v][syst1]+= pow(fHistSpectrumErrorSistUnCorr[m][syst][syst1]->GetBinContent(v+1),2); 
	  }
	}
	/*
	  cout << "\n\ntry 0   per pT V0 = " << v << " e molt " << m << endl;
	  cout << "alcuni valori possono essere zero come conseguenza dello smooth " << endl;
	  for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	  if (syst==9 || syst == avoidthissyst || syst == avoidthissystbis || syst==avoidthissysttris) continue;
	  cout << "syst n. " << syst << endl;	
	  if ((Correlated[m][v][syst]==kFALSE &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){	  
	  cout << "errore sist scorrelato in delta phi " << NSpectrumErrorSistUnCorr[m][v][syst][syst1]<< endl;
	  }
	  if ((Correlated[m][v][syst] &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){
	  cout << "errore sist correlato in delta phi " << NSpectrumErrorSistCorr[m][v][syst][syst1]<< endl;
	  }
	  }
	*/


	NSpectrumFinal[m][v][syst1]= NSpectrum[m][v][0];
	NSpectrumErrorSoloStatFinal[m][v][syst1]=NSpectrumErrorSoloStat[m][v][0];
	NSpectrumErrorFinal[m][v][syst1]= sqrt(NSpectrumErrorSoloStatFinal[m][v][syst1] +NSpectrumErrorSistCorrSI[m][v][syst1] + NSpectrumErrorSistUnCorrSI[m][v][syst1]);
	NSpectrumErrorSoloStatFinal[m][v][syst1]=sqrt(NSpectrumErrorSoloStat[m][v][0]);
	NSpectrumErrorSistCorrSI[m][v][syst1]=sqrt( NSpectrumErrorSistCorrSI[m][v][syst1]);
	NSpectrumErrorSistUnCorrSI[m][v][syst1]=sqrt( NSpectrumErrorSistUnCorrSI[m][v][syst1]);
	fHistSpectrumErrorSistCorrSI[m][syst1]->SetBinContent(v+1,NSpectrumErrorSistCorrSI[m][v][syst1]); 
	fHistSpectrumErrorSistUnCorrSI[m][syst1]->SetBinContent(v+1,NSpectrumErrorSistUnCorrSI[m][v][syst1]); 
	//***************************************************************

	for(Int_t i =0; i <2; i++){
	  DeltaPhiWidth[i]=TMath::Abs(fHistPhiDistr_Corr[m][v][0]->GetXaxis()->GetBinLowEdge(fHistPhiDistr_Corr[m][v][0]->FindBin(ALowBinFit[syst1])) - fHistPhiDistr_Corr[m][v][0]->GetXaxis()->GetBinUpEdge(fHistPhiDistr_Corr[m][v][0]->FindBin(AUpBinFit[syst1])));
	  DeltaPhiWidthApprox[i]=TMath::Abs(AUpBinFit[i]-ALowBinFit[i]);
	  if (TypeAnalysis==0)     	  DeltaPhiWidth[i]=TMath::Abs(fHistPhiDistr_Corr[m][v][0]->GetXaxis()->GetBinLowEdge(fHistPhiDistr_Corr[m][v][0]->FindBin(ALowBinFit[0])) - fHistPhiDistr_Corr[m][v][0]->GetXaxis()->GetBinUpEdge(fHistPhiDistr_Corr[m][v][0]->FindBin(AUpBinFit[0])));
	}

	
	cout << "\n\ntry 1 " <<  "per m " << m << " v " << v << " syst1 " << syst1<<endl;
	cout << "errore statistico qui calcolato" <<     NSpectrumErrorSoloStatFinal[m][v][syst1]<< endl;
	cout << "errore sist scorrelato in delta phi " << NSpectrumErrorSistUnCorrSI[m][v][syst1]<< endl;
	cout << "errore sist correlato in delta phi " << NSpectrumErrorSistCorrSI[m][v][syst1]<< endl; 
	cout << "errore sist totale(da phi) + stat" << NSpectrumErrorFinal[m][v][syst1]<< endl; 

	
	NSpectrumErrorSistCorrSI[m][v][syst1]=NSpectrumErrorSistCorrSI[m][v][syst1]/NTrigger[m];
	NSpectrumFinal[m][v][syst1]= NSpectrumFinal[m][v][syst1]/NTrigger[m];
	NSpectrumErrorFinal[m][v][syst1]= NSpectrumErrorFinal[m][v][syst1]/NTrigger[m];
	NSpectrumErrorSistUnCorrSI[m][v][syst1]= NSpectrumErrorSistUnCorrSI[m][v][syst1]/NTrigger[m];
	NSpectrumErrorSoloStatFinal[m][v][syst1]=	NSpectrumErrorSoloStatFinal[m][v][syst1]/NTrigger[m];
	
	
	//**********************************************************************************
	
	for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	  if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	  if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  NSpectrumErrorSistUnCorr[m][v][syst][syst1]=NSpectrumErrorSistUnCorr[m][v][syst][syst1]/NTrigger[m];
	  NSpectrumErrorSistCorr[m][v][syst][syst1]=NSpectrumErrorSistCorr[m][v][syst][syst1]/NTrigger[m];
	  NSpectrumErrorSistUnCorr[m][v][syst][syst1]=NSpectrumErrorSistUnCorr[m][v][syst][syst1]/fHistSpectrum[m][0]->GetBinWidth(v+1);
	  NSpectrumErrorSistCorr[m][v][syst][syst1]=NSpectrumErrorSistCorr[m][v][syst][syst1]/fHistSpectrum[m][0]->GetBinWidth(v+1);;
	  if (TypeAnalysis !=3){
	    NSpectrumErrorSistUnCorr[m][v][syst][syst1]=NSpectrumErrorSistUnCorr[m][v][syst][syst1]/DeltaPhiWidth[syst1];
	    NSpectrumErrorSistCorr[m][v][syst][syst1]= NSpectrumErrorSistCorr[m][v][syst][syst1]/DeltaPhiWidth[syst1];
	  }
	   
	}
	
	/*
	  cout << "\n\ntry intermediate " << endl;
	  for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	  if (syst==9 || syst == avoidthissyst || syst == avoidthissystbis || syst==avoidthissysttris) continue;
	  cout << "syst n. " << syst << endl;	
	  if ((Correlated[m][v][syst]==kFALSE &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){	  
	  cout << "errore sist scorrelato in delta phi " << NSpectrumErrorSistUnCorr[m][v][syst][syst1]<< endl;
	  }
	  if ((Correlated[m][v][syst] &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){
	  cout << "errore sist correlato in delta phi " << NSpectrumErrorSistCorr[m][v][syst][syst1]<< endl;
	  }
	  }
	*/

	//**********************************************************************************
	if (fHistSpectrum[m][0]->GetBinWidth(v+1)!=0){ //added
	  NSpectrumFinal[m][v][syst1]= NSpectrumFinal[m][v][syst1]/fHistSpectrum[m][0]->GetBinWidth(v+1);
	  NSpectrumErrorFinal[m][v][syst1]= NSpectrumErrorFinal[m][v][syst1]/fHistSpectrum[m][0]->GetBinWidth(v+1);
	  NSpectrumErrorSistUnCorrSI[m][v][syst1]= 	NSpectrumErrorSistUnCorrSI[m][v][syst1]/fHistSpectrum[m][0]->GetBinWidth(v+1);
	  NSpectrumErrorSoloStatFinal[m][v][syst1]= 	NSpectrumErrorSoloStatFinal[m][v][syst1]/fHistSpectrum[m][0]->GetBinWidth(v+1);
	  NSpectrumErrorSistCorrSI[m][v][syst1]=NSpectrumErrorSistCorrSI[m][v][syst1]/fHistSpectrum[m][0]->GetBinWidth(v+1);
	} //added
	if (TypeAnalysis !=3){
	  NSpectrumFinal[m][v][syst1]=  NSpectrumFinal[m][v][syst1]/DeltaPhiWidth[syst1];
	  NSpectrumErrorFinal[m][v][syst1]=  NSpectrumErrorFinal[m][v][syst1]/DeltaPhiWidth[syst1];
	  NSpectrumErrorSistUnCorrSI[m][v][syst1]= 	NSpectrumErrorSistUnCorrSI[m][v][syst1]/DeltaPhiWidth[syst1];
	  NSpectrumErrorSoloStatFinal[m][v][syst1]= 	NSpectrumErrorSoloStatFinal[m][v][syst1]/DeltaPhiWidth[syst1];
	  NSpectrumErrorSistCorrSI[m][v][syst1]=NSpectrumErrorSistCorrSI[m][v][syst1]/DeltaPhiWidth[syst1];
	}
	//faccio spettri in Pt per i diversi sistematici che hanno un effetto sullo spettro in pt
	cout << syst1 << endl;
	fHistSpectrum[m][syst1]->SetBinContent(v+1, NSpectrumFinal[m][v][syst1]);
	fHistSpectrum[m][syst1]->SetBinError(v+1, NSpectrumErrorFinal[m][v][syst1]);
	if(syst1==0){
	  fHistSpectrum_masterOnlyStat[m]->SetBinContent(v+1, NSpectrumFinal[m][v][syst1]);
	  fHistSpectrum_masterOnlyStat[m]->SetBinError(v+1, NSpectrumErrorSoloStatFinal[m][v][syst1]);
	}

	cout << "\nvalore spettro " << NSpectrumFinal[m][v][syst1] << endl;
	cout << "valore spettro da isto " << 	fHistSpectrum_masterOnlyStat[m]->GetBinContent(v+1)<< endl;
	cout << "m " << m << " v " << v << " syst1 " << syst1<<endl;
	cout << "errore statistico " <<    fHistSpectrum_masterOnlyStat[m]->GetBinError(v+1)/fHistSpectrum_masterOnlyStat[m]->GetBinContent(v+1)<< endl;
	cout << "errore statistico qui calcolato" <<     NSpectrumErrorSoloStatFinal[m][v][syst1]/ fHistSpectrum_masterOnlyStat[m]->GetBinContent(v+1)<< endl;
	cout << "errore sist scorrelato in delta phi " << NSpectrumErrorSistUnCorrSI[m][v][syst1]/ fHistSpectrum_masterOnlyStat[m]->GetBinContent(v+1)<< endl;
	cout << "errore sist correlato in delta phi " << NSpectrumErrorSistCorrSI[m][v][syst1]/ fHistSpectrum_masterOnlyStat[m]->GetBinContent(v+1)<< endl; 
	cout << "errore sist totale(da phi) + stat" << NSpectrumErrorFinal[m][v][syst1]/ fHistSpectrum_masterOnlyStat[m]->GetBinContent(v+1)<< endl; 
	for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	  if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	  if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	  cout << "syst n. " << syst << endl;	
	  if ((Correlated[m][v][syst]==kFALSE &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){	  
	    cout << "errore sist scorrelato in delta phi " << NSpectrumErrorSistUnCorr[m][v][syst][syst1]<< endl;
	  }
	  if ((Correlated[m][v][syst] &&BarlowPassed[m][v][syst]==kTRUE && syst< 12) ){
	    cout << "errore sist correlato in delta phi " << NSpectrumErrorSistCorr[m][v][syst][syst1]<< endl;
	  }

	}

      }//end loop on v
    }//end loop syst 1

    //*****************************************************************
    fHistSpectrum_master[m]=(TH1D*)    fHistSpectrum[m][0]->Clone(Form("fHistSpectrum_master_%i",m));
    fHistSpectrum_masterSystCorr[m]=      (TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_masterSystCorr"+Smolt[m]);
    fHistSpectrum_masterSystUnCorr[m]=      (TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_masterSystUnCorr"+Smolt[m]);
    fHistSpectrum_masterTotalUnCorr[m]=      (TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_masterTotalUnCorr"+Smolt[m]);
    fHistSpectrum_masterTotal[m]=      (TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_masterTotal"+Smolt[m]);

    for(Int_t syst1=0; syst1<syst1_limit; syst1++){
      fHistSpectrum_ratio[m][syst1]=      (TH1D*)fHistSpectrum[m][syst1]->Clone("fHistSpectrum_ratio"+Smolt[m]+Form("_SysBC%i",syst1));
      fHistSpectrum_Barlow[m][syst1]=      (TH1D*)fHistSpectrum[m][syst1]->Clone("fHistSpectrum_Barlow"+Smolt[m]+Form("_SysBC%i",syst1));
    }

   
    //Barlow check for spectrum
    for(Int_t syst=0; syst< syst1_limit; syst++){    
      //      cout << "\n\n\n ************************************syst "<<syst << endl;
      if(syst==0){	
	fHistSpectrum_ratio[m][syst]={0};
	fHistSpectrum_Barlow[m][syst]={0};
      }
      if(syst!=0){	
	for(Int_t j=2; j<= numPtV0Max; j++){

	  if(fHistSpectrum_master[m]->GetBinContent(j)!=0){
	    fHistSpectrum_ratio[m][syst]->SetBinContent(j,fHistSpectrum[m][syst]->GetBinContent(j)/fHistSpectrum_master[m]->GetBinContent(j));
	    fHistSpectrum_ratio[m][syst]->SetBinError(j,0);
	  }
	  if(sqrt(TMath::Abs(pow(fHistSpectrum[m][syst]->GetBinError(j),2)-pow(fHistSpectrum_master[m]->GetBinError(j),2)))!=0){
	    fHistSpectrum_Barlow[m][syst]->SetBinContent(j,(fHistSpectrum[m][syst]->GetBinContent(j)-fHistSpectrum_master[m]->GetBinContent(j))/sqrt(TMath::Abs(pow(fHistSpectrum[m][syst]->GetBinError(j),2)-pow(fHistSpectrum_master[m]->GetBinError(j),2))));
	    fHistSpectrum_Barlow[m][syst]->SetBinError(j,0);
	  }
	  else{
	    fHistSpectrum_Barlow[m][syst]->SetBinContent(j,0);
	    fHistSpectrum_Barlow[m][syst]->SetBinError(j,0);
	  }

	}
      }
	
      Int_t count=0;
      CorrelatedSpectrum[m][0]=kTRUE;
      fHistSpectrum_Corr[m][syst]=(TH1D*)fHistSpectrum[m][syst]->Clone(Form("fHistSpectrum_Corr_m%i_syst%i",m,syst));
      //	fHistSpectrum_Corr[m][syst]->SetBinContent(j,0);
      if (syst!=0) {
	for(Int_t j=2; j<= numPtV0Max; j++){
	  if (TMath::Abs(fHistSpectrum_Barlow[m][syst]->GetBinContent(j)) >= 2) count++;
	}
	if (count >=2) BarlowPassedSpectrum[m][syst]=kTRUE;
	if ( BarlowPassedSpectrum[m][syst]==kTRUE) {

	  NumberBarlowPassedSpectrum[syst]+=1;
	  for(Int_t j=1; j<= fHistSpectrum[m][syst]->GetNbinsX(); j++){
	    MeanSpectrum[m][syst] += fHistSpectrum_ratio[m][syst]->GetBinContent(j);
	  }
	  MeanSpectrum[m][syst]=MeanSpectrum[m][syst]/fHistSpectrum[m][syst]->GetNbinsX();
	  for(Int_t j=1; j<= fHistSpectrum[m][syst]->GetNbinsX(); j++){
	    SigmaSpectrum[m][syst] += pow(fHistSpectrum_ratio[m][syst]->GetBinContent(j)-MeanSpectrum[m][syst],2);
	  }
	  SigmaSpectrum[m][syst]=sqrt(SigmaSpectrum[m][syst]/(fHistSpectrum[m][syst]->GetNbinsX())/(fHistSpectrum[m][syst]->GetNbinsX()+1));
	  if(TMath::Abs((MeanSpectrum[m][syst]-1))>numSigmaCorr*SigmaSpectrum[m][syst]) CorrelatedSpectrum[m][syst]=kTRUE;
	  
	  //	    if (CorrelatedSpectrum[m][syst]==kFALSE){
	  //	  cout << "  m " << m<< "  syst " << syst << endl;
	  fHistSigmaBarlowSpectrum[m][syst]=(TH1D*)fHistSpectrum[m][syst]->Clone(Form("fHistSpectrum_SigmaBarlow_m%i_syst%i",m,syst));

	  for(Int_t j=1; j<= fHistSpectrum[m][syst]->GetNbinsX(); j++){
	    fHistSigmaBarlowSpectrum[m][syst]->SetBinContent(j,(fHistSpectrum[m][syst]->GetBinContent(j)-fHistSpectrum_master[m]->GetBinContent(j))/sqrt(12.));
	    fHistSpectrum_Corr[m][syst]->SetBinContent(j,0);
	  }
	  fHistSigmaBarlowSpectrum[m][syst]->Smooth();	    
	  //}
	  if (CorrelatedSpectrum[m][syst]==kTRUE){
	    NumberCorrSpectrum[m]+=1;
	  }
	}
      }
    }//end loop syst


    //***************************************************************************
    cout << "setting sist errors to spectrum " << endl;
    cout << "\n\nmult "<< m << endl;
    for(Int_t j=1; j< numPtV0Max; j++){
      SigmaSystSpectrumUnCorr[m][j]={0};
      SigmaTotalSpectrumUnCorr[m][j]={0};
      SigmaTotalSpectrum[m][j]={0};
      SigmaSystSpectrumCorr[m][j]={0};
      //cout << "************\n\n" << endl;
      for(Int_t syst=1; syst< syst1_limit; syst++){    
	// cout << " j " <<j << " syst  " <<syst << endl;
	if (CorrelatedSpectrum[m][syst]==kTRUE) continue;
	if ( BarlowPassedSpectrum[m][syst]==kFALSE) continue;
	//	cout << " componente non correlata m" << m << " v " << j << " syst " << syst << endl;
	//cout << fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1)        	<< endl;
	SigmaSystSpectrumUnCorr[m][j]+= pow(fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1),2);
      }

      SigmaSystSpectrumUnCorrBC[m][j]=      sqrt(SigmaSystSpectrumUnCorr[m][j]) ; //uncorr from BC

      SigmaSystSpectrumUnCorr[m][j]+= pow(NSpectrumErrorFinal[m][j][0],2); 

      SigmaTotalSpectrumUnCorr[m][j]=sqrt(SigmaSystSpectrumUnCorr[m][j]); //total uncorr (con stat)

      SigmaSystSpectrumUnCorr[m][j]=       SigmaSystSpectrumUnCorr[m][j] -pow(  NSpectrumErrorSoloStatFinal[m][j][0],2); //- stat

      SigmaSystSpectrumUnCorr[m][j]= sqrt(SigmaSystSpectrumUnCorr[m][j]);    
      cout << " m " << m << " j " <<  j <<  " syst uncorr in pt from BC only: " <<  SigmaSystSpectrumUnCorrBC[m][j] << " syst uncorr in pt: " <<  SigmaSystSpectrumUnCorr[m][j] << " syst uncorr in pt + stat: " <<       SigmaTotalSpectrumUnCorr[m][j]<< endl;
      cout << "\n\n" << endl;

      for(Int_t syst=1; syst< syst1_limit; syst++){    
	//cout <<m << "  " << j << "  " <<syst << endl;
	if (CorrelatedSpectrum[m][syst]==kFALSE) continue;
	if ( BarlowPassedSpectrum[m][syst]==kFALSE) continue;
	//	cout << " componente correlata " << m << " v " << j << " syst " << syst << endl;
	//cout << fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1)        	<< endl;
	SigmaSystSpectrumCorr[m][j]+= pow(fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1),2);
	// for(Int_t v =0; v < numPtV0; v++){
	//   YieldPerErrore[m][syst]+=   (fHistSpectrum_master[m]->GetBinContent(v+1) + fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(v+1));
	//   YieldDefault[m][syst]+=   (fHistSpectrum_master[m]->GetBinContent(v+1));
	// }
	// SigmaSystSpectrumCorrOK[m]+=pow(TMath::Abs(YieldPerErrore[m][syst] - YieldDefault[m][syst]),2);
	
      }
      //      SigmaSystSpectrumCorrOK[m]=sqrt( SigmaSystSpectrumCorrOK[m]);
      SigmaSystSpectrumCorrBC[m][j]=      sqrt(SigmaSystSpectrumCorr[m][j]) ; //corr from BC (no other errors correlated in pT are studied)
      SigmaTotalSpectrum[m][j]=    sqrt( pow( SigmaTotalSpectrumUnCorr[m][j],2) +pow(SigmaSystSpectrumCorr[m][j],2)); //total 
      SigmaSystSpectrumCorr[m][j]=sqrt(SigmaSystSpectrumCorr[m][j]);
      fHistSpectrum_masterSystCorr[m]->SetBinError(j+1,SigmaSystSpectrumCorr[m][j]);
      fHistSpectrum_masterSystUnCorr[m]->SetBinError(j+1,SigmaSystSpectrumUnCorr[m][j]);
      fHistSpectrum_masterTotalUnCorr[m]->SetBinError(j+1,SigmaTotalSpectrumUnCorr[m][j]);
      fHistSpectrum_masterTotal[m]->SetBinError(j+1,SigmaTotalSpectrum[m][j]);
      //      cout << " m " << m << " j " <<  j << " sist corr in pt " <<      SigmaSystSpectrumCorr[m][j] << endl;
    }


    //propago correttamente la componente correlata in pt dell'errore sull'istogramma yield vs mult
    for(Int_t syst=1; syst< syst1_limit; syst++){    
      //cout <<m << "  " << j << "  " <<syst << endl;
      if (CorrelatedSpectrum[m][syst]==kFALSE) continue;
      if ( BarlowPassedSpectrum[m][syst]==kFALSE) continue;
      for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
	YieldPerErrore[m][syst]+= (  (fHistSpectrum_master[m]->GetBinContent(v+1) + fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(v+1))*fHistSpectrum_master[m]->GetBinWidth(v+1));
	YieldDefault[m][syst]+=   (fHistSpectrum_master[m]->GetBinContent(v+1)*fHistSpectrum_master[m]->GetBinWidth(v+1));
      }
      SigmaSystSpectrumCorrOK[m]+=pow(TMath::Abs(YieldPerErrore[m][syst] - YieldDefault[m][syst]),2);
      
    }
    SigmaSystSpectrumCorrOK[m]=sqrt( SigmaSystSpectrumCorrOK[m]);    
  }//end loop m

  cout << "**********************************************************************"<< endl;
  cout << " end loop on m...starting to calculate yield... " << endl;
  cout << "**********************************************************************"<< endl;

  // cout << " i draw canvas " << endl;
  TCanvas *canvasSys[4];
  for (Int_t i=0; i < 4; i++){
    canvasSys[i]=new TCanvas(Form("canvas%i",i),Form("canvas%i",i), 1300, 1000);
    canvasSys[i]->Divide(3,2);
  }

  //yield vs multiplicity

  Double_t Yield[nummolt+1]={0};
  Double_t YieldErrTotal[nummolt+1]={0};
  Double_t YieldErrStat[nummolt+1]={0};
  Double_t YieldErrSist[nummolt+1]={0};
  Double_t YieldErrSistUnCorr[nummolt+1]={0};
  Double_t YieldErrSistCorr[nummolt+1]={0};
  TH1D* fHistYieldvsErrTotal=new TH1D ("fHistYieldvsErrTot","fHistYieldvsErrTot",300,0,30);
  TH1D* fHistYieldvsErrSoloStat=new TH1D ("fHistYieldvsErrSoloStat","fHistYieldvsErrSoloStat",300,0,30);
  TH1D* fHistYieldFitvsErrSoloStat=new TH1D ("fHistYieldFitvsErrSoloStat","fHistYieldFitvsErrSoloStat",300,0,30);
  TH1D* fHistYieldvsErrSoloSist=new TH1D ("fHistYieldvsErrSoloSist","fHistYieldvsErrSoloSist",300, 0, 30);
  TH1D* fHistYieldvsErrSoloStatMB=new TH1D ("fHistYieldvsErrSoloStatMB","fHistYieldvsErrSoloStatMB",300,0,30);
  TH1D* fHistYieldvsErrSoloSistMB=new TH1D ("fHistYieldvsErrSoloSistMB","fHistYieldvsErrSoloSistMB",300, 0, 30);
  TH1D* fHistYieldvsErrErrTotal=new TH1D ("fHistYieldvsErrErrTot","fHistYieldvsErrErrTot",300,0,30);
  TH1D* fHistYieldvsErrErrSoloStat=new TH1D ("fHistYieldvsErrErrSoloStat","fHistYieldvsErrErrSoloStat",300,0,30);
  TH1D* fHistYieldvsErrErrSoloSist=new TH1D ("fHistYieldvsErrErrSoloSist","fHistYieldvsErrErrSoloSist",300, 0, 30);
  TH1D* fHistYieldvsErrErrTotalMB   =new TH1D ("fHistYieldvsErrErrTotMB","fHistYieldvsErrErrTotMB",300,0,30);
  TH1D* fHistYieldvsErrErrSoloStatMB=new TH1D ("fHistYieldvsErrErrSoloStatMB","fHistYieldvsErrErrSoloStatMB",300,0,30);
  TH1D* fHistYieldvsErrErrSoloSistMB=new TH1D ("fHistYieldvsErrErrSoloSistMB","fHistYieldvsErrErrSoloSistMB",300, 0, 30);
  TH1D* fHistYieldvsErrTotalRatio=new TH1D ("fHistYieldvsErrTotRatio","fHistYieldvsErrTotRatio",300,0,30);
  TH1D* fHistYieldvsErrSoloStatRatio=new TH1D ("fHistYieldvsErrSoloStatRatio","fHistYieldvsErrSoloStatRatio",300,0,30);
  TH1D* fHistYieldvsErrSoloSistRatio=new TH1D ("fHistYieldvsErrSoloSistRatio","fHistYieldvsErrSoloSistRatio",300, 0, 30);

  TH1D* fHistYieldErrSistUnCorr=new TH1D ("fHistYieldErrSistUnCorr","fHistYieldErrSistUnCorr",300, 0, 30);
  TH1D* fHistYieldErrSistCorr=new TH1D ("fHistYieldErrSistCorr","fHistYieldErrSistCorr",300, 0, 30);

  TH1D* fHistYieldMultClass=new TH1D ("fHistYieldMultClass","fHistYieldMultClass",nummolt, Nmolt);
  
  Double_t multUniform[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};
  Double_t multChiara[nummolt+1]={0};
  Double_t multdistr[11]={0};
  Double_t multdistrtotal=0;
  Double_t limits[11]={0,1,5,10,15,20,30,40,50,70,100};
  TH1D* fHist_mult;
  fHist_mult=(TH1D*)list->FindObject("fHist_multiplicity");
  for(Int_t lim=0; lim< 10; lim++){
    for(Int_t l=fHist_mult->FindBin(limits[lim]+0.01); l<=fHist_mult->FindBin(limits[lim+1]-0.01); l++) {
      multdistr[lim]+=fHist_mult->GetBinContent(l);
    }  
    multdistrtotal+=multdistr[lim];
  }


  multChiara[0]= (26.02*multdistr[0]+20.02*multdistr[1])/(multdistr[0]+multdistr[1]);
  multChiara[1]= 16.17;
  multChiara[2]= (13.77*multdistr[3]+12.04*multdistr[4]+10.02*multdistr[5])/(multdistr[3]+multdistr[4]+multdistr[5]);
  multChiara[3]= (7.95*multdistr[6]+6.32*multdistr[7])/(multdistr[6]+multdistr[7]);
  multChiara[4]= (4.5*multdistr[8]+2.55*multdistr[9])/(multdistr[8]+multdistr[9]);
  multChiara[5]= (26.02*multdistr[0]+20.02*multdistr[1]+16.17*multdistr[2]+13.77*multdistr[3]+12.04*multdistr[4]+10.02*multdistr[5]+7.95*multdistr[6]+6.32*multdistr[7]+4.5*multdistr[8]+2.55*multdistr[9])/(multdistrtotal);

  //*********************rimetto mult tabulata
  Double_t   mult0[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};
  Double_t   mult1[nummolt+1]={26.02, 20.02, 14.97, 10.69, 4.4};

  Double_t   mult[nummolt+1]={26.02, 20.02, 14.97, 10.69, 4.4};
  for (Int_t m=0; m<nummolt; m++){
    if (MultBinning==0)mult[m] = mult0[m];
    else    if (MultBinning==1)mult[m] = mult1[m];
  }
  //*****************************************
  // for(Int_t j =1; j <=300; j++){
  //   fHistYieldvsErrSoloSist->SetBinContent(j, 0.00000001);
  //   fHistYieldvsErrSoloStat->SetBinContent(j, 0.00000001);
  //   fHistYieldvsErrTotal->SetBinContent(j, 0.00000001);
  
  // }



  //fit to the pt spectra (I fit the spectrum with errors stat + syst coming from all selections apart from choice of BC region ********************************+
  AliPWGFunc pwgfunc;
  Int_t numfittipo=5;
  Double_t YieldExtr[nummolt+1][numfittipo]={0};
  Double_t YieldExtrErrStat[nummolt+1][numfittipo]={0};
  Int_t ColorFit[numfittipo]={860, 881, 868, 628, 419};
  TLegend *legendfit=new TLegend(0.6, 0.6, 0.9, 0.9);
  TString	nameFit[numfittipo]={"mT-scaling", "pT-scaling", "Boltzmann", "Fermi-Dirac", "Levi"};//"Bose-Einstein"};
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);


  TFitResultPtr fFitResultPtr0[nummolt+1][syst1_limit][numfittipo]; 
  TF1* fit_MTscaling[nummolt+1][syst1_limit][numfittipo]; 
  Double_t YieldFitIntegralTypeFit[nummolt+1][numfittipo]={0};
  Double_t YieldFitIntegralErrorTypeFit[nummolt+1][numfittipo]={0};
  Double_t YieldFitIntegral[nummolt+1]={0};
  Double_t YieldFitIntegralError[nummolt+1]={0};

  //    LevyTsallis(kFitFunctionNames[0].data(), kParticleMass),
  //    pwgfunc.GetBoltzmann(kParticleMass, 0.1, 1, kFitFunctionNames[1].data()),
  //    pwgfunc.GetPTExp(0.1, 1, kFitFunctionNames[3].data()),
  TString   nameMTscaling[nummolt+1][syst1_limit][numfittipo];
  gStyle->SetOptFit(1111);
  for(Int_t m=0; m< nummolt +1; m++){
    //    if (m==0) continue;
    // if (m>=3) continue;     //added  
    legendCorrBis[m]= new TLegend(0.5, 0.7, 0.9, 0.9);
    legendCorrBis[m]->SetHeader("Selezioni applicate");
    canvasSys[0]->cd(m+1);
    gPad->SetLeftMargin(0.15);
    //       cout << "m " << m << endl;
    Int_t count4=0;
    for(Int_t syst=0; syst< syst1_limit; syst++){
      for (Int_t typefit =0; typefit<numfittipo; typefit++){
	if (    nameFit[typefit]== "pT-scaling") continue;
	pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
	nameMTscaling[m][syst][typefit] = Form("fitMTscaling_m%i_syst%i_fit%i",m,syst, typefit);
	cout << nameFit[typefit]<< endl;
	if (typefit==0)      fit_MTscaling[m][syst][typefit]=    pwgfunc.GetMTExp(massParticle[type], 0.1, 0.04, nameMTscaling[m][syst][typefit]); //mass, T, norm, name
	if (typefit==1)      fit_MTscaling[m][syst][typefit]=    pwgfunc.GetPTExp(0.1, 0.04, nameMTscaling[m][syst][typefit]); //mass, T, norm, name
	if (typefit==2)      fit_MTscaling[m][syst][typefit]=    pwgfunc.GetBoltzmann(massParticle[type],0.1, 0.04, nameMTscaling[m][syst][typefit]); 
	if (typefit==3)      fit_MTscaling[m][syst][typefit]=    pwgfunc.GetFermiDirac(massParticle[type],0.1, 0.04, nameMTscaling[m][syst][typefit]);
 	if (typefit==4)  {
	  fit_MTscaling[m][syst][typefit]=    pwgfunc.GetLevi(massParticle[type],0.1, 0.03, 0.04, nameMTscaling[m][syst][typefit]);	  //norm, n, T, mass (but the function must be called with these parameters in inverse order)

	  if (type==8) fit_MTscaling[m][syst][typefit]->SetParLimits(0, 0,fHistSpectrum[m][syst]->GetBinContent(fHistSpectrum[m][syst]->GetMaximumBin())*0.5*10);
	  else if (type==0)  fit_MTscaling[m][syst][typefit]->SetParLimits(0, 0,fHistSpectrum[m][syst]->GetBinContent(fHistSpectrum[m][syst]->GetMaximumBin())*0.5*10); 
	  fit_MTscaling[m][syst][typefit]->SetParLimits(2, 0.1, 10);
	  fit_MTscaling[m][syst][typefit]->SetParLimits(1, 2, 30);
	  fit_MTscaling[m][syst][typefit]->SetParameter(2, 0.7);
	}
	if (m==0 && syst==0)	legendfit->AddEntry( fit_MTscaling[m][syst][typefit], nameFit[typefit], "l"); 
	fit_MTscaling[m][syst][typefit]->SetLineColor(ColorFit[typefit]);
	fit_MTscaling[m][syst][typefit]->SetRange(LowRange[m], UpRange[m]);
	fFitResultPtr0[m][syst][typefit]=	fHistSpectrum[m][syst]->Fit(	fit_MTscaling[m][syst][typefit],"SR0");

	YieldFitIntegralTypeFit[m][typefit]=   fit_MTscaling[m][syst][typefit]->Integral(0,8);
	YieldFitIntegralErrorTypeFit[m][typefit]=   fit_MTscaling[m][syst][typefit]->IntegralError(0,8, fFitResultPtr0[m][syst][typefit] ->GetParams(),(fFitResultPtr0[m][syst][typefit]->GetCovarianceMatrix()).GetMatrixArray());

	YieldExtr[m][typefit]=   fit_MTscaling[m][0][typefit]->Integral(0,LowRangeSpectrumPart[m]);
	YieldExtrErrStat[m][typefit]=     fit_MTscaling[m][0][typefit]->IntegralError(0,LowRangeSpectrumPart[m], fFitResultPtr0[m][0][typefit] ->GetParams(),(fFitResultPtr0[m][0][typefit]->GetCovarianceMatrix()).GetMatrixArray());
    

	  cout << "\n\n\n******* hey I am there !! " << endl;
	if (nameFit[typefit]==FitFixed)  {
	  cout << "\n\n\n******* hey I am there !! " << endl;
	  YieldFitIntegral[m]=   	YieldFitIntegralTypeFit[m][typefit];
	  YieldFitIntegralError[m]=	YieldFitIntegralErrorTypeFit[m][typefit];
	}
	if (typefit==4) cout << "maximum of normalization " << fHistSpectrum[m][syst]->GetBinContent(fHistSpectrum[m][syst]->GetMaximumBin())*0.5*8 << endl;
	//	fHistSpectrum[m][syst]->Fit(	fit_MTscaling[m][syst][typefit],"", "",1,4);
	fit_MTscaling[m][syst][typefit]->SetRange(0,8);

	//fit_MTscaling[m][syst][typefit]->SetRange(0, 80);
	/*
	  cout << "m " << m << "chisquare/NDF " << fit_MTscaling[m][syst][typefit]->GetChisquare() << "/" << fit_MTscaling[m][syst][typefit]->GetNDF() <<"\n" <<  endl;
	  cout << "integral percentage for pT < 0.5 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 0.5)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
	  cout << "integral percentage for pT < 1.0 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 1)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
	  cout << "integral percentage for pT < 1.5 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 1.5)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
	  cout << "Integral up to 30 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
	*/
      }  

      //      fHistSpectrum[m][syst]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistSpectrum[m][syst]->GetXaxis()->SetTitle("");
      //      fHistSpectrum[m][syst]->GetYaxis()->SetTitle("1/N_{Trigg} dN/dp_{T} 1/#Delta#eta #Delta#phi");
      fHistSpectrum[m][syst]->GetYaxis()->SetTitle("");
      fHistSpectrum[m][syst]->GetYaxis()->SetTitleOffset(1.8);
      fHistSpectrum[m][syst]->SetTitle("p_{T} spectrum of K^{0}_{S} produced in jet in multiplicity range "+Smolt[m]+ "%");
      if (isBulk)      fHistSpectrum[m][syst]->SetTitle("p_{T} spectrum of K^{0}_{S} produced in UE in multiplicity range "+Smolt[m]+ "%");
      fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,0.04);
      fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,SpectrumSup[type][TypeAnalysis]);
      if (ishhCorr){ 
	if (!isBulk)       fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,0.3);
	if (isBulk)       fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,0.3);
	if (isTotal)       fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,0.3);
      }
      fHistSpectrum[m][syst]->SetMarkerStyle(Marker[syst]);
      fHistSpectrum[m][syst]->SetLineColor(Color[syst]);
      fHistSpectrum[m][syst]->SetMarkerColor(Color[syst]);
      //      fHistSpectrum[m][syst]->Draw("samee");
      fHistSpectrumPart[m][syst] =     (TH1D*)  fHistSpectrum[m][syst]->Clone(Form("fHistSpectrumPart_m%i_syst%i",m,syst));

      for(Int_t b=1;  b<fHistSpectrumPart[m][syst]->FindBin(LowRangeSpectrumPart[m]+0.001); b++){
	fHistSpectrumPart[m][syst]->SetBinContent(b,0);
	fHistSpectrumPart[m][syst]->SetBinError(b,0);
      }


      if (type==8 && m==1 && TypeAnalysis==0) {
	for(Int_t b=1;  b<fHistSpectrumPart[m][syst]->FindBin(8); b++){
	  if(  b==fHistSpectrumPart[m][syst]->FindBin(2.2)){
	    fHistSpectrumPart[m][syst]->SetBinContent(b,0);
	    fHistSpectrumPart[m][syst]->SetBinError(b,0);
	  }
	}
      }
 
      fHistSpectrumPart[m][syst]->SetTitle("Multiplicity class "+SmoltLegend[m]);
      if (type>0)      fHistSpectrumPart[m][syst]->GetXaxis()->SetTitle(Form("p^{#%s}_{T} (GeV/c)", tipoTitle[type].Data())); 
      else     fHistSpectrumPart[m][syst]->GetXaxis()->SetTitle(Form("p^{%s}_{T} (GeV/c)", tipoTitle[type].Data())); 
      fHistSpectrumPart[m][syst]->GetXaxis()->SetTitleSize(0.04);
      fHistSpectrumPart[m][syst]->GetXaxis()->SetTitleOffset(1);
      fHistSpectrumPart[m][syst]->Draw("samee");
      for (Int_t typefit =0; typefit<numfittipo; typefit++){
	if (    nameFit[typefit]== "pT-scaling") continue;
	//	if (    nameFit[typefit]!= "Fermi-Dirac") continue;
	fit_MTscaling[m][syst][typefit]->Draw("same");
      }
      legendfit->Draw("");
      if (!isBulk)	legendCorrBis[m]->AddEntry(fHistSpectrum[m][syst],SSystJet[syst],"pel");   
      if (isBulk)	legendCorrBis[m]->AddEntry(fHistSpectrum[m][syst],SSystBulk[syst],"pel");   
      count4++;
      if(syst==1) legendCorrBis[m]->Draw();
    }
  }


  TFile *fileout = new TFile(stringout ,"RECREATE");
  for(Int_t m=0; m< nummolt +1; m++){
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (    nameFit[typefit]== "pT-scaling") continue;
      cout << " saving the TF1 " << endl;
      fileout->WriteTObject(fHistSpectrumPart[m][0]);
      fileout  ->WriteTObject(fit_MTscaling[m][0][typefit]);
    }
  }
  fileout-> WriteTObject(  fHistBinCountingRegionAndPtFit);
  fileout->WriteTObject( canvasSys[0]); 

  TH1F *  fHistSpectrumRatioFit[nummolt+1][syst1_limit][numfittipo];
  TCanvas *canvasRatioFit = new TCanvas("canvasRatioFit", "canvasRatioFit", 800, 500);
  canvasRatioFit->Divide(3,2);
  TF1 * rettaUno= new TF1("rettaUno", "pol0",0,8);
  rettaUno->FixParameter(0,1);
  for(Int_t m=0; m< nummolt +1; m++){
    //    if (m==0) continue;
    canvasRatioFit->cd(m+1);
    gPad->SetLeftMargin(0.15);
    for(Int_t syst=0; syst< syst1_limit; syst++){
      for (Int_t typefit =0; typefit<numfittipo; typefit++){
	if (    nameFit[typefit]== "pT-scaling") continue;
	//	fHistSpectrumRatioFit[m][syst][typefit] = (TH1F*) fHistSpectrum[m][syst]->Clone(Form("HistRatioFit_m%i_syst%i_typefit%i", m, syst, typefit));
	fHistSpectrumRatioFit[m][syst][typefit] = (TH1F*) fHistSpectrumPart[m][syst]->Clone(Form("HistRatioFit_m%i_syst%i_typefit%i", m, syst, typefit));
	fHistSpectrumRatioFit[m][syst][typefit]->Divide(fit_MTscaling[m][syst][typefit]);
	fHistSpectrumRatioFit[m][syst][typefit]->GetXaxis()->SetTitle("");
	fHistSpectrumRatioFit[m][syst][typefit]->GetYaxis()->SetTitle("");
	fHistSpectrumRatioFit[m][syst][typefit]->GetYaxis()->SetTitleOffset(1.8);
	fHistSpectrumRatioFit[m][syst][typefit]->GetYaxis()->SetRangeUser(0.4,2);
	if (TypeAnalysis==0)      fHistSpectrumRatioFit[m][syst][typefit]->GetYaxis()->SetRangeUser(0.4,2);
	fHistSpectrumRatioFit[m][syst][typefit]->SetLineColor(ColorFit[typefit]);
	fHistSpectrumRatioFit[m][syst][typefit]->SetMarkerColor(ColorFit[typefit]);
	fHistSpectrumRatioFit[m][syst][typefit]->Draw("samee");
	rettaUno->SetLineColor(kBlack);
	if (typefit==0)      rettaUno->Draw("same");
      }
      legendfit->Draw("");
    }
  }
  fileout->WriteTObject( canvasRatioFit); 

  cout << "\n\n Yield vs multiplicity *******************************\n"<< endl;
  for(Int_t m=0; m < nummolt+1; m++){
    //    if (m==0) continue;
    // if (m>=3) continue;     //added  
    Int_t typefitFixed=0;
    for (Int_t typefit=0; typefit<numfittipo; typefit++){
      if (nameFit[typefit]==FitFixed)    typefitFixed=typefit;
    }
    Yield[m]+=   fit_MTscaling[m][0][typefitFixed]->Integral(0,LowRangeSpectrumPart[m]);
    YieldErrStat[m]+= pow(    fit_MTscaling[m][0][typefitFixed]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][0][typefitFixed] ->GetParams(),(fFitResultPtr0[m][0][typefitFixed]->GetCovarianceMatrix()).GetMatrixArray() ), 2);
    
    if (type==8 && m==1 && TypeAnalysis==0) {
      Yield[m]+=   fit_MTscaling[m][0][typefitFixed]->Integral(2, 2.5);
    }  
    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
    if (type==8 && m==1 && v!=fHistSpectrumPart[m][0]->FindBin(2.2) && TypeAnalysis==0) {
      Yield[m]+=  ( fHistSpectrum_master[m]->GetBinContent(v+1)*fHistSpectrum_master[m]->GetBinWidth(v+1));
      YieldErrStat[m]+= pow(   fHistSpectrum_masterOnlyStat[m]->GetBinError(v+1)*fHistSpectrum_master[m]->GetBinWidth(v+1),2);
      YieldErrSistUnCorr[m]+=pow(    fHistSpectrum_masterSystUnCorr[m]->GetBinError(v+1)*fHistSpectrum_master[m]->GetBinWidth(v+1),2);
    }
    else {
      Yield[m]+=  ( fHistSpectrum_master[m]->GetBinContent(v+1)*fHistSpectrum_master[m]->GetBinWidth(v+1));
      YieldErrStat[m]+= pow(   fHistSpectrum_masterOnlyStat[m]->GetBinError(v+1)*fHistSpectrum_master[m]->GetBinWidth(v+1),2);
      YieldErrSistUnCorr[m]+=pow(    fHistSpectrum_masterSystUnCorr[m]->GetBinError(v+1)*fHistSpectrum_master[m]->GetBinWidth(v+1),2);
    }
    }
    YieldErrSistCorr[m]= pow(SigmaSystSpectrumCorrOK[m],2);
    YieldErrStat[m]=sqrt( YieldErrStat[m]);
    YieldErrSistUnCorr[m]=sqrt( YieldErrSistUnCorr[m]);
    YieldErrSistCorr[m]=sqrt( YieldErrSistCorr[m]);
    //smooth degli errori sistematici con histo intermedi
    fHistYieldErrSistUnCorr->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]),YieldErrSistUnCorr[m]);
    fHistYieldErrSistCorr->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]),YieldErrSistCorr[m]);
    //*****************  


  }
  //se inserisco smooth errore sistematico va a zero (perché tanti bin valgono zero)
  // fHistYieldErrSistUnCorr->Smooth();
  // fHistYieldErrSistCorr->Smooth();
  for(Int_t m=0; m < nummolt+1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    YieldErrSistUnCorr[m]=  fHistYieldErrSistUnCorr->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]));
    YieldErrSistCorr[m]=  fHistYieldErrSistCorr->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]));
    YieldErrSist[m]=sqrt(pow(YieldErrSistUnCorr[m],2) + pow(YieldErrSistCorr[m],2));
    YieldErrTotal[m]=sqrt( pow(YieldErrStat[m],2)+ pow(YieldErrSistUnCorr[m],2)+ pow(YieldErrSistCorr[m],2));  
  }



  cout << "***************************"<< endl;
  for(Int_t m=0; m < nummolt+1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    if (m < nummolt){
      fHistYieldvsErrSoloSist->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);
      fHistYieldvsErrSoloStat->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);  
      fHistYieldvsErrTotal->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);
      fHistYieldvsErrSoloSist->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrSistUnCorr[m],2) + pow(YieldErrSistCorr[m],2)));
      fHistYieldvsErrSoloStat->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrStat[m]);  
      fHistYieldvsErrTotal->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrTotal[m]);

      //      fHistYieldFitvsErrSoloSist->SetBinContent(fHistYieldFitvsErrTotal->FindBin(mult[m]), YieldFitIntegral[m]);
      fHistYieldFitvsErrSoloStat->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), YieldFitIntegral[m]);  
      //      fHistYieldFitvsErrTotal->SetBinContent(fHistYieldFitvsErrTotal->FindBin(mult[m]), YieldFitIntegral[m]);
      fHistYieldFitvsErrSoloStat->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), YieldFitIntegralError[m]);  

      fHistYieldvsErrSoloSistRatio->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]/Yield[5]);
      fHistYieldvsErrSoloStatRatio->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]/Yield[5]);  
      fHistYieldvsErrTotalRatio->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]/Yield[5]);
      fHistYieldvsErrSoloSistRatio->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrSist[m]/Yield[5],2)+pow(YieldErrSist[5]/pow(Yield[5],2)*Yield[m],2)));
      fHistYieldvsErrSoloStatRatio->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrStat[m]/Yield[5],2)+pow(YieldErrStat[5]/pow(Yield[5],2)*Yield[m],2)));
      fHistYieldvsErrTotalRatio->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrTotal[m]/Yield[5],2)+pow(YieldErrTotal[5]/pow(Yield[5],2)*Yield[m],2)));

      fHistYieldvsErrErrSoloSist->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrSistUnCorr[m],2) + pow(YieldErrSistCorr[m],2))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m])));
      fHistYieldvsErrErrSoloStat->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrStat[m]/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m])));  
      fHistYieldvsErrErrTotal->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrTotal[m]/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m])));
      fHistYieldvsErrErrSoloSist->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), 0);
      fHistYieldvsErrErrSoloStat->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), 0);
      fHistYieldvsErrErrTotal->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), 0);

    }
    else{
      fHistYieldvsErrSoloSist->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), 0);
      fHistYieldvsErrSoloStat->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), 0);  
      fHistYieldFitvsErrSoloStat->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), 0);  
      fHistYieldvsErrTotal->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), 0);
      fHistYieldvsErrSoloSistMB->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);
      fHistYieldvsErrSoloStatMB->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);  
      fHistYieldvsErrSoloSistMB->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrSistUnCorr[m],2) + pow(YieldErrSistCorr[m],2)));
      fHistYieldvsErrSoloStatMB->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrStat[m]);  
      fHistYieldvsErrErrSoloSistMB->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrSistUnCorr[m],2) + pow(YieldErrSistCorr[m],2))/fHistYieldvsErrSoloSistMB->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m])));
      fHistYieldvsErrErrSoloStatMB->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrStat[m]/fHistYieldvsErrSoloSistMB->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m])));  
      fHistYieldvsErrErrTotalMB->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrTotal[m]/fHistYieldvsErrSoloSistMB->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m])));
      fHistYieldvsErrErrSoloSistMB->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), 0);
      fHistYieldvsErrErrSoloStatMB->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), 0);
      fHistYieldvsErrErrTotalMB->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), 0);
    
    }


    cout << "\n************  m=" << m << endl;
    cout << "valori dei vettori, errori assoluti"  << endl;
    cout << "yield           "<<     Yield[m] << endl;
    cout << "error stat      "<<     YieldErrStat[m]<< endl;
    cout << "error uncorr pt " <<    YieldErrSistUnCorr[m]<< endl;
    cout << "error   corr pt " <<    YieldErrSistCorr[m]<< endl;
    cout << "error sist      "<<     YieldErrSist[m]<< endl;
    cout << "error total     "<<     YieldErrTotal[m]<< endl;
    if (m< nummolt) {
      cout << "valori degli istogrammi, errori relativi " << endl;
      cout << "yield           "<<     fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "error stat rel  "<<     fHistYieldvsErrSoloStat->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "error stat rel yield from fit  "<<     fHistYieldFitvsErrSoloStat->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "error sist rel  "<<     fHistYieldvsErrSoloSist->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "error total rel "<<     fHistYieldvsErrTotal->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "*********" << endl;
      // cout << "gli stessi errori ma assoluti" << endl;
      // cout <<  fHistYieldvsErrErrSoloStat->GetBinContent( fHistYieldvsErrSoloStat->FindBin(mult[m])) << endl;
      // cout <<  fHistYieldvsErrErrSoloSist->GetBinContent( fHistYieldvsErrSoloSist->FindBin(mult[m])) << endl;
      // cout <<  fHistYieldvsErrErrTotal->GetBinContent( fHistYieldvsErrTotal->FindBin(mult[m])) << endl;
    }
  }

  //  cout << "\n\n*************fit functions****************"<< endl;

  TF1 *expo = new TF1("expo", "expo", 0,30);
  fHistYieldvsErrTotal->Fit(expo,"R+");
  //  cout << "\nchi quadro/NDF expo" << expo->GetChisquare() << "/" << expo->GetNDF()<< endl;
  //  cout << expo->GetParameter(0) << " exp^" << expo->GetParameter(1) << endl;
  // TF1 *expo1 = new TF1("expo1", "[0]*TMath::Exp([1]*x) + [2]", 0,30);
  //  expo1->SetLineColor(3);
  //  fHistYieldvsErrTotal->Fit(expo1,"R+");
  //  cout << "\nchi quadro/NDF expo1" << expo->GetChisquare() << "/" << expo->GetNDF()<< endl;
  //  cout << expo1->GetParameter(0) << " exp^" << expo1->GetParameter(1)<< " + " <<expo1->GetParameter(2)<< endl;
  TF1 *power = new TF1("", "[1]*x**[0]",1,30);
  power->SetLineColor(1);
  fHistYieldvsErrTotal->Fit(power,"R+");
  //  cout << "\nchi quadro/NDF power law" << power->GetChisquare() << "/" << power->GetNDF()<< endl;
  //  cout << power->GetParameter(1) << " x^" << power->GetParameter(0) << endl;

  TF1 *rettayield = new TF1("rettayield", "pol1", 0, 30);
  rettayield->SetLineColor(3);
  fHistYieldvsErrTotal->Fit(rettayield,"R+");
  //  cout << "\nchi quadro/NDF rettayield" << rettayield->GetChisquare() << "/" << rettayield->GetNDF()<< endl;


  //canvas: yield con errori statistici e sistematici*****************************************************
  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 1300,1000);
  canvasYield->cd();
  canvasYield->SetLogy();
  canvasYield->SetLogx();
  fHistYieldvsErrSoloStat->SetTitle ("Yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
  fHistYieldvsErrSoloStat->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
  fHistYieldvsErrSoloStat->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi");
  //  if (isBulk || isTotal)  fHistYieldvsErrSoloStat->GetYaxis()->SetTitle ("C x N_{K^{0}_{S}}/N_{Trigg}");
  fHistYieldvsErrSoloStat->GetYaxis()->SetTitleOffset(1.1);
  fHistYieldvsErrSoloStat->GetXaxis()->SetTitleOffset(1.1);
  fHistYieldvsErrSoloStat->GetYaxis()->SetTitleSize(0.042);

  fHistYieldFitvsErrSoloStat->SetTitle ("YieldFit of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
  fHistYieldFitvsErrSoloStat->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
  fHistYieldFitvsErrSoloStat->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi");
  //  if (isBulk || isTotal)  fHistYieldFitvsErrSoloStat->GetYaxis()->SetTitle ("C x N_{K^{0}_{S}}/N_{Trigg}");
  fHistYieldFitvsErrSoloStat->GetYaxis()->SetTitleOffset(1.1);
  fHistYieldFitvsErrSoloStat->GetXaxis()->SetTitleOffset(1.1);
  fHistYieldFitvsErrSoloStat->GetYaxis()->SetTitleSize(0.042);

  fHistYieldvsErrSoloSist->SetTitle ("Yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
  fHistYieldvsErrSoloSist->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
  fHistYieldvsErrSoloSist->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi");
  fHistYieldvsErrSoloSist->GetYaxis()->SetTitleOffset(1.1);
  fHistYieldvsErrSoloSist->GetXaxis()->SetTitleOffset(1.1);
  fHistYieldvsErrSoloSist->GetYaxis()->SetTitleSize(0.042);

  fHistYieldvsErrSoloStatRatio->SetTitle ("Relative yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
  fHistYieldvsErrSoloStatRatio->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
  fHistYieldvsErrSoloStatRatio->GetYaxis()->SetTitle ("(N_{K^{0}_{S}}/N_{Trigg}) / (N_{K^{0}_{S}}/N_{Trigg})_{0-100 %}");
  fHistYieldvsErrSoloStatRatio->GetYaxis()->SetTitleOffset(1.1);
  fHistYieldvsErrSoloStatRatio->GetYaxis()->SetTitleSize(0.042);
  fHistYieldvsErrSoloStatRatio->GetXaxis()->SetTitleOffset(1.1);

  fHistYieldvsErrSoloSistRatio->SetTitle ("Relative yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
  fHistYieldvsErrSoloSistRatio->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
  fHistYieldvsErrSoloSistRatio->GetYaxis()->SetTitle ("(N_{K^{0}_{S}}/N_{Trigg}) / (N_{K^{0}_{S}}/N_{Trigg})_{0-100 %}");
  fHistYieldvsErrSoloSistRatio->GetYaxis()->SetTitleOffset(1.1);
  fHistYieldvsErrSoloSistRatio->GetYaxis()->SetTitleSize(0.042);
  fHistYieldvsErrSoloSistRatio->GetXaxis()->SetTitleOffset(1.1);

  fHistYieldvsErrSoloStat->GetYaxis()->SetRangeUser(0.02,0.4);
  fHistYieldvsErrSoloStat->GetXaxis()->SetRangeUser(1,50);
  fHistYieldvsErrSoloStat->SetMarkerStyle(1);
  fHistYieldvsErrSoloStat->SetLineColor(1);
  fHistYieldvsErrSoloStat->SetMarkerColor(1);
  fHistYieldvsErrSoloStat->Draw("samee");
  fHistYieldFitvsErrSoloStat->GetYaxis()->SetRangeUser(0.02,0.4);
  fHistYieldFitvsErrSoloStat->GetXaxis()->SetRangeUser(1,50);
  fHistYieldFitvsErrSoloStat->SetMarkerStyle(1);
  fHistYieldFitvsErrSoloStat->SetLineColor(1);
  fHistYieldFitvsErrSoloStat->SetMarkerColor(1);
  fHistYieldFitvsErrSoloStat->Draw("samee");

  fHistYieldvsErrSoloStatMB->GetYaxis()->SetRangeUser(0.02,0.4);
  fHistYieldvsErrSoloStatMB->GetXaxis()->SetRangeUser(1,50);
  fHistYieldvsErrSoloStatMB->SetMarkerStyle(1);
  fHistYieldvsErrSoloStatMB->SetLineColor(1);
  fHistYieldvsErrSoloStatMB->SetMarkerColor(1);

  fHistYieldvsErrSoloSistMB->GetYaxis()->SetRangeUser(0.02,0.4);
  fHistYieldvsErrSoloSistMB->GetXaxis()->SetRangeUser(1,50);
  fHistYieldvsErrSoloSistMB->SetMarkerStyle(1);
  fHistYieldvsErrSoloSistMB->SetLineColor(1);
  fHistYieldvsErrSoloSistMB->SetMarkerColor(1);

  fHistYieldvsErrSoloSist->GetYaxis()->SetRangeUser(0.02,0.4);
  fHistYieldvsErrSoloSist->GetXaxis()->SetRangeUser(1,50);
  fHistYieldvsErrSoloSist->SetMarkerStyle(1);
  fHistYieldvsErrSoloSist->SetLineColor(1);
  fHistYieldvsErrSoloSist->SetMarkerColor(1);
  fHistYieldvsErrSoloSist->SetFillStyle(0);
  
  fHistYieldvsErrSoloSist->Draw("samee2");
  // fHistYieldvsErrTotal->GetYaxis()->SetRangeUser(0,0.2);
  // fHistYieldvsErrTotal->SetMarkerStyle(20);
  // fHistYieldvsErrTotal->SetLineColor(1);
  // fHistYieldvsErrTotal->SetMarkerColor(3);
  // fHistYieldvsErrTotal->SetFillStyle(0);
  // fHistYieldvsErrTotal->Draw("samee2");
  expo->Draw("same");
  //    expo1->Draw("same");
  power->Draw("same");

  legendYield->AddEntry(expo, SFit[0], "pl");
  legendYield->AddEntry(power, SFit[1], "pl");
  legendYield->Draw();  

  //  cout << "\n\n\ndrawing canvas " << endl;

  //canvas: errori relativi dello yield, statistici, sistematici e totali*********************************
  TCanvas *canvasYieldError = new TCanvas("canvasYieldError", "canvasYieldError", 1300,1000);

  canvasYieldError->cd();

  fHistYieldvsErrErrSoloStat->SetTitle ("Relative uncertainty vs V0M multiplicity");
  fHistYieldvsErrErrSoloStat->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
  fHistYieldvsErrErrSoloStat->GetYaxis()->SetTitle ("Relative uncertainty");
  fHistYieldvsErrErrSoloStat->GetYaxis()->SetRangeUser(0,0.165);
  fHistYieldvsErrErrSoloStat->GetXaxis()->SetRangeUser(1,50);
  fHistYieldvsErrErrSoloStat->SetMarkerStyle(23);
  fHistYieldvsErrErrSoloStat->SetLineColor(1);
  fHistYieldvsErrErrSoloStat->SetMarkerColor(1);
  fHistYieldvsErrErrSoloStat->Draw("samep");
  legendYieldError->AddEntry(    fHistYieldvsErrErrSoloStat, "stat.", "pel");

  fHistYieldvsErrErrSoloSist->SetTitle ("Relative uncertainty vs V0M multiplicity");
  fHistYieldvsErrErrSoloSist->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
  fHistYieldvsErrErrSoloSist->GetYaxis()->SetTitle ("Relative uncertainty");
  fHistYieldvsErrErrSoloSist->GetYaxis()->SetRangeUser(0,0.165);
  fHistYieldvsErrErrSoloSist->GetXaxis()->SetRangeUser(1,50);
  fHistYieldvsErrErrSoloSist->SetMarkerStyle(20);
  fHistYieldvsErrErrSoloSist->SetLineColor(2);
  fHistYieldvsErrErrSoloSist->SetMarkerColor(2);
  fHistYieldvsErrErrSoloSist->SetFillStyle(0);  
  fHistYieldvsErrErrSoloSist->Draw("samep");
  legendYieldError->AddEntry(    fHistYieldvsErrErrSoloSist, "syst.", "pel");

  fHistYieldvsErrErrTotal->GetYaxis()->SetRangeUser(0,0.165);
  fHistYieldvsErrErrTotal->GetXaxis()->SetRangeUser(1,50);
  fHistYieldvsErrErrTotal->SetMarkerStyle(24);
  fHistYieldvsErrErrTotal->SetLineColor(4);
  fHistYieldvsErrErrTotal->SetMarkerColor(4);
  fHistYieldvsErrErrTotal->SetFillStyle(0);
  fHistYieldvsErrErrTotal->Draw("samep");
  legendYieldError->AddEntry(    fHistYieldvsErrErrTotal, "total", "pel");
   
  legendYieldError->Draw();  

  //  cout << "\n\n\ndrawing canvas " << endl;

  //  TFile *fileout = new TFile(stringout ,"RECREATE");
  /*
    for (Int_t m =0; m<nummolt+1; m++){
    for(Int_t syst=0; syst<12; syst++){
    fHistSpectrumErrorSistUnCorr[m][syst][0]->Write();
    fHistSpectrumErrorSistCorr[m][syst][0]->Write();
    }
    }
  */
  // fileout->WriteTObject( fHistYieldvsErrSoloStat);
  //fileout->WriteTObject( fHistYieldvsErrSoloSist);
  fHistYieldvsErrSoloSist->Write();
  fHistYieldvsErrSoloStat->Write();
  fHistYieldFitvsErrSoloStat->Write();
  fHistYieldvsErrSoloSistMB->Write();
  fHistYieldvsErrSoloStatMB->Write();
  fHistYieldvsErrTotal->Write();
  fHistYieldvsErrSoloSistRatio->Write();
  fHistYieldvsErrSoloStatRatio->Write();
  fHistYieldvsErrTotalRatio->Write();
  fHistYieldvsErrErrSoloSist->Write();
  fHistYieldvsErrErrSoloStat->Write();
  fHistYieldvsErrErrTotal->Write();
  fHistYieldvsErrErrSoloSistMB->Write();
  fHistYieldvsErrErrSoloStatMB->Write();
  fHistYieldvsErrErrTotalMB->Write();
  fileout->WriteTObject( canvasYield);
  fileout->WriteTObject( canvasYieldError);
  TF1* retta = new TF1("retta", "pol0", -10, 10);
  retta->FixParameter(0,1);
  retta->SetLineColor(1);
  retta->SetLineWidth(0.1);
  TF1* retta1 = new TF1("retta1", "pol0", -10, 10);
  retta1->FixParameter(0,0);
  retta1->SetLineColor(1);
  retta1->SetLineWidth(0.1);

  legendErrorSpectrum= new TLegend(0.5, 0.65, 0.9, 0.9);
  //legendErrorSpectrum->SetHeader("Type of error");
  
  for(Int_t m=nummolt; m>=0; m--){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    canvasSys[1]->cd(m+1);
    gPad->SetLeftMargin(0.15);

    fHistSpectrum_masterOnlyStat[m]->GetXaxis()->SetTitle("p^{Assoc}_{T} (GeV/c)");
    //    fHistSpectrum_masterOnlyStat[m]->GetYaxis()->SetTitle("1/#Delta#eta #Delta#phi 1/N_{Trigg} dN/dp_{T}");
    fHistSpectrum_masterOnlyStat[m]->GetYaxis()->SetTitleOffset(1.6);
    fHistSpectrum_masterOnlyStat[m]->GetXaxis()->SetTitleOffset(1.05);
    fHistSpectrum_masterOnlyStat[m]->GetYaxis()->SetTitleSize(0.045);
    fHistSpectrum_masterOnlyStat[m]->GetXaxis()->SetTitleSize(0.045);
    fHistSpectrum_masterOnlyStat[m]->SetTitle("Multiplicity class "+SmoltLegend[m]);
    fHistSpectrum_masterOnlyStat[m]->GetYaxis()->SetRangeUser(0,SpectrumSup[type][TypeAnalysis]);
    if (ishhCorr){
      if (!isBulk)    fHistSpectrum_masterOnlyStat[m]->GetYaxis()->SetRangeUser(0,0.3);
      if (isBulk)    fHistSpectrum_masterOnlyStat[m]->GetYaxis()->SetRangeUser(0,3);
      if (isTotal)    fHistSpectrum_masterOnlyStat[m]->GetYaxis()->SetRangeUser(0,3);
    }
    fHistSpectrum_masterOnlyStat[m]->SetMarkerStyle(1);
    fHistSpectrum_masterOnlyStat[m]->SetLineColor(1);
    fHistSpectrum_masterOnlyStat[m]->SetMarkerColor(1);
   
    fHistSpectrum_masterSystCorr[m]->GetXaxis()->SetTitle("p^{Assoc}_{T} (GeV/c)");
    //    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetTitle("1/#Delta#eta #Delta#phi 1/N_{Trigg} dN/dp_{T}");
    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetTitleOffset(1.6);
    fHistSpectrum_masterSystCorr[m]->GetXaxis()->SetTitleOffset(1.05);
    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetTitleSize(0.045);
    fHistSpectrum_masterSystCorr[m]->GetXaxis()->SetTitleSize(0.045);
    fHistSpectrum_masterSystCorr[m]->SetTitle("Multiplicity class "+SmoltLegend[m]);
    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,SpectrumSup[type][TypeAnalysis]);
    if (ishhCorr){
      if (!isBulk)    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,0.3);
      if (isBulk)    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,3);
      if (isTotal)    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,3);
    }
    fHistSpectrum_masterSystCorr[m]->SetMarkerStyle(1);
    fHistSpectrum_masterSystCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystCorr[m]->SetFillStyle(9);
    fHistSpectrum_masterSystCorr[m]->SetFillStyle(1001);
    if (SigmaSystSpectrumCorr[m][0]!=0)    fHistSpectrum_masterSystCorr[m]->SetFillColorAlpha(kAzure-4,0.003);

    fHistSpectrum_masterSystUnCorr[m]->GetXaxis()->SetTitle("p^{Assoc}_{T} (GeV/c)");
    //    fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetTitle("1/#Delta#eta #Delta#phi 1/N_{Trigg} dN/dp_{T}");
    fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetTitleOffset(1.6);
    fHistSpectrum_masterSystUnCorr[m]->GetXaxis()->SetTitleOffset(1.05);
    fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetTitleSize(0.045);
    fHistSpectrum_masterSystUnCorr[m]->GetXaxis()->SetTitleSize(0.045);
    fHistSpectrum_masterSystUnCorr[m]->SetTitle("Multiplicity class "+SmoltLegend[m]);
    fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetRangeUser(0,SpectrumSup[type][TypeAnalysis]);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerStyle(1);
    fHistSpectrum_masterSystUnCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerColor(1);
    fHistSpectrum_masterSystUnCorr[m]->SetFillStyle(0);
    //    fHistSpectrum_masterSystUnCorr[m]->SetFillColorAlpha(kWhite,0.003);
    //    fHistSpectrum_masterSystUnCorr[m]->SetFillStyle(9);

    fHistSpectrum_masterSystCorr[m]->Draw("samee2 p");
    fHistSpectrum_masterSystUnCorr[m]->Draw("samee2 p");  
    fHistSpectrum_masterOnlyStat[m]->Draw("samee p");
    fileout->WriteTObject(fHistSpectrum_master[m]);


    if(m==5){
      legendErrorSpectrum->AddEntry(fHistSpectrum_masterOnlyStat[m],SErrorSpectrum[0], "el");  
      //      if (!(isBulk==0 && isMC==0 && isTotal==0) && !(isBulk==0 && isMC==1 && isTotal==0))      legendErrorSpectrum->AddEntry(fHistSpectrum_masterSystCorr[m],SErrorSpectrum[2], "f");   
      if (isBulk)      legendErrorSpectrum->AddEntry(fHistSpectrum_masterSystCorr[m],SErrorSpectrum[2], "f");   
      legendErrorSpectrum->AddEntry(fHistSpectrum_masterSystUnCorr[m],SErrorSpectrum[1], "f");   
    }
    //    legendErrorSpectrum->Draw();
  }

  canvasSys[1]->SaveAs("FinalOutput/DATA"+year0 +"/PtSpectrum"+JetOrBulk[jet]+DataOrMC[isMC]+".pdf");
  fileout->WriteTObject( canvasSys[1]);


  canvasSpectrumSysBis[0]=new TCanvas("canvasSpectrumSysBis0","canvasSpectrumSysBis0", 1300, 1000);
  canvasSpectrumSysBis[1]=new TCanvas("canvasSpectrumSysBis1","canvasSpectrumSysBis1", 1300, 1000);
  canvasSpectrumSysBis[0]->Divide(3,3);
  canvasSpectrumSysBis[1]->Divide(3,3);

  for(Int_t m=0; m< nummolt +1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    canvasSpectrumSys[2*m]=new TCanvas(Form("canvasSpectrumSys0_%i",m),Form("canvasSpectrumSys0_%i",m), 1300, 1000);
    canvasSpectrumSys[2*m+1]=new TCanvas(Form("canvasSpectrumSys1_%i",m),Form("canvasSpectrumSys1_%i",m), 1300, 1000);
    canvasSpectrumSys[2*m]->Divide(4,3);
    canvasSpectrumSys[2*m+1]->Divide(4,3);

    canvasPhiSys[4*m]=new TCanvas(Form("canvasSys0_%i",m),Form("canvasSys0_%i",m), 1300, 1000);
    canvasPhiSys[4*m+1]=new TCanvas(Form("canvasSys1_%i",m),Form("canvasSys1_%i",m), 1300, 1000);
    canvasPhiSys[4*m+2]=new TCanvas(Form("canvasSys2_%i",m),Form("canvasSys2_%i",m), 1300, 1000);
    canvasPhiSys[4*m+3]=new TCanvas(Form("canvasSys3_%i",m),Form("canvasSys3_%i",m), 1300, 1000);

    canvasPhiSys[4*m]->Divide(2,3);
    canvasPhiSys[4*m+1]->Divide(2,3);
    canvasPhiSys[4*m+2]->Divide(2,3);
    canvasPhiSys[4*m+3]->Divide(2,3);


  }
  for(Int_t m=0; m< nummolt +1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    legend_Bpassed_Spectrum[m]= new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_Corr_Spectrum[m]= new TLegend(0.7, 0.7, 0.9, 0.9);

    //
    Int_t counter=0;
    Int_t counter2=0;
    for(Int_t syst=0; syst<syst1_limit; syst++){
	  
      // if (syst==9 || syst == avoidthissyst) continue;
      // if (CorrelatedBis[m][syst]==kFALSE) continue;
      // if (syst>0 && syst < 12) continue; 
      if (m<3)	  canvasSpectrumSysBis[0]->cd(m+1);
      else 	  canvasSpectrumSysBis[1]->cd(m+1-3);
      gPad->SetLeftMargin(0.15);
      fHistSpectrum[m][syst]->GetXaxis()->SetTitle("p^{Assoc}_{T} (GeV/c)");
      fHistSpectrum[m][syst]->GetXaxis()->SetTitleOffset(1.1);
      fHistSpectrum[m][syst]->GetYaxis()->SetTitleOffset(1.6);
      fHistSpectrum[m][syst]->GetYaxis()->SetTitleSize(0.045);
      fHistSpectrum[m][syst]->GetXaxis()->SetTitleSize(0.045);
      //      fHistSpectrum[m][syst]->GetYaxis()->SetTitle("1/N_{Trigg} dN/dp_{T}");
      //      if (isBulk) fHistSpectrum[m][syst]->GetYaxis()->SetTitle("C x 1/N_{Trigg} dN/dp_{T}");
      fHistSpectrum[m][syst]->SetTitle("Multiplicity class " + SmoltLegend[m]);
      fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,0.04);
      fHistSpectrum[m][syst]->SetMarkerStyle(Marker[syst]);
      fHistSpectrum[m][syst]->SetLineColor(Color[syst]);
      fHistSpectrum[m][syst]->SetMarkerColor(Color[syst]);
      fHistSpectrum[m][syst]->Draw("samee");
      //      if(syst==syst1_limit -1 ) legendCorrBis[m]->Draw();

      if(BarlowPassedSpectrum[m][syst]==kTRUE) {
	counter++;
	if(syst!=0){
	  if (m<3) canvasSpectrumSysBis[0]->cd(6+m+1);
	  else canvasSpectrumSysBis[1]->cd(6+m+1-3);
	  fHistSpectrum_Barlow[m][syst]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	  fHistSpectrum_Barlow[m][syst]->GetYaxis()->SetTitle("Barlow variable");
	  fHistSpectrum_Barlow[m][syst]->SetTitle("Barlow variable distribution for significant systematic effects");
	  fHistSpectrum_Barlow[m][syst]->GetYaxis()->SetRangeUser(-30,30);
	  fHistSpectrum_Barlow[m][syst]->SetMarkerStyle(Marker[syst]);
	  fHistSpectrum_Barlow[m][syst]->SetLineColor(Color[syst]);
	  fHistSpectrum_Barlow[m][syst]->SetMarkerColor(Color[syst]);
	  fHistSpectrum_Barlow[m][syst]->Draw("samep");
	  retta1->SetLineColor(1);
	  retta1->SetLineWidth(0.1);
	  retta1->Draw("same");
	  if (!ishhCorr)	  legend_Bpassed_Spectrum[m]->AddEntry(fHistSpectrum_Barlow[m][syst],SSyst[syst],"pel");   
	  if (ishhCorr)	  legend_Bpassed_Spectrum[m]->AddEntry(fHistSpectrum_Barlow[m][syst],SSysthh[syst],"pel");   
	  //	  if(counter==NumberBarlowPassedSpectrum[m]) legend_Bpassed_Spectrum[m]->Draw();

	  if (CorrelatedSpectrum[m][syst]==kTRUE){	  
	    counter2++;
	    if (m<3) canvasSpectrumSysBis[0]->cd(3+m+1);
	    else canvasSpectrumSysBis[1]->cd(3+m+1-3);
	    fHistSpectrum_ratio[m][syst]->GetYaxis()->SetRangeUser(-2,4);
	    fHistSpectrum_ratio[m][syst]->SetMarkerStyle(Marker[syst]);
	    fHistSpectrum_ratio[m][syst]->SetLineColor(Color[syst]);
	    fHistSpectrum_ratio[m][syst]->SetMarkerColor(Color[syst]);
	    fHistSpectrum_ratio[m][syst]->Draw("samep");
	    retta->SetLineColor(1);
	    retta->SetLineWidth(0.1);
	    retta->Draw("same");
	    if (!ishhCorr) 	    legend_Corr_Spectrum[m]->AddEntry(fHistSpectrum_ratio[m][syst],SSyst[syst],"pel");   
	    if (ishhCorr) 	    legend_Corr_Spectrum[m]->AddEntry(fHistSpectrum_ratio[m][syst],SSysthh[syst],"pel");   
	    //  cout << counter2 << " "<<NumberCorr[m] << endl;
	    //	    if(counter2==NumberCorrSpectrum[m]) legend_Corr_Spectrum[m]->Draw();
	  }
	  
	}
      }
    }
  }


  for(Int_t m=0; m< nummolt +1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      legend_corr[m][v]= new TLegend(0.6, 0.7, 0.9, 0.9);
      legend_corr[m][v]->SetHeader("Selezioni applicate");
      legend_Bpassed[m][v]= new TLegend(0.6, 0.7, 0.9, 0.9);
      legend_Bpassed[m][v]->SetHeader("Selezioni applicate");
      legend_UnCorr[m][v]= new TLegend(0.6, 0.7, 0.9, 0.9);
      legend_UnCorr[m][v]->SetHeader("Selezioni applicate");
	
      Int_t counter=0;
      Int_t counter2=0;
      Int_t counter_unc=0;
      //      cout << "******************** " << m << "  " << v << endl;
      for(Int_t syst=0; syst<numSysV0Global+numsystPhiGlobal; syst++){
	if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	if (v<2)	  canvasPhiSys[4*m]->cd(v+1);
	else 	  if (v<4)	  canvasPhiSys[4*m+1]->cd(v+1-2);
	else 	  if (v<6)	  canvasPhiSys[4*m+2]->cd(v+1-4);
	else	  canvasPhiSys[4*m+3]->cd(v+1-6);
	if (!isBulk)	fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,8000);
	if (isBulk)	fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(0,7000);
	if (m==5)	fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,8000);
	if (type==4|| type==5 || type==8){
	  if (!isBulk)	fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,12000);
	  if (isBulk)	fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(0,12000);
	}
	fHistPhiDistr[m][v][syst]->GetXaxis()->SetTitle("#Delta#phi");
	fHistPhiDistr[m][v][syst]->GetYaxis()->SetTitle("Counts");
	fHistPhiDistr[m][v][syst]->SetMarkerStyle(Marker[syst]);
	fHistPhiDistr[m][v][syst]->SetLineColor(Color[syst]);
	fHistPhiDistr[m][v][syst]->SetMarkerColor(Color[syst]);
	if (m==0 && v==0 && !ishhCorr) legend3->AddEntry(fHistPhiDistr[m][v][syst],SSyst[syst],"pel");   
	if (m==0 && v==0 && ishhCorr) legend3->AddEntry(fHistPhiDistr[m][v][syst],SSysthh[syst],"pel");   

	fHistPhiDistr[m][v][syst]->Draw("samee");
	if(syst == numSysV0Global+numsystPhiGlobal-1) legend3->Draw();

	if(BarlowPassed[m][v][syst]==kTRUE) {
	  counter++;
	  if(syst!=0){
	    if (v<2)	  canvasPhiSys[4*m]->cd(4+v+1);
	    else 	  if (v<4)	  canvasPhiSys[4*m+1]->cd(4+v+1-2);
	    else 	  if (v<6)	  canvasPhiSys[4*m+2]->cd(4+v+1-4);
	    else	  canvasPhiSys[4*m+3]->cd(4+v+1-6);
	    fHistPhiDistr_Barlow[m][v][syst]->GetYaxis()->SetRangeUser(-30,30);
	    fHistPhiDistr_Barlow[m][v][syst]->SetTitle("Barlow variable distribution for significant systematic effects");	   
	    fHistPhiDistr_Barlow[m][v][syst]->GetXaxis()->SetTitle("#Delta#phi");	   
	    fHistPhiDistr_Barlow[m][v][syst]->GetYaxis()->SetTitle("Barlow variable");
	    fHistPhiDistr_Barlow[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistPhiDistr_Barlow[m][v][syst]->SetLineColor(Color[syst]);
	    fHistPhiDistr_Barlow[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistPhiDistr_Barlow[m][v][syst]->Draw("samep");
	    retta1->Draw("same"); 
	    if (!ishhCorr)	    legend_Bpassed[m][v]->AddEntry(fHistPhiDistr_Barlow[m][v][syst],SSyst[syst],"pel");   
	    if (ishhCorr)	    legend_Bpassed[m][v]->AddEntry(fHistPhiDistr_Barlow[m][v][syst],SSysthh[syst],"pel");   
	    if(counter==NumberBarlowPassed[m][v]) legend_Bpassed[m][v]->Draw();
	    //cout << counter << " "<<NumberBarlowPassed[m][v] << endl;

	    if (Correlated[m][v][syst]==kTRUE){	  
	      counter2++;
	      if (v<2)	  canvasPhiSys[4*m]->cd(2+v+1);
	      else 	  if (v<4)	  canvasPhiSys[4*m+1]->cd(2+v+1-2);
	      else 	  if (v<6)	  canvasPhiSys[4*m+2]->cd(2+v+1-4);
	      else	  canvasPhiSys[4*m+3]->cd(2+v+1-6);
	      if (isBulk) fHistPhiDistr_ratio[m][v][syst]->GetYaxis()->SetRangeUser(0.5,1.5);
	      fHistPhiDistr_ratio[m][v][syst]->SetTitle("Ratio of #Delta#phi distribution to the default one for systematic effects correlated across #Delta #phi");	   
	      fHistPhiDistr_ratio[m][v][syst]->GetXaxis()->SetTitle("#Delta#phi");	   
	      fHistPhiDistr_ratio[m][v][syst]->GetYaxis()->SetTitle("Ratio");
	      
	      fHistPhiDistr_ratio[m][v][syst]->GetYaxis()->SetRangeUser(-2,4);
	      fHistPhiDistr_ratio[m][v][syst]->SetMarkerStyle(Marker[syst]);
	      fHistPhiDistr_ratio[m][v][syst]->SetLineColor(Color[syst]);
	      fHistPhiDistr_ratio[m][v][syst]->SetMarkerColor(Color[syst]);
	      fHistPhiDistr_ratio[m][v][syst]->Draw("samep");
	      retta->SetLineColor(1);
	      retta->SetLineWidth(0.1);
	      retta->Draw("same");
	      if (!ishhCorr) 	      legend_corr[m][v]->AddEntry(fHistPhiDistr_ratio[m][v][syst],SSyst[syst],"pel");   
	      if (ishhCorr) 	      legend_corr[m][v]->AddEntry(fHistPhiDistr_ratio[m][v][syst],SSysthh[syst],"pel");   
	      //  cout << counter2 << " "<<NumberCorr[m][v] << endl;
	      if(counter2==NumberCorr[m][v]) legend_corr[m][v]->Draw();
	    }
	  
	  }
	

	  if (Correlated[m][v][syst]==kTRUE){ //remember: for syst=0 was put kTRUE
	    // cout << "\n\ndrawing canvas " << endl;
	    // cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	    if (v<4)  canvasSpectrumSys[2*m]->cd(v+1);
	    else  canvasSpectrumSys[2*m+1]->cd(v+1-4);
	    fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,8000);
	    if (m==5)	  fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,20000);
	    fHistPhiDistr[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistPhiDistr[m][v][syst]->SetLineColor(Color[syst]);
	    fHistPhiDistr[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistPhiDistr[m][v][syst]->Draw("samee");
	    //	  if(counter2==NumberCorr[m][v]) legend_corr[m][v]->Draw();
	    if(syst==0) legend_corr[m][v]->Draw();
	 
	    /* 
	       if(syst!=0){
	       if (v<4)	  canvasSpectrumSys[2*m]->cd(4+v+1);
	       else	  canvasSpectrumSys[2*m+1]->cd(4+v+1-4);
	       fHistPhiDistr_ratio[m][v][syst]->GetYaxis()->SetRangeUser(-5,5);
	       fHistPhiDistr_ratio[m][v][syst]->SetMarkerStyle(Marker[syst]);
	       fHistPhiDistr_ratio[m][v][syst]->SetLineColor(Color[syst]);
	       fHistPhiDistr_ratio[m][v][syst]->SetMarkerColor(Color[syst]);
	       fHistPhiDistr_ratio[m][v][syst]->Draw("samee");

	       }
	    */
	    /*
	      if(syst==0){
	      //	    cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	      if (v<4)	  canvasSpectrumSys[2*m]->cd(4+v+1);
	      else	  canvasSpectrumSys[2*m+1]->cd(4+v+1-4);
	      fHistPhiDistr_master[m][v]->GetYaxis()->SetRangeUser(-1000,8000);
	      if (m==5) fHistPhiDistr_master[m][v]->GetYaxis()->SetRangeUser(-1000,20000);
	      fHistPhiDistr_master[m][v]->SetMarkerStyle(Marker[syst]);
	      fHistPhiDistr_master[m][v]->SetLineColor(Color[syst]);
	      fHistPhiDistr_master[m][v]->SetMarkerColor(Color[syst]);
	      fHistPhiDistr_master[m][v]->Draw("samee");

	      }
	    */


  
	  }

	  //	  if (syst==0 || Correlated[m][v][syst]==kFALSE){ //remember: for syst=0 was put kTRUE
	  if (kTRUE){
	    //  cout << "syst " << syst << endl;
	    counter_unc++;
	    if (v<4)	  canvasSpectrumSys[2*m]->cd(8+v+1);
	    else	  canvasSpectrumSys[2*m+1]->cd(8+v+1-4);
	    fHistSigmaSystNSmoothed[m][v][syst]->GetYaxis()->SetRangeUser(0,2);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetLineColor(Color[syst]);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetTitle("Relative uncorrelated systematic uncertainty");
	    fHistSigmaSystNSmoothed[m][v][syst]->GetXaxis()->SetTitle("#Delta#phi");
	    fHistSigmaSystNSmoothed[m][v][syst]->GetYaxis()->SetTitle("Relative uncertainty");
	    fHistSigmaSystNSmoothed[m][v][syst]->Draw("same");
	    if (!ishhCorr && syst !=0 &&  Correlated[m][v][syst]==kFALSE)	    legend_UnCorr[m][v]->AddEntry(fHistSigmaSystNSmoothed[m][v][syst],SSyst[syst],"pel");   
	    if (ishhCorr && syst !=0 &&  Correlated[m][v][syst]==kFALSE)	    legend_UnCorr[m][v]->AddEntry(fHistSigmaSystNSmoothed[m][v][syst],SSysthh[syst],"pel");   
	    else  legend_UnCorr[m][v]->AddEntry(fHistSigmaSystNSmoothed[m][v][syst],"total uncorrelated uncertainty","pel");   
	    if(counter_unc==(NumberBarlowPassed[m][v]-NumberCorr[m][v])) legend_UnCorr[m][v]->Draw();

	    //	cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	    if (v<4)	  canvasSpectrumSys[2*m]->cd(4+v+1);
	    else	  canvasSpectrumSys[2*m+1]->cd(4+v+1-4);
	    fHistSigmaSyst[m][v][syst]->GetYaxis()->SetRangeUser(0,400);
	    fHistSigmaSyst[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistSigmaSyst[m][v][syst]->SetLineColor(Color[syst]);
	    fHistSigmaSyst[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistSigmaSyst[m][v][syst]->Draw("sameep");
	    //	cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	  }
	}

      }
    }
  }//end loop on mult


  cout << "\n\n ********************** Values of spectrum + errors ; If you are interested, decomment me!! " << endl;
  cout << "master val, master_bis val, master_bis error(stat), errore stat , errore stat + sist uncorr(deltaphi), errore stat +sist(delta phi) , errore uncorr (sist), error uncorr sist + stat, error corr (sist) "<< endl;


  for(Int_t m=0; m< nummolt +1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    fHistSpectrumErrorSistBC[m]=(TH1D*)fHistSpectrum[m][0]->Clone(Form("fHistSpectrumErrorSistBC_m%i",m));
    fHistSpectrumErrorSistPhi[m]=(TH1D*)fHistSpectrum[m][0]->Clone(Form("fHistSpectrumErrorSistPhi_m%i",m));
    cout << "*********************************mult " <<m << endl;
    for(Int_t j=1; j < numPtV0Max; j++){
      for(Int_t syst=0; syst<syst1_limit; syst++){
	if(syst==9) continue;
	//	cout << "syst " << syst << "  " << fHistSpectrum_master[m]->GetBinContent(j+1) << "  " << fHistSpectrum[m][syst]->GetBinContent(j+1) <<  "  " <<  ALowBinFit[syst] << "  " << AUpBinFit[syst] << "  " << endl;//<<fHistPhiDistr_Corr[m][j][syst]->FindBin(ALowBinFit[syst])<<  endl;

      }
      //        cout << "^^^^" << j<< "        " << fHistSpectrum_master[m]->GetBinContent(j+1)<< "  " <<  "     " << fHistSpectrum_masterOnlyStat[m]->GetBinContent(j+1)<< "     " <<fHistSpectrum_masterOnlyStat[m]->GetBinError(j+1) << "  " << NSpectrumErrorSoloStatFinal[m][j][0]<< "  <<   " <<NSpectrumErrorSistUnCorr[m][j][0] << "  " <<NSpectrumErrorFinal[m][j][0]<< " <= " << SigmaSystSpectrumUnCorr[m][j] << " < " << SigmaTotalSpectrumUnCorr[m][j] << " !=0 but not consid.     " << SigmaSystSpectrumCorr[m][j] << "          "   << endl;
    
      //      cout << "\n\n ********************** Values of spectrum + errors ; If you are interested, decomment me!! " << endl;
      /*
      cout << "^^^^" << j<< " spectrum value: " << fHistSpectrum_master[m]->GetBinContent(j+1)<< ", stat rel: " <<fHistSpectrum_masterOnlyStat[m]->GetBinError(j+1)/fHistSpectrum_master[m]->GetBinContent(j+1) <<", sist uncorr phi (rel): " <<NSpectrumErrorSistUnCorrSI[m][j][0]/fHistSpectrum_master[m]->GetBinContent(j+1) << ", sist corr phi (rel): "<<NSpectrumErrorSistCorrSI[m][j][0]/ fHistSpectrum_masterOnlyStat[m]->GetBinContent(j+1) << ", sist phi (rel): " <<sqrt(pow(NSpectrumErrorFinal[m][j][0],2) - pow(fHistSpectrum_masterOnlyStat[m]->GetBinError(j+1),2))/fHistSpectrum_master[m]->GetBinContent(j+1)<<", sist uncorr pt (BC only): " << SigmaSystSpectrumUnCorrBC[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1)<<  ", sist uncorr pt (rel): " << SigmaSystSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1) << ", tot uncorr (sist+stat)  (rel): " << SigmaTotalSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1) << ", tot (sist + stat)(rel): " << SigmaTotalSpectrum[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1) << ", sist corr pt (rel): " << SigmaSystSpectrumCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1)  << endl;
      */
      for(Int_t syst=0; syst< numsystPhiGlobal + numSysV0Global; syst++){    
	//cout << "\n\n\n ************************************syst "<<syst << endl
	if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
     	//inserisco errori sist correlati in delta phi nello spettro relativo a syst =0
	if ((Correlated[m][j][syst]==kTRUE && syst< 12) ){
	  //	  cout << "sist corr pt x syst " << syst << "  " << SSyst[syst] << "  " <<	    NSpectrumErrorSistCorr[m][j][syst][0]<< "  "<< endl;
	}
      }
      //      cout << "*************" << endl;
      fHistSpectrumErrorSistBC[m]->SetBinContent(j+1,sqrt(pow(SigmaSystSpectrumUnCorrBC[m][j],2) + pow(SigmaSystSpectrumCorr[m][j],2))/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrumErrorSistBC[m]->SetBinError(j+1,0);

      fHistSpectrumErrorSistPhi[m]->SetBinContent(j+1,sqrt(pow(NSpectrumErrorFinal[m][j][0],2) - pow(fHistSpectrum_masterOnlyStat[m]->GetBinError(j+1),2))/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrumErrorSistPhi[m]->SetBinError(j+1,0);

    }
    //   fHistSpectrumErrorSistBC[m]->Smooth();
  }
  //  cout << "end loop on m " << endl;

  TLegend *legenderrorspectrum = new TLegend(0.15, 0.7,0.45,0.9);
  //    legenderrorspectrum->SetHeader("");


  for(Int_t m=0; m < nummolt+1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    //cout << "\n\n"<< m << endl;
    for(Int_t syst=1; syst< numsystPhiGlobal + numSysV0Global; syst++){  
      // cout << syst << endl;  
      if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
      if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
      if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
      fHistSpectrumErrorSist[m][syst]=(TH1D*)fHistSpectrum[m][syst]->Clone(Form("fHistSpectrumErrorSistUnCorr_m%i_syst%i",m,syst));
    }
  
    for(Int_t j =1; j < numPtV0Max; j++){
      // cout << j << endl;
      fHistSpectrum_masterOnlyStatRel[m]= (TH1D*)      fHistSpectrum_masterOnlyStat[m]->Clone("fHistSpectrum_"+Smolt[m]+"_Rel_eta<0.5");
      fHistSpectrum_masterOnlyStatRel[m]->SetBinContent(j+1,NSpectrumErrorSoloStatFinal[m][j][0]/fHistSpectrum_master[m]->GetBinContent(j+1));
      //      fHistSpectrum_masterOnlyStat[m]->SetBinError(j+1,sqrt(pow(NSpectrumErrorFinal[m][j][0],2) - pow(fHistSpectrum_masterOnlyStat[m]->GetBinError(j+1),2))/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrum_masterSystCorr[m]->SetBinContent(j+1,SigmaSystSpectrumCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrum_masterSystUnCorr[m]->SetBinContent(j+1,SigmaSystSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrum_masterTotalUnCorr[m]->SetBinContent(j+1,SigmaTotalSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrum_masterTotal[m]->SetBinContent(j+1,sqrt(pow(SigmaTotalSpectrumUnCorr[m][j],2) + pow(SigmaSystSpectrumCorr[m][j],2))/fHistSpectrum_master[m]->GetBinContent(j+1));

      for(Int_t syst=1; syst< numsystPhiGlobal + numSysV0Global; syst++){  
	// cout << syst << endl;  
	if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
	if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
	if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
     
	fHistSpectrumErrorSist[m][syst]->SetBinContent(j+1,0);
	fHistSpectrumErrorSist[m][syst]->SetBinError(j+1,0);
	if ((BarlowPassed[m][j][syst]==kTRUE && syst< 12) ){
	  if (!Correlated[m][j][syst])fHistSpectrumErrorSist[m][syst]->SetBinContent(j+1,NSpectrumErrorSistUnCorr[m][j][syst][0]/fHistSpectrum_master[m]->GetBinContent(j+1));
	  if (Correlated[m][j][syst])fHistSpectrumErrorSist[m][syst]->SetBinContent(j+1,NSpectrumErrorSistCorr[m][j][syst][0]/fHistSpectrum_master[m]->GetBinContent(j+1));
	  fHistSpectrumErrorSist[m][syst]->SetBinError(j+1,0);
	  //	  cout << " m " << m << " syst " << syst << " j " << j << "  " << NSpectrumErrorSist[m][j][syst][0]/fHistSpectrum_master[m]->GetBinContent(j+1) << endl;  
	  //	  cout << " m " << m << " syst " << syst << " j " << j << "  " << fHistSpectrumErrorSist[m][syst]->GetBinContent(j+1) << endl;  
	}
      }
    }
  }

  cout << "ecco le differenze tra molteplicità" << endl;
  for(Int_t m=0; m< nummolt+1; m++){
    //    cout << m << "   " <<multUniform[m] << "  " << mult[m]<< endl;
  }


  TCanvas *canvasSysError=new TCanvas("canvasSysError", "canvasSysError",1300,1000);
  canvasSysError->Divide(3,2);

  for(Int_t m=0; m < nummolt+1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    canvasSys[2]->cd(m+1);
    gPad->SetLeftMargin(0.15);
    fHistSpectrum_masterOnlyStatRel[m]->GetXaxis()->SetTitle("p^{Assoc}_{T} (GeV/c)");
    fHistSpectrum_masterOnlyStatRel[m]->GetYaxis()->SetTitle("Relative uncertainty");
    fHistSpectrum_masterOnlyStatRel[m]->GetYaxis()->SetTitleOffset(1.3);
    fHistSpectrum_masterOnlyStatRel[m]->GetYaxis()->SetTitleSize(0.045);
    //    fHistSpectrum_masterOnlyStatRel[m]->GetXaxis()->SetTitleOffset(1.3);
    fHistSpectrum_masterOnlyStatRel[m]->GetXaxis()->SetTitleSize(0.048);
    fHistSpectrum_masterOnlyStatRel[m]->SetTitle("Relative systematic and statistic uncertainties ");
    fHistSpectrum_masterOnlyStatRel[m]->SetTitle("Multiplicity class "+SmoltLegend[m]);
    fHistSpectrum_masterOnlyStatRel[m]->GetYaxis()->SetRangeUser(0,0.45);
    if (isBulk )  fHistSpectrum_masterOnlyStatRel[m]->GetYaxis()->SetRangeUser(0,0.8);
    if (isBulk && isMC)  fHistSpectrum_masterOnlyStatRel[m]->GetYaxis()->SetRangeUser(0,2);
    fHistSpectrum_masterOnlyStatRel[m]->SetMarkerStyle(20);
    fHistSpectrum_masterOnlyStatRel[m]->SetLineColor(1);
    fHistSpectrum_masterOnlyStatRel[m]->SetMarkerColor(2);
    fHistSpectrum_masterOnlyStatRel[m]->Draw("same");
    if (m==0)legenderrorspectrum->AddEntry(    fHistSpectrum_masterOnlyStatRel[m], "stat. " , "pel");
    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,0.45);
    if (isBulk)    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,0.8);
    if (isBulk && isMC)    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,2);
    fHistSpectrum_masterSystCorr[m]->SetMarkerStyle(25);
    fHistSpectrum_masterSystCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystCorr[m]->SetMarkerColor(9);
    fHistSpectrum_masterSystCorr[m]->SetFillStyle(0);
    fHistSpectrum_masterSystCorr[m]->Draw("same");
    if (isBulk){
      if (m==0) legenderrorspectrum->AddEntry(    fHistSpectrum_masterSystCorr[m], "syst. corr. " , "pel");
    }
    fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetRangeUser(0,0.45);
    if (isBulk)     fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetRangeUser(0,0.8);
    if (isBulk && isMC)     fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetRangeUser(0,2);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerStyle(20);
    fHistSpectrum_masterSystUnCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerColor(3);
    fHistSpectrum_masterSystUnCorr[m]->SetFillStyle(0);
    fHistSpectrum_masterSystUnCorr[m]->Draw("same");
    if (m==0)legenderrorspectrum->AddEntry(    fHistSpectrum_masterSystUnCorr[m], "syst. uncorr. " , "pel");
    fHistSpectrum_masterTotal[m]->GetYaxis()->SetRangeUser(0,0.45);
    fHistSpectrum_masterTotal[m]->SetMarkerStyle(20);
    fHistSpectrum_masterTotal[m]->SetLineColor(1);
    fHistSpectrum_masterTotal[m]->SetMarkerColor(4);
    fHistSpectrum_masterTotal[m]->SetFillStyle(0);
    fHistSpectrum_masterTotal[m]->Draw("same");
    if (m==0) legenderrorspectrum->AddEntry(    fHistSpectrum_masterTotal[m], "total " , "pel");
    legenderrorspectrum->Draw();
  }
   
  fileout->WriteTObject( canvasSys[2]);
  canvasSys[2]->SaveAs("FinalOutput/DATA"+year0 +"/RelSysPt"+JetOrBulk[jet]+DataOrMC[isMC]+"Bis.pdf");

  //***************************+saving canvas to file *********************************

  // for(Int_t i =0; i<4; i++){
  // fileout->WriteTObject( canvasSys[i]);
  // }

  fileout->WriteTObject( canvasSpectrumSysBis[0]);
  fileout->WriteTObject( canvasSpectrumSysBis[1]);

  for(Int_t m=0; m< nummolt +1; m++){
    //if (m==0) continue;

    // if (m>=3) continue;     //added  
    fileout->WriteTObject( canvasSpectrumSys[2*m]);
    fileout->WriteTObject( canvasSpectrumSys[2*m+1]);
  }
  for(Int_t m=0; m< nummolt +1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    fileout->WriteTObject( canvasPhiSys[4*m]);
    fileout->WriteTObject( canvasPhiSys[4*m+1]);
    fileout->WriteTObject( canvasPhiSys[4*m+2]);
    fileout->WriteTObject( canvasPhiSys[4*m+3]);
  }



  /*
    for(Int_t m=0; m< nummolt +1; m++){
    cout << "*********************************mult " <<m << endl;
    for(Int_t j=0; j < numPtV0; j++){
    cout << " j " << j << endl;
    for(Int_t syst=1; syst< numsystPhiGlobal + numSysV0Global; syst++){   
    cout << "\n\n\n ************************************syst "<<syst << endl;
    if (syst==9 || syst == avoidthissyst || syst == avoidthissystbis) continue;
    //inserisco errori sist correlati in delta phi nello spettro relativo a syst =0
    if ((Correlated[m][j][syst]==kTRUE && syst< 12) ){
    cout << "sist corr pt x syst " << syst << "  " << SSyst[syst] << "  " <<	    NSpectrumErrorSistCorr[m][j][syst][0]<< "  "<<endl;
    //	  cout << fHistSpectrumErrorSist[m][syst]->GetBinContent(j+1) << endl;
    }
    }
    //cout << "*************" << endl;
    }
    }
  */

  cout << "\n\n\n facciamo tabella per la tesi!!!!" << endl;
  for(Int_t syst=1; syst< numsystPhiGlobal + numSysV0Global; syst++){   
    SigmaBarlowMin[syst]=100000; 
  }

  Double_t MinValue=0;
  Double_t MaxValue=0;
  if (isBulk || isTotal ){
    MinValue=ALowBinFit[0];
    MaxValue=AUpBinFit[0];
  }
  else{
    MinValue=-0.2;
    MaxValue=0.2;
  }
  for(Int_t syst=1; syst< numsystPhiGlobal + numSysV0Global; syst++){   
    if (ishhCorr && (syst==3 || syst==4 || syst==5)) continue; //isto non definiti per questi sistematici
    if (!ishhCorr && (syst==9 || syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
    if (ishhCorr && (syst ==avoidthissyst || syst==avoidthissystbis || syst==avoidthissysttris)) continue;
    for(Int_t m=0; m< nummolt +1; m++){
      //if (m==0) continue;
      // if (m>=3) continue;     //added  
      //cout << "*********************************mult " <<m << endl;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if ( BarlowPassed[m][v][syst]==kFALSE) continue;
	for(Int_t j=fHistPhiDistr_Corr[m][v][syst]->FindBin(MinValue); j<= fHistPhiDistr_Corr[m][v][syst]->FindBin(MaxValue); j++){
	  if ((fHistPhiDistr[m][v][syst]->GetBinContent(j+1)) ==0) continue;
	  //cout << "\n\n\n ************************************syst "<<syst << endl;

	  //if (m==0 && v ==0 && j==fHistPhiDistr_master_Corr[m][v][syst]->FindBin(ALowBinFit[0])){
	  //SigmaBarlowMax[syst]=	 TMath::Abs( SigmaBarlow[m][v][syst][j]/fHistPhiDistr_master[m][v][syst]->GetBinContent(j+1));
	  //  SigmaBarlowMin[syst]=	1000000000;
	  // }
	  //	  if (SigmaBarlow[m][v][syst][j] ==0) cout <<m << " " << v << " " << syst << "  " << j << " " << SigmaBarlow[m][v][syst][j]<< endl;
	  if (SigmaBarlowMax[syst]< TMath::Abs(SigmaBarlow[m][v][syst][j]/fHistPhiDistr_master[m][v]->GetBinContent(j+1))) SigmaBarlowMax[syst]=TMath::Abs(SigmaBarlow[m][v][syst][j]/fHistPhiDistr_master[m][v]->GetBinContent(j+1));
	  //cout <<m << " " << v << " " << syst << "  " << j << " " << SigmaBarlowMin[syst]<< endl;
	  if (SigmaBarlow[m][v][syst][j]!=0 && SigmaBarlowMin[syst]> TMath::Abs(SigmaBarlow[m][v][syst][j]/fHistPhiDistr_master[m][v]->GetBinContent(j+1))) SigmaBarlowMin[syst]=TMath::Abs(SigmaBarlow[m][v][syst][j]/fHistPhiDistr_master[m][v]->GetBinContent(j+1));
	}
      }
      //cout << "*************" << endl;
    }
  }

  for(Int_t syst=1; syst< numsystPhiGlobal + numSysV0Global; syst++){   
    //    cout << syst << "  " << SSyst[syst]<<  "  " <<SigmaBarlowMin[syst] << "  " << SigmaBarlowMax[syst]<< endl;
  }


 

  cout << "\n\ndelta phi width " << DeltaPhiWidth[0]<< endl;
  cout << "delta phi width (approx) " << DeltaPhiWidthApprox[0]<< endl;
  cout << "minimo dell'intervallo di BC " <<ALowBinFit[0] << " " <<  fHistPhiDistr_Corr[1][1][0]->GetXaxis()->GetBinLowEdge(fHistPhiDistr_Corr[1][1][0]->FindBin(ALowBinFit[0]))<< endl;
  cout << "massimo dell'intervallo di BC " <<AUpBinFit[0] << " " <<  fHistPhiDistr_Corr[1][1][0]->GetXaxis()->GetBinUpEdge(fHistPhiDistr_Corr[1][1][0]->FindBin(AUpBinFit[0]))<< endl;
  cout << " bin width " <<  fHistPhiDistr_Corr[1][1][0]->GetXaxis()->GetBinWidth(1) << endl;
  cout << " minimo bin dell'intervallo di BC " << fHistPhiDistr_Corr[1][1][0]->FindBin(ALowBinFit[0])<< " massimo " <<  fHistPhiDistr_Corr[1][1][0]->FindBin(AUpBinFit[0])<< endl;
  for (Int_t m=0; m< nummolt+1; m++){
    //if (m==0) continue;
    // if (m>=3) continue;     //added  
    for (Int_t v=PtV0Min; v< numPtV0Max; v++){
      //      cout << " m " << m << " PtV0 " << v << " n spectrum final " << NSpectrumFinal[m][v][0] << "  " << NSpectrumFinal[m][v][1]<< " e dall isto " << 	  fHistSpectrum[m][0]->GetBinContent(v+1)<< endl;
    }
  }

  TCanvas *canvasFit[5];
  for (Int_t i=0; i < 5; i++){
    canvasFit[i]=new TCanvas(Form("canvasFit%i",i),Form("canvasFit%i",i), 1300, 1000);
    canvasFit[i]->Divide(3,2);
  }
  TH1F* histFit[7][numfittipo];
  for (Int_t i=0; i < 7; i++){
    for (Int_t typefit=0; typefit<numfittipo; typefit++){
      if (    nameFit[typefit]== "pT-scaling") continue;
      histFit[i][typefit]=new TH1F(Form("histFit%i_%i",i, typefit),Form("histFit%i_%i",i, typefit),nummolt, Nmolt );
      histFit[i][typefit]->SetLineColor(ColorFit[typefit]);
      histFit[i][typefit]->SetMarkerColor(ColorFit[typefit]);
    }
  }
    cout << "\n\n********************************* info about integrals (integral fraction below fixed pt, parameter values, chisquare. \n If interested decomment me! " << endl;
  for (Int_t m=0; m< nummolt+1; m++){
    //if (m==0) continue;
    Int_t syst=0;


    cout << " multiplcity " << m << endl;
    for (Int_t typefit=0; typefit<numfittipo; typefit++){
      if (    nameFit[typefit]== "pT-scaling") continue;
    /*    cout << "\n\n********************************* info about integrals (integral fraction below fixed pt, parameter values, chisquare. \n If interested decomment me! " << endl;
      cout << "\n" << nameFit[typefit]<< endl;   
      cout << "m " << m << "chisquare/NDF " << fit_MTscaling[m][syst][typefit]->GetChisquare() << "/" << fit_MTscaling[m][syst][typefit]->GetNDF() <<"\n" <<  endl;
      cout << " T " << fit_MTscaling[m][syst][typefit]->GetParameter(1)<< endl;
      if (typefit==numfittipo-1)    cout << " n " << fit_MTscaling[m][syst][typefit]->GetParameter(2)<< endl;
      if (    nameFit[typefit]== "Levi")    cout << " n " << fit_MTscaling[m][syst][typefit]->GetParameter(1)<< endl;
      cout << "integral percentage for pT < 0.5 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 0.5)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
      cout << "integral percentage for pT < 1.0 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 1)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
      cout << "integral percentage for pT < 1.5 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 1.5)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
      cout << "integral percentage for pT >8  GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(8, 30)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
      cout << "integral percentage for pT <30 (over 80)   GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 30)/fit_MTscaling[m][syst][typefit]->Integral(0, 80)<< endl;
      cout << "Integral up to 30 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
      cout << "Integral up to 0.5 GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 0.5)<< endl;
      cout << "integral percentage for pT >8  GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(8, 30)/fit_MTscaling[m][syst][typefit]->Integral(0, 30)<< endl;
      cout << "integral percentage for pT <30 (over 80)   GeV/c " << fit_MTscaling[m][syst][typefit]->Integral(0, 30)/fit_MTscaling[m][syst][typefit]->Integral(0, 80)<< endl;
      cout << "adjacent bin integral " <<       fHistSpectrum[m][syst]->GetBinContent(2) *       fHistSpectrum[m][syst]->GetBinWidth(1) << endl;
    */
      //fill histFit
      if( m==nummolt) continue;
      if (typefit!=numfittipo-1){
	histFit[0][typefit]->SetBinContent(m+1,fit_MTscaling[m][syst][typefit]->GetParameter(1));
	histFit[0][typefit]->SetBinError(m+1,fit_MTscaling[m][syst][typefit]->GetParError(1));
      }
      else {
	histFit[0][typefit]->SetBinContent(m+1,fit_MTscaling[m][syst][typefit]->GetParameter(2));
	histFit[0][typefit]->SetBinError(m+1,fit_MTscaling[m][syst][typefit]->GetParError(2));
      }
      //      histFit[0][typefit]->GetYaxis()->SetTitle("parameter T");
      histFit[0][typefit]->GetYaxis()->SetTitle(" ");
      histFit[0][typefit]->GetYaxis()->SetTitleOffset(0.05);
      histFit[0][typefit]->GetYaxis()->SetRangeUser(0,2.);

      histFit[1][typefit]->SetBinContent(m+1,YieldExtr[m][typefit]/fit_MTscaling[m][syst][typefit]->Integral(0, 50));
      histFit[1][typefit]->SetBinError(m+1,YieldExtrErrStat[m][typefit]/fit_MTscaling[m][syst][typefit]->Integral(0, 50));
      histFit[1][typefit]->GetYaxis()->SetTitle(""); //extrapolated fraction
      histFit[1][typefit]->GetYaxis()->SetRangeUser(0,0.3);

      histFit[2][typefit]->SetBinContent(m+1,fit_MTscaling[m][syst][typefit]->Integral(0, 0.5)/fit_MTscaling[m][syst][typefit]->Integral(0, 50));
      //      histFit[2][typefit]->GetYaxis()->SetTitle("integral fraction below 0.5 GeV/c");
      histFit[2][typefit]->GetYaxis()->SetTitle(" ");
      histFit[2][typefit]->GetYaxis()->SetTitleOffset(1.5);
      histFit[2][typefit]->GetYaxis()->SetRangeUser(0,1);

      histFit[3][typefit]->SetBinContent(m+1,fit_MTscaling[m][syst][typefit]->Integral(0, 1)/fit_MTscaling[m][syst][typefit]->Integral(0, 50));
      //      histFit[3][typefit]->GetYaxis()->SetTitle("integral fraction below 2 GeV/c");
      histFit[3][typefit]->GetYaxis()->SetTitle(" ");
      histFit[3][typefit]->GetYaxis()->SetTitleOffset(1.5);
      histFit[3][typefit]->GetYaxis()->SetRangeUser(0,1);

      histFit[4][typefit]->GetYaxis()->SetTitle("#Chi^{2}");
      histFit[4][typefit]->SetBinContent(m+1, fit_MTscaling[m][syst][typefit]->GetChisquare());

      histFit[5][typefit]->SetBinContent(m+1,fit_MTscaling[m][syst][typefit]->Integral(8, 30)/fit_MTscaling[m][syst][typefit]->Integral(0, 50));
      histFit[5][typefit]->GetYaxis()->SetTitle("fraction above 8 GeV/c");
      histFit[5][typefit]->GetYaxis()->SetRangeUser(0,0.01);

      histFit[6][typefit]->SetBinContent(m+1,fit_MTscaling[m][syst][typefit]->GetParameter(0));
      histFit[6][typefit]->SetBinError(m+1,fit_MTscaling[m][syst][typefit]->GetParError(0));
      histFit[6][typefit]->GetYaxis()->SetTitle("norm");

    }
  }

  for (Int_t i=0; i < 6; i++){
    canvasFit[0]->cd(i+1);
    for (Int_t typefit=0; typefit<numfittipo; typefit++){
      if (    nameFit[typefit]== "pT-scaling") continue;
      histFit[i][typefit]->Draw("same");
    }
    legendfit->Draw("");
  }


  fileout->WriteTObject( canvasFit[0]);
  fileout->Close();

  cout << "\n\n I've produced the file " << stringout << endl; 
  cout << "check if number of trigger particles is calculated in the correct way! Different histograms should be used if trigger particles are the highest pt charged particles or if trigger particles are all the charged particle with pT > pT min" << endl;


  cout <<  "number of topo syst " << numSysV0Global << " total number of syst (topo + inv mass + DeltaEta + DeltaPhi " << numSystGlobal <<  " n. of inv mass region + DeltaEta region syst " << numsystPhiGlobal  <<endl;

  Int_t typefitFixed=0;
  for (Int_t typefit=0; typefit<numfittipo; typefit++){
    if (nameFit[typefit]==FitFixed)    typefitFixed=typefit;
  }

  for (Int_t m=0; m<nummolt+1; m++){
    cout << m << endl;
    cout << "integral for 0 <pT < 8 GeV/c " << fit_MTscaling[m][0][typefitFixed]->Integral(0, 8) << endl;
    cout << " low range spectrum (for smaller pT the integral of fit function is used to extrapolate " << LowRangeSpectrumPart[m] << endl;
    cout << " low range fit " <<LowRange[m] << endl; 
    cout << " up range fit " <<UpRange[m] << endl; 
    cout << "fraction extrapolated " <<    fit_MTscaling[m][0][typefitFixed]->Integral(0,LowRangeSpectrumPart[m]) / fit_MTscaling[m][0][typefitFixed]->Integral(0, 8) << endl;
  }
  /*just a check
    for (Int_t m=0; m<nummolt+1; m++){
    cout << "\n\n m " << m << endl;
    for(Int_t v=0; v < numPtV0; v++){
    cout << "before; LowRange " << LowRange[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
    if (LowRange[m]<=NPtV0[v]) break;
    cout << "LowRange " << LowRange[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
    }
    cout <<   PtBinMin[m]<< endl;
    }
  */
  for (Int_t m=0; m<nummolt+1; m++){
    if (m< nummolt) {
      cout << "valori degli istogrammi, errori relativi " << endl;
      cout << "yield           "<<     fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "error stat rel  "<<     fHistYieldvsErrSoloStat->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "error stat rel yield of fit with "<< FitFixed<< " " <<   fHistYieldFitvsErrSoloStat->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      //      cout << YieldFitIntegralErrorTypeFit[m][typefitFixed] <<"  "  << YieldFitIntegralError[m] << "  " <<  fit_MTscaling[m][0][typefitFixed]->IntegralError(0,8, fFitResultPtr0[m][0][typefitFixed] ->GetParams(),(fFitResultPtr0[m][0][typefitFixed]->GetCovarianceMatrix()).GetMatrixArray())<< " " << fit_MTscaling[m][0][typefitFixed]->Integral(0,8)<< endl;
      //      cout << "rel error on extrapolated yield " <<  YieldExtrErrStat[m] << " " <<YieldExtr[m] << " " << YieldExtrErrStat[m]/YieldExtr[m]<< endl;
      cout << "error sist rel  "<<     fHistYieldvsErrSoloSist->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "error total rel "<<     fHistYieldvsErrTotal->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
      cout << "*********" << endl;
      // cout << "gli stessi errori ma assoluti" << endl;
      // cout <<  fHistYieldvsErrErrSoloStat->GetBinContent( fHistYieldvsErrSoloStat->FindBin(mult[m])) << endl;
      // cout <<  fHistYieldvsErrErrSoloSist->GetBinContent( fHistYieldvsErrSoloSist->FindBin(mult[m])) << endl;
      // cout <<  fHistYieldvsErrErrTotal->GetBinContent( fHistYieldvsErrTotal->FindBin(mult[m])) << endl;
    }
  }
  cout << numPtV0Max << endl;
  cout << NPtV0[numPtV0Max-1] << endl;
}


