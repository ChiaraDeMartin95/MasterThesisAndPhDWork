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

void StyleHisto(TH1D *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title){
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
  //
}

void PtSpectraNew( Int_t type=0,  Int_t TypeAnalysis=1,Int_t numsysV0index=500, Int_t numsysTriggerindex=99, Int_t sysPhi = 0, Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=1,   Int_t israp=0,TString year=""/*"1617_AOD234_hK0s"/*"1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"*//*"Run2DataRed_hXi"/"2016kehjl_hK0s"*/, TString yearDCAzTrigger = ""/*"1617_hK0s"*/,  TString Path1 =""/*"_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, Bool_t isEfficiency=1,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=0, Int_t MultBinning=1, Int_t PtBinning=0, Bool_t SysTuvaWay=0,  Bool_t isNewInputPath=1,  Bool_t isHM =1, Bool_t ispp5TeV=0, Bool_t isEffCorr=1){

  Bool_t isGenOnTheFly = 0;
  if (isMC) isGenOnTheFly = 1; 
  //isGenOnTheFly --> events were generated on the fly and only the kinematic part is saved; the multiplicity distribution in percentile classes is not abvailable, instead classes based on the number of particles in the V0 acceptance are used
  if (!isMC && isGenOnTheFly) return;
  if (isGenOnTheFly) { //these variabes have no meaning for the MCtruth analysis -- they are set to zero in order not to appear in output file name
    MultBinning = 0;
    isHM =0;
    isEffCorr=0;
    isNewInputPath = 1;
    isEfficiency = 0;
  }

  //isEffCorr = 1 --> for K0s analysis. Efficiency was found to have a bug, the correctd efficiency is applied if isEffCorr = 1
  if (type==8) isEffCorr=0;

  Float_t PtTrigMin2=0.15;
  TString yearLowPtTrig = "";
  if (ispp5TeV) MultBinning =3;
  Int_t sys=0;
  //if I wnat to reproduce preliminary results I hae to set isErrDPhiCorrectlyProp=0
  if (type!=0 && type!=8) {cout << "type==0 (K0s) or type==8 (Xi) " << endl; return;}
  if (TypeAnalysis>2) {cout << "sys errors not yet implemented for these regions " << endl; return;}

  if (TypeAnalysis==0 && sysPhi>2) return;
  if (TypeAnalysis==1 && sysPhi>3) return;
  if (TypeAnalysis==2 && sysPhi>0) return;

  if (type==0){
    PtBinning=1;
    isEffCorr =1; //fix of efficiency
    yearDCAzTrigger="1617_hK0s";
    year= "1617_AOD234_hK0s";
    isNewInputPath=1;
    numsysV0index=200;
    MultBinning=0;
    if (isHM) {
      year= "AllhK0sHM_RedNo16k";
      isNewInputPath=1;
      numsysV0index=200;
      MultBinning=1;
    }
    if (ispp5TeV) {
      year= "17pq_hK0s";
      isNewInputPath=1;
      numsysV0index=186;
      MultBinning=3;
    }
    if (isMC){
      MultBinning = 0;
      isHM =0;
      isEffCorr =0;
      year= "_PythiaRopes_Test1";
    }
  }
  else if (type==8){
    PtBinning=0;
    yearDCAzTrigger="Run2DataRed_MECorr_hXi";
    year = "161718Full_AOD234_hXi";
    yearLowPtTrig = "_161718Full_AOD234_hXi";
    isNewInputPath=1;
    numsysV0index=400; //was 500
    PtTrigMin2=0.15;
    MultBinning=0;
    if (isHM) {
      year= "161718_HM_hXi_WithFlat16k_No18p";
      yearLowPtTrig = "_161718_HM_hXi_WithFlat16k_No18p";
      isNewInputPath=1;
      //      numsysV0index=360;
      numsysV0index=390;
      MultBinning=1;
      PtTrigMin2=3;
    }
    if (ispp5TeV) {
      year= "17pq_hXi";
      yearLowPtTrig = "_17pq_pp5TeV_hXi_pttrig0.15";
      isNewInputPath=1;
      numsysV0index=400;
      MultBinning=3;
    }
    if (isMC){
      MultBinning = 0;
      isHM =0;
      year= "_PythiaRopes_Test1";
    }
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

  TString hhCorr[2]={"", "_hhCorr"};
  TFile *fileinbis;
  TFile *filein;

  TFile *fileInDefault;
  TFile *fileInDefaultOOJ;
  TFile *fileSignalSys;
  TFile *fileSysTuva;
  TFile *fileDCAzTrigger;
  TFile *fileinDeltaEtaSys;

  TString PathInOOJ;
  TString PathInDefault;
  TString PathInDeltaEtaSys;
  TString PathInSignalSys;
  TString PathInDCAzSys;
  TString PurityPath;

  TString file = year;
  if (PtBinning!=0)  file+=Form("_PtBinning%i", PtBinning);
  file+=Path1;
  TString PathInBis =  "FinalOutput/AnalysisResults" + year/* + Path1*/  + ".root";
  if (isMC && isEfficiency)  PathInBis =  "FinalOutput/AnalysisResults" + year  + "_MCEff" /*+ Path1*/ +".root";
  if (isMC && !isEfficiency) PathInBis =  "FinalOutput/AnalysisResults" + year  + "_MCTruth.root";
  if (ishhCorr && !isMC)  PathInBis =  "FinalOutput/AnalysisResults" + year  +"_hhCorr" +Path1 + ".root";
  if (ishhCorr && isMC)  PathInBis =  "FinalOutput/AnalysisResults" + year  + "_hhCorr_MCEff" + Path1 + ".root";
  cout << "\npath in (task output) "<< PathInBis << endl;
  fileinbis=new TFile(PathInBis,"");
  if (!fileinbis){cout << PathInBis << " does not exist " << endl; return;}

  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  const Int_t nummolt=10;
  Int_t nummoltMax =10;
  if (!isGenOnTheFly) nummoltMax = 5;
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
  Double_t Nmoltpp5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};
  Double_t NmoltGenOnTheFly[nummolt+1]={0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300};
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  TString SmoltGenOnTheFly[nummolt+1]={"0-30", "30-60", "60-90", "90-120", "120-150", "150-180", "180-210", "210-240", "240-270", "270-300", "0-300"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegendpp5TeV[nummolt+1]={"0-10 %", "10-100 %", "100-100 %", "100-100 %", "100-100 %", "0-100 %"};
  TString SmoltLegendGenOnTheFly[nummolt+1]={"0-30", "30-60", "60-90", "90-120", "120-150", "150-180", "180-210", "210-240", "240-270", "270-300", "0-300"};

  for (Int_t m=0; m<nummoltMax+1; m++){
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
    else if (MultBinning==3){
      Nmolt[m] = Nmoltpp5TeV[m];
      Smolt[m] = Smoltpp5TeV[m];
      SmoltLegend[m] = SmoltLegendpp5TeV[m];
    }
  }
  if (isHM){
    Nmolt[1] = 0.001;
    Nmolt[2] = 0.005;
    Nmolt[3] = 0.01;
    Nmolt[4] = 0.05;
    Nmolt[5] = 0.1;
    Smolt[0] = "0-0.001";
    Smolt[1] = "0.001-0.005";
    Smolt[2] = "0.005-0.01";
    Smolt[3] = "0.01-0.05";
    Smolt[4] = "0.05-0.1";
    Smolt[5] = "0-0.1";
    SmoltLegend[0] = "0-0.001 %";
    SmoltLegend[1] = "0.001-0.005 %";
    SmoltLegend[2] = "0.005-0.01 %";
    SmoltLegend[3] = "0.01-0.05 %";
    SmoltLegend[4] = "0.05-0.1 %";
    SmoltLegend[5] = "0-0.1 %";
    if (MultBinning==1){
      Nmolt[1] = 0;
      Nmolt[2] = 0;
      Smolt[0] = "0-0a";
      Smolt[1] = "0-0b";
      Smolt[2] = "0-0.01";
      SmoltLegend[2] = "0-0.01 %";
    }
  }

  if (isGenOnTheFly){
    for (Int_t m=0; m<nummoltMax+1; m++){
      Smolt[m] = SmoltGenOnTheFly[m];
      Nmolt[m] = NmoltGenOnTheFly[m];
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

  TString  SmoltDCAzTrigger[nummolt+1] = {""};
  for (Int_t m=0; m<nummoltMax+1; m++){
    SmoltDCAzTrigger[m] = Smolt[m];
    if (isHM)     SmoltDCAzTrigger[m] = "_all";
    if (ispp5TeV)     SmoltDCAzTrigger[m] = "_all";
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

  if (isMC && !isEfficiency){
    SPtV01[1]= "0-0.5";
    NPtV01[1]= 0.;
  }

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
  TString SFit[3]={"exponential", "power law", "retta"};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};
  Int_t inutile[numSyst]=     {0,1,    1, 0,  1,  0,  1,  1,  1,  0,  1,  1,  0};
  Int_t Marker[numSyst]={20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  Int_t MarkerBetter[numSyst]={1,22, 32, 1, 29, 1,   3,  34, 33, 1, 20, 21, 22};
  //  Int_t Color[numSyst]= {1,  2,  3,  4,  5,  6,  7,  7,  4, 10,  6,  1,  2};
  //  Int_t Color[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 628,868,1, 909, 1,801,418,860,  1, 881, 7, 1};
  Int_t ColorBetter[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};
  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Color[3] ={628, 418, 600};
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

  TLegend *legendError = new TLegend(0.6, 0.6, 0.9, 0.9);

  Float_t YSup=0.012;
  Float_t YInf=-0.001;
  if (type==0 && TypeAnalysis==0) YSup=0.008;
  if (type==0 && isHM){
    if (TypeAnalysis!=0) YSup = 0.02;
  }
  if (type==8){
    YSup=0.0006;
    YInf = -0.0001;
    if (isHM){
      if (TypeAnalysis!=0)  YSup=0.0015;
    }
  }

  TH1D* fHistSpectrum_master[nummolt+1];
  TH1D* fHistSpectrum_masterOnlyStat[nummolt+1];
  TH1D* fHistSpectrum_masterOnlyStatRel[nummolt+1];
  TH1D* fHistSpectrum_masterSystCorr[nummolt+1];
  TH1D* fHistSpectrum_masterSystUnCorr[nummolt+1];
  TH1D* fHistSpectrum[nummolt+1][numSyst];  

  TH1D* fHistPhiDistr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosist[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistUnCorr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistCorr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistRelError[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistRelErrorCorr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistRelErrorUnCorr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistRelErrSE[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistSE[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistSERelErr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistSETuva[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistDCAz[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistDCAzRelErr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistDeltaEta[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solostat[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_master[nummolt+1][numPtV0];
  TH1F* hRelErrorPurity[nummolt+1];

  TF1*  gauss[nummolt+1][numPtV0][numSyst];
  TF1*  gaussint[nummolt+1][numPtV0][numSyst];

  Double_t YieldPerErrore[nummolt+1][numSyst]={0};
  Double_t YieldDefault[nummolt+1][numSyst]= {0};
  
  TH1F*   fHistV0EfficiencyPtBins[nummolt];
  TH1F*   HistContV0PtBins[nummolt];

  TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Lambda", "Xi","Xi", "Omega", "Omega", "Xi", "Omega"};
  TString TaskName = "";
  if (isNewInputPath){
    //    if (isMC)  TaskName = "_MCTruth_PtTrigMin3.0_PtTrigMax15.0";
    if (isGenOnTheFly) TaskName = "_PtTrigMin3.0_PtTrigMax15.0";
    else TaskName = "_PtTrigMin3.0_PtTrigMax15.0";
  }
  cout << "TaskName " << TaskName << endl;
  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask"+dirinputtype[type]+TaskName);
  if (!dir) {cout << " input dir not found " << endl; return;}

  TString ContainerName ="";
  if (isNewInputPath) {
    if (year.Index("AOD234")!=-1) ContainerName= "_hK0s_Task_";
    else if (year.Index("AOD235")!=-1) ContainerName= "_hK0s_Task_RecoAndEfficiency";
    else  ContainerName= "_hK0s_Task_";
    if (type==8) {
      ContainerName= "_hXi_Task_Default";
      if (ispp5TeV)       ContainerName= "_hXi_Task_";
    }
    if (isGenOnTheFly) {
      if (type == 0) ContainerName = "_hK0s_Task_K0s";
      else ContainerName = "_hK0s_Task_Xi";
    }
  }
  cout << "Container name " << ContainerName << endl;
  TList *list = (TList*)dir->Get("MyOutputContainer"+ ContainerName);
  if (!list) {cout << " input list not found " << endl; return;}


  TString stringout;
  stringout = Dir+"/DATA"+year0+"/PtSpectraNew" +hhCorr[ishhCorr];
  stringout+="_"+ year;
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+=   Form("_SysPhi%i_PtMin%.1f_", sysPhi, PtTrigMin);
  stringout+= RegionType[TypeAnalysis];
  if (SysTuvaWay) stringout+= "_SysTuvaWay";
  //  stringout += "_Prova";
  if (MultBinning!=0) stringout += Form("_MultBinning%i", MultBinning);
  //  stringout += "_NewTopoSelSys";
  // stringout += "_Boh";
  TString stringoutpdf = stringout;
  //  stringout += "_OldDPhiRange";
  if (isEffCorr) stringout += "_EffCorr";
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  //********************************************************************* 
  //definition of integral regions in deltaPhi distribution
  //********************************************************************* 
  const Int_t numSysDPhi =4;
  Float_t ALowBin[numSysDPhi]={-1}; 
  Float_t AUpBin[numSysDPhi]={1};
  Float_t ALowBinFit[numSysDPhi]={-1};
  Float_t AUpBinFit[numSysDPhi]={1};
  Float_t DeltaPhiWidth[numSysDPhi]={0};
  Float_t DeltaPhiWidthApprox[numSysDPhi]={0};

  if (TypeAnalysis==0 || TypeAnalysis==10){
    if (!ishhCorr && (type==0 || type==8)){
      ALowBinFit[0]=	ALowBin[0]={-1.05};
      AUpBinFit[0]=	AUpBin[0]={1.05};

      ALowBinFit[1]=	ALowBin[1]={-0.8};
      AUpBinFit[1]=	AUpBin[1]={0.8};

      ALowBinFit[2]=	ALowBin[2]={-1.2};
      AUpBinFit[2]=	AUpBin[2]={1.2};
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
    if (type==0){
      //OLD      ALowBinFit[0]=	ALowBin[0]={1.1};
      ALowBinFit[0]=	ALowBin[0]={1}; //NEW
      AUpBinFit[0]=	AUpBin[0]={1.8};

      ALowBinFit[1]=	ALowBin[1]={1.1}; //was 1 
      AUpBinFit[1]=	AUpBin[1]={1.8};
    
      ALowBinFit[2]=	ALowBin[2]={1.1};
      AUpBinFit[2]=	AUpBin[2]={2};

      ALowBinFit[3]=	ALowBin[3]={1};
      AUpBinFit[3]=	AUpBin[3]={2};
    }
    else {
      ALowBinFit[0]=	ALowBin[0]={1};
      //   ALowBinFit[0]=	ALowBin[0]={0.85};
      AUpBinFit[0]=	AUpBin[0]={1.8};

      ALowBinFit[1]=	ALowBin[1]={1.1};
      AUpBinFit[1]=	AUpBin[1]={1.8};
    
      ALowBinFit[2]=	ALowBin[2]={1.1};
      AUpBinFit[2]=	AUpBin[2]={2};

      //ALowBinFit[0]=	ALowBin[0]={0.85};
      ALowBinFit[3]=	ALowBin[3]={1};
      AUpBinFit[3]=	AUpBin[3]={2};
      /*
      ALowBinFit[0]=	ALowBin[0]={1.1};
      AUpBinFit[0]=	AUpBin[0]={2};

      ALowBinFit[1]=	ALowBin[1]={1};
      AUpBinFit[1]=	AUpBin[1]={2};
    
      ALowBinFit[2]=	ALowBin[2]={1};
      AUpBinFit[2]=	AUpBin[2]={1.8};

      ALowBinFit[3]=	ALowBin[3]={1.1};
      AUpBinFit[3]=	AUpBin[3]={1.8};
      */
    }
    /*
      ALowBinFit[1]=	ALowBin[1]={2};
      AUpBinFit[1]=	AUpBin[1]={4.28};

      ALowBinFit[2]=	ALowBin[2]={1};
      AUpBinFit[2]=	AUpBin[2]={4.28};
    */
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

  Double_t LowRangeJet[nummolt+1]= {1, 1, 1, 1, 1, 1}; 
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

  if (type==0){
    if (year=="1617_AOD234_hK0s") {
      for (Int_t m =0; m<nummoltMax+1; m++){
	if (PtTrigMin==3){
	  LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.5;
	  UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=3;
	  //	  if (m==4) UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=3;
	  LowRangeBulk[m]= 0.5;
	  LowRangeAll[m]= 0.5;
	  UpRangeAll[m]= 2;
	  UpRangeBulk[m]= 2;
	  LowRangeASBFit[m]=      LowRangeASZYAM[m]= 0.1;
	  UpRangeASBFit[m]=       UpRangeASZYAM[m]= 2;
	}
      }
    }
    else if (isHM){
      for (Int_t m=0; m<nummoltMax+1; m++){
	LowRangeJet[m] = 0.5;
	//LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.8;
	//      if (m==2) LowRangeJet[m] = 0.8;
	UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
	LowRangeBulk[m]= 0.5;
	LowRangeAll[m]= 0.5;
	UpRangeAll[m]= 3; //2
	UpRangeBulk[m]= 3; //2
      }
    }
    else if (ispp5TeV){
      for (Int_t m =0; m<nummoltMax+1; m++){
        LowRangeJet[m] = 0.5;
        UpRangeJet[m] = 4;
	if (m==0) LowRangeJet[m] = 0.8;
        LowRangeBulk[m]= 0.1;
        LowRangeAll[m]= 0.1;
        UpRangeAll[m]= 2.0;
        UpRangeBulk[m]= 2;
      }
    }
  }

  else if (type==8){
    if (year == "161718Full_AOD234_hXi"){
      for (Int_t m=0; m<nummoltMax+1; m++){
	LowRangeJet[m] = 1;
	if (m==1 || m ==0 || m==4) LowRangeJet[m]=  1.5;
	UpRangeJet[m] =4;
	LowRangeBulk[m]= 0.5;
	LowRangeAll[m]= 0.5;
	UpRangeBulk[m]= 3;
	UpRangeAll[m]= 3;
      }
    }
    else if (year == "161718_HM_hXi_WithFlat16k_No18p"){
      for (Int_t m=0; m<nummoltMax+1; m++){
	if (MultBinning ==1) {
	  LowRangeJet[m] = 1.5; //1.5
	  //        LowRangeBulkWithFull[m]= 1;
	  LowRangeBulk[m] = 1.;
	  LowRangeAll[m] = 1.;
	  UpRangeJet[m] = 8;
	  //        UpRangeBulkWithFull[m] = 3;
	  UpRangeBulk[m] = 3;
	  UpRangeAll[m] = 3;
	  LowRangeJetZYAM[m]= 1.0;
	  UpRangeJetZYAM[m]= 8;
	}
      }
    }
    else if (ispp5TeV){
      for (Int_t m=0; m<nummoltMax+1; m++){
	LowRangeJet[m] = 1.5;
	//      if (m!=0 && TypeOOJSub==1)  LowRangeJet[m] = 1;
	LowRangeBulk[m] = 1.;
	LowRangeAll[m] = 1.;
	UpRangeJet[m] = 4;
	if (m==0) UpRangeBulk[m] = 4;
	else  UpRangeBulk[m] = 3;
	if (m==0)      UpRangeAll[m] = 4;
	else      UpRangeAll[m] = 3;
      }
    }
  }

  //end of low and upper ranges for the fit****************************
  //************************************************************

  Int_t PtBinMin[nummolt+1]={0};
  for (Int_t m=0; m<nummoltMax+1; m++){
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

    //defining region where bin counting is performed (LowRangeSpectrumPart < pt < UpRangeSpectrumPart)
    //LowRangeSpectrumPart:

    if  (TypeAnalysis==1 || TypeAnalysis==2){
      if (type>0 && !isMC)    LowRangeSpectrumPart[m] = 0.5;
      else if (type>0 && isMC && m!=5)    LowRangeSpectrumPart[m] = 1.;
      else     LowRangeSpectrumPart[m] = 0.1;
    }
    else  LowRangeSpectrumPart[m] = LowRange[m];

    UpRangeSpectrumPart[m] = 8.;

    if (year== "1617_AOD234_hK0s"){
      if (SkipAssoc)     if (m==4) UpRangeSpectrumPart[m] = 4; //not enough stat above                                           
    }
    else if (year== "161718Full_AOD234_hXi"){
      if (m==4 || m ==3) {
	if (SkipAssoc)      UpRangeSpectrumPart[m] = 4; //not enough stat above
        LowRangeSpectrumPart[m] = 1; //not enough stat below
      }
      if (m==0 || m ==1 || m ==4) {
        if (TypeAnalysis==0)    LowRangeSpectrumPart[m] = 1.5;
      }
      if (m==4){
	if (SkipAssoc){
	  if (TypeAnalysis==0 || TypeAnalysis==1) UpRangeSpectrumPart[m] = 4.;
	}
      }
    }
    else if (year == "161718_HM_hXi_WithFlat16k_No18p"){
      if (TypeAnalysis!=0) LowRangeSpectrumPart[m] = 1.;
    }
    else if (year == "17pq_hXi"){
      LowRangeSpectrumPart[m] = 1.;
      //      UpRangeSpectrumPart[m] = 4;
      if (TypeAnalysis==0)   LowRangeSpectrumPart[m] = 1.5;
    }
    for(Int_t v=0; v < numPtV0Max; v++){
      if (LowRangeSpectrumPart[m]<=NPtV0[v]) break;
      PtBinMin[m]++;
      //      cout << "LowRange " << LowRangeSpectrumPart[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
    }
  }

  Int_t sysV0=0;    
  if(isMC && isEfficiency) file = year + "_MCEff" +Path1;
  //  file+=Form("_PtBinning%i", PtBinning);
  //  file+= Path1;

  TString titleX=  "p_{T} (GeV/c)";
  TString titleY=  "1/#Delta #eta #Delta #phi 1/N_{trigg} dN/dp_{T}";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 

  Float_t  YieldPt[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrStat[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSist[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistCorr[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistUnCorr[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistSE[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistSETuva[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistDCAz[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistPurity[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistDeltaEta[nummolt+1][numPtV0] ={0};
  TH1D* fHistSpectrumStat[nummolt+1];
  TH1D* fHistSpectrumSist[nummolt+1];
  TH1D* fHistSpectrumSistSE[nummolt+1];
  TH1D* fHistSpectrumSistSETuva[nummolt+1];
  TH1D* fHistSpectrumSistDCAz[nummolt+1];
  TH1D* fHistSpectrumSistPurity[nummolt+1];
  TH1D* fHistSpectrumSistDeltaEta[nummolt+1];
  TH1D* fHistSpectrumStatRelError[nummolt+1];
  TH1D* fHistSpectrumSistRelError[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorSE[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorSETuva[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorDeltaEta[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorDCAz[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorPurity[nummolt+1];

  TF1 * lineat1 = new TF1("pol0","pol0", -1./2*TMath::Pi(), 3./2*TMath::Pi());
  lineat1->SetLineColor(kBlack);
  lineat1->FixParameter(0,0);

  TCanvas* canvasPlotProj[nummolt+1];
  TCanvas* canvasPtSpectra = new TCanvas ("canvasPtSpectra", "canvasPtSpectra", 1300, 800);
  if (isHM || MultBinning==3)   canvasPtSpectra->Divide(2,2);
  else  canvasPtSpectra->Divide(3,2);
  TCanvas* canvasPtSpectraRelError = new TCanvas ("canvasPtSpectraRelError", "canvasPtSpectraRelError", 1300, 800);
  if (isHM || MultBinning==3)   canvasPtSpectraRelError->Divide(2,2);
  else  canvasPtSpectraRelError->Divide(3,2);
  TCanvas* canvasPtSpectraRelErrorAll = new TCanvas ("canvasPtSpectraRelErrorAll", "canvasPtSpectraRelErrorAll", 1300, 800);
  if (isHM || MultBinning==3)   canvasPtSpectraRelErrorAll->Divide(2,2);
  else  canvasPtSpectraRelErrorAll->Divide(3,2);

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
    else   if (TypeAnalysis==1){
      if (isHM) LimSupError =0.04;
      else LimSupError =0.1;
    }
    else   if (TypeAnalysis==2) {
      LimSupError =0.04;
      if (isHM)  LimSupError =0.04;
    }
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSupError =0.4;
    else   if (TypeAnalysis==1) LimSupError =0.25;
    else   if (TypeAnalysis==2) LimSupError =0.1; //0.05
    if (isHM){
      if (TypeAnalysis==1) LimSupError =0.1;
    }
  }

  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //    cout << " m " << m << endl;
    canvasPlotProj[m] = new TCanvas (Form("canvasPlot_m%i", m), Form("canvasPlot_m%i", m), 1300, 800);
    if (type==8)    canvasPlotProj[m]-> Divide(4, 2);
    else     canvasPlotProj[m]-> Divide((float)numPtV0/2+1, 2);
    fHistSpectrumStat[m]=new TH1D ("fHistSpectrum_"+Smolt[m],"fHistSpectrumStat_"+Smolt[m], numPtV0Max, NPtV0) ;
    fHistSpectrumStat[m]->Sumw2();
    fHistSpectrumSist[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSist_"+Smolt[m]);
    fHistSpectrumSistSE[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistSE_"+Smolt[m]);
    fHistSpectrumSistSETuva[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistSETuva_"+Smolt[m]);
    fHistSpectrumSistDCAz[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDCAz_"+Smolt[m]);
    fHistSpectrumSistPurity[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistPurity_"+Smolt[m]);
    if (TypeAnalysis!=2)    fHistSpectrumSistDeltaEta[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDeltaEta_"+Smolt[m]);

    //PathInDefault definition
    if (isMC && !isEfficiency) {
      PathInDefault=Dir+"/DATA"+year0+"/histo/AngularCorrelation" + year;
      PathInDefault += "_MCTruth";
      if (PtBinning!=0) PathInDefault += Form("_PtBinning%i", PtBinning);
    }
    else {
      PathInDefault = Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file;
    }
    if(type>=0){
      PathInDefault +="_"+tipo[type];
      PathInDefault +=Srap[israp];
      PathInDefault +=SSkipAssoc[SkipAssoc];
    }
    if (type==0 && isHM) PathInDefault += "_isBkgParab";
    PathInDefault+= hhCorr[ishhCorr]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", 0, sysV0, sys, PtTrigMin)+"_Output";
    if (!(isMC && !isEfficiency)){
      if (type==0 || (type==8 && isHM))      PathInDefault += "_Sidebands";
      if (ispp5TeV && type==8) PathInDefault += "_IsMEFrom13TeV";
      PathInDefault += "_IsEtaEff";
    }
    if (MultBinning!=0) PathInDefault += Form("_MultBinning%i", MultBinning);
    if (isMC && !isEfficiency) PathInDefault += "_IsMEFromCorrectCentrality";
    if (type==0) PathInDefault +="_NewdEtaChoice";
    if (isEffCorr) PathInDefault += "_EffCorr";
    if (isMC && !isEfficiency) PathInDefault += "_MCPrediction";
    PathInDefault += ".root";

    cout << "\nDefault file: " << PathInDefault << endl;
    fileInDefault = new TFile(PathInDefault, "");

    //PathInOOj definition (for jet production of Xi)
    PathInOOJ= "OOJComparison"+year;
    PathInOOJ +=  yearLowPtTrig;
    if (PtBinning>0) PathInOOJ+=Form("_PtBinning%i",PtBinning);
    PathInOOJ+= Path1;
    //  PathInOOJ+="_Jet0.75";
    if(type>=0){
      PathInOOJ +="_"+tipo[type];
      PathInOOJ +=Srap[israp];
      PathInOOJ +=SSkipAssoc[SkipAssoc];
    }

    PathInOOJ+=Form("_sys0_PtTrigMin%.1f_PtTrigMin%.1f_Output", PtTrigMin, PtTrigMin2);
    PathInOOJ += "_IsEtaEff";
    if (ispp5TeV && type==8) PathInOOJ += "_isMEFrom13TeV";
    if (MultBinning!=0) PathInOOJ += Form("_MultBinning%i", MultBinning);
    if (isHM)     PathInOOJ += "_Sidebands";
    //if (isHM) PathInOOJ += "_isOOJFromAllMult";
    if (!isHM)    PathInOOJ += "_NewDEtaChoice";
    PathInOOJ += ".root";

    if (type==8 && TypeAnalysis==0){
    fileInDefaultOOJ = new TFile(PathInOOJ,"");
    }

    //PathInSignalSys definition
    TString PathInSysDCAz="";
    TString PathInSys=Dir+"/DATA2016/histo/SignalExtractionStudy";
    TString PathInSysTuva=Dir+"/DATA2016/LoosestTightestTopoSel";
    PathInSysDCAz=PathInSys;
    PathInSysDCAz+=yearDCAzTrigger ;
    PathInSys+=year ;
    PathInSysTuva+=year ;
    TString PathInSysDeltaEta=Dir+"/DATA2016/SysDeltaEta";
    PathInSysDeltaEta+=year ;
    if (PtBinning>0) {
      PathInSys +=Form("_PtBinning%i",PtBinning);
      PathInSysDCAz +=Form("_PtBinning%i",PtBinning);
      PathInSysTuva +=Form("_PtBinning%i",PtBinning);
      PathInSysDeltaEta +=Form("_PtBinning%i",PtBinning);
    }
    if(type>=0){
      if (!ishhCorr) {
	PathInSys +="_"+tipo[type];
	PathInSysDCAz +="_"+tipo[type];
	PathInSysTuva +="_"+tipo[type];
	PathInSysDeltaEta +="_"+tipo[type];
      }
      PathInSys +=Srap[israp];
      //PathInSys +=SSkipAssoc[SkipAssoc];
      PathInSysDCAz +=Srap[israp];
      //PathInSysDCAz +=SSkipAssoc[SkipAssoc];
      PathInSysTuva +=Srap[israp];
      PathInSysTuva +=SSkipAssoc[SkipAssoc];
      PathInSysDeltaEta +=Srap[israp];
      PathInSysDeltaEta +=SSkipAssoc[SkipAssoc];
    }

    PathInSignalSys= PathInSys+  Form("_SysT%i_SysV0Num%i_Sys%i_PtMin%.1f", 0, numsysV0index, sys, PtTrigMin);
    if (MultBinning!=0) PathInSignalSys += Form("_MultBinning%i", MultBinning);
    if (!ispp5TeV)   PathInSignalSys += "_Bis";
    PathInSignalSys += ".root";
    PathInSysTuva= PathInSysTuva+  Form("_PtMin%.1f_", PtTrigMin)+RegionTypeOld[TypeAnalysis]+".root" ;
    PathInDCAzSys= PathInSysDCAz +  Form("_SysTNum%i_SysV0%i_Sys%i_PtMin%.1f.root", numsysTriggerindex, 0, sys, PtTrigMin);

    //file from where I get the syst uncertainty on pt spectra associated to fit procedure
    PurityPath = "FinalOutput/DATA2016/SysPurity";
    PurityPath +=year;
    if (PtBinning>0)      PurityPath+="_PtBinning1";
    PurityPath +="_"+tipo[type];
    PurityPath +=Srap[israp];
    //    PurityPath +=SSkipAssoc[SkipAssoc];
    PurityPath +=Form("_PtMin%.1f",PtTrigMin);
    if (MultBinning!=0) PurityPath += Form("_MultBinning%i", MultBinning);
    PurityPath +=".root";

    cout <<"Syst. assoc to signal extraction from file: " <<  PathInSignalSys << endl;
    fileSignalSys = new TFile(PathInSignalSys, "");
    if (SysTuvaWay)    fileSysTuva = new TFile(PathInSysTuva, "");
    fileDCAzTrigger = new TFile(PathInDCAzSys, "");

    //    if (!filein) {cout << filein << "does not exist " << endl ; return; }
    if (!fileSignalSys) {cout << PathInSignalSys << "does not exist " << endl ; return; }
    if (!fileSysTuva && SysTuvaWay) {cout << PathInSysTuva << "does not exist " << endl ; return; }
    if (!fileDCAzTrigger) {cout << PathInDCAzSys << "does not exist " << endl ; return; }

    if(TypeAnalysis!=2){
      PathInDeltaEtaSys = PathInSysDeltaEta + hhCorr[ishhCorr]+Form("_PtMin%.1f_", PtTrigMin)+RegionType[TypeAnalysis];
      if (MultBinning!=0) PathInDeltaEtaSys += Form("_MultBinning%i", MultBinning);
      if ((type==0 || isHM) && !(isMC && !isEfficiency)) PathInDeltaEtaSys += "_Sidebands";
      if (type==0) {
	PathInDeltaEtaSys +="_NewdEtaChoice";
      }
      if (type==8 && TypeAnalysis==0) {
	if (!isHM) PathInDeltaEtaSys +="_NewDEtaChoice";
      }
      if (type==8 && TypeAnalysis==1 && ispp5TeV) {
	PathInDeltaEtaSys +="_NewDEtaChoice";
      }
      if (isEffCorr) PathInDeltaEtaSys += "_EffCorr";
      PathInDeltaEtaSys += ".root";
      cout <<"Syst. assoc to choice of DeltaEta: " <<  PathInSysDeltaEta << endl;
      fileinDeltaEtaSys = new TFile(PathInDeltaEtaSys, "");
      cout << "File with DeltaEta syst.: " << PathInDeltaEtaSys << endl;
      if (!fileinDeltaEtaSys) {cout << PathInDeltaEtaSys << "does not exist " << endl ; return; }
      //      cout << "\nSyst. associated to DeltaEta region from file " <<       PathInDeltaEtaSys<<endl;
    }

    //   cout <<     "I start the loop on pT " << endl;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      //      cout << v << endl;
      if (TypeAnalysis==0){
	if (type==8 && NPtV0[v]<2.5-0.001)  {
	  fHistPhiDistr_master[m][v]=(TH1D*)fileInDefaultOOJ->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmooth");
	}
	else fHistPhiDistr_master[m][v]=(TH1D*)fileInDefault->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_TrCorr");
      }
      if (TypeAnalysis==1) fHistPhiDistr_master[m][v]=(TH1D*)fileInDefault->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr_TrCorr");
      if (TypeAnalysis==2 || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_master[m][v]=(TH1D*)fileInDefault->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulk_EffCorr_TrCorr");
      if (TypeAnalysis==3) fHistPhiDistr_master[m][v]=(TH1D*)fileInDefault->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulk_EffCorrNotScaled_TrCorr");
      cout << "ME_m"<<Smolt[m]<<"_v"<<SPtV0[v]<<"_AC_phi_V0Sub_JetBulk_EffCorr_TrCorr" << endl;
      if(!fHistPhiDistr_master[m][v]) {cout << "histo master not present " << endl; return;}
      fHistPhiDistr_master[m][v]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]);

      //getting phi distributions with systematic errors from signal extraction
      if (!(isMC && !isEfficiency)) {
	if (TypeAnalysis==0) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_DefSys");
	if (TypeAnalysis==1) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr_DefSys");
	if (TypeAnalysis==2 || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr_DefSys");
	if (TypeAnalysis==3) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled_DefSys");
	if(!fHistPhiDistr_solosistSE[m][v]){ cout << " histo with signal extraction error not presetn " << endl; return;}
	if (SysTuvaWay){ 
	  fHistPhiDistr_solosistSETuva[m][v]=(TH1D*)fileSysTuva->Get("PhiDistr_solosistSETuvaWay_m"+Smolt[m]+"_v"+SPtV0[v]);
	  if(!fHistPhiDistr_solosistSETuva[m][v])  {cout << "histo SETuva not present " << endl; return;}
	} 
	
	if (TypeAnalysis==0) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+SmoltDCAzTrigger[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_DefSys");
	if (TypeAnalysis==1) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+SmoltDCAzTrigger[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr_DefSys");
	if (TypeAnalysis==2 || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+SmoltDCAzTrigger[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr_DefSys");
	if (TypeAnalysis==3) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+SmoltDCAzTrigger[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled_DefSys");

	if(!fHistPhiDistr_solosistDCAz[m][v]){ cout << " histo with DCAzTrigger error not present" << endl; return;}
      }
      else {
	cout << "***Define fake histograms for systematic studies for MC truth***" << endl;
	fHistPhiDistr_solosistSE[m][v]=(TH1D*)fHistPhiDistr_master[m][v]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_FakeSEForMCTruth");
	fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fHistPhiDistr_master[m][v]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_FakeDCAzForMCTruth");
      }

      if (kTRUE){ //SMOOTH
	fHistPhiDistr_solosistDCAzRelErr[m][v] = (TH1D*) fHistPhiDistr_solosistDCAz[m][v]->Clone(Form("fHistPhiDistr_solosistDCAzRelErr_m%i_v%i",m, v));
	for (Int_t dphi=1; dphi<= fHistPhiDistr_solosistDCAz[m][v]->GetNbinsX(); dphi++){
	  fHistPhiDistr_solosistDCAzRelErr[m][v]->SetBinContent(dphi, fHistPhiDistr_solosistDCAz[m][v]->GetBinError(dphi)/fHistPhiDistr_solosistDCAz[m][v]->GetBinContent(dphi));
	}
	fHistPhiDistr_solosistDCAzRelErr[m][v]->Smooth(3);

	fHistPhiDistr_solosistSERelErr[m][v] = (TH1D*) fHistPhiDistr_solosistSE[m][v]->Clone(Form("fHistPhiDistr_solosistSERelErr_m%i_v%i",m, v));
	for (Int_t dphi=1; dphi<= fHistPhiDistr_solosistSE[m][v]->GetNbinsX(); dphi++){
	  fHistPhiDistr_solosistSERelErr[m][v]->SetBinContent(dphi, fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi)/fHistPhiDistr_solosistSE[m][v]->GetBinContent(dphi));
	}
	if (TypeAnalysis==0){
	  fHistPhiDistr_solosistSERelErr[m][v]->Smooth(3);
	}
	for (Int_t dphi=1; dphi<= fHistPhiDistr_solosistSE[m][v]->GetNbinsX(); dphi++){
	  fHistPhiDistr_solosistSE[m][v]->SetBinContent(dphi, fHistPhiDistr_master[m][v]->GetBinContent(dphi));
	  fHistPhiDistr_solosistSE[m][v]->SetBinError(dphi, fHistPhiDistr_solosistSERelErr[m][v]->GetBinContent(dphi)*fHistPhiDistr_master[m][v]->GetBinContent(dphi));
	}
      }

      if (TypeAnalysis!=2){
	fHistPhiDistr_solosistDeltaEta[m][v]=(TH1D*)fileinDeltaEtaSys->Get("PhiDistr_solosistDeltaEta_m"+Smolt[m]+"_v"+SPtV0[v]);
	if(!fHistPhiDistr_solosistDeltaEta[m][v])  {cout << "histo DeltaEta not present " << endl; return;}
	fHistPhiDistr_solosistDeltaEta[m][v]->SetName("PhiDistrSysDEta_m"+Smolt[m]+"_v"+SPtV0[v]);
	//	if ((type==8 && TypeAnalysis==1 && !isHM) || (type==8 && TypeAnalysis==0 && ispp5TeV)) {
	if ((type==8 && TypeAnalysis==1 && !isHM)) {
	  for (Int_t b=1; b<= fHistPhiDistr_solosistDeltaEta[m][v]->GetNbinsX(); b++){
	    fHistPhiDistr_solosistDeltaEta[m][v]->SetBinError(b,0); //I checked and the variation is not Barlow significant in the 1 < dphi < 2 range
	  }
	}
      }

      //*************************************************************************
      //Systematic sources correlated in dphi
      //A) Invariant mass fit procedure
      //	cout << " Syst associated to fit procedure: " << endl;
      if (isMC && !isEfficiency) PurityPath = "FinalOutput/DATA2016/SysPurity161718Full_AOD234_hXi_Xi_Eta0.8_PtMin3.0_Try.root"; //take as input a fake histogram from another file and set everythin to zero
      cout << "Purity Path: " << PurityPath << endl;
      TFile *PurityFile = new TFile(PurityPath, "");
      if (!PurityFile) {cout << "Purity file not available " << endl;  return;}
      if (isMC && !isEfficiency) hRelErrorPurity[m] = (TH1F*) PurityFile->Get(Form("hRelError_m%i", 5));
      else  hRelErrorPurity[m] = (TH1F*) PurityFile->Get(Form("hRelError_m%i", m));
      if (!hRelErrorPurity[m]) {cout << "Input histogram for purity not present " << endl; return;}
      if (type==8){
	for (Int_t b=1; b<= hRelErrorPurity[m]->GetNbinsX(); b++){
	  hRelErrorPurity[m]->SetBinContent(b, 0.005);
	}
      }
      if (isMC && !isEfficiency){
	for (Int_t b=1; b<= hRelErrorPurity[m]->GetNbinsX(); b++){
	  hRelErrorPurity[m]->SetBinContent(b, 0);
	}
      }
      //*************************************************************************

      fHistPhiDistr[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solostat[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solostat_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosist[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosist_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistUnCorr[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistUnCorr_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistCorr[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistCorr_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistRelError[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistRelError_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistRelErrorCorr[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistRelErrorCorr_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistRelErrorUnCorr[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistRelErrorUnCorr_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistRelErrSE[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistRelErrSE_m"+Smolt[m]+"_v"+SPtV0[v]);

      for(Int_t i =0; i <numSysDPhi; i++){
	DeltaPhiWidth[i]=TMath::Abs(fHistPhiDistr_master[m][v]->GetXaxis()->GetBinLowEdge(fHistPhiDistr_master[m][v]->FindBin(ALowBinFit[sysPhi])) - fHistPhiDistr_master[m][v]->GetXaxis()->GetBinUpEdge(fHistPhiDistr_master[m][v]->FindBin(AUpBinFit[sysPhi])));
	DeltaPhiWidthApprox[i]=TMath::Abs(AUpBinFit[i]-ALowBinFit[i]);
	if (TypeAnalysis==0){
	  DeltaPhiWidth[i]=TMath::Abs(fHistPhiDistr_master[m][v]->GetXaxis()->GetBinLowEdge(fHistPhiDistr_master[m][v]->FindBin(ALowBinFit[0])) - fHistPhiDistr_master[m][v]->GetXaxis()->GetBinUpEdge(fHistPhiDistr_master[m][v]->FindBin(AUpBinFit[0])));
	}
	if (TypeAnalysis==8 || TypeAnalysis==9 || TypeAnalysis==10) DeltaPhiWidth[i]=1; //not scaled for AS and JetNS                      
      }
      //      cout << " v " << v << endl;
    } //end loop v
    //cout << " end loop " << endl;
  } //end loop m

  Float_t DCAzSystError =0;
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    cout << "\nMultiplicity " << m << endl;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      //      cout << "I've set errors " << endl;
      for(Int_t dphi=1; dphi<=fHistPhiDistr_solostat[m][v]->GetNbinsX(); dphi++){
	//if (isHM || ispp5TeV) DCAzSystError =	fHistPhiDistr_solosistDCAzRelErr[m][v]->GetBinContent(dphi)*fHistPhiDistr_solosist[m][v]->GetBinContent(dphi);
	if (kTRUE) DCAzSystError =	fHistPhiDistr_solosistDCAzRelErr[m][v]->GetBinContent(dphi)*fHistPhiDistr_solosist[m][v]->GetBinContent(dphi);
	else DCAzSystError = fHistPhiDistr_solosistDCAz[m][v]->GetBinError(dphi);

	if (TypeAnalysis==2) fHistPhiDistr_solosist[m][v]->SetBinError(dphi,    sqrt(pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)+ pow(DCAzSystError,2) ));
	else 	if (TypeAnalysis==1 || TypeAnalysis==0) fHistPhiDistr_solosist[m][v]->SetBinError(dphi,    sqrt(pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)+ pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi),2) + pow(DCAzSystError,2)));

	//****SYSTEMATIC UNCERTAINTIES CORRELATED and UNCORRELATED IN DPHI******
	//TOPO sel: correlated in dphi, DCAz trigger not correlated in dphi
	if (TypeAnalysis==2) {
	  fHistPhiDistr_solosistCorr[m][v]->SetBinError(dphi, sqrt(pow( hRelErrorPurity[m] ->GetBinContent( hRelErrorPurity[m] ->FindBin(NPtV0[v]+0.001)) *fHistPhiDistr_solosistCorr[m][v]->GetBinContent(dphi),2) + pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)));
	  fHistPhiDistr_solosistUnCorr[m][v]->SetBinError(dphi, DCAzSystError);
	}
	else { //how to treat DeltaEta selection?
	  //13 TeV HM: TOPO sel: correlated in dphi, DCAz trigger not correlated in dphi, choice of dEta uncorrelated in dphi for JET and OOJ
	  if (isHM){
	    fHistPhiDistr_solosistCorr[m][v]->SetBinError(dphi, sqrt(pow( hRelErrorPurity[m] ->GetBinContent( hRelErrorPurity[m] ->FindBin(NPtV0[v]+0.001)) *fHistPhiDistr_solosistCorr[m][v]->GetBinContent(dphi),2) + pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)));

	    fHistPhiDistr_solosistUnCorr[m][v]->SetBinError(dphi, sqrt(pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi),2) + pow(DCAzSystError,2)));
	  }
	  //5 TeV MB: TOPO sel: correlated in dphi, DCAz trigger not correlated in dphi, choice of dEta uncorrelated in dphi for JET and OOJ 
	  else if (ispp5TeV){ //for both Xi and K0s the DeltaEta variation is not correlated in dphi
	    fHistPhiDistr_solosistUnCorr[m][v]->SetBinError(dphi, sqrt(pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi),2) + pow(DCAzSystError,2)));
	    fHistPhiDistr_solosistCorr[m][v]->SetBinError(dphi, sqrt(pow( hRelErrorPurity[m] ->GetBinContent( hRelErrorPurity[m] ->FindBin(NPtV0[v]+0.001)) *fHistPhiDistr_solosistCorr[m][v]->GetBinContent(dphi),2) + pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)));
	  }
	  //13 TeV MB: TOPO sel: correlated in dphi, DCAz trigger not correlated in dphi, choice of dEta not correlated in dphi
	  else {
	    fHistPhiDistr_solosistCorr[m][v]->SetBinError(dphi, sqrt(pow( hRelErrorPurity[m] ->GetBinContent( hRelErrorPurity[m] ->FindBin(NPtV0[v]+0.001)) *fHistPhiDistr_solosistCorr[m][v]->GetBinContent(dphi),2) + pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)));

	    if (TypeAnalysis!=2) fHistPhiDistr_solosistUnCorr[m][v]->SetBinError(dphi, sqrt(pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi),2) + pow(DCAzSystError,2)));
	    else fHistPhiDistr_solosistUnCorr[m][v]->SetBinError(dphi, DCAzSystError);
	  }
	}
	//***********************************************************************
	fHistPhiDistr_solosistRelError[m][v]->SetBinContent( dphi,   TMath::Abs(  fHistPhiDistr_solosist[m][v]->GetBinError(dphi)/fHistPhiDistr_solosist[m][v]->GetBinContent(dphi)));
	fHistPhiDistr_solosistRelErrorCorr[m][v]->SetBinContent( dphi,   TMath::Abs(  fHistPhiDistr_solosistCorr[m][v]->GetBinError(dphi)/fHistPhiDistr_solosist[m][v]->GetBinContent(dphi)));
	fHistPhiDistr_solosistRelErrorUnCorr[m][v]->SetBinContent( dphi,   TMath::Abs(  fHistPhiDistr_solosistUnCorr[m][v]->GetBinError(dphi)/fHistPhiDistr_solosist[m][v]->GetBinContent(dphi)));
	fHistPhiDistr_solosistRelErrSE[m][v]->SetBinContent( dphi,   TMath::Abs(  fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi)/fHistPhiDistr_solosist[m][v]->GetBinContent(dphi)));
      }

      /*
      fHistPhiDistr_solosistRelError[m][v]->Smooth();
      fHistPhiDistr_solosistRelErrorCorr[m][v]->Smooth();
      fHistPhiDistr_solosistRelErrorUnCorr[m][v]->Smooth();
      fHistPhiDistr_solosistRelErrSE[m][v]->Smooth();
      for(Int_t dphi=1; dphi<=fHistPhiDistr_solostat[m][v]->GetNbinsX(); dphi++){
	if (!isHM && !ispp5TeV && type==8 && TypeAnalysis==0 && m==4 && v==PtV0Min+1){
	  cout << "\nnon smoothed " << 	  fHistPhiDistr_solosist[m][v]->GetBinError(dphi);
	  cout << " rel error " << 	  fHistPhiDistr_solosistRelError[m][v]->GetBinContent(dphi) << " content " << fHistPhiDistr_solosist[m][v]->GetBinContent(dphi);
	  fHistPhiDistr_solosist[m][v]->SetBinError(dphi,      fHistPhiDistr_solosistRelError[m][v]->GetBinContent(dphi)*      fHistPhiDistr_solosist[m][v]->GetBinContent(dphi));
	  fHistPhiDistr_solosistCorr[m][v]->SetBinError(dphi,      fHistPhiDistr_solosistRelErrorCorr[m][v]->GetBinContent(dphi)*      fHistPhiDistr_solosist[m][v]->GetBinContent(dphi));
	  fHistPhiDistr_solosistUnCorr[m][v]->SetBinError(dphi,      fHistPhiDistr_solosistRelErrorUnCorr[m][v]->GetBinContent(dphi)*      fHistPhiDistr_solosist[m][v]->GetBinContent(dphi));
	  fHistPhiDistr_solosistSE[m][v]->SetBinError(dphi,      fHistPhiDistr_solosistRelErrSE[m][v]->GetBinContent(dphi)*      fHistPhiDistr_solosist[m][v]->GetBinContent(dphi));
	  cout << " smoothed " << 	  fHistPhiDistr_solosist[m][v]->GetBinError(dphi)<< endl;
	}
      }
      */
      YieldPt[m][v]=0;
      YieldPtErrStat[m][v]=0;
      YieldPtErrSist[m][v]=0;  
      YieldPtErrSistCorr[m][v]=0;
      YieldPtErrSistUnCorr[m][v]=0;

      YieldPtErrSistSE[m][v]=0;
      YieldPtErrSistSETuva[m][v]=0;
      YieldPtErrSistDCAz[m][v]=0;
      YieldPtErrSistPurity[m][v]=0;
      YieldPtErrSistDeltaEta[m][v]=0;

      for(Int_t dphi=fHistPhiDistr_solostat[m][v]->FindBin(ALowBinFit[sysPhi]); dphi<=fHistPhiDistr_solostat[m][v]->FindBin(AUpBinFit[sysPhi]); dphi++){
	YieldPt[m][v]+= fHistPhiDistr_solostat[m][v]->GetBinContent(dphi);
	YieldPtErrStat[m][v]+= pow(fHistPhiDistr_solostat[m][v]->GetBinError(dphi),2);
	YieldPtErrSistUnCorr[m][v]+= pow(fHistPhiDistr_solosistUnCorr[m][v]->GetBinError(dphi),2); 
	YieldPtErrSistCorr[m][v]+= fHistPhiDistr_solosistCorr[m][v]->GetBinError(dphi);

	//	if (isHM || ispp5TeV)  DCAzSystError =	fHistPhiDistr_solosistDCAzRelErr[m][v]->GetBinContent(dphi)*fHistPhiDistr_solosist[m][v]->GetBinContent(dphi);
	if (kTRUE)  DCAzSystError =	fHistPhiDistr_solosistDCAzRelErr[m][v]->GetBinContent(dphi)*fHistPhiDistr_solosist[m][v]->GetBinContent(dphi);
	else DCAzSystError = fHistPhiDistr_solosistDCAz[m][v]->GetBinError(dphi)/fHistPhiDistr_solosistDCAz[m][v]->GetBinContent(dphi) * fHistPhiDistr_solosist[m][v]->GetBinContent(dphi);

	//****propagate error according to their (un)correlation in dphi***************
	YieldPtErrSistSE[m][v]+= fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi); //phi corr
	YieldPtErrSistDCAz[m][v]+= pow(DCAzSystError,2); //phi uncorr
	if (SysTuvaWay)	YieldPtErrSistSETuva[m][v]+= fHistPhiDistr_solosistSETuva[m][v]->GetBinError(dphi);
	else	YieldPtErrSistSETuva[m][v]+= fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi);

	//DeltaEta part
	if (TypeAnalysis!=2){
	if (isHM){
	  if (type==0){
	    YieldPtErrSistDeltaEta[m][v]+= pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi), 2); //phi uncorr
	  } else {
	  if (TypeAnalysis==0) YieldPtErrSistDeltaEta[m][v]+= pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi), 2); //phi uncorr
	  else if (TypeAnalysis==1) YieldPtErrSistDeltaEta[m][v]+= fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi); //phi corr
	  }
	}
	else if (ispp5TeV){
	  if (TypeAnalysis!=2)	  YieldPtErrSistDeltaEta[m][v]+= pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi), 2); //phi uncor
	}
	else {
	  if (TypeAnalysis==0) YieldPtErrSistDeltaEta[m][v]+= pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi), 2); //phi uncorr
	  if (TypeAnalysis==1 && type==0) YieldPtErrSistDeltaEta[m][v]+= pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi), 2); //phi uncorr
	  //not significant for OOJ
	}
	}
      } //end of dphi loop

      YieldPtErrStat[m][v]= sqrt(YieldPtErrStat[m][v]);
      YieldPtErrSistDCAz[m][v] = sqrt(YieldPtErrSistDCAz[m][v]); //phi uncorr

      if (TypeAnalysis!=2){
      if (isHM){
	if (type==0){
	  YieldPtErrSistDeltaEta[m][v]= sqrt(YieldPtErrSistDeltaEta[m][v]); //phi uncorr
	}
	else {
	if (TypeAnalysis==0)	YieldPtErrSistDeltaEta[m][v]= sqrt(YieldPtErrSistDeltaEta[m][v]); //phi uncorr
	}
      }
      else if (ispp5TeV){ //delta eta choice not barlow significant
	  if (TypeAnalysis!=2)	YieldPtErrSistDeltaEta[m][v]= sqrt(YieldPtErrSistDeltaEta[m][v]); //phi uncorr
      }
      else {
	if (TypeAnalysis==0) YieldPtErrSistDeltaEta[m][v]= sqrt(YieldPtErrSistDeltaEta[m][v]); //phi uncorr
	else if (TypeAnalysis==1 && type==0) YieldPtErrSistDeltaEta[m][v]= sqrt(YieldPtErrSistDeltaEta[m][v]); //phi uncorr
      }
      }
      YieldPtErrSistUnCorr[m][v]= sqrt(YieldPtErrSistUnCorr[m][v]);

      //Below I add systematic effect associated to purity (fit to bkg of invariant mass distribution)
      YieldPtErrSistPurity[m][v]= hRelErrorPurity[m]->GetBinContent(hRelErrorPurity[m]->FindBin(NPtV0[v]+0.0001))*YieldPt[m][v];

      YieldPtErrSistCorr[m][v]= sqrt(pow(YieldPtErrSistCorr[m][v],2) + pow(YieldPtErrSistPurity[m][v], 2));
      YieldPtErrSist[m][v]= sqrt(pow(YieldPtErrSistCorr[m][v], 2) + pow(YieldPtErrSistUnCorr[m][v], 2));

      cout << "\e[35m " << SmoltLegend[m] << "\e[36m " << NPtV0[v] << " < pt < "<<  NPtV0[v+1] << " GeV/c, \e[39m  Yield: " <<  YieldPt[m][v] << " Rel. uncertainties on yield: "<< endl; 
      cout << "Purity: " <<     YieldPtErrSistPurity[m][v]/YieldPt[m][v] << " DCAz: " <<       YieldPtErrSistDCAz[m][v]/YieldPt[m][v] << " DeltaEta: " <<       YieldPtErrSistDeltaEta[m][v]/YieldPt[m][v] << " Topo sel: " <<       YieldPtErrSistSE[m][v]/YieldPt[m][v] << " Tot: " <<       YieldPtErrSist[m][v]/YieldPt[m][v]<<endl;
      cout << "Relative uncorr: " << YieldPtErrSistUnCorr[m][v]/YieldPt[m][v]<< endl;
      cout << "Relative corr: " << YieldPtErrSistCorr[m][v]/YieldPt[m][v]<< endl;
      //      for(Int_t v=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]-0.001); v<fHistSpectrumStat[m]->GetNbinsX() ; v++ ){
      for(Int_t v=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]-0.001); v<fHistSpectrumStat[m]->FindBin(UpRangeSpectrumPart[m]-0.001); v++){
	fHistSpectrumStat[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumStat[m]->SetBinError(v+1, YieldPtErrStat[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumSist[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumSist[m]->SetBinError(v+1, YieldPtErrSist[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

	fHistSpectrumSistSE[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumSistSE[m]->SetBinError(v+1, YieldPtErrSistSE[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

	fHistSpectrumSistSETuva[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumSistSETuva[m]->SetBinError(v+1, YieldPtErrSistSETuva[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

	fHistSpectrumSistDCAz[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumSistDCAz[m]->SetBinError(v+1, YieldPtErrSistDCAz[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

	fHistSpectrumSistPurity[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumSistPurity[m]->SetBinError(v+1, YieldPtErrSistPurity[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

	if (TypeAnalysis!=2){
	  fHistSpectrumSistDeltaEta[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	  fHistSpectrumSistDeltaEta[m]->SetBinError(v+1, YieldPtErrSistDeltaEta[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	}

      }
      /*
	for(Int_t b=1;  b<fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); b++){
	fHistSpectrumStat[m]->SetBinContent(b,0);
	fHistSpectrumStat[m]->SetBinError(b,0);
	fHistSpectrumSist[m]->SetBinContent(b,0);
	fHistSpectrumSist[m]->SetBinError(b,0);
	fHistSpectrumSistSE[m]->SetBinContent(b,0);
	fHistSpectrumSistSE[m]->SetBinError(b,0);
	if (TypeAnalysis!=2){
	fHistSpectrumSistDeltaEta[m]->SetBinContent(b,0);
	fHistSpectrumSistDeltaEta[m]->SetBinError(b,0);
	}
	}
      */
      if (TypeAnalysis!=0) YInf=0;
      canvasPlotProj[m]->cd(v+1);
      fHistPhiDistr_solostat[m][v]->GetYaxis()->SetRangeUser(YInf, YSup);
      fHistPhiDistr_solosist[m][v]->GetYaxis()->SetRangeUser(YInf, YSup);
      fHistPhiDistr_solosistCorr[m][v]->GetYaxis()->SetRangeUser(YInf, YSup);
      fHistPhiDistr_solostat[m][v]->Draw("same");
      fHistPhiDistr_solosist[m][v]->SetFillStyle(0);
      fHistPhiDistr_solosist[m][v]->Draw("same e2");
      fHistPhiDistr_solosistCorr[m][v]->SetFillStyle(3001);
      fHistPhiDistr_solosistCorr[m][v]->SetFillColorAlpha(kGreen+3,1);
      fHistPhiDistr_solosistCorr[m][v]->Draw("same e2");
      if (TypeAnalysis==0) lineat1->Draw("same");
      
    } //end loop on pt v0
    //    fHistSpectrumStat[m]->Scale(1./NTrigger[m]);
    //    fHistSpectrumSist[m]->Scale(1./NTrigger[m]);
    fHistSpectrumStat[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSist[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistSE[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistSETuva[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistDCAz[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistPurity[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    if (TypeAnalysis!=2)    fHistSpectrumSistDeltaEta[m]->Scale(1./DeltaPhiWidth[sysPhi]);

    fHistSpectrumStatRelError[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumStatRelError_"+Smolt[m]);
    fHistSpectrumSistRelError[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelError_"+Smolt[m]);
    fHistSpectrumSistRelErrorSE[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorSE_"+Smolt[m]);
    fHistSpectrumSistRelErrorSETuva[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorSETuva_"+Smolt[m]);
    fHistSpectrumSistRelErrorDeltaEta[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorDeltaEta_"+Smolt[m]);
    fHistSpectrumSistRelErrorDCAz[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorDCAz_"+Smolt[m]);
    fHistSpectrumSistRelErrorPurity[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorPurity_"+Smolt[m]);

    //    for(Int_t  b=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
    for(Int_t  b=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); b<=fHistSpectrumStat[m]->FindBin(UpRangeSpectrumPart[m]-0.001) ; b++){
      if (fHistSpectrumStat[m]->GetBinContent(b)!=0){
	fHistSpectrumStatRelError[m]->SetBinContent(b,    fHistSpectrumStat[m]->GetBinError(b)/    fHistSpectrumStat[m] ->GetBinContent(b));
	fHistSpectrumSistRelError[m]->SetBinContent(b,    fHistSpectrumSist[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
	fHistSpectrumSistRelErrorSE[m]->SetBinContent(b,    fHistSpectrumSistSE[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
	fHistSpectrumSistRelErrorSETuva[m]->SetBinContent(b,    fHistSpectrumSistSETuva[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
	fHistSpectrumSistRelErrorDCAz[m]->SetBinContent(b,    fHistSpectrumSistDCAz[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
	fHistSpectrumSistRelErrorPurity[m]->SetBinContent(b,    fHistSpectrumSistPurity[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
	if (TypeAnalysis!=2)    fHistSpectrumSistRelErrorDeltaEta[m]->SetBinContent(b,    fHistSpectrumSistDeltaEta[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));

	fHistSpectrumStatRelError[m]->SetBinError(b,0);
	fHistSpectrumSistRelError[m]->SetBinError(b,0);
	fHistSpectrumSistRelErrorSE[m]->SetBinError(b,0);
	fHistSpectrumSistRelErrorSETuva[m]->SetBinError(b,0);
	fHistSpectrumSistRelErrorDCAz[m]->SetBinError(b,0);
	fHistSpectrumSistRelErrorPurity[m]->SetBinError(b,0);
	fHistSpectrumSistRelErrorDeltaEta[m]->SetBinError(b,0);
      }
      else {
	fHistSpectrumStatRelError[m]->SetBinContent(b,0);
	fHistSpectrumSistRelError[m]->SetBinContent(b,0);
	fHistSpectrumSistRelErrorSE[m]->SetBinContent(b,0);
	fHistSpectrumSistRelErrorSETuva[m]->SetBinContent(b,0);
	fHistSpectrumSistRelErrorDCAz[m]->SetBinContent(b,0);
	fHistSpectrumSistRelErrorPurity[m]->SetBinContent(b,0);
	fHistSpectrumSistRelErrorDeltaEta[m]->SetBinContent(b,0);
      }
    }

    fHistSpectrumSistRelErrorSE[m]->GetXaxis()->SetRangeUser(LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorSETuva[m]->GetXaxis()->SetRangeUser(LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorDCAz[m]->GetXaxis()->SetRangeUser(LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorPurity[m]->GetXaxis()->SetRangeUser(LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorDeltaEta[m]->GetXaxis()->SetRangeUser(LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);

    if (isHM && type==8 && TypeAnalysis==0){
      fHistSpectrumSistRelErrorSE[m]->Smooth(1, "R");
      fHistSpectrumSistRelErrorSETuva[m]->Smooth(1, "R");
      fHistSpectrumSistRelErrorDCAz[m]->Smooth(1, "R");
      fHistSpectrumSistRelErrorPurity[m]->Smooth(1, "R");
      fHistSpectrumSistRelErrorDeltaEta[m]->Smooth(1, "R");
    }
    else if (type==0 && TypeAnalysis!=2){
      fHistSpectrumSistRelErrorDeltaEta[m]->Smooth(1, "R");
    }
    else if (ispp5TeV && TypeAnalysis!=2){
      fHistSpectrumSistRelErrorDeltaEta[m]->Smooth(1, "R");
    }
    if (ispp5TeV && type==0){
      fHistSpectrumSistRelErrorSE[m]->Smooth(1, "R");
    }

    for(Int_t  v=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); v<=fHistSpectrumStat[m]->FindBin(UpRangeSpectrumPart[m]-0.001) ; v++){
      if (fHistSpectrumStat[m]->GetBinContent(v)==0){
	fHistSpectrumSistRelError[m]->SetBinContent(v,0);
	fHistSpectrumSistRelErrorSE[m]->SetBinContent(v,0);
	fHistSpectrumSistRelErrorSETuva[m]->SetBinContent(v,0);
	fHistSpectrumSistRelErrorDCAz[m]->SetBinContent(v,0);
	fHistSpectrumSistRelErrorPurity[m]->SetBinContent(v,0);
	fHistSpectrumSistRelErrorDeltaEta[m]->SetBinContent(v,0);
      }
      fHistSpectrumSistSE[m]->SetBinError(v, fHistSpectrumSistRelErrorSE[m]->GetBinContent(v)*fHistSpectrumSist[m]->GetBinContent(v));
      fHistSpectrumSistSETuva[m]->SetBinError(v,fHistSpectrumSistRelErrorSETuva[m]->GetBinContent(v)*fHistSpectrumSist[m]->GetBinContent(v));
      fHistSpectrumSistDCAz[m]->SetBinError(v, fHistSpectrumSistRelErrorDCAz[m]->GetBinContent(v)*fHistSpectrumSist[m]->GetBinContent(v));
      fHistSpectrumSistPurity[m]->SetBinError(v, fHistSpectrumSistRelErrorPurity[m]->GetBinContent(v)*fHistSpectrumSist[m]->GetBinContent(v));
      if (TypeAnalysis!=2){
	fHistSpectrumSistDeltaEta[m]->SetBinError(v, fHistSpectrumSistRelErrorDeltaEta[m]->GetBinContent(v)*fHistSpectrumSist[m]->GetBinContent(v));
	fHistSpectrumSist[m]->SetBinError(v, sqrt(pow(fHistSpectrumSistSE[m]->GetBinError(v),2) + pow(fHistSpectrumSistDCAz[m]->GetBinError(v),2)+ pow(fHistSpectrumSistPurity[m]->GetBinError(v),2)+ pow(fHistSpectrumSistDeltaEta[m]->GetBinError(v),2)));
      }
      else {
	fHistSpectrumSist[m]->SetBinError(v, sqrt(pow(fHistSpectrumSistSE[m]->GetBinError(v),2) + pow(fHistSpectrumSistDCAz[m]->GetBinError(v),2)+ pow(fHistSpectrumSistPurity[m]->GetBinError(v),2)));
      }
      fHistSpectrumSistRelError[m]->SetBinContent(v,  fHistSpectrumSist[m]->GetBinError(v)/fHistSpectrumSist[m]->GetBinContent(v));
    }

    if (isHM)     canvasPtSpectra->cd(m+1-2); 
    else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectra->cd(3); 
    else     canvasPtSpectra->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSist[m], 0, LimSup, Color[TypeAnalysis], 6, titleX, titleY, title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumStat[m], 0, LimSup, Color[TypeAnalysis], 6, titleX, titleY,  title+SmoltLegend[m]);
    fHistSpectrumSist[m]    ->SetFillStyle(0);
    fHistSpectrumSist[m]->Draw("same e2");
    fHistSpectrumStat[m]->Draw("same e");

    if (isHM)     canvasPtSpectraRelErrorAll->cd(m+1-2); 
    else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraRelErrorAll->cd(3); 
    else     canvasPtSpectraRelErrorAll->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumStatRelError[m], 0, LimSupError, Color[TypeAnalysis], 33, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelError[m], 0, LimSupError, Color[TypeAnalysis], 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorSE[m], 0, LimSupError, 807, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorSETuva[m], 0, LimSupError, 630, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorDCAz[m], 0, LimSupError, 867, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorPurity[m], 0, LimSupError, 909, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorDeltaEta[m], 0, LimSupError, 881, 27, titleX, titleYRel,  title+SmoltLegend[m]);

    if(m==nummoltMax)    legendError->AddEntry(fHistSpectrumStatRelError[m], "stat.", "pl");
    if(m==nummoltMax && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelError[m], "syst.", "pl");
    if(m==nummoltMax)    legendError->AddEntry(fHistSpectrumSistRelErrorSE[m], "syst. topol. sel. M1", "l");
    if(m==nummoltMax && SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorSETuva[m], "syst topol. sel. M2", "l");
    if(m==nummoltMax && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorDCAz[m], "syst. DCAzTrigger", "l");
    if(m==nummoltMax && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorPurity[m], "syst. purity", "l");
    if(m==nummoltMax && TypeAnalysis!=2 && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorDeltaEta[m], "syst. #Delta #eta region", "l");
    if (!SysTuvaWay)    fHistSpectrumSistRelError[m]->Draw("same p");
    fHistSpectrumSistRelErrorSE[m]->Draw("same");
    if (SysTuvaWay)    fHistSpectrumSistRelErrorSETuva[m]->Draw("same");
    if (!SysTuvaWay){
      fHistSpectrumSistRelErrorDCAz[m]->Draw("same");
      fHistSpectrumSistRelErrorPurity[m]->Draw("same");
      if (TypeAnalysis!=2)  fHistSpectrumSistRelErrorDeltaEta[m]->Draw("same");
    }
    fHistSpectrumStatRelError[m]->Draw("same p");
    legendError->Draw("");

    if (isHM)     canvasPtSpectraRelError->cd(m+1-2); 
    else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraRelError->cd(3); 
    else   canvasPtSpectraRelError->cd(m+1);
    gPad->SetLeftMargin(0.15);
    fHistSpectrumSistRelError[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");
  }//end loop on m

  if (TypeAnalysis==0) {
    cout << "\n\e[36m***** Significance of jet peak***** \e[39m" << endl;
    for (Int_t m=0; m<nummoltMax+1; m++){
      if (isHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      cout <<"\n\e[35mMoltiplicity: " << Smolt[m] << " % \e[39m" << endl;
      for(Int_t v=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); v<=fHistSpectrumStat[m]->FindBin(UpRangeSpectrumPart[m]-0.001); v++){
	cout << SPtV0[v-1] << " significance: " << fHistSpectrumStat[m]->GetBinContent(v)/sqrt(pow(fHistSpectrumSist[m]->GetBinError(v),2) + pow(fHistSpectrumStat[m]->GetBinError(v),2))<< endl;
      }
    }
  }

  //  cout << "\nGoing to write on file " << endl;   
  fileout->WriteTObject(canvasPtSpectra);
  canvasPtSpectra->SaveAs(stringoutpdf + "_PtSpectra.pdf");
  fileout->WriteTObject(canvasPtSpectraRelError);
  canvasPtSpectraRelError->SaveAs(stringoutpdf + "_PtSpectraRelError.pdf");
  fileout->WriteTObject(canvasPtSpectraRelErrorAll);
  canvasPtSpectraRelErrorAll->SaveAs(stringoutpdf + "_PtSpectraAllRelErrors.pdf");
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    fileout->WriteTObject(canvasPlotProj[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorSE[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorDeltaEta[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorDCAz[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorPurity[m]);
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fileout->WriteTObject(fHistPhiDistr_solostat[m][v]);
      fileout->WriteTObject(fHistPhiDistr_solosist[m][v]);
      fileout->WriteTObject(fHistPhiDistr_solosistUnCorr[m][v]);
      fileout->WriteTObject(fHistPhiDistr_solosistCorr[m][v]);
    }
    fileout->WriteTObject(      fHistSpectrumSist[m]);
    fileout->WriteTObject(      fHistSpectrumStat[m]);

  }

  TH1F * DeltaPhiLimit = new TH1F ("DeltaPhiLimit", "DeltaPhiLimit", 2,0,2);
  DeltaPhiLimit->SetBinContent(1, fHistPhiDistr_master[nummoltMax][1]->GetXaxis()->GetBinLowEdge( fHistPhiDistr_master[nummoltMax][1]->FindBin(ALowBin[sysPhi] )));
  DeltaPhiLimit->SetBinContent(2, fHistPhiDistr_master[nummoltMax][1]->GetXaxis()->GetBinUpEdge( fHistPhiDistr_master[nummoltMax][1]->FindBin(AUpBin[sysPhi] )));
  DeltaPhiLimit->GetXaxis()->SetBinLabel(1,"Low DeltaPhi");
  DeltaPhiLimit->GetXaxis()->SetBinLabel(2,"Up DeltaPhi");
  fileout->WriteTObject(  DeltaPhiLimit);
  fileout->Close();

  //  cout << "pt spectra value for m==5 " << endl;
  for(Int_t v=PtV0Min; v < numPtV0Max; v++){
    //    cout << "v: " << v << " " << fHistSpectrumStat[5]->GetBinContent(v+1);
  }

  cout << "\n------------------------------------------------"<< endl;
  cout << "DeltaPhi width of one histo bin " << fHistPhiDistr_master[nummoltMax][1]->GetBinWidth(1)<< endl;
  cout <<  "DPhi range: "<< ALowBin[sysPhi]<<" - "<<AUpBin[sysPhi] << ", effective range: "<<  fHistPhiDistr_master[nummoltMax][1]->GetXaxis()->GetBinLowEdge( fHistPhiDistr_master[nummoltMax][1]->FindBin(ALowBin[sysPhi] ))<<" - "<< fHistPhiDistr_master[nummoltMax][1]->GetXaxis()->GetBinUpEdge( fHistPhiDistr_master[nummoltMax][1]->FindBin(AUpBin[sysPhi] ))<<endl;
  cout << "DeltaPhiWidth " << DeltaPhiWidth[sysPhi] << endl;

  cout << "\n------------------------------------------------"<< endl;
  cout << "\e[35mdefault plots from file:\e[39m "<<endl;
  cout << PathInDefault << endl;
  if (TypeAnalysis==0 && type==8)  cout << PathInOOJ << endl;

  if (!(isMC && !isEfficiency)){
    cout << "\n***Summary of syst uncertainty:*** " << endl;
    cout << "Purity: input file is relative uncertainty on spectrum " << endl;
    cout << "Systematics on topological selections: relative uncertainty is computed and assigned to default dphi projections" << endl;
    cout << "Systematics on DCAz selections: relative uncertainty is computed and assigned to the default dphi projections" << endl;
    if (TypeAnalysis!=2) cout << "Systematics on DeltaEta choice: comparison of dphi projections is done in this macro" << endl;

    cout << "\nStarting from the files:\nFor syst. associated to topo sel:\n "  << PathInSignalSys << "\nFor syst. associated to DeltaEta region:\n " << PathInDeltaEtaSys << "\nFor syst. associated to DCAz trigger:\n " << PathInDCAzSys << endl;
    cout << "For syst. associated to fit procedure:\n " <<     PurityPath <<endl;
  }
  cout << "\nI have created the file " << stringout << "\n" << endl;

}

