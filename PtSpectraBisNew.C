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
#include <TLegendEntry.h>
#include <TFile.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include <TRandom.h>
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include </data/dataalice/cdemart/AliPhysicsChiara/AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C>
#include "Macros/constants.h"
//#include </data/dataalice/cdemart/ALICE_analysis_tutorial/YieldMean.C>

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

void PtSpectraBisNew(Int_t type=8,  Int_t TypeAnalysis=0, Bool_t isppHM =0,Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=0,   Int_t israp=0, TString year=""/*"1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"/"Run2DataRed_hXi"/*"2016kehjl_hK0s"*/, Bool_t isEfficiency=0,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=0, Int_t MultBinning=1, Int_t PtBinning=1, TString FitFixed="Boltzmann"/*"Fermi-Dirac"*/,  Bool_t TwoFitFunctions=0, Bool_t isNormCorrFullyComputed=1, Bool_t isMeanMacro=0, Bool_t ispp5TeV=0,  Bool_t isErrorAssumedPtCorr=1, Bool_t isFitForPlot=0, Bool_t isEffCorr=0, Bool_t isdNdEtaTriggered=1, Int_t sys=0, Int_t MaterialBudgetCorr=2, Int_t isOnlyPlottingRelError=0, Int_t isTopoSel = 0){

  //isTopoSel
  //=0 --> default one
  //=1 --> default one but starting from file named "default"
  //=2 --> isTightest
  //=3 --> isLoosest
  TString Decision = "";
  if (isTopoSel!=0) {
    cout << "Are you sure tou want to run with the non-default selections of topo variables? (type y if you are)" << endl;
    cin>> Decision;
    if (Decision!="y") return;
  }

  Bool_t isGenOnTheFly = 0;
  if (isMC) isGenOnTheFly = 1;
  Bool_t isWingsCorrectionApplied =0;
  Int_t MonashTune =0;
  //isGenOnTheFly --> events were generated on the fly and only the kinematic part is saved; the multiplicity distribution in percentile classes is not abvailable, instead classes based on the number of particles in the V0 acceptance are used
  if (!isMC && isGenOnTheFly) return;
  if (isGenOnTheFly) { //these variabes have no meaning for the MCtruth analysis -- they are set to zero in order not to appear in output file name
    ispp5TeV=0;
    MultBinning = 0;
    isppHM =0;
    isEffCorr=0;
    isEfficiency = 0;
    isNormCorrFullyComputed = 0; //normalisation factor not needed for MC truth
    MaterialBudgetCorr=0;
  }
  if (type==0 || isGenOnTheFly){
  cout <<"Do you want to analyse the files with the wings correction applied? Type 1 if you DO want" << endl;
  cin >> isWingsCorrectionApplied;
  }

  if (TypeAnalysis>3) {cout << "sys errors not yet implemented for these regions " << endl; return;}
  if (type!=0 && type!=8) {cout << "Macro not working for particles different from K0s (type=0) and Xi (type=8)" << endl; return;}

  //NB: "ChangesIncluded" labels the new file produced by choosing Sidebands for K0s and Xi in HM, !SkipAssoc, and larger dEta choice for K0s

  //files from where norm factor is taken
  TString SfileNormCorrFC ="";

  if (type==0){
    MultBinning=0;
    PtBinning=1;
    if (!isMC) year= "1617_AOD234_hK0s";
    else  year= "";
    SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0_FIXED.root";
    isEffCorr=1; //fix in efficiency for K0s --- introduced in MAY 2022
    MaterialBudgetCorr=1;
    if (isGenOnTheFly){
      MaterialBudgetCorr=0;
      MultBinning = 0;
      isppHM =0;
      isEffCorr =0;
      cout << "Should be analyse Ropes (=1) or MonashDefault (=2) or EPOSLHC (=3)? " << endl;
      cin >> MonashTune;
      if (MonashTune == 1)    year = "PythiaRopes";
      else if (MonashTune ==2)        year = "PythiaMonash";
      else if (MonashTune ==3)        year = "EPOSLHC_3BEvForhK0s";
      else {cout << "option not valid " << endl; return;}
    }
  }
  else if (type==8){
    isEffCorr=0;
    MultBinning=0;
    PtBinning=0;
    if (!isMC) year="161718Full_AOD234_hXi";
    else  year=  "";
    SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue161718_hXi_Hybrid_vs_161718_hXi_Xi_Eta0.8_PtMin3.0";
    if (TypeAnalysis==3) SfileNormCorrFC += "_isBulkBlue";
    SfileNormCorrFC += "_FIXED.root";
    MaterialBudgetCorr=2;
    if (isGenOnTheFly){
      MaterialBudgetCorr=0;
      MultBinning = 0;
      isppHM =0;
      cout << "Should be analyse Ropes (=1) or MonashDefault (=2) or EPOSLHC (=3)? " << endl;
      cin >> MonashTune;
      if (MonashTune == 1)    year = "PythiaRopes_IncreasedStatXi";
      else if (MonashTune ==2)        year = "PythiaMonash_IncreasedStatXi";
      else if (MonashTune ==3)        year = "EPOSLHC_7BEvForhXi";
      else {cout << "option not valid " << endl; return;}
    }
  }

  if (isppHM) {
    MultBinning=1;
    MaterialBudgetCorr=2;
    if (type==0){
      PtBinning=1;
      if (!isMC) year= "AllhK0sHM_RedNo16k";
      else  year= "";
      SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0_FIXED.root";
    }
    else if (type==8){
      PtBinning=0;
      if (!isMC) year="161718_HM_hXi_WithFlat16k_No18p";
      else  year=  "";
      SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue161718_hXi_Hybrid_vs_161718_hXi_Xi_Eta0.8_PtMin3.0";
      if (TypeAnalysis==3) SfileNormCorrFC += "_isBulkBlue";
      SfileNormCorrFC += "_FIXED.root";
    }
  }
  else if (ispp5TeV){
    MaterialBudgetCorr=2;
    MultBinning=3;
    if (type==0){
      PtBinning=1;
      if (!isMC) year= "17pq_hK0s";
      else  year= "";
      SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue17pq_pp5TeV_Hybrid_vs_17pq_pp5TeV_PtBinning1_K0s_Eta0.8_PtMin3.0_MultBinning3_FIXED.root";
      if (TypeAnalysis==1 || TypeAnalysis==0) SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0_FIXED.root";
      //use the 13 TeV one for jet and ooj
    }
    else if (type==8){
      PtBinning=0;
      if (!isMC) year="17pq_hXi";
      else  year=  "";
      SfileNormCorrFC = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue161718_hXi_Hybrid_vs_161718_hXi_Xi_Eta0.8_PtMin3.0_FIXED.root";
    }
  }

  TString RegionType[4] = {"Jet", "Bulk", "Inclusive", "Bulk"};
  TString RegionTypeOld[4] = {"Jet", "Bulk", "All", "Bulk"};
  TString RegionTypeNew[4] = {"Jet", "Bulk", "Full", "Bulk"};

  const   Int_t NumberTypeAnalysis=11;
  cout << "Here's the meaning of different values of TypeAnalysis:" << endl;
  cout << 0 << "in-jet production " << endl;
  cout << 1 << "out-of-jet production  (Delta Phi between 1 and 2)" << endl;
  cout << 2 << "inclusive production (from JetBulkEffCorr)" << endl;
  cout << 3 << "inclusive production not scaled by DeltaEta and DeltaPhi (for comparison with published data) (from JetBulkEffCorr)" << endl;

  gStyle->SetOptStat(0);

  TF1 * lineat2= new TF1 ("lineat2", "pol0", 0,8);
  lineat2->FixParameter(0,2);
  lineat2->SetLineColor(kBlack);
  lineat2->SetLineWidth(1);
  TF1 * lineatm2= new TF1 ("lineatm2", "pol0", 0,8);
  lineatm2->FixParameter(0,-2);
  lineatm2->SetLineColor(1);
  lineatm2->SetLineWidth(1);
  TF1 * lineat0= new TF1 ("lineat0", "pol0", 0,8);
  lineat0->FixParameter(0,0);
  lineat0->SetLineColor(kBlack);
  lineat0->SetLineWidth(1);
  TF1 * lineat1= new TF1 ("lineat1", "pol0", 0,8);
  lineat1->FixParameter(0,1);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineWidth(1);

  TString hhCorr[2]={"", "_hhCorr"};
  TFile *fileinbis;
  TFile *filein;

  TString PathIn1;
  TString file = year;
  if (PtBinning!=0)  file+=Form("_PtBinning%i", PtBinning);

  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  Int_t nummoltMax = nummolt;
  if (!isGenOnTheFly) nummoltMax = 5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9
  const Int_t numPtTrigger=1;
  const Int_t numsysDPhi=4;
  const Int_t numtipo=10;

  Float_t LowPtLimitForAvgPtFS[nummolt+1] = {0};

  Int_t PtV0Min = 0; //0 el
  if (type>0)   PtV0Min =1;
  if (type==0 && PtBinning!=0) PtV0Min=1;
  if (isMC) {
    PtV0Min =0;
    if (type==0 && PtBinning!=0) PtV0Min=1;
  }

  TString DataOrMC[2]={"Data", "MC"};
  Int_t jet=TypeAnalysis;

  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  Float_t SpectrumSup[numtipo][NumberTypeAnalysis]={{0.015,0.2, 0.2,0.8,0.2,0.2,0.015, 0.015, 0.03, 0.03, 0.03},{0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001, 0.002, 0.002, 0.002},{0.015,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001}};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100}; 

  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  Double_t Nmoltpp5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};
  Double_t NmoltHM[nummolt+1]={0, 0, 0, 0.01, 0.05, 0.1}; 
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  TString SmoltHM[nummolt+1]={"0-0a", "0-0b", "0-0.01", "0.01-0.05", "0.05-0.1", "0-0.1"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegendHM[nummolt+1]={"0-0a %", "0-0b %", "0-0.01 %", "0.01-0.05 %", "0.05-0.1 %", "0-0.1 %"};
  TString SmoltLegendpp5TeV[nummolt+1]={"0-10 %", "10-100 %", "100-100 %", "100-100 %", "100-100 %", "0-100 %"};

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
    if (isppHM){
      Nmolt[m] = NmoltHM[m];
      Smolt[m] = SmoltHM[m];
      SmoltLegend[m] = SmoltLegendHM[m];
    }
  }
  if (isGenOnTheFly){
    for (Int_t m=0; m<nummoltMax+1; m++){
      Smolt[m] = SmoltGenOnTheFly[m];
      Nmolt[m] = NmoltGenOnTheFly[m];
      SmoltLegend[m] = SmoltGenOnTheFly[m];
    }
  }

  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", "", ""};  
  //TString SPtV0[numPtV0]={"", "0.5-1", "0.5-1",  "1-1.5","1.5-2", "2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (isGenOnTheFly && type==8) SPtV0[0]={"0-0.5"};

  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100, 100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  NPtV0[1]=0.5;
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", "", ""};
  SNPtV0[1]={"0.5"};

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
  TString SSystJet[2]={"BC [-1.0, 1.0]", "BC [-1.2, 1.2]"};
  TString SSystBulk[3]={"BC [1.0, 2.0]", "BC [2.0, 4.28]", "BC [1.0, 4.28]"};
  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Color[4] ={628, 418, 600, 418};

  TString stringout;
  TString PathIn0;
  TString PathIn;
  PathIn0 = Dir+"/DATA"+year0+"/PtSpectra";
  //  if (isppHM || ispp5TeV) PathIn0 +="New"; 
  PathIn0 +="New"; 
  PathIn0 += hhCorr[0];
  PathIn0 +="_" + year;
  if (PtBinning>0) PathIn0 +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    PathIn0 +="_"+tipo[type];
    PathIn0 +=Srap[israp];
    PathIn0 +=SSkipAssoc[SkipAssoc];
  }

  stringout = Dir+"/DATA"+year0+"/";
  stringout += "PtSpectraBisNew" +hhCorr[0];
  if (isppHM)   stringout += "_pp13TeVHM";
  if (ispp5TeV)   stringout += "_pp5TeV";
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  //  if (year!="1617_hK0s" && year!="Run2DataRed_MECorr_hXi") stringout+= "_"+year;
  stringout+= "_"+year;
  if(type>=0){
    stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+=   Form("_PtMin%.1f_", PtTrigMin);
  if (TypeAnalysis==3)  stringout+= "BulkBlue";
  else   stringout+= RegionType[TypeAnalysis];
  //  stringout += "_Try";
  //stringout += "_LowPtExtr0.5";
  if (TwoFitFunctions) stringout += "_TwoFitFunctions";
  if (isNormCorrFullyComputed) stringout+="_isNormCorrFullyComputed";
  //stringout += "_Fixed";
  if (isMeanMacro)  stringout += "_YieldMeanMacro"; 
  if (isErrorAssumedPtCorr)  stringout += "_isErrorAssumedPtCorr"; 
  if (!isGenOnTheFly)  stringout += "_ChangesIncluded";
  if (isdNdEtaTriggered) stringout += "_isdNdEtaTriggered";
  //  stringout += "_NewTopoSelSys";
  //stringout += "_PicAN";
  //  stringout +="_Ter";
  //  stringout +="_Prova";
  //  if (type==8)  PathInFakeSBDef += "_SBMEFromPeak";
  if (isFitForPlot)  stringout += "_isFitForPlot";
  //  stringout += "_NewDPhiRangeBis";
  //  stringout += "_OldDPhiRange";
  if (MultBinning!=0) stringout += Form("_MultBinning%i", MultBinning);
  if (isEffCorr) stringout += "_EffCorr";
  if (isWingsCorrectionApplied) {
    stringout += "_isWingsCorrectionApplied";
    if (isGenOnTheFly){
      if (type==8) stringout += "New";
    }
    else stringout += "New";
  }
  if (sys!=0) stringout += Form("_dEtaSys%i", sys);
  if (MaterialBudgetCorr==1) stringout += "_MatBudgetCorr";
  else if (MaterialBudgetCorr==2) stringout += "_MatBudgetCorrFAST";
  if (isTopoSel!=0){
    if (isTopoSel==1) stringout += "_SysV0Default";
    else if (isTopoSel==2) stringout += "_SysV0Tightest";
    else if (isTopoSel==3) stringout += "_SysV0Loosest";
  }
  TString PathOutPictures = stringout;
  //  stringout += "_PForT";
  stringout += "_MultCorr";
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

  if (type==0){
    if (year=="1617_AOD234_hK0s") {    
      for (Int_t m=0; m<nummoltMax+1; m++){
	LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=3;
	LowRangeBulk[m]= 0.5; 
	LowRangeAll[m]= 0.5; 
	UpRangeAll[m]= 2.5; 
	UpRangeBulk[m]= 2.5; 
	if (m==4) {
	  LowRangeBulk[m]= 0.1; 
	  UpRangeBulk[m]= 2.0; 
	}
	LowPtLimitForAvgPtFS[m] = 0.5;
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

	LowPtLimitForAvgPtFS[m] = 1.5;
      }
    }
  }

  if (isppHM){
    for (Int_t m=0; m<nummoltMax+1; m++){
      if (type==0){
	LowRangeJet[m] = 0.5;    
	UpRangeJet[m] = 4;
	LowRangeBulk[m]= 0.5;
	LowRangeAll[m]= 0.5;
	if (MaterialBudgetCorr==0){
	  UpRangeAll[m]= 3.; 
	  UpRangeBulk[m]= 3.; 
	}
	else {
	  UpRangeAll[m]= 2.5; 
	  UpRangeBulk[m]= 2.5; 
	}
	LowPtLimitForAvgPtFS[m] = 0.8;
	if (isFitForPlot) {
	  UpRangeBulk[m] = 8;
	  UpRangeAll[m] = 8;
	  UpRangeJet[m] = 8;
	}
      }
      else if (type==8){
	LowRangeJet[m] = 1.5;
	//	LowRangeJet[m] = 2;
	//	if (m==2) LowRangeJet[m] = 2; //check
	LowRangeBulk[m] = 1.;
	LowRangeAll[m] = 1.;
	UpRangeJet[m] = 8;
	UpRangeBulk[m] = 3;
	UpRangeAll[m] = 3;
	LowPtLimitForAvgPtFS[m] = 2;
	if (isFitForPlot) {
	  UpRangeBulk[m] = 8;
	  UpRangeAll[m] = 8;
	  UpRangeJet[m] = 8;
	}
      }
    }
  }
  else if (ispp5TeV){
    if (type==0){
      for (Int_t m =0; m<nummoltMax+1; m++){
        LowRangeJet[m] = 0.5;
        UpRangeJet[m] = 4;
	if (m==0) LowRangeJet[m] = 0.8;
        LowRangeBulk[m]= 0.1;
        LowRangeAll[m]= 0.1;
        UpRangeAll[m]= 2.0;
        UpRangeBulk[m]= 2;
	if (isFitForPlot) {
	  UpRangeBulk[m] = 8;
	  UpRangeAll[m] = 8;
	  UpRangeJet[m] = 8;
	}
      }
    }
    else if (type==8){
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
	if (isFitForPlot) {
	  UpRangeBulk[m] = 8;
	  UpRangeAll[m] = 8;
	  UpRangeJet[m] = 8;
	}
      }
    }
  }

  if (isMC && !isEfficiency){ //MCPrediction
    for (Int_t m=0; m<nummoltMax+1; m++){
      LowRangeAll[m] = 3.; //we are only interested in high-pt extrapolation (which is negligible however) 
      LowRangeBulk[m] = 3.;
      if (type==0){
	LowRangeAll[m] = 1.;
	LowRangeBulk[m] = 1.;
	LowRangeJet[m] = 0.5;
	if (TypeAnalysis==0)	LowPtLimitForAvgPtFS[m] = 0.5;
      }
      else {
	/*	
	if (m==nummoltMax-1 || m==nummoltMax-2 || m==nummoltMax-3) LowRangeJet[m] = 1.5;
	else    LowRangeJet[m] = 1.;
	*/
	if (MonashTune==1){
	  if (m==nummoltMax-1 || m == nummoltMax-2)  LowRangeJet[m] = 1.5;
          else LowRangeJet[m] = 1.;
	  //	  LowRangeJet[m] = 1.5;
	}
	else if (MonashTune==2){
	  if (m==nummoltMax-1) LowRangeJet[m] = 2;
	  else if (m==nummoltMax-4 || m==nummoltMax-3 || m==nummoltMax-2) LowRangeJet[m] = 1.5;
	  else LowRangeJet[m] = 1.5;
	}
	else if (MonashTune==3){
	  if (m==nummoltMax-1) LowRangeJet[m] = 1.5;
          //          else if (m==nummoltMax-4 || m==nummoltMax-3 || m==nummoltMax-2) LowRangeJet[m] = 1.5;
          else LowRangeJet[m] = 1.;
	}
	if (TypeAnalysis==0)	LowPtLimitForAvgPtFS[m] = LowRangeJet[m];
      }
      if (type==0)       UpRangeJet[m] = 3;
      else     {
	//UpRangeJet[m] = 4;
	if (m==nummoltMax-1 || m==nummoltMax-2 || m==nummoltMax-3 || m==nummoltMax-4) UpRangeJet[m] = 8;
	else    UpRangeJet[m] = 4;	
      }
      UpRangeAll[m] = 8;
      UpRangeBulk[m] = 8;
      cout <<  LowRangeJet[m] << "- " <<  UpRangeJet[m] << endl;
    }
  }
  Float_t LowPtLimitForAvgPtFSAllMult = LowPtLimitForAvgPtFS[nummolt];
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
      if (type>0)    LowRangeSpectrumPart[m] = 0.5;
      else LowRangeSpectrumPart[m] = 0.1;
      if (isMC & !isEfficiency) {
        LowRangeSpectrumPart[m] = 0;
      }
    }
    else  {
      LowRangeSpectrumPart[m] = LowRange[m];
    }

    if (year== "161718Full_AOD234_hXi"){
      if (m==4 || m ==3)  LowRangeSpectrumPart[m] = 1; //not enough stat below
      if (m==0 || m ==1 || m ==4) {
        if (TypeAnalysis==0)    LowRangeSpectrumPart[m] = 1.5;
      }
    }
    else if (year == "161718_HM_hXi_WithFlat16k_No18p"){
      if (TypeAnalysis!=0) LowRangeSpectrumPart[m] = 1.;
    }
    else if (year == "17pq_hXi"){
      LowRangeSpectrumPart[m] = 1.;
      if (TypeAnalysis==0)   LowRangeSpectrumPart[m] = 1.5;
    }

    for(Int_t v=0; v < numPtV0Max; v++){
      if (LowRangeSpectrumPart[m]<=NPtV0[v]) break;
      PtBinMin[m]++;
      //      cout << "LowRange " << LowRangeSpectrumPart[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
      //      cout << "PtBinMin " << PtBinMin[m] << endl;
    }

    if (isGenOnTheFly && type==0 && PtBinMin[m]==0) PtBinMin[m] = 1; //to skip first bin, which is 0-0

    //UpRangeSpectrumPart:
    UpRangeSpectrumPart[m] = 8;
    if (year== "161718Full_AOD234_hXi"){
      //if (type==8 && m==4 && (TypeAnalysis==0 || TypeAnalysis==1))    UpRangeSpectrumPart[m] = 4;
      if (m==4 || m ==3) {
        //UpRangeSpectrumPart[m] = 4; //not enough stat above
      }
    }
    else if (year == "17pq_hXi"){
      if (SkipAssoc) UpRangeSpectrumPart[m] = 4;
      else UpRangeSpectrumPart[m] = 8;
    }
  }

  for (Int_t m=0; m<nummoltMax+1; m++){
    cout << "\n\e[32m*******Multiplicity: " << SmoltLegend[m] << " *************\e[39m" << endl;
    cout << "fit range " << LowRange[m] << "- " << UpRange[m] << endl;
    cout << "!extrapolation range " << LowRangeSpectrumPart[m] << "- " << UpRangeSpectrumPart[m] << endl;
  }

  if(isMC && isEfficiency) file = year + "_MCEff";
  //  file+=Form("_PtBinning%i", PtBinning);

  TString titleX=  "#it{p}_{T} (GeV/#it{c})";
  TString titleY=  "1/#Delta#eta #Delta#phi 1/N_{trigg} dN/dp_{T}";
  titleY = "";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 

  TLegend *legendPhi = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legendError = new TLegend(0.6, 0.7, 0.9, 0.9);
  TLegend *legendError2 = new TLegend(0.6, 0.7, 0.9, 0.9);

  TH1F* fHistSpectrumStat[nummolt+1];
  TH1F* fHistSpectrumStatNotNorm[nummolt+1];
  TH1F* fHistSpectrumStatUp[nummolt+1];
  TH1F* fHistSpectrumStatDown[nummolt+1];
  TH1F* fHistSpectrumStatHard[nummolt+1];
  TH1F* fHistSpectrumStatSoft[nummolt+1];
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
  TH1F* fHistSpectrumSistMultUnCorr[nummolt+1];
  TH1F* fHistSpectrumSistLeadTrack[nummolt+1];
  TH1F* fHistSpectrumSistMCChoice[nummolt+1];
  TH1F* fHistSpectrumSistFakeSB[nummolt+1];
  TH1F* fHistSpectrumSistOOJSubDef[nummolt+1];

  //histos for relative uncertainty
  TH1F* fHistSpectrumStatRelError[nummolt+1]; //stat
  TH1F* fHistSpectrumSistRelError[nummolt+1]; //sist on the DeltaPhi projections
  TH1F* fHistSpectrumSistRelErrorDPhi[nummolt+1]; //sist assoc to choice of DeltaPhi interval
  TH1F* fHistSpectrumSistRelErrorAll[nummolt+1]; //total sist
  TH1F* fHistSpectrumSistRelErrorMultUnCorr[nummolt+1]; //sist uncorrelated in multiplicity
  TH1F* fHistSpectrumSistRelErrorDCAz[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorPurity[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorSE[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorDeltaEta[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMB[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMBCorrection[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorOOBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorIBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorLeadTrack[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMCChoice[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorFakeSB[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorOOJSubDef[nummolt+1]; 

  TCanvas* canvasYield = new TCanvas ("canvasYield", "canvasYield", 1300, 800);
  TCanvas* canvasYieldNormCorr = new TCanvas ("canvasYieldNormCorr", "canvasYieldNormCorr", 1300, 800);
  TCanvas* canvasYieldNormCorrRatio = new TCanvas ("canvasYieldNormCorrRatio", "canvasYieldNormCorrRatio", 1300, 800);
  TCanvas* canvasPtvsMult = new TCanvas ("canvasPtvsMult", "canvasPtvsMult", 1300, 800);
  TCanvas* canvasYieldErr = new TCanvas ("canvasYieldErr", "canvasYieldErr", 1300, 800);
  TCanvas* canvasYieldSistFracMultUnCorr = new TCanvas ("canvasYieldSistFracMultUnCorr", "canvasYieldSistFracMultUnCorr", 1300, 800);
  TCanvas* canvasPtErr = new TCanvas ("canvasPtErr", "canvasPtErr", 1300, 800);
  TCanvas* canvasNormFactor = new TCanvas ("canvasNormFactor", "canvasNormFactor", 1300, 800);
  if (isppHM || MultBinning==3)   canvasNormFactor->Divide(2,2);
  else   canvasNormFactor->Divide(3,2);
  TCanvas* canvasPtSpectra = new TCanvas ("canvasPtSpectra", "canvasPtSpectra", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectra->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectra->Divide(2,2);
  else   canvasPtSpectra->Divide(3,2);
  TCanvas* canvaspol0 = new TCanvas ("canvaspol0", "canvaspol0", 1300, 800);
  if (isppHM || MultBinning==3)   canvaspol0->Divide(2,2);
  else canvaspol0->Divide(3,2);
  TCanvas* canvasOOJSubDef = new TCanvas ("canvasOOJSubDef", "canvasOOJSubDef", 1300, 800);
  canvasOOJSubDef->Divide(3,2);
  TCanvas* canvasLeadTrackDef = new TCanvas ("canvasLeadTrackDef", "canvasLeadTrackDef", 1300, 800);
  canvasLeadTrackDef->Divide(3,2);
  TCanvas* canvasMCChoiceDef = new TCanvas ("canvasMCChoiceDef", "canvasMCChoiceDef", 1300, 800);
  canvasMCChoiceDef->Divide(3,2);
  TCanvas* canvasFakeSBDef = new TCanvas ("canvasFakeSBDef", "canvasFakeSBDef", 1300, 800);
  canvasFakeSBDef->Divide(3,2);

  TCanvas* canvasNormCorr = new TCanvas ("canvasNormCorr", "canvasNormCorr", 1300, 800);
  canvasNormCorr->Divide(3,2);
  TCanvas* canvasNormCorrRatio = new TCanvas ("canvasNormCorrRatio", "canvasNormCorrRatio", 1300, 800);
  canvasNormCorrRatio->Divide(3,2);

  TCanvas* canvasPtSpectraFit = new TCanvas ("canvasPtSpectraFit", "canvasPtSpectraFit", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraFit->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraFit->Divide(2,2);
  else   canvasPtSpectraFit->Divide(3,2);
  TCanvas* canvasPtSpectraFitUp = new TCanvas ("canvasPtSpectraFitUp", "canvasPtSpectraFitUp", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraFitUp->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraFitUp->Divide(2,2);
  else  canvasPtSpectraFitUp->Divide(3,2);
  TCanvas* canvasPtSpectraFitDown = new TCanvas ("canvasPtSpectraFitDown", "canvasPtSpectraFitDown", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraFitDown->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraFitDown->Divide(2,2);
  else  canvasPtSpectraFitDown->Divide(3,2);
  TCanvas* canvasPtSpectraFitHard = new TCanvas ("canvasPtSpectraFitHard", "canvasPtSpectraFitHard", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraFitHard->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraFitHard->Divide(2,2);
  else  canvasPtSpectraFitHard->Divide(3,2);
  TCanvas* canvasPtSpectraFitSoft = new TCanvas ("canvasPtSpectraFitSoft", "canvasPtSpectraFitSoft", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraFitSoft->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraFitSoft->Divide(2,2);
  else  canvasPtSpectraFitSoft->Divide(3,2);
  TCanvas* canvasDummy = new TCanvas ("canvasDummy", "canvasDummy", 1300, 800);
  TCanvas* canvasPtSpectraFitRatio = new TCanvas ("canvasPtSpectraFitRatio", "canvasPtSpectraFitRatio", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraFitRatio->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraFitRatio->Divide(2,2);
  else  canvasPtSpectraFitRatio->Divide(3,2);

  TCanvas *canvasExtrFractionLowPt = new TCanvas("canvasExtrFractionLowPt", "canvasExtrFractionLowPt", 1300, 800);
  TCanvas *canvasExtrFractionHighPt = new TCanvas("canvasExtrFractionHighPt", "canvasExtrFractionHighPt", 1300, 800);

  TCanvas* canvasPtSpectraFitBis = new TCanvas ("canvasPtSpectraFitBis", "canvasPtSpectraFitBis", 1300, 800);
  canvasPtSpectraFitBis->Divide(3,2);

  TCanvas* canvasPtSpectraAll = new TCanvas ("canvasPtSpectraAll", "canvasPtSpectraAll", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraAll->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraAll->Divide(2,2);
  else   canvasPtSpectraAll->Divide(3,2);
  TCanvas* canvasBarlow = new TCanvas ("canvasBarlow", "canvasBarlow", 1300, 800);
  if (isGenOnTheFly) canvasBarlow->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasBarlow->Divide(2,2);
  else  canvasBarlow->Divide(3,2);
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

  TCanvas* canvasRatioFakeSBDef = new TCanvas ("canvasRatioFakeSBDef", "canvasRatioFakeSBDef", 1300, 800);
  canvasRatioFakeSBDef->Divide(3,2);

  TCanvas* canvasPtSpectraRelError = new TCanvas ("canvasPtSpectraRelError", "canvasPtSpectraRelError", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraRelError->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraRelError->Divide(2,2);  
  else  canvasPtSpectraRelError->Divide(3,2);
  TCanvas* canvasPtSpectraRelErrorAll = new TCanvas ("canvasPtSpectraRelErrorAll", "canvasPtSpectraRelErrorAll", 1300, 800);
  if (isGenOnTheFly) canvasPtSpectraRelErrorAll->Divide((nummoltMax+2)/2, 2);
  else if (isppHM || MultBinning==3)   canvasPtSpectraRelErrorAll->Divide(2,2);  
  else  canvasPtSpectraRelErrorAll->Divide(3,2);

  Float_t HighPtRatio[numtipo] = {1.15, 0, 0, 0, 0, 0, 0, 0, 1.2, 0};
  Float_t LowPtRatio[numtipo] = {0.85, 0, 0, 0, 0, 0, 0, 0, 0.8, 0};
  if (type==8 && TypeAnalysis==0) {
    LowPtRatio[type] = 0.6; 
    HighPtRatio[type] = 1.4;
  }

  cout << "\nprendo histo per confronto con dati pubblicati " << endl;                                        
  TString PathDatiPubblicati ="";                                                                               
  if (type==0) PathDatiPubblicati = "HEPData-ins1748157-v1-Table_1.root";                                       
  else if (type==8) PathDatiPubblicati = "HEPData-1583750454-v1-Table_3.root";                                  
  TFile *filedatipubblicati = new TFile(PathDatiPubblicati, "");                                                
  if (!filedatipubblicati) {cout << "file dati pubblicati not there " << endl; return;}                            TString STable = "Table 3";                                                                                      if (type==0) STable = "Table 1";                                                                              
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
  TH1F* hRatioFakeSBDef[nummolt+1];
  Float_t BarlowVarFakeSBDef[nummolt+1][numPtV0]={0};

  Float_t ErrDPhi[nummolt+1][numPtV0]={0};
  Int_t ErrDPhiCounter[nummolt+1]={0};

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
  Float_t LimInfPtvsMultTight=10e-5;
  Float_t LimSupPtvsMultTight=3;

  Int_t nrebin[nummolt+1] ={8, 8, 8, 8, 8, 8};
  if (type==0){
    if (TypeAnalysis==0) LimSupPtvsMult=4; //2.3
    if (TypeAnalysis==1) LimSupPtvsMult=2; //1.2
    if (TypeAnalysis==2) LimSupPtvsMult=2; //1.2
    if (TypeAnalysis==0) LimSupPtvsMultTight=2.7;
    if (TypeAnalysis==1) LimSupPtvsMultTight=1.5;
    if (TypeAnalysis==2) LimSupPtvsMultTight=1.5;
    for (Int_t m=0; m<=nummoltMax+1; m++){
      if (isppHM) nrebin[m] = 8;
    }
  }

  if (type==8){
    if (TypeAnalysis==0) LimSupPtvsMult=8; //5
    if (TypeAnalysis==1) LimSupPtvsMult=4; //2.4
    if (TypeAnalysis==2) LimSupPtvsMult=4;
    if (TypeAnalysis==0) LimSupPtvsMultTight=4;
    if (TypeAnalysis==1) LimSupPtvsMultTight=2.5;
    if (TypeAnalysis==2) LimSupPtvsMultTight=2.5;
    for (Int_t m=0; m<=nummoltMax+1; m++){
      if (isppHM) nrebin[m] = 8;
      else if (ispp5TeV) nrebin[m] = 8;
      else {
	if (TypeAnalysis==0)	nrebin[m] = 16;
	else if (m==3 || m==4) 	nrebin[m] = 16;
      }
    }
  }
  if (type==0){
    if (TypeAnalysis==0) LimInfPtvsMult=0.8;
    if (TypeAnalysis==1) LimInfPtvsMult=0.8;
    if (TypeAnalysis==2) LimInfPtvsMult=0.8;
    if (TypeAnalysis==0) LimInfPtvsMultTight=1.5;
    if (TypeAnalysis==1) LimInfPtvsMultTight=0.5;
    if (TypeAnalysis==2) LimInfPtvsMultTight=0.5;
  }
  if (type==8){
    if (TypeAnalysis==0) LimInfPtvsMult=1.2;
    if (TypeAnalysis==1) LimInfPtvsMult=1.2;
    if (TypeAnalysis==2) LimInfPtvsMult=1.2;
  }

  Float_t LimSupYieldErr=0.04;
  Float_t LimSupPtErr=0.04;
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
    if (ispp5TeV) LimSup = 0.12;
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
    LimSupError=0.7-10e-4;
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
  Float_t  UpBinDPhi[numsysDPhi]={0};


  //first part: I get default spectra with stat uncertainty
  cout << "\n**************************"<<endl;
  cout << "First part: I get default spectra with stat uncertainty";
  TString  PathInDef=PathIn0;
  PathInDef+=   Form("_SysPhi%i_PtMin%.1f_", 0, PtTrigMin);
  PathInDef+= RegionType[TypeAnalysis];
  //  PathInDef += "_OldDPhiRange";
  if (MultBinning!=0) PathInDef += Form("_MultBinning%i", MultBinning);
  if (isEffCorr) PathInDef += "_EffCorr";
  if (isWingsCorrectionApplied) {
    PathInDef += "_isWingsCorrectionApplied";
    if (isGenOnTheFly){
      if (type==8) PathInDef += "New";
    }
    else PathInDef += "New";
  }
  if (sys!=0) PathInDef += Form("_dEtaSys%i", sys);
  if (MaterialBudgetCorr==1) PathInDef += "_MatBudgetCorr";
  if (isTopoSel!=0){
    if (isTopoSel==1) PathInDef += "_SysV0Default";
    else if (isTopoSel==2) PathInDef += "_SysV0Tightest";
    else if (isTopoSel==3) PathInDef += "_SysV0Loosest";
  }
  PathInDef += ".root";
  cout << " from the file: " << PathInDef << endl;
  TFile *  fileinDef = new TFile(PathInDef, "");
  if (!fileinDef) {cout << PathInDef << " does not exist" << endl; return;}
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && m<2) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //    cout << Smolt[m] << endl;
    fHistSpectrumSistRelErrorDCAz[m]    =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorDCAz_"+Smolt[m]); 
    fHistSpectrumSistRelErrorPurity[m]    =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorPurity_"+Smolt[m]); 
    fHistSpectrumSistRelErrorSE[m]      =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorSE_"+Smolt[m]);
    fHistSpectrumSistRelErrorDeltaEta[m]=(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorDeltaEta_"+Smolt[m]);
    if (!fHistSpectrumSistRelErrorDCAz[m]) {cout << "Rel syst uncertainty DCAz " << endl; return;}
    if (!fHistSpectrumSistRelErrorPurity[m]) {cout << "Rel syst uncertainty purity " << endl; return;}
    if (!fHistSpectrumSistRelErrorSE[m]) {cout << "Rel syst uncertainty topological selections " << endl; return;}
    if (!fHistSpectrumSistRelErrorDeltaEta[m]){cout << "Rel syst uncertainty DeltaEta " << endl; return;}
    fHistSpectrumStat[m]=(TH1F*)fileinDef->Get("fHistSpectrum_"+Smolt[m]);
    fHistSpectrumStatUp[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumUp_"+Smolt[m]);
    fHistSpectrumStatDown[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumDown_"+Smolt[m]);
    fHistSpectrumStatHard[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumHard_"+Smolt[m]);
    fHistSpectrumStatSoft[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSoft_"+Smolt[m]);
    fHistSpectrumTemp[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumTemp_"+Smolt[m]);
    fHistSpectrumSist[m]=(TH1F*)fileinDef->Get("fHistSpectrumSist_"+Smolt[m]);
    hDeltaPhiLimit= (TH1F*) fileinDef-> Get("DeltaPhiLimit");
    if (!hDeltaPhiLimit) {cout << "DetaPhiLimit not there " << endl;return;}
    LowBinDPhi[0] =	hDeltaPhiLimit->GetBinContent(1);
    UpBinDPhi[0] =	hDeltaPhiLimit->GetBinContent(2);

    if (!fHistSpectrumStat[m]) {cout << "fHistSpectrumStat does not exist" << endl; return;}
    if (!fHistSpectrumSist[m]) {cout << "fHistSpectrumSist does not exist" << endl; return;}

    fHistSpectrumSistRelErrorAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorAll_"+Smolt[m]);
    fHistSpectrumSistRelErrorMultUnCorr[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorMultUnCorr_"+Smolt[m]);
    fHistSpectrumStatRelError[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumStatRelError_"+Smolt[m]);
    fHistSpectrumSistRelErrorDPhi[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorDPhi_"+Smolt[m]);
    fHistSpectrumSistRelErrorMB[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMB_"+Smolt[m]);
    fHistSpectrumSistRelErrorMBCorrection[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMBCorrection_"+Smolt[m]);
    fHistSpectrumSistRelErrorOOBPU[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorOOBPU_"+Smolt[m]);
    fHistSpectrumSistRelErrorIBPU[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorIBPU_"+Smolt[m]);
    fHistSpectrumSistRelErrorLeadTrack[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorLeadTrack_"+Smolt[m]);
    fHistSpectrumSistRelErrorMCChoice[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMCChoice_"+Smolt[m]);
    fHistSpectrumSistRelErrorFakeSB[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorFakeSB_"+Smolt[m]);
    fHistSpectrumSistRelErrorOOJSubDef[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorOOJSubDef_"+Smolt[m]);

    fHistSpectrum_max[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMax_"+Smolt[m]);
    fHistSpectrum_min[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMin_"+Smolt[m]);
    fHistSpectrumSistDPhi[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDPhi_"+Smolt[m]);
    fHistSpectrumSistOOJ[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistOOJ_"+Smolt[m]);
    fHistSpectrumSistAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistAll_"+Smolt[m]);
    fHistSpectrumSistMultUnCorr[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistMultUnCorr_"+Smolt[m]);
    fHistSpectrumSistLeadTrack[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistLeadTrack_"+Smolt[m]);
    fHistSpectrumSistMCChoice[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistMCChoice_"+Smolt[m]);
    fHistSpectrumSistFakeSB[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistFakeSB_"+Smolt[m]);
    fHistSpectrumSistOOJSubDef[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistOOJSubDef_"+Smolt[m]);

    for(Int_t v=PtV0Min; v <    fHistSpectrumSistAll[m]->GetNbinsX() ; v++){
      //      cout << "v " << v << endl;
      if (fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSist[m]->GetBinError(v+1)/ fHistSpectrumSist[m]->GetBinContent(v+1)); //this works for inclusive, for jet and OOJ I will change it
	fHistSpectrumStatRelError[m]->SetBinContent(v+1, fHistSpectrumStat[m]->GetBinError(v+1)/ fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      else       if (type==0 && NPtV0[v] == 0.1){ //?
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

  //Rel uncertainties taken from Fiorella's analysis
  Float_t  Sigma2OOBPU[nummolt+1] ={0};
  Float_t  Sigma2IBPU[nummolt+1] ={0};
  Float_t  Sigma2MB[nummolt+1]={0};   
  Float_t OOBPU = 0.012; //relative syst. uncertainty on pT spectrum associated to OOB PileUp taken from Fiorella's paper and "averaged" over pT
  Float_t IBPU=0.02; //same but for Inbunch pileup
  Float_t MB = 0.01; //same but for material budget
  Float_t MBCorrection = 0; //when I apply MaterialBudgetCorr == 2 I apply a FAST correction  =  the ratio between the efficiency w the correct Material Budget description and the old efficiency. 
  //Now this uncertainty is set to zero since all the uncertainty related to the material budget is inlcuded in MB

  if (type==8){
    OOBPU = IBPU = MB = 0.02;
  }
  if (isMC && !isEfficiency){
    OOBPU = IBPU = MB = 0;
  }

  //Pt-dependent material budget
  Float_t MBPtXi[numPtV0+1] = {0, 0.05, 0.03, 0.022, 0.018, 0.015, 0.01, 0.006 };
  Float_t MBPtK0s[numPtV0+1] = {0, 0.03, 0.017, 0.01, 0.007, 0.006, 0.005, 0.005, 0.005, 0.005};

  //second part: I evaluate syst uncertainty associated to choice of DeltaPhi (for jet and out of jet)
  Int_t NSign=3;

  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistDPhi[m]->SetBinError(v+1, 0);
    }
  }

  if (TypeAnalysis!=2){
    cout << "\n**************************"<<endl;
    cout << "Second part: I evaluate syst uncertainty associated to choice of DeltaPhi (for jet and out of jet)\n " << endl;
    cout << "2a. Evaluating Barlow significance of DeltaPhi change" << endl;
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      cout << " m " << m << endl;
      for (Int_t sysDPhi=1; sysDPhi<numsysDPhi; sysDPhi++){
	if (TypeAnalysis==0 && sysDPhi>2) continue;
	if (TypeAnalysis==0 && sysDPhi==1) continue; //|dphi| < 0.86 does not include all of the jet, I cannot take it as systematic uncertainty
	if (TypeAnalysis==0 && type==0 && MonashTune==3) continue;
	PathIn=PathIn0;
	PathIn+=   Form("_SysPhi%i_PtMin%.1f_", sysDPhi, PtTrigMin);
	PathIn+= RegionType[TypeAnalysis];
	if (MultBinning!=0) PathIn += Form("_MultBinning%i", MultBinning);
	if (isEffCorr) PathIn += "_EffCorr";
	if (isWingsCorrectionApplied) {
	  PathIn += "_isWingsCorrectionApplied";
	  if (isGenOnTheFly){
	    if (type==8) PathIn += "New";
	  }
	  else PathIn += "New";
	}
	if (sys!=0) PathIn += Form("_dEtaSys%i", sys);
	if (MaterialBudgetCorr==1) PathIn += "_MatBudgetCorr";
	/*
	if (isTopoSel!=0){
	  if (isTopoSel==1) PathIn += "_SysV0Default";
	  else if (isTopoSel==2) PathIn += "_SysV0Tightest";
	  else if (isTopoSel==3) PathIn += "_SysV0Loosest";
	}
	*/
	PathIn+= ".root";
	cout << "" << PathIn << endl;
	filein = new TFile(PathIn, "");
	if (!filein) {cout << "file is missing " << endl; return;}
	hDeltaPhiLimit= (TH1F*) filein-> Get("DeltaPhiLimit");
	if (!hDeltaPhiLimit) {cout << "DetaPhiLimit not there " << endl;return;}
	LowBinDPhi[sysDPhi] =	hDeltaPhiLimit->GetBinContent(1);
	UpBinDPhi[sysDPhi] =	hDeltaPhiLimit->GetBinContent(2);
	fHistSpectrumStatDPhi[m][sysDPhi]=(TH1F*)filein->Get("fHistSpectrum_"+Smolt[m]);
	fHistSpectrumStatDPhi[m][sysDPhi]->SetName(Form("fHistSpectrumStat_sysDPhi%i", sysDPhi)+Smolt[m]);

	StyleHisto(fHistSpectrumStatDPhi[m][sysDPhi], 0, LimSup, ColorsysDPhi[sysDPhi], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
	StyleHisto(fHistSpectrumStat[m], 0, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m],  0, 0, 0);

	if (isppHM)     canvasPtSpectraAll->cd(m+1-2);
	else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraAll->cd(3);
	else     canvasPtSpectraAll->cd(m+1);
	//      cout << "I'm drawing on the canvas " << endl;
	gPad->SetLeftMargin(0.15);
	if(m==nummoltMax && sysDPhi==1) legendPhi->AddEntry(fHistSpectrumStat[m],Form("%.2f-%.2f", LowBinDPhi[0], UpBinDPhi[0]) ,"l");
	if(m==nummoltMax) legendPhi->AddEntry(fHistSpectrumStatDPhi[m][sysDPhi],Form("%.2f-%.2f", LowBinDPhi[sysDPhi], UpBinDPhi[sysDPhi]) ,"l");
	//	if (sysDPhi==numsysDPhi-1) fHistSpectrumStat[m]->DrawClone("");
	fHistSpectrumStat[m]->DrawClone("same");
	fHistSpectrumStatDPhi[m][sysDPhi]->Draw("same");
	if (sysDPhi==numsysDPhi-1) legendPhi->Draw("same");
	if (TypeAnalysis==0 && sysDPhi==numsysDPhi-2) legendPhi->Draw("same");

	hBarlowVar[m][sysDPhi]=(TH1F*)    fHistSpectrumStatDPhi[m][sysDPhi]->Clone("fHistBarlowVar_"+Smolt[m]);
	for ( Int_t b=1; b<=hBarlowVar[m][sysDPhi]->GetNbinsX(); b++){
	  hBarlowVar[m][sysDPhi]->SetBinContent(b,0);
	  hBarlowVar[m][sysDPhi]->SetBinError(b,0);
	}
	for(Int_t v=PtBinMin[m]; v < numPtV0Max; v++){
	  //	  cout << " v " << v << " " << PtBinMin[m] << " " << NPtV0[v] << endl;
	  if (v==1 && TypeAnalysis==0 && type==8) continue;
	  if (NPtV0[v] == 4 && TypeAnalysis==0 && type==8 && m==4) continue;
	  if (NPtV0[v] == 4 && type==0 && m==4) continue;
	  //	cout << "\nv: " << v << endl;
	  //	cout << "...histo...dphi " << sysDPhi << " " << fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1)<< endl;
	  BarlowVar[m][v][sysDPhi] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatDPhi[m][sysDPhi]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) - pow(fHistSpectrumStatDPhi[m][sysDPhi]->GetBinError(v+1),2)));
	  if (fHistSpectrumStat[m]->GetBinContent(v+1) ==0)  BarlowVar[m][v][sysDPhi] =0;
	  if (TMath::Abs(BarlowVar[m][v][sysDPhi])>2) BarlowSign[m][sysDPhi]++;
	  hBarlowVar[m][sysDPhi] ->SetBinContent(v+1, BarlowVar[m][v][sysDPhi]) ;
	  hBarlowVar[m][sysDPhi] ->SetBinError(v+1, 0) ;
	  //cout << "barlow var " << NPtV0[v] << " "  << 	  hBarlowVar[m][sysDPhi] ->GetBinContent(v+1)<<" " <<	  hBarlowVar[m][sysDPhi] ->GetBinError(v+1) << endl;
	}//end loop v
	if (BarlowSign[m][sysDPhi]>= NSign) BarlowPassed[m][sysDPhi]=1;
	Int_t LimDPhi = 5;
	if (isppHM && type==0){
	  LimDPhi =6;
	}
	StyleHisto(hBarlowVar[m][sysDPhi], -LimDPhi, LimDPhi, ColorsysDPhi[sysDPhi], 33, titleX, "N_{#sigma}(Barlow) = (Def - Var)/Err",  title+SmoltLegend[m],  0, 0, 0);
	if ( BarlowPassed[m][sysDPhi])       StyleHisto(hBarlowVar[m][sysDPhi], -LimDPhi, LimDPhi, ColorsysDPhi[sysDPhi], 27,  titleX, "N_{#sigma}(Barlow) = (Def - Var)/Err",  title+SmoltLegend[m],  0, 0, 0);

	//	cout << "Barlow passed " << BarlowPassed[m][sysDPhi] << endl;
	if (isppHM)     canvasBarlow->cd(m+1-2);
	else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasBarlow->cd(3);
	else     canvasBarlow->cd(m+1);
	gPad->SetLeftMargin(0.15);
	hBarlowVar[m][sysDPhi] ->Draw("same p");
	lineat0->Draw("same");

	if (!BarlowPassed[m][sysDPhi]) continue;
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
	fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1,0);
	fHistSpectrumSistRelErrorDPhi[m]->SetBinError(v+1,0);
	if (v==1 && TypeAnalysis==0 && type==0) continue;
	if (NPtV0[v]==4 && type==0 && m==4) continue;
	ErrDPhi[m][v] = TMath::Abs(fHistSpectrum_min[m]->GetBinContent(v+1) - fHistSpectrum_max[m]->GetBinContent(v+1))/2;
	fHistSpectrumSistDPhi[m]->SetBinError(v+1, ErrDPhi[m][v]);

	if (fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  if (ErrDPhi[m][v]!=0)	  ErrDPhiCounter[m]++;
	  fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1, ErrDPhi[m][v]/fHistSpectrumStat[m]->GetBinContent(v+1));
	}
	fHistSpectrumSistRelErrorDPhi[m]->SetBinError(v+1,0);
      } //end new
    }//end loop on m

    for(Int_t m=0; m<=nummoltMax; m++){
      cout << "Multiplicity " << m << endl;
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	//	if (v==1 && TypeAnalysis==0 && type==0) continue; //?
	//	if (NPtV0[v]==4 && type==0 && m==4) continue; //?
	if (TypeAnalysis==0){
	  //	  cout << "ErrDPhi " << ErrDPhi[nummolt][v] << endl;
	  //	  cout << "stat " << fHistSpectrumStat[nummolt]->GetBinContent(v+1) << endl;
	  fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1, ErrDPhi[nummoltMax][v]/fHistSpectrumStat[nummoltMax]->GetBinContent(v+1));
	  fHistSpectrumSistDPhi[m]->SetBinError(v+1, ErrDPhi[nummoltMax][v]/fHistSpectrumStat[nummoltMax]->GetBinContent(v+1)* fHistSpectrumStat[m]->GetBinContent(v+1));
	  cout << fHistSpectrumSistRelErrorDPhi[m]->GetBinContent(v) << endl;
	}
      }

      //SMOOTH the relative uncertainty and compute the absolute uncertainty
      fHistSpectrumSistRelErrorDPhi[m]->GetXaxis()->SetRangeUser(LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
      fHistSpectrumSistRelErrorDPhi[m]->Smooth(1, "R");
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (fHistSpectrumStat[m]->GetBinContent(v+1)==0){
	  fHistSpectrumSistRelErrorDPhi[m]->SetBinContent(v+1,0);
	  fHistSpectrumSistDPhi[m]->SetBinError(v+1,0);
	}
	else fHistSpectrumSistDPhi[m]->SetBinError(v+1, fHistSpectrumSistRelErrorDPhi[m]->GetBinContent(v+1) * fHistSpectrumStat[m]->GetBinContent(v+1));
      }
    }
  } //end of type analysis

  //diagnostic cout:
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (m!=0) continue; //choose the multiplicity for diagnosis
    cout << "\nSyst. uncert. associated to choice of DeltaPhi,  (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " +- " << fHistSpectrumSistDPhi[m]->GetBinError(b) << " (" << fHistSpectrumSistDPhi[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ") " << endl;
    }
  }

  //here syst related to pt < pt,trigg 
  //cout << "\n*******************************************"<<endl;
  //cout<< "hXi: systematic effect associated to selection pt,assoc < pt,Trigg  " << endl;
  Float_t   YieldSpectrumErrLeadTrack[nummolt+1]={0};
  TString    PathInLeadTrackDef;
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistLeadTrack[m]->SetBinError(v+1,0);
    }
  }
  if (type==8 && kFALSE){
    cout << "\n*******************************************"<<endl;
    cout<< "hXi: systematic effect associated to selection pt,assoc < pt,Trigg  " << endl;
    PathInLeadTrackDef = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    if (type==8 && TypeAnalysis==0){
      PathInLeadTrackDef += "_OOJNoTriggerSmoothed";
    }
    if(type>=0){
      PathInLeadTrackDef +="_"+tipo[type];
      PathInLeadTrackDef +=Srap[israp];
      PathInLeadTrackDef +="_AllAssoc";
    }
    //    if (IsHighPtExtr) PathInLeadTrackDef+="_HighPtExtr";
    PathInLeadTrackDef+= hhCorr[0]+"_"+RegionTypeOld[TypeAnalysis] +DataOrMC[isMC] + Form("_PtMin%.1f", PtTrigMin);
    PathInLeadTrackDef+="_IsEtaEff";
    if (MultBinning!=0) PathInLeadTrackDef += Form("_MultBinning%i", MultBinning);
    PathInLeadTrackDef += ".root";
    cout << "\n\n" << PathInLeadTrackDef << endl;
    TFile *  fileinLeadTrackDef = new TFile(PathInLeadTrackDef, "");
    if (!fileinLeadTrackDef) {cout << PathInLeadTrackDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
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
      StyleHisto(hBarlowVarLeadTrackDef[m], -10, 10, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);
      if (BarlowPassedLeadTrackDef[m]) StyleHisto(hBarlowVarLeadTrackDef[m], -10, 10, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);

      //      cout << "Barlow passed " << BarlowPassedLeadTrackDef[m] << endl;
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
	//	if (NPtV0[v] <3 || TypeAnalysis==2){
	if (NPtV0[v] <3){
	  fHistSpectrumSistLeadTrack[m]->SetBinError(v+1,0);
	  continue;
	}
	if (v==numPtV0Max-1 && m==4 && TypeAnalysis==0 && type==8) continue;
	if (fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1) ==0) continue;
	//	if(TMath::Abs(BarlowVarLeadTrackDef[m][v])>2){
	if (BarlowPassedLeadTrackDef[m]){
	  fHistSpectrumSistLeadTrack[m]->SetBinError(v+1, TMath::Abs(fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1))/2);
	  //YieldSpectrumErrLeadTrack[m]+= TMath::Abs(   fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatLeadTrackDef[m]->GetBinContent(v+1));
	}
      }
    }//end loop m
  }//end loop syst associated to choice pt <pt,trigg
  //diagnostic cout:
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (m!=0) continue; //choose the multiplicity for diagnosis
    cout << "\nSyst. uncert. associated to choice pt < pt,trigg (not inlcuded now, should be zero) ,  (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " +- " << fHistSpectrumSistLeadTrack[m]->GetBinError(b) << " (" << fHistSpectrumSistLeadTrack[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ") " << endl;
    }
  }


  cout << "\n*******************************************"<<endl;
  cout<< "Systematic effect associated to MC used to calculate efficiency  " << endl;
  Float_t   YieldSpectrumErrMCChoice[nummolt+1]={0};
  TString    PathInMCChoiceDef;
  TLegend * legendMCChoice= new TLegend(0.6, 0.6, 0.9, 0.9);
  Int_t IsSignMCChoice=0;
  Int_t Varm = 0;
  TString PathMCChoiceSyst="";

  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistMCChoice[m]->SetBinError(v+1,0);
    }
  }

  if (type==0 && kFALSE){
    // cout << "\n*******************************************"<<endl;
    // cout<< "Systematic effect associated to choice of MC used to calculat K0s/Xi efficiency  " << endl;
    PathInMCChoiceDef = Dir+"/DATA"+year0+"/SystematicAnalysis" + year;
    if (PtBinning>0) PathInMCChoiceDef +=Form("_PtBinning%i",PtBinning);
    if (type==8 && TypeAnalysis==0) PathInMCChoiceDef += "_OOJSmoothedBis";
    if(type>=0){
      PathInMCChoiceDef +="_"+tipo[type];
      PathInMCChoiceDef +=Srap[israp];
      PathInMCChoiceDef +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInMCChoiceDef+="_HighPtExtr";
    PathInMCChoiceDef+= hhCorr[0]+"_"+RegionTypeOld[TypeAnalysis] +DataOrMC[isMC] + Form("_PtMin%.1f_EPOS_IsEtaEff.root", PtTrigMin);
    cout << "\n\n" << PathInMCChoiceDef << endl;
    TFile *  fileinMCChoiceDef = new TFile(PathInMCChoiceDef, "");
    if (!fileinMCChoiceDef) {cout << PathInMCChoiceDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
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
      if (m==nummoltMax)      legendMCChoice->AddEntry(fHistSpectrumStat[m], "PYTHIA8", "ple");
      fHistSpectrumStat[m]->DrawClone("ep");
      StyleHisto(fHistSpectrumStatMCChoiceDef[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==nummoltMax)      legendMCChoice->AddEntry(fHistSpectrumStatMCChoiceDef[m], "EPOS", "ple");
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
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (IsSignMCChoice>0){
	  fHistSpectrumSistMCChoice[m]->SetBinError(v+1,	  (fHistSpectrumStat[Varm]->GetBinContent(v+1) -    fHistSpectrumStatMCChoiceDef[Varm]->GetBinContent(v+1))/2/fHistSpectrumSistMCChoice[Varm]->GetBinContent(v+1)*fHistSpectrumSistMCChoice[m]->GetBinContent(v+1));
	}
      }
    }
  }//end loop syst associated to MC used to calculate efficiency
  else if (type==0 && !isGenOnTheFly){
    //set pt-independent relative uncertainty of 3%
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistMCChoice[m]->SetBinError(v+1, 0.01*fHistSpectrumSistMCChoice[m]->GetBinContent(v+1));
      }
    }
  }
  else if (type==8){
    //set pt-independent relative uncertainty of 3%
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	//	fHistSpectrumSistMCChoice[m]->SetBinError(v+1, 0.03*fHistSpectrumSistMCChoice[m]->GetBinContent(v+1));
	fHistSpectrumSistMCChoice[m]->SetBinError(v+1, 0);
      }
    }
  }

  //diagnostic cout:
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (m!=nummoltMax) continue; //choose the multiplicity for diagnosis
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
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistFakeSB[m]->SetBinError(v+1,0);
    }
  }
  cout << "\n*******************************************"<<endl;
  cout<< "Systematic effect associated to how fake K0s/Xi are sub  " << endl;
  if (type==8 && !isppHM && !ispp5TeV && !isGenOnTheFly){
    PathInFakeSBDef = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    if (PtBinning>0) PathInFakeSBDef +=Form("_PtBinning%i",PtBinning);
    if (type==8 && isppHM && TypeAnalysis==0) PathInFakeSBDef +="_OOJAllMult";
    if(type>=0){
      PathInFakeSBDef +="_"+tipo[type];
      PathInFakeSBDef +=Srap[israp];
      PathInFakeSBDef +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInFakeSBDef+="_HighPtExtr";
    PathInFakeSBDef+= hhCorr[0]+"_"+RegionTypeOld[TypeAnalysis] +DataOrMC[isMC] + Form("_PtMin%.1f", PtTrigMin);
    PathInFakeSBDef += "_Sidebands";
    //    if (type==8)  PathInFakeSBDef += "_SBMEFromPeak";
    PathInFakeSBDef+="_IsEtaEff";
    if (MultBinning!=0) PathInFakeSBDef += Form("_MultBinning%i", MultBinning);
    if (type==0 && TypeAnalysis==0 && isppHM) PathInFakeSBDef += "_NewdEtaChoice";
    PathInFakeSBDef += ".root";

    //    cout << "\nFrom file: " << PathInFakeSBDef << endl;
    TFile *  fileinFakeSBDef = new TFile(PathInFakeSBDef, "");
    if (!fileinFakeSBDef) {cout << PathInFakeSBDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      //    cout << " m " << m << endl;
      fHistSpectrumStatFakeSBDef[m]    =(TH1F*)fileinFakeSBDef->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      if (!fHistSpectrumStatFakeSBDef[m]) {cout << " I was looking for histo in " << PathInFakeSBDef  << endl; return;}

      hBarlowVarFakeSBDef[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistBarlowVarFakeSBDef_"+Smolt[m]);
      hRatioFakeSBDef[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistRatioFakeSBDef_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	//	if (NPtV0[v] <3) continue;
	if (fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1) ==0) continue;
	BarlowVarFakeSBDef[m][v] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) -pow(fHistSpectrumStatFakeSBDef[m]->GetBinError(v+1),2)));
	if (TMath::Abs(BarlowVarFakeSBDef[m][v])>2) BarlowSignFakeSBDef[m]++;
	hBarlowVarFakeSBDef[m] ->SetBinContent(v+1, BarlowVarFakeSBDef[m][v]) ;
	hBarlowVarFakeSBDef[m] ->SetBinError(v+1, 0) ;
	//	cout << "barlow var v" << v << " "  << 	  hBarlowVarFakeSBDef[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarFakeSBDef[m] ->GetBinError(v+1) << endl;
      }//end loop v

      hRatioFakeSBDef[m]->Divide(fHistSpectrumStatFakeSBDef[m]);

      Int_t      LimFakeSB = 5;
      if (isppHM) {
	LimFakeSB = 10;
	if (TypeAnalysis==0) LimFakeSB = 25;
      }
      
      if (BarlowSignFakeSBDef[m]>= 3) BarlowPassedFakeSBDef[m]=1;
      StyleHisto(hBarlowVarFakeSBDef[m], -LimFakeSB, LimFakeSB, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow) (Def - Var)",  title+SmoltLegend[m], 0, 0, 0);
      if (BarlowPassedFakeSBDef[m]) StyleHisto(hBarlowVarFakeSBDef[m], -LimFakeSB, LimFakeSB, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow) (Def - Var)",  title+SmoltLegend[m], 0, 0, 0);

      //      cout << "Barlow passed " << BarlowPassedFakeSBDef[m] << endl;
      canvasBarlowFakeSBDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      hBarlowVarFakeSBDef[m] ->Draw("same p");
      lineat2->Draw("same");
      lineatm2->Draw("same");

      canvasRatioFakeSBDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(hRatioFakeSBDef[m], 0.8, 1.2, Color[TypeAnalysis], 27,  titleX, "Ratio Def/Var",  title+SmoltLegend[m], 0, 0, 0);
      hRatioFakeSBDef[m]->GetYaxis()->SetRangeUser(0.8, 1.2);
      hRatioFakeSBDef[m] ->Draw("same p");
      lineat1->Draw("same");

      canvasFakeSBDef->cd(m+1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==nummoltMax)      legendFakeSB->AddEntry(fHistSpectrumStat[m], "Default", "ple");
      fHistSpectrumStat[m]->DrawClone("ep");
      StyleHisto(fHistSpectrumStatFakeSBDef[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      if (m==nummoltMax)      legendFakeSB->AddEntry(fHistSpectrumStatFakeSBDef[m], "Sidebands", "ple");
      fHistSpectrumStatFakeSBDef[m]->DrawClone("same ep");
      legendFakeSB->Draw("");

      //      cout << " going to set errors " << endl;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistFakeSB[m]->SetBinError(v+1,0);
	if (fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1) ==0) continue;
	//	if(TMath::Abs(BarlowVarFakeSBDef[m][v])>2){
	//	if (BarlowPassedFakeSBDef[m]){
	if (m==nummoltMax)	{
	  fHistSpectrumSistFakeSB[m]->SetBinError(v+1,	  (fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatFakeSBDef[m]->GetBinContent(v+1))/2);
	}
	//	}
      }
    }//end loop m
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	if (m!=nummoltMax)	fHistSpectrumSistFakeSB[m]->SetBinError(v+1,	  (fHistSpectrumStat[5]->GetBinContent(v+1) -    fHistSpectrumStatFakeSBDef[5]->GetBinContent(v+1))/2/ fHistSpectrumSistFakeSB[5]->GetBinContent(v+1)* fHistSpectrumSistFakeSB[m]->GetBinContent(v+1));
	//TO SMOOTH
	fHistSpectrumSistRelErrorFakeSB[m]->SetBinContent(v+1,	 fHistSpectrumSistFakeSB[m]->GetBinError(v+1)/ fHistSpectrumSistFakeSB[m]->GetBinContent(v+1));
      }
    }
    for(Int_t m=0; m<nummoltMax +1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      //      cout << "Multiplicity " << m << endl;
      fHistSpectrumSistRelErrorFakeSB[m]->GetXaxis()->SetRangeUser(LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
      fHistSpectrumSistRelErrorFakeSB[m]->Smooth(1, "R"); //SMOOTH
      if (ispp5TeV && type==0 && TypeAnalysis==0)       fHistSpectrumSistRelErrorFakeSB[m]->Smooth(1, "R"); //SMOOTH
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistFakeSB[m]->SetBinError(v+1, fHistSpectrumSistRelErrorFakeSB[m]->GetBinContent(v+1)*fHistSpectrumSistFakeSB[m]->GetBinContent(v+1));
	/*
	cout << " Rel Error v " << v << " " << fHistSpectrumSistRelErrorFakeSB[m]->GetBinContent(v+1) << endl;
	cout << "Content v " << v << " " << fHistSpectrumSistFakeSB[m]->GetBinContent(v+1) << endl;
	cout << "Error v " << v << " " << fHistSpectrumSistFakeSB[m]->GetBinError(v+1) << endl;
	*/
      }
    }
  }//end loop syst associated to how fake K0s/Xi are subtracted

  TString PathInSBRelErr = "FinalOutput/DATA2016/PtSpectraBisNew_161718Full_AOD234_hXi_Xi_Eta0.8_AllAssoc_PtMin3.0_";
  PathInSBRelErr += RegionType[TypeAnalysis];
  PathInSBRelErr += "_isNormCorrFullyComputed_isErrorAssumedPtCorr_ChangesIncluded.root";
  if (type==8 && ispp5TeV){
    TFile * fileinFakeSBRelErr = new TFile(PathInSBRelErr);
    for(Int_t m=0; m<nummoltMax +1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      fHistSpectrumSistRelErrorFakeSB[m]    =(TH1F*)fileinFakeSBRelErr->Get("fHistSpectrumSistRelErrorFakeSB__all"); 
      if (!fHistSpectrumSistRelErrorFakeSB[m]) {cout << " I was looking for histo in " << PathInSBRelErr  << endl; return;}
      for (Int_t v=0; v< fHistSpectrumSistRelErrorFakeSB[m]->GetNbinsX(); v++){
	if (fHistSpectrumStat[m]->GetBinContent(v+1)==0)       fHistSpectrumSistRelErrorFakeSB[m]->SetBinContent(v+1, 0);
      }
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumSistFakeSB[m]->SetBinError(v+1, fHistSpectrumSistRelErrorFakeSB[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinContent(v+1));
      }
    }
  }
  
  //diagnostic cout:
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (m!=nummoltMax) continue; //choose the multiplicity for diagnosis
    cout << "\nSyst. fake K0s/Xi uncert. associated to spectrum (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " +- " << fHistSpectrumSistFakeSB[m]->GetBinError(b) << " (" << fHistSpectrumSistFakeSB[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ") " << endl;
    }
  }

  //*************************************************************************
  //for Xi in-jet production: systematic effect associated to ooj subtraction
  //*************************************************************************

  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fHistSpectrumSistOOJSubDef[m]->SetBinError(v+1,0);
    }
  }
  Float_t   YieldSpectrumErrOOJSub[nummolt+1]={0};
  TString    PathInOOJSubDef;
  if (TypeAnalysis==0 && type==8 && !isGenOnTheFly){
    cout << "\n*******************************************"<<endl;
    cout<< "hXi: systematic effect associated to OOJ subtraction " << endl;
    PathInOOJSubDef = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    if (isppHM) PathInOOJSubDef += "_OOJAllMult";/*"_OOJAllMultAllPt"; /*"_OOJAllMultAtHighPt";*/ //for Xi HM I compare the spectra for pt > 2.5 GeV/c, where both methods can be used
    if(type>=0){
      PathInOOJSubDef +="_"+tipo[type];
      PathInOOJSubDef +=Srap[israp];
      PathInOOJSubDef +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInOOJSubDef+="_HighPtExtr";
    PathInOOJSubDef+= hhCorr[0]+"_Jet" +DataOrMC[isMC] + Form("_PtMin%.1f", PtTrigMin);
    if (ispp5TeV && type==8) PathInOOJSubDef+="_IsMEFrom13TeV";
    PathInOOJSubDef+= "_IsEtaEff";
    if (MultBinning!=0) PathInOOJSubDef += Form("_MultBinning%i", MultBinning);
    PathInOOJSubDef+=".root";
    cout << "\n\n" << PathInOOJSubDef << endl;
    TFile *  fileinOOJSubDef = new TFile(PathInOOJSubDef, "");
    if (!fileinOOJSubDef) {cout << PathInOOJSubDef << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      //    cout << " m " << m << endl;
      fHistSpectrumStatOOJSubDef[m]    =(TH1F*)fileinOOJSubDef->Get(Form("fHistSpectrumPart_m%i_syst0", m)); 
      if (!fHistSpectrumStatOOJSubDef[m]) {cout << " I was looking for histo in " << PathInOOJSubDef  << endl; return;}

      hBarlowVarOOJSubDef[m]=(TH1F*)    fHistSpectrumStat[m]->Clone("fHistBarlowVarOOJSubDef_"+Smolt[m]);
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	cout << "spectrum " <<  fHistSpectrumStatOOJSubDef[m]->GetBinContent(v) << endl;
	if ((NPtV0[v] >2.5 || NPtV0[v]<2) && !isppHM)continue;
	//	if (isppHM && NPtV0[v] <2.5) continue;
	if (fHistSpectrumStatOOJSubDef[m]->GetBinContent(v+1) ==0) continue;
	BarlowVarOOJSubDef[m][v] = (    fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatOOJSubDef[m]->GetBinContent(v+1))/sqrt(TMath::Abs(pow( fHistSpectrumStat[m]->GetBinError(v+1),2) -pow(fHistSpectrumStatOOJSubDef[m]->GetBinError(v+1),2)));
	if (TMath::Abs(BarlowVarOOJSubDef[m][v])>2) BarlowSignOOJSubDef[m]++;
	hBarlowVarOOJSubDef[m] ->SetBinContent(v+1, BarlowVarOOJSubDef[m][v]) ;
	hBarlowVarOOJSubDef[m] ->SetBinError(v+1, 0) ;
	//cout << "barlow var v" << v << " "  << 	  hBarlowVarOOJSubDef[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarOOJSubDef[m] ->GetBinError(v+1) << endl;
      }//end loop v

      if (BarlowSignOOJSubDef[m]>= 1) BarlowPassedOOJSubDef[m]=1;
      if (isppHM) BarlowPassedOOJSubDef[m] = 1; //the only non significant multiplciity intevral is the 0.01-0.05% one, but it makes no sense to skip just one mult interval
      StyleHisto(hBarlowVarOOJSubDef[m], -5, 5, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);
      if (BarlowPassedOOJSubDef[m]) StyleHisto(hBarlowVarOOJSubDef[m], -5, 5, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0, 0, 0);

      //      cout << "Barlow passed " << BarlowPassedOOJSubDef[m] << endl;
      canvasBarlowOOJSubDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      hBarlowVarOOJSubDef[m] ->DrawClone("same p");

      canvasOOJSubDef->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStat[m]->DrawClone("ep");
      StyleHisto(fHistSpectrumStatOOJSubDef[m], LimInf, LimSup, 1, 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
      fHistSpectrumStatOOJSubDef[m]->DrawClone("same ep");

      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	cout << "spectrum " <<  fHistSpectrumStatOOJSubDef[m]->GetBinContent(v) << endl;
	fHistSpectrumSistOOJSubDef[m]->SetBinError(v+1,0);
	if (NPtV0[v] >=2.5 && !isppHM) continue;
	if (BarlowPassedOOJSubDef[m]){
	  if (isppHM ) {
	  }
	  else 	  fHistSpectrumSistOOJSubDef[m]->SetBinError(v+1,	  TMath::Abs(fHistSpectrumStat[m]->GetBinContent(fHistSpectrumStat[m]->FindBin(2.01)) -    fHistSpectrumStatOOJSubDef[m]->GetBinContent(fHistSpectrumStatOOJSubDef[m]->FindBin(2.01)))/2/ fHistSpectrumStat[m]->GetBinContent(fHistSpectrumStat[m]->FindBin(2.01))* fHistSpectrumStat[m]->GetBinContent(v+1));
	  //	  YieldSpectrumErrOOJSub[m]+= TMath::Abs(   fHistSpectrumStat[m]->GetBinContent(v+1) -    fHistSpectrumStatOOJSubDef[m]->GetBinContent(v+1));
	}
      }
    }//end loop m

    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue; 
      if (isppHM) {
	for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	  if (NPtV0[v]>=1.5 && NPtV0[v]< 2.5) {
	    fHistSpectrumSistOOJSubDef[m]->SetBinError(v+1,	  TMath::Abs(fHistSpectrumStat[3]->GetBinContent(v+1) -    fHistSpectrumStatOOJSubDef[3]->GetBinContent(v+1))/2/ fHistSpectrumStat[3]->GetBinContent(v+1)* fHistSpectrumStat[m]->GetBinContent(v+1)); //take uncertainty from 0.01-0.05% mult class
	    fHistSpectrumSistRelErrorOOJSubDef[m]->SetBinContent(v+1, fHistSpectrumSistOOJSubDef[m]->GetBinError(v+1)/fHistSpectrumSistOOJSubDef[m]->GetBinContent(v+1));
	  }
	  else {
	    fHistSpectrumSistOOJSubDef[m]->SetBinError(v+1,0);
	    fHistSpectrumSistRelErrorOOJSubDef[m]->SetBinContent(v+1,0);
	  }
	}
      }
    }
  } //end loop ooj sub for hXi

  //diagnostic cout:
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (type==0) continue;
    //    if (m!=nummoltMax && isppHM) continue; //choose the multiplicity for diagnosis
    if (m!=0 && !isppHM) continue; //choose the multiplicity for diagnosis
    cout << "\n*******************************************"<<endl;
    cout << "\nSyst. uncert. associated to how OOJ is subtracted (only for Jet) ,  (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " +- " << fHistSpectrumSistOOJSubDef[m]->GetBinError(b) << " (" << fHistSpectrumSistOOJSubDef[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ") " << endl;
    }
  }

  //third part: for the jet, I evaluate sys associated to type of OOJ subtraction (for the jet)
  TString    PathInpol0;
  if (kFALSE){
  //  if (TypeAnalysis==0 && type==0 && !isppHM && !ispp5TeV){
    cout << "\n*******************************************"<<endl;
    cout << "3. Syst. associated to pol0 OOJ subtraction for hK0s"<< endl;

    PathInpol0 = Dir+"/DATA"+year0+"/SystematicAnalysis" +year;
    if (PtBinning>0) PathInpol0 += "_PtBinning1";
    if(type>=0){
      PathInpol0 +="_"+tipo[type];
      PathInpol0 +=Srap[israp];
      PathInpol0 +=SSkipAssoc[SkipAssoc];
    }
    //    if (IsHighPtExtr) PathInpol0+="_HighPtExtr";
    PathInpol0+= hhCorr[0]+"_JetBFit" +DataOrMC[isMC] + Form("_PtMin%.1f", PtTrigMin);
    PathInpol0 += "_IsEtaEff";
    if (MultBinning!=0) PathInpol0 += Form("_MultBinning%i", MultBinning);
    PathInpol0 +=".root";
    cout << "\nFrom file: " << PathInpol0 << endl;
    TFile *  fileinpol0 = new TFile(PathInpol0, "");
    if (!fileinpol0) {cout << PathInpol0 << " does not exist" << endl; return;}

    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
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
	//	cout << "barlow var v" << v << " "  << 	  hBarlowVarpol0[m] ->GetBinContent(v+1)<<" " <<	  hBarlowVarpol0[m] ->GetBinError(v+1) << endl;

      }//end loop v

      if (BarlowSignpol0[m]>= NSign) BarlowPassedpol0[m]=1;
      StyleHisto(hBarlowVarpol0[m], -5, 5, Color[TypeAnalysis], 33, titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0,0,0);
      if (BarlowPassedpol0[m]) StyleHisto(hBarlowVarpol0[m], -5, 5, Color[TypeAnalysis], 27,  titleX, "N_{#sigma}(Barlow)",  title+SmoltLegend[m], 0,0,0);

      //      cout << "Barlow passed " << BarlowPassedpol0[m] << endl;
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

  cout << "************UNCERTAINTIES UNCORRELATED WITH MULTIPLICITY**************" << endl;

  //Topological selections
  TString sfileTopoSelUncorr = "UncertaintiesUncorrMult_" + tipo[type]+"_"+ RegionType[2] + "_TopoSel_Fixed.root";
  TH1F * hFractionUnCorrTopoSel[nummoltMax+1];
  TFile * fileTopoSelUncorr = new TFile(sfileTopoSelUncorr, "");
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;  
    if (m!=nummoltMax) {
      if (isppHM) hFractionUnCorrTopoSel[m] = (TH1F*)fileTopoSelUncorr->Get("hSpectrumFracUncorr_0-5");
      else if (ispp5TeV) hFractionUnCorrTopoSel[m] = (TH1F*)fileTopoSelUncorr->Get("hSpectrumFracUncorr_0-5");
      else {
	hFractionUnCorrTopoSel[m] = (TH1F*)fileTopoSelUncorr->Get("hSpectrumFracUncorr_"+Smolt[m]);
	if (type==8 && m==1) 	hFractionUnCorrTopoSel[m] = (TH1F*)fileTopoSelUncorr->Get("hSpectrumFracUncorr_0-5");
      }
      if (!hFractionUnCorrTopoSel[m]) {cout << "Missing histo: hFractionUnCorrTopoSel for mult: " << m << endl; return;} 
      hFractionUnCorrTopoSel[m]->SetName(Form("hSpectrumFracUncorr_TopoSel_m%i", m));
    }
    else { //define histo for 0-100% (and 0-0.1%) but fill it with zero
      hFractionUnCorrTopoSel[m] = (TH1F*)fileTopoSelUncorr->Get("hSpectrumFracUncorr_0-5"); 
      if (!hFractionUnCorrTopoSel[m]) {cout << "Missing histo: hFractionUnCorrTopoSel for mult: " << m << endl; return;} 
      hFractionUnCorrTopoSel[m]->SetName(Form("hSpectrumFracUncorr_TopoSel_m%i", m));
      for (Int_t b=1; b<=  hFractionUnCorrTopoSel[m]->GetNbinsX(); b++){
	hFractionUnCorrTopoSel[m]->SetBinContent(b, 0);
      }
    }
  }

  //DPhi selections
  TString sfileDPhiUncorr = "UncertaintiesUncorrMult_" + tipo[type]+"_"+ RegionType[TypeAnalysis] + "_dPhiSel_Fixed.root";
  TString sfileDPhiUncorrK0s = "";
  if (type==8 && TypeAnalysis==1) {
    sfileDPhiUncorrK0s = "UncertaintiesUncorrMult_" + tipo[0]+"_"+ RegionType[TypeAnalysis] + "_dPhiSel_Fixed.root";
  }
  TH1F * hFractionUnCorrDPhi[nummoltMax+1];
  TH1F * hFractionUnCorrDPhiK0s[nummoltMax+1];
  if (TypeAnalysis!=2){
    TFile * fileDPhiUncorr = new TFile(sfileDPhiUncorr, "");
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;  
      cout << "\nMult: " << Smolt[m] << endl;
      if (m!=nummoltMax) {
	if (ispp5TeV){
	  if (m==0)	  hFractionUnCorrDPhi[m] = (TH1F*)fileDPhiUncorr->Get("hSpectrumFracUncorr_10-30");
	  else if (m==1)  hFractionUnCorrDPhi[m] = (TH1F*)fileDPhiUncorr->Get("hSpectrumFracUncorr_30-50");
	}
	else {
	  hFractionUnCorrDPhi[m] = (TH1F*)fileDPhiUncorr->Get("hSpectrumFracUncorr_"+Smolt[m]);
	  if (type==8 && m==1 && !isppHM) hFractionUnCorrTopoSel[m] = (TH1F*)fileTopoSelUncorr->Get("hSpectrumFracUncorr_10-30");
	}
	if (!hFractionUnCorrDPhi[m]) {cout << "Missing histo: hFractionUnCorrDPhi for mult: " << m << endl; return;} 
	hFractionUnCorrDPhi[m]->SetName(Form("hSpectrumFracUncorr_DPhi_m%i", m));
	if (type==8 && TypeAnalysis==1) {
	  TFile * fileDPhiUncorrK0s = new TFile(sfileDPhiUncorrK0s, "");
	  if (ispp5TeV){
	    if (m==0)	  hFractionUnCorrDPhiK0s[m] = (TH1F*)fileDPhiUncorrK0s->Get("hSpectrumFracUncorr_10-30");
	    else if (m==1)  hFractionUnCorrDPhiK0s[m] = (TH1F*)fileDPhiUncorrK0s->Get("hSpectrumFracUncorr_30-50");
	  }
	  else hFractionUnCorrDPhiK0s[m] = (TH1F*)fileDPhiUncorrK0s->Get("hSpectrumFracUncorr_"+Smolt[m]);
	  if (!hFractionUnCorrDPhiK0s[m]) {cout << "Missing histo: hFractionUnCorrDPhiK0s for mult: " << m << endl; return;} 
	  hFractionUnCorrDPhiK0s[m]->SetName(Form("hSpectrumFracUncorr_DPhi_K0s_m%i", m));
	  for (Int_t b=1; b<=  hFractionUnCorrDPhiK0s[m]->GetNbinsX(); b++){
	    Float_t pt = hFractionUnCorrDPhiK0s[m]->GetBinCenter(b);
	    cout << "pt K0s: " << pt << " " <<  hFractionUnCorrDPhiK0s[m]->GetBinContent(b) << endl;
	  }
	}
      }
      else {//define histo for 0-100% (and 0-0.1%) but fill it with zero
	hFractionUnCorrDPhi[m] = (TH1F*)fileDPhiUncorr->Get("hSpectrumFracUncorr_0-5");
	if (!hFractionUnCorrDPhi[m]) {cout << "Missing histo: hFractionUnCorrDPhi for mult: " << m << endl; return;} 
	hFractionUnCorrDPhi[m]->SetName(Form("hSpectrumFracUncorr_DPhi_m%i", m));
	for (Int_t b=1; b<=  hFractionUnCorrDPhi[m]->GetNbinsX(); b++){
	  hFractionUnCorrDPhi[m]->SetBinContent(b, 0);
	}
      }
    }
  }

  if (type==8 && TypeAnalysis==1) {
    for(Int_t m=0; m<nummoltMax; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;  
      cout << "\nMult: " << Smolt[m] << endl;
      for (Int_t b=1; b<=  hFractionUnCorrDPhi[m]->GetNbinsX(); b++){
	Float_t pt = hFractionUnCorrDPhi[m]->GetBinCenter(b);
	hFractionUnCorrDPhi[m]->SetBinContent(b, hFractionUnCorrDPhiK0s[m]->GetBinContent(hFractionUnCorrDPhiK0s[m]->FindBin(pt)));
	cout << "pt: " << pt << " " <<  hFractionUnCorrDPhi[m]->GetBinContent(b) << endl;
      }
    }
  }

  //    cout << "Correction by efficiency ratio between 22e1 and old efficiency; The MC 22e1 is characterised by a R-dependent material budget" << endl;
  TString PathInMBCorr = "";
  PathInMBCorr = "FinalOutput/DATA2016/Efficiency/RatioBetweenEficiency22e1ToEfficiencyOld";
  if (type==0) PathInMBCorr += "_K0s.root";
  else  PathInMBCorr += "_Xi.root";
  TFile * fileMBCorr;
  TH1F * histoMBCorr[nummoltMax+1];
  TH1F * histoMBCorrForError[nummoltMax+1];

  if (MaterialBudgetCorr!=0){
    cout << "Correction by efficiency ratio between 22e1 and old efficiency; The MC 22e1 is characterised by a R-dependent material budget" << endl;
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      if (m!=0 && !isppHM) continue; //choose the multiplicity for diagnosis
      if (isppHM && m!=5) continue;
      cout << "\n*******************************************"<<endl;
      cout << "\nStat and syst uncertainty BEFORE correction by efficiency scaling,  (mult "<< m << ")" << endl;
      for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
	cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " (stat. rel. error: " << fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ", syst. rel. error:  "<< fHistSpectrumSistAll[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << " )" << endl;
	//      cout << fHistSpectrumSistMultUnCorr[m]->GetBinError(b+1) << " " << fHistSpectrumSistAll[m]->GetBinError(b+1) << endl;
      }
    }

    fileMBCorr = new TFile(PathInMBCorr);

    for (Int_t m=0; m<=nummoltMax; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      if (ispp5TeV) {
	histoMBCorr[m] = (TH1F*)fileMBCorr->Get(Form("fHistCorrFactor_%i", 5)); //22e1/old efficiency
	histoMBCorr[m]->SetName(Form("fHistCorrFactorValue_%i", m));
	histoMBCorrForError[m] = (TH1F*)fileMBCorr->Get(Form("fHistCorrFactor_%i", 0)); //uncertainty from 0-5% (the largest one)
	histoMBCorrForError[m]->SetName(Form("fHistCorrFactorForError_%i", m));
      }
      else if (isppHM){
	histoMBCorr[m] = (TH1F*)fileMBCorr->Get(Form("fHistCorrFactor_%i", 0)); //22e1/old efficiency
	histoMBCorr[m]->SetName(Form("fHistCorrFactorValue_%i", m));
	histoMBCorrForError[m] = (TH1F*)fileMBCorr->Get(Form("fHistCorrFactor_%i", 0)); //22e1/old efficiency
	histoMBCorrForError[m]->SetName(Form("fHistCorrFactorForError_%i", m));
      }
      else {
	histoMBCorr[m] = (TH1F*)fileMBCorr->Get(Form("fHistCorrFactor_%i", m)); //22e1/old efficiency
	histoMBCorr[m]->SetName(Form("fHistCorrFactorValue_%i", m));
	histoMBCorrForError[m] = (TH1F*)fileMBCorr->Get(Form("fHistCorrFactor_%i", m)); //22e1/old efficiency
	histoMBCorrForError[m]->SetName(Form("fHistCorrFactorForError_%i", m));
      }

      if (type==8){
	for (Int_t b=1; b<= histoMBCorr[m]->GetNbinsX(); b++){
	  if (histoMBCorr[m]->GetBinContent(b)!=0){
	    histoMBCorr[m]->SetBinContent(b, 1.02); //FIXED VALUES
	    histoMBCorr[m]->SetBinError(b, 0); 
	  }
	}
      }
    }
  }

  TString titleEff = "";

  //Summing errors altogether********************

  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      Sigma2OOBPU[m] =pow(OOBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
      Sigma2IBPU[m] =pow(IBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;

      if (type==0) {
	if (MaterialBudgetCorr!=0){
	  MB = TMath::Abs(histoMBCorrForError[m]->GetBinContent(v)-1)/2;
	}
	else MB = MBPtK0s[v];
      }
      else{
	if (MaterialBudgetCorr!=0)	MB = 0.02; //fixed 2%
	else  MB = MBPtXi[v];
      }

      if (isMC && !isEfficiency) MB = 0;

      Sigma2MB[m]= pow(MB *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;

      /*
      if (!isMC && ispp5TeV && type==0 && MaterialBudgetCorr==2){
	Sigma2MB[m] += pow(0.02*fHistSpectrumStat[m]->GetBinContent(v+1),2); //to take into account multiplicity dependence of correction factor
	MBCorrection = 0.02; 
      }
      else if (!isMC && type==8 && MaterialBudgetCorr==2){
	Sigma2MB[m] += pow(0.02*fHistSpectrumStat[m]->GetBinContent(v+1),2); //to take into account multiplicity dependence of correction factor
	MBCorrection = 0.02; 
      }
      */

      fHistSpectrumSistAll[m]->SetBinError(v+1,sqrt(pow(fHistSpectrumSistDPhi[m]->GetBinError(v+1),2) + pow(fHistSpectrumSist[m]->GetBinError(v+1),2) + Sigma2OOBPU[m]+ Sigma2IBPU[m] + Sigma2MB[m] + pow(fHistSpectrumSistLeadTrack[m]->GetBinError(v+1),2)+  pow(fHistSpectrumSistMCChoice[m]->GetBinError(v+1),2)+ pow(fHistSpectrumSistFakeSB[m]->GetBinError(v+1),2) + pow(fHistSpectrumSistOOJSubDef[m]->GetBinError(v+1),2)));

      if (TypeAnalysis==2) {
	fHistSpectrumSistMultUnCorr[m]->SetBinError(v+1,fHistSpectrumSistRelErrorSE[m]->GetBinContent(v+1) * fHistSpectrumSist[m]->GetBinContent(v+1)*hFractionUnCorrTopoSel[m]->GetBinContent(v+1));
      }
      else {
	if (pow(fHistSpectrumSistDPhi[m]->GetBinError(v+1)*hFractionUnCorrDPhi[m]->GetBinContent(v+1),2) + pow(fHistSpectrumSistRelErrorSE[m]->GetBinContent(v+1) * fHistSpectrumSist[m]->GetBinContent(v+1)*hFractionUnCorrTopoSel[m]->GetBinContent(v+1),2) > 0 ) {
	  fHistSpectrumSistMultUnCorr[m]->SetBinError(v+1,sqrt(pow(fHistSpectrumSistDPhi[m]->GetBinError(v+1)*hFractionUnCorrDPhi[m]->GetBinContent(v+1),2) + pow(fHistSpectrumSistRelErrorSE[m]->GetBinContent(v+1) * fHistSpectrumSist[m]->GetBinContent(v+1)*hFractionUnCorrTopoSel[m]->GetBinContent(v+1),2)));
	}
	else       fHistSpectrumSistMultUnCorr[m]->SetBinError(v+1,0);
      }

      cout << "\nMult: " << Smolt[m] << " pt: " <<  fHistSpectrumSistAll[m]->GetBinCenter(v+1) << " uncorrelated fraction on the total: " <<  fHistSpectrumSistMultUnCorr[m]->GetBinError(v+1) / fHistSpectrumSistAll[m]->GetBinError(v+1) << endl;
      cout << "Topo sel uncorr fraction: " << hFractionUnCorrTopoSel[m]->GetBinContent(v+1) << endl;
      cout << "Fraction of syt. related to topo sel: " << fHistSpectrumSistRelErrorSE[m]->GetBinContent(v+1) *  fHistSpectrumSistAll[m]->GetBinContent(v+1) / fHistSpectrumSistAll[m]->GetBinError(v+1)<< endl;
      if (TypeAnalysis!=2){
	cout << "DPhi uncorr fraction: " << hFractionUnCorrDPhi[m]->GetBinContent(v+1) << endl;
	cout << "Fraction of syt. related to dPhi choice: " << fHistSpectrumSistDPhi[m]->GetBinError(v+1) / fHistSpectrumSistAll[m]->GetBinError(v+1)<< endl;
      }

      if (fHistSpectrumSistAll[m]->GetBinContent(v+1) == 0) fHistSpectrumSistAll[m]->SetBinError(v+1, 0);
      if (fHistSpectrumSistMultUnCorr[m]->GetBinContent(v+1) == 0) fHistSpectrumSistMultUnCorr[m]->SetBinError(v+1, 0);

      if (fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSistAll[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumSistRelErrorMultUnCorr[m]->SetBinContent(v+1, fHistSpectrumSistMultUnCorr[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumSistRelErrorMB[m]->SetBinContent(v+1,MB);
	fHistSpectrumSistRelErrorMBCorrection[m]->SetBinContent(v+1,MBCorrection);
	fHistSpectrumSistRelErrorIBPU[m]->SetBinContent(v+1,IBPU);
	fHistSpectrumSistRelErrorOOBPU[m]->SetBinContent(v+1,OOBPU);
	fHistSpectrumSistRelErrorLeadTrack[m]->SetBinContent(v+1,fHistSpectrumSistLeadTrack[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
       	fHistSpectrumSistRelErrorMCChoice[m]->SetBinContent(v+1,fHistSpectrumSistMCChoice[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumSistRelErrorFakeSB[m]->SetBinContent(v+1,fHistSpectrumSistFakeSB[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
       	fHistSpectrumSistRelErrorOOJSubDef[m]->SetBinContent(v+1,fHistSpectrumSistOOJSubDef[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
      }

      fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMultUnCorr[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMB[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMBCorrection[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorDPhi[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorIBPU[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorOOBPU[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorLeadTrack[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorMCChoice[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorFakeSB[m]->SetBinError(v+1,0);
      fHistSpectrumSistRelErrorOOJSubDef[m]->SetBinError(v+1,0);
      //      cout << " stat rel error " <<       fHistSpectrumStatRelError[m]->GetBinContent(v+1)<< endl;
      //      cout << " sist rel error assoc to DeltaPhi choice " <<       fHistSpectrumSistRelErrorDPhi[m]->GetBinContent(v+1)<< endl;
    }
  
    //SMOOTH (just for K0s MB 13 TeV and just for visualisation purposes to avoid drops in uncertainty)
    if (type==0 && !isppHM && !ispp5TeV){
      fHistSpectrumSistRelErrorMB[m]->GetXaxis()->SetRangeUser(2, 8);
      fHistSpectrumSistRelErrorMB[m]->Smooth(1, "R");
      fHistSpectrumSistRelErrorMB[m]->GetXaxis()->SetRangeUser(0, 8);
    }

    //  cout << "Drawing sist rel error " << endl;
    if (isppHM)     canvasPtSpectraRelErrorAll->cd(m+1-2);
    else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraRelErrorAll->cd(3);
    else     canvasPtSpectraRelErrorAll->cd(m+1);
    gPad->SetLeftMargin(0.15);

    if (isOnlyPlottingRelError) titleEff = "";
    else titleEff = title+SmoltLegend[m];

    StyleHisto(fHistSpectrumStatRelError[m], LimInfError, LimSupError, Color[TypeAnalysis], 33, titleX, titleYRel ,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumStatRelError[m]->SetMarkerSize(2);
    StyleHisto(fHistSpectrumSistRelErrorAll[m], LimInfError, LimSupError, Color[TypeAnalysis], 27, titleX, titleYRel ,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorAll[m]->SetMarkerSize(2);
    StyleHisto(fHistSpectrumSistRelErrorMultUnCorr[m], LimInfError, LimSupError, Color[TypeAnalysis], 30, titleX, titleYRel ,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorMultUnCorr[m]->SetMarkerSize(2);
    StyleHisto(fHistSpectrumSistRelErrorSE[m], LimInfError, LimSupError, 807, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorSE[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorDCAz[m], LimInfError, LimSupError, 867, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorDCAz[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorPurity[m], LimInfError, LimSupError, 834, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorPurity[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorDeltaEta[m], LimInfError, LimSupError, 881, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorDeltaEta[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorDPhi[m], LimInfError, LimSupError, 922, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorDPhi[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorMB[m], LimInfError, LimSupError, 401, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorMB[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorMBCorrection[m], LimInfError, LimSupError, 401, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorMBCorrection[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorOOBPU[m], LimInfError, LimSupError, 825, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorOOBPU[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorIBPU[m], LimInfError, LimSupError, 631, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorIBPU[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorLeadTrack[m], LimInfError, LimSupError, 630, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    StyleHisto(fHistSpectrumSistRelErrorMCChoice[m], LimInfError, LimSupError, 907, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorMCChoice[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorFakeSB[m], LimInfError, LimSupError, 600, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorFakeSB[m]->SetLineWidth(2);
    StyleHisto(fHistSpectrumSistRelErrorOOJSubDef[m], LimInfError, LimSupError, 907, 27, titleX, titleYRel,  titleEff,1, LowRangeSpectrumPart[m], UpRangeSpectrumPart[m]);
    fHistSpectrumSistRelErrorOOJSubDef[m]->SetLineWidth(2);
    gPad->SetLogy();

    TLegend *legendErrorAll= new TLegend(0.6, 0.1, 0.9, 0.4);
    legendErrorAll->AddEntry(fHistSpectrumStatRelError[m], "stat.", "pl");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMultUnCorr[m], "syst. mult. uncorr.", "pl");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorSE[m], "syst. signal extr.", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDCAz[m], "syst. DCAzTrigger", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorPurity[m], "syst. purity", "l");
    if(TypeAnalysis!=2)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDeltaEta[m], "syst. #Delta#eta region", "l");
    if(TypeAnalysis!=2)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDPhi[m], "syst. #Delta#phi region", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMB[m], "Material Budget", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorOOBPU[m], "Out-of-bunch PU", "l");
    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorIBPU[m], "In-bunch PU", "l");
    if(BarlowPassedLeadTrackDef[m])   legendErrorAll->AddEntry(fHistSpectrumSistRelErrorLeadTrack[m], "p_{T} < p_{T}^{Trigg}", "l");
    //    if(BarlowPassedMCChoiceDef[m])  
    if (type==0)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMCChoice[m], "MCChoice", "l");
    if (type==8 && !isppHM)    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorFakeSB[m], "Fake "+tipo[type], "l");
    if (BarlowPassedOOJSubDef[m])    legendErrorAll->AddEntry(fHistSpectrumSistRelErrorOOJSubDef[m], "OOJ sub", "l");

    fHistSpectrumSistRelErrorAll[m]->DrawClone("same p");
    fHistSpectrumSistRelErrorMultUnCorr[m]->DrawClone("same p");
    fHistSpectrumSistRelErrorMB[m]->DrawClone("same");
    fHistSpectrumSistRelErrorMBCorrection[m]->DrawClone("same");
    fHistSpectrumSistRelErrorIBPU[m]->DrawClone("same");
    fHistSpectrumSistRelErrorOOBPU[m]->DrawClone("same");
    fHistSpectrumSistRelErrorSE[m]->DrawClone("same");
    fHistSpectrumSistRelErrorDCAz[m]->DrawClone("same");
    fHistSpectrumSistRelErrorPurity[m]->DrawClone("same");
    fHistSpectrumSistRelErrorLeadTrack[m]->DrawClone("same");
    if (type==0)    fHistSpectrumSistRelErrorMCChoice[m]->DrawClone("same");
    fHistSpectrumSistRelErrorFakeSB[m]->DrawClone("same");
    fHistSpectrumSistRelErrorOOJSubDef[m]->DrawClone("same");
    if (TypeAnalysis!=2) {
      fHistSpectrumSistRelErrorDeltaEta[m]->DrawClone("same");
      fHistSpectrumSistRelErrorDPhi[m]->DrawClone("same");
    }
    fHistSpectrumStatRelError[m]->Draw("same p");
    legendErrorAll->Draw("");
    cout << "Hola ! " << endl;
  }
    
  if (MaterialBudgetCorr==2){
    cout << "Correction by efficiency ratio between 22e1 and old efficiency; The MC 22e1 is characterised by a R-dependent material budget" << endl;
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      if (m!=0 && !isppHM) continue; //choose the multiplicity for diagnosis
      if (isppHM && m!=5) continue;
      cout << "\n*******************************************"<<endl;
      cout << "\nStat and syst uncertainty BEFORE correction by efficiency scaling,  (mult "<< m << ")" << endl;
      for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
	cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " (stat. rel. error: " << fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ", syst. rel. error:  "<< fHistSpectrumSistAll[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << " )" << endl;
	//	cout << fHistSpectrumSistMultUnCorr[m]->GetBinError(b+1) << " " << fHistSpectrumSistAll[m]->GetBinError(b+1) << endl;
      }
    }

    Int_t CorrespondingBin=0;
    Float_t RelError[nummoltMax+1] = {0};
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;

      for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ //I only have to *scale* sist error
	//	cout <<	fHistSpectrumSistAll[m]->GetBinLowEdge(b) << " vs " << histoMBCorr[m]->GetBinLowEdge(histoMBCorr[m]->FindBin(fHistSpectrumStat[m]->GetBinCenter(b)))<< endl;
	CorrespondingBin = histoMBCorr[m]->FindBin(fHistSpectrumStat[m]->GetBinCenter(b));
	RelError[m] = sqrt(pow(fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b), 2) + pow(histoMBCorr[m]->GetBinError(CorrespondingBin)/histoMBCorr[m]->GetBinContent(CorrespondingBin),2));

	fHistSpectrumStat[m]->SetBinContent(b,fHistSpectrumStat[m]->GetBinContent(b)/histoMBCorr[m]->GetBinContent(CorrespondingBin));
	fHistSpectrumStat[m]->SetBinError(b, RelError[m]*fHistSpectrumStat[m]->GetBinContent(b));

	fHistSpectrumSistAll[m]->SetBinContent(b,fHistSpectrumSistAll[m]->GetBinContent(b)/histoMBCorr[m]->GetBinContent(CorrespondingBin));
	fHistSpectrumSistAll[m]->SetBinError(b,fHistSpectrumSistAll[m]->GetBinError(b)/histoMBCorr[m]->GetBinContent(CorrespondingBin));
	fHistSpectrumSistMultUnCorr[m]->SetBinContent(b,fHistSpectrumSistMultUnCorr[m]->GetBinContent(b)/histoMBCorr[m]->GetBinContent(CorrespondingBin));
	fHistSpectrumSistMultUnCorr[m]->SetBinError(b,fHistSpectrumSistMultUnCorr[m]->GetBinError(b)/histoMBCorr[m]->GetBinContent(CorrespondingBin));

      }
    }
  }
  //diagnostic cout:
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (m!=0 && !isppHM) continue; //choose the multiplicity for diagnosis
    if (isppHM && m!=5) continue;
    cout << "\n*******************************************"<<endl;
    cout << "\nStat and syst uncertainty after correction by efficiency scaling,  (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      //      cout << "sist content " << fHistSpectrumSistMultUnCorr[m]->GetBinContent(b) << endl;
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " (stat. rel. error: " << fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ", syst. rel. error:  "<< fHistSpectrumSistAll[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << " ) [mult uncorr fraction = " << fHistSpectrumSistMultUnCorr[m]->GetBinError(b)/fHistSpectrumSistAll[m]->GetBinError(b) << "]" << endl;
      //      cout << fHistSpectrumSistMultUnCorr[m]->GetBinError(b+1) << " " << fHistSpectrumSistAll[m]->GetBinError(b+1) << endl;
    }
  }

  cout << "\n //*************spectra normalization ******************" << endl;

  TFile* fileNormCorrFC;
  TH1F * fHistNormCorrFC[nummolt+1];
  TH1F * fHistNormCorrAllMultFC[nummolt+1];
  TF1 * pol0NormFactor[nummolt+1];
  Float_t pol0RelError=0;
  Float_t SpectrumStatValueBeforeScaling[nummolt+1][numPtV0] = {0};
  Float_t SpectrumStatErrorBeforeScaling[nummolt+1][numPtV0] = {0};

  if (isNormCorrFullyComputed==1){ // In this case I perform a correction with the *comnplete* normalisation factor.
    cout << "Normalization factor for comparison with event generators ---- fully computed!" << endl;
    cout << "From file: " << SfileNormCorrFC << endl;
    fileNormCorrFC=new TFile(SfileNormCorrFC,"");
    if (!fileNormCorrFC) {cout << "file norm not there " << endl; return;}
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;

      fHistSpectrumStatNotNorm[m] = (TH1F*) fHistSpectrumStat[m]->Clone(""+ Smolt[m] + "_NotNorm");
      if (isppHM)     canvasNormFactor->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasNormFactor->cd(3);
      else     canvasNormFactor->cd(m+1);

      if (TypeAnalysis==3) RegionTypeOld[TypeAnalysis] = "BulkBlue";

      if (isppHM)      fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",0));
      else if (ispp5TeV) {
	if (type==8 || (type==0 && TypeAnalysis!=2)) fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",5));
	else 	fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));
      }
      else fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));

      if (TypeAnalysis==3) RegionTypeOld[TypeAnalysis] = "Bulk";
      pol0NormFactor[m] = new TF1(Form("pol0NormFactor_m%i", m), "pol0", 0, 8);
      if (TypeAnalysis == 0 && type==8) {
	cout << "pol0 fit to normalization factor (only for Xi in jet)" << endl;
	fHistNormCorrFC[m]-> Fit(pol0NormFactor[m], "R+");
      }

      fHistNormCorrFC[m]-> Draw("");
      if (TypeAnalysis == 0 && type==8)  {
	for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ 
	  SpectrumStatValueBeforeScaling[m][b-1] = 	fHistSpectrumStat[m]->GetBinContent(b);
	  SpectrumStatErrorBeforeScaling[m][b-1] = 	fHistSpectrumStat[m]->GetBinError(b);
	}
 	fHistSpectrumStat[m]->Scale(1./pol0NormFactor[m]->GetParameter(0));
	for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ 
	  pol0RelError = pol0NormFactor[m]->GetParError(0)/pol0NormFactor[m]->GetParameter(0);
	  fHistSpectrumStat[m]->SetBinError(b, fHistSpectrumStat[m]->GetBinContent(b) * sqrt(pow(SpectrumStatErrorBeforeScaling[m][b-1]/SpectrumStatValueBeforeScaling[m][b-1], 2) + pow(pol0RelError, 2)));
	  cout << "pol0Rel error " << pol0RelError << endl;
	  cout << fHistSpectrumStat[m]->GetBinError(b) / fHistSpectrumStat[m]->GetBinContent(b) << endl;
	}
      }
      else {
	if (m==4){
	  fHistNormCorrFC[m]->SetBinContent(fHistNormCorrFC[m]->FindBin(6), fHistNormCorrFC[3]->GetBinContent(fHistNormCorrFC[m]->FindBin(6))); //this bin is empty for 50-100% norm factor
	}
	fHistSpectrumStat[m]->Divide(fHistNormCorrFC[m]);
      }
      
      for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ //I only have to *scale* sist error
	//	cout <<       fHistNormCorrFC[m]->GetBinContent(b) << endl;
	if (TypeAnalysis == 0 && type==8)  {
	  fHistSpectrumSistAll[m]->SetBinContent(b,fHistSpectrumSistAll[m]->GetBinContent(b)/pol0NormFactor[m]->GetParameter(0));
	  fHistSpectrumSistAll[m]->SetBinError(b,fHistSpectrumSistAll[m]->GetBinError(b)/pol0NormFactor[m]->GetParameter(0));
	  fHistSpectrumSistMultUnCorr[m]->SetBinContent(b,fHistSpectrumSistMultUnCorr[m]->GetBinContent(b)/pol0NormFactor[m]->GetParameter(0));
	  fHistSpectrumSistMultUnCorr[m]->SetBinError(b,fHistSpectrumSistMultUnCorr[m]->GetBinError(b)/pol0NormFactor[m]->GetParameter(0));
	}
	else {
	  fHistSpectrumSistAll[m]->SetBinContent(b,fHistSpectrumSistAll[m]->GetBinContent(b)/fHistNormCorrFC[m]->GetBinContent(b));
	  fHistSpectrumSistAll[m]->SetBinError(b,fHistSpectrumSistAll[m]->GetBinError(b)/fHistNormCorrFC[m]->GetBinContent(b));
	  fHistSpectrumSistMultUnCorr[m]->SetBinContent(b,fHistSpectrumSistMultUnCorr[m]->GetBinContent(b)/fHistNormCorrFC[m]->GetBinContent(b));
	  fHistSpectrumSistMultUnCorr[m]->SetBinError(b,fHistSpectrumSistMultUnCorr[m]->GetBinError(b)/fHistNormCorrFC[m]->GetBinContent(b));

	  /*
	  if (fHistNormCorrFC[m]->GetBinContent(b) ==0 ) { //this should be avoided if possible
	    fHistSpectrumSistAll[m]->SetBinContent(b,0);
	    fHistSpectrumSistAll[m]->SetBinError(b,0);
	  }
	  */
	}
      }
    }
  }

  cout << "...end of spectra normalization *********************** \n" << endl;

  //diagnostic cout:
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (m!=0 && !isppHM) continue; //choose the multiplicity for diagnosis
    if (isppHM && m!=5) continue;
    cout << "\n*******************************************"<<endl;
    cout << "\nStat and syst uncertainty after normalization factor correction,  (mult "<< m << ")" << endl;
    for (Int_t b=PtV0Min; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      //      cout << "sist content " << fHistSpectrumSistMultUnCorr[m]->GetBinContent(b) << endl;
      cout << NPtV0[b-1]<<"-"<<NPtV0[b] << ":    " << fHistSpectrumStat[m]->GetBinContent(b) << " (stat. rel. error: " << fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << ", syst. rel. error:  "<< fHistSpectrumSistAll[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b) << " ) [mult uncorr fraction = " << fHistSpectrumSistMultUnCorr[m]->GetBinError(b)/fHistSpectrumSistAll[m]->GetBinError(b) << "]" << endl; 
      //      cout << fHistSpectrumSistMultUnCorr[m]->GetBinError(b+1) << " " << fHistSpectrumSistAll[m]->GetBinError(b+1) << endl;
    }
  }

  cout << "UPDATE of STAT REL ERROR " << endl;
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v <    fHistSpectrumSistAll[m]->GetNbinsX() ; v++){
      if (fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	fHistSpectrumStatRelError[m]->SetBinContent(v+1, fHistSpectrumStat[m]->GetBinError(v+1)/ fHistSpectrumStat[m]->GetBinContent(v+1));
      }
    }
  }
  cout << "END OF UPDATE of STAT REL ERROR " << endl;

  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;

    if (isppHM)     canvasPtSpectra->cd(m+1-2);
    else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectra->cd(3);
    else     canvasPtSpectra->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistAll[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[TypeAnalysis], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
    fHistSpectrumStat[m]->DrawClone("same");
    fHistSpectrumSistAll[m]->SetFillStyle(0);
    fHistSpectrumSistAll[m]->DrawClone("same e2");
    if (m==nummoltMax){
      legendError2->AddEntry(fHistSpectrumStat[m], "stat.", "ple");
      legendError2->AddEntry(fHistSpectrumSistAll[m], "syst.", "fe");
    }
    legendError2->Draw("");
    cout << "\n\e[32mMultiplicity: " << SmoltLegend[m]  << "\e[39m" << endl;
    for (Int_t b= PtBinMin[m] ; b< fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << " " << NPtV0[b] << " < pt < " << NPtV0[b+1] << " GeV/c , Yield: " << fHistSpectrumStat[m]->GetBinContent(b+1) << ", rel. stat. error: " << fHistSpectrumStat[m]->GetBinError(b+1)/ fHistSpectrumStat[m]->GetBinContent(b+1) << endl;
      cout << " " << NPtV0[b] << "< pt < " << NPtV0[b+1] << " GeV/c, Yield: " << fHistSpectrumSistAll[m]->GetBinContent(b+1) << ", rel. sist. error: " << fHistSpectrumSistAll[m]->GetBinError(b+1)/ fHistSpectrumSistAll[m]->GetBinContent(b+1) << endl;
      cout << " " << NPtV0[b] << "< pt < " << NPtV0[b+1] << " GeV/c, Yield: " << fHistSpectrumSistMultUnCorr[m]->GetBinContent(b+1) << ", rel. sist. error (MULT UNCORR): " << fHistSpectrumSistMultUnCorr[m]->GetBinError(b+1)/ fHistSpectrumSistMultUnCorr[m]->GetBinContent(b+1) << " = " << fHistSpectrumSistMultUnCorr[m]->GetBinError(b+1)/ fHistSpectrumSistAll[m]->GetBinError(b+1)<< " of the total syst. error " << endl;
      //      cout << fHistSpectrumSistMultUnCorr[m]->GetBinError(b+1) << " " << fHistSpectrumSistAll[m]->GetBinError(b+1) << endl;
    }

    if (isppHM)     canvasPtSpectraRelErrorAll->cd(m+1-2);
    else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraRelErrorAll->cd(3);
    else     canvasPtSpectraRelErrorAll->cd(m+1);

    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistRelErrorAll[m], LimInfError, LimSupError,Color[TypeAnalysis] , 27, titleX, titleYRel, titleEff, 0, 0, 0);
    StyleHisto(fHistSpectrumSistRelErrorMultUnCorr[m], LimInfError, LimSupError,Color[TypeAnalysis] , 30, titleX, titleYRel, titleEff, 0, 0, 0);
    StyleHisto(fHistSpectrumStatRelError[m], LimInfError, LimSupError,Color[TypeAnalysis] , 33, titleX, titleYRel, titleEff, 0, 0, 0);
    StyleHisto(fHistSpectrumSistRelErrorPublished, LimInfError, LimSupError, 922, 33, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);

    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumSistRelErrorMultUnCorr[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");

    if (isppHM)     canvasPtSpectraRelError->cd(m+1-2);
    else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraRelError->cd(3);
    else     canvasPtSpectraRelError->cd(m+1);
    gPad->SetLeftMargin(0.15);
    if (m==nummoltMax) fHistSpectrumSistRelErrorPublished->Draw("same");
    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumSistRelErrorMultUnCorr[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");
    if (m==nummoltMax){
      legendError->AddEntry(    fHistSpectrumStatRelError[m], "stat.", "pl");
      legendError->AddEntry(    fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
      legendError->AddEntry(    fHistSpectrumSistRelErrorMultUnCorr[m], "syst. mult. uncorr", "pl");
    }
    if (m==nummoltMax)      legendError->AddEntry(fHistSpectrumSistRelErrorPublished, "syst. published", "l");
    legendError->Draw("");
  } //end loop m

  cout << "Drawing plot rel error for thesis"  << endl;
  //PLOTS
  TCanvas * canvasNR = new TCanvas("canvasNR", "canvasNR", 1000, 1200); //only one pad
  canvasNR->SetFillColor(0);
  canvasNR->SetTickx(1);
  canvasNR->SetTicky(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  canvasNR->cd();
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.02);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLogy();
  Int_t multChosen = 5;

  //*********************
  fHistSpectrumSistRelErrorDCAz[multChosen]->Smooth(1, "R");
  if (TypeAnalysis!=2) {
    //    if (type==8 && TypeAnalysis==0)    fHistSpectrumSistRelErrorDeltaEta[multChosen]->Smooth(1, "R");
    //if (type==0)     fHistSpectrumSistRelErrorDeltaEta[multChosen]->Smooth(1, "R");
    if (type==0) fHistSpectrumSistRelErrorDPhi[multChosen]->Smooth(1, "R");
  }
  fHistSpectrumSistRelErrorMB[multChosen]->Smooth(1, "R");

  fHistSpectrumSistRelErrorAll[multChosen]->Draw("same p");
  if (multChosen!=5)  fHistSpectrumSistRelErrorMultUnCorr[multChosen]->Draw("same p"); //it is set to zero in the multiplicity class 0-100%
  fHistSpectrumStatRelError[multChosen]->DrawClone("same pl");
  fHistSpectrumSistRelErrorMB[multChosen]->Draw("same");
  fHistSpectrumSistRelErrorIBPU[multChosen]->Draw("same");
  fHistSpectrumSistRelErrorOOBPU[multChosen]->Draw("same");
  fHistSpectrumSistRelErrorSE[multChosen]->Draw("same");
  fHistSpectrumSistRelErrorDCAz[multChosen]->Draw("same");
  fHistSpectrumSistRelErrorPurity[multChosen]->Draw("same");
  if (type==0)    fHistSpectrumSistRelErrorMCChoice[multChosen]->Draw("same");
  if (type==8 && !isppHM)  fHistSpectrumSistRelErrorFakeSB[multChosen]->Draw("same");
  //  if (type==8)  fHistSpectrumSistRelErrorOOJSubDef[multChosen]->Draw("same");
  if (TypeAnalysis!=2) {
    if (type==8 && TypeAnalysis==0)    fHistSpectrumSistRelErrorDeltaEta[multChosen]->Draw("same");
    if (type==0) fHistSpectrumSistRelErrorDeltaEta[multChosen]->Draw("same");
    if (type==0) fHistSpectrumSistRelErrorDPhi[multChosen]->Draw("same");
    if (type==8) {
      if (!isppHM && !ispp5TeV) fHistSpectrumSistRelErrorDPhi[0]->Draw("same");
      else fHistSpectrumSistRelErrorDPhi[multChosen]->Draw("same");
    }
  }

  TString NameP1[2]={"K_{S}^{0}", "#Xi"};
  TLegend *legendAllErrors;
  if (type==0) legendAllErrors= new TLegend(0.6, 0.21, 0.9, 0.46);
  else legendAllErrors= new TLegend(0.6, 0.17, 0.9, 0.42);
  legendAllErrors->SetFillStyle(0);
  legendAllErrors->SetTextAlign(12);
  legendAllErrors->SetTextSize(0.027);
  legendAllErrors->AddEntry(fHistSpectrumStatRelError[multChosen], "stat. error", "pl");
  legendAllErrors->AddEntry(fHistSpectrumSistRelErrorAll[multChosen], "syst. error", "pl");
  if (multChosen!=5) legendAllErrors->AddEntry(fHistSpectrumSistRelErrorMultUnCorr[multChosen], "syst. error (mult. uncorr.)", "pl");
  if (type==0)  legendAllErrors->AddEntry(fHistSpectrumSistRelErrorSE[multChosen], NameP1[type]+" candidates sel.", "l");
  else   legendAllErrors->AddEntry(fHistSpectrumSistRelErrorSE[multChosen], NameP1[1]+" candidates sel.", "l");
  legendAllErrors->AddEntry(fHistSpectrumSistRelErrorDCAz[multChosen], "DCAz trigger", "l");
  legendAllErrors->AddEntry(fHistSpectrumSistRelErrorPurity[multChosen], "signal extr.", "l");
  legendAllErrors->AddEntry(fHistSpectrumSistRelErrorMB[multChosen], "material budget", "l");
  legendAllErrors->AddEntry(fHistSpectrumSistRelErrorOOBPU[multChosen], "out-of-bunch pile-up", "l");
  legendAllErrors->AddEntry(fHistSpectrumSistRelErrorIBPU[multChosen], "in-bunch pile-up", "l");
  if (type==0)    legendAllErrors->AddEntry(fHistSpectrumSistRelErrorMCChoice[multChosen], "Monte Carlo", "l");
  if (type==8 && !isppHM)     legendAllErrors->AddEntry(fHistSpectrumSistRelErrorFakeSB[multChosen], "fake "+NameP1[1], "l");
  if (BarlowPassedOOJSubDef[multChosen])    legendAllErrors->AddEntry(fHistSpectrumSistRelErrorOOJSubDef[multChosen], "out-of-jet subtr.", "l");

  TLegend *legendDeltaErrors= new TLegend(0.7, 0.81, 1, 0.96);
  legendDeltaErrors->SetFillStyle(0);
  legendDeltaErrors->SetTextAlign(12);
  legendDeltaErrors->SetTextSize(0.027);
  legendDeltaErrors->AddEntry(fHistSpectrumStatRelError[multChosen], "stat. error", "pl");
  legendDeltaErrors->AddEntry(fHistSpectrumSistRelErrorAll[multChosen], "syst. error", "pl");
  if (multChosen!=5)  legendDeltaErrors->AddEntry(fHistSpectrumSistRelErrorMultUnCorr[multChosen], "syst. error (mult. uncorr.)", "pl");
  if (type==0){
    if(TypeAnalysis!=2)    legendDeltaErrors->AddEntry(fHistSpectrumSistRelErrorDeltaEta[multChosen], "#Delta#eta region", "l");
    if(TypeAnalysis!=2)    legendDeltaErrors->AddEntry(fHistSpectrumSistRelErrorDPhi[multChosen], "#Delta#varphi region", "l");
  }
  else {
    if(TypeAnalysis==0)    legendDeltaErrors->AddEntry(fHistSpectrumSistRelErrorDeltaEta[multChosen], "#Delta#eta region", "l");
  }

  if (TypeAnalysis==2)  legendAllErrors->Draw("");
  //  else   if (TypeAnalysis==1)  legendDeltaErrors->Draw("");
  else legendDeltaErrors->Draw("");

  TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
  TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
  TString sRegion1[3]={"|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 1.1", "0.86 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
  //  TString sRegion1[3]={" #color[628]{:|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 1.1}", " #color[418]{:0.86 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8}", " #color[600]{:|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2}"};
  TLegend *legendRegionF;
  if (TypeAnalysis==2) legendRegionF=new TLegend(0.08, 0.74, 0.7, 0.81);
  else legendRegionF=new TLegend(0.08, 0.88, 0.7, 0.95);
  legendRegionF->SetFillStyle(0);
  TLegendEntry * lRe1 = legendRegionF->AddEntry("", sRegion[TypeAnalysis], "");
  lRe1->SetTextSize(0.034);
  lRe1->SetTextAlign(12);
  TLegendEntry *leR2=      legendRegionF->AddEntry("", sRegion1[TypeAnalysis], "");
  leR2->SetTextSize(0.03);
  leR2->SetTextAlign(12);
  legendRegionF->Draw("");

  TLegend *Legend1=new TLegend(0.08,0.84,0.7,0.96);
  Legend1->SetFillStyle(0);
  Legend1->SetTextAlign(12);
  Legend1->SetTextSize(0.034);
  Legend1->AddEntry("", "#bf{This work}", "");
  if (ispp5TeV)      Legend1->AddEntry("", "pp, #sqrt{#it{s}} = 5.02 TeV", "");
  else 	Legend1->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  if (type==0)  Legend1->AddEntry("", NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  else   Legend1->AddEntry("", NameP[1]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
  //  Legend1->AddEntry("", "V0M Multiplicity Percentile "+Smolt[multChosen], "");
  Legend1->AddEntry("", "V0M Multiplicity Percentile 0-100%", "");
  //  Legend1->AddEntry("", sRegion[TypeAnalysis] + sRegion1[TypeAnalysis], "");
  if (TypeAnalysis==2)  Legend1->Draw("");

  canvasNR->SaveAs("RelError"+year+"_"+RegionTypeNew[TypeAnalysis]+".pdf");
  TFile * temp = new TFile("temp.root", "RECREATE");
  canvasNR->Write();
  temp->Close();
  if (isOnlyPlottingRelError) return;

  cout << "(end of)Drawing plot rel error for thesis"  << endl;

  //fourth part: fit to obtain pt-integrated yield vs mult
  AliPWGFunc pwgfunc;
  const Int_t numfittipo=4;
  Int_t numfittipoEff = 0; //effective number of fit function used
  TLegend *legendfit=new TLegend(0.6, 0.6, 0.9, 0.9);
  TString   nameMTscaling[nummolt+1][numfittipo];

  TF1 * rettaUno= new TF1("rettaUno", "pol0",0,8);
  rettaUno->FixParameter(0,1);

  TH1F* fHistSpectrumRatioFit[nummolt+1][numfittipo];
  TH1F* hhoutYield[nummolt+1];
  TH1F* hhoutYieldMy[nummolt+1];
  TH1F* hhoutYieldRatioToMine[nummolt+1];
  TH1F* hhoutRelStat[nummolt+1];
  TH1F* hhoutRelStatMy[nummolt+1];
  TH1F* hhoutRelSystHigh[nummolt+1];
  TH1F* hhoutRelSystHighMy[nummolt+1];
  TH1F* hhoutRelSystLow[nummolt+1];
  TH1F* hhoutRelSystLowMy[nummolt+1];

  TH1F* hhoutAvgPt[nummolt+1];
  TH1F* hhoutAvgPtMy[nummolt+1];
  TH1F* hhoutAvgPtRatioToMine[nummolt+1];
  TH1F* hhoutAvgPtRelStat[nummolt+1];
  TH1F* hhoutAvgPtRelStatMy[nummolt+1];
  TH1F* hhoutAvgPtRelSystHigh[nummolt+1];
  TH1F* hhoutAvgPtRelSystHighMy[nummolt+1];
  TH1F* hhoutAvgPtRelSystLow[nummolt+1];
  TH1F* hhoutAvgPtRelSystLowMy[nummolt+1];

  Float_t hhoutYieldAvg[nummolt+1]={0};
  Float_t hhoutRelStatAvg[nummolt+1]={0};
  Float_t hhoutRelSystHighAvg[nummolt+1]={0};
  Float_t hhoutRelSystLowAvg[nummolt+1]={0};

  Float_t hhoutAvgPtAvg[nummolt+1]={0};
  Float_t hhoutAvgPtRelStatAvg[nummolt+1]={0};
  Float_t hhoutAvgPtRelSystHighAvg[nummolt+1]={0};
  Float_t hhoutAvgPtRelSystLowAvg[nummolt+1]={0};

  TH1 *hhout[numfittipo][nummolt+1];
  TString  Titlehhout[9] = {"kYield",
    "kYieldStat",
    "kYieldSysHi",
    "kYieldSysLo",
    "kMean",
    "kMeanStat",
    "kMeanSysHi",
    "kMeanSysLo",
    "kExtra"};

  Int_t ColorFit[numfittipo+1]={860, 881, 868, 628, 419};
  TFitResultPtr fFitResultPtr0[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrUp[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrDown[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrHard[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrSoft[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtr1[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrTemp[nummolt+1][numfittipo];
  TF1* fit_MTscaling[nummolt+1][numfittipo];
  TF1* fit_MTscalingUp[nummolt+1][numfittipo];
  TF1* fit_MTscalingDown[nummolt+1][numfittipo];
  TF1* fit_MTscalingHard[nummolt+1][numfittipo];
  TF1* fit_MTscalingSoft[nummolt+1][numfittipo];
  TF1* fit_MTscalingBis[nummolt+1][numfittipo];
  TF1* fit_MTscalingTemp[nummolt+1][numfittipo];
  TString nameFit[numfittipo+1]={"mT-scaling", "Boltzmann", "Fermi-Dirac", "Levi", ""};
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  Int_t factor=1;

  Float_t   AvgPt[nummolt+1]={0};
  Float_t   AvgPtFS[nummolt+1]={0};
  Float_t   AvgPtMeanMacro[nummolt+1]={0};

  Float_t   AvgPtDiscrErr[nummolt+1] = {0};
  Float_t   AvgPtFSNum[nummolt+1]={0};
  Float_t   AvgPtFSDenom[nummolt+1]={0};
  Float_t   AvgPtFSTemp[nummolt+1]={0};
  Float_t   AvgPtFSNumTemp[nummolt+1]={0};
  Float_t   AvgPtFSDenomTemp[nummolt+1]={0};
  Float_t   AvgPtFSHard[nummolt+1]={0};
  Float_t   AvgPtFSNumHard[nummolt+1]={0};
  Float_t   AvgPtFSDenomHard[nummolt+1]={0};
  Float_t   AvgPtFSSoft[nummolt+1]={0};
  Float_t   AvgPtFSNumSoft[nummolt+1]={0};
  Float_t   AvgPtFSDenomSoft[nummolt+1]={0};
  Float_t   AvgPtHard[nummolt+1]={0};
  Float_t   AvgPtSoft[nummolt+1]={0};
  Float_t   AvgPtMax[nummolt+1]={0};
  Float_t   AvgPtMin[nummolt+1]={0};
  Float_t   AvgTemp[3] = {0};
  Float_t   AvgPtFit[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitMinFit[nummolt+1]={0};
  Float_t   AvgPtFitMaxFit[nummolt+1]={0};
  Float_t   AvgPtMeanMacroMinFit[nummolt+1]={0};
  Float_t   AvgPtMeanMacroMaxFit[nummolt+1]={0};
  Float_t   AvgPtFitUp[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitDown[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitHard[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitSoft[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitTemp[nummolt+1][numfittipo]={0};

  Float_t   AvgPtSistErr[nummolt+1]={0};
  Float_t   AvgPtStatErr[nummolt+1]={0};
  Float_t   AvgPtSistErrFS[nummolt+1]={0};
  Float_t   AvgPtStatErrFS[nummolt+1]={0};
  Float_t   AvgPtSistErrMeanMacro[nummolt+1]={0};
  Float_t   AvgPtStatErrMeanMacro[nummolt+1]={0};
  Float_t   AvgPtSistErrFit[nummolt+1]={0};
  Float_t   AvgPtSistErrFitMeanMacro[nummolt+1]={0};
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
  Float_t    Yield4FitErrSistLowPt[nummolt+1]={0};
  Float_t    Yield4FitErrSistHighPt[nummolt+1]={0};
  Float_t    YieldExtrErrSistLow[nummolt+1]={0};
  Float_t    YieldErrSystExtrLowPt[nummolt+1]={0};
  Float_t    YieldErrSystExtrHighPt[nummolt+1]={0};
  Float_t    YieldSpectrum[nummolt+1]={0};
  Float_t    YieldSpectrumErrStat[nummolt+1]={0};
  Float_t    YieldSpectrumErrSist[nummolt+1]={0};
  Float_t    YieldSpectrumErrSistMultUnCorr[nummolt+1]={0};
  Float_t    Yield[nummolt+1]={0};
  Float_t    FracErrorUnCorr[nummolt+1]={0};
  Float_t    YieldErrStat[nummolt+1]={0};
  Float_t    YieldErrSist[nummolt+1]={0};
  Float_t    YieldErrSistMultUnCorr[nummolt+1]={0};
  Float_t    YieldErrSistMy[nummolt+1]={0};
  Float_t    YieldErrSistLow[nummolt+1]={0};
  Float_t    YieldErrRelNormFactor[nummolt+1]={0, 0, 0, 0.015, 0.138, 0.08};

  TCanvas* canvasFitResult = new TCanvas ("canvasFitResult", "canvasFitResult", 1300, 800);
  TH1F * hFitResult[numfittipo];
  TH1F * hExtrFractionLowPt[numfittipo+1];
  TH1F * hExtrFractionHighPt[numfittipo+1];
  TH1F*  fHistAvgPtDistr[nummolt+1][numfittipo];
  TH1F*  fHistAvgPtFSDistr[nummolt+1];

  cout << "\n\e[35m*******************************************************\n" << endl;
  cout << "Fitting the pt spectra... " <<endl; 
  cout << "\n*******************************************************\n\e[39m" << endl;
  for (Int_t typefit =0; typefit<numfittipo; typefit++){
    cout << "ciao " << endl;
    hFitResult[typefit] = new TH1F ("hFitResult"+nameFit[typefit], "hFitResult",nummoltMax, Nmolt);
  }
  for (Int_t typefit =0; typefit<numfittipo+1; typefit++){
    hExtrFractionHighPt[typefit] = new TH1F ("hExtrFractionHighPt"+nameFit[typefit], "hExtrFractionHighPt",nummoltMax, Nmolt);
    hExtrFractionLowPt[typefit] = new TH1F ("hExtrFractionLowPt"+nameFit[typefit], "hExtrFractionLowPt",nummoltMax, Nmolt);
  }
  for(Int_t m=0; m<nummoltMax+1; m++){
    numfittipoEff=0;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && m<2) continue;
    cout << "\n\e[35m*** Multiplicity: " << SmoltLegend[m] << " ***\e[39m" << endl;
    hhoutYield[m] = new TH1F(Form("hhoutYield_m%i", m),Form("hhoutYield_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutYieldMy[m] = new TH1F(Form("hhoutYieldMy_m%i", m),Form("hhoutYieldMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelStat[m] = new TH1F(Form("hhoutRelStat_m%i", m),Form("hhoutRelStat_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelStatMy[m] = new TH1F(Form("hhoutRelStatMy_m%i", m),Form("hhoutRelStatMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystHigh[m] = new TH1F(Form("hhoutRelSystHigh_m%i", m),Form("hhoutRelSystHigh_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystHighMy[m] = new TH1F(Form("hhoutRelSystHighMy_m%i", m),Form("hhoutRelSystHighMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystLow[m] = new TH1F(Form("hhoutRelSystLow_m%i", m),Form("hhoutRelSystLow_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystLowMy[m] = new TH1F(Form("hhoutRelSystLowMy_m%i", m),Form("hhoutRelSystLowMy_m%i", m) , numfittipo+1, 0, numfittipo+1);

    hhoutAvgPt[m] = new TH1F(Form("hhoutAvgPt_m%i", m),Form("hhoutAvgPt_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtMy[m] = new TH1F(Form("hhoutAvgPtMy_m%i", m),Form("hhoutAvgPtMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelStat[m] = new TH1F(Form("hhoutAvgPtRelStat_m%i", m),Form("hhoutAvgPtRelStat_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelStatMy[m] = new TH1F(Form("hhoutAvgPtRelStatMy_m%i", m),Form("hhoutAvgPtRelStatMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystHigh[m] = new TH1F(Form("hhoutAvgPtRelSystHigh_m%i", m),Form("hhoutAvgPtRelSystHigh_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystHighMy[m] = new TH1F(Form("hhoutAvgPtRelSystHighMy_m%i", m),Form("hhoutAvgPtRelSystHighMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystLow[m] = new TH1F(Form("hhoutAvgPtRelSystLow_m%i", m),Form("hhoutAvgPtRelSystLow_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystLowMy[m] = new TH1F(Form("hhoutAvgPtRelSystLowMy_m%i", m),Form("hhoutAvgPtRelSystLowMy_m%i", m) , numfittipo+1, 0, numfittipo+1);

    fHistAvgPtFSDistr[m] = new TH1F (Form("fHistAvgPtFSDistr_m%i", m), Form("fHistAvgPtFSDistr_m%i", m), 10000, LimInfPtvsMult, LimSupPtvsMult);

    AvgPtMeanMacroMaxFit[m] = 0;
    AvgPtMeanMacroMinFit[m] = 1000;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
      if (type==8 && TypeAnalysis==10) factor=2; //otherwise fit is not properly done                 
      pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
      nameMTscaling[m][typefit] = Form("fitMTscaling_m%i_fit%i",m, typefit);
      cout << "\n***********Fitting pt spectra with: " << nameFit[typefit]<< " (Quiet mode selected) "<< endl;
      if (typefit==0)      fit_MTscaling[m][typefit]=    pwgfunc.GetMTExp(massParticle[type], 0.1, 0.04*factor, nameMTscaling[m][typefit]); //mass, T, norm, name                                
      if (typefit==1) {
	//0: normalisation
	//1: temperature
	fit_MTscaling[m][typefit]=    pwgfunc.GetBoltzmann(massParticle[type],0.1, 0.04*factor, nameMTscaling[m][typefit]);
	if (MaterialBudgetCorr==2 && type==0 && isppHM && TypeAnalysis!=0) {
	  //	  fit_MTscaling[m][typefit]->SetParLimits(0, 2, 10);
	  //	  fit_MTscaling[m][typefit]->SetParameter(0, 4);
	}
	if (type==8 && TypeAnalysis==0 && isppHM)	fit_MTscaling[m][typefit]->SetParLimits(1, 0.1, 10);
	if (type==8 && TypeAnalysis==0 && isGenOnTheFly && MonashTune==2)	fit_MTscaling[m][typefit]->SetParLimits(1, 0.1, 10);
      }
      if (typefit==2)      fit_MTscaling[m][typefit]=    pwgfunc.GetFermiDirac(massParticle[type],0.1, 0.04*factor, nameMTscaling[m][typefit]);
      if (typefit==3)  {
	fit_MTscaling[m][typefit]=    pwgfunc.GetLevi(massParticle[type],0.1, 0.03, 0.04*factor, nameMTscaling[m][typefit]);  //norm, n, T, mass (but the function must be called with these parameters in inverse order) 
	fit_MTscaling[m][typefit]->SetParLimits(0, 0,fHistSpectrumStat[m]->GetBinContent(fHistSpectrumStat[m]->GetMaximumBin())*0.5*10); //norm
	fit_MTscaling[m][typefit]->SetParLimits(2, 0.1, 10); //T
	fit_MTscaling[m][typefit]->SetParLimits(1, 2, 30); //n
	if (MaterialBudgetCorr==2 && type==0 && isppHM && TypeAnalysis!=0) {
	  fit_MTscaling[m][typefit]->SetParLimits(2, 0.4, 10); //T
	}
	if (type==8) {
	  fit_MTscaling[m][typefit]->SetParLimits(1, 15,30);
	  if (isGenOnTheFly && TypeAnalysis==0) fit_MTscaling[m][typefit]->SetParLimits(1, 2, 100000); 
	}
	fit_MTscaling[m][typefit]->SetParameter(2, 0.7);
	if (ispp5TeV && type==0 && TypeAnalysis==2)   fit_MTscaling[m][typefit]->SetParLimits(2, 0.25, 10);
      }
      if (m==nummoltMax)    legendfit->AddEntry( fit_MTscaling[m][typefit], nameFit[typefit],"l");

      fit_MTscaling[m][typefit]->SetLineColor(ColorFit[typefit]);

      fit_MTscalingBis[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Bis");

      fit_MTscaling[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fit_MTscalingBis[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      if (TwoFitFunctions) {
	fit_MTscaling[m][typefit]->SetRange(LowRange[m], 1.5);
	fit_MTscalingBis[m][typefit]->SetRange(2, UpRange[m]);
      }

      if (isMeanMacro){
	hhout[typefit][m] = YieldMean(fHistSpectrumStat[m], fHistSpectrumSistAll[m],  fit_MTscalingBis[m][typefit], 0, 20, 0.01, 0.1, "0qI", "log.root", LowRange[m], UpRange[m]);
	hhout[typefit][m] -> SetLineColor(ColorFit[typefit]);
	hhout[typefit][m] -> SetName("hhout_"+ nameFit[typefit]+Form("_m%i", m));
	hhout[typefit][m]->GetYaxis()->SetRangeUser(0,2);
	for (Int_t b=1; b<= hhout[typefit][m]->GetNbinsX(); b++){
	  hhout[typefit][m] ->GetXaxis() ->SetBinLabel(b,Titlehhout[b-1] );
	}
	cout << "m " << m << " typefit " << typefit << endl;
	hhoutYield[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutRelStat[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutRelSystHigh[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutRelSystLow[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

	hhoutAvgPt[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutAvgPtRelStat[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutAvgPtRelSystHigh[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutAvgPtRelSystLow[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

	cout << "m " << m << " typefit " << typefit << endl;
	hhoutYield[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(9)); //extrapolated fraction
	hhoutYieldAvg[m] +=       hhoutYield[m] ->GetBinContent(typefit+1); 

	hhoutRelStat[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(2)); // stat error
	hhoutRelStatAvg[m] +=       hhoutRelStat[m] ->GetBinContent(typefit+1); 

	hhoutRelSystHigh[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(3)); // syst error
	hhoutRelSystHighAvg[m] +=       hhoutRelSystHigh[m] ->GetBinContent(typefit+1); 

	hhoutRelSystLow[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(4)); // syst error
	hhoutRelSystLowAvg[m] +=       hhoutRelSystLow[m] ->GetBinContent(typefit+1); 

	hhoutAvgPt[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(5));  //avg pt
	hhoutAvgPtAvg[m] +=       hhoutAvgPt[m] ->GetBinContent(typefit+1); 

	hhoutAvgPtRelStat[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(6)); // stat error of avg pt
	hhoutAvgPtRelStatAvg[m] +=       hhoutAvgPtRelStat[m] ->GetBinContent(typefit+1); 

	hhoutAvgPtRelSystHigh[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(7)); // syst error of avg pt
	hhoutAvgPtRelSystHighAvg[m] +=       hhoutAvgPtRelSystHigh[m] ->GetBinContent(typefit+1); 

	hhoutAvgPtRelSystLow[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(8)); // syst error of avg pt
	hhoutAvgPtRelSystLowAvg[m] +=       hhoutAvgPtRelSystLow[m] ->GetBinContent(typefit+1); 

	//2.C SYSTEMATIC RELATED TO CHOICE OF FIT FUNCTION
	if(AvgPtMeanMacroMaxFit[m] < hhout[typefit][m]->GetBinContent(5)){
	  AvgPtMeanMacroMaxFit[m] = hhout[typefit][m]->GetBinContent(5);
	}
	if(AvgPtMeanMacroMinFit[m] > hhout[typefit][m]->GetBinContent(5)){
	  AvgPtMeanMacroMinFit[m] = hhout[typefit][m]->GetBinContent(5);
	}
	cout <<  AvgPtMeanMacroMaxFit[m] << " " <<   AvgPtMeanMacroMinFit[m]<< endl;

      }

      fFitResultPtr0[m][typefit] = fHistSpectrumStat[m]->Fit(fit_MTscaling[m][typefit],"SR0I");
      fit_MTscaling[m][typefit]->SetRange(0,20);

      if (isppHM)     canvasPtSpectraFit->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraFit->cd(3);
      else     canvasPtSpectraFit->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStat[m]->Draw("same");
      fit_MTscaling[m][typefit]->Draw("same");
      legendfit->Draw("");

      fFitResultPtr1[m][typefit]=       fHistSpectrumStat[m]->Fit(    fit_MTscalingBis[m][typefit],"SR0IQ");
      fit_MTscalingBis[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitBis->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraFitBis->cd(3);
      else     canvasPtSpectraFitBis->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStat[m]->Draw("same");
      fit_MTscalingBis[m][typefit]->Draw("same");
      legendfit->Draw("");

      fHistSpectrumRatioFit[m][typefit] = (TH1F*) fHistSpectrumStat[m]->Clone(Form("HistRatioFit_m%i_typefit%i", m, typefit));
      fHistSpectrumRatioFit[m][typefit]->Divide(fit_MTscaling[m][typefit]);
      if (isppHM)     canvasPtSpectraFitRatio->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraFitRatio->cd(3);
      else     canvasPtSpectraFitRatio->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumRatioFit[m][typefit], 0.4, 2, ColorFit[typefit], 33, "p_{T} (GeV/c)", "", "Spectra/Fit ratio", 0, 0, 0);
      rettaUno->SetLineColor(kBlack);
      fHistSpectrumRatioFit[m][typefit]->Draw("samee");
      if (typefit==0)      rettaUno->Draw("same");
      legendfit->Draw("");
      hFitResult[typefit] ->SetBinContent(m+1,fit_MTscaling[m][typefit]->GetChisquare()/fit_MTscaling[m][typefit]->GetNDF());

      //DEFINE + FIT SPECTRA TO COMPUTE SYST UNCERTAINTY
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumStatUp[m]->SetBinContent(v+1,fHistSpectrumSistAll[m]->GetBinError(v+1)+fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumStatDown[m]->SetBinContent(v+1,-fHistSpectrumSistAll[m]->GetBinError(v+1)+fHistSpectrumStat[m]->GetBinContent(v+1));
      }
      if (isMC && !isEfficiency){
	if (type == 8 && (TypeAnalysis == 2 || TypeAnalysis == 1 || TypeAnalysis == 0)){
	  fHistSpectrumStatHard[m] = (TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumHard_"+Smolt[m]);
	  fHistSpectrumStatSoft[m] = (TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSoft_"+Smolt[m]);
	}
      }
      else { //only when systematic error is non-zero!
	cout << "If I break, it's because the systematic error of the spectra is zero! " << endl;
	fHistSpectrumStatHard[m] = (TH1F*)YieldMean_ReturnExtremeHardHisto(fHistSpectrumSist[m]); //I should use only the pt-uncorr part of the systematic uncertainty (I am including some pt-correlated part as well...)
	fHistSpectrumStatHard[m]->SetName("fHistSpectrumHard_"+Smolt[m]);
	fHistSpectrumStatSoft[m] = (TH1F*)YieldMean_ReturnExtremeSoftHisto(fHistSpectrumSist[m]);
	fHistSpectrumStatSoft[m]->SetName("fHistSpectrumSoft_"+Smolt[m]);
      }

      //fit +1sigma sistematica
      fit_MTscalingUp[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Up");    
      fit_MTscalingUp[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrUp[m][typefit]=       fHistSpectrumStatUp[m]->Fit(    fit_MTscalingUp[m][typefit],"SR0IQ");
      fit_MTscalingUp[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitUp->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraFitUp->cd(3);
      else     canvasPtSpectraFitUp->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatUp[m]->Draw("same");
      fit_MTscalingUp[m][typefit]->Draw("same");
      legendfit->Draw("");

      //fit -1sigma sistematica
      fit_MTscalingDown[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Down");
      fit_MTscalingDown[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrDown[m][typefit]=       fHistSpectrumStatDown[m]->Fit(    fit_MTscalingDown[m][typefit],"SR0IQ");
      fit_MTscalingDown[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitDown->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraFitDown->cd(3);
      else     canvasPtSpectraFitDown->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatDown[m]->Draw("same");
      fit_MTscalingDown[m][typefit]->Draw("same");
      legendfit->Draw("");

      //fit of softened spectrum
      fit_MTscalingSoft[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Soft");
      fit_MTscalingSoft[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrSoft[m][typefit]=       fHistSpectrumStatSoft[m]->Fit(    fit_MTscalingSoft[m][typefit],"SR0IQ");
      fit_MTscalingSoft[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitSoft->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraFitSoft->cd(3);
      else     canvasPtSpectraFitSoft->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatSoft[m]->Draw("same");
      fit_MTscalingSoft[m][typefit]->Draw("same");
      legendfit->Draw("");

      //fit of hardened spectrum
      fit_MTscalingHard[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Hard");
      fit_MTscalingHard[m][typefit]->SetRange(LowRange[m], UpRange[m]);
     fFitResultPtrHard[m][typefit]=       fHistSpectrumStatHard[m]->Fit(    fit_MTscalingHard[m][typefit],"SR0IQ");
      fit_MTscalingHard[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitHard->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummoltMax)  canvasPtSpectraFitHard->cd(3);
      else     canvasPtSpectraFitHard->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatHard[m]->Draw("same");
      fit_MTscalingHard[m][typefit]->Draw("same");
      legendfit->Draw("");

      canvasDummy->cd();

      //1A. CALCULATE YIELDS
      YieldExtrLowPt[m][typefit]= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPt[m][typefit] = fit_MTscalingBis[m][typefit]->Integral(UpRangeSpectrumPart[m],20);
      YieldExtrLowPtAvg[m]+= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvg[m]+= fit_MTscalingBis[m][typefit]->Integral(UpRangeSpectrumPart[m],20);     

      hhoutYieldMy[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutRelStatMy[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutRelSystHighMy[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutRelSystLowMy[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

      hhoutAvgPtMy[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutAvgPtRelStatMy[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutAvgPtRelSystHighMy[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutAvgPtRelSystLowMy[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

      hhoutYieldMy[m] ->SetBinContent(typefit+1, YieldExtrLowPt[m][typefit] + YieldExtrHighPt[m][typefit]);
      //      hhoutRelStatMy[m] ->SetBinContent(typefit+1, sqrt(pow(    fit_MTscaling[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2) + pow ( fit_MTscalingBis[m][typefit]->IntegralError(UpRangeSpectrumPart[m],20,fFitResultPtr1[m][typefit] ->GetParams(),(fFitResultPtr1[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ),2)) ); //this is only the stat uncertainty related to the extrapolated part!
      hhoutRelStatMy[m] ->SetBinContent(typefit+1, 0);
      hhoutRelSystHighMy[m] ->SetBinContent(typefit+1, 0);
      hhoutRelSystLowMy[m] ->SetBinContent(typefit+1, 0);

      hhoutAvgPtRelStatMy[m] ->SetBinContent(typefit+1, 0);
      hhoutAvgPtRelSystHighMy[m] ->SetBinContent(typefit+1, 0);
      hhoutAvgPtRelSystLowMy[m] ->SetBinContent(typefit+1, 0);

      hExtrFractionLowPt[typefit] ->SetBinContent(m+1, YieldExtrLowPt[m][typefit]);
      hExtrFractionHighPt[typefit] ->SetBinContent(m+1, YieldExtrHighPt[m][typefit]);

      //1.B CALCULATE STATISTICAL ERROR OF EXTRAPOLATED YIELD
      /* WRONG: I am doing 4 measurements of the same quantity, but 4 fits with 4 different functions. The statistical error on the extraplated part should not be sigma / sqrt(4) but ~sigma (I can take an average sigma of the 4 functions) */
      YieldErrStatLowPtAvg[m]+= fit_MTscaling[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() );
      YieldErrStatHighPtAvg[m]+= fit_MTscalingBis[m][typefit]->IntegralError(UpRangeSpectrumPart[m],20,fFitResultPtr1[m][typefit] ->GetParams(),(fFitResultPtr1[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() );

      //1.C CALCULATE SYST UNCERTAIANTY OF EXTRAPOLATED YIELD
      YieldExtrLowPtAvgUp[m]+= fit_MTscalingUp[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvgUp[m]+= fit_MTscalingUp[m][typefit]->Integral(UpRangeSpectrumPart[m],20);     
      YieldExtrLowPtAvgDown[m]+= fit_MTscalingDown[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvgDown[m]+= fit_MTscalingDown[m][typefit]->Integral(UpRangeSpectrumPart[m],20);

      //1.D SYSTEMATIC RELATED TO CHOICE OF FIT FUNCTION
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

      //2.A CALCULATE AVERAGE PT VS MULT FROM FIT
      Int_t numInt = 600; //number of intervals in [0, UpRangePtInterval]
      Float_t Pti=0;
      Float_t UpRangePtInterval = 20;
      Float_t DeltaPt=	UpRangePtInterval/numInt; //interval width

      for (Int_t i=0; i<= numInt; i++){
	Pti = UpRangePtInterval/numInt*i;
	AvgPtFit[m][typefit] += fit_MTscalingBis[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingBis[m][typefit]->Integral(0,UpRangePtInterval);
	//2.C CALCULATE SYST UNCERTAIANTY OF AVG PT FROM FIT
	AvgPtFitHard[m][typefit] += fit_MTscalingHard[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingHard[m][typefit]->Integral(0,UpRangePtInterval);
        AvgPtFitSoft[m][typefit] += fit_MTscalingSoft[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingSoft[m][typefit]->Integral(0,UpRangePtInterval);
      } 

      AvgPt[m] += AvgPtFit[m][typefit];
      AvgPtHard[m] += AvgPtFitHard[m][typefit];
      AvgPtSoft[m] += AvgPtFitSoft[m][typefit];
      hhoutAvgPtMy[m] ->SetBinContent(typefit+1,  AvgPtFit[m][typefit]);

      //2.B CALCULATE STATISTICAL ERROR OF AVG PT FROM FIT
      fit_MTscalingTemp[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Temp");
      Int_t IterNum =200; //500
      fHistAvgPtDistr[m][typefit] = new TH1F (Form("fHistAvgPtDistr_m%i_typefit%i", m, typefit), Form("fHistAvgPtDistr_m%i_typefit%i", m, typefit), 10000, LimInfPtvsMult, LimSupPtvsMult);

      for (Int_t i= 0; i<IterNum;i++ ){
	Float_t Temp=0;
	for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
	  //extract a random number
	  gRandom->SetSeed(i*fHistSpectrumStat[m]->GetNbinsX()+b);
	  Temp=	  gRandom->Gaus(0, fHistSpectrumStat[m]->GetBinError(b));
	  fHistSpectrumTemp[m]->SetBinContent(b,fHistSpectrumStat[m]->GetBinContent(b) + Temp);
	  fHistSpectrumTemp[m]->SetBinError(b,fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b)*fHistSpectrumTemp[m]->GetBinContent(b));
	  //	  cout <<" stat " <<  fHistSpectrumStat[m]->GetBinContent(b) << " temp " <<  fHistSpectrumTemp[m]->GetBinContent(b) << endl;
	}
	fit_MTscalingTemp[m][typefit]->SetRange(LowRange[m], UpRange[m]);
	fFitResultPtrTemp[m][typefit]=fHistSpectrumTemp[m]->Fit(fit_MTscalingTemp[m][typefit],"SR0QI");
	fit_MTscalingTemp[m][typefit]->SetRange(0,20);

	if (typefit==0){ //--------------------------------
	  //to compute avg pt from spectrum and stat error of avg pt obtained from spectrum!
	  for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
	    if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowPtLimitForAvgPtFS[m]) continue;
	    AvgPtFSNum[m] +=  fHistSpectrumStat[m]->GetBinContent(b) *  fHistSpectrumStat[m]->GetBinCenter(b)* fHistSpectrumStat[m]->GetBinWidth(b);
	    AvgPtFSDenom[m] +=  fHistSpectrumStat[m]->GetBinContent(b) *  fHistSpectrumStat[m]->GetBinWidth(b);
	    AvgPtFSNumTemp[m] +=  fHistSpectrumTemp[m]->GetBinContent(b) *  fHistSpectrumTemp[m]->GetBinCenter(b)* fHistSpectrumTemp[m]->GetBinWidth(b);
	    AvgPtFSDenomTemp[m] +=  fHistSpectrumTemp[m]->GetBinContent(b) *  fHistSpectrumTemp[m]->GetBinWidth(b);
	  }
	  AvgPtFS[m] = AvgPtFSNum[m]/AvgPtFSDenom[m];
	  AvgPtFSTemp[m] = AvgPtFSNumTemp[m]/AvgPtFSDenomTemp[m];
	}//------------------------------------------------

	AvgPtFitTemp[m][typefit]=0;
	for (Int_t l=0; l<= numInt; l++){
	  Pti = UpRangePtInterval/numInt*l;
	  AvgPtFitTemp[m][typefit] += fit_MTscalingTemp[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingTemp[m][typefit]->Integral(0,UpRangePtInterval);
	  //	cout << " AVERAGE PT at iteration " << i << " of "<< numInt << " (pt = "  << Pti << "): " << 	AvgPtFit[m][typefit] << endl;
	}

	fHistAvgPtDistr[m][typefit] ->Fill(AvgPtFitTemp[m][typefit]);
	if (typefit ==0) fHistAvgPtFSDistr[m] ->Fill(AvgPtFSTemp[m]);
      } //end number of iterations

      TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
      if (type==0 && TypeAnalysis==0 && isppHM && m==0) nrebin[m] *=2;
      if (type==8 && TypeAnalysis==0 && isppHM) nrebin[m] = nrebin[m]*4;
      fHistAvgPtDistr[m][typefit]->Rebin(nrebin[m]);
      gaus->SetParameter(2, fHistAvgPtDistr[m][typefit]->GetRMS());
      fHistAvgPtDistr[m][typefit]->Fit(gaus, "q");
      AvgPtStatErrFit[m][typefit] = 	fHistAvgPtDistr[m][typefit]->GetRMS()/	fHistAvgPtDistr[m][typefit]->GetMean()* AvgPtFit[m][typefit];
      //AvgPtStatErrFit[m][typefit] = 	gaus->GetParameter(2)/gaus->GetParameter(1)* AvgPtFit[m][typefit];
      AvgPtStatErr[m]+= AvgPtStatErrFit[m][typefit];
    } //end loop typefit
    
    //2.C SYSTEMATIC RELATED TO CHOICE OF FIT FUNCTION
    AvgPtFitMaxFit[m] = 0;
    AvgPtFitMinFit[m] = 1000;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
      numfittipoEff++;
      if(AvgPtFitMaxFit[m] < AvgPtFit[m][typefit]){
	AvgPtFitMaxFit[m] = AvgPtFit[m][typefit];
      }
      if(AvgPtFitMinFit[m] > AvgPtFit[m][typefit]){
	AvgPtFitMinFit[m] = AvgPtFit[m][typefit];
      }
    }

    YieldExtrHighPtAvg[m]=YieldExtrHighPtAvg[m]/numfittipoEff;
    YieldExtrLowPtAvg[m]=YieldExtrLowPtAvg[m]/numfittipoEff;
    YieldExtrHighPtAvgUp[m]=YieldExtrHighPtAvgUp[m]/numfittipoEff;
    YieldExtrHighPtAvgDown[m]=YieldExtrHighPtAvgDown[m]/numfittipoEff;
    YieldExtrLowPtAvgUp[m]=YieldExtrLowPtAvgUp[m]/numfittipoEff;
    YieldExtrLowPtAvgDown[m]=YieldExtrLowPtAvgDown[m]/numfittipoEff;
    YieldErrStatHighPtAvg[m]=YieldErrStatHighPtAvg[m]/numfittipoEff;
    YieldErrStatLowPtAvg[m]=YieldErrStatLowPtAvg[m]/numfittipoEff;
    /*WRONG
    YieldErrStatHighPtAvg[m]=sqrt(YieldErrStatHighPtAvg[m])/numfittipoEff;
    YieldErrStatLowPtAvg[m]=sqrt(YieldErrStatLowPtAvg[m])/numfittipoEff;
    */

    //1.A EXTRAPOLATED YIELD
    YieldExtr[m] =     YieldExtrHighPtAvg[m]+    YieldExtrLowPtAvg[m]; 

    //1.B STATISTICAL ERROR OF EXTRAPOLATED YIELD
    YieldExtrErrStat[m] = sqrt(    pow(YieldErrStatHighPtAvg[m],2)+    pow(YieldErrStatLowPtAvg[m],2)); 

    //1.C SYSTEMATIC ERROR OF EXTRAPOLATED YIELD
    YieldErrSystExtrLowPt[m] = (YieldExtrLowPtAvgUp[m]- YieldExtrLowPtAvgDown[m])/2;
    YieldErrSystExtrHighPt[m] = (YieldExtrHighPtAvgUp[m]- YieldExtrHighPtAvgDown[m])/2;
    if (!isErrorAssumedPtCorr){
      YieldExtrErrSist[m]=sqrt(pow(  YieldErrSystExtrLowPt[m],2) + pow(  YieldErrSystExtrHighPt[m] ,2)) ;
    }
    else {
      YieldExtrErrSist[m]= YieldErrSystExtrLowPt[m]  + YieldErrSystExtrHighPt[m] ;
    }

    //1.D SYSTEMATIC ERROR OF YIELD RELATED TO CHOICE OF FIT FUNCTION
    Yield4FitErrSistHighPt[m] = (YieldExtrMaxHighPt[m]-YieldExtrMinHighPt[m])/2; 
    Yield4FitErrSistLowPt[m] = (YieldExtrMaxLowPt[m]-YieldExtrMinLowPt[m])/2; 
    YieldExtrErrSistFourFit[m]=sqrt(pow(    Yield4FitErrSistHighPt[m],2) + pow(    Yield4FitErrSistLowPt[m],2));

    //STATISTICAL AND SYSTEMATIC ERROR OF MEASURED YIELD
    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      if (NPtV0[v] >= UpRangeSpectrumPart[m]) continue; 
      YieldSpectrum[m] +=  ( fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1));
      YieldSpectrumErrStat[m] += pow ( fHistSpectrumStat[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1),2);
      if (!isErrorAssumedPtCorr){
	YieldSpectrumErrSist[m] += pow ( fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1),2);
	YieldSpectrumErrSistMultUnCorr[m] += pow ( fHistSpectrumSistMultUnCorr[m]->GetBinError(v+1)*fHistSpectrumSistMultUnCorr[m]->GetBinWidth(v+1),2);
      }
      else {
	YieldSpectrumErrSist[m] += fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1);
	YieldSpectrumErrSistMultUnCorr[m] += fHistSpectrumSistMultUnCorr[m]->GetBinError(v+1)*fHistSpectrumSistMultUnCorr[m]->GetBinWidth(v+1);
      }
    }

    YieldSpectrumErrStat[m]=sqrt( YieldSpectrumErrStat[m]);
    if (!isErrorAssumedPtCorr){
      YieldSpectrumErrSist[m]=sqrt( YieldSpectrumErrSist[m]);
      YieldSpectrumErrSistMultUnCorr[m]=sqrt( YieldSpectrumErrSistMultUnCorr[m]);
    }

    //FINAL RESULTS FOR YIELD
    Yield[m] =  YieldSpectrum[m]+YieldExtr[m];
    YieldErrStat[m] = sqrt(pow(YieldSpectrumErrStat[m],2) + pow(YieldExtrErrStat[m],2));
    FracErrorUnCorr[m] = YieldSpectrumErrSistMultUnCorr[m]/ YieldSpectrumErrSist[m];

    if (!isErrorAssumedPtCorr){
      YieldErrSist[m] = sqrt(pow(YieldSpectrumErrSist[m],2) +  pow(YieldErrSystExtrLowPt[m],2) + pow( YieldErrSystExtrHighPt[m],2));
      YieldErrSistMultUnCorr[m] = sqrt(pow(YieldSpectrumErrSistMultUnCorr[m],2) +  pow(FracErrorUnCorr[m]*YieldErrSystExtrLowPt[m],2) + pow(FracErrorUnCorr[m]* YieldErrSystExtrHighPt[m],2));
    }
    else {
      YieldErrSist[m] = YieldSpectrumErrSist[m] +  YieldErrSystExtrLowPt[m] + YieldErrSystExtrHighPt[m];
      YieldErrSistMultUnCorr[m] = YieldSpectrumErrSistMultUnCorr[m] +  FracErrorUnCorr[m]*YieldErrSystExtrLowPt[m] + FracErrorUnCorr[m]*YieldErrSystExtrHighPt[m];
    }
    YieldErrSistMy[m] = YieldErrSist[m];
    YieldErrSist[m]= sqrt(pow(YieldErrSist[m],2) + pow(Yield4FitErrSistLowPt[m],2) + pow(Yield4FitErrSistHighPt[m],2) );

    canvasExtrFractionLowPt->cd();
    gPad->SetLeftMargin(0.15);
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
      hExtrFractionLowPt[typefit]->SetBinContent(m+1, hExtrFractionLowPt[typefit]->GetBinContent(m+1)/Yield[m]);
      if (TypeAnalysis==0) {
	StyleHisto(hExtrFractionLowPt[typefit], 0, 0.5, ColorFit[typefit], 33, "p_{T} (GeV/c)", "", "ExtrapolatedFraction", 0, 0, 0);
      }
      else {
	if (type==0)	StyleHisto(hExtrFractionLowPt[typefit], 0, 0.1, ColorFit[typefit], 33, "p_{T} (GeV/c)", "", "ExtrapolatedFraction", 0, 0, 0);
	else 	StyleHisto(hExtrFractionLowPt[typefit], 0, 1, ColorFit[typefit], 33, "p_{T} (GeV/c)", "", "ExtrapolatedFraction", 0, 0, 0);
      }
      hExtrFractionLowPt[typefit]->Draw("same");
    }
    hExtrFractionLowPt[numfittipo]->SetBinContent(m+1, YieldExtrLowPtAvg[m]/Yield[m]);
    hExtrFractionLowPt[numfittipo]->SetLineColor(kBlack);
    hExtrFractionLowPt[numfittipo]->SetLineWidth(3);
    hExtrFractionLowPt[numfittipo]->Draw("same");

    canvasExtrFractionHighPt->cd();
    gPad->SetLeftMargin(0.15);
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
      hExtrFractionHighPt[typefit]->SetBinContent(m+1, hExtrFractionHighPt[typefit]->GetBinContent(m+1)/Yield[m]);
      if (TypeAnalysis==0)      StyleHisto(hExtrFractionHighPt[typefit], 0, 0.1, ColorFit[typefit], 33, "p_{T} (GeV/c)", "", "ExtrapolatedFraction", 0, 0, 0);
      else StyleHisto(hExtrFractionHighPt[typefit], 0, 0.005, ColorFit[typefit], 33, "p_{T} (GeV/c)", "", "ExtrapolatedFraction", 0, 0, 0);
      hExtrFractionHighPt[typefit]->Draw("same");
    }

    hExtrFractionHighPt[numfittipo]->SetBinContent(m+1, YieldExtrHighPtAvg[m]/Yield[m]);
    hExtrFractionHighPt[numfittipo]->SetLineColor(kBlack);
    hExtrFractionHighPt[numfittipo]->SetLineWidth(3);
    hExtrFractionHighPt[numfittipo]->Draw("same");

    //******************************************************************
    //FINAL VALUES TO AVERAGE PT VS MULT OBTAINED FROM FIT
    AvgPt[m] = AvgPt[m]/numfittipoEff; //AVERAGE PT
    AvgPtStatErr[m] = AvgPtStatErr[m]/numfittipoEff; //STATISTICAL ERROR OF AVG PT
    AvgPtSistErrFit[m] =  (AvgPtFitMaxFit[m] - AvgPtFitMinFit[m])/2; //SYST ERROR OF AVG PT RELATED TO CHOICE OF FIT FUNCTION
    if (isMeanMacro)     AvgPtSistErrFitMeanMacro[m] =  (AvgPtMeanMacroMaxFit[m] - AvgPtMeanMacroMinFit[m])/2; //SYST ERROR OF AVG PT CALCULATED BY YIELDMEAN MACRO RELATED TO CHOICE OF FIT FUNCTION

    //SYST ERROR OF AVG PT:
    AvgPtSoft[m] = AvgPtSoft[m]/numfittipoEff; 
    AvgPtHard[m] = AvgPtHard[m]/numfittipoEff;
    AvgPtSistErr[m] = TMath::Abs((AvgPtHard[m]- AvgPtSoft[m]))/2;
    //******************************************************************

    //******************************************************************
    //CALCULATE AVERAGE PT VS MULT FROM SPECTRUM
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowPtLimitForAvgPtFS[m]) continue;
      AvgPtFSNumHard[m] +=  fHistSpectrumStatHard[m]->GetBinContent(b) *  fHistSpectrumStatHard[m]->GetBinCenter(b)* fHistSpectrumStatHard[m]->GetBinWidth(b);
      AvgPtFSDenomHard[m] +=  fHistSpectrumStatHard[m]->GetBinContent(b) *  fHistSpectrumStatHard[m]->GetBinWidth(b);
      AvgPtFSNumSoft[m] +=  fHistSpectrumStatSoft[m]->GetBinContent(b) *  fHistSpectrumStatSoft[m]->GetBinCenter(b)* fHistSpectrumStatSoft[m]->GetBinWidth(b);
      AvgPtFSDenomSoft[m] +=  fHistSpectrumStatSoft[m]->GetBinContent(b) *  fHistSpectrumStatSoft[m]->GetBinWidth(b);
    }

    Float_t DeltaP =0;
    Float_t AvP =0;
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowPtLimitForAvgPtFS[m]) continue;
      DeltaP = fHistSpectrumStat[m]->GetBinWidth(b);
      AvP = fHistSpectrumStat[m]->GetBinCenter(b);
      //...not sure how I computed the errors below....:)
      /*
      AvgPtStatErrFS[m] += pow((DeltaP*AvP*AvgPtFSDenom[m]-pow(DeltaP,2)*AvP*fHistSpectrumStat[m]->GetBinContent(b))/pow(AvgPtFSDenom[m],2) *fHistSpectrumStat[m]->GetBinError(b) ,2);
      AvgPtSistErrFS[m] += pow((DeltaP*AvP*AvgPtFSDenom[m]-pow(DeltaP,2)*AvP*fHistSpectrumSist[m]->GetBinContent(b))/pow(AvgPtFSDenom[m],2) *fHistSpectrumSist[m]->GetBinError(b) ,2);
      */
    }

    AvgPtFSHard[m] = AvgPtFSNumHard[m]/AvgPtFSDenomHard[m];
    AvgPtFSSoft[m] = AvgPtFSNumSoft[m]/AvgPtFSDenomSoft[m];
    AvgPtSistErrFS[m] = TMath::Abs(AvgPtFSHard[m]-AvgPtFSSoft[m])/2;
    AvgPtStatErrFS[m] = fHistAvgPtFSDistr[m]->GetRMS()/	fHistAvgPtFSDistr[m]->GetMean()* AvgPtFS[m];
    //***************************************************

    //CALCULATE DISCRETIZATION ERROR
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowRange[m]) continue;
      //      AvgPtDiscrErr[m] += pow(fHistSpectrumStat[m]->GetBinContent(b),2) * pow(fHistSpectrumStat[m]->GetBinWidth(b),4)/12; 
      AvgPtDiscrErr[m] += fHistSpectrumStat[m]->GetBinContent(b) * pow(fHistSpectrumStat[m]->GetBinWidth(b),2); 
    }
    //    AvgPtDiscrErr[m] = sqrt(AvgPtDiscrErr[m])/YieldSpectrum[m];
    AvgPtDiscrErr[m] = AvgPtDiscrErr[m]/YieldSpectrum[m]/sqrt(12);

    //FILLING HISTOGRAMS
    hhoutYield[m] ->SetBinContent(numfittipo+1, hhoutYieldAvg[m]/numfittipoEff);
    hhoutYieldMy[m] ->SetBinContent(numfittipo+1, YieldExtr[m]);
    hhoutAvgPt[m] ->SetBinContent(numfittipo+1, hhoutAvgPtAvg[m]/numfittipoEff);
    hhoutAvgPtMy[m] ->SetBinContent(numfittipo+1, AvgPt[m]);

    hhoutYieldRatioToMine[m] = (TH1F*)    hhoutYield[m]->Clone(Form("hhoutYieldRatioToMine_m%i", m));
    hhoutYieldRatioToMine[m] ->Divide(hhoutYieldMy[m]);
    hhoutYieldRatioToMine[m] ->GetYaxis()->SetRangeUser(0.9,1.1);
    hhoutYieldRatioToMine[m]->SetLineColor(ColorMult[m]);

    hhoutAvgPtRatioToMine[m] = (TH1F*)    hhoutAvgPt[m]->Clone(Form("hhoutAvgPtRatioToMine_m%i", m));
    hhoutAvgPtRatioToMine[m] ->Divide(hhoutAvgPtMy[m]);
    hhoutAvgPtRatioToMine[m] ->GetYaxis()->SetRangeUser(0.9,1.1);
    hhoutAvgPtRatioToMine[m]->SetLineColor(ColorMult[m]);

    hhoutRelStat[m] ->SetBinContent(numfittipo+1, hhoutRelStatAvg[m]/numfittipoEff);
    hhoutRelSystHigh[m] ->SetBinContent(numfittipo+1, hhoutRelSystHighAvg[m]/numfittipoEff);
    hhoutRelSystLow[m] ->SetBinContent(numfittipo+1, hhoutRelSystLowAvg[m]/numfittipoEff);
    hhoutAvgPtRelStat[m] ->SetBinContent(numfittipo+1, hhoutAvgPtRelStatAvg[m]/numfittipoEff);
    hhoutAvgPtRelSystHigh[m] ->SetBinContent(numfittipo+1, hhoutAvgPtRelSystHighAvg[m]/numfittipoEff);
    hhoutAvgPtRelSystLow[m] ->SetBinContent(numfittipo+1, hhoutAvgPtRelSystLowAvg[m]/numfittipoEff);

    hhoutRelStat[m] ->Scale(1./Yield[m]);
    hhoutRelSystHigh[m] ->Scale(1./Yield[m]);
    hhoutRelSystLow[m] ->Scale(1./Yield[m]);
    hhoutAvgPtRelStat[m] ->Scale(1./AvgPt[m]);
    hhoutAvgPtRelSystHigh[m] ->Scale(1./AvgPt[m]);
    hhoutAvgPtRelSystLow[m] ->Scale(1./AvgPt[m]);

    hhoutRelStat[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutRelSystHigh[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutRelSystLow[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelStat[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutAvgPtRelSystHigh[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelSystLow[m]->GetYaxis()->SetRangeUser(0,0.1);

    hhoutRelStat[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystHigh[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystLow[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelStat[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystHigh[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystLow[m]->SetLineColor(ColorMult[m]);

    hhoutRelStatMy[m] ->SetBinContent(numfittipo+1, YieldErrStat[m]);
    hhoutRelStatMy[m] ->Scale(1./Yield[m]);
    hhoutRelSystHighMy[m] ->SetBinContent(numfittipo+1, YieldErrSistMy[m]/Yield[m]);
    hhoutRelSystLowMy[m] ->SetBinContent(numfittipo+1, YieldErrSistMy[m]/ Yield[m]);

    hhoutRelStatMy[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutRelSystHighMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutRelSystLowMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutRelStatMy[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystHighMy[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystLowMy[m]->SetLineColor(ColorMult[m]);

    hhoutAvgPtRelStatMy[m] ->SetBinContent(numfittipo+1, AvgPtStatErr[m]/AvgPt[m]);
    hhoutAvgPtRelSystHighMy[m] ->SetBinContent(numfittipo+1, AvgPtSistErr[m]/AvgPt[m]);
    hhoutAvgPtRelSystLowMy[m] ->SetBinContent(numfittipo+1, AvgPtSistErr[m]/ AvgPt[m]);

    hhoutAvgPtRelStatMy[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutAvgPtRelSystHighMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelSystLowMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelStatMy[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystHighMy[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystLowMy[m]->SetLineColor(ColorMult[m]);

    cout << "\n\n\e[35m*** Yields *** " << SmoltLegend[m] << "\e[39m " << endl;
    cout << "\nSpectra yield: " << endl;
    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      cout << NPtV0[v] << " < pt < " << NPtV0[v+1] << ": " << fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1) <<  " +- " << fHistSpectrumStat[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1) << " (stat.) +- " <<  fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1) << " (sist.) [mult. uncorr. fraction = "<< fHistSpectrumSistMultUnCorr[m]->GetBinError(v+1)/fHistSpectrumSistAll[m]->GetBinError(v+1)<< "]" << endl;
      //      cout << fHistSpectrumSistMultUnCorr[m]->GetBinError(v+1) << " " << fHistSpectrumSistAll[m]->GetBinError(v+1) << endl;
    }

    cout << "\npT integrated yield: " << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) +- "<< YieldErrSist[m] << " (sist) "<< endl;
    cout << " (stat rel: " <<  YieldErrStat[m]/ Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldErrSist[m]/ Yield[m] << ") [mult. uncorr. fraction = " << YieldErrSistMultUnCorr[m]/YieldErrSist[m]<< "]" << endl;
    cout << " (sist mult. uncorr. rel: " <<  YieldErrSistMultUnCorr[m]/ Yield[m] << ") " << endl;

    cout << "\nYield from spectrum: " <<   YieldSpectrum[m] << " +- " << YieldSpectrumErrStat[m]<<" (stat) +- " <<YieldSpectrumErrSist[m]<< " (sist)"  << endl;
    cout << " (stat rel: " <<  YieldSpectrumErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldSpectrumErrSist[m]/Yield[m] << ") [mult. uncorr. fraction = " << YieldSpectrumErrSistMultUnCorr[m]/YieldSpectrumErrSist[m] << "]" << endl;
    cout << " (sist mult uncorr rel: " <<  YieldSpectrumErrSistMultUnCorr[m]/Yield[m] << ") " << endl; 
    cout << "Errors relative to yield from spectrum:" << endl;
    cout << " (stat rel: " <<  YieldSpectrumErrStat[m]/YieldSpectrum[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldSpectrumErrSist[m]/YieldSpectrum[m] << ") [mult. uncorr. fraction = " << YieldSpectrumErrSistMultUnCorr[m]/YieldSpectrumErrSist[m] << "]" << endl;
    cout << " (sist mult uncorr rel: " <<  YieldSpectrumErrSistMultUnCorr[m]/YieldSpectrum[m] << ") " << endl; 

    cout << "\nYield from extrapolation: " <<     YieldExtr[m] << " +- " <<     YieldExtrErrStat[m]<< " (stat)  +- " << YieldExtrErrSistFourFit[m]<< " (sist 4 fit) " << YieldExtrErrSist[m] << " (sist extr) "<< endl;
    cout << " (stat rel: " <<  YieldExtrErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist 4 fit rel: " <<  YieldExtrErrSistFourFit[m]/Yield[m] << ") "<<endl;
    cout << " (sist rel: " <<  YieldExtrErrSist[m]/Yield[m] << ") "<<endl;
    cout << "Errors relative to extrapolated yield:" << endl;
    cout << " (stat rel: " <<  YieldExtrErrStat[m]/YieldExtr[m] << ") " << endl;
    cout << " (sist 4 fit rel: " <<  YieldExtrErrSistFourFit[m]/YieldExtr[m] << ") "<<endl;
    cout << " (sist rel: " <<  YieldExtrErrSist[m]/YieldExtr[m] << ") "<<endl;

    cout << "\nFraction of yield from low pt extrapolation " <<  YieldExtrLowPtAvg[m]/Yield[m] << endl;
    cout << "Fraction of yield from high pt extrapolation " <<  YieldExtrHighPtAvg[m]/Yield[m] << endl;
    cout << "\nFit range: " << LowRange[m] << " - " << UpRange[m] << endl;
    cout << "Bin content range: " << LowRangeSpectrumPart[m] << " - " << UpRangeSpectrumPart[m] << endl;

    cout << "\e[35m**** Avg pt obtained with different fit functions **** " << SmoltLegend[m] << "\e[39m" << endl;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
      cout << "\n ***" << endl;
      cout << "Avg pt (from fit) " << nameFit[typefit]<< ": " << AvgPtFit[m][typefit]<< " GeV/c" <<endl;
      cout << "Avg pt (from distribution of <pt> obtained with random variation of spectrum within statistical uncertainty): " << endl;
      //      cout << "Avg pt (from last random variation of spectrum) " << AvgPtFitTemp[m][typefit]<< " GeV/c " << endl;
      cout << " MEAN: " << fHistAvgPtDistr[m][typefit]->GetMean()<< ", RMS: " << fHistAvgPtDistr[m][typefit]->GetRMS() << endl;
      cout << " bin width of <pt> distribution (check it is much smaller than RMS): " << fHistAvgPtDistr[m][typefit]->GetXaxis()->GetBinWidth(1)<< endl;
    }
    cout << "\n\n\e[35m**** Avg pt from avg of different fit functions **** " << SmoltLegend[m] <<"\e[39m" << endl;
    cout << "Average pt: " << AvgPt[m]<< " +- " << AvgPtStatErr[m] << " (stat.) +- " << AvgPtSistErr[m] << " (syst. no function choice) +- "<<AvgPtSistErrFit[m] <<" (syst. function choice)"<<  endl;
    cout << "Rel errors:  " <<  AvgPtStatErr[m]/AvgPt[m] << " (stat.) +- " << AvgPtSistErr[m]/AvgPt[m] << " (syst. no function choice) +- "<< AvgPtSistErrFit[m]/AvgPt[m]<< " (syst. function choice)"<<  endl;
    
    cout << "\nAvg pt from spectrum (no fit, calculated for pt > " << 	LowPtLimitForAvgPtFS[m]  << "): " << AvgPtFS[m]<< " +- " << AvgPtStatErrFS[m] << " (stat.) +- " << AvgPtSistErrFS[m] << " (syst. no function choice) "<< endl; 
    cout << "From distribution of mean pt obtained by varying spectra within stat. uncertainty: \n MEAN: " << fHistAvgPtFSDistr[m]->GetMean()<< ", RMS: " << fHistAvgPtFSDistr[m]->GetRMS() << endl;  
    cout << "Rel errors:  " <<  AvgPtStatErrFS[m]/AvgPtFS[m] << " (stat.) +- " << AvgPtSistErrFS[m]/AvgPtFS[m] << " (syst. no function choice) +- "<< endl;
    cout << "Mean of fHistSpectrumStat: " << fHistSpectrumStat[m]->GetMean() << endl;

    if (isMeanMacro){
      AvgPtMeanMacro[m] =  hhoutAvgPt[m]->GetBinContent(numfittipo+1);
      AvgPtStatErrMeanMacro[m] = hhoutAvgPtRelStat[m]->GetBinContent(numfittipo+1);
      AvgPtSistErrMeanMacro[m] = (hhoutAvgPtRelSystHigh[m]->GetBinContent(numfittipo+1) + hhoutAvgPtRelSystLow[m]->GetBinContent(numfittipo+1))/2;
    cout << "\nAverage pt from Mean Macro: " << AvgPtMeanMacro[m] << " +- " <<  AvgPtStatErrMeanMacro[m] << " (stat.) +- " << AvgPtSistErrMeanMacro[m]<< " (syst. no function choice) +- "<< "-" <<" (syst. function choice)"<<  endl;
    cout << "Rel errors:  " <<  AvgPtStatErrMeanMacro[m]/AvgPtMeanMacro[m] << " (stat.) +- " << AvgPtSistErrMeanMacro[m]/AvgPtMeanMacro[m] << " (syst. no function choice) " <<  AvgPtSistErrFitMeanMacro[m]/AvgPtMeanMacro[m]<< " (syst. function choice)"<<  endl;
    }

    cout << "Discretization error: " << AvgPtDiscrErr[m] << " (Rel to AvgPt from mean macro: " << AvgPtDiscrErr[m]/AvgPtMeanMacro[m] << ") " << endl;

  }//end loop m


  cout << "\n\nYields : " << endl;
  for (Int_t m=0; m<nummoltMax +1; m++){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && m<2) continue;
    cout <<"\n\e[32m************ Multiplicity: " << SmoltLegend[m] << " *****************\e[39m"  << endl;

    cout << "pT integrated yield: " << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) +- "<< YieldErrSist[m] << " (sist) "<< endl;
    cout << " (stat rel: " <<  YieldErrStat[m]/ Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldErrSist[m]/ Yield[m] << ") " << endl;
    cout << " (sist mult. uncorr. rel: " <<  YieldErrSistMultUnCorr[m]/ Yield[m] << " = " << YieldErrSistMultUnCorr[m]/YieldErrSist[m] << " of the total) " << endl;

    cout << "\nYield from spectrum: " <<   YieldSpectrum[m] << " +- " << YieldSpectrumErrStat[m]<<" (stat) +- " <<YieldSpectrumErrSist[m]<< " (sist)"  << endl;
    cout << " (stat rel: " <<  YieldSpectrumErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldSpectrumErrSist[m]/Yield[m] << ") "<<endl;
    cout << " (sist mult uncorr rel: " <<  YieldSpectrumErrSistMultUnCorr[m]/Yield[m] << " = " << YieldSpectrumErrSistMultUnCorr[m]/YieldSpectrumErrSist[m] << " of the total) "<<endl;

    cout << "\nYield from extrapolation: " <<     YieldExtr[m] << " +- " <<     YieldExtrErrStat[m]<< " (stat)  +- " << YieldExtrErrSistFourFit[m]<< " (sist 4 fit) " << YieldExtrErrSist[m] << " (sist extr) "<< endl;
    cout << " (stat rel: " <<  YieldExtrErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist 4 fit rel: " <<  YieldExtrErrSistFourFit[m]/Yield[m] << ") "<<endl;
    cout << " (sist rel: " <<  YieldExtrErrSist[m]/Yield[m] << ") "<<endl;

    cout << "\nFraction of yield from low pt extrapolation " <<  YieldExtrLowPtAvg[m]/Yield[m] << endl;
    cout << "Fraction of yield from high pt extrapolation " <<  YieldExtrHighPtAvg[m]/Yield[m] << endl;
    cout << "\nFit range: " << LowRange[m] << " - " << UpRange[m] << endl;
    cout << "Bin content range: " << LowRangeSpectrumPart[m] << " - " << UpRangeSpectrumPart[m] << endl;
  }

  TLegend *legendYield=new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend *legendPtvsMult=new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend *legendYieldErr=new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legendPtErr=new TLegend(0.6, 0.6, 0.9, 0.9);

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
  TH1F* fHistYieldStat=new TH1F ("fHistYieldStatD","fHistYieldStatD",200,0,40);
  TH1F* fHistYieldSist=new TH1F ("fHistYieldSist","fHistYieldSist",200,0,40);

  TH1F* fHistPtvsMultStat = new TH1F ("fHistPtvsMultStat","fHistPtvsMultStat",200,0,40);
  TH1F* fHistPtvsMultSist = new TH1F ("fHistPtvsMultSist","fHistPtvsMultSist",200,0,40);
  TH1F* fHistPtvsMultStatFromSpectrum = new TH1F ("fHistPtvsMultStatFromSpectrum","fHistPtvsMultStatFromSpectrum",200,0,40);
  TH1F* fHistPtvsMultSistFromSpectrum = new TH1F ("fHistPtvsMultSistFromSpectrum","fHistPtvsMultSistFromSpectrum",200,0,40);
  TH1F* fHistPtvsMultStatMeanMacro = new TH1F ("fHistPtvsMultStatMeanMacro","fHistPtvsMultStatMeanMacro",200,0,40);
  TH1F* fHistPtvsMultSistMeanMacro = new TH1F ("fHistPtvsMultSistMeanMacro","fHistPtvsMultSistMeanMacro",200,0,40);

  TH1F* fHistPtStatRelErr = new TH1F ("fHistPtStatRelErr","fHistPtStatRelErr",200,0,40);
  TH1F* fHistPtSistRelErr = new TH1F ("fHistPtSistRelErr","fHistPtSistRelErr",200,0,40);
  TH1F* fHistPtStatRelErrFromSpectrum = new TH1F ("fHistPtStatRelErrFromSpectrum","fHistPtStatRelErrFromSpectrum",200,0,40);
  TH1F* fHistPtSistRelErrFromSpectrum = new TH1F ("fHistPtSistRelErrFromSpectrum","fHistPtSistRelErrFromSpectrum",200,0,40);
  TH1F* fHistPtStatRelErrMeanMacro = new TH1F ("fHistPtStatRelErrMeanMacro","fHistPtStatRelErrMeanMacro",200,0,40);
  TH1F* fHistPtSistRelErrMeanMacro = new TH1F ("fHistPtSistRelErrMeanMacro","fHistPtSistRelErrMeanMacro",200,0,40);
  TH1F* fHistPtSistRelErrFit = new TH1F ("fHistPtSistRelErrFit","fHistPtSistRelErrFit",200,0,40);

  TH1F* fHistYieldSistNoExtr=new TH1F ("fHistYieldSistNoExtr","fHistYieldSistNoExtr",200,0,40);
  TH1F* fHistYieldSistMultUnCorr=new TH1F ("fHistYieldSistMultUnCorr","fHistYieldSistMultUnCorr",200,0,40);
  TH1F* fHistYieldStatRelErr=new TH1F ("fHistYieldStatRelErr","fHistYieldStatRelErr",200,0,40);
  TH1F* fHistYieldSistRelErr=new TH1F ("fHistYieldSistRelErr","fHistYieldSistRelErr",200,0,40);
  TH1F* fHistYieldSistNoExtrRelErr=new TH1F ("fHistYieldSistNoExtrRelErr","fHistYieldSistNoExtrRelErr",200,0,40);
  TH1F* fHistYieldSist4FitRelErr=new TH1F ("fHistYieldSist4FitRelErr","fHistYieldSist4FitRelErr",200,0,40);
  TH1F* fHistYieldSistMultUnCorrRelErr=new TH1F ("fHistYieldSistMultUnCorrRelErr","fHistYieldSistMultUnCorrRelErr",200,0,40);
  TH1F* fHistYieldFractionSistMultUnCorr=new TH1F ("fHistYieldFractionSistMultUnCorr","fHistYieldFractionSistMultUnCorr",200,0,40);
  TH1F* fHistYieldStatEtaEff=new TH1F ("fHistYieldStatEtaEff","fHistYieldStatEtaEff",200,0,40);
  TH1F* fHistYieldStatEtaEffRatio;
  TH1F* fHistYieldStatEtaEffRatioRef;
  TH1F* fHistYieldStatNormCorr=new TH1F ("fHistYieldStatNormCorr","fHistYieldStatNormCorr",200,0,40);
  TH1F* fHistYieldStatNormCorrRatio;
  TH1F* fHistYieldStatNormCorrRatioRef;

  Float_t   dNdEta[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};
  if (isdNdEtaTriggered){ //values estimated with official task
    dNdEta[0] = 24.039;
    dNdEta[1] = 18.930;
    dNdEta[2] = 14.799;
    dNdEta[3] = 10.703;
    dNdEta[4] = 7.396;
    dNdEta[5] = 15.853;
  }
  if (isppHM) {
    dNdEta[0] = 39.40;
    dNdEta[1] = 36.89;
    dNdEta[2] = 35.16;
    dNdEta[3] = 32.57;
    dNdEta[4] = 30.43;
    dNdEta[5] = 31.5;
    if (MultBinning==1){
      dNdEta[0] =0;
      dNdEta[1] =0;
      dNdEta[2] =36.29; //values from 16l with 18d8 MC
      dNdEta[3] =32.57;
      dNdEta[4] =30.43;
      dNdEta[5] = 31.5;
      if (isdNdEtaTriggered){ //values estimated with official task
	dNdEta[0] =0;
	dNdEta[1] =0;
	dNdEta[2] =37.6; 
	dNdEta[3] =34.063;
	dNdEta[4] =32.047;
	dNdEta[5] =33.479;
      }
    }
  }
  if (ispp5TeV){
    dNdEta[0] = 15.27;
    dNdEta[1] = 11.91;
    dNdEta[2] = 8.73;
    dNdEta[3] = 5.78;
    dNdEta[4] = 3.03;
    dNdEta[5] = 5.49;
    if (MultBinning==3){
      /*      
      dNdEta[0] = 13.595; //estimated by me = 13.89, tabulated = 13.595
      dNdEta[1] = 4.91;//estimated by me = 6.95, tabulated= 4.91
      dNdEta[5] = 5.49;
      */
      dNdEta[0] = 13.89; //these values were estimated by me
      dNdEta[1] = 6.95;
      dNdEta[5] = 5.49;
      if (isdNdEtaTriggered){ //values estimated with official task
	dNdEta[0] = 16.93;
	dNdEta[1] = 10.22;
	dNdEta[5] = 12.51;
      }
    }
  }

  if (isGenOnTheFly){
    if (isdNdEtaTriggered){
      if (MonashTune==1) { //ropes
	dNdEta[0] = 8.45;
	dNdEta[1] = 12.42;
	dNdEta[2] = 15.11;
	dNdEta[3] = 17.15;
	dNdEta[4] = 19.12;
	dNdEta[5] = 20.99;
	dNdEta[6] = 22.75;
	dNdEta[7] = 24.84;
	dNdEta[8] = 27.36;
	dNdEta[9] = 30.35;
	dNdEta[10] = 16.5;
      } else if (MonashTune==2){ //default
	dNdEta[0] = 8.27;
	dNdEta[1] = 12.48;
	dNdEta[2] = 15.23;
	dNdEta[3] = 17.26;
	dNdEta[4] = 19.22;
	dNdEta[5] = 21.1;
	dNdEta[6] = 22.9;
	dNdEta[7] = 25.1;
	dNdEta[8] = 27.81;
	dNdEta[9] = 31.26;
	dNdEta[10] = 16.78;
      } else if (MonashTune==3){ //default
	dNdEta[0] = 7.97;
	dNdEta[1] = 11.63;
	dNdEta[2] = 14.26;
	dNdEta[3] = 16.29;
	dNdEta[4] = 18.29;
	dNdEta[5] = 20.25;
	dNdEta[6] = 22.16;
	dNdEta[7] = 24.55;
	dNdEta[8] = 27.49;
	dNdEta[9] = 31.61;
	dNdEta[10] = 17.29;
      }
    }
    else { 
      //    Monash Default CR, all events INEL > 0 (to simulate tabulated values)
      dNdEta[0] = 4.01;
      dNdEta[1] = 9.04;
      dNdEta[2] = 12.1;
      dNdEta[3] = 14.29;
      dNdEta[4] = 16.4;
      dNdEta[5] = 18.46;
      dNdEta[6] = 20.37;
      dNdEta[7] = 22.68;
      dNdEta[8] = 25.58;
      dNdEta[9] = 29.09;
      dNdEta[10] = 7.00;
    }
  }

  Int_t LimSupChi=120;
  if (type==1) LimSupChi=100;
  canvasFitResult->cd();
  for (Int_t typefit=0; typefit<numfittipo; typefit++){
    if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
    StyleHisto( hFitResult[typefit], 0,LimSupChi, ColorFit[typefit], 33, "Multiplicity class", "#chi^{2}/NDF", "#chi^{2} vs multiplicity", 0, 0, 0);
    hFitResult[typefit]->GetYaxis()->SetTitleOffset(1.2);
    hFitResult[typefit]->Draw("same");
    if (typefit==0)    legendfit->Draw("");
  }

  for(Int_t m=0; m<nummoltMax; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    fHistPtvsMultStat->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPt[m]);
    fHistPtvsMultSist->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPt[m]);
    fHistPtvsMultStat->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErr[m]);
    fHistPtvsMultSist->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),sqrt(pow(AvgPtSistErr[m],2) + pow(AvgPtSistErrFit[m],2)));
    fHistPtvsMultStatFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtFS[m]);
    fHistPtvsMultSistFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtFS[m]);
    fHistPtvsMultStatFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErrFS[m]);
    fHistPtvsMultSistFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtSistErrFS[m]);
    if (isMeanMacro){
    fHistPtvsMultStatMeanMacro->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtMeanMacro[m]);
    fHistPtvsMultSistMeanMacro->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtMeanMacro[m]);
    fHistPtvsMultStatMeanMacro->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]), AvgPtStatErrMeanMacro[m]);
    fHistPtvsMultSistMeanMacro->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]), sqrt(pow(AvgPtSistErrMeanMacro[m],2) + pow(AvgPtSistErrFitMeanMacro[m]/AvgPt[m]*AvgPtMeanMacro[m],2)));
    }

    fHistPtStatRelErr->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErr[m]/AvgPt[m]);
    fHistPtStatRelErr->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),0);
    fHistPtSistRelErr->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErr[m]/AvgPt[m]);
    fHistPtSistRelErr->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);

    fHistPtStatRelErrFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErrFS[m]/AvgPtFS[m]);
    fHistPtStatRelErrFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),0);
    fHistPtSistRelErrFromSpectrum->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErrFS[m]/AvgPtFS[m]);
    fHistPtSistRelErrFromSpectrum->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);

    if (isMeanMacro){
      fHistPtStatRelErrMeanMacro->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErrMeanMacro[m]/AvgPtMeanMacro[m]);
      fHistPtStatRelErrMeanMacro->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),0);
      fHistPtSistRelErrMeanMacro->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErrMeanMacro[m]/AvgPtMeanMacro[m]);
      fHistPtSistRelErrMeanMacro->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);
    }

    fHistPtSistRelErrFit->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErrFitMeanMacro[m]/AvgPtMeanMacro[m]);
    fHistPtSistRelErrFit->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);
  }

  for(Int_t m=0; m<nummoltMax; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    fHistYieldStat->SetBinContent(fHistYieldStat->FindBin(dNdEta[m]),Yield[m]);
    fHistYieldStat->SetBinError(fHistYieldStat->FindBin(dNdEta[m]),YieldErrStat[m]);
    fHistYieldSist->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),Yield[m]);
    fHistYieldSist->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),YieldErrSist[m]);
    fHistYieldSistNoExtr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),Yield[m]);
    fHistYieldSistNoExtr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]), YieldSpectrumErrSist[m]);
    fHistYieldSistMultUnCorr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),Yield[m]);
    fHistYieldSistMultUnCorr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]), YieldErrSistMultUnCorr[m]);

    fHistYieldStatRelErr->SetBinContent(fHistYieldStat->FindBin(dNdEta[m]),YieldErrStat[m]/Yield[m]);
    fHistYieldStatRelErr->SetBinError(fHistYieldStat->FindBin(dNdEta[m]),0);
    fHistYieldSistRelErr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),YieldErrSist[m]/Yield[m]);
    fHistYieldSistRelErr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),0);
    fHistYieldSistNoExtrRelErr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]), YieldSpectrumErrSist[m]/Yield[m]);
    fHistYieldSistNoExtrRelErr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),0);
    fHistYieldSist4FitRelErr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),  YieldExtrErrSistFourFit[m]/Yield[m]);
    fHistYieldSist4FitRelErr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),0);
    fHistYieldSistMultUnCorrRelErr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]), YieldErrSistMultUnCorr[m]/Yield[m]);
    fHistYieldSistMultUnCorrRelErr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),0);
    fHistYieldFractionSistMultUnCorr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]), YieldErrSistMultUnCorr[m]/YieldErrSist[m]);
    fHistYieldFractionSistMultUnCorr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),0);

    fHistYieldStat->SetMarkerSize(1.5);
    fHistYieldSist->SetMarkerSize(1.5);
    fHistYieldSistNoExtr->SetMarkerSize(1.5);
    fHistYieldSistMultUnCorr->SetMarkerSize(1.5);
    fHistYieldStatRelErr->SetMarkerSize(2);
    fHistYieldSistRelErr->SetMarkerSize(2);
    fHistYieldSistNoExtrRelErr->SetMarkerSize(2);
    fHistYieldSist4FitRelErr->SetMarkerSize(2);
    fHistYieldSistMultUnCorrRelErr->SetMarkerSize(2);
    fHistYieldFractionSistMultUnCorr->SetMarkerSize(2);

  }


  Float_t Xl[nummolt+1]= {0};
  Float_t Xh[nummolt+1]= {0};
  Float_t multctrbin[nummolt+1] ={0} ;

  Float_t LimSupdNdEtaAxis = 25;
  if (isppHM) LimSupdNdEtaAxis = 40;
  if (isGenOnTheFly) LimSupdNdEtaAxis = 40;

  for (Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    multctrbin[m] =   fHistYieldStat->GetXaxis()->GetBinCenter(  fHistYieldStat->FindBin(dNdEta[m]));
  }

  canvasYield->cd();
  StyleHisto(fHistYieldStat, LimInfYield, LimSupYield, Color[TypeAnalysis], 1, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity",1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSist, LimInfYield, LimSupYield, Color[TypeAnalysis], 1, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSistNoExtr, LimInfYield, LimSupYield, Color[TypeAnalysis], 1, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSistMultUnCorr, LimInfYield, LimSupYield, Color[TypeAnalysis], 1, titleYieldX, titleYieldY, titleYield[TypeAnalysis] + " yield vs multiplicity", 1, 0, LimSupdNdEtaAxis);
     
  fHistYieldStat->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldSist->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldSistNoExtr->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldSistMultUnCorr->GetYaxis()->SetTitleOffset(1.2);
  if (TypeAnalysis==0){
    fHistYieldStat->Fit(pol0, "R0");
  }

  fHistYieldStat->DrawClone("e");
  //  fHistYieldSistNoExtr->SetFillStyle(9);
  fHistYieldSistNoExtr->SetFillStyle(3001);
  fHistYieldSistNoExtr->SetFillColorAlpha(Color[TypeAnalysis], 1);
  fHistYieldSistNoExtr->Draw("same e2");
  fHistYieldSistMultUnCorr->SetFillStyle(1001);
  fHistYieldSistMultUnCorr->SetFillColorAlpha(Color[TypeAnalysis], 1);
  fHistYieldSistMultUnCorr->Draw("same e2");
  fHistYieldSist->SetFillStyle(0);
  //  if (!(TypeAnalysis==0 && type==8))  fHistYieldSist->Draw("same e2");
  fHistYieldSist->Draw("same e2");
  legendYield->AddEntry(fHistYieldSist, "syst.", "ef");
  legendYield->AddEntry(fHistYieldSistNoExtr, "syst. (no extr)", "ef");
  legendYield->AddEntry(fHistYieldSistMultUnCorr, "syst. (mult uncorr)", "ef");
  legendYield->AddEntry(fHistYieldStat, "stat.", "pel");
  legendYield->Draw("");
  canvasYield->SaveAs("provay.pdf");

  fHistYieldStat->SetName("fHistYieldStat");
  fHistYieldStat->SetTitle("fHistYieldStat");

  canvasPtvsMult->cd();
  StyleHisto(fHistPtvsMultStat, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 33, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity "+ RegionTypeNew[TypeAnalysis],1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultSist, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 33, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" + RegionTypeNew[TypeAnalysis], 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultStatFromSpectrum, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 4, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity "+ RegionTypeNew[TypeAnalysis],1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultSistFromSpectrum, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis], 4, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" + RegionTypeNew[TypeAnalysis], 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultStatMeanMacro, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis]+2, 22, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity "+ RegionTypeNew[TypeAnalysis],1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultSistMeanMacro, LimInfPtvsMult, LimSupPtvsMult, Color[TypeAnalysis]+2, 22, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" + RegionTypeNew[TypeAnalysis], 1, 0, LimSupdNdEtaAxis);
     
  fHistPtvsMultStat->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSist->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultStatFromSpectrum->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSistFromSpectrum->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultStatMeanMacro->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSistMeanMacro->GetYaxis()->SetTitleOffset(1.2);

  if (TypeAnalysis==0){
    fHistPtvsMultStat->Fit(pol0, "R0");
  }

  fHistPtvsMultStat->Draw("e");
  fHistPtvsMultSist->SetFillStyle(0);
  fHistPtvsMultSist->Draw("same e2");
  fHistPtvsMultStatFromSpectrum->Draw("same e");
  fHistPtvsMultSistFromSpectrum->SetFillStyle(0);
  fHistPtvsMultSistFromSpectrum->Draw("same e2");
  fHistPtvsMultStatMeanMacro->Draw("same e");
  fHistPtvsMultSistMeanMacro->SetFillStyle(0);
  fHistPtvsMultSistMeanMacro->Draw("same e2");
  legendPtvsMult->AddEntry(fHistPtvsMultSist, "syst. from fit", "ef");
  legendPtvsMult->AddEntry(fHistPtvsMultStat, "stat. from fit", "pel");
  legendPtvsMult->AddEntry(fHistPtvsMultSistFromSpectrum, Form("syst. from spectrum for pt>%.1f", LowPtLimitForAvgPtFSAllMult), "ef");
  legendPtvsMult->AddEntry(fHistPtvsMultStatFromSpectrum, Form("stat. from spectrum for pt>%.1f", LowPtLimitForAvgPtFSAllMult), "pel");
  legendPtvsMult->AddEntry(fHistPtvsMultSistMeanMacro, "syst. from MeanMacro", "ef");
  legendPtvsMult->AddEntry(fHistPtvsMultStatMeanMacro, "stat. from MeanMacro", "pel");
  legendPtvsMult->Draw("");

  canvasYieldErr->cd();
  StyleHisto(fHistYieldStatRelErr, 10e-5, LimSupYieldErr, Color[TypeAnalysis], 33, titleYieldX, titleYRel,titleYRel+" of "+ titleYield[TypeAnalysis] + " yield vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  fHistYieldStatRelErr->GetYaxis()->SetTitleOffset(1.2);
  StyleHisto(fHistYieldSistRelErr, 10e-5, LimSupYieldErr, Color[TypeAnalysis], 27, titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSistNoExtrRelErr, 10e-5, LimSupYieldErr, 881, 27,  titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSist4FitRelErr, 10e-5, LimSupYieldErr, 797, 20,  titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSistMultUnCorrRelErr, 10e-5, LimSupYieldErr, 921, 30,  titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldFractionSistMultUnCorr, 10e-5, 0.5, 921, 33,  titleYieldX, titleYRel, "Fraction of uncorrelated syst. uncertainty",  1, 0, LimSupdNdEtaAxis);
  fHistYieldSistMultUnCorrRelErr->SetMarkerSize(2.1);
  legendYieldErr->AddEntry(fHistYieldStatRelErr, "stat.", "pl");
  //  if (!(TypeAnalysis==0 && type==8))   legendYieldErr->AddEntry(fHistYieldSistRelErr, "syst.", "pl");
  legendYieldErr->AddEntry(fHistYieldSistRelErr, "syst.", "pl");
  legendYieldErr->AddEntry(fHistYieldSistNoExtrRelErr, "syst. (no extr)", "pl");
  legendYieldErr->AddEntry(fHistYieldSist4FitRelErr, "syst. 4 fit choice", "pl");
  legendYieldErr->AddEntry(fHistYieldSistMultUnCorrRelErr, "syst. (mult uncorr)", "pl");
  fHistYieldSistRelErr->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldStatRelErr->Draw("e p");
  fHistYieldSistNoExtrRelErr->Draw("same p");
  fHistYieldSist4FitRelErr->Draw("same p");
  fHistYieldSistMultUnCorrRelErr->Draw("same p");
  fHistYieldSistRelErr->Draw("same p");
  legendYieldErr->Draw("");

  canvasYieldSistFracMultUnCorr->cd();
  fHistYieldFractionSistMultUnCorr->Draw("same hist");

  canvasPtErr->cd();
  StyleHisto(fHistPtStatRelErr, 10e-5, LimSupPtErr, 628, 33, titleYieldX, titleYRel,titleYRel+" of "+ titleYield[TypeAnalysis] + " <pt> vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtStatRelErrFromSpectrum, 10e-5, LimSupPtErr, 630, 27, titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " <pt> vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtStatRelErrMeanMacro, 10e-5, LimSupPtErr, 881, 21,  titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);

  StyleHisto(fHistPtSistRelErr, 10e-5, LimSupPtErr, 418, 33, titleYieldX, titleYRel,titleYRel+" of "+ titleYield[TypeAnalysis] + " <pt> vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtSistRelErrFromSpectrum, 10e-5, LimSupPtErr, 420, 27, titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " <pt> vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtSistRelErrMeanMacro, 10e-5, LimSupPtErr, 424, 21,  titleYieldX, titleYRel, titleYRel+" of "+titleYield[TypeAnalysis] + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);

  StyleHisto(fHistPtSistRelErrFit, 10e-5, LimSupPtErr, 1, 22, titleYieldX, titleYRel,titleYRel+" of "+ titleYield[TypeAnalysis] + " <pt> vs multiplicity", 1, 0, LimSupdNdEtaAxis);

  legendPtErr->AddEntry(fHistPtStatRelErr, "stat. from fit", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErr, "syst. from fit", "pl");
  legendPtErr->AddEntry(fHistPtStatRelErrFromSpectrum, "stat. from spectrum", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErrFromSpectrum, "syst. from spectrum", "pl");
  legendPtErr->AddEntry(fHistPtStatRelErrMeanMacro, "stat. from mean macro", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErrMeanMacro, "syst. from mean macro", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErrFit, "syst. from fit function choice (Mean Macro)", "pl");

  fHistPtStatRelErr->Draw("e p");
  fHistPtSistRelErr->Draw("same p");
  fHistPtStatRelErrFromSpectrum->Draw("same e p");
  fHistPtSistRelErrFromSpectrum->Draw("same p");
  fHistPtStatRelErrMeanMacro->Draw("same e p");
  fHistPtSistRelErrMeanMacro->Draw("same p");
  fHistPtSistRelErrFit->Draw("same p");
  legendPtErr->Draw("");

  cout << "\n\n going to write on file " << endl;   
  if (isMeanMacro){
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      fileout->WriteTObject(hhoutYield[m]);
      fileout->WriteTObject(hhoutYieldMy[m]);
      fileout->WriteTObject(hhoutYieldRatioToMine[m]);
      fileout->WriteTObject(hhoutRelStat[m]);
      fileout->WriteTObject(hhoutRelStatMy[m]);
      fileout->WriteTObject(hhoutRelSystHigh[m]);
      fileout->WriteTObject(hhoutRelSystHighMy[m]);
      fileout->WriteTObject(hhoutRelSystLow[m]);
      fileout->WriteTObject(hhoutRelSystLowMy[m]);

      fileout->WriteTObject(hhoutAvgPt[m]);
      fileout->WriteTObject(hhoutAvgPtMy[m]);
      fileout->WriteTObject(hhoutAvgPtRatioToMine[m]);
      fileout->WriteTObject(hhoutAvgPtRelStat[m]);
      fileout->WriteTObject(hhoutAvgPtRelStatMy[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystHigh[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystHighMy[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystLow[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystLowMy[m]);

      for (Int_t typefit=0; typefit<numfittipo; typefit++){
	if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
	fileout->WriteTObject(hhout[typefit][m]);
	fileout->WriteTObject(fit_MTscaling[m][typefit]);
	//	fileout->WriteTObject(fHistAvgPtDistr[m][typefit]);
      }
    }
  }
  fileout->WriteTObject(canvasYield);
  fileout->WriteTObject(canvasYieldErr);
  fileout->WriteTObject(canvasYieldSistFracMultUnCorr);
  fileout->WriteTObject(canvasPtvsMult);
  fileout->WriteTObject(canvasPtErr);
  fileout->WriteTObject(fHistYieldStat);
  fileout->WriteTObject(fHistYieldSist);
  fileout->WriteTObject(fHistYieldSistNoExtr);
  fileout->WriteTObject(fHistYieldSistMultUnCorr);
  fileout->WriteTObject(fHistPtvsMultStat);
  fileout->WriteTObject(fHistPtvsMultSist);
  fileout->WriteTObject(fHistPtvsMultStatFromSpectrum);
  fileout->WriteTObject(fHistPtvsMultSistFromSpectrum);
  if (isMeanMacro){
    fileout->WriteTObject(fHistPtvsMultStatMeanMacro);
    fileout->WriteTObject(fHistPtvsMultSistMeanMacro);
  }
  fileout->WriteTObject(fHistYieldStatRelErr);
  fileout->WriteTObject(fHistYieldSistRelErr);
  fileout->WriteTObject(fHistYieldSistNoExtrRelErr);
  fileout->WriteTObject(fHistYieldSist4FitRelErr);
  fileout->WriteTObject(fHistYieldSistMultUnCorrRelErr);
  fileout->WriteTObject(fHistYieldFractionSistMultUnCorr);
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    fHistAvgPtFSDistr[m] ->GetXaxis()->SetRangeUser(LimInfPtvsMultTight, LimSupPtvsMultTight);
    fHistAvgPtDistr[m][0] ->GetXaxis()->SetRangeUser(LimInfPtvsMultTight, LimSupPtvsMultTight);
    fileout->WriteTObject( fHistAvgPtDistr[m][0]);
    fileout->WriteTObject( fHistAvgPtFSDistr[m]);
  }
  if (!(isMC && !isEfficiency)){
    if (TypeAnalysis==0 && type==8) {
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
    fileout->WriteTObject(canvasRatioFakeSBDef);
    fileout->WriteTObject(canvaspol0);
    fileout->WriteTObject(canvasBarlowpol0);
  }
  fileout->WriteTObject(canvasPtSpectraAll);
  fileout->WriteTObject(canvasBarlow);
  fileout->WriteTObject(canvasPtSpectra);
  fileout->WriteTObject(canvasPtSpectraFit);
  fileout->WriteTObject(canvasPtSpectraFitRatio);
  fileout->WriteTObject(canvasPtSpectraFitUp);
  fileout->WriteTObject(canvasPtSpectraFitDown);
  fileout->WriteTObject(canvasPtSpectraFitHard);
  fileout->WriteTObject(canvasPtSpectraFitSoft);
  fileout->WriteTObject(canvasPtSpectraFitBis);
  fileout->WriteTObject(canvasPtSpectraRelError);
  fileout->WriteTObject(canvasPtSpectraRelErrorAll);
  fileout->WriteTObject(canvasFitResult);
  fileout->WriteTObject(canvasExtrFractionLowPt);
  fileout->WriteTObject(canvasExtrFractionHighPt);
  if (isNormCorrFullyComputed==1) fileout->WriteTObject(canvasNormFactor);
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      for (Int_t typefit=0; typefit<numfittipo; typefit++){
	if (type==0 && TypeAnalysis==0 && MonashTune==3 && nameFit[typefit]=="Boltzmann" ) continue;
	fileout->WriteTObject(fit_MTscaling[m][typefit]);
      }
    fileout->WriteTObject(fHistSpectrumSistFakeSB[m]);
    fileout->WriteTObject(fHistSpectrumSistOOJSubDef[m]);
    fileout->WriteTObject(fHistSpectrumSistAll[m]);
    fileout->WriteTObject(fHistSpectrumSistMultUnCorr[m]);
    if (isNormCorrFullyComputed)    fileout->WriteTObject(fHistSpectrumStatNotNorm[m]);
    fileout->WriteTObject(fHistSpectrumStat[m]);
    fileout->WriteTObject(fHistSpectrumStatUp[m]);
    fileout->WriteTObject(fHistSpectrumStatDown[m]);
    fileout->WriteTObject(fHistSpectrumStatHard[m]);
    fileout->WriteTObject(fHistSpectrumStatSoft[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorMCChoice[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorFakeSB[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorOOJSubDef[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorAll[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorMultUnCorr[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorMB[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorMBCorrection[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorIBPU[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorOOBPU[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorSE[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorDCAz[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorPurity[m]);
    if (type==0)        fileout->WriteTObject(fHistSpectrumSistRelErrorMCChoice[m]);
    if (TypeAnalysis!=2) {
      fileout->WriteTObject(fHistSpectrumSistRelErrorDeltaEta[m]);
      fileout->WriteTObject(fHistSpectrumSistRelErrorDPhi[m]);
    }
  }

  //  TString DirPicture = "PictureForNote/";
  if (TypeAnalysis==3) RegionType[TypeAnalysis] = "BulkBlue";
  TString DirPicture = PathOutPictures;
  TString TypeAnalysisPic = "";
  if (isppHM) TypeAnalysisPic ="_isppHM";
  else if (ispp5TeV) TypeAnalysisPic ="_ispp5TeV";
  canvasYield->SaveAs(DirPicture+"YieldvsMult"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasYieldErr->SaveAs(DirPicture+"YieldvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasYieldSistFracMultUnCorr->SaveAs(DirPicture+"FracSistMultUnCorr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtErr->SaveAs(DirPicture+"PtvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtvsMult->SaveAs(DirPicture+"AvgPtvsMult"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectra->SaveAs(DirPicture+"PtSpectra"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraFit->SaveAs(DirPicture+"PtSpectraFit"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraFitRatio->SaveAs(DirPicture+"PtSpectraFitRatio"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraFitBis->SaveAs(DirPicture+"PtSpectraFitBis"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraRelError->SaveAs(DirPicture+"PtSpectraRelErr"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasPtSpectraRelErrorAll->SaveAs(DirPicture+"PtSpectraRelErrAll"+tipo[type]+RegionType[TypeAnalysis]+".pdf");
  canvasFitResult->SaveAs(DirPicture+"ChiSquare"+tipo[type]+RegionType[TypeAnalysis]+".pdf");

  TString SisMCTruth = "";
  if (isMC && !isEfficiency) SisMCTruth = "MCTruth/"; 
  if (MonashTune==1) SisMCTruth += "PythiaRopes/";
  else   if (MonashTune==2) SisMCTruth += "PythiaMonash/";
  else   if (MonashTune==3) SisMCTruth += "EPOSLHC/";
  canvasYield->SaveAs("PicturePaperProposal/"+SisMCTruth+"YieldvsMult"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasYieldErr->SaveAs("PicturePaperProposal/"+SisMCTruth+"YieldvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasYieldSistFracMultUnCorr->SaveAs("PicturePaperProposal/"+SisMCTruth+"FracSistMultUnCorr"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtErr->SaveAs("PicturePaperProposal/"+SisMCTruth+"PtvsMultErr"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtvsMult->SaveAs("PicturePaperProposal/"+SisMCTruth+"AvgPtvsMult"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtSpectra->SaveAs("PicturePaperProposal/"+SisMCTruth+"PtSpectra"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtSpectraFit->SaveAs("PicturePaperProposal/"+SisMCTruth+"PtSpectraFit"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtSpectraFitRatio->SaveAs("PicturePaperProposal/"+SisMCTruth+"PtSpectraRatioFit"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtSpectraFitBis->SaveAs("PicturePaperProposal/"+SisMCTruth+"PtSpectraFitBis"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtSpectraRelError->SaveAs("PicturePaperProposal/"+SisMCTruth+"PtSpectraRelErr"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasPtSpectraRelErrorAll->SaveAs("PicturePaperProposal/"+SisMCTruth+"PtSpectraRelErrAll"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasFitResult->SaveAs("PicturePaperProposal/"+SisMCTruth+"ChiSquare"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasExtrFractionLowPt->SaveAs("PicturePaperProposal/"+SisMCTruth+"ExtrFractionLowPt"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");
  canvasExtrFractionHighPt->SaveAs("PicturePaperProposal/"+SisMCTruth+"ExtrFractionHighPt"+tipo[type]+RegionType[TypeAnalysis]+TypeAnalysisPic+".pdf");

  fileout->Close();

  if (TypeAnalysis==0)  cout <<" q " <<  pol0->GetParameter(0) << " red chisq " << pol0->GetChisquare()/pol0->GetNDF() << endl;

  cout << "\n\e[35mStarting from the file(s):"<< endl;
  cout  << "->for default spectra: \e[39m";
  cout << PathInDef << endl;
  if (TypeAnalysis==0) cout << "Rel uncertainty associated to dphi choice is taken from the whole multiplicity interval " << endl;
  if (TypeAnalysis!=2) cout << "Rel uncertainty associated to dphi choice is smoothed " << endl;
  if (!(isMC && !isEfficiency)){
    if (type==0 && TypeAnalysis==0)  cout << "\n\e[35m->to get uncertainty related to OOJ subtraction:\e[39m " << PathInpol0<< endl
				       ;
    if (type==8) cout << "\n\e[35m->to get uncertainty related to OOJ subtraction:\e[39m " << PathInOOJSubDef << endl;
    if (type==8 && SkipAssoc) cout << "\n\e[35m->to get uncertainty related to pt < pt,trig:\e[39m " <<   PathInLeadTrackDef << endl;
    if (type==8) {
      cout << "To get uncertainty related to the use of purity instead of sidebands" << endl;
      if (ispp5TeV) cout << PathInSBRelErr  << endl;
      else cout <<  PathInFakeSBDef << endl;
    }
  cout << "\e[35m->to get uncertainty related to sidebands:\e[39m ";
  if (type ==0) cout << PathInFakeSBDef<< endl;
  if (isNormCorrFullyComputed) cout << "\n\e[35mFile from where normalisation factor is taken:\e[39m " << SfileNormCorrFC << endl;

  cout << "\nWith respect to preliminaries: " << endl;
  cout << "1) The final spectra are obtained using PYTHIA8 efficiency, and are not an average of the spectra obtaiend with EPOS and PYTHIA8, becasue I couldn't compute 2D efficiency using EPOS. The relative uncertainty associated to the choice of MC is however computed from 0-100% class, starting from the file:\n" <<PathInMCChoiceDef << endl;
  }

  cout << "\nI have created the file:\n" << stringout << "\n" <<endl;

  if (TypeAnalysis==3) cout << "Syst. errors taken from bulk! (to be done for bulkblue) " << endl;
  cout << "Remember to change iternum!!!!!!!!" << endl;

  cout <<"Number of functions with which fit is performed " << numfittipoEff << endl;
}

