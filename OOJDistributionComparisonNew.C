#include <TLine.h>
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
#include <TFile.h>
#include <TLegend.h>

void OOJDistributionComparisonNew(Int_t indexSysV0=0,  Int_t sysTrigger=0,  Int_t indexsysTrigger=0,  Int_t sys=0,Bool_t ishhCorr=0, Float_t PtTrigMin1=3, Float_t PtTrigMax1 =15, Float_t PtTrigMin2=0.15, Float_t PtTrigMax2 =2.5, Bool_t SkipAssoc=1,Int_t israp=0, Int_t sysV0=0,Int_t type=8, Int_t PtIntervalShown=1,   TString year0 = "2016",TString yearPtTrig3=/*"161718_HM_hXi"/*/"161718Full_AOD234_hXi"/*"17pq_hXi"/*"1617_hK0s"/*"Run2DataRed_MECorr_hXi"*/,TString year=/*"161718_HM_hXi"*/"161718Full_AOD234_hXi"/*"1617_AOD234_hK0s"*"17pq_pp5TeV_hXi_pttrig0.15"/*"AllMC_hXi"/*"2016kehjl_hK0s"*"2016k_hK0s"/*"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"*/, TString yearMCPtTrig3=""/*"1617MC_hK0s"*/,TString yearMC=""/*/"AllMC_hXi"/*"2018f1_extra_hK0s"/*"2016kl_hXi"/*2018f1_extra_hK0s_30runs_150MeV"*/,  TString Path1 =""/*"_Jet0.75"/*"_PtTrigMax2.5_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, TString Path2 =""/*"_PtTrigMax2.5"/*"_NewMultClassBis_Jet0.75"*/, TString Dir ="FinalOutput/", Bool_t isEnlargedDeltaEta=0, Int_t isMC=0, Int_t MultBinning=0, Int_t PtBinning=0,  Bool_t isSysDef=1, Bool_t isDefaultSel=0, Bool_t isEPOS=0, Bool_t isEtaEff=1, Bool_t isNewInputPath=1, Bool_t isppHM=0, Bool_t isOOJFromK0s=0, Bool_t isBkgParab=0, Bool_t isOOJFromAllMult=0, Bool_t isMEFrom13TeV=0){

  //------ INFO ABOUT MAIN HISTOGRAMS (main ways to subtract OOJ)--------
  //---> Jet distr is projection of |dEta| < 0.75
  //---> Jet+OJ distr is projection of |dEta| < 1.17
  //  DEFAULT METHOD: taken from AngularCorrelation_fisrtCasc.C output (Jet distr - BULK Def)
  //  DEFAULT METHOD SMOOTH: (Jet distr - BULK Def Rebinned and Smoothed)
  //  REBSMOOTH: (Jet distr - BULK New 0-100% scaled M1, Rebinned and Smoothed) M1=scaled to Jet distr in 1 < dphi < 2
  //  REBSMOOTHFIT: (Jet distr - Pol 0 fit to BULK New 0-100% scaled M1, Rebinned and Smoothed) M1=scaled to Jet distr in 1 < dphi < 2
  //  REBSMOOTHBIS: (Jet distr - BULK New 0-100% scaled M2, Rebinned and Smoothed) M2=scaled to default BULK distribution in all dphi interval
  //  INCLUSIVEREBSMOOTH: (Jet+OJ distr - BULK New 0-100% scaled M2, Rebinned and Smoothed) M2=scaled to default BULK distribution in all dphi interval
  //----------------------------------------------------------------------

  if (year == "2016k_hK0s") Path1 = "_PtTrigMax2.5";
  if (yearPtTrig3 == "17pq_hXi") MultBinning=3;
  Bool_t IsSpecial=0;
  if (yearPtTrig3 == "1617_hK0s") IsSpecial=1;

  //if isEPOS==1, the MC used to correct for efficiency is EPOS (but just for the default ooj distirbutions, not the new ones

  if (isSysDef && sys!=0) return;
  if (isDefaultSel && sysTrigger!=0) return;
  //isMC = 0 if only data are studied, 1 if only MC is studied, 2 if both are studied 
  Int_t LimInfMC=0;
  Int_t LimSupMC=0;
  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};
  if (isMC==1) {
    LimInfMC=1; 
    LimSupMC=1; 
  }
  else  if (isMC==2) {
    LimInfMC=0; 
    LimSupMC=1; 
  }

  if (PtIntervalShown>6){
    cout << " Thera are not so many Pt assoc intervals " << endl;
    return;
  }
  if (ishhCorr==1){ year = "2016k"; Path1="_New"; Path2="_New"; yearMC = "2018f1_extra";}

  if (!ishhCorr && type==0 && (((year!="2016k_onlyTriggerWithHighestPt"&& year!="2018f1_extra_onlyTriggerWithHighestPt") || yearMC!="2018f1_extra_onlyTriggerWithHighestPt"))) {
    //    cout << "for hV0 correlation you should use year = 2016k_onlyTriggerWithHighestPt together with the MC yearMC = 2018f1_extra_onlyTriggerWithHighestPt" << endl;
    // return; 
  }

  Float_t PtTrigMin=0;
  gStyle->SetOptStat(0);
  //  gStyle->SetOptFit(1111);

  //lista degli effetti  sistematici studiati in questa macro
  if (sys==3 || sys>6) return;
  if (sysV0>6)return;
  if (ishhCorr && sys<4 && sys!=0) return; //se faccio correlazione hh non posso studiare sistematico associato a scelta regione Sidebands e peak nella distribuzione di massa invariante
  if (sysV0>2 && ishhCorr) return;
  //  if (sysTrigger!=0) return;
  if (sysV0!=0 && sys!=0) return;

  Double_t LimSupY[2] = {0.0006, 0.05}; 
  Double_t LimInfY[2] = {-0.0001, -0.004}; 
  TString hhCorr[2]= {"", "_hhCorr"};
  Int_t sysang=sys;
  Int_t sysOOJ = 0;
  if (sys==5){ sysOOJ=5; sysang=0;}

  TLegend *legendPt = new TLegend(0.7, 0.7, 0.9, 0.9);
  legendPt->SetHeader("p_T^{Assoc} intervals");

  TLegend *legendmult = new TLegend(0.5, 0.1, 0.9, 0.4);
  legendmult->SetHeader("Multiplicity classes");

  Dir+="DATA"+year0;
  TString file[2];
  TString fileOOJ[2];

  file[0] = yearPtTrig3;
  //  if (sysang==0)  file[0] += Path1;
  if (PtBinning>0) file[0]+= Form("_PtBinning%i", PtBinning);

  fileOOJ[0] = year;
  if (PtBinning>0) fileOOJ[0]+= Form("_PtBinning%i", PtBinning);
  if (sysOOJ==0)  fileOOJ[0] += Path1;
 
  file[1] = yearMCPtTrig3 + "_MCEff" ;  
  if (sysang==0)  file[1] += Path1;
  if (PtBinning>0) file[1]+= Form("_PtBinning%i", PtBinning);
  file[1]+=Path2; 

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9;
  Int_t numPtV0Max =numPtV0;
  if (PtBinning==0)numPtV0Max =numPtV0-1;
  Int_t   numPtV0MaxOOJ=5;
  if (type==0) numPtV0MaxOOJ=6;
  if (isOOJFromK0s) numPtV0MaxOOJ=numPtV0Max;
  if (isOOJFromAllMult) numPtV0MaxOOJ=numPtV0Max;
  Int_t   numPtV0MaxUniversal=0;
  const Int_t numPtTrig=10;
  const Int_t numtipo=10;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
 
  Int_t PtV0Min=1;
  if (!ishhCorr && type==0) PtV0Min=0;
 
  Int_t ColorPt[numPtV0]= {401,801,628,909,881,860,868,842, 921};

  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  //  TString SPtV0[numPtV0]={"", "", "0.5-1", "1-1.5","1.5-2","2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8,100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", ""};
  SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }

  TString SPtV01[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};

  if (PtBinning==1){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }

  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString BkgType[2]={"", "_isBkgParab"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  TLegend * legendDataMC = new TLegend(0.6, 0.6, 0.9, 0.9);
  Int_t ColorWidth[2]={868, 418};
  TString MCOrNot[2]={" Data" , " Pythia8" };
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  TString Smolt5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  Double_t Nmolt5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};
  TString Smolt[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,2,7,15,30,100}; 

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Nmolt[m] = Nmolt0[m];
      Smolt[m] = Smolt0[m];
    }
    else if (MultBinning==1){
      Nmolt[m] = Nmolt1[m];
      Smolt[m] = Smolt1[m];
    }
    else if (MultBinning==2){
      Nmolt[m] = Nmolt2[m];
      Smolt[m] = Smolt2[m];
    }
    else if (MultBinning==3){
      Smolt[m] = Smolt5TeV[m];
      Nmolt[m] = Nmolt5TeV[m];
    }
  }
  if (isppHM){
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
    if (MultBinning==1){
      Nmolt[1] = 0;
      Nmolt[2] = 0;
      Smolt[0] = "0-0a";
      Smolt[1] = "0-0b";
      Smolt[2] = "0-0.01";
    }
  }

  TString Ssideband[numSB]={"", "_SB"};
  TString Ssidebandtitle[numSB]={"", " (sidebands)"};
  TString Szeta[numzeta]={""};

  TString nameSE[nummolt+1][numPtV0];
  TString nameME[nummolt+1][numPtV0];
  TString nameMEEtaProj[nummolt+1][numPtV0];
  TString nameMEPhiProj[nummolt+1][numPtV0];
  TString namePhiProjJet[nummolt+1][numPtV0];
  TString namePhiProjInclusive[nummolt+1][numPtV0];
  TString namePhiProjJetZYAM[nummolt+1][numPtV0];
  TString namePhiProjJetFromBulkFit[nummolt+1][numPtV0];
  TString namePhiProjJetNotBulkSub[nummolt+1][numPtV0];
  TString namePhiProjBulk[nummolt+1][numPtV0];
  TString namePhiProjBulkRatio[nummolt+1][numPtV0];
  TString namePhiProjJetBulk[nummolt+1][numPtV0];
  TString nameEtaProj[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[nummolt+1][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[nummolt+1][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJet[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetInclusive[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRebSmoothCorrMult[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetCorrMultBis[nummolt+1][numPtV0][numPtTrig];
  TF1*  pol0JetRebSmoothBis[nummolt+1][numPtV0];
  TF1*  pol0JetRebSmooth[nummolt+1][numPtV0];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetZYAM[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[nummolt+1][numPtV0][numPtTrig];

  TF1*	  pol0ZYAM[nummolt+1][numPtV0][numPtTrig];
  TF1*	  pol0BulkBis[nummolt+1][numPtV0][numPtTrig];
  TF1*	  pol0BulkSmooth[nummolt+1][numPtV0][numPtTrig];

  TH1F *hDeltaEtaDeltaPhi_PhiProjInclusive[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubRebin[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNewMethod[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNDRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNRebinDRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethod[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkReb[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkRebSmoothed[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk_mall[numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkDefault[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaled[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledToInclusiveAllMult[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultRebin[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultRebinSmoothed[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkScaledRebin[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkDenom[nummolt+1][numPtV0][numPtTrig];
  TF1  *BulkRatioPol0[nummolt+1][numPtV0];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProjPHalf[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProjNHalf[nummolt+1][numPtV0][numPtTrig];
  TString PathIn[numPtTrig];
  TString PathInDefaultEta[numPtTrig];

  TF1* Baseline[nummolt+1][numPtV0][numPtTrig];
  TF1* Gaussian[nummolt+1][numPtV0][numPtTrig];
  TF1* GaussianAS[nummolt+1][numPtV0][numPtTrig];
  TF1* GaussianEta[nummolt+1][numPtV0][numPtTrig];

  //  TLine *tlineEtaSx=new TLine(-1.04, -TMath::Pi()/2, -1.04,  3*TMath::Pi()/2);
  //  TLine *tlineEtaDx=new TLine(1.04, -TMath::Pi()/2, 1.04,  3*TMath::Pi()/2);
  TLine *tlineEtaSx=new TLine(-0.78, -TMath::Pi()/2, -0.78,  3*TMath::Pi()/2);
  TLine *tlineEtaDx=new TLine(0.78, -TMath::Pi()/2, 0.78,  3*TMath::Pi()/2);
  TLine *tlineEtaInclSx=new TLine(-1.14, -TMath::Pi()/2, -1.14,  3*TMath::Pi()/2);
  TLine *tlineEtaInclDx=new TLine(1.14, -TMath::Pi()/2, 1.14,  3*TMath::Pi()/2);
  //  TLine *tlinePhiSx=new TLine(-1.5, -1.06, 1.5,  -1.06); //was 1.32
  //  TLine *tlinePhiDx=new TLine(-1.5, +1.06, 1.5,  +1.06);
  TLine *tlinePhiSx=new TLine(-1.5, -0.85, 1.5,  -0.85); //was 1.32
  TLine *tlinePhiDx=new TLine(-1.5, 0.85, 1.5,  0.85);
  TLine *tlinePhiBase=new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2,0);
  TLine *tlineAtOne=new TLine(-1, 1, 1, 1);
  TLine *tlineAtOneAllDeltaPhi=new TLine(-TMath::Pi()/2, 1, 3*TMath::Pi()/2, 1);

  TCanvas *  canvasWidthGaussian[2];
  TCanvas *  canvasWidthGaussianEta[2];
  TCanvas *  canvasWidthGaussianAS[2];
  for (Int_t i =0; i<2; i++){  
    canvasWidthGaussian[i] = new TCanvas (Form("canvasWidthGaussian%i",i), Form("canvasWidthGaussian%i",i), 800, 500);
    canvasWidthGaussian[i]->Divide(3,2);
    canvasWidthGaussianEta[i] = new TCanvas (Form("canvasWidthGaussianEta%i",i), Form("canvasWidthGaussianEta%i",i), 800, 500);
    canvasWidthGaussianEta[i]->Divide(3,2);
    canvasWidthGaussianAS[i] = new TCanvas (Form("canvasWidthGaussianAS%i",i), Form("canvasWidthGaussianAS%i",i), 800, 500);
    canvasWidthGaussianAS[i]->Divide(3,2);
  }

  TCanvas *canvasPlot[nummolt+1][2];
  TCanvas *canvasPlotME[nummolt+1][2];
  TCanvas *canvasPlotMEEtaProj[2];
  TCanvas *canvasPlotMEPhiProj[2];
  TCanvas *canvasPlotProj[nummolt+1][2];
  TCanvas *canvasPlotProjBis[nummolt+1][2];
  TCanvas *canvasPlotProjAllPt[nummolt+1][2];
  TCanvas *canvasPlotProjSmoothed[nummolt+1][2];
  TCanvas *canvasPlotProjRatioJet[nummolt+1][2];
  TCanvas *canvasPlotProjEta[nummolt+1][3][2];

  TCanvas *canvasPlotOOJDistr[nummolt+1][2][2];
  TCanvas *canvasPlotOOJDistrMultComp[2];

  TH1F *HistoWidthGaussian[nummolt+1][2][numPtTrig];
  TH1F *HistoWidthGaussianAS[nummolt+1][2][numPtTrig];
  TH1F *HistoWidthGaussianEta[nummolt+1][2][numPtTrig];
  TH1F *RatioHistoWidthGaussian[nummolt+1][numPtTrig];
  TH1F *RatioHistoWidthGaussianAS[nummolt+1][numPtTrig];
  TH1F *RatioHistoWidthGaussianEta[nummolt+1][numPtTrig];
  TString PhiRegion[3]={"All", "DeltaPhi>Pi/2", "DeltaPhi < Pi/2"};


  TString nomefileoutput= "OOJComparison"+yearPtTrig3+"_"+year;
  if (PtBinning>0) nomefileoutput+=Form("_PtBinning%i",PtBinning);
  //  if (sys==0) nomefileoutput +=Path1;
  TString nomepdffile;
  if(type>=0){
    nomefileoutput +="_"+tipo[type];
    nomefileoutput +=Srap[israp];
    nomefileoutput +=SSkipAssoc[SkipAssoc];
  }

  if(isSysDef && isDefaultSel)      nomefileoutput += hhCorr[ishhCorr] +  Form("_SysV0Default");
  else   if(isSysDef && !isDefaultSel && sysTrigger==0)      nomefileoutput += hhCorr[ishhCorr] +  Form("_SysV0index%i", indexSysV0);
  else   if(isSysDef && !isDefaultSel && sysTrigger==1)      nomefileoutput += hhCorr[ishhCorr] +  Form("_SysTindex%i", indexsysTrigger);
  else  nomefileoutput += hhCorr[ishhCorr]+ Form("_sys%i", sys);

  nomefileoutput += Form("_PtTrigMin%.1f_PtTrigMin%.1f_Output", PtTrigMin1, PtTrigMin2);
  if (isEPOS)   nomepdffile+= "_EPOS";
  if (isEPOS)   nomefileoutput+= "_EPOS";
  if (isEtaEff) nomefileoutput+= "_IsEtaEff";
  if (isMEFrom13TeV) nomefileoutput += "_isMEFrom13TeV";
  if (MultBinning!=0) nomefileoutput+= Form("_MultBinning%i", MultBinning);
  if (isOOJFromK0s) nomefileoutput+= "_isOOJFromK0s";
  if (isOOJFromAllMult) nomefileoutput+= "_isOOJFromAllMult";
  nomepdffile = nomefileoutput;
  nomefileoutput+= ".root";
   
  TFile *fileout = new TFile(nomefileoutput, "RECREATE");

  Float_t PtTrigChosen=0;
  Float_t PtTrigMax=0;

  Int_t numC =4;
  if (PtBinning>0) numC=6;
  if (isOOJFromAllMult) numC = numPtV0MaxOOJ-1;
  cout << " creating canvases " << endl;

  for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
    canvasPlotOOJDistrMultComp[IntisMC]=new TCanvas(Form("canvasPlotOOJDistrMultComp_MC%i", IntisMC), "canvasPlotOOJDistrMultComp", 1300, 800);
    canvasPlotOOJDistrMultComp[IntisMC]->Divide(numC,2);
  }

  for(Int_t m=nummolt; m>=0; m--){
    //      if (m==0) continue;
    cout << "\n\n m " << m << endl;
    for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
      for (Int_t i=0; i<2; i++){
	if(i==0)	canvasPlotOOJDistr[m][IntisMC][i]=new TCanvas(Form("canvasPlotOOJDistr_m%i_MC%i",m, IntisMC), "canvasPlotOOJDistr"+Smolt[m], 1300, 800);
	else 	canvasPlotOOJDistr[m][IntisMC][i]=new TCanvas(Form("canvasPlotOOJDistrRatio_m%i_MC%i",m, IntisMC), "canvasPlotOOJDistrRatio"+Smolt[m], 1300, 800);
	canvasPlotOOJDistr[m][IntisMC][i]->Divide(numC,2);
      }
      canvasPlot[m][IntisMC]=new TCanvas(Form("canvasPlot_m%i_MC%i",m, IntisMC), "canvasPlot"+Smolt[m], 1300, 800);
      canvasPlot[m][IntisMC]->Divide(numC,2);
      canvasPlotME[m][IntisMC]=new TCanvas(Form("canvasPlotME_m%i_MC%i",m,  IntisMC), "canvasPlotME"+Smolt[m], 1300, 800);
      canvasPlotME[m][IntisMC]->Divide(numC,2);
      canvasPlotProj[m][IntisMC]=new TCanvas(Form("canvasPlotProj_m%i_MC%i",m,  IntisMC), "canvasPlotProj"+Smolt[m], 1300, 800);
      canvasPlotProj[m][IntisMC]->Divide(numC,4);
      canvasPlotProjBis[m][IntisMC]=new TCanvas(Form("canvasPlotProjBis_m%i_MC%i",m,  IntisMC), "canvasPlotProjBis"+Smolt[m], 1300, 800);
      canvasPlotProjBis[m][IntisMC]->Divide(numC,2);
      canvasPlotProjAllPt[m][IntisMC]=new TCanvas(Form("canvasPlotProjAllPt_m%i_MC%i",m,  IntisMC), "canvasPlotProj"+Smolt[m], 1300, 800);
      canvasPlotProjAllPt[m][IntisMC]->Divide(numC,2);
      
      canvasPlotProjRatioJet[m][IntisMC]=new TCanvas(Form("canvasPlotProjRatioJet_m%i_MC%i",m,  IntisMC), "canvasPlotProjRatioJet"+Smolt[m], 1300, 800);
      canvasPlotProjRatioJet[m][IntisMC]->Divide(numC,2);
      for (Int_t t=0; t<3; t++){
	canvasPlotProjEta[m][t][IntisMC]=new TCanvas(Form("canvasPlotProjEta_m%i_%i_MC%i",m,t, IntisMC), "canvasPlotProjEta"+Smolt[m]+"_"+PhiRegion[t], 1300, 800);
	canvasPlotProjEta[m][t][IntisMC]->Divide(numC,2);
      }
    }
  }

  TCanvas*   canvasSummedPt=new TCanvas ("canvasSummedPt", "canvasSummedPt", 800, 500);
  canvasSummedPt->Divide(nummolt+1, 2);
  for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
    canvasPlotMEEtaProj[IntisMC]=new TCanvas(Form("canvasPlotMEEtaProj_%i", IntisMC), Form("canvasPlotMEEtaProj_%i", IntisMC), 1300, 800);
    canvasPlotMEEtaProj[IntisMC]->Divide(8,2);
    canvasPlotMEPhiProj[IntisMC]=new TCanvas(Form("canvasPlotMEPhiProj_%i", IntisMC), Form("canvasPlotMEPhiProj_%i", IntisMC), 1300, 800);
    canvasPlotMEPhiProj[IntisMC]->Divide(8,2);
  }

  //I start the loop overt the two differetn files ******************
  Float_t IntegralFactorBulkNew[nummolt+1][numPtV0][numPtTrig]={0};
  Float_t IntegralFactorBulkDef[nummolt+1][numPtV0][numPtTrig]={0};
  Float_t IntegralFactorBulkNewAllMult[nummolt+1][numPtV0][numPtTrig]={0};
  TLegend * OOJlegend = new TLegend (0.6, 0.7, 0.9, 0.9);
  cout << "\n\n\n********+ looping over the two files " << endl;
  Float_t ScalingFactorJetNotBulkSub[nummolt+1][numPtV0]={0};
  Float_t ScalingFactorBulk[nummolt+1][numPtV0]={   0};
  Float_t ScalingFactorInclusive[nummolt+1][numPtV0]={   0};

  TH1F* hDeltaEtaLimits;
  TH1F* hDeltaEtaLimitsDef;
  Float_t  DeltaEtaInclusive=0;
  Float_t  DeltaEtaJet=0;
  Float_t  DeltaEtaJetDef=0;

  for (Int_t LoopFile=0; LoopFile<=1; LoopFile++){
    if (LoopFile==0)   numPtV0MaxUniversal=numPtV0Max;
    else  numPtV0MaxUniversal=numPtV0MaxOOJ;
    //if (LoopFile==1) continue;
    cout << "\n\n\n Looping on the two files... " << LoopFile << endl;
    //**************calcolo numero particelle di trigger*******************                                                                 
    //*********************************************************************                                                                
    if (LoopFile==0)    PtTrigChosen=   PtTrigMin1;
    else     PtTrigChosen=   PtTrigMin2;
    if (LoopFile==0)    PtTrigMax=   PtTrigMax1;
    else     PtTrigMax=   PtTrigMax2;
                  
    TString yearFinal;             
    TString yearMCFinal;             
    if (LoopFile==0){
      yearFinal=yearPtTrig3;
      yearMCFinal=yearMCPtTrig3;
    }
    else {
      yearFinal=year;
      yearMCFinal=yearMC;
    }
    Int_t NTrigger[nummolt+1][2];
    TFile *fileinbis[2]; //Data and MC
    TString PathInBis[2];
    PathInBis[0] =  "FinalOutput/AnalysisResults" + yearFinal +".root"; //+ Path1  + ".root"; change
    PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMCFinal  + "_MCEff.root";// + Path1 +".root"; change
    if (ishhCorr){  
      PathInBis[0] =  "FinalOutput/AnalysisResults" + year +Path1+ "_hhCorr.root";
      PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMC  + Path1+"_hhCorr_MCEff.root";
    }

    TFile *fileinbisPart1;
    TFile *fileinbisPart2;
    TString PathInBisPart1;
    TString PathInBisPart2;
    TDirectoryFile *dirPart1;
    TDirectoryFile *dirPart2;
    TList *listPart1;
    TList *listPart2;

    for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){ //change
      if (!IsSpecial){
	fileinbis[IntisMC]=new TFile(PathInBis[IntisMC],"");
	cout <<"file from task: " << fileinbis[IntisMC] << endl;
	cout << "Path in " << PathInBis[IntisMC] << endl;
	if (!fileinbis[IntisMC]) return;
      }
      else {
	PathInBisPart1 =  "FinalOutput/AnalysisResults2016kehjl_hK0s.root"; //+ Path1  + ".root"; change      
	PathInBisPart2 =  "FinalOutput/AnalysisResultsLHC17_hK0s.root"; //+ Path1  + ".root"; change          
	fileinbisPart1=new TFile(PathInBisPart1,"");
	fileinbisPart2=new TFile(PathInBisPart2,"");
	cout << "Path in " << PathInBisPart1 << " and " << PathInBisPart2 << endl;
      }

      TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Lambda", "Xi",   "Xi", "Omega","Omega", "Xi", "Omega"};
      TDirectoryFile *dir;
      if (isNewInputPath){
	if (yearPtTrig3=="17pq_hXi"){
	  if (LoopFile==0)	  dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTaskXi_PtTrigMin3.0_PtTrigMax15.0");
	  else 	  dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTaskXi_PtTrigMin0.2_PtTrigMax2.5");
	}
	else  {
	  dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTaskXi_PtTrigMin3.0_PtTrigMax15.0");
	  if (isOOJFromK0s && LoopFile==1)  dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTask_PtTrigMin3.0_PtTrigMax15.0");
	}
	if (!dir) {cout << " input dir not present " << endl; return;}
      } else {
	if (IsSpecial){
	  dirPart1 = (TDirectoryFile*)fileinbisPart1->Get("MyTask"+dirinputtype[type]);
	  dirPart2 = (TDirectoryFile*)fileinbisPart2->Get("MyTask"+dirinputtype[type]);
	}
	else {
	  dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTask"+dirinputtype[type]);
	  if (!dir) return;
	}
      }

      TString NameContainer= "";
      if (isNewInputPath) {
	if (yearPtTrig3=="17pq_hXi"){
	  NameContainer = "_hXi_Task_";
	}
	else {
	  if (LoopFile==0) NameContainer = "_hXi_Task_Default";
	  else {
	    if (isOOJFromK0s)  NameContainer = "_hK0s_Task_";
	    else  NameContainer = "_hXi_Task_LowPtTrig";
	  }
	}
      }

      TList *list;
      if (!IsSpecial){
	list = (TList*)dir->Get("MyOutputContainer" + NameContainer);
	if (!list) {cout << "input list not present " << endl; return;}
      }
      else {
	listPart1 = (TList*)dirPart1->Get("MyOutputContainer"+NameContainer);
	listPart2 = (TList*)dirPart2->Get("MyOutputContainer"+NameContainer);
	if (!listPart1 || !listPart2) return;
      }

      TH2D *fHistTriggervsMult;
      TH2D *fHistTriggervsMult2;
      if (!IsSpecial){
	fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtMaxvsMultBefAll");
	if (!fHistTriggervsMult ) return;
      }
      else   {
	cout << " I'm special ! " << endl;
	fHistTriggervsMult2         = (TH2D*)listPart2->FindObject("fHistPtMaxvsMultBefAll");
	if (!      fHistTriggervsMult2 ) return;
	fHistTriggervsMult2         ->SetName("fHistPtMaxvsMultBefAllDenom");
	cout << "entries 1; " << fHistTriggervsMult2->GetEntries() << endl;
	fHistTriggervsMult         = (TH2D*)listPart1->FindObject("fHistPtMaxvsMultBefAll");
	if (!      fHistTriggervsMult ) return;
	cout << "entries 2 " << fHistTriggervsMult->GetEntries() << endl;
	fHistTriggervsMult ->Add(      fHistTriggervsMult2) ;
	cout << "sum ; " << fHistTriggervsMult->GetEntries() << endl;

      }

      TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigChosen+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMax-0.00001) );

      //*********************************************************************                                                                                                
 
      for(Int_t m=0; m<nummolt+1; m++){
	if (LoopFile==1 && year == "161718Full_AOD234_hXi"){
	  Nmolt[m] = Nmolt0[m];
	  Smolt[m] = Smolt0[m];
	}
	if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
	if (isppHM && (m==0 || m==1)) continue;
	//      if (m==0) continue;
	NTrigger[m][IntisMC]=0;
	if(m<nummolt){
	  for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.00001); j<=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.00001); j++ ){
	    NTrigger[m][IntisMC]+= fHistTriggervsMult_MultProj->GetBinContent(j);
	  }
	}
	else {
	  for(Int_t j=1; j<=fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	    NTrigger[m][IntisMC]+= fHistTriggervsMult_MultProj->GetBinContent(j);
	  }
	}
	cout << "n trigger in mult range (all triggers)    " << m << "  " <<  " is MC " << IntisMC << "  " << NTrigger[m][IntisMC] <<   endl;
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////              

    gStyle->SetOptStat(0);
 
    TFile *filein[numPtTrig];
    TFile *fileinDefaultEta[numPtTrig];
    Int_t counter=0;

    //tentative  for (Int_t fileinit=0; fileinit< 2; fileinit++){
    for(Int_t m=nummolt; m>=0; m--){
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      if (isppHM && (m==0 || m==1)) continue;

      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){

	PtTrigMin=PtTrig+PtTrigChosen;
	if (PtTrigMin!=PtTrigChosen) continue;
	counter++;
	cout << "PtTrig " << PtTrig << endl;

	HistoWidthGaussian[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
	HistoWidthGaussian[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
	HistoWidthGaussianAS[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianASm%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
	HistoWidthGaussianAS[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianASm%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
	HistoWidthGaussianEta[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianEtam%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
	HistoWidthGaussianEta[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianEtam%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);

      }


      for(Int_t v=PtV0Min; v<numPtV0Max; v++){ 
	cout << " Loop file " << LoopFile << " v " << v << " " << NPtV0[v] <<endl; 
	if (v>= numPtV0MaxOOJ && LoopFile==1) {
	  cout << "v is " << v << " and I am going to continue " << endl;
	  continue;
	}
	//    if (v!=PtIntervalShown) continue;
	//	if (v>4) continue;
	nameSE[m][v]="ME_";
	nameSE[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_AC";
	nameME[m][v]="ME_m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_norm";
	nameMEEtaProj[m][v] =nameME[m][v]+"_EtaProj";
	nameMEPhiProj[m][v] =nameME[m][v]+"_PhiProj"; 
	namePhiProjJet[m][v]= nameSE[m][v] + "_phi_V0Sub_BulkSub_EffCorr";
	namePhiProjJetZYAM[m][v]= nameSE[m][v] + "_phi_V0Sub_EffCorr_BulkSubZYAM";
	namePhiProjJetFromBulkFit[m][v]= nameSE[m][v] + "_phi_V0Sub_EffCorr_BulkSubBulkFit";
	namePhiProjJetNotBulkSub[m][v]= nameSE[m][v] + "_phi_V0Sub_EffCorr";
	namePhiProjBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_Bulk_EffCorr";
	namePhiProjJetBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_JetBulkEffCorr";
	//      cout << nameSE[m][v]<< endl;

	for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	  PtTrigMin=PtTrig+PtTrigChosen;
	  if (PtTrigMin!=PtTrigChosen) continue;
	  cout << "PtTrig " << PtTrig << endl;
	  cout << "PtTrigMin " << PtTrigMin << endl;
	  cout << "PtTrigChosen " << PtTrigChosen << endl;
	  cout << "PtTrigMax " << PtTrigMax << endl;
	  nameEtaProj[m][v][PtTrig]= nameSE[m][v] + Form("_EtaProj_PtTrig%i", PtTrig);
	  //	if (PtTrigMin==4 || PtTrigMin==5 || PtTrigMin>10)continue;
	  if (PtTrigMin>7) continue;
	  for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
	    PathIn[PtTrig]= Dir+"/histo/AngularCorrelation";
	    if (LoopFile==0) PathIn[PtTrig]  += file[IntisMC];
	    //if (!(LoopFile==0 && NPtV0[v]>2.5))  PathIn[PtTrig]  += Path1;
	    //	    if (!(LoopFile==0))  PathIn[PtTrig]  += Path1;
	    if (LoopFile==1)	    PathIn[PtTrig]  += fileOOJ[IntisMC];
	    cout << "PathIn " << PathIn[PtTrig] << endl;
	    if(type>=0){
	      //	      if (LoopFile==1)	      PathIn[PtTrig] += Form("_PtTrigMax%.1f",PtTrigMax);
	      if (LoopFile==1 && isOOJFromK0s){
		PathIn[PtTrig] +="_"+tipo[0];
	      }
	      else PathIn[PtTrig] +="_"+tipo[type];
	      PathIn[PtTrig] +=Srap[israp];
	      if (LoopFile==1 && isOOJFromK0s){
		PathIn[PtTrig] +=BkgType[isBkgParab];
	      }
	      //	      if (LoopFile==0 && NPtV0[v]>2.5)	      PathIn[PtTrig] +=SSkipAssoc[SkipAssoc];
	      //if (LoopFile==0)	      PathIn[PtTrig] +=SSkipAssoc[SkipAssoc];
	      PathIn[PtTrig] +=SSkipAssoc[SkipAssoc];
	    }
	    PathInDefaultEta[PtTrig]= PathIn[PtTrig];
	    if (LoopFile==0){
	      PathInDefaultEta[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, 0, PtTrigChosen)+"_Output";
	      if (isEtaEff) PathInDefaultEta[PtTrig]+="_IsEtaEff";
	      if (MultBinning!=0) PathInDefaultEta[PtTrig]+=Form("_MultBinning%i", MultBinning);
	      PathInDefaultEta[PtTrig] += ".root";

	      if(isSysDef && isDefaultSel)   PathIn[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysT%i_SysV0Default_Sys%i_PtMin%.1f", sysTrigger, sysang, PtTrigChosen)+"_Output";   
	      else  if(isSysDef && !isDefaultSel && sysTrigger==0)   PathIn[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysT%i_SysV0index%i_Sys%i_PtMin%.1f", sysTrigger, indexSysV0, sysang, PtTrigChosen)+"_Output";   
	      else  if(isSysDef && !isDefaultSel && sysTrigger==1)   PathIn[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysTindex%i_SysV0%i_Sys%i_PtMin%.1f", indexsysTrigger, 0, sysang, PtTrigChosen)+"_Output";   
	      else  {
		PathIn[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigChosen)+"_Output";
		if (isEPOS) PathIn[PtTrig]+="_EPOS";
	      }
	      if (isMEFrom13TeV) PathIn[PtTrig]+= "_IsMEFrom13TeV";
	      if (isEtaEff) PathIn[PtTrig]+="_IsEtaEff";
	      if (MultBinning!=0) PathIn[PtTrig]+=Form("_MultBinning%i", MultBinning);
	      PathIn[PtTrig] += ".root";
	    }
	    else{
	      PathIn[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", 0, sysV0, sysOOJ, PtTrigChosen)+"_Output";
	      if (isEtaEff) PathIn[PtTrig]+="_IsEtaEff";
	      if (MultBinning!=0) PathIn[PtTrig]+=Form("_MultBinning%i", MultBinning);
	      PathIn[PtTrig]+=".root";
	    }
	    if (!ishhCorr && isEnlargedDeltaEta)	  PathIn[PtTrig]= Dir+"/histo/AngularCorrelation" + file[IntisMC] + hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin)+"_DeltaEtaPhiEnlarged_Output.root";
	    cout <<"path in : " <<  PathIn[PtTrig] << endl;
	    cout <<"path in default : " <<  PathInDefaultEta[PtTrig] << endl;
	    //	    return;
	    filein[PtTrig]= new TFile(PathIn[PtTrig], "");
	    if (!filein[PtTrig]) {cout << PathIn[PtTrig]  << " is not there! " << endl; return;}

	    if (LoopFile==0){
	      fileinDefaultEta[PtTrig]= new TFile(PathInDefaultEta[PtTrig], "");
	      if (!fileinDefaultEta[PtTrig]) {cout << PathInDefaultEta[PtTrig]  << " is not there! " << endl; return;}
	      hDeltaEtaLimits=(TH1F*)	    filein[PtTrig]->Get("fHistEtaLimitsOfRegion");
	      hDeltaEtaLimitsDef=(TH1F*)	    fileinDefaultEta[PtTrig]->Get("fHistEtaLimitsOfRegion");
	      DeltaEtaInclusive = 	  2*  hDeltaEtaLimits->GetBinContent(6);
	      DeltaEtaJet = 	  2*  hDeltaEtaLimits->GetBinContent(2);
	      DeltaEtaJetDef = 	  2*  hDeltaEtaLimitsDef->GetBinContent(2);
	    }
	    hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameSE[m][v]);
	    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameME[m][v]);

	    if (v==PtV0Min)  hDeltaEtaDeltaPhi_MEbins[m][numPtV0MaxUniversal-1][PtTrig]= (TH2F*)filein[PtTrig]->Get("ME_m"+ Smolt[m]+"_v"+SPtV0[numPtV0MaxUniversal-1]+Ssideband[0]+"_norm");
	    if (!hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]){cout << "no SE 2D histo for v = "<<v << " name of the histo " << nameSE[m][v] << endl;  return;}
	    if (!hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]){cout << "no ME 2D histo for v = " << v << endl;  return;}
	    if (!hDeltaEtaDeltaPhi_MEbins[m][numPtV0MaxUniversal-1][PtTrig]){cout << "no ME 2D histo for v = " <<numPtV0Max-1 <<  endl;  return;}

	    //projection along eta of the Mixed Event distribution
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->ProjectionX(nameMEEtaProj[m][v],0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetMarkerColor(ColorPt[v]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetNbinsY());

	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Clone(nameMEEtaProj[m][v]+"_Ratio");
	    //	  if (v==1)	{
	    hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][numPtV0MaxUniversal-1][PtTrig]->ProjectionX(nameMEEtaProj[m][numPtV0MaxUniversal-1]+ "_Master",0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetNbinsY());
	    //	  }
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]);

	    for (Int_t i=1; i<=   hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetNbinsX(); i++){
	      hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->SetBinError(i, sqrt(TMath::Abs(pow( hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetBinError(i),2)-pow( hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->GetBinError(i),2)))/ hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->GetBinContent(i));
	    }


	    cout << "ho preso isto" << endl;

	    //projection along phi of the Mixed Event distribution
	    hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->ProjectionY(nameMEPhiProj[m][v],0,-1, "E");

	    hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	    hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->SetMarkerColor(ColorPt[v]);

	    hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->Clone(nameMEPhiProj[m][v]+"_Ratio");
	    if (v==PtV0Min)	  hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]= 	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][numPtV0MaxUniversal-1][PtTrig]->ProjectionY(nameMEPhiProj[m][numPtV0MaxUniversal-1]+"_Master",0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]);

	    for (Int_t i=1; i<=   hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->GetNbinsX(); i++){
	      hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->SetBinError(i,sqrt( TMath::Abs( pow(hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->GetBinError(i),2)- pow(hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]->GetBinError(i),2)))/hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]->GetBinContent(i));
	    }
	    cout << "ho preso isto" << endl;

	    //other projections
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetZYAM[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetFromBulkFit[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetNotBulkSub[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]);

	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]+"_RelError");
	    hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]+"_RelError");
	    hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]+"_RelError");
	    hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig], 0,-1, "E");
	    hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig],hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(+TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(+3./2*TMath::Pi()- 0.001)  , "E");
	    hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig],hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(-TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(+1./2*TMath::Pi()- 0.001)  , "E");


	    cout << "histos were searched for in "<< filein[PtTrig]->GetName() << endl;
	    if(! hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]) {cout << " missing histo jet" << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]) {cout << " missing histo jet zyam " << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]) {cout << " missing histo bulk from fit" << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]) {cout << " missing histo jetnotbulksub " << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]) {cout << " missing histo bulk" << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]) {cout << " missing histo jetbulk" << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]) {cout << " missing histo jet rel error" << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]) {cout << " missing histo bulk rel error" << endl; return;}
	    if(! hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]) {cout << " missing histo jetbulkrel error " << endl; return;}

	    cout << "ho preso isto" << endl;

	    //	    cout << "this rebin is performed for visualization purposes only " << endl;
	    if (m!=nummolt){
	      //	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Rebin(2);
	      //	      hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Rebin(2);
	      //	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Rebin(2);
	    }
	    
	    //canvasPlot[m]->cd(PtTrig+1);
	    canvasPlot[m][IntisMC]->cd(v+1);
	    cout << "scelto cd " << endl;
	    //	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetLogz();
	    //	gPad->SetLogz();
	    TString TitleHisto=Form("%.1f < p_{T} < %.1f GeV/c ", NPtV0[v], NPtV0[v+1]) + Smolt[m] +"%";
	    //"Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin) +MCOrNot[IntisMC]
	    hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetTitle("AC v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->SetTitle(TitleHisto);
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->SetTitle(TitleHisto);
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->SetTitle(TitleHisto);
	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetTitle(TitleHisto);
	    hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->Draw("colz");
	    tlinePhiDx->Draw();
	    tlinePhiSx->Draw();
	    tlineEtaDx->Draw();
	    tlineEtaSx->Draw();
	    tlineEtaInclDx->Draw();
	    tlineEtaInclSx->Draw();

	    canvasPlotME[m][IntisMC]->cd(v+1);
	    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetTitle("AC v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->Draw("colz");

	    cout << "m " << m << " v " << v << " "  << counter << endl;
	    if (counter==1)	  legendPt->AddEntry(	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig], SPtV0[v], "pl");
	    canvasPlotMEEtaProj[IntisMC]->cd(m+1);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Draw("same");
	    if (v==numPtV0-1) legendPt->Draw("same");
	    canvasPlotMEEtaProj[IntisMC]->cd(m+1+6);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->Draw("same");
	    if (v==numPtV0-1) legendPt->Draw("same");

	    canvasPlotMEPhiProj[IntisMC]->cd(m+1);
	    hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,20);
	    hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->Draw("same");
	    if (v==numPtV0-1) legendPt->Draw("same");
	    canvasPlotMEPhiProj[IntisMC]->cd(m+1+6);
	    hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	    hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->Draw("same");
	    if (v==numPtV0-1) legendPt->Draw("same");

	    //	canvasPlot[m]->cd(PtTrig+1+7);
	    cout << "setting to zero last bin " << endl;
	    //	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	    //	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	    //	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);


	    hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->SetLineColor(860);
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetLineColor(kBlue);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->SetLineColor(kGreen+3);

	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetLineColor(628);
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetLineColor(418);

	    hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->SetMarkerColor(860);
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetMarkerColor(628);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->SetMarkerColor(kGreen+3);

	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetMarkerColor(628);
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetMarkerColor(418);
	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]    ->SetLineColor(kBlue-3);
	    //	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetMarkerColor(868);
	    cout << "setting to zero last bin 1" << endl;
	    if (hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig] && hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]&& hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]){
	      cout << "setting to zero last bin 2" << endl;
	      hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->SetLineColor(628);
	      hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->SetLineColor(418);
	      hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->SetLineColor(868);
	      hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->SetMarkerColor(628);
	      hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->SetMarkerColor(418);
	      hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->SetMarkerColor(868);
	      hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.07);
	      hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->GetYaxis()->SetLabelSize(0.07);
	      hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->GetYaxis()->SetLabelSize(0.07);
	    }
	    //division by number of trigger particles
	    cout << "division by trigger particles " << endl;
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Scale(1./NTrigger[m][IntisMC]);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->Scale(1./NTrigger[m][IntisMC]);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->Scale(1./NTrigger[m][IntisMC]);
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->Scale(1./NTrigger[m][IntisMC]);
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Scale(1./NTrigger[m][IntisMC]);
	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Scale(1./NTrigger[m][IntisMC]);

	    //sum of different ptV0 bins
	    if (v==PtV0Min){
	      cout << "Adding histos to make ptv0 sum" << endl;
	      hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    =	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Clone(Form("PhiProjJet_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	      hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   =	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Clone(Form("PhiProjBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	      hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]=	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Clone(Form("PhiProjJetBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
    
	    }

	    hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->SetTitle(Form("PhiProjJet_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->SetTitle(Form("PhiProjBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));	  
	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->SetTitle(Form("PhiProjJetBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));

	    if (v!=1){
	      cout << "Adding histos to make ptv0 sum" << endl;
	      hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]);
	    }

	    gPad->SetLeftMargin(0.2);

	    Gaussian[m][v][PtTrig]= new TF1 ( Form("GaussianFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  -0.8, 0.8);
	    GaussianAS[m][v][PtTrig]= new TF1 ( Form("GaussianASFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  TMath::Pi()-1.2,   TMath::Pi()+1.2);
	    GaussianEta[m][v][PtTrig]= new TF1 ( Form("GaussianEtaFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  -0.8,   0.8);
	    Baseline[m][v][PtTrig]= new TF1 ( Form("Baseline_m%i_v%i_PtTrig%i", m,v,PtTrig),"pol0",  1.5,  4);
	    Baseline[m][v][PtTrig]->SetLineColor(881);
	    GaussianAS[m][v][PtTrig]->SetLineColor(419);
	    GaussianEta[m][v][PtTrig]->SetLineColor(860);

	    cout << "m " << m << " v " << v << " PtTrig " << PtTrig << "  " << MCOrNot[IntisMC] << endl;
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Fit(Gaussian[m][v][PtTrig], "R");
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Fit(Baseline[m][v][PtTrig], "R0");
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Fit(GaussianAS[m][v][PtTrig], "R0");
	    if (v>1 && v <4) GaussianEta[m][v][PtTrig]->SetRange(-0.5, 0.5);
	    else if (v>=4)	  GaussianEta[m][v][PtTrig]->SetRange(-0.3, 0.3);
	    hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Fit(GaussianEta[m][v][PtTrig], "R");

	    cout << "\nbaseline value" <<  Baseline[m][v][PtTrig]->GetParameter(0)<< " +- " <<  Baseline[m][v][PtTrig]->GetParError(0)<< endl;
	    cout << "\njet fit: chi square " << Gaussian[m][v][PtTrig]->GetChisquare() << " NDF " << Gaussian[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << Gaussian[m][v][PtTrig]->GetChisquare()/Gaussian[m][v][PtTrig]->GetNDF() << endl;
	    cout << "\nAwaySide fit: chi square " << GaussianAS[m][v][PtTrig]->GetChisquare() << " NDF " << GaussianAS[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << GaussianAS[m][v][PtTrig]->GetChisquare()/GaussianAS[m][v][PtTrig]->GetNDF() << endl;

	    //	  if (	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetMaximum() > 	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetMaximum()) 	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Draw("");

	    if (m<=1 && MultBinning==1) LimSupY[0]=0.001;
	    if (type==0)	  LimSupY[0]=0.03;
	    if (type==0)	  LimInfY[0]=-0.002;
	    if (type==8){
	      if (year.Index("pp5TeV")!=-1){
		if (m==1) LimSupY[0]=0.0005;
		//		LimSupY[0]=0.0001;
		if (m==0) LimSupY[0]=0.0008;
	      }
	      else {
		if (m==0 || m==1)  LimSupY[0]=0.0012;
	      }
	      if (isppHM) {
		LimSupY[0]=0.0012;
	      }
	    }
	    cout << " ciao ciao " << endl;
	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]); //-0.004 for hh and hK0s
	    cout << " ciao ciao " << endl;
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(LimInfY[ishhCorr],LimSupY[ishhCorr]);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]);
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]);
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]);

	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetLabelSize(0.04);
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->GetYaxis()->SetLabelSize(0.04);
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.04);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.04);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.04);
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetLabelSize(0.04);

	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetTitle("N/N_{Trigg} per #Delta#eta"); 
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetTitle("N/N_{Trigg} per #Delta#eta"); 
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetTitle("N/N_{Trigg} per #Delta#eta"); 
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]   ->GetYaxis()->SetTitle("N/N_{Trigg} per #Delta#eta"); 

	    hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.5); 
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetTitleOffset(1.5); 
	    hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetTitleOffset(1.5); 
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]   ->GetYaxis()->SetTitleOffset(1.5); 

	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetMarkerColor(kBlue);
	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetLineColor(kBlue);
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetMarkerColor(628);
	    hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetLineColor(628);

	    pol0ZYAM[m][v][PtTrig]= new TF1("pol0",Form("pol0ZYAM_m%i_v%i",m,v), 1,2);
	    pol0ZYAM[m][v][PtTrig]->SetLineColor(kRed);
	    //	    hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Fit(pol0ZYAM[m][v][PtTrig], "R+");

	    pol0BulkBis[m][v][PtTrig]= new TF1("pol0",Form("pol0BulkBis_m%i_v%i",m,v), -1,1);
	    pol0BulkSmooth[m][v][PtTrig]= new TF1("pol0",Form("pol0BulkSmooth_m%i_v%i",m,v), -TMath::Pi()/2,3/2*TMath::Pi());
	    pol0BulkSmooth[m][v][PtTrig]->SetLineColor(418);
	    if (LoopFile==0)	    pol0BulkBis[m][v][PtTrig]->SetLineColor(418);
	    else pol0BulkBis[m][v][PtTrig]->SetLineColor(420);
	    
	    cout << " ciao ciao " << endl;
	    if (LoopFile==0){
	      ScalingFactorInclusive[m][v]=	      hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Integral( hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->FindBin(1),  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->FindBin(2));
	      IntegralFactorBulkDef[m][v][PtTrig]=	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Integral();
	      ScalingFactorJetNotBulkSub[m][v]=	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Integral( hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->FindBin(1),  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->FindBin(2));
	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethod[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Clone(namePhiProjJetNotBulkSub[m][v] + "_DefaultMethod");
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_NewMethod");
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_NewMethodRebinSub");
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]=	(TH1F*)      	      hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_DefaultMethod");
	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig] =(TH1F*)	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Clone(namePhiProjJetNotBulkSub[m][v]+"_Rebin");
	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubRebin[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Clone(namePhiProjJetNotBulkSub[m][v] + "_Rebinned");
	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubRebin[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig] ->Rebin(2);

	      if (PtBinning==0)	      canvasPlotProj[m][IntisMC]->cd(v);
	      else       canvasPlotProj[m][IntisMC]->cd(v+1);
	      gPad->SetLeftMargin(0.15);
	      if (v<numPtV0MaxOOJ){
		hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Draw("hist p e");
		//	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethod[m][v][PtTrig]->Draw("");
		hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]    ->SetMarkerStyle(33);
		hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]    ->Draw("same");
		hDeltaEtaDeltaPhi_PhiProjBulkReb[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Clone(namePhiProjBulk[m][v]+"_Rebin");
		hDeltaEtaDeltaPhi_PhiProjBulkReb[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjBulkReb[m][v][PtTrig]->Draw("same");
		tlinePhiBase->Draw("same");
	      }
	      if (PtBinning==0)	      canvasPlotProjBis[m][IntisMC]->cd(v);
	      else       canvasPlotProjBis[m][IntisMC]->cd(v+1);
	      if (v<numPtV0MaxOOJ){
		hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8*hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetBinContent( hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetMinimumBin()), 1.2*hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetBinContent( hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetMaximumBin()));
		hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Draw("hist p e");
		hDeltaEtaDeltaPhi_PhiProjBulkReb[m][v][PtTrig]->Draw("same");
	      }
	      if (PtBinning==0)	      canvasPlotProjBis[m][IntisMC]->cd(v+numC);
	      else       canvasPlotProjBis[m][IntisMC]->cd(v+1+numC);
	      if (v<numPtV0MaxOOJ){
		//		hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetYaxis()->SetRangeUser(3*hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetBinContent( hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetMinimumBin()), 3*hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetBinContent( hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetMaximumBin()));
		hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetYaxis()->SetRangeUser(-0.01*10e-3, 0.02*10e-3);
		cout <<  hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetBinContent( hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetMaximumBin())<< endl;
		hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]    ->DrawClone("same");
		tlinePhiBase->Draw("same");
	      }

	      if (PtBinning==0)	      canvasPlotProj[m][IntisMC]->cd(v+numC);
	      else	      canvasPlotProj[m][IntisMC]->cd(v+numC+1);
	      gPad->SetLeftMargin(0.15);
	      if (v<numPtV0MaxOOJ){
		hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethod[m][v][PtTrig]->Draw("");
		tlinePhiBase->Draw("same");
	      }

	      if (PtBinning==0)	      canvasPlotProj[m][IntisMC]->cd(v+numC);
	      else	      canvasPlotProj[m][IntisMC]->cd(v+2*numC+1);
	      gPad->SetLeftMargin(0.15);
	      if (v<numPtV0MaxOOJ){
		hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Draw("");
		hDeltaEtaDeltaPhi_PhiProjJetInclusive[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_InclusiveRebSmooth");	 
		tlinePhiBase->Draw("same");
	      }

	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_RebSmooth");	 
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothCorrMult[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_RebSmoothCorrMult");	 
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_RebSmoothBis");	 
	      hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_RebCorrMultBis");	 
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_RebSmoothFit");	 
	
	      //subtraction using rebin + smooth OOJ distribution from events with pT,Trig > 3 geV 
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Clone(namePhiProjJet[m][v] + "_DefaultMethodSmooth");
	      hDeltaEtaDeltaPhi_PhiProjBulkRebSmoothed[m][v][PtTrig]= (TH1F*)	      	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Clone(namePhiProjBulk[m][v] + "_RebSmooth");	 
	      hDeltaEtaDeltaPhi_PhiProjBulkRebSmoothed[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkRebSmoothed[m][v][PtTrig]->Smooth();
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->Add(	      hDeltaEtaDeltaPhi_PhiProjBulkRebSmoothed[m][v][PtTrig],-1);

	      for (Int_t b=1; b<= hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->GetNbinsX(); b++ ){
		//		hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->SetBinError(b, TMath::Abs(	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetBinError(b)-hDeltaEtaDeltaPhi_PhiProjBulkRebSmoothed[m][v][PtTrig]->GetBinError(b)) );

	      }

	      canvasPlotProjAllPt[m][IntisMC]->cd(v+1);
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->SetLineColor(881);
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->SetMarkerColor(881);
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjBulkRebSmoothed[m][v][PtTrig]->Draw("");
	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Draw("same");
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->Draw("same");
	      tlinePhiBase->Draw("same");
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]    ->Draw("same");
	      if (PtBinning==0)	      canvasPlotProj[m][IntisMC]->cd(v);
	      else 	      canvasPlotProj[m][IntisMC]->cd(v+1);
	      if (v<numPtV0MaxOOJ){
		//c	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->DrawClone("same");
	      }

	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]->Scale(NTrigger[m][IntisMC]);
	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetDefaultMethodSmooth[m][v][PtTrig]);
	    }
	    //	    }	    
	    //I rescale the histogram of the OOJ distribution

	    if (LoopFile==1)	 {
	      if (PtBinning==0)	      canvasPlotProj[m][IntisMC]->cd(v+numC);
	      else	      canvasPlotProj[m][IntisMC]->cd(v+numC+1);
	      ScalingFactorBulk[m][v]=	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Integral( hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->FindBin(1),  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->FindBin(2));
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Sumw2();
	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "Scaled");
	      if (m==nummolt)  hDeltaEtaDeltaPhi_PhiProjBulk_mall[v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulk[nummolt][v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "_mall");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulk_mall[v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "ScaledAllMult");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledToInclusiveAllMult[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulk_mall[v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "ScaledToInclusive");

	      hDeltaEtaDeltaPhi_PhiProjBulkScaledToInclusiveAllMult[m][v][PtTrig]->Scale(ScalingFactorInclusive[m][v]/ScalingFactorBulk[nummolt][v]);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]->Scale(ScalingFactorJetNotBulkSub[m][v]/ScalingFactorBulk[m][v]);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->Scale(ScalingFactorJetNotBulkSub[m][v]/ScalingFactorBulk[nummolt][v]);

	      tlineAtOne->Draw("same"); 
	      //cout <<"\\n\n\n****** m " << m << " v " << v <<" " <<  	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->GetBinContent(1)<< endl;
	      //	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->Scale(0.02/ScalingFactorBulk[nummolt][v]);
	      //	      cout << "scaling factor " << ScalingFactorBulk[nummolt][v]<< endl;
	      //	      cout <<  	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->GetBinContent(1)<< endl;

	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]->SetLineColor(420);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]->SetMarkerColor(420);

	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->SetLineColor(kGreen-6);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->SetMarkerColor(kGreen-6);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->SetMarkerStyle(33);

	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]->Draw("same");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->Draw("same");

	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig],-1);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->SetLineColor(kRed+2);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->SetMarkerColor(kRed+2);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->SetMarkerStyle(33);
	      cout << " \n\n\n rebinning " << endl;
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->Draw("same");
	
	      //the histogram above is equal to the histogram 'JetNewMethodRebinSub' since rebinning -> subtracting == subtracting --> rebinning (only a difference in the error if histograms involve in subtraction are treated as fully correlated. In this case I treat them as uncorrelated (...not correct...) and therefore there are no differences

	      hDeltaEtaDeltaPhi_PhiProjBulkScaledRebin[m][v][PtTrig] = (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "ScaledRebin");
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledRebin[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaledRebin[m][v][PtTrig],-1);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->SetLineColor(kRed-7);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->SetMarkerColor(kRed-7);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->Draw("same");

	      hDeltaEtaDeltaPhi_PhiProjJetNDRatio[m][v][PtTrig]=(TH1F*)	      hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->Clone(namePhiProjJet[m][v]+"_NDRatio");
	      hDeltaEtaDeltaPhi_PhiProjJetNRebinDRatio[m][v][PtTrig]=(TH1F*)	      hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->Clone(namePhiProjJet[m][v]+"_NRebinDRatio");
	      hDeltaEtaDeltaPhi_PhiProjJetNDRatio[m][v][PtTrig]->Divide(hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjJetNRebinDRatio[m][v][PtTrig]->Divide(hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]);

	      for (Int_t b=1; b<=		hDeltaEtaDeltaPhi_PhiProjJetNDRatio[m][v][PtTrig]->GetNbinsX() ; b++){
		hDeltaEtaDeltaPhi_PhiProjJetNDRatio[m][v][PtTrig]->SetBinError(b,TMath::Abs(hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetBinError(b)-hDeltaEtaDeltaPhi_PhiProjJetNewMethod[m][v][PtTrig]->GetBinError(b)));
		hDeltaEtaDeltaPhi_PhiProjJetNRebinDRatio[m][v][PtTrig]->SetBinError(b,TMath::Abs(hDeltaEtaDeltaPhi_PhiProjJetDefaultMethod[m][v][PtTrig]->GetBinError(b)-hDeltaEtaDeltaPhi_PhiProjJetNewMethodRebinSub[m][v][PtTrig]->GetBinError(b)));
	      }

	      if (PtBinning==0)	      canvasPlotProj[m][IntisMC]->cd(v+2*numC);
	      else	      canvasPlotProj[m][IntisMC]->cd(v+2*numC+1);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultRebin[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMult[m][v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "ScaledAllMultRebin");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultRebin[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultRebin[m][v][PtTrig]->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]=(TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultRebin[m][v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "ScaledAllMultSmoothed");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]-> Smooth();
	      pol0BulkSmooth[m][v][PtTrig]->SetRange(-1, 1);
	      //c	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]->Fit(pol0BulkSmooth[m][v][PtTrig], "R+");
	      pol0BulkSmooth[m][v][PtTrig]->SetRange(-1./2*TMath::Pi(), 3./2*TMath::Pi());
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig],-1);

	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]->Smooth();
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothCorrMult[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig],-1);

	      for (Int_t b=1; b<= hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->GetNbinsX(); b++ ){
		//		hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->SetBinError(b, TMath::Abs(	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetBinError(b)-hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]->GetBinError(b)) );
		//		hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]->SetBinError(b, TMath::Abs(	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->GetBinError(b)-pol0BulkSmooth[m][v][PtTrig]->GetParError(0)));

	      }
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->SetLineColor(628);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->SetMarkerColor(628);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]->SetLineColor(kRed-7);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]->SetMarkerColor(kRed-7);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]->SetMarkerStyle(33);

	      //OOJ distribution for NoTriggerEvents and m0-100 scaled by integral of OOJ distribution default
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulk_mall[v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "ScaledBisAllMult");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]= (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "ScaledBis");
	      IntegralFactorBulkNewAllMult[m][v][PtTrig]=	      hDeltaEtaDeltaPhi_PhiProjBulk_mall[v][PtTrig]   ->Integral();
	      IntegralFactorBulkNew[m][v][PtTrig]=	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Integral();
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->Scale(IntegralFactorBulkDef[m][v][PtTrig]/IntegralFactorBulkNewAllMult[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]->Scale(IntegralFactorBulkDef[m][v][PtTrig]/IntegralFactorBulkNew[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->SetLineColor(kGreen+3);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->SetMarkerColor(kGreen+3);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]->SetLineColor(kViolet+5);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]->SetMarkerColor(kViolet+5);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]->SetMarkerStyle(33);
	      //	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->Smooth();
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]->Smooth();

	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig],-1);
	      hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig],-1);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]->SetLineColor(kAzure+7);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]->SetMarkerColor(kAzure+7);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]->SetLineColor(kBlue +2);
	      hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]->SetMarkerColor(kBlue+2);
	      hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]->SetMarkerStyle(33);

	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]->Scale(DeltaEtaJet/DeltaEtaJetDef); //for a correct sub of OOJ, JetNotBulkSUb was divided by the DeltaEta associated to that systematic choice. Here, the scaling by the default DeltaEta is recovered.
	      hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]->Scale(DeltaEtaJet/DeltaEtaJetDef);

	      pol0JetRebSmooth[m][v]= new TF1(Form("pol0JetRebSmooth_m%i_v%i", m, v), "pol0", 1.5, 3.5);
	      pol0JetRebSmoothBis[m][v]= new TF1(Form("pol0JetRebSmoothBis_m%i_v%i", m, v), "pol0", 1.5, 3.5);
	      pol0JetRebSmooth[m][v]->SetLineColor(628);
	      pol0JetRebSmoothBis[m][v]->SetLineColor(kAzure+7);
	      //	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->Fit(pol0JetRebSmooth[m][v], "R+");
	      //	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]->Fit(pol0JetRebSmoothBis[m][v], "R+");

	      //	      pol0BulkSmooth[m][v][PtTrig]->SetRange(-3/2*TMath::Pi(), TMath::Pi());
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]->Add(pol0BulkSmooth[m][v][PtTrig],-1, "i");
	      hDeltaEtaDeltaPhi_PhiProjJetNotBulkSubDefaultMethodRebin[m][v][PtTrig]->Draw("hist ep");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]->Draw("same hist ep");

	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->DrawClone("same ep");
	      //	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]->DrawClone("same hist ep");
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]->Draw("same ep");
	      //	      hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]->Draw("same ep");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->Draw("same hist ep");
	      //	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBis[m][v][PtTrig]->DrawClone("same hist ep");
	      tlinePhiBase->Draw("same"); 

	      //	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->Scale(NTrigger[m][IntisMC]);
	      //	      hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]->Scale(NTrigger[m][IntisMC]);
	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]);
	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetRebSmoothCorrMult[m][v][PtTrig]);
	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetRebSmoothBis[m][v][PtTrig]);
	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetRebSmoothFit[m][v][PtTrig]);
	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetRebCorrMultBis[m][v][PtTrig]);

	      if (PtBinning==0)	      canvasPlotProjBis[m][IntisMC]->cd(v);
	      else       canvasPlotProjBis[m][IntisMC]->cd(v+1);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]->SetMarkerColor(kViolet);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]->SetLineColor(kViolet);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultSmoothed[m][v][PtTrig]->Draw("same hist ep");
	      //hDeltaEtaDeltaPhi_PhiProjBulkScaledAllMultRebin[m][v][PtTrig]->Draw("same hist ep");

	      if (PtBinning==0)	      canvasPlotProjBis[m][IntisMC]->cd(v+numC);
	      else	      canvasPlotProjBis[m][IntisMC]->cd(v+numC+1);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->SetLineColor(kAzure+5);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->SetMarkerColor(kAzure+5);
	      hDeltaEtaDeltaPhi_PhiProjJetRebSmooth[m][v][PtTrig]->DrawClone("same ep");


	      if (PtBinning==0)	      canvasPlotProj[m][IntisMC]->cd(v+3*numC);
	      else	      canvasPlotProj[m][IntisMC]->cd(v+3*numC+1);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledToInclusiveAllMult[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledToInclusiveAllMult[m][v][PtTrig]->Smooth();
	     
	      //	      hDeltaEtaDeltaPhi_PhiProjJetInclusive[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaledToInclusiveAllMult[m][v][PtTrig],-1);
	      hDeltaEtaDeltaPhi_PhiProjJetInclusive[m][v][PtTrig]->Add(hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig],-1);
	      hDeltaEtaDeltaPhi_PhiProjJetInclusive[m][v][PtTrig]->Scale(DeltaEtaInclusive/DeltaEtaJet);
	      hDeltaEtaDeltaPhi_PhiProjJetInclusive[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjJetInclusive[m][v][PtTrig]->Draw("same");
	      // hDeltaEtaDeltaPhi_PhiProjBulkScaledToInclusiveAllMult[m][v][PtTrig]->Draw("same");	      
	      hDeltaEtaDeltaPhi_PhiProjBulkScaledBisAllMult[m][v][PtTrig]->Draw("same");
	      fileout->WriteTObject(hDeltaEtaDeltaPhi_PhiProjJetInclusive[m][v][PtTrig]);
	    }

	    //	  pol0ZYAM[m][v][PtTrig]->Draw("same");
	    //	  pol0BulkBis[m][v][PtTrig]->Draw("same");
	    
	    //	  hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->Draw("same");
	    //	  hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->Draw("same");
       

	    if (numC==4)	    canvasPlotOOJDistr[m][IntisMC][0]->cd(v);
	    else if (numC==6) 	    canvasPlotOOJDistr[m][IntisMC][0]->cd(v+1);
	    if (LoopFile==0)	{
	      //	      IntegralFactorBulkDef[m][v][PtTrig]=	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Integral();
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetLineColor(418);
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetMarkerColor(418);
	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig] = (TH1F*)hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig] ->Clone("PhiProjBulkRatioOOJ");
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetRangeUser(0,LimSupY[ishhCorr]);
	      hDeltaEtaDeltaPhi_PhiProjBulkDefault[m][v][PtTrig] = (TH1F*)hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig] ->Clone("PhiProjBulkDefault");
	      if (m==nummolt && v==PtV0Min) OOJlegend  ->AddEntry(hDeltaEtaDeltaPhi_PhiProjBulkDefault[m][v][PtTrig], "OOJ (w Trigg)", "pl");
	      hDeltaEtaDeltaPhi_PhiProjBulkDefault[m][v][PtTrig]   ->Draw("");
	    }
	    else{
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Fit(pol0BulkBis[m][v][PtTrig], "R+");
	      IntegralFactorBulkNew[m][v][PtTrig]=	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Integral();
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetLineColor(420);
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetMarkerColor(420);
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig] =(TH1F*) hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Clone(namePhiProjBulk[m][v]+ "_IntegralScaled");

	      hDeltaEtaDeltaPhi_PhiProjBulkDenom[m][v][PtTrig] = (TH1F*) hDeltaEtaDeltaPhi_PhiProjBulk[nummolt][v][PtTrig] ->Clone(namePhiProjBulk[m][v]+ "_DenomRatio");
	      //hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Scale(IntegralFactorBulkDef[m][v][PtTrig]/IntegralFactorBulkNew[m][v][PtTrig]);
	      //hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Scale(ScalingFactorJetNotBulkSub[m][v]/ScalingFactorBulk[m][v]);
	      //hDeltaEtaDeltaPhi_PhiProjBulkDenom[m][v][PtTrig]   ->Scale(IntegralFactorBulkDef[m][v][PtTrig]/IntegralFactorBulkNew[nummolt][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjBulkDenom[m][v][PtTrig]   ->Scale(ScalingFactorJetNotBulkSub[m][v]/ScalingFactorBulk[nummolt][v]);
 
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetRangeUser(0,LimSupY[ishhCorr]);
	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig] ->Divide(hDeltaEtaDeltaPhi_PhiProjBulkDenom[m][v][PtTrig]);
	      BulkRatioPol0[m][v] = new TF1(Form("bulkratiopol0_m%i_v%i", m, v), "pol0", -1, 1);
	      BulkRatioPol0[m][v]->SetLineColor(kGreen-3);
	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig]->Fit(BulkRatioPol0[m][v], "R+");
	      cout <<  BulkRatioPol0[m][v]->GetParameter(0) << " +- " << BulkRatioPol0[m][v]->GetParError(0)<< endl;
	      if (m==nummolt && v==PtV0Min) OOJlegend  ->AddEntry(hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig], "OOJ (w/o Trig) scaled to all OOJ","pl");
	      if (m==nummolt && v==PtV0Min) OOJlegend  ->AddEntry(hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig], "OOJ (w/o Trig) scaled to J+OJ","pl");
	      hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Draw("same");
	      hDeltaEtaDeltaPhi_PhiProjBulkScaled[m][v][PtTrig]   ->Draw("same");
	      OOJlegend->Draw("");
	      
	    }


	    if (LoopFile==1)	 {
	      if (numC==4)	      canvasPlotOOJDistr[m][IntisMC][1]->cd(v);
	      else 	      canvasPlotOOJDistr[m][IntisMC][1]->cd(v+1);
	      gStyle->SetOptStat(0);
	      gPad->SetLeftMargin(0.15);
	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig]->GetYaxis()->SetTitle("Default/New OOJ distribution");
	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.2);
	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 2);
	      if (type==0)	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.95, 1.05);
	      hDeltaEtaDeltaPhi_PhiProjBulkRatio[m][v][PtTrig]->Draw("");
	      tlineAtOneAllDeltaPhi->Draw("same"); 
	      if (numC==4)	      canvasPlotOOJDistr[m][IntisMC][1]->cd(v+4);
	      else 	      canvasPlotOOJDistr[m][IntisMC][1]->cd(v+7);
	      hDeltaEtaDeltaPhi_PhiProjJetNDRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(-10, 10);
	      hDeltaEtaDeltaPhi_PhiProjJetNRebinDRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(-10, 10);
	      hDeltaEtaDeltaPhi_PhiProjJetNDRatio[m][v][PtTrig]->Draw("e");
	      hDeltaEtaDeltaPhi_PhiProjJetNRebinDRatio[m][v][PtTrig]->Draw("e same");
	      tlineAtOneAllDeltaPhi->Draw("same"); 

	    }

	    //multiplicity comparison
	    if (LoopFile==1){
	      if (numC==4)	      canvasPlotOOJDistrMultComp[IntisMC]->cd(v);
	      else 	      canvasPlotOOJDistrMultComp[IntisMC]->cd(v+1);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->Scale(1./IntegralFactorBulkNew[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->SetMarkerColor(ColorPt[m]);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->SetLineColor(ColorPt[m]);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->Rebin(2);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,0.06);
	      if (v==PtV0Min)	      legendmult->AddEntry(hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig], Smolt[m], "pl");
	      if (m==nummolt)	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->Draw("");
	      else 	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->Draw("same");
	      if (m==0) legendmult->Draw("");

	      if (numC==4)	      canvasPlotOOJDistrMultComp[IntisMC]->cd(v+4);
	      else 	      canvasPlotOOJDistrMultComp[IntisMC]->cd(v+7);
	      gPad->SetLeftMargin(0.15);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig] = (TH1F*) 	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[m][v][PtTrig]->Clone(namePhiProjBulk[m][v]+ "_IntegralScaledRatio");
	      if (m!=nummolt) 	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig]->Divide(	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaled[nummolt][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.5,1.5);
	      if (type==0)	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.95,1.05);
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig]->GetYaxis()->SetTitle("OOJ_{mult} / OOJ_{0-100%}");
	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.2);
	      if (m==nummolt)	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig]->Draw("");
	      else 	      hDeltaEtaDeltaPhi_PhiProjBulkIntegralScaledRatio[m][v][PtTrig]->Draw("same");
	      tlineAtOneAllDeltaPhi->Draw("same"); 
	      if (m==0 && v==0) legendmult->Draw("");
	      //	      if (m==0) legendmult->Draw("");

	    }


	    canvasPlotProjRatioJet[m][IntisMC]->cd(v+1);
	    cout << "cd chosen " << endl;
	    gPad->SetLeftMargin(0.2);

	    hDeltaEtaDeltaPhi_PhiProjJetRatio[m][v][PtTrig]            =(TH1F*)	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]           ->Clone(namePhiProjJet[m][v]+"_Ratio");
	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]        =(TH1F*)	  hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]       ->Clone(namePhiProjJetZYAM[m][v]+"_Ratio");
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig] =(TH1F*)	  hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]->Clone(namePhiProjJetFromBulkFit[m][v]+"_Ratio");

	    hDeltaEtaDeltaPhi_PhiProjJetRatio[m][v][PtTrig] ->Divide( hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig] ->Divide( hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig] ->Divide( hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]);

	    for (Int_t b = 1; b<=   hDeltaEtaDeltaPhi_PhiProjJetRatio[m][v][PtTrig] ->GetNbinsX(); b++){
	      hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]->SetBinError(b,TMath::Abs(hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetBinError(b) - hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]->GetBinError(b)));
	      hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig]->SetBinError(b,TMath::Abs(hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetBinError(b) - hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]->GetBinError(b)));
	    }

	    hDeltaEtaDeltaPhi_PhiProjJetRatio[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(0.5, 1.5);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(0, 3);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(0, 3);

	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]    ->GetXaxis()->SetRangeUser(-1.0, 1.0);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig]    ->GetXaxis()->SetRangeUser(-1.0, 1.0);

	    hDeltaEtaDeltaPhi_PhiProjJetRatio[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.04);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.04);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.04);

	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]          -> SetMarkerColor(628);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig]   -> SetMarkerColor(kGreen+3);
	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]          -> SetLineColor(628);
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig]   -> SetLineColor(kGreen+3);

	    //	  hDeltaEtaDeltaPhi_PhiProjJetRatio[m][v][PtTrig]    ->Draw("");
	    hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[m][v][PtTrig]    ->Draw("same");
	    hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[m][v][PtTrig]    ->Draw("same");
	    tlineAtOne->Draw(""); 

	    if (IntisMC==1){
	      //	  canvasPlotProj[m][IntisMC]->cd(v+8);
	      if (hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig] && hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]&& hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]){
		/*
		  hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->Draw("");
		  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->Draw("same");
		  hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->Draw("same");
		*/
		cout << "here an average value of the relative error of the jet deltaPhi projection in the range ~[-1; 1]" << endl;
		Int_t Counter=0;
		Float_t Average=0;
		for (Int_t b=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()    ->FindBin(-1); b<=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()->FindBin(1); b++){
		  Counter++;
		  Average+=  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->GetBinError(b);
		}
		Average = Average/Counter;
		cout << "average " << Average << endl;
	      }
	    }

	    for (Int_t t=0; t<3; t++){
	      canvasPlotProjEta[m][t][IntisMC]->cd(v+1);
	      gPad->SetLeftMargin(0.15);
	      if (t==0){
		hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
		hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Draw("");
	      }
	      else if (t==1){
		hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
		hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Draw("");
	      }
	      else if (t==2){
		hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
		hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->Draw("");
	      }

	    }
	    HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinContent(v+1, Gaussian[m][v][PtTrig]->GetParameter(2));
	    HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinError(v+1, Gaussian[m][v][PtTrig]->GetParError(2));
	    HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinContent(v+1, GaussianAS[m][v][PtTrig]->GetParameter(2));
	    HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinError(v+1, GaussianAS[m][v][PtTrig]->GetParError(2));
	    HistoWidthGaussianEta[m][IntisMC][PtTrig]->SetBinContent(v+1, GaussianEta[m][v][PtTrig]->GetParameter(2));
	    HistoWidthGaussianEta[m][IntisMC][PtTrig]->SetBinError(v+1, GaussianEta[m][v][PtTrig]->GetParError(2));

	  } //end of IntisMC loop
	} //end of PtTrigMin loop
      } //end of v loop


      //I draw delta-phi projection done summing over pT assoc bins*****
      for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
	if (!IntisMC)      canvasSummedPt->cd(1+m);
	else       canvasSummedPt->cd(1+m+nummolt+1);
	for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	  if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	  for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
	    if (v>= numPtV0MaxOOJ && LoopFile==1) continue;
	    cout <<" m " << m <<  " v " << v << " "  <<  	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetMaximum() << endl;
	  }
	  cout <<" m " << m <<  " " <<	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetMaximum() << endl;
	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Scale(2.28);
	  hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->Scale(2.28);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->Scale(2.28);
	  cout <<" m " << m <<  " " <<	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetMaximum() << endl;
	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetYaxis()->SetRangeUser(-0.004, 0.007);//0.5 for hh
	  hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->GetYaxis()->SetRangeUser(-0.004, 0.007);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(-0.004, 0.007);

	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetYaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->GetYaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Draw("");
	  hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->Draw("same");
	}
      }

      //**************************************************************
      for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
	for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	  if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;

	  HistoWidthGaussian[m][IntisMC][PtTrig]->GetXaxis()->SetTitle("p_{T, assoc}");
	  HistoWidthGaussian[m][IntisMC][PtTrig]->GetYaxis()->SetTitle("#sigma (rad)");
	  HistoWidthGaussian[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0,0.7);
	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetLineColor(ColorWidth[IntisMC]);
	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetMarkerColor(ColorWidth[IntisMC]);
	  if (m==0)	legendDataMC->AddEntry(	HistoWidthGaussian[m][IntisMC][PtTrig],MCOrNot[IntisMC], "pl"); 

	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetXaxis()->SetTitle("p_{T, assoc}");
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetYaxis()->SetTitle("#sigma (rad)");
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0,5);
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetLineColor(ColorWidth[IntisMC]);
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetMarkerColor(ColorWidth[IntisMC]);

	  canvasWidthGaussian[0]->cd(m+1);
	  HistoWidthGaussian[m][IntisMC][PtTrig]->Draw("same");
	  if (IntisMC==1)	legendDataMC->Draw("same");

	  canvasWidthGaussianAS[0]->cd(m+1);
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->Draw("same");
	  if (IntisMC==1)	legendDataMC->Draw("same");

	  canvasWidthGaussianEta[0]->cd(m+1);
	  HistoWidthGaussianEta[m][IntisMC][PtTrig]->Draw("same");
	  if (IntisMC==1)	legendDataMC->Draw("same");
	} 
      } //end second intisMC loop

      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	PtTrigMin=PtTrig+PtTrigChosen;
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	RatioHistoWidthGaussian[m][PtTrig]=(TH1F*)	HistoWidthGaussian[m][0][PtTrig]->Clone(Form("RatioDataPythiaWidthm%i_PtTrig%.0f", m, PtTrigMin));
	RatioHistoWidthGaussian[m][PtTrig]->SetTitle(Form("RatioDataPythiaWidthm%i_PtTrig%.0f", m, PtTrigMin));
	RatioHistoWidthGaussian[m][PtTrig]->Sumw2();
	RatioHistoWidthGaussian[m][PtTrig]->Divide(HistoWidthGaussian[m][1][PtTrig]);
	RatioHistoWidthGaussian[m][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);

	RatioHistoWidthGaussianAS[m][PtTrig]=(TH1F*)	HistoWidthGaussianAS[m][0][PtTrig]->Clone(Form("RatioDataPythiaWidthASm%i_PtTrig%.0f", m, PtTrigMin));
	RatioHistoWidthGaussianAS[m][PtTrig]->SetTitle(Form("RatioDataPythiaWidthASm%i_PtTrig%.0f", m, PtTrigMin));
	RatioHistoWidthGaussianAS[m][PtTrig]->Sumw2();
	RatioHistoWidthGaussianAS[m][PtTrig]->Divide(HistoWidthGaussianAS[m][1][PtTrig]);
	RatioHistoWidthGaussianAS[m][PtTrig]->GetYaxis()->SetRangeUser(0,2);

	canvasWidthGaussian[1]->cd(m+1);
	RatioHistoWidthGaussian[m][PtTrig]->Draw("same");

	canvasWidthGaussianAS[1]->cd(m+1);
	RatioHistoWidthGaussianAS[m][PtTrig]->Draw("same");	  
      }

      //    if (m==0 || m>3){
   

      //    }

      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	/*
	  fileout->WriteTObject(RatioHistoWidthGaussian[m][PtTrig]);
	  fileout->WriteTObject(HistoWidthGaussian[m][0][PtTrig]);
	  fileout->WriteTObject(HistoWidthGaussian[m][1][PtTrig]);
	  fileout->WriteTObject(RatioHistoWidthGaussianAS[m][PtTrig]);
	  fileout->WriteTObject(HistoWidthGaussianAS[m][0][PtTrig]);
	  fileout->WriteTObject(HistoWidthGaussianAS[m][1][PtTrig]);
	  fileout->WriteTObject(HistoWidthGaussianEta[m][0][PtTrig]);
	  fileout->WriteTObject(HistoWidthGaussianEta[m][1][PtTrig]);
	*/
      }
    } //end of mult loop
  } //end loop on input file for OOJ comparison

  for(Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
    fileout->WriteTObject(canvasPlotOOJDistrMultComp[IntisMC]); 
  }
  for(Int_t m=nummolt; m>=0; m--){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && (m==0 || m==1)) continue;
    for(Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
      fileout->WriteTObject(canvasPlot[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotOOJDistr[m][IntisMC][0]); 
      fileout->WriteTObject(canvasPlotOOJDistr[m][IntisMC][1]); 
      fileout->WriteTObject(canvasPlotME[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotProj[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotProjBis[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotProjAllPt[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotProjRatioJet[m][IntisMC]); 

      if (m==nummolt && IntisMC == LimInfMC) {
	canvasPlotProj[m][IntisMC]->SaveAs(nomepdffile+"_PlotProj.pdf(");
	canvasPlotProjBis[m][IntisMC]->SaveAs(nomepdffile+"_PlotProjBis.pdf(");
	canvasPlotProjAllPt[m][IntisMC]->SaveAs(nomepdffile+"_PlotProjAllPt.pdf(");
      }
      else if ((m==0 || (isppHM && m==2)) && IntisMC == LimSupMC) {
	canvasPlotProj[m][IntisMC]->SaveAs(nomepdffile+"_PlotProj.pdf)");  
	canvasPlotProjBis[m][IntisMC]->SaveAs(nomepdffile+"_PlotProjBis.pdf)");  
	canvasPlotProjAllPt[m][IntisMC]->SaveAs(nomepdffile+"_PlotProjAllPt.pdf)");  
      }
      else {
	canvasPlotProj[m][IntisMC]->SaveAs(nomepdffile+"_PlotProj.pdf"); 
	canvasPlotProjBis[m][IntisMC]->SaveAs(nomepdffile+"_PlotProjBis.pdf"); 
	canvasPlotProjAllPt[m][IntisMC]->SaveAs(nomepdffile+"_PlotProjAllPt.pdf"); 
      }

      canvasPlot[m][IntisMC]->Close(); 
      canvasPlotME[m][IntisMC]->Close(); 
      canvasPlotProj[m][IntisMC]->Close(); 
      canvasPlotProjBis[m][IntisMC]->Close(); 
      canvasPlotProjRatioJet[m][IntisMC]->Close(); 

      for (Int_t t=0; t<3; t++){
	fileout->WriteTObject(canvasPlotProjEta[m][t][IntisMC]); 

      }
    }
  }

  for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){ 
    fileout->WriteTObject(canvasPlotMEEtaProj[IntisMC]); 
    fileout->WriteTObject(canvasPlotMEPhiProj[IntisMC]); 
  }
  fileout->WriteTObject(canvasSummedPt); 
  fileout->WriteTObject(canvasWidthGaussian[0]);
  fileout->WriteTObject(canvasWidthGaussianAS[0]);
  fileout->WriteTObject(canvasWidthGaussian[1]);
  fileout->WriteTObject(canvasWidthGaussianAS[1]);
  fileout->WriteTObject(canvasWidthGaussianEta[0]);
  fileout->WriteTObject(canvasWidthGaussianEta[1]);


  cout << " end loop on two files " << endl;

  fileout->Close();

  cout << "baseline fits " << endl;
  for(Int_t m=nummolt; m>=0; m--){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && (m==0 || m==1)) continue;
    //      if (m==0) continue;
    cout << "\n\n" << endl;
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	cout << " m " << m << " v " << v << " PtTrig " << PtTrig  << endl;
	cout << "baseline value" <<  Baseline[m][v][PtTrig]->GetParameter(0)<< " +- " <<  Baseline[m][v][PtTrig]->GetParError(0)<< " number of sigmas from zero " <<  TMath::Abs(Baseline[m][v][PtTrig]->GetParameter(0)/ Baseline[m][v][PtTrig]->GetParError(0)) << endl;
      }
    }
  }

  for(Int_t m=nummolt; m>=0; m--){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && (m==0 || m==1)) continue;
    //      if (m==0) continue;
    cout << "\n\n" << endl;
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
 
	if (hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig] && hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]&& hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]){
	  cout << " m " << m << " v " << v << " PtTrig " << PtTrig  << endl;
	  //	  cout << "here an average value of the relative error of the jet deltaPhi projection in the range ~[-1; 1]" << endl;
	  Int_t Counter=0;
	  Float_t Average=0;
	  for (Int_t b=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()    ->FindBin(-1); b<=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()->FindBin(1); b++){
	    Counter++;
	    Average+= TMath::Abs( hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->GetBinContent(b));
	  }
	  Average = Average/Counter;
	  cout << "average " << Average << endl;
	}
      }
    }
  }

  for(Int_t m=nummolt; m>=0; m--){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && (m==0 || m==1)) continue;
    //      if (m==0) continue;
    cout << "\n\n" << endl;
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	cout << " v " << endl;
	//	cout << "JET+BULK normalization low limit " << hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->FindBin(1) << " upper limit " <<   hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->FindBin(2)<< " integral: " << 	 hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Integral( hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->FindBin(1),  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->FindBin(2))<< endl;
	//	cout << "BULK normalization low limit " << hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->FindBin(1) << " upper limit " <<   hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->FindBin(2) << " integral: " << 	hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Integral( hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->FindBin(1),  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->FindBin(2)) << endl;
	cout <<   ScalingFactorJetNotBulkSub[m][v] << " " <<    ScalingFactorBulk[m][v] << endl;
      }
    }
  }

  cout << "pol0 fit to default/new OOJ distribution ratio " << endl;
  for(Int_t m=nummolt; m>=0; m--){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && (m==0 || m==1)) continue;
    cout << "\nm " << m << endl;
    for(Int_t v=PtV0Min; v<numPtV0MaxOOJ; v++){  
      cout << " v " << v <<"\n" << BulkRatioPol0[m][v]->GetParameter(0) << " +- " << BulkRatioPol0[m][v]->GetParError(0)<< endl;
      cout << "scaling factor M2 " << ScalingFactorJetNotBulkSub[m][v]/ScalingFactorBulk[nummolt][v] << endl;
      cout << "scaling factor M1 " << IntegralFactorBulkDef[m][v][0]/IntegralFactorBulkNew[nummolt][v][0] << endl;
      cout << "scaling factor M1 (mult corr) " << IntegralFactorBulkDef[m][v][0]/IntegralFactorBulkNew[m][v][0] << endl;
      //      cout << "pol0 to JetRebSmooth in 1 < dphi < 4 " <<     pol0JetRebSmooth[m][v]->GetParameter(0) <<" +- " << pol0JetRebSmooth[m][v]->GetParError(0)<< endl;
      //      cout << "pol0 to JetRebSmoothBis in 1 < dphi < 4 " <<     pol0JetRebSmoothBis[m][v]->GetParameter(0)<< " +- " << pol0JetRebSmoothBis[m][v]->GetParError(0)<< endl;
    }
  }

  cout << "DeltaEtaJet (default selection, or sys==0) " << DeltaEtaJetDef << endl;
  cout << "DeltaEtaJet " << DeltaEtaJet << "DeltaEtaInclusive " << DeltaEtaInclusive << endl;
  cout << "\n\n ho creato il file " << nomefileoutput << endl;
  cout << " ho creato le canvas " << nomepdffile +"_PlotProjAllPt.pdf e " << nomepdffile+"_PlotProj.pdf"<< endl; 
}

