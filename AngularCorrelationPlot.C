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

void AngularCorrelationPlot(Bool_t isTrigEff =0, Bool_t isTrigEffComp=0, Bool_t ishhCorr=0, Float_t PtTrigChosen=3, Float_t PtTrigMax =15, Bool_t SkipAssoc=0,Int_t israp=0, Int_t sysV0=0, Int_t sysTrigger=0,Int_t sys=0, Int_t type=0, Int_t PtIntervalShown=1,   TString year0 = "2016",TString year="16kl_hK0s"/*_Hybrid"/*"2019h11_HM_hK0s"/*"161718_hXi"/*"1617_GP_hK0s"/"17pq_hK0s"/"161718_HM_hXi_WithFlat16k_No18p"/*"161718_HM_hXi"/"161718Full_AOD234_hXi"/*"17pq_pp5TeV_hXi_pttrig0.15"/*"17pq_hXi"/*"17pq_pp5TeV_Hybrid"/"1617_AOD234_hK0s"/*"17pq_hXi"/*"LHC16kl_pass2_GP_Fio"/*"1617GP_hK0s"/*"1617_AOD234_hK0s"/*"161718_hXi"/*"161718_MD_hXi_New"/*"2018f1_extra_hK0s_Fio"/"17pq_hK0s"/*"LHC17_AOD234_Red"/"AllhK0sHM_RedNo16k"/*"2016k_HM_hK0s"/*"1617GP_hK0s_Hybrid_New"/*"1617_hK0s"/*"2018g4_extra_hXi_SelTrigger"/*"161718_MD_hXi"/"2018f1g4_extra_EtaEff_hXi"/*"2018g4_extra_EtaEff_Hybrid_hK0s"/""/*"161718_MD_EtaEff_hXi"/*"AllMC_hXi"/*"2016kehjl_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*/, TString yearMC="16kl_hK0s"/*"16kl_hK0s_Hybrid"/*"2019h11_HM_hK0s"/*"161718_hXi"/"161718_HM_hXi"/*"17pq_hXi"/*"17pq_pp5TeV_Hybrid"/"1617_AOD234_hK0s"/*"17pq_hK0s"/*"LHC16kl_pass2_GP_Fio"/*"1617GP_hK0s"/*"161718_hXi"/*"161718_MD_New_hXi"/*"2018f1_extra_hK0s_Fio"/*"1617GP_hK0s_Hybrid_New"/*"1617GP_hK0s"/*"1617MC_hK0s"/*"161718_MD_EtaEff_hXi"/"2018g4_extra_EtaEff_Hybrid_hK0s"/*"2018f1g4_extra_EtaEff_hXi"/"2018g4_extra_EtaEff_hK0s"/*"2018f1_extra_Reco_hK0s"/"1617MC_hK0s"/"AllMC_hXi"/*"161718_MD_hXi_Hybrid"*/,  TString Path1 =""/*"_Jet0.75"/*"_PtTrigMax2.5_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, TString Path2 =""/*"_PtTrigMax2.5"/*"_NewMultClassBis_Jet0.75"*/, TString Dir ="FinalOutput/", Bool_t isEnlargedDeltaEta=0, Int_t isMC=1, Int_t isEfficiency=1, Int_t MultBinning=0, Int_t PtBinning=1, Bool_t isSidebands=0, Bool_t isSidebandsAnalysis =0, Bool_t IsMEFromHybrid=0, Bool_t isMEFromCorrectCentrality=0, Bool_t isCompWithMEFromHybrid=0, TString yearHybrid =""/*"2018g4_extra_EtaEff_Hybrid_hK0s"/* "161718_MD_hXi_Hybrid_MCTruth"*/, Bool_t IsParticleTrue=0, Bool_t isEtaEff=1, Bool_t isEtaEffComp=0, Bool_t isEta05=0 , Bool_t isCompWithMEFromXi=0, TString yearXiComp= "161718_MD_EtaEff_hXi",  Int_t TriggerPDGChoice=0, Bool_t isNewInputPath=1, Bool_t isHM=0, Bool_t isEPOSEff=0, Bool_t isMEFromPeak=0, Bool_t isMEFromK0s=0, Bool_t isBkgParab=0 , Bool_t isSysDef =0, Int_t isDefaultSel = 0, Bool_t isMCForNorm=0, Bool_t VarRange=0){

  Int_t SidebandsSide =0;
  if (isSidebands){
  cout << "Do you want to use the left + right sidebands (0), the left one only (1), or the right one only (2)?" << endl;
  cin >> SidebandsSide;
  }
  Bool_t isNewDEtaJet=1;

  if (isMC && !isEfficiency) isEtaEff=0;
  Bool_t isMEFrom13TeV = 0;
  if ((isDefaultSel==0 || isDefaultSel > 3) && isSysDef) return;
  //isDefaultSel=3 -> isLoosest, isDefaultSel=2 -> isTightest

  if (isMEFromK0s && type==0) {cout << "the option isMEFromK0s is meant to be used when Xi is being analyzed" <<endl; return;}
  if (year=="17pq_hXi") MultBinning=3;
  if (year=="17pq_pp5TeV_hXi_pttrig0.15") MultBinning=3;
  //isSidebands = 1 :to display also the ME and the SE/ME obtained from sdebands of invariant mass (to check for example if there is enough statistics to use the sidebands)
  //isSidebandsAnalysis = 1 :to take dphi projections obtained by subtracting fake K0s/Xi contribution using the sidbands distribution

  //isEtaEffComp=1: comparison between ME obtained with pt-dependent eff and ME obtained with pt and eta-dependent efficiency is done
  if (isEtaEffComp==1 && isEtaEff==0) isEtaEffComp=0;

  //isCompWithMEFromXi=1: comparison between K0s and Xi deltaEta proj of ME (only available when running macro for K0s) 
  if (isCompWithMEFromXi && type==8) {cout << " comparison between Xi and K0s ME is done only when running the macro for K0s and we need to use the same pt binning! " << endl; return; }

  TString sTriggerPDGChoice[3] = {"", "_IsOnlypiKpemu", "_IsNotSigmaOnly"};
  if (!isMC) TriggerPDGChoice=0;

  Bool_t IsSpecial=0;
  if (year == "1617_hK0s" || yearMC=="1617MC_hK0s") IsSpecial=1;
  if (isMC==2 && year == "1617_hK0s" && yearMC!="1617MC_hK0s") return;
  if (isMC ==2 && year != "1617_hK0s" && yearMC=="1617MC_hK0s") return;

  if (type==0) {
    //    year = "1617_hK0s";
    //    PtBinning=1;
  }
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
  gStyle->SetOptFit(1111);

  //lista degli effetti  sistematici studiati in questa macro
  if (sys==3 || sys>8) return;
  if (sysV0>6)return;
  if (ishhCorr && sys<4 && sys!=0) return; //se faccio correlazione hh non posso studiare sistematico associato a scelta regione Sidebands e peak nella distribuzione di massa invariante
  if (sysV0>2 && ishhCorr) return;
  if (sysTrigger!=0) return;
  if (sysV0!=0 && sys!=0) return;

  Double_t LimSupY[2] = {0.0006, 0.05}; 
  Double_t LimInfY[2] = {-0.0001, -0.004}; 
  TString hhCorr[2]= {"", "_hhCorr"};
  Int_t sysang=0;
  if (sys==1) sysang=1;
  if (sys==2) sysang=2;
  if (sys==4) sysang =4;
  if (sys==7) sysang =7;
  if (sys==8) sysang =8;
  //sysang =4; //new

  TLegend *legendPt = new TLegend(0.6, 0.7, 0.9, 0.9);
  legendPt->SetHeader("p_T^{Assoc} intervals");
  TLegend *legendMult = new TLegend(0.6, 0.7, 0.9, 0.9);
  legendMult->SetHeader("Multiplicity classes");
  //  Int_t Colormult[21]={1, 401, 801, 628, 909, 881, 860, 868, 841, 418, 628, 909, 881, 867, 921, 401,  841, 862, 866, 865, 864};
  Int_t Colormult[]={1, 801, 628, 867, 909,  881, 860, 868, 841, 418, 628, 909, 881, 867, 921, 401,  841, 862, 866, 865, 864};

  Dir+="DATA"+year0;
  TString file[2];
  file[0] = year;
  if (PtBinning>0) file[0]+= Form("_PtBinning%i", PtBinning);
  file[0]+= Path1;
  file[1] = yearMC;
  if (isEfficiency) file[1]+= "_MCEff" ;
  else file[1]+= "_MCTruth" ;
  file[1]+= Path1;
  if (PtBinning>0) file[1]+= Form("_PtBinning%i", PtBinning);
  file[1]+=Path2; 

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9;
  Int_t numPtV0Max =numPtV0;
  if (PtBinning==0)numPtV0Max =numPtV0-1;
  const Int_t numPtTrig=10;
  const Int_t numtipo=10;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
 
  Int_t PtV0Min=1;
  if (!ishhCorr && type==0) PtV0Min=0;
  if (isMC && !isEfficiency && !isMCForNorm) PtV0Min =0;
  Int_t numPtV0Chosen = 4; 
  Int_t multchosen =5;

  Int_t ColorPt[numPtV0]= {401,801,628,909,881,860,868,842, 921};

  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  //  TString SPtV0[numPtV0]={"", "", "0.5-1", "1-1.5","1.5-2","2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  if (isMC && !isEfficiency && !isMCForNorm) SPtV0[0]={"0-0.5"};
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
  if (isEta05) Srap[0] = "_Eta0.5";
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  TLegend * legendDataMC = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend * legendSB = new TLegend(0.5, 0.7, 0.9, 0.9);
  Int_t ColorWidth[2]={868, 418};
  TString MCOrNot[2]={" Data" , " Pythia8" };
  TString SmoltBis[nummolt+1]={"0-5%", "5-10%", "10-30%", "30-50%", "50-100%", "0-100%"};
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  TString SmoltBispp5TeV[nummolt+1]={"0-10%", "10-100%", "100-100", "100-100", "100-100", "0-100%"};
  Double_t Nmoltpp5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};
  TString Smolt[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,2,7,15,30,100}; 

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Nmolt[m] = Nmolt0[m];
      Smolt[m] = Smolt0[m];
    }
    if (MultBinning==1){
      Nmolt[m] = Nmolt1[m];
      Smolt[m] = Smolt1[m];
    }
    if (MultBinning==2){
      Nmolt[m] = Nmolt2[m];
      Smolt[m] = Smolt2[m];
    }
    else if (MultBinning==3){
      Smolt[m] = Smoltpp5TeV[m];
      SmoltBis[m] = SmoltBispp5TeV[m];
      Nmolt[m] = Nmoltpp5TeV[m];
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
    SmoltBis[0] = "0-0.001 %";
    SmoltBis[1] = "0.001-0.005 %";
    SmoltBis[2] = "0.005-0.01 %";
    SmoltBis[3] = "0.01-0.05 %";
    SmoltBis[4] = "0.05-0.1 %";
    SmoltBis[5] = "0-0.1 %";
    if (MultBinning == 1) {
      Nmolt[1] = 0;
      Nmolt[2] = 0;
      Smolt[0] = "0-0a";
      Smolt[1] = "0-0b";
      Smolt[2] = "0-0.01";
      SmoltBis[0] = "0-0a %";
      SmoltBis[1] = "0-0b %";
      SmoltBis[2] = "0-0.01 %";
    }
  }

  TString Ssideband[numSB]={"", "_SB"};
  TString Ssidebandtitle[numSB]={"", " (sidebands)"};
  TString Szeta[numzeta]={""};
  TString PtTitle = "";

  TString nameisEtaEff="";
  if (isEtaEff) nameisEtaEff="_Effw";
  TString nameRawSE[nummolt+1][numPtV0];
  TString nameSE[nummolt+1][numPtV0];
  TString nameME[nummolt+1][numPtV0];
  TString nameSESB[nummolt+1][numPtV0];
  TString nameMESB[nummolt+1][numPtV0];
  TString nameMEEtaProj[nummolt+1][numPtV0];
  TString nameMEPhiProj[nummolt+1][numPtV0];
  TString      namePhiProjFakeJet[nummolt+1][numPtV0];
  TString      namePhiProjFakeBulk[nummolt+1][numPtV0];
  TString      namePhiProjFakeAll[nummolt+1][numPtV0];
  TString      namePhiProjSBJet[nummolt+1][numPtV0];
  TString      namePhiProjSBBulk[nummolt+1][numPtV0];
  TString      namePhiProjSBAll[nummolt+1][numPtV0];

  TString namePhiProjJet[nummolt+1][numPtV0];
  TString namePhiProjJetZYAM[nummolt+1][numPtV0];
  TString namePhiProjJetFromBulkFit[nummolt+1][numPtV0];
  TString namePhiProjJetNotBulkSub[nummolt+1][numPtV0];
  TString namePhiProjBulk[nummolt+1][numPtV0];
  TString namePhiProjJetBulk[nummolt+1][numPtV0];
  TString nameEtaProj[nummolt+1][numPtV0][numPtTrig];
  TString nameRawEtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F * hMEEntriesvsMult= new TH1F ("hMEEntriesvsMult", "hMEEntriesvsMult", nummolt, Nmolt);
  TH1F * hSEEntriesvsMult= new TH1F ("hSEEntriesvsMult", "hSEEntriesvsMult", nummolt, Nmolt);
  TH1F * hMEtoSEEntriesvsMult;
  TH1F * hMEEntriesvsPt[nummolt+1];
  TH1F * hSEEntriesvsPt[nummolt+1];
  TH1F * hMEtoSEEntriesvsPt[nummolt+1];
  Int_t MEEntries[nummolt]={0};
  Int_t SEEntries[nummolt]={0};
  TH2F *hDeltaEtaDeltaPhi_RawSEbins[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_MEbinsHybrid[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_MEbinsXi[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_MEbinsComp[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_SEbinsSB[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_MEbinsSB[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjXi[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjComp[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjCompRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[nummolt+1][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjMasterMolt[numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[nummolt+1][numPtTrig];

  TH1F *hDeltaEtaDeltaPhi_PhiProjJetSB[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetFakeSB[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkSB[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjAllSB[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjAllFakeSB[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[nummolt+1][numPtV0][numPtTrig];

  TF1 * polJetSB[nummolt+1][numPtV0][numPtTrig];
  TF1 * polBulkSB[nummolt+1][numPtV0][numPtTrig];
  TF1 * polAllSB[nummolt+1][numPtV0][numPtTrig];
  TF1 * pol0Fit[nummolt+1][numPtV0][numPtTrig];

  TH1F *hDeltaEtaDeltaPhi_PhiProjJet[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetZYAM[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetFromBulkFitRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetZYAMRatio[nummolt+1][numPtV0][numPtTrig];

  TF1*	  pol0ZYAM[nummolt+1][numPtV0][numPtTrig];
  TF1*	  pol0ZYAMOrigin[nummolt+1][numPtV0][numPtTrig];
  TF1*	  pol0BulkBis[nummolt+1][numPtV0][numPtTrig];

  TH1F *hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_RawSEEtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_RawSEEtaProjPHalfRatioToME[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProjPHalf[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProjNHalf[nummolt+1][numPtV0][numPtTrig];

  TF1* pol0Centr[nummolt+1][numPtV0][numPtTrig];
  TF1 * pol0SideDX[nummolt+1][numPtV0][numPtTrig];
  TF1 * pol0SideSX[nummolt+1][numPtV0][numPtTrig];

  TString PathIn[numPtTrig];
  TString PathInComp[numPtTrig];
  TString PathInMEHybrid[numPtTrig];
  TString PathInMEXi;

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
  TLine *tlinePhiSx=new TLine(-1.15, -0.85, 1.15,  -0.85); //was 1.32
  TLine *tlinePhiDx=new TLine(-1.15, 0.85, 1.15,  0.85);
  TLine *tlinePhiBase=new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2,0);
  TLine *tlineAtOne=new TLine(-1, 1, 1, 1);
  TLine *tlineAtOneLong=new TLine(-1.2, 1, 1.2, 1);
  TLine *tlineAtOnePt=new TLine(0, 1, 8, 1);
  tlineAtOnePt->SetLineColor(881);
  TLine*   lineEta075;
  TLine*   lineEtam075;

  TCanvas*  canvasMEtoSEMult = new TCanvas("canvasMEtoSEMult", "canvasMEtoSEMult", 1300, 800);
  TCanvas*  canvasMEtoSE = new TCanvas("canvasMEtoSE", "canvasMEtoSE", 1300, 800);
  canvasMEtoSE->Divide(3,1);
  TCanvas *  canvasWidthGaussian[2];
  TCanvas *  canvasWidthGaussianEta[2];
  TCanvas *  canvasWidthGaussianAS[2];
  TCanvas * canvasWings[2];

  TCanvas*    canvasSBJetFitResult = new TCanvas("canvasSBJetFitResult", "canvasSBJetFitResult", 1300,800);
  canvasSBJetFitResult->Divide(3,3);
  TCanvas*    canvasSBBulkFitResult = new TCanvas("canvasSBBulkFitResult", "canvasSBBulkFitResult", 1300,800);
  canvasSBBulkFitResult->Divide(3,3);
  TCanvas*    canvasSBAllFitResult = new TCanvas("canvasSBAllFitResult", "canvasSBAllFitResult", 1300,800);
  canvasSBAllFitResult->Divide(3,3);

  for (Int_t i =0; i<2; i++){  
    canvasWings[i] = new TCanvas (Form("canvasWings%i",i), Form("canvasWings%i",i), 800, 500);
    canvasWings[i]->Divide(3,2);

    canvasWidthGaussian[i] = new TCanvas (Form("canvasWidthGaussian%i",i), Form("canvasWidthGaussian%i",i), 800, 500);
    canvasWidthGaussian[i]->Divide(3,2);
    canvasWidthGaussianEta[i] = new TCanvas (Form("canvasWidthGaussianEta%i",i), Form("canvasWidthGaussianEta%i",i), 800, 500);
    canvasWidthGaussianEta[i]->Divide(3,2);
    canvasWidthGaussianAS[i] = new TCanvas (Form("canvasWidthGaussianAS%i",i), Form("canvasWidthGaussianAS%i",i), 800, 500);
    canvasWidthGaussianAS[i]->Divide(3,2);
  }

  TCanvas *canvasPlot[nummolt+1][2];
  TCanvas *canvasPlotSB[nummolt+1][2];
  TCanvas *canvasPlotSE[nummolt+1][2];
  TCanvas *canvasPlotME[nummolt+1][2];
  TCanvas *canvasPlotMEEtaProjComp[nummolt+1][2];
  TCanvas *canvasPlotMEEtaProjCompRatio[nummolt+1][2];
  TCanvas *canvasPlotMESB[nummolt+1][2];
  TCanvas *canvasPlotMEEtaProj[2];
  TCanvas *canvasPlotMEEtaProjCompHybrid[nummolt+1][2];
  TCanvas *canvasPlotMEEtaProjCompHybridRatio[nummolt+1][2];
  TCanvas *canvasPlotMEEtaProjCompXi[nummolt+1][2];
  TCanvas *canvasPlotMEEtaProjCompXiRatio[nummolt+1][2];
  TCanvas *canvasPlotSEEtaProj[2];
  TCanvas *canvasPlotMEEtaProjMolt[2];
  TCanvas *canvasPlotMEPhiProj[2];
  TCanvas *canvasPlotProj[nummolt+1][2];
  TCanvas *canvasPlotProjSB[nummolt+1][2][3];
  TCanvas *canvasPlotProjRatioSB[nummolt+1][2][3];
  TCanvas *canvasPlotProjRatioJet[nummolt+1][2];
  TCanvas *canvasPlotProjEta[nummolt+1][3][2];
  TCanvas *canvasPlotRawSEProjEta[nummolt+1][4][2];

  TH1F *HistoWidthGaussian[nummolt+1][2][numPtTrig];
  TH1F *HistoWidthGaussianAS[nummolt+1][2][numPtTrig];
  TH1F *HistoWidthGaussianEta[nummolt+1][2][numPtTrig];
  TH1F *RatioHistoWidthGaussian[nummolt+1][numPtTrig];
  TH1F *RatioHistoWidthGaussianAS[nummolt+1][numPtTrig];
  TH1F *RatioHistoWidthGaussianEta[nummolt+1][numPtTrig];
  TH1F * hWings[nummolt+1][2][numPtTrig];
  TH1F * hWingsLeft[nummolt+1][2][numPtTrig];
  TH1F * hWingsRight[nummolt+1][2][numPtTrig];
  TH1F*  hJetSBFitResult[nummolt+1];
  TH1F*  hBulkSBFitResult[nummolt+1];
  TH1F*  hAllSBFitResult[nummolt+1];

  TString PhiRegion[3]={"All", "DeltaPhi>Pi/2", "DeltaPhi < Pi/2"};
  TString Region[3]={"Jet", "Bulk", "All"};

  //*********************************************************************                                                                                                  
  //**************calcolo numero particelle di trigger*******************                                                                                                
  Float_t NTrigger[nummolt+1][2];
  TH1F * histoNTrigger;                           
  TH1F * histoNTriggerMult;                           

  TFile *fileinbis[2]; //Data and MC
  TFile *fileinbisPart1[2]; //Data and MC
  TFile *fileinbisPart2[2]; //Data and MC
  TString PathInBis[2];
  TString PathInBisPart1[2];
  TString PathInBisPart2[2];
  PathInBis[0] =  "FinalOutput/AnalysisResults" + year +".root"; //+ Path1  + ".root"; change
  PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMC;
  if (isEfficiency)  PathInBis[1] += "_MCEff";
  else  PathInBis[1] += "_MCTruth";
  PathInBis[1] += ".root";
  if (ishhCorr){  
    PathInBis[0] =  "FinalOutput/AnalysisResults" + year +Path1+ "_hhCorr.root";
    PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMC  + Path1+"_hhCorr_MCEff.root";
  }

  if (IsSpecial){
    PathInBisPart1[0] =  "FinalOutput/AnalysisResults2016kehjl_hK0s.root"; //+ Path1  + ".root"; change
    PathInBisPart1[1] =  "FinalOutput/AnalysisResults2018f1d8_extra_hK0s_MCEff.root";// + Path1 +".root"; change

    PathInBisPart2[0] =  "FinalOutput/AnalysisResultsLHC17_hK0s.root"; //+ Path1  + ".root"; change
    PathInBisPart2[1] =  "FinalOutput/AnalysisResultsLHC17anch17_hK0s_MCEff.root"; //+ Path1  + ".root"; change


  }

  for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){ //change
    fileinbis[IntisMC]=new TFile(PathInBis[IntisMC],"");
    if (IsSpecial){
      fileinbisPart1[IntisMC]=new TFile(PathInBisPart1[IntisMC],"");
      fileinbisPart2[IntisMC]=new TFile(PathInBisPart2[IntisMC],"");
    }
    cout <<"file from task: " << fileinbis[IntisMC] << endl;
    if (!fileinbis[IntisMC]) return;

    TString ContainerName = "";
    if (isNewInputPath) {
      //      if (year.Index("Fio")!=-1) ContainerName="_hK0s_Task_suffix";
      if (year == "2016k_HM_hK0s") ContainerName="_hK0s_Task_suffix";
      else if (isMC && !isEfficiency){
	if (type==8)	ContainerName="_h"+ tipo[type] +"_Task_MCTruth";
	else ContainerName="_h"+ tipo[type] +"_Task_Truth";
	if (year.Index("Hybrid")!=-1) ContainerName="_h"+tipo[type] +"_Task_Hybrid";
	else if (year.Index("Fio")!=-1) ContainerName="_hK0s_Task_Truth";
      }
      else if (isMC && isEfficiency){
	ContainerName="_h"+ tipo[type] +"_Task_RecoAndEfficiency";
      }
      else {
	if (type==0) ContainerName="_hK0s_Task_";
	else ContainerName="_hXi_Task_Default";
	//	else ContainerName="_hXi_Task_";
      }
    }
    TString TaskName ="";
    if (isNewInputPath) {
      if (year == "2016k_HM_hK0s") TaskName = "_PtTrigMin3.0_PtTrigMax30.0";
      //      else if (year.Index("Fio")!=-1) TaskName = "_MCTruth_PtTrigMin0.2_PtTrigMax15.0";
      else if (isMC) TaskName = "_MCTruth_PtTrigMin3.0_PtTrigMax15.0";
      else if (year == "17pq_pp5TeV_hXi_pttrig0.15") TaskName = "_PtTrigMin0.2_PtTrigMax2.5";
      else TaskName = "_PtTrigMin3.0_PtTrigMax15.0";
    }

    TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Lambda", "Xi",   "Xi", "Omega","Omega", "Xi", "Omega"};
    cout << "Drectory name " << "MyTask"<<dirinputtype[type] << TaskName<< endl;
    TDirectoryFile *dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTask"+dirinputtype[type]+ TaskName);
    if (!dir) {cout << " directory not present " << endl; return;}
    TDirectoryFile *dirPart1;
    TDirectoryFile *dirPart2;
    TList *listPart1;
    TList *listPart2;
    if (IsSpecial){
      dirPart1 = (TDirectoryFile*)fileinbisPart1[IntisMC]->Get("MyTask"+dirinputtype[type]);
      dirPart2 = (TDirectoryFile*)fileinbisPart2[IntisMC]->Get("MyTask"+dirinputtype[type]);
      listPart1 = (TList*)dirPart1->Get("MyOutputContainer"+ ContainerName);
      listPart2 = (TList*)dirPart2->Get("MyOutputContainer"+ ContainerName);
    }
    TList *list = (TList*)dir->Get("MyOutputContainer"+ ContainerName);
    cout << "Container name " << ContainerName<< endl;
    if (!list) {cout << " list not there " << endl; return;}
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


    cout << "ok up to here " << endl;
    TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigChosen+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMax-0.00001) );

    //*********************************************************************                                                                     
    for(Int_t m=0; m<nummolt+1; m++){
      if (isHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      //      if (m==0) continue;
      NTrigger[m][IntisMC]=0;
      if(m<nummolt){
	for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.0001); j<=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.0001); j++ ){
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
  TFile *fileinComp[numPtTrig];
  TFile *fileinCompMEHybrid[numPtTrig];
  TFile *fileinCompMEXi;
  TString nomefileoutput= "AngularCorrelationPlot";
  if (isMC==0)  nomefileoutput+=year;
  else if (isMC==1) nomefileoutput+=yearMC;
  else nomefileoutput+=year + "_" +yearMC;
  if (isMC && !isEfficiency) nomefileoutput+="_MCTruth";
  if (PtBinning!=0)  nomefileoutput +=Form("_PtBinning%i",PtBinning) +Path1;
  if(type>=0){
    nomefileoutput +="_"+tipo[type];
    nomefileoutput +=Srap[israp];
    nomefileoutput +=SSkipAssoc[SkipAssoc];
  }
  if (isSysDef){
    if (isDefaultSel==1) 	 nomefileoutput += "_SysV0Default";
    else if (isDefaultSel==2) 	  nomefileoutput += "_SysV0Tightest";
    else if (isDefaultSel==3) 	 nomefileoutput += "_SysV0Loosest";
  }
  nomefileoutput += hhCorr[ishhCorr] +  isMCOrData[isMC] + Form("_PtTrigMin%.1f_Output", PtTrigChosen);
  if (IsParticleTrue)  nomefileoutput+= "_IsParticleTrue";
  if (IsMEFromHybrid) nomefileoutput+= "_IsMEFromHybrid";
  if (isMEFromK0s) nomefileoutput+= "_IsMEFromK0s";
  if (isMEFromPeak)  nomefileoutput+= "_IsMEFromPeak";
  if (isMEFromCorrectCentrality) nomefileoutput+= "_IsMEFromCorrectCentrality";
  if (isMEFrom13TeV) nomefileoutput += "_IsMEFrom13TeV";
  if (isEPOSEff)  nomefileoutput+= "_EPOS"; 
  if (SidebandsSide ==1) nomefileoutput+= "_LeftSB";
  else if (SidebandsSide ==2) nomefileoutput+= "_RightSB";
  if (isEtaEff) nomefileoutput+= "_IsEtaEff";
  if (isTrigEff) nomefileoutput+= "_IsTrigEff";
  if (isTrigEffComp) nomefileoutput+= "_IsTrigEffComp";
  if (isSidebandsAnalysis) nomefileoutput+= "_SB";
  nomefileoutput+= sTriggerPDGChoice[TriggerPDGChoice];
  if (MultBinning !=0) nomefileoutput+= Form("_MultBinning%i", MultBinning);
  if (isMC && !isEfficiency && !isMCForNorm) nomefileoutput += "_MCPrediction";
  if (isNewDEtaJet)  nomefileoutput += "_NewdEtaChoice";
  if (sysang!=0)   nomefileoutput += Form("_sys%i", sysang);
  //  nomefileoutput += "_Try";
  //  nomefileoutput+= "_LargerJetWidthLargerBulkWidth";
  if (VarRange!=0) nomefileoutput+= Form("_VarRange%i", VarRange);
  TString nomefileoutputpdf=nomefileoutput;
  //  nomefileoutput +="_thinptbins";
  nomefileoutput += ".root";
  //  nomefileoutputpdf += ".pdf";

  if (!ishhCorr && isEnlargedDeltaEta)  nomefileoutput= "AngularCorrelationPlot" + hhCorr[ishhCorr] + Form("_PtTrigMin%i_DeltaEtaPhiEnlarged_Output.root", PtTrigChosen);
   
  TFile *fileout = new TFile(nomefileoutput, "RECREATE");

  //*****************************************************
  Int_t counter=0;
  TCanvas*   canvasSummedPt=new TCanvas ("canvasSummedPt", "canvasSummedPt", 800, 500);
  canvasSummedPt->Divide(nummolt+1, 2);
  for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
    canvasPlotMEEtaProj[IntisMC]=new TCanvas(Form("canvasPlotMEEtaProj_%i", IntisMC), Form("canvasPlotMEEtaProj_%i", IntisMC), 1300, 800);
    canvasPlotMEEtaProj[IntisMC]->Divide(3,2);
    canvasPlotSEEtaProj[IntisMC]=new TCanvas(Form("canvasPlotSEEtaProj_%i", IntisMC), Form("canvasPlotSEEtaProj_%i", IntisMC), 1300, 800);
    canvasPlotSEEtaProj[IntisMC]->Divide(6,2);
    canvasPlotMEEtaProjMolt[IntisMC]=new TCanvas(Form("canvasPlotMEEtaProjMolt_%i", IntisMC), Form("canvasPlotMEEtaProjMolt_%i", IntisMC), 1300, 800);
    canvasPlotMEEtaProjMolt[IntisMC]->Divide(3,3);
    canvasPlotMEPhiProj[IntisMC]=new TCanvas(Form("canvasPlotMEPhiProj_%i", IntisMC), Form("canvasPlotMEPhiProj_%i", IntisMC), 1300, 800);
    canvasPlotMEPhiProj[IntisMC]->Divide(6,2);
  }

  Int_t numC =4;
  if (PtBinning>0) numC=5;
  for(Int_t m=nummolt; m>=0; m--){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //      if (m==0) continue;
    cout << "\n\n m " << m << endl;
    for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
      canvasPlot[m][IntisMC]=new TCanvas(Form("canvasPlot_m%i_MC%i",m, IntisMC), "canvasPlot"+Smolt[m], 1300, 800);
      canvasPlot[m][IntisMC]->Divide(numC,2);
      canvasPlotSB[m][IntisMC]=new TCanvas(Form("canvasPlotSB_m%i_MC%i",m, IntisMC), "canvasPlotSB"+Smolt[m], 1300, 800);
      canvasPlotSB[m][IntisMC]->Divide(numC,2);
      canvasPlotME[m][IntisMC]=new TCanvas(Form("canvasPlotME_m%i_MC%i",m,  IntisMC), "canvasPlotME"+Smolt[m], 1300, 800);
      canvasPlotME[m][IntisMC]->Divide(numC,2);
      canvasPlotSE[m][IntisMC]=new TCanvas(Form("canvasPlotSE_m%i_MC%i",m,  IntisMC), "canvasPlotSE"+Smolt[m], 1300, 800);
      canvasPlotSE[m][IntisMC]->Divide(numC,2);

      canvasPlotMEEtaProjCompHybrid[m][IntisMC]=new TCanvas(Form("canvasPlotMEEtaProjCompHybrid_m%i_MC%i",m,  IntisMC), "canvasPlotMEEtaProjCompHybrid"+Smolt[m], 1300, 800);
      canvasPlotMEEtaProjCompHybrid[m][IntisMC]->Divide(numC,2);

      canvasPlotMEEtaProjCompHybridRatio[m][IntisMC]=new TCanvas(Form("canvasPlotMEEtaProjCompHybridRatio_m%i_MC%i",m,  IntisMC), "canvasPlotMEEtaProjCompHybridRatio"+Smolt[m], 1300, 800);
      canvasPlotMEEtaProjCompHybridRatio[m][IntisMC]->Divide(numC,2);

      canvasPlotMEEtaProjCompXi[m][IntisMC]=new TCanvas(Form("canvasPlotMEEtaProjCompXi_m%i_MC%i",m,  IntisMC), "canvasPlotMEEtaProjCompXi"+Smolt[m], 1300, 800);
      canvasPlotMEEtaProjCompXi[m][IntisMC]->Divide(numC,2);

      canvasPlotMEEtaProjCompXiRatio[m][IntisMC]=new TCanvas(Form("canvasPlotMEEtaProjCompXiRatio_m%i_MC%i",m,  IntisMC), "canvasPlotMEEtaProjCompXiRatio"+Smolt[m], 1300, 800);
      canvasPlotMEEtaProjCompXiRatio[m][IntisMC]->Divide(numC,2);

      canvasPlotMEEtaProjComp[m][IntisMC]=new TCanvas(Form("canvasPlotMEEtaProjComp_m%i_MC%i",m,  IntisMC), "canvasPlotMEEtaProjComp"+Smolt[m], 1300, 800);
      canvasPlotMEEtaProjComp[m][IntisMC]->Divide(numC,2);
      canvasPlotMEEtaProjCompRatio[m][IntisMC]=new TCanvas(Form("canvasPlotMEEtaProjCompRatio_m%i_MC%i",m,  IntisMC), "canvasPlotMEEtaProjCompRatio"+Smolt[m], 1300, 800);
      canvasPlotMEEtaProjCompRatio[m][IntisMC]->Divide(numC,2);
      canvasPlotMESB[m][IntisMC]=new TCanvas(Form("canvasPlotMESB_m%i_MC%i",m,  IntisMC), "canvasPlotMESB"+Smolt[m], 1300, 800);
      canvasPlotMESB[m][IntisMC]->Divide(numC,2);

      //oldcanvasPlotProj[m][IntisMC]=new TCanvas(Form("canvasPlotProj_m%i_MC%i",m,  IntisMC), "canvasPlotProj"+Smolt[m], 1300, 800);
      canvasPlotProj[m][IntisMC]=new TCanvas(Form("canvasPlotProj_m%i_MC%i",m,  IntisMC), "canvasPlotProj"+Smolt[m],1300, 800);
      //old      canvasPlotProj[m][IntisMC]->Divide(numC,2);
      //      if (isMC && !isEfficiency)      canvasPlotProj[m][IntisMC]->Divide(numC,3);
      canvasPlotProj[m][IntisMC]->Divide(3,3);
      canvasPlotProjRatioJet[m][IntisMC]=new TCanvas(Form("canvasPlotProjRatioJet_m%i_MC%i",m,  IntisMC), "canvasPlotProjRatioJet"+Smolt[m], 1300, 800);
      canvasPlotProjRatioJet[m][IntisMC]->Divide(numC,2);


      for (Int_t t=0; t<4; t++){
	if (t<3){
	  canvasPlotProjSB[m][IntisMC][t]=new TCanvas(Form("canvasPlotProjSB_m%i_MC%i",m,  IntisMC)+Region[t], "canvasPlotProjSB"+Smolt[m]+Region[t], 1300, 800);
	  canvasPlotProjSB[m][IntisMC][t]->Divide(3,3);
	  canvasPlotProjRatioSB[m][IntisMC][t]=new TCanvas(Form("canvasPlotProjRatioSB_m%i_MC%i",m,  IntisMC)+Region[t], "canvasPlotProjRatioSB"+Smolt[m]+Region[t],  1300, 800);
	  canvasPlotProjRatioSB[m][IntisMC][t]->Divide(3,3);
	  //	canvasPlotProjEta[m][t][IntisMC]=new TCanvas(Form("canvasPlotProjEta_m%i_%i_MC%i",m,t, IntisMC), "canvasPlotProjEta"+Smolt[m]+"_"+PhiRegion[t], 1300, 800);
	  canvasPlotProjEta[m][t][IntisMC]=new TCanvas(Form("canvasPlotProjEta_m%i_%i_MC%i",m,t, IntisMC), "canvasPlotProjEta"+Smolt[m]+"_"+PhiRegion[t], 1300, 800);
	  //	canvasPlotProjEta[m][t][IntisMC]->Divide(numC,2);
	  canvasPlotProjEta[m][t][IntisMC]->Divide(3,3);
	}
	canvasPlotRawSEProjEta[m][t][IntisMC]=new TCanvas(Form("canvasPlotRawSEProjEta_m%i_%i_MC%i",m,t, IntisMC), "canvasPlotRawSEProjEta"+Smolt[m]+ Region[t], 1300, 800);
	canvasPlotRawSEProjEta[m][t][IntisMC]->Divide(3,3);

      }
    }

    hMEEntriesvsPt[m] = new TH1F (Form("hMEEntriesvsPt_m%i", m), Form("hMEEntriesvsPt_m%i", m), numPtV0Max, NPtV0);
    hSEEntriesvsPt[m] = new TH1F (Form("hSEEntriesvsPt_m%i", m), Form("hSEEntriesvsPt_m%i", m), numPtV0Max, NPtV0);
  
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
      hWings[m][0][PtTrig] = new TH1F(Form("hWings_m%i_Data", m), Form("hWings_m%i_Data", m), numPtV0, NPtV0);
      hWings[m][1][PtTrig] = new TH1F(Form("hWings_m%i_Pythia8", m), Form("hWings_m%i_Pythia", m), numPtV0, NPtV0);
      hWingsLeft[m][0][PtTrig] = new TH1F(Form("hWingsLeft_m%i_Data", m), Form("hWingsLeft_m%i_Data", m), numPtV0, NPtV0);
      hWingsLeft[m][1][PtTrig] = new TH1F(Form("hWingsLeft_m%i_Pythia8", m), Form("hWingsLeft_m%i_Pythia", m), numPtV0, NPtV0);
      hWingsRight[m][0][PtTrig] = new TH1F(Form("hWingsRight_m%i_Data", m), Form("hWingsRight_m%i_Data", m), numPtV0, NPtV0);
      hWingsRight[m][1][PtTrig] = new TH1F(Form("hWingsRight_m%i_Pythia8", m), Form("hWingsRight_m%i_Pythia", m), numPtV0, NPtV0);

      hJetSBFitResult[m]=new TH1F(Form("hJetSBFitResult_m%i", m), Form("hJetSBFitResult_m%i", m), numPtV0, NPtV0);
      hBulkSBFitResult[m]=new TH1F(Form("hBulkSBFitResult_m%i", m), Form("hBulkSBFitResult_m%i", m), numPtV0, NPtV0);
      hAllSBFitResult[m]=new TH1F(Form("hAllSBFitResult_m%i", m), Form("hAllSBFitResult_m%i", m), numPtV0, NPtV0);
    }


    for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
      //    if (v!=PtIntervalShown) continue;
cout << "\n\n " << endl;
      if (TMath::Abs(PtTrigMax-2.5) < 0.001 && type==8 && v>4 ) continue;
      nameRawSE[m][v]="SE_";
      nameRawSE[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0];
      if (isEtaEff)       nameRawSE[m][v]+=nameisEtaEff;
      nameSE[m][v]="ME_";
      nameSE[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_AC";
      nameSESB[m][v]="ME_m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[1]+"_AC";
      nameME[m][v]="ME_m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_norm";
      nameMESB[m][v]="ME_m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[1]+"_norm";

      nameMEEtaProj[m][v] =nameME[m][v]+"_EtaProj";
      nameMEPhiProj[m][v] =nameME[m][v]+"_PhiProj"; 
      
      namePhiProjFakeJet[m][v] = nameSE[m][v] + "_phi_etaJet_FakeSB";
      namePhiProjFakeBulk[m][v] = nameSE[m][v] + "_phi_etaBI_FakeSB";
      namePhiProjFakeAll[m][v] = nameSE[m][v] + "_phi_etaAll_FakeSB";
      /*
	namePhiProjFakeJet[m][v] = nameSE[m][v] + "_phi_etaJet";
	namePhiProjFakeBulk[m][v] = nameSE[m][v] + "_phi_etaBI";
	namePhiProjFakeAll[m][v] = nameSE[m][v] + "_phi_etaAll";
      */
      namePhiProjSBJet[m][v] = nameSESB[m][v] + "_phi_etaJet";
      namePhiProjSBBulk[m][v] = nameSESB[m][v] + "_phi_etaBI";
      namePhiProjSBAll[m][v] = nameSESB[m][v] + "_phi_etaAll";

      namePhiProjJet[m][v]= nameSE[m][v] + "_phi_V0Sub_BulkSub_EffCorr";
      namePhiProjJetZYAM[m][v]= nameSE[m][v] + "_phi_V0Sub_EffCorr_BulkSubZYAM";
      namePhiProjJetFromBulkFit[m][v]= nameSE[m][v] + "_phi_V0Sub_EffCorr_BulkSubBulkFit";
      namePhiProjJetNotBulkSub[m][v]= nameSE[m][v] + "_phi_V0Sub_EffCorr";
      namePhiProjBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_Bulk_EffCorr";
      namePhiProjJetBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_JetBulkEffCorr";

      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	PtTrigMin=PtTrig+PtTrigChosen;
	if (PtTrigMin!=PtTrigChosen) continue;
	cout << "PtTrig " << PtTrig << endl;
	cout << "PtTrigMin " << PtTrigMin << endl;
	cout << "PtTrigChosen " << PtTrigChosen << endl;
	nameEtaProj[m][v][PtTrig]= nameSE[m][v] + Form("_EtaProj_PtTrig%i", PtTrig);
	nameRawEtaProj[m][v][PtTrig]= nameRawSE[m][v] + Form("_EtaProj_PtTrig%i", PtTrig);
	//	if (PtTrigMin==4 || PtTrigMin==5 || PtTrigMin>10)continue;
	if (PtTrigMin>7) continue;
	for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
	
	  PathIn[PtTrig]= Dir+"/histo/AngularCorrelation" + file[IntisMC];
	  //	  PathIn[PtTrig] += "_Jet0.75";
	  PathInMEHybrid[PtTrig] = Dir+"/histo/AngularCorrelation" + yearHybrid+"_MCTruth";
	  PathInMEXi = Dir+"/histo/AngularCorrelation" + yearXiComp+ "_MCEff";
	  if (PtBinning==1) PathInMEHybrid[PtTrig]+="_PtBinning1";
	  if(type>=0){
	    PathIn[PtTrig] +="_"+tipo[type];
	    PathIn[PtTrig] +=Srap[israp];
	    PathIn[PtTrig] +=SSkipAssoc[SkipAssoc];
	    PathInMEHybrid[PtTrig] += "_"+tipo[type]+ Srap[israp]+SSkipAssoc[0];
	    PathInMEXi+= "_"+tipo[8]+ Srap[israp]+SSkipAssoc[SkipAssoc];
	  }
	  PathIn[PtTrig]+= hhCorr[ishhCorr];
	  if (isBkgParab)	  PathIn[PtTrig]+= "_isBkgParab";
	  if (isSysDef){
	    if (isDefaultSel==1) 	  PathIn[PtTrig]+= Form("_SysT%i_SysV0Default_Sys%i_PtMin%.1f", sysTrigger, sysang, PtTrigMin);
	    else if (isDefaultSel==2) 	  PathIn[PtTrig]+= Form("_SysT%i_SysV0Tightest_Sys%i_PtMin%.1f", sysTrigger, sysang, PtTrigMin);
	    else if (isDefaultSel==3) 	  PathIn[PtTrig]+= Form("_SysT%i_SysV0Loosest_Sys%i_PtMin%.1f", sysTrigger, sysang, PtTrigMin);
	  }
	  else 	  PathIn[PtTrig]+= Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin);
	  PathIn[PtTrig] +="_Output";
	  PathInMEHybrid[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin)+"_Output.root";
	  PathInMEXi += hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, 2.)+"_Output_IsParticleTrue_IsEtaEff.root";
	  if (IsParticleTrue)  PathIn[PtTrig]+= "_IsParticleTrue";
	  PathIn[PtTrig]+= sTriggerPDGChoice[TriggerPDGChoice];
	  if  (IsMEFromHybrid) PathIn[PtTrig]+= "_IsMEFromHybrid";
	  if (isMEFromK0s) PathIn[PtTrig]+= "_IsMEFromK0s";
	  PathInComp[PtTrig] = PathIn[PtTrig];
	  if (isMEFromCorrectCentrality) PathIn[PtTrig]+= "_IsMEFromCorrectCentrality";
	  if (isEPOSEff) PathIn[PtTrig]+= "_EPOS";
	  if (isSidebandsAnalysis) PathIn[PtTrig]+= "_Sidebands";
	  if (isMEFromPeak)  PathIn[PtTrig]+= "_IsMEFromPeak";
	  if (isMEFrom13TeV) PathIn[PtTrig] += "_IsMEFrom13TeV";
	  if (SidebandsSide ==1) PathIn[PtTrig]+= "_LeftSB";
	  else if (SidebandsSide ==2) PathIn[PtTrig]+= "_RightSB";
	  if (isEtaEff) PathIn[PtTrig]+= "_IsEtaEff";
	  if (isTrigEff) PathIn[PtTrig]+= "_IsTrigEff";
	  if (MultBinning!=0) PathIn[PtTrig]+=Form("_MultBinning%i", MultBinning);
	  //	  PathIn[PtTrig] += "_thinptbins";
	  if (isNewDEtaJet)	  PathIn[PtTrig] += "_NewdEtaChoice";
	  if (isMC && !isEfficiency && !isMCForNorm) PathIn[PtTrig]+= "_MCPrediction";
	  if (isMC && !isEfficiency && !isMCForNorm) PathInComp[PtTrig]+= "_MCPrediction";
	  //	  PathIn[PtTrig] += "_Try1";
	  //	  PathInComp[PtTrig]+= "_Try";
	  if (VarRange!=0) PathIn[PtTrig]+= Form("_VarRange%i", VarRange);
	  if (VarRange!=0) PathInComp[PtTrig]+= Form("_VarRange%i", VarRange);
	  PathIn[PtTrig]+= ".root";
	  PathInComp[PtTrig]+= ".root";
	  if (!ishhCorr && isEnlargedDeltaEta)	  PathIn[PtTrig]= Dir+"/histo/AngularCorrelation" + file[IntisMC] + hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin)+"_DeltaEtaPhiEnlarged_Output.root";
	  cout <<"path in : " <<  PathIn[PtTrig] << endl;
	  filein[PtTrig]= new TFile(PathIn[PtTrig], "");

	  if (isTrigEff){
	    histoNTrigger = (TH1F*)filein[PtTrig]->Get("fHistNTrigger");
	    histoNTriggerMult = (TH1F*)filein[PtTrig]->Get("fHistNTriggerMult");
	    if (!histoNTrigger || !histoNTriggerMult) {cout << "No Ntrigger info saved " << endl; return;}
	    if (m<nummolt)	NTrigger[m][IntisMC] = histoNTrigger->GetBinContent(m+1);
	    else NTrigger[m][IntisMC] = histoNTriggerMult->GetBinContent(1);
	    cout << NTrigger[m][IntisMC]<< endl;
	    //	    return;
	  }

	  if(isEtaEffComp){
	    cout <<"path in of file with histos obtained from neglecting eta-dependence of efficiency: " <<  PathInComp[PtTrig] << endl;
	    fileinComp[PtTrig]= new TFile(PathInComp[PtTrig], "");
	  }

	  if (isCompWithMEFromHybrid){
	    //	    PathInMEHybrid[PtTrig]= "FinalOutput/DATA2016/histo/AngularCorrelation161718_MD_hXi_Hybrid_MCTruth_Xi_Eta0.8_AllAssoc_SysT0_SysV00_Sys0_PtMin2.0_Output.root";
	    cout <<"path in of file with histos obtained from hybrid run: " <<  PathInMEHybrid[PtTrig] << endl;
	    fileinCompMEHybrid[PtTrig]= new TFile(PathInMEHybrid[PtTrig], "");
	    if (! fileinCompMEHybrid[PtTrig]) {cout << " no file " << PathInMEHybrid[PtTrig] << endl; return;}
       	  }

	  if (isCompWithMEFromXi){
	    cout <<"path in of file with Xi histos: " <<  PathInMEXi << endl;
	    fileinCompMEXi = new TFile(PathInMEXi, "");
	    if (!fileinCompMEXi) {cout << " no file " << PathInMEXi << endl; return;}
	  }

	  hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameRawSE[m][v]);
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameSE[m][v]);
	  hDeltaEtaDeltaPhi_SEbinsSB[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameSESB[m][v]);
	  //	  else 	  hDeltaEtaDeltaPhi_SEbinsSB[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameSE[m][v]);

	  if (isEtaEffComp){
	    hDeltaEtaDeltaPhi_MEbinsComp[m][v][PtTrig]= (TH2F*)fileinComp[PtTrig]->Get(nameME[m][v]);
	    hDeltaEtaDeltaPhi_MEbinsComp[m][v][PtTrig]->SetName(nameME[m][v]+"_Comp");
	  }

	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameME[m][v]);
	  hDeltaEtaDeltaPhi_MEbins[multchosen][v][PtTrig]= (TH2F*)filein[PtTrig]->Get("ME_m"+ Smolt[multchosen]+"_v"+SPtV0[v]+Ssideband[0]+"_norm");

	  if (isCompWithMEFromHybrid){
	    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetName("Dummy");
	    hDeltaEtaDeltaPhi_MEbinsHybrid[m][v][PtTrig]= (TH2F*)fileinCompMEHybrid[PtTrig]->Get(nameME[m][v]);
	    if (!hDeltaEtaDeltaPhi_MEbinsHybrid[m][v][PtTrig]) {cout << "m " << m << " v " << v << " name " << nameME[m][v] << " no histo from hybrid" << endl; return;}
	    hDeltaEtaDeltaPhi_MEbinsHybrid[m][v][PtTrig]->SetName(nameME[m][v]+"_Hybrid");
	    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetName(nameME[m][v]);
	  }

	  if (isCompWithMEFromXi && NPtV0[v]>=0.5){
	  //	  if (isCompWithMEFromXi){
	    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetName("Dummy");
	    hDeltaEtaDeltaPhi_MEbinsXi[m][v][PtTrig]= (TH2F*)fileinCompMEXi->Get(nameME[m][v]);
	    if (!hDeltaEtaDeltaPhi_MEbinsXi[m][v][PtTrig]) {cout << "m " << m << " v " << v << " name " << nameME[m][v] << " no histo from Xi" << endl; return;}
	    hDeltaEtaDeltaPhi_MEbinsXi[m][v][PtTrig]->SetName(nameME[m][v]+"_Xi");
	    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetName(nameME[m][v]);
	  }
	  if (isSidebands)  hDeltaEtaDeltaPhi_MEbinsSB[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameMESB[m][v]);
	  //	  else 	  hDeltaEtaDeltaPhi_MEbinsSB[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameME[m][v]);
	  if (v==PtV0Min)  hDeltaEtaDeltaPhi_MEbins[m][numPtV0Chosen][PtTrig]= (TH2F*)filein[PtTrig]->Get("ME_m"+ Smolt[m]+"_v"+SPtV0[numPtV0Chosen]+Ssideband[0]+"_norm");
	  if (!hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]){cout << "no Raw SE 2D histo for v = "<<v << " name of the histo " << nameRawSE[m][v] << endl;  return;}
	  if (!hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]){cout << "no SE 2D histo for v = "<<v << " name of the histo " << nameSE[m][v] << endl;  return;}
	  if (!hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]){cout << "no ME 2D histo for v = " << v << " name of the histo " << nameME[m][v] <<  endl;  return;}
	  if (isSidebands){
	    if (!hDeltaEtaDeltaPhi_SEbinsSB[m][v][PtTrig]){cout << "no SE SB 2D histo for v = "<<v << " name of the histo " << nameSE[m][v] << endl;  return;}
	    if (!hDeltaEtaDeltaPhi_MEbinsSB[m][v][PtTrig]){cout << "no ME SB 2D histo for v = " << v << endl;  return;}
	  }
	  if (!hDeltaEtaDeltaPhi_MEbins[m][numPtV0Chosen][PtTrig]){cout << "no ME 2D histo (fixed Pt for comparison) for v = " <<numPtV0Chosen<< " name " << "ME_m"<< Smolt[m]<<"_v"<<SPtV0[numPtV0Chosen]<<Ssideband[0]<<"_norm"<< endl;  return;}

	  //projection along eta of the Mixed Event distribution
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->ProjectionX(nameMEEtaProj[m][v],0,-1, "E");
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetMarkerColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetMarkerStyle(33);
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetNbinsY());

	  hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Clone(nameMEEtaProj[m][v]+"_Ratio");
	  //	  if (v==1)	{
	  hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][numPtV0Chosen][PtTrig]->ProjectionX(nameMEEtaProj[m][numPtV0Chosen]+ "_Master",0,-1, "E");
	  hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetNbinsY());
	  //	  }
	  hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]);

	  for (Int_t i=1; i<=   hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->SetBinError(i, sqrt(TMath::Abs(pow( hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetBinError(i),2)-pow( hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->GetBinError(i),2)))/ hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->GetBinContent(i));
	  }

	  if (m==nummolt){
	    hDeltaEtaDeltaPhi_MEbins_EtaProjMasterMolt[v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[multchosen][v][PtTrig]->ProjectionX(nameMEEtaProj[multchosen][v]+ "_MasterMolt",0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjMasterMolt[v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[multchosen][v][PtTrig]->GetNbinsY());
	  }
	  if (m!=multchosen) {
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->ProjectionX(nameMEEtaProj[m][v]+ "_MoltRatio",0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetNbinsY());
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_EtaProjMasterMolt[v][PtTrig]);
	  }
	  cout << "ho preso isto" << endl;

	  //projection along phi of the Mixed Event distribution
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->ProjectionY(nameMEPhiProj[m][v],0,-1, "E");

	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->SetMarkerColor(ColorPt[v]);
	  cout << "ho preso isto" << endl;
	  if (v==PtV0Min){
	    hDeltaEtaDeltaPhi_MEbins_PhiProj[m][numPtV0Chosen][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][numPtV0Chosen][PtTrig]->ProjectionY(nameMEPhiProj[m][numPtV0Chosen],0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]=(TH1F*) hDeltaEtaDeltaPhi_MEbins_PhiProj[m][numPtV0Chosen][PtTrig]->Clone(nameMEPhiProj[m][numPtV0Chosen]+"_Master");
	    hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]->Rebin(2);
	  }
	  cout << "ho preso isto" << endl;
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->Rebin(2);
	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->Clone(nameMEPhiProj[m][v]+"_Ratio");

	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]);
	  cout << "ho preso isto" << endl;
	  for (Int_t i=1; i<=   hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->GetNbinsX(); i++){
	    //	    hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->SetBinError(i,sqrt( TMath::Abs( pow(hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->GetBinError(i),2)- pow(hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]->GetBinError(i),2)))/hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]->GetBinContent(i));
	  }
	  cout << "ho preso isto" << endl;

	  //other projections
	  if (isSidebands){
	    hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjSBJet[m][v]);
	    //	    hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjFakeJet[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjSBBulk[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjFakeBulk[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjSBAll[m][v]);
	    hDeltaEtaDeltaPhi_PhiProjAllFakeSB[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjFakeAll[m][v]);

	    polJetSB[m][v][PtTrig] = new TF1 (Form("pol0Jet_m%i_v%i", m, v),Form("pol0Jet_m%i_v%i", m, v), -TMath::Pi()/2, 3./2*TMath::Pi());
	    polBulkSB[m][v][PtTrig] = new TF1 (Form("pol0Bulk_m%i_v%i", m, v), Form("pol0Bulk_m%i_v%i", m, v),-TMath::Pi()/2, 3./2*TMath::Pi());
	    polAllSB[m][v][PtTrig] = new TF1 (Form("pol0All_m%i_v%i", m, v),Form("pol0All_m%i_v%i", m, v), -TMath::Pi()/2, 3./2*TMath::Pi());

	    polJetSB[m][v][PtTrig]->SetLineWidth(0.2);
	    polJetSB[m][v][PtTrig]->SetLineColor(1);
	    polAllSB[m][v][PtTrig]->SetLineWidth(0.2);
	    polAllSB[m][v][PtTrig]->SetLineColor(1);
	    polBulkSB[m][v][PtTrig]->SetLineWidth(0.2);
	    polBulkSB[m][v][PtTrig]->SetLineColor(1);

	  }

	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetZYAM[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetFromBulkFit[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetNotBulkSub[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]);

	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]+"_RelError");
	  hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]+"_RelError");
	  hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]+"_RelError");

	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->ProjectionX(nameRawEtaProj[m][v][PtTrig], 0,-1, "E");
	  hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->ProjectionX(nameRawEtaProj[m][v][PtTrig]+"PHalf",hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetYaxis()->FindBin(+TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetYaxis()->FindBin(+3./2*TMath::Pi()- 0.001)  , "E");
	  hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->ProjectionX(nameRawEtaProj[m][v][PtTrig]+"NHalf",hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetYaxis()->FindBin(-TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetYaxis()->FindBin(+1./2*TMath::Pi()- 0.001)  , "E");

	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig], 0,-1, "E");
	  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig]+"PHalf",hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetYaxis()->FindBin(+TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetYaxis()->FindBin(+3./2*TMath::Pi()- 0.001)  , "E");
	  //	  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig]+"PHalf",hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetYaxis()->FindBin(3+0.001) ,hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetYaxis()->FindBin(3./2 * TMath::Pi()- 0.001)  , "E");
	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig]+"NHalf",hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetYaxis()->FindBin(-TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetYaxis()->FindBin(+1./2*TMath::Pi()- 0.001)  , "E");

	  pol0SideDX[m][v][PtTrig] = new TF1(Form("pol0SideDX_m%i_v%i", m, v),"pol0", 0.75, 1.15);
	  pol0SideSX[m][v][PtTrig] = new TF1(Form("pol0SideSX_m%i_v%i", m, v),"pol0", -1.15, -0.75);
	  //	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->FindBin(-1.2), 	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->FindBin(-0.75));
	  pol0Centr[m][v][PtTrig] = new TF1( Form("pol0Centr_m%i_v%i", m, v), "pol0",-0.75, 0.75);//,	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->FindBin(-0.75), 	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->FindBin(0.75));
	  pol0Centr[m][v][PtTrig]->SetLineColor(628);
	  pol0SideSX[m][v][PtTrig]->SetLineColor(kGreen+3);
	  pol0SideDX[m][v][PtTrig]->SetLineColor(kGreen+3);

	  cout << "ho preso isto" << endl;

	  cout << "this rebin is performed for visualization purposes only " << endl;
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Rebin(2);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Rebin(2);
	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Rebin(2);

	  canvasPlotSB[m][IntisMC]->cd(v+1);
	  if (isSidebands){
	    hDeltaEtaDeltaPhi_SEbinsSB[m][v][PtTrig]->SetTitle(SmoltBis[m] + " AC v"+SPtV0[v]+ " sidebands");
	    hDeltaEtaDeltaPhi_SEbinsSB[m][v][PtTrig]->Draw("colz");
	  }

	  //canvasPlot[m]->cd(PtTrig+1);
	  canvasPlot[m][IntisMC]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.15);
	  cout << "scelto cd " << endl;
	  //	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetLogz();
	  //	gPad->SetLogz();
	  //	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetTitle("AC v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetTitle(SmoltBis[m] + " AC v"+SPtV0[v]);
	  TString TitleOld="Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin) +MCOrNot[IntisMC];
	  TString TitleNew =SmoltBis[m] + " p^{K^{0}_{S}}_{T}  [" + SPtV0[v]+") GeV/c" ;
	  if (type==8) TitleNew ="p^{#Xi}_{T}  [" + SPtV0[v]+") GeV/c" ;
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->SetTitle(TitleNew);
	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->SetTitle(TitleNew);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->SetTitle(TitleNew);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetTitle(TitleNew);
	  //	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.5, 1.5);
	  if (year.Index("Fio")!=-1) 	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.5, 1.5);
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->Draw("colz");
	  tlinePhiDx->Draw();
	  tlinePhiSx->Draw();
	  tlineEtaDx->Draw();
	  tlineEtaSx->Draw();
	  tlineEtaInclDx->Draw();
	  tlineEtaInclSx->Draw();

	  canvasPlotME[m][IntisMC]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.15);
	  //	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetTitle("AC v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetTitle("ME " + SmoltBis[m]+ " v"+SPtV0[v]);
	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->Draw("colz");

	  canvasPlotSE[m][IntisMC]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.15);
	  //	  hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->SetTitle("Raw SE v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->SetTitle("Raw SE "+ SmoltBis[m]+ " v"+SPtV0[v]);
	  hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->Draw("colz");

	  if (isEtaEffComp){
	    canvasPlotMEEtaProjComp[m][IntisMC]->cd(v+1);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" ME proj %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Draw("");

	    hDeltaEtaDeltaPhi_MEbins_EtaProjComp[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbinsComp[m][v][PtTrig]->ProjectionX(nameMEEtaProj[m][v]+"_Comp",0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjComp[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbinsComp[m][v][PtTrig]->GetNbinsY());
	    hDeltaEtaDeltaPhi_MEbins_EtaProjComp[m][v][PtTrig]->SetLineColor(kBlack);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjComp[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjComp[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" ME proj %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProjComp[m][v][PtTrig]->Draw("same");

	    canvasPlotMEEtaProjCompRatio[m][IntisMC]->cd(v+1);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjCompRatio[m][v][PtTrig]= (TH1F*) hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Clone(nameMEEtaProj[m][v]+"_CompRatio");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjCompRatio[m][v][PtTrig]->Divide(	  hDeltaEtaDeltaPhi_MEbins_EtaProjComp[m][v][PtTrig]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjCompRatio[m][v][PtTrig]->SetLineColor(kBlack);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjCompRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjCompRatio[m][v][PtTrig]->SetTitle(Form("ME proj ratio correct/old  %.1f < p_{T} < %.1f",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProjCompRatio[m][v][PtTrig]->Draw("same");

	  }

	  canvasPlotMESB[m][IntisMC]->cd(v+1);
	  if (isSidebands){
	    hDeltaEtaDeltaPhi_MEbinsSB[m][v][PtTrig]->SetTitle(SmoltBis[m] + "ME sidebands v"+SPtV0[v]);
	    hDeltaEtaDeltaPhi_MEbinsSB[m][v][PtTrig]->Draw("colz");
	  }

	  cout << "m " << m << " v " << v << " "  << counter << endl;
	  if (counter==1)	  legendPt->AddEntry(	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig], SPtV0[v], "pl");

	  canvasPlotSEEtaProj[IntisMC]->cd(m+1);
	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->SetMarkerColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetNbinsY());
	  //	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->Integral());
	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2* hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetMaximum());
	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" Raw SE proj %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->DrawClone("same");
	  if (v==numPtV0-1) legendPt->Draw("same");
	  /*
	    canvasPlotSEEtaProj[IntisMC]->cd(m+1+6);
	    hDeltaEtaDeltaPhi_RawSEEtaProjRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	    hDeltaEtaDeltaPhi_RawSEEtaProjRatio[m][v][PtTrig]->Draw("same");
	    if (v==numPtV0-1) legendPt->Draw("same");
	  */

	  if (isCompWithMEFromHybrid){
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbinsHybrid[m][v][PtTrig]->ProjectionX(nameMEEtaProj[m][v]+"_Hybrid",0,-1, "E");
	    cout << " proj done " << endl;
	    //	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->SetLineColor(kBlack);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->SetMarkerColor(kBlack);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbinsHybrid[m][v][PtTrig]->GetNbinsY());
	  
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetTitle("Pair acceptance");	  
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->GetYaxis()->SetTitle("Pair acceptance");	 
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->GetYaxis()->SetTitleSize(0.04);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetTitleSize(0.04);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->GetXaxis()->SetTitleSize(0.04);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetXaxis()->SetTitleSize(0.04);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.2);

	    canvasPlotMEEtaProjCompHybrid[m][IntisMC]->cd(v+1);
	    gPad->SetLeftMargin(0.15);
	    gPad->SetBottomMargin(0.15);
	    cout << " here " << endl;
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->SetMarkerStyle(27);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" ME proj %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->DrawClone("");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->SetMarkerStyle(33);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" ME proj %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]->DrawClone("same");

	    canvasPlotMEEtaProjCompHybridRatio[m][IntisMC]->cd(v+1);
	    cout << " here " << endl;
	    gPad->SetLeftMargin(0.15);
	    gPad->SetBottomMargin(0.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]=(TH1F*) 	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Clone(nameMEEtaProj[m][v]+"_RatioToHybrid");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" %.1f < p_{T} < %.1f",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->Divide(hDeltaEtaDeltaPhi_MEbins_EtaProjHybrid[m][v][PtTrig]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->Draw("");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->GetYaxis()->SetTitle("Pair acceptance ratio");	  
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->GetYaxis()->SetTitleSize(0.05);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.2);
	    pol0Fit[m][v][PtTrig] = new TF1 (Form("pol0_m%i_v%i", m, v), "pol0", -1.15, 1.15);
	    //	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToHybrid[m][v][PtTrig]->Fit(pol0Fit[m][v][PtTrig], "R0");
	    //	    TLegend * legendMEratio = new TLegend (0.7, 0.7, 0.9, 0.9);
	    //	    legendMEratio->AddEntry(pol0Fit[m][v][PtTrig], Form("q= %.2f", pol0Fit[m][v][PtTrig]->GetParameter(0)), "l");
	    tlineAtOneLong->Draw("same"); 
	  }


	  if (isCompWithMEFromXi  && NPtV0[v]>=0.5){
	    //if (isCompWithMEFromXi){
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbinsXi[m][v][PtTrig]->ProjectionX(nameMEEtaProj[m][v]+"_Xi",0,-1, "E");
	    cout << " proj done " << endl;
	    //	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->SetLineColor(kBlack);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->SetMarkerColor(kBlack);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbinsXi[m][v][PtTrig]->GetNbinsY());
	  
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->GetYaxis()->SetTitle("Pair acceptance Xi");	 
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->GetYaxis()->SetTitleSize(0.04);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->GetXaxis()->SetTitleSize(0.04);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.2);

	    canvasPlotMEEtaProjCompXi[m][IntisMC]->cd(v+1);
	    gPad->SetLeftMargin(0.15);
	    gPad->SetBottomMargin(0.15);
	    cout << " here " << endl;
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->DrawClone("");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->SetMarkerStyle(33);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" %.1f < p_{T} < %.1f",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]->DrawClone("same");

	    canvasPlotMEEtaProjCompXiRatio[m][IntisMC]->cd(v+1);
	    cout << " here " << endl;
	    gPad->SetLeftMargin(0.15);
	    gPad->SetBottomMargin(0.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]=(TH1F*) 	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Clone(nameMEEtaProj[m][v]+"_RatioToXi");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" %.1f < p_{T} < %.1f",NPtV0[v], NPtV0[v+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->Divide(hDeltaEtaDeltaPhi_MEbins_EtaProjXi[m][v][PtTrig]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->Draw("");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->GetYaxis()->SetTitle("Pair acceptance ratio K0s to Xi");	  
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->GetYaxis()->SetTitleSize(0.05);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioToXi[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.2);
	    tlineAtOneLong->Draw("same"); 
	  }

	  if (m>2){
	    canvasPlotMEEtaProj[IntisMC]->cd(m-3+1); //-3
	    canvasPlotMEEtaProj[IntisMC]->SetLeftMargin(0.15);
	    canvasPlotMEEtaProj[IntisMC]->SetBottomMargin(0.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetTitle("ME " + SmoltBis[m]+"%" + " v"+ SPtV0[v]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Draw("same");
	    if (v==numPtV0-1) legendPt->Draw("same");
	    canvasPlotMEEtaProj[IntisMC]->cd(m+1); //+6
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->SetTitle("ME ratio " + Smolt[m]+ Form(" to %.1f < p_{T} < %.1f", NPtV0[numPtV0Chosen],  NPtV0[numPtV0Chosen+1]));
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->Draw("same");
	    if (v==numPtV0-1) legendPt->Draw("same");
	  }

	  canvasPlotMEEtaProjMolt[IntisMC]->SetLeftMargin(0.15);
	  canvasPlotMEEtaProjMolt[IntisMC]->SetBottomMargin(0.25);
	  canvasPlotMEEtaProjMolt[IntisMC]->cd(v+1);
	  //	  if (m!=nummolt ){
	  if (m!=multchosen){
	    cout << "here molt " << endl;
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->SetLineColor(Colormult[m]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->SetMarkerColor(Colormult[m]);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->SetMarkerStyle(33);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->SetTitle("ME ratio to " + Smolt[multchosen]+Form(" %.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));	
	    if (v==PtV0Min)	    legendMult->AddEntry(	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig], SmoltBis[m], "pl");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.9, 1.1);
	    if (type==0 && year=="1617_hK0s") 	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.9, 1.1);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1);
	    //	    if (m==0 || m ==1 || m==2)	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->Draw("same");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatioMolt[m][v][PtTrig]->Draw("same");
	  }
	  if (m==0) tlineAtOneLong->Draw("same");
	  if (m==0) legendMult->Draw("same");
	  if (isHM && m==2) legendMult->Draw("same");
	  

	  canvasPlotMEPhiProj[IntisMC]->cd(m+1);
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,40);
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->Draw("same");
	  if (v==numPtV0-1) legendPt->Draw("same");
	  canvasPlotMEPhiProj[IntisMC]->cd(m+1+6);
	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->Draw("same");
	  if (v==numPtV0-1) legendPt->Draw("same");

	  //	canvasPlot[m]->cd(PtTrig+1+7);
	  cout << "setting to zero last bin " << endl;
	  //	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  //	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  //	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->SetLineColor(860);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetLineColor(kBlue);
	  hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->SetLineColor(kGreen+3);

	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetLineColor(628);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetLineColor(418);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetLineColor(868);

	  hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->SetMarkerColor(860);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetMarkerColor(628);
	  hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->SetMarkerColor(kGreen+3);

	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetMarkerColor(628);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetMarkerColor(418);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetMarkerColor(868);

	  if (hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig] && hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]&& hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]){
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

	  for (Int_t t =0; t<3; t++){
	    if (isSidebands){
	      canvasPlotProjSB[m][IntisMC][t]->cd(v+1);
	      if (type==8 && v!=numPtV0Max-1) 	  canvasPlotProjSB[m][IntisMC][t]->cd(v);
	      gStyle->SetOptStat(0);
	      gPad->SetLeftMargin(0.2);
	      if (m==nummolt && v==PtV0Min && t==0){
		legendSB->AddEntry(	  hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig], "Sidebands", "pl");
		legendSB->AddEntry(	  hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig], "Peak region", "pl");
	      }
	      PtTitle =Form("%.1f < p_{T} < %.1f GeV/c" , NPtV0[v], NPtV0[v+1]);
	      if (t==0){
		hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]->SetTitle(PtTitle);	 
		hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig]->SetTitle(PtTitle);	 
		hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]->GetYaxis()->SetTitle("Arbitrary units");
		hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig]->GetYaxis()->SetTitle("Arbitrary units");
		hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.4);
		hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.4);
		hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]->Draw("same");
		hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig]->Draw("same");
	      }
	      else if (t==1){
		hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]->SetTitle(PtTitle);	 
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[m][v][PtTrig]->SetTitle(PtTitle);	 
		hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]->GetYaxis()->SetTitle("Arbitrary units");
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[m][v][PtTrig]->GetYaxis()->SetTitle("Arbitrary units");
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.4);
		hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.4);
		hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]->Draw("same");
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[m][v][PtTrig]->Draw("same");
	      }
	      else {
		hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]->SetTitle(PtTitle);	 
		hDeltaEtaDeltaPhi_PhiProjAllFakeSB[m][v][PtTrig]->SetTitle(PtTitle);	 
		hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]->GetYaxis()->SetTitle("Arbitrary units");
		hDeltaEtaDeltaPhi_PhiProjAllFakeSB[m][v][PtTrig]->GetYaxis()->SetTitle("Arbitrary units");
		hDeltaEtaDeltaPhi_PhiProjAllFakeSB[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.4);
		hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.4);
		hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjAllFakeSB[m][v][PtTrig]->Rebin(2);
		hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]->Draw("same");
		hDeltaEtaDeltaPhi_PhiProjAllFakeSB[m][v][PtTrig]->Draw("same");
	      }
	      legendSB->Draw("");

	      hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]= (TH1F*) 	    hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]->Clone(namePhiProjFakeAll[m][v]+"_Ratio");
	      hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]= (TH1F*) 	    hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]->Clone(namePhiProjFakeBulk[m][v]+"_Ratio");
	      hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]= (TH1F*) 	    hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]->Clone(namePhiProjFakeJet[m][v]+"_Ratio");

	      hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->Sumw2();
	      hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->Sumw2();
	      hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->Sumw2();

	      hDeltaEtaDeltaPhi_PhiProjAllSB[m][v][PtTrig]->Sumw2();
	      hDeltaEtaDeltaPhi_PhiProjBulkSB[m][v][PtTrig]->Sumw2();
	      hDeltaEtaDeltaPhi_PhiProjJetSB[m][v][PtTrig]->Sumw2();

	      hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->Divide(	    hDeltaEtaDeltaPhi_PhiProjAllFakeSB[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->Divide(	    hDeltaEtaDeltaPhi_PhiProjBulkFakeSB[m][v][PtTrig]);
	      hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->Divide(	    hDeltaEtaDeltaPhi_PhiProjJetFakeSB[m][v][PtTrig]);


	      canvasPlotProjRatioSB[m][IntisMC][t]->cd(v+1);
	      if (type==8 && v!=numPtV0Max-1) 	  canvasPlotProjRatioSB[m][IntisMC][t]->cd(v);
	      gStyle->SetOptFit(0);
	      gStyle->SetOptStat(0);
	      gPad->SetLeftMargin(0.2);
	      PtTitle =Form("%.1f < p_{T} < %.1f GeV/c" , NPtV0[v], NPtV0[v+1]);
	      if (t==0){
		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->SetLineColor(628);
		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->SetMarkerColor(628);
		//		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.5, 1.5);
		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5);
		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetTitle("Sidebands/Fake SB");
		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->SetTitle(PtTitle);
		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->Fit(polJetSB[m][v][PtTrig], "R+");
		hDeltaEtaDeltaPhi_PhiProjJetFakeSBRatio[m][v][PtTrig]->Draw("same");
		hJetSBFitResult[m]->SetBinContent(v+1, polJetSB[m][v][PtTrig]->GetParameter(0));
		hJetSBFitResult[m]->SetBinError(v+1, polJetSB[m][v][PtTrig]->GetParError(0));
	      }
	      else if (t==1){
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->SetLineColor(kGreen+2);
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->SetMarkerColor(kGreen+2);

		//		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.5, 1.5);
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.5);
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetTitle("Sidebands/Fake SB");
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->SetTitle(PtTitle);
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->Fit(polBulkSB[m][v][PtTrig], "R+");
		hDeltaEtaDeltaPhi_PhiProjBulkFakeSBRatio[m][v][PtTrig]->Draw("same");
		hBulkSBFitResult[m]->SetBinContent(v+1, polBulkSB[m][v][PtTrig]->GetParameter(0));
		hBulkSBFitResult[m]->SetBinError(v+1, polBulkSB[m][v][PtTrig]->GetParError(0));

	      }
	      else {
		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->SetLineColor(kBlue);
		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->SetMarkerColor(kBlue);

		//		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.5, 1.5);
		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5);
		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->GetYaxis()->SetTitle("Sidebands/Fake SB");
		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->SetTitle(PtTitle);
		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->Fit(polAllSB[m][v][PtTrig], "R+");
		hDeltaEtaDeltaPhi_PhiProjAllFakeSBRatio[m][v][PtTrig]->Draw("same");
		hAllSBFitResult[m]->SetBinContent(v+1, polAllSB[m][v][PtTrig]->GetParameter(0));
		hAllSBFitResult[m]->SetBinError(v+1, polAllSB[m][v][PtTrig]->GetParError(0));

	      }
	    }
	  }

	  canvasPlotProj[m][IntisMC]->cd(v+1);
	  if (type==8 && v!=numPtV0Max-1) 	  canvasPlotProj[m][IntisMC]->cd(v);
	  if (isMC && !isEfficiency) 	  canvasPlotProj[m][IntisMC]->cd(v+1);
	  cout << "cd chosen " << endl;

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
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Sumw2();
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Scale(1./	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Integral());
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Fit(GaussianEta[m][v][PtTrig], "R");

	  cout << "\nbaseline value" <<  Baseline[m][v][PtTrig]->GetParameter(0)<< " +- " <<  Baseline[m][v][PtTrig]->GetParError(0)<< endl;
	  cout << "\njet fit: chi square " << Gaussian[m][v][PtTrig]->GetChisquare() << " NDF " << Gaussian[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << Gaussian[m][v][PtTrig]->GetChisquare()/Gaussian[m][v][PtTrig]->GetNDF() << endl;
	  cout << "\nAwaySide fit: chi square " << GaussianAS[m][v][PtTrig]->GetChisquare() << " NDF " << GaussianAS[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << GaussianAS[m][v][PtTrig]->GetChisquare()/GaussianAS[m][v][PtTrig]->GetNDF() << endl;

	  //	  if (	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetMaximum() > 	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetMaximum()) 	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Draw("");

	  if (m<=1 && MultBinning==1) LimSupY[0]=0.001;
	  if (type==0)	{
	    LimSupY[0]=0.014; //was 0.014
	    LimInfY[0]=-0.001; //was 0.002
	    if (m>=3 && m<5){
	      LimSupY[0]=0.014;
	      LimInfY[0]=-0.001;
	    }
	    if (NPtV0[v] >=1.2 && NPtV0[v] < 2.5)  LimSupY[0]=0.01;
	    else if (NPtV0[v] >= 2.5)  LimSupY[0]=0.008;
	    if (isHM) {
	      LimSupY[0] = 0.02;
	      if (NPtV0[v] >=1.2 && NPtV0[v] < 2.5)  LimSupY[0]=0.015;
	      else if (NPtV0[v] >= 2.5)  LimSupY[0]=0.008;
	    }
	    if (year.Index("Fio")!=-1){
	      if (NPtV0[v] <1.2)  LimSupY[0]=0.012;
	      if (NPtV0[v] >=1.2 && NPtV0[v] < 2.5)  LimSupY[0]=0.008;
	      else if (NPtV0[v] >= 2.5)  LimSupY[0]=0.003;
	      if (m>2) {
		if (NPtV0[v] <1.2)  LimSupY[0]=0.004;
		if (NPtV0[v] >=1.2 && NPtV0[v] < 2.5)  LimSupY[0]=0.002;
		else if (NPtV0[v] >= 2.5)  LimSupY[0]=0.001;
	      }
	    }
	  }
	  else {
	    if (isHM){
	      //	      LimSupY[0] = 0.005;
	      LimSupY[0] = 0.002;
	      //	      if (NPtV0[v] >=2 && NPtV0[v] < 4)  LimSupY[0]=0.0035;
	      if (NPtV0[v] >=2 && NPtV0[v] < 4)  LimSupY[0]=0.001;
	      //	      else if (NPtV0[v] >= 2.5)  LimSupY[0]=0.0025;
	      else if (NPtV0[v] >= 2.5)  LimSupY[0]=0.001;
	    }
	    else if (m==0) LimSupY[0] = 0.0012;
	  }

	  cout << LimSupY[0] << endl;
	  //	  return;
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]); //-0.004 for hh and hK0s
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(LimInfY[ishhCorr], LimSupY[ishhCorr]);
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

	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetTitleOffset(2.2); //was 1.8
	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->GetYaxis()->SetTitleOffset(2.2); 
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetTitleOffset(2.2); 
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetTitleOffset(2.2); 

	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetStats(0);
	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->SetStats(0);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->SetStats(0);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->SetStats(0);

	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetMarkerColor(kBlue);
	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]    ->SetLineColor(kBlue);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetMarkerColor(628);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetLineColor(628);

	  //	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Draw("");
	  pol0ZYAMOrigin[m][v][PtTrig]=(TF1*)filein[PtTrig]->Get(Form("pol0ZYAM_m%i_v%i",m,v));
	  if (!pol0ZYAMOrigin[m][v][PtTrig]) {cout << "pol0Zyam not saved in input file " << endl;}
	  if (pol0ZYAMOrigin[m][v][PtTrig]){
	    pol0ZYAM[m][v][PtTrig]= new TF1("pol0",Form("pol0ZYAM_m%i_v%i",m,v), -TMath::Pi()/2,3./2*TMath::Pi());
	    pol0ZYAM[m][v][PtTrig]->SetLineColor(kRed);
	    pol0ZYAM[m][v][PtTrig]->SetParameter(0, pol0ZYAMOrigin[m][v][PtTrig]->GetParameter(0)/NTrigger[m][IntisMC]*2); //2 for rebin 
	    pol0ZYAMOrigin[m][v][PtTrig]->SetParameter(0, pol0ZYAMOrigin[m][v][PtTrig]->GetParameter(0)/NTrigger[m][IntisMC]*2); //2 for rebin 
	    pol0ZYAM[m][v][PtTrig]->SetLineColor(860);
	    pol0ZYAM[m][v][PtTrig]->SetLineWidth(0.3);
	    //	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Fit(pol0ZYAM[m][v][PtTrig], "R+");
	  }

	  pol0BulkBis[m][v][PtTrig]= new TF1("pol0",Form("pol0BulkBis_m%i_v%i",m,v), -1,1);
	  pol0BulkBis[m][v][PtTrig]->SetLineColor(kGreen+3);
	  //	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Fit(pol0BulkBis[m][v][PtTrig], "R+");

	  hDeltaEtaDeltaPhi_PhiProjJetNotBulkSub[m][v][PtTrig]->Draw("");
	  //	  if (pol0ZYAMOrigin[m][v][PtTrig])	  pol0ZYAM[m][v][PtTrig]->Draw("same");
	  //	  pol0ZYAMOrigin[m][v][PtTrig]->Draw("same");
	  //	  pol0BulkBis[m][v][PtTrig]->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Draw("same");
	  //	  hDeltaEtaDeltaPhi_PhiProjJetZYAM[m][v][PtTrig]    ->Draw("same");
	  //	  hDeltaEtaDeltaPhi_PhiProjJetFromBulkFit[m][v][PtTrig]    ->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Draw("same");
	  tlinePhiBase->Draw("");

	  canvasPlotProjRatioJet[m][IntisMC]->cd(v+1);
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
	      hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" AC proj %.1f <p_{T} <%.1f GeV/c", NPtV0[v], NPtV0[v+1])); 
	      hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1); 
	      hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	      hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->DrawClone("");
	      lineEtam075=new TLine(-0.75, 0, -0.75, 1* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
	      lineEta075=new TLine(0.75, 0, 0.75, 1* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
							    
	    }
	    else if (t==1){
	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Sumw2();
	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Scale(1./	  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Integral());
	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2* hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetMaximum());
	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" AC proj(+) %.1f <p_{T} <%.1f GeV/c", NPtV0[v], NPtV0[v+1])); 
	      //	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form("%.1f <p_{T} <%.1f GeV/c", NPtV0[v], NPtV0[v+1])); 
	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1); 
	      //	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Fit(pol0Centr[m][v][PtTrig], "R+");
	      //	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Fit(pol0SideSX[m][v][PtTrig], "R+");
	      //	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Fit(pol0SideDX[m][v][PtTrig], "R+");
	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	      hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->DrawClone("");
	      lineEtam075=new TLine(-0.75, 0, -0.75,  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetMaximum());
	      lineEta075=new TLine(0.75, 0,  0.75,  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetMaximum());

	    }
	    else if (t==2){
	      hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->Sumw2();
	      hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->Scale(1./	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->Integral());
	      hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2* hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetMaximum());
	      cout << "\n\n****** fit to dEta proj in the away side " << endl;
	      hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" AC proj(-) %.1f <p_{T} <%.1f GeV/c", NPtV0[v], NPtV0[v+1])); 
	      //hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form("%.1f <p_{T} <%.1f GeV/c", NPtV0[v], NPtV0[v+1])); 
	      hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1); 
	      hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	      hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->DrawClone("same");
	      lineEtam075=new TLine(-0.75, 0, -0.75, hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetMaximum());
	      lineEta075=new TLine(0.75, 0, 0.75,  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetMaximum());

	    }

	    lineEta075->SetLineColor(628);
	    lineEtam075->SetLineColor(628);
	    lineEta075->Draw("same");
	    lineEtam075->Draw("same");
	  }


	  for (Int_t t=0; t<3; t++){
	    canvasPlotRawSEProjEta[m][t][IntisMC]->cd(v+1);
	    gPad->SetLeftMargin(0.15);

	    if (t==0){
	      hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->SetLineColor(kBlack);
	      hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->SetMarkerColor(kBlack);
	      hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1); 
	      hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	      hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->Integral());
	      hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2* hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetMaximum());
	      hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->Draw("");						
	      lineEtam075=new TLine(-0.75, 0, -0.75, hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetMaximum());
	      lineEta075=new TLine(0.75, 0, 0.75,  hDeltaEtaDeltaPhi_RawSEEtaProj[m][v][PtTrig]->GetMaximum());

							    
	    }
	    else if (t==1){
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->SetLineColor(kBlack);
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->SetMarkerColor(kBlack);
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->SetMarkerStyle(33);
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1); 
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	      //	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2* hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->GetMaximum());
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" Raw SE proj(+) %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->Scale(2./(hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->GetBinContent(hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->FindBin(0.05)) + hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->GetBinContent(hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->FindBin(-0.05))));
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2);
	      hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	      hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->Draw("");				     
	      hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->DrawClone("same");				    
	      lineEtam075=new TLine(-0.75, 0, -0.75, 1.2);
	      lineEta075=new TLine(0.75, 0, 0.75, 1.2);

	    }
	    else if (t==2){
	      hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->SetLineColor(kBlack);
	      hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->SetMarkerColor(kBlack);
	      hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->GetXaxis()->SetTitleOffset(1); 
	      hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->GetXaxis()->SetRangeUser(-1.15, 1.15);
	      hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2* hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->GetMaximum());
	      hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" Raw SE proj(-) %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	      hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->Draw("");			      		
	      lineEtam075=new TLine(-0.75, 0, -0.75, hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->GetMaximum());
	      lineEta075=new TLine(0.75, 0, 0.75,  hDeltaEtaDeltaPhi_RawSEEtaProjNHalf[m][v][PtTrig]->GetMaximum());

	    }

	    lineEta075->SetLineColor(628);
	    lineEtam075->SetLineColor(628);
	    lineEta075->Draw("same");
	    lineEtam075->Draw("same");
	  }

	  canvasPlotRawSEProjEta[m][3][IntisMC]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  hDeltaEtaDeltaPhi_RawSEEtaProjPHalfRatioToME[m][v][PtTrig]=(TH1F*) hDeltaEtaDeltaPhi_RawSEEtaProjPHalf[m][v][PtTrig]->Clone(nameRawEtaProj[m][v][PtTrig]+"PHalfRatioToME");
	  hDeltaEtaDeltaPhi_RawSEEtaProjPHalfRatioToME[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]);
	  hDeltaEtaDeltaPhi_RawSEEtaProjPHalfRatioToME[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.9, 1.1);
	  hDeltaEtaDeltaPhi_RawSEEtaProjPHalfRatioToME[m][v][PtTrig]->SetTitle(SmoltBis[m] + Form(" SE(+)/ME proj %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]));
	  hDeltaEtaDeltaPhi_RawSEEtaProjPHalfRatioToME[m][v][PtTrig]->Draw("");				     
	  lineEtam075=new TLine(-0.75, 0.9, -0.75, 1.1);
	  lineEta075=new TLine(0.75, 0.9, 0.75, 1.1);
	  lineEta075->SetLineColor(628);
	  lineEtam075->SetLineColor(628);
	  lineEta075->Draw("same");
	  lineEtam075->Draw("same");
	  tlineAtOneLong->Draw("same"); 

	  cout << " hWings m: " <<m << " v: " << v <<  endl;
	  hWings[m][IntisMC][PtTrig]->SetBinContent(v+1,(pol0SideDX[m][v][PtTrig]->GetParameter(0) + pol0SideSX[m][v][PtTrig]->GetParameter(0))/ (2*pol0Centr[m][v][PtTrig]->GetParameter(0)));
	  cout << pol0Centr[m][v][PtTrig]->GetParameter(0) << " dx " << pol0SideDX[m][v][PtTrig]->GetParameter(0) << " sx " << pol0SideSX[m][v][PtTrig]->GetParameter(0)<< endl;
	  cout <<  hWings[m][IntisMC][PtTrig]->GetBinContent(v+1)<< endl;
	  hWings[m][IntisMC][PtTrig]->SetBinError(v+1,sqrt(pow(pol0SideDX[m][v][PtTrig]->GetParError(0)/2/pol0Centr[m][v][PtTrig]->GetParameter(0),2) + pow(pol0SideSX[m][v][PtTrig]->GetParError(0)/2/pol0Centr[m][v][PtTrig]->GetParameter(0),2)+ pow(hWings[m][IntisMC][PtTrig]->GetBinContent(v+1)/pol0Centr[m][v][PtTrig]->GetParameter(0) * pol0Centr[m][v][PtTrig]->GetParError(0),2) ));
	  //	  hWings[m][IntisMC][PtTrig]->SetBinError(v+1,0);

	  hWingsRight[m][IntisMC][PtTrig]->SetBinContent(v+1,pol0SideDX[m][v][PtTrig]->GetParameter(0)/pol0Centr[m][v][PtTrig]->GetParameter(0));
	  hWingsRight[m][IntisMC][PtTrig]->SetBinError(v+1, 	  hWingsRight[m][IntisMC][PtTrig]->GetBinContent(v+1) * sqrt(pow(pol0SideDX[m][v][PtTrig]->GetParError(0)/pol0SideDX[m][v][PtTrig]->GetParameter(0) ,2) + pow( pol0Centr[m][v][PtTrig]->GetParError(0)/pol0Centr[m][v][PtTrig]->GetParameter(0),2) ));

	  hWingsLeft[m][IntisMC][PtTrig]->SetBinContent(v+1,pol0SideSX[m][v][PtTrig]->GetParameter(0)/pol0Centr[m][v][PtTrig]->GetParameter(0));
	  hWingsLeft[m][IntisMC][PtTrig]->SetBinError(v+1, 	  hWingsLeft[m][IntisMC][PtTrig]->GetBinContent(v+1) * sqrt(pow(pol0SideSX[m][v][PtTrig]->GetParError(0)/pol0SideSX[m][v][PtTrig]->GetParameter(0) ,2) + pow( pol0Centr[m][v][PtTrig]->GetParError(0)/pol0Centr[m][v][PtTrig]->GetParameter(0),2) ));

	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinContent(v+1, Gaussian[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinError(v+1, Gaussian[m][v][PtTrig]->GetParError(2));
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinContent(v+1, GaussianAS[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinError(v+1, GaussianAS[m][v][PtTrig]->GetParError(2));
	  HistoWidthGaussianEta[m][IntisMC][PtTrig]->SetBinContent(v+1, GaussianEta[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussianEta[m][IntisMC][PtTrig]->SetBinError(v+1, GaussianEta[m][v][PtTrig]->GetParError(2));

	} //end of IntisMC loop
      } //end of PtTrigMin loop
    } //end of v loop

    if (isSidebands){
      canvasSBJetFitResult->cd(m+1);
      hJetSBFitResult[m]->GetYaxis()->SetRangeUser(0.8,1.2);
      hJetSBFitResult[m]->Draw("");
      canvasSBBulkFitResult->cd(m+1);
      hBulkSBFitResult[m]->GetYaxis()->SetRangeUser(0.8,1.2);
      hBulkSBFitResult[m]->Draw("");
      canvasSBAllFitResult->cd(m+1);
      hAllSBFitResult[m]->GetYaxis()->SetRangeUser(0.8,1.2);
      hAllSBFitResult[m]->Draw("");
    }

    //I draw delta-phi projection done summing over pT assoc bins*****
    for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
      if (!IntisMC)      canvasSummedPt->cd(1+m);
      else       canvasSummedPt->cd(1+m+nummolt+1);
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
	  if (TMath::Abs(PtTrigMax-2.5) < 0.001 && type==8 && v>4 ) continue;
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

	canvasWings[0]->cd(m+1);
	hWings[m][IntisMC][PtTrig]->SetLineColor(kBlack);
	hWings[m][IntisMC][PtTrig]->SetMarkerColor(kBlack);
	hWings[m][IntisMC][PtTrig]->SetMarkerStyle(33);
	hWings[m][IntisMC][PtTrig]->GetXaxis()->SetRangeUser(0.01,7);
	hWings[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0.9,1.1);
	hWingsLeft[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0.9,1.1);
	hWingsLeft[m][IntisMC][PtTrig]->SetLineColor(kRed);
	hWingsLeft[m][IntisMC][PtTrig]->SetMarkerColor(kRed);
	hWingsLeft[m][IntisMC][PtTrig]->SetMarkerStyle(33);
	hWingsLeft[m][IntisMC][PtTrig]->GetXaxis()->SetRangeUser(0.01,7);
	hWingsRight[m][IntisMC][PtTrig]->SetLineColor(kBlue);
	hWingsRight[m][IntisMC][PtTrig]->SetMarkerColor(kBlue);
	hWingsRight[m][IntisMC][PtTrig]->SetMarkerStyle(33);
	hWingsRight[m][IntisMC][PtTrig]->GetXaxis()->SetRangeUser(0.01,7);
	hWingsRight[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0.9,1.1);
	hWings[m][IntisMC][PtTrig]->Draw("e");
	hWingsLeft[m][IntisMC][PtTrig]->Draw("same e");
	hWingsRight[m][IntisMC][PtTrig]->Draw("same e");
	tlineAtOnePt->Draw("same");
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

      //Ratio between ME and SE entries 
      MEEntries[m]=0;
      SEEntries[m]=0;
      for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
	if (TMath::Abs(PtTrigMax-2.5) < 0.001 && type==8 && v>4 ) continue;
	MEEntries[m] +=    hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetEntries();
	SEEntries[m] +=    hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetEntries();
	hMEEntriesvsPt[m]->SetBinContent(v+1, hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetEntries());
	hSEEntriesvsPt[m]->SetBinContent(v+1, hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetEntries());
	cout << " v " << v << endl;
	cout << "ME entries " << hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetEntries() << endl;
	cout << "SE entries " << hDeltaEtaDeltaPhi_RawSEbins[m][v][PtTrig]->GetEntries() << endl;
      }
    }
    hMEEntriesvsMult->SetBinContent(m+1, MEEntries[m]);
    hMEEntriesvsMult->SetBinError(m+1, sqrt(MEEntries[m]));
    hSEEntriesvsMult->SetBinContent(m+1, SEEntries[m]);
    hSEEntriesvsMult->SetBinError(m+1, sqrt(SEEntries[m]));
    hMEEntriesvsPt[m]->SetMarkerColor(Colormult[m]);
    hMEEntriesvsPt[m]->SetMarkerStyle(33);
    hMEEntriesvsPt[m]->SetLineColor(Colormult[m]);
    hSEEntriesvsPt[m]->SetMarkerColor(Colormult[m]);
    hSEEntriesvsPt[m]->SetMarkerStyle(33);
    hSEEntriesvsPt[m]->SetLineColor(Colormult[m]);

    hMEEntriesvsPt[m]->Sumw2();
    hSEEntriesvsPt[m]->Sumw2();
    hMEtoSEEntriesvsPt[m]= (TH1F*)hMEEntriesvsPt[m]->Clone("hMEtoSEEntriesvsPt"+ Smolt[m]);
    hMEtoSEEntriesvsPt[m]->Divide(hSEEntriesvsPt[m]);
    
    canvasMEtoSEMult->cd();
    hMEtoSEEntriesvsPt[m]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMEtoSEEntriesvsPt[m]->GetYaxis()->SetTitle("ME entries / SE entries");
    hMEtoSEEntriesvsPt[m]->GetYaxis()->SetRangeUser(0,12);
    hMEtoSEEntriesvsPt[m]->SetTitle("");
    hMEtoSEEntriesvsPt[m]->Draw("same");
    if (m==0)    legendMult->Draw("");

    //**************************************************************
    //*** Write on file *************
    //**************************************************************
    for(Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){
      for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
	if (TMath::Abs(PtTrigMax-2.5) < 0.001 && type==8 && v>4 ) continue;
	fileout->WriteTObject(hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][0]);
      }
      //      fileout->WriteTObject(hWings[m][IntisMC][0]);
      //      fileout->WriteTObject(hWingsRight[m][IntisMC][0]);
      //      fileout->WriteTObject(hWingsLeft[m][IntisMC][0]);
      fileout->WriteTObject(canvasPlotSE[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotME[m][IntisMC]); 
      fileout->WriteTObject(canvasPlot[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotProj[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotProjRatioJet[m][IntisMC]);
      if (m==nummolt && IntisMC==LimInfMC) canvasPlotSE[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf("); 
      else canvasPlotSE[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      canvasPlotME[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      canvasPlot[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      canvasPlotProj[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      if (isSidebandsAnalysis){
	for (Int_t t=0; t<3; t++){
	  canvasPlotProjSB[m][IntisMC][t]->SaveAs(nomefileoutputpdf+".pdf");
	  canvasPlotProjRatioSB[m][IntisMC][t]->SaveAs(nomefileoutputpdf+".pdf");
	}
      }
      fileout->WriteTObject(canvasPlotSB[m][IntisMC]); 
      fileout->WriteTObject(canvasPlotMESB[m][IntisMC]); 
      canvasPlotSB[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      canvasPlotMESB[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      if (isEtaEffComp){
	fileout->WriteTObject(canvasPlotMEEtaProjComp[m][IntisMC]);
	fileout->WriteTObject(canvasPlotMEEtaProjCompRatio[m][IntisMC]);
	canvasPlotMEEtaProjComp[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
	canvasPlotMEEtaProjCompRatio[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      }
      if (isCompWithMEFromHybrid){
	fileout->WriteTObject(canvasPlotMEEtaProjCompHybrid[m][IntisMC]);
	fileout->WriteTObject(canvasPlotMEEtaProjCompHybridRatio[m][IntisMC]);
	canvasPlotMEEtaProjCompHybrid[m][IntisMC]->SaveAs(nomefileoutputpdf + "_"+ Smolt[m]+"_MECompHybrid.pdf");
	canvasPlotMEEtaProjCompHybridRatio[m][IntisMC]->SaveAs(nomefileoutputpdf + "_" + Smolt[m] +"_MECompHybridRatio.pdf");
	canvasPlotMEEtaProjCompHybrid[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
	canvasPlotMEEtaProjCompHybridRatio[m][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      }
      if (isCompWithMEFromXi){
	fileout->WriteTObject(canvasPlotMEEtaProjCompXi[m][IntisMC]);
	fileout->WriteTObject(canvasPlotMEEtaProjCompXiRatio[m][IntisMC]);
	canvasPlotMEEtaProjCompXi[m][IntisMC]->SaveAs(nomefileoutputpdf + "_"+ Smolt[m]+"_MECompXi.pdf");
	canvasPlotMEEtaProjCompXiRatio[m][IntisMC]->SaveAs(nomefileoutputpdf + "_" + Smolt[m] +"_MECompXiRatio.pdf");
	canvasPlotMEEtaProjCompXi[m][IntisMC]->SaveAs(nomefileoutputpdf + ".pdf");
	canvasPlotMEEtaProjCompXiRatio[m][IntisMC]->SaveAs(nomefileoutputpdf + ".pdf");
      }

      for (Int_t t=0; t<3; t++){
	if (isSidebandsAnalysis){
	  canvasPlotProjSB[m][IntisMC][t]->SaveAs(nomefileoutputpdf+"_SBcomp_"+Smolt[m]+Region[t]+".pdf");
	  canvasPlotProjRatioSB[m][IntisMC][t]->SaveAs(nomefileoutputpdf+"_SBratio_"+Smolt[m]+Region[t]+".pdf");
	  fileout->WriteTObject( canvasPlotProjSB[m][IntisMC][t]);
	  fileout->WriteTObject( canvasPlotProjRatioSB[m][IntisMC][t]);
	}
	fileout->WriteTObject(canvasPlotProjEta[m][t][IntisMC]); 
	fileout->WriteTObject(canvasPlotRawSEProjEta[m][t][IntisMC]); 
	canvasPlotProjEta[m][t][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
	canvasPlotRawSEProjEta[m][t][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
      }
      fileout->WriteTObject(canvasPlotRawSEProjEta[m][3][IntisMC]); 
      canvasPlotRawSEProjEta[m][3][IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
    }

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
      
      fileout->WriteTObject(RatioHistoWidthGaussian[m][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussian[m][0][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussian[m][1][PtTrig]);
      fileout->WriteTObject(RatioHistoWidthGaussianAS[m][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussianAS[m][0][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussianAS[m][1][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussianEta[m][0][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussianEta[m][1][PtTrig]);
      
    }
  } //end of mult loop


  hMEEntriesvsMult->Sumw2();
  hSEEntriesvsMult->Sumw2();
  hMEtoSEEntriesvsMult= (TH1F*)     hMEEntriesvsMult->Clone("hMEtoSEEntriesvsMult");
  hMEtoSEEntriesvsMult->Divide(    hSEEntriesvsMult);
  hSEEntriesvsMult->SetTitle("SE entries");
  hMEEntriesvsMult->SetTitle("ME entries");
  hMEtoSEEntriesvsMult->SetTitle("ME entries/SE entries");
  hSEEntriesvsMult->GetXaxis()->SetTitle("Multiplicity class");
  hMEEntriesvsMult->GetXaxis()->SetTitle("Multiplicity class");
  hMEtoSEEntriesvsMult->GetXaxis()->SetTitle("Multiplicity class");
  hMEtoSEEntriesvsMult->GetYaxis()->SetRangeUser(0, 12);
  canvasMEtoSE->cd(1);
  hMEEntriesvsMult->Draw("");
  canvasMEtoSE->cd(2);
  hSEEntriesvsMult->Draw("");
  canvasMEtoSE->cd(3);
  hMEtoSEEntriesvsMult->Draw("");

  for (Int_t IntisMC=LimInfMC; IntisMC<=LimSupMC; IntisMC++){ 
    if (isSidebands){
      fileout->WriteTObject(canvasSBJetFitResult);
      fileout->WriteTObject(canvasSBAllFitResult);
      fileout->WriteTObject(canvasSBBulkFitResult);
    }
    fileout->WriteTObject(canvasPlotSEEtaProj[IntisMC]); 
    fileout->WriteTObject(canvasPlotMEEtaProj[IntisMC]); 
    fileout->WriteTObject(canvasPlotMEPhiProj[IntisMC]); 
    fileout->WriteTObject(canvasPlotMEEtaProjMolt[IntisMC]); 
    canvasPlotSEEtaProj[IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
    canvasPlotMEEtaProj[IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
    canvasPlotMEPhiProj[IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
    canvasPlotMEEtaProjMolt[IntisMC]->SaveAs(nomefileoutputpdf+".pdf");
  }
  fileout->WriteTObject(canvasSummedPt); 
  fileout->WriteTObject(canvasWings[0]);
  fileout->WriteTObject(canvasWidthGaussian[0]);
  fileout->WriteTObject(canvasWidthGaussianAS[0]);
  fileout->WriteTObject(canvasWidthGaussian[1]);
  fileout->WriteTObject(canvasWidthGaussianAS[1]);
  fileout->WriteTObject(canvasWidthGaussianEta[0]);
  fileout->WriteTObject(canvasWidthGaussianEta[1]);
  fileout->WriteTObject(canvasMEtoSE);
  fileout->WriteTObject(canvasMEtoSEMult);
  fileout->Close();

  canvasMEtoSEMult->SaveAs(nomefileoutputpdf+".pdf");
  canvasMEtoSE->SaveAs(nomefileoutputpdf+".pdf");
  canvasWings[0]->SaveAs(nomefileoutputpdf+".pdf)");

  cout << "baseline fits " << endl;
  for(Int_t m=nummolt; m>=0; m--){
    if (isHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //      if (m==0) continue;
    cout << "\n\n" << endl;
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
      if (TMath::Abs(PtTrigMax-2.5) < 0.001 && type==8 && v>4 ) continue;
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	cout << " m " << m << " v " << v << " PtTrig " << PtTrig  << endl;
	cout << "baseline value" <<  Baseline[m][v][PtTrig]->GetParameter(0)<< " +- " <<  Baseline[m][v][PtTrig]->GetParError(0)<< " number of sigmas from zero " <<  TMath::Abs(Baseline[m][v][PtTrig]->GetParameter(0)/ Baseline[m][v][PtTrig]->GetParError(0)) << endl;
      }
    }
  }

  for(Int_t m=nummolt; m>=0; m--){
    if (isHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //      if (m==0) continue;
    cout << "\n\n" << endl;
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){  
      if (TMath::Abs(PtTrigMax-2.5) < 0.001 && type==8 && v>4 ) continue;
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

  if (isTrigEff){
    for(Int_t m=0; m<nummolt+1; m++){
      if (isHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      if (isMC)      cout << "m " << m << " NTrigger " << NTrigger[m][1] << endl;
      else       cout << "m " << m << " NTrigger " << NTrigger[m][0] << endl;
    }
  }

  cout << "\n a partire dai file: \n" << PathIn[0] << endl; 
  if (isEtaEffComp) cout <<PathInComp[0] << "\n" ;
  if (isCompWithMEFromHybrid) cout <<PathInMEHybrid << "\n" ;
  if (isCompWithMEFromXi) cout <<PathInMEXi << "\n" ;

  cout << "\n\n ho creato il file " << nomefileoutput << " (also pdf)" << endl;

  for(Int_t m=nummolt; m>=0; m--){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    cout << "n trigger in mult range (all triggers)    " << m << "  " <<  " is MC " << 0 << "  " << NTrigger[m][0] <<   endl;
  }
}

