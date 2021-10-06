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

void AngularCorrelation_firstCasc(Int_t indexSysV0=0,  Int_t sysTrigger=0,  Int_t indexsysTrigger=0,Bool_t ishhCorr=0, Int_t sys=0, Bool_t SkipAssoc=1, Float_t ptjmin=3,  Float_t PtTrigMinFit=3, Int_t sysV0=0,Int_t type=0,  Int_t israp=0, Bool_t isMC=0, Bool_t isEfficiency=0, TString yearME=/*"161718_HM_hXi"/"161718Full_AOD234_hXi"/*"17pq_pp5TeV_hXi_pttrig0.15"/"17pq_hK0s"/*"17pq_pp5TeV_Hybrid"*/"1617_AOD234_hK0s"/*"17pq_hK0s"/*"LHC16kl_pass2_GP_Fio_Hybrid"/"1617GP_hK0s_Hybrid_New"/*"1617_AOD234_hK0s"/*"161718_hXi"/*"161718_MD_New_hXi_Hybrid"/*"2018f1_extra_hK0s_Fio"/"17pq_hK0s"/*"LHC17_AOD234_Red"/"AllhK0sHM_RedNo16k"/*"1617GP_hK0s"/*"2018f1_extra_15runs"/*"2018f1_extra_MylabelBis_15runs_hK0s_Hybrid"/*"2018f1_extra_Mylabel_15runs_hK0s_Hybrid"/*"2016k_HM_hK0s"/"1617_hK0s"/*"161718_MD_hXi_Hybrid_MCTruth"/*"2018f1g4_extra_hXi_Hybrid_MCTruth"/"2018f1_extra_15runs_NohDaughtersofK0s_hK0s_Hybrid"/*"2018g4_extra_EtaEff_hK0s_MCEff"*/, TString year= /*"161718_HM_hXi" /"161718Full_AOD234_hXi"/*"17pq_pp5TeV_hXi_pttrig0.15"/"17pq_hXi"/*"17pq_pp5TeV_Hybrid"*/"1617_AOD234_hK0s"/*"17pq_hK0s"/*"LHC16kl_pass2_GP_Fio_Hybrid"/"1617GP_hK0s_Hybrid_New"/*"1617_AOD234_hK0s"/*"161718_hXi"/*"161718_MD_New_hXi_Hybrid"/*"2018f1_extra_hK0s_Fio"/"17pq_hK0s"/*"LHC17_AOD234_Red"/"AllhK0sHM_RedNo16k"/*"1617GP_hK0s"/*"2018f1_extra_15runs"/*"2018f1_extra_MylabelBis_15runs_hK0s_Hybrid"//"2016k_HM_hK0s"/*"2018f1_extra_15runs_NohDaughtersofK0s_hK0s_Hybrid"/*"2018f1_extra_5runs_label_Hybrid_hK0s"/"1617_hK0s"/*"2018g4_extra_EtaEff_hK0s"/*"2018g4_extra_EtaEff_Hybrid_hK0s"/"161718_MD_EtaEff_hXi"/*"LHC16_17_GP_Hybrid_hXi"/*"2018f1g4_extra_hXi"/"2018g4_extra_hXi_SelTrigger"/*_15runsBis"/"1617_hK0s"/*AllMC_hXi"/"2018f1_extra_hK0s"/*"2016k_hK0s"/*"2016k_MECorr"/"Run2DataRed_MECorr_hXi"/*/, TString yearMC=""/*"161718Full_AOD235_hXi"/*"1617GP_hK0s_Hybrid_New" /"17pq_pp5TeV_Hybrid"/*"1617_GP_AOD235_With18c12b"/*"17pq_hK0s"/*"LHC16kl_pass2_GP_Fio_Hybrid"/*"1617GP_hK0s"/*"1617_GP_AOD235"/*"1617_GP_AOD235_With18c12b"/*"161718_hXi"/*"17pq_hK0s"/*"161718_MD_New_hXi_Hybrid"/*"2018f1_extra_hK0s_Fio"/*""//*"1617GP_hK0s"/*"2018f1_extra_15runs"/*"2018f1_extra_MylabelBis_15runs_hK0s_Hybrid"/"2019h11c_extra_HM_hK0s"/*"2018f1_extra_15runs_NohDaughtersofK0s_hK0s_Hybrid"/*"2018f1_extra_5runs_label_Hybrid_hK0s"/"1617MC_hK0s"/"2018g4_extra_EtaEff_hK0s"/*"161718_MD_EtaEff_hXi"/"2018f1g4_extra_hXi"/*"2018g4_extra_EtaEff_Hybrid_hK0s"/*"161718_MD_hXi_Hybrid"/*"LHC16_17_GP_Hybrid_hXi"/*"17d20bEPOS_hK0s"/"2018g4_extra_hXi_SelTrigger"/*"1617MC_hK0s"/"AllMC_hXi"/*"1617GP_hXi"/*"2016kl_hXi"/*"2018f1_extra_MECorr"/"2018f1_extra_Hybrid_hK0s"/*"17d20bEPOS_hXi"*/,  TString Path1 =""/*"_PtTrigMax2.5"/*"_NewMultClassBis"*/, TString Dir ="FinalOutput/",  Float_t ptjmax=15, Bool_t isBkgParab=1, Bool_t isMeanFixedPDG=1, Int_t MultBinning=0, Int_t PtBinning=0, Bool_t isSysDef=0, Bool_t isDefaultSel=0, Bool_t IsPtTrigMultDep=0, Bool_t isLoosest=0, Bool_t isTightest=0, Bool_t IsParticleTrue=0, Bool_t IsEfficiencyMassSel=0, Bool_t isSidebands=0,  Bool_t isMEFromHybrid=0, Bool_t isMEFromCorrectCentrality = 0, Bool_t isEtaEff=1, Bool_t isMEFromK0s=0, Bool_t isNewInputPath=1, Bool_t isHM=0, Bool_t isAllDeltaEta=0, Bool_t isMEFromPeak=0, Bool_t isMEFrom13TeV=0){ 

  if (isMEFromK0s && type==0) {cout << "the option isMEFromK0s is meant to be used when Xi is being analyzed" << endl; return;}
  if (year=="17pq_hXi") MultBinning=3;
  if (year=="17pq_pp5TeV_hXi_pttrig0.15") MultBinning=3;

  Int_t rebin=2;
  Int_t rebinx=2;  
  Int_t rebiny=2;
  TString year0 = "2016";
  Bool_t isEta05=0;
  Bool_t MasterThesisAnalysis=0;
  Bool_t isEnlargedDeltaEtaPhi=0;
  Bool_t isNotSigmaTrigger=0;
  Bool_t IsSpecial=0;
  Int_t TriggerPDGChoice=0;
  if (year == "1617_hK0s" || yearMC=="1617MC_hK0s") IsSpecial=1;

  TString sTriggerPDGChoice[3] = {"", "_IsOnlypiKpemu", "_IsNotSigmaOnly"};
  if (!isMC) TriggerPDGChoice=0;

  //isEtaEff==1 : eta dependence of efficiency considered
  if (isMC && !isEfficiency) {
    isEtaEff=0; //no efficiency correction needed in this case                                                
    IsParticleTrue=0; //if I deal with MCTruth/MCHybrid, all associated particles are true by default        
  }

  if (isSysDef && sys!=0) return; 
  //masterthesis analysis è usato quando devo fare analisi presentata per la tesi
  if (ishhCorr && type!=0) {cout << " for hh correlation type should be set to zero in order to get the correct files " << endl; return;}

  if (isDefaultSel && (isLoosest || isTightest) ) return;
  if (isLoosest && isTightest) return;
  if (isDefaultSel && sysTrigger==1) {cout << "to run the default selections you should put sysTrigger==0) " << endl;  return; }

  const  Float_t PtTrigMin=ptjmin;
  //  gStyle->SetOptStat(0);

  //lista degli effetti  sistematici studiati in questa macro
  if (type>0 && sys==5 && ptjmin==3) {cout << "enlargening of DeltaEta region for bulk is posible only for K0s " << endl; return;}
  if (sys==3 || sys>6) return;
  if (sysV0>6)return;
  if (ishhCorr && sys<4 && sys!=0) return; //se faccio correlazione hh non posso studiare sistematico associato a scelta regione Sidebands e peak nella distribuzione di massa invariante
  if (sysV0>2 && ishhCorr) return;
  //  if (sysTrigger!=0) return;
  if (sysV0!=0 && sys!=0) return;

  Int_t sysang=0;
  if (sys==1) sysang=1;
  if (sys==2) sysang=2;

  const Int_t numtipo=10;
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  if(isEta05) Srap[0] = "_Eta0.5";
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
 
  Float_t JetValue=0.84;
  Float_t JetValueDefault=0.84; //questo valore è indipendente dalla scelta del valore di sys ed è il valore per cui viene diviso lo yield nel jet, in modo da essere confrontatato con lo yield OJ e JOJ

  if (ishhCorr || (!ishhCorr && isEnlargedDeltaEtaPhi)) {//rendo ampiezza Jet uguale per hh e hK0s (l'output di hK0s è nel file DeltaEtaPhiEnlarged)
    JetValueDefault =1.;
    if ((sys==0 || sys==5)) JetValue=1.;
    else  if (sys==4) JetValue=0.9;
  }
  else {
    if (sys==4) JetValue=0.9;
    if (sys==6) JetValue=0.74;
  }


  Float_t BulkLowValue=0.84;
  Float_t BulkUpValue=1.1;
  Float_t InclusiveUpValue=1.1;
  if (isMC && !isEfficiency && isAllDeltaEta){
    InclusiveUpValue = 10;
  }
  Float_t InclusiveLowValue= -InclusiveUpValue;

  //InclusiveLowValue=0; // for asymmetric DeltaEta interval
  //  InclusiveUpValue=0-0.001; // for asymmetric DeltaEta interval

  if (sys==4 && !ishhCorr){
    BulkLowValue=0.9;
    BulkUpValue=1.1;
    // BulkUpValue=1.25;
  }
  if (sys==5 && !ishhCorr){
    BulkUpValue=1.25;
  }
  if (sys==6 && !ishhCorr){
    BulkLowValue=0.74;
  }

  if ((sys==0 || sys==4) && ishhCorr){
    BulkLowValue=1.05;
    BulkUpValue=1.3;
  }
  if (sys==5 && ishhCorr){
    BulkLowValue=1.15;
    BulkUpValue=1.4;
  }

  Dir+="DATA"+year0;
  TString file = year;
  TString file2 = year;
  if (isMC && isEfficiency) file += "_MCEff";
  //  file +=Path1; 
  if (ishhCorr) {
    file2+="_hhCorr";
    //    file+="_hhCorr";
    if (isMC)     file2+="_MCEff";
  }
  if (PtBinning>0)  file2+=Form("_PtBinning%i",PtBinning);
  file2+=Path1;
  //  if(isMC && isEfficiency) file = year + "_MCEff";
  if(isMC && !isEfficiency) file = yearMC + "_MCTruth";
  TString fileMC = yearMC;
  Int_t MC=0;
  if (isMC) MC=1;

  TString PathIn;
  TString PathME;
  TString PathInME;
  PathIn= Dir+"/histo/AngularCorrelation" + file;
  PathME= Dir+"/histo/AngularCorrelation" + yearME;
  if (PtBinning>0)  PathME+=Form("_PtBinning%i",PtBinning);
  //  if (isMC && !isEfficiency)PathIn+="_MCTruth";
  if (PtBinning>0)    PathIn+=Form("_PtBinning%i",PtBinning);
  PathIn+=Path1;
  if(type>=0){
    PathIn +="_";
    if (!ishhCorr)	  PathIn+=tipo[type];
    PathIn +=Srap[israp];
    PathIn +=SSkipAssoc[SkipAssoc];
    if (ishhCorr)      PathIn +="_hhCorr";
    PathME +="_";
    if (!ishhCorr){
      if (isMEFromK0s) PathME+=tipo[0];
      else PathME+=tipo[type];
    }
    PathME +=Srap[israp];
    if (isMEFromHybrid)    PathME +=SSkipAssoc[0];
    else     PathME +=SSkipAssoc[SkipAssoc];
    if (ishhCorr)      PathME +="_hhCorr";

  }
  if (!isMEFromHybrid && !isMEFromK0s) {
    if (isMEFrom13TeV) PathInME=PathME;
    else  PathInME=PathIn;
    PathInME += Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", 0 /*sysTrigger*/, sysV0, sysang, PtTrigMin);
    if (IsPtTrigMultDep) PathInME+="_IsPtTrigMultDep";
    if (IsParticleTrue) PathInME+= "_IsParticleTrue"; 
    PathInME+= sTriggerPDGChoice[TriggerPDGChoice];
    if (isNotSigmaTrigger) PathInME+= "_IsNotSigmaTrigger";
  }
  else if (isMEFromK0s){
    PathInME = PathME +  Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", 0 /*sysTrigger*/, sysV0, sysang, 3.);
    if (IsParticleTrue) PathInME+= "_IsParticleTrue"; 
    //PathInME+= "_isEtaEff";
  }
  else  PathInME = PathME +  Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", 0 /*sysTrigger*/, sysV0, sysang, PtTrigMin);
  //  PathInME += "_thinptbins";
  //  PathInME+= "_DCAz0.5";
  //  PathInME+= "_isPrimaryTrigger";

  if (isEtaEff) PathInME += "_isEtaEff";
  //PathInME += "_isEtaEff";
  if (MultBinning!=0) PathInME += Form("_MultBinning%i", MultBinning);
  //  PathInME += "_AllEff2019";
  //  PathInME+="_EPOS";
  PathInME += ".root";

  if(isSysDef && isDefaultSel)    PathIn+= Form("_SysT%i_SysV0Default_Sys%i_PtMin%.1f", sysTrigger,sysang, PtTrigMin);
  else  if(isSysDef && isLoosest)    PathIn+= Form("_SysT%i_SysV0Loosest_Sys%i_PtMin%.1f", sysTrigger,sysang, PtTrigMin);
  else  if(isSysDef && isTightest)    PathIn+= Form("_SysT%i_SysV0Tightest_Sys%i_PtMin%.1f", sysTrigger,sysang, PtTrigMin);
  else if(isSysDef && !isDefaultSel &&sysTrigger==0)    PathIn+= Form("_SysT%i_SysV0index%i_Sys%i_PtMin%.1f", sysTrigger, indexSysV0, sysang, PtTrigMin);
  else if(isSysDef && !isDefaultSel &&sysTrigger==1)    PathIn+= Form("_SysTindex%i_SysV0%i_Sys%i_PtMin%.1f", indexsysTrigger, 0, sysang, PtTrigMin);
  else  {
    PathIn+= Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin);
    if ( IsPtTrigMultDep)  PathIn+="_IsPtTrigMultDep";
    if (IsParticleTrue) PathIn+= "_IsParticleTrue"; 
    PathIn+= sTriggerPDGChoice[TriggerPDGChoice];
    if (isNotSigmaTrigger) PathIn+= "_IsNotSigmaTrigger";
    //    PathIn+="_thinptbins";
    //    PathIn+= "_DCAz0.5";
    //    PathIn+= "_isPrimaryTrigger";
  }
  if (isEtaEff) PathIn += "_isEtaEff";
  //PathIn += "_isEtaEff";
  if (MultBinning!=0) PathIn += Form("_MultBinning%i", MultBinning);
  //    PathIn+= "_AllEff2019";
  //  PathIn+="_EPOS";
  PathIn+=".root";

  if (MasterThesisAnalysis){
    PathIn= Dir+"/histo/AngularCorrelation" + file + Form("_SysT%i_SysV0%i_Sys%i", sysTrigger, sysV0, sysang)+".root";
  }
 
  TString PathInBis =  "FinalOutput/AnalysisResults" + file  + ".root";
  if (ishhCorr) PathInBis="FinalOutput/AnalysisResults" + file2  + ".root"; 

  TString PathInBisPart1;
  TString PathInBisPart2;
  if (IsSpecial){
    if (year == "1617_hK0s"){
      PathInBisPart1 =  "FinalOutput/AnalysisResults2016kehjl_hK0s.root"; 
      PathInBisPart2 =  "FinalOutput/AnalysisResultsLHC17_hK0s.root"; 
    }
    else if (year == "1617MC_hK0s"){
      PathInBisPart1 =  "FinalOutput/AnalysisResults2018f1d8_extra_hK0s_MCEff.root";
      PathInBisPart2 =  "FinalOutput/AnalysisResultsLHC17anch17_hK0s_MCEff.root"; 
    }
  }

  TString PathInEfficiency = Dir+ "/Efficiency/Efficiency" +fileMC;
  if (PtBinning>0)    PathInEfficiency+=Form("_PtBinning%i",PtBinning);
  //   PathInEfficiency+=Form("_PtBinning%i",1);
  PathInEfficiency+=Path1;// change
  if(type>=0){
    PathInEfficiency +="_";
    if (!ishhCorr)             PathInEfficiency +=tipo[type];
    PathInEfficiency +=Srap[israp];
    PathInEfficiency +=SSkipAssoc[SkipAssoc];
    if (ishhCorr)     PathInEfficiency +="_hhCorr";
  }

  if (isSysDef && isDefaultSel)   PathInEfficiency+= Form("_SysT%i_SysV0Default_PtMin%.1f.root", sysTrigger, PtTrigMin);
  else   if (isSysDef && isLoosest)   PathInEfficiency+= Form("_SysT%i_SysV0Loosest_PtMin%.1f.root", sysTrigger, PtTrigMin);
  else   if (isSysDef && isTightest)   PathInEfficiency+= Form("_SysT%i_SysV0Tightest_PtMin%.1f.root", sysTrigger, PtTrigMin);
  else if(isSysDef && !isDefaultSel && sysTrigger==0)   PathInEfficiency+= Form("_SysT%i_SysV0index%i_PtMin%.1f.root", sysTrigger, indexSysV0, PtTrigMin);
  else if(isSysDef && !isDefaultSel && sysTrigger==1)   PathInEfficiency+= Form("_SysTindex%i_SysV0%i_PtMin%.1f.root", indexsysTrigger, 0, PtTrigMin);
  else {
    PathInEfficiency+= Form("_SysT%i_SysV0%i_PtMin%.1f", sysTrigger, sysV0, PtTrigMin);
    if (IsPtTrigMultDep)    PathInEfficiency+= "_IsPtTrigMultDep";
    if (IsEfficiencyMassSel)    PathInEfficiency+= "_isEfficiencyMassSel";
    if (MultBinning!=0) PathInEfficiency += Form("_MultBinning%i", MultBinning);
    PathInEfficiency+= ".root";
  }
  //  if (ishhCorr) PathInEfficiency = Dir+ "/Efficiency/Efficiency" +fileMC  +"_hhCorr" + Path1+ Form("_SysT%i_SysV0%i_PtMin%.1f", sysTrigger, sysV0, PtTrigMin)+".root";
  if (MasterThesisAnalysis) PathInEfficiency = Dir+ "/Efficiency/Efficiency" +fileMC  +Form("_MCEff_Efficiency_SysT%i_SysV0%i", sysTrigger, sysV0)+".root";
 
  TString PathInMass= Dir+"/invmass_distribution_thesis/invmass_distribution";
  if (isMC && isEfficiency) PathInMass+="_MCEff";
  TString Title;
  if (PtBinning>0)  PathInMass+= Form("_PtBinning%i",PtBinning);
  PathInMass+= Path1; //+"_" ;
  //  PathInMass+= "Eta05_" ; //change
  PathInMass+="_";
  PathInMass+= year;

  cout << "path in: " << PathInBis << endl;
  cout << "path in: " << PathInME << endl;
  cout << "path in angular correlation: " << PathIn << endl;
  cout << "path in angular correlation (from where I take ME): " << PathInME << endl;
  if (!(isMC && !isEfficiency)){
    cout << "path efficiency: " <<PathInEfficiency << endl;
    cout << "path in mass first part: " << PathInMass << endl;
  }
  TString PathInMassDef;

  TFile *filepurezza;
  TFile *filein = new TFile(PathIn);
  if (!filein) return;
  TFile *fileinME = new TFile(PathInME);
  if (!fileinME) return;
  TFile *fileinbis = new TFile(PathInBis);
  TFile *fileinbisPart1;
  TFile *fileinbisPart2;
  if (IsSpecial){
    fileinbisPart1=new TFile(PathInBisPart1,"");
    fileinbisPart2=new TFile(PathInBisPart2,"");
  }
  if (!fileinbis) return;
  TFile *fileinEfficiency = new TFile(PathInEfficiency);
  if (!isEtaEff){
    if (!fileinEfficiency) return;
  }
  else if  (!(isMC && !isEfficiency)){
    if (!fileinEfficiency) return;
  }

  TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Xi", "Lambda","Xi", "Omega", "Omega", "Xi", "Omega"};

  TString ContainerName = "";
  if (isNewInputPath) {
    if (type==0){
      if (year.Index("Hybrid")!=-1)  ContainerName = "_hK0s_Task_Hybrid"; 
      else if (isMC && !isEfficiency) {
	if (year.Index("Fio")!=-1) ContainerName = "_hK0s_Task_Truth"; 
	//	else ContainerName = "_hK0s_Task_MCTruth"; 
	else ContainerName = "_hK0s_Task_Truth"; 
      }
      else if (isMC && isEfficiency){
	ContainerName = "_hK0s_Task_RecoAndEfficiency"; 
      }
      //      else ContainerName = "_hK0s_Task_Default"; 
      else ContainerName = "_hK0s_Task_"; 
    }
    else {
      //      if (year.Index("New")!=-1)  	ContainerName="_h" + tipo[type] +"_Task_MCTruth"; 
      if (year.Index("New")!=-1)  	ContainerName="_h" + tipo[type] +"_Task_Hybrid"; 
      else if (year.Index("Hybrid")!=-1)  ContainerName = "_hXi_Task_Hybrid"; 
      //      else if (isMC && !isEfficiency) ContainerName = "_hXi_Task_MCTruth"; 
      else if (isMC && !isEfficiency) ContainerName = "_hXi_Task_Truth"; 
      //else ContainerName = "_hXi_Task"; 
      else  {
	ContainerName="_hXi_Task_";
	if (type==8) 	ContainerName="_hXi_Task_Default";
	if (isHM) ContainerName="_hXi_Task_Default";
	if (TMath::Abs(PtTrigMin - 0.15) < 0.001) 	ContainerName="_hXi_Task_LowPtTrig";
      }
    }
  }
  TString TaskName ="";
  if (isNewInputPath) {
    if (!(year.Index("Hybrid")==-1)){
      if (year.Index("New")!=-1) TaskName = "_PtTrigMin3.0_PtTrigMax15.0";     
      else TaskName = "_MCTruth_PtTrigMin3.0_PtTrigMax15.0"; 
    }
    else if (isMC) {
      //      if (year.Index("Fio")!=-1) TaskName = "_MCTruth_PtTrigMin0.2_PtTrigMax15.0"; 
      if (year.Index("New")!=-1) TaskName = "_MCTruth_PtTrigMin3.0_PtTrigMax15.0"; 
      else  TaskName = "_MCTruth_PtTrigMin3.0_PtTrigMax15.0"; 
    }
    //    else TaskName = "_PtTrigMin3.0_PtTrigMax30.0";
    else {
      if (year == "2016k_HM_hK0s") {
	TaskName = "_PtTrigMin3.0_PtTrigMax30.0";      
      }
      else if (year == "17pq_pp5TeV_hXi_pttrig0.15"){
	TaskName = "_PtTrigMin0.2_PtTrigMax2.5";      
      }
      else  TaskName = "_PtTrigMin3.0_PtTrigMax15.0";
    }
  }

  cout << "Taskname: " << TaskName << endl;
  cout << "Container name: " << ContainerName << endl;
  TDirectoryFile *dir;
  if (type==0)    dir= (TDirectoryFile*)fileinbis->Get("MyTask"+ TaskName);
  else dir= (TDirectoryFile*)fileinbis->Get("MyTask"+dirinputtype[type]+ TaskName);
  if (!dir) {cout << "directory not present " << endl; return;}
  TDirectoryFile *dirPart1;
  TDirectoryFile *dirPart2;
  TList *listPart1;
  TList *listPart2;
  if (IsSpecial){
    dirPart1 = (TDirectoryFile*)fileinbisPart1->Get("MyTask"+dirinputtype[type]);
    dirPart2 = (TDirectoryFile*)fileinbisPart2->Get("MyTask"+dirinputtype[type]);
    listPart1 = (TList*)dirPart1->Get("MyOutputContainer"+ ContainerName);
    listPart2 = (TList*)dirPart2->Get("MyOutputContainer"+ ContainerName);
  }

  TList *list = (TList*)dir->Get("MyOutputContainer" + ContainerName);
  TList *list2 = (TList*)dir->Get("MyOutputContainer3"+ ContainerName);
  if (!list) {cout << "list not there " << endl; return;}

  cout << "\n********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "REMEMBER TO RUN READTREEPLCHIARA_first.C and READTREEPLCHIARA_second.C FIRST! "  << endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9;
  const Int_t numPtTrigger=1;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
  Int_t PtV0Min=1; //el 1
  if (!ishhCorr && type==0) PtV0Min=0;

  Int_t Marker[nummolt+1]={7,4,20,22,29,28};
  Int_t Color[nummolt+1]={2,3,4,6,7, 9};
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};

  Float_t ScaleFactor[nummolt+1][numzeta][numPtV0][numPtTrigger]={0};
  Float_t ScaleFactorBulk[nummolt+1][numzeta][numPtV0][numPtTrigger]={0}; 
  Float_t ScaleFactorJet[nummolt+1][numzeta][numPtV0][numPtTrigger]={0}; 
  Float_t ScaleFactorJetNotDef[nummolt+1][numzeta][numPtV0][numPtTrigger]={0}; 
  Float_t ScaleFactorJetBulk[nummolt+1][numzeta][numPtV0][numPtTrigger]={0}; 

  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 

  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  Double_t Nmoltpp5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Smolt[m] = Smolt0[m];
      Nmolt[m] = Nmolt0[m];
    }
    if (MultBinning==1){
      Smolt[m] = Smolt1[m];
      Nmolt[m] = Nmolt1[m];
    }
    if (MultBinning==2){
      Smolt[m] = Smolt2[m];
      Nmolt[m] = Nmolt2[m];
    }
    else if (MultBinning==3){
      Smolt[m] = Smoltpp5TeV[m];
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
    if (MultBinning==1){
      Nmolt[1] = 0;
      Nmolt[2] = 0;
      Smolt[0] = "0-0a";
      Smolt[1] = "0-0b";
      Smolt[2] = "0-0.01";
    }
  }

  Int_t numMultBins=100;
  Float_t UpperLimitMult = 100;
  if (isHM) UpperLimitMult = 0.1;

  TString Ssideband[numSB]={"", "_SB"};
  TString Ssidebandtitle[numSB]={"", " (sidebands)"};
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  //  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  //  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};


  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  //  TString SPtV0[numPtV0]={"", "","0.5-1", "1-1.5","1.5-2", "2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0.1-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", ""};
  //TString SNPtV0[numPtV0+1]={"0.0","0.0", "0.5","1.0","1.5","2.0","3.0","4.0","8.0"};
  SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }

  TString SPtV01[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  TString SNPtV01[numPtV0]={"0.1", "0.5", "0.8", "1.2","1.6", "2.0","2.5", "3.0", "4.0"};
  Double_t NPtV01[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};
  Int_t numPtV0Max=numPtV0;
  if (PtBinning==1) numPtV0Max = numPtV0;
  else numPtV0Max = numPtV0-1;
  if (PtBinning==1){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)   SPtV0[v] = SPtV01[v];
      if (v<numPtV0Max)   SNPtV0[v] = SNPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }

  TString SPtTrigger[numPtTrigger]={"3-30"};
  Double_t NPtTrigger[numPtTrigger+1]={ptjmin,ptjmax};
  TString Seta[numetabis]={"_eta<0.5", "_allEta"};
  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Int_t numbin=2; //questo numero deve essere pari

  Float_t binwx;
  Float_t binwy;
  Bool_t Bool;
 
  TH1F *fHistEtaLimitsOfRegion = new TH1F("fHistEtaLimitsOfRegion", "fHistEtaLimitsOfRegion",6, 0, 6);

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



  /*
    TH2D *fHistTriggervsMultPtRange =(TH2D*) fHistTriggervsMult->Clone("fHistTriggervsMultPtRange");

    fHistTriggervsMultPtRange->GetXaxis()->SetRange(fHistTriggervsMultPtRange->GetXaxis()->FindBin(PtTrigMultDep[m]+0.0001),fHistTriggervsMultPtRange->GetXaxis()->FindBin(PtTrigMax-0.0001) );
  */

  TH1D *fHistTriggervsMult_MultProj;

  fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMin+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(ptjmax-0.00001) );

  TH1D *HistoTriggerEfficiency;
  TH1D *HistContTriggerMolt;

  if (!isEtaEff){
    if (!(isMC && !isEfficiency)){
      HistoTriggerEfficiency     = (TH1D*)fileinEfficiency->Get("HistoTriggerEfficiency"); //trigger efficiency vs mult
      if (!HistoTriggerEfficiency) {cout << " no histotriggerefficiency" << endl; return;}
      HistContTriggerMolt        = (TH1D*)fileinEfficiency->Get("HistContTriggerMolt"); //trigger contamination factor vs mult
      if (!HistContTriggerMolt) {cout << " no histoConttriggermolt" << endl; return;}
    }
  }

  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_SEbins_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins_rapMolt[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MErap[nummolt+1][numzeta][numPtV0][numPtTrigger];
  // Float_t ssbFactor[nummolt+1][numPtV0];    
  TH1F * histoSSB[nummolt+1];
  TH1F * histo_Bcentral[nummolt+1];
  TH1F * histo_Bside[nummolt+1];
  TH1F*  fHistPtvsMult_mult[nummolt+1];
  TH1D* HistBI[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D* HistBII[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];

  TString RegionType[3]={"Jet", "BI", "All"};

  TH2D *hDeltaEtaDeltaPhi_ME_normbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_ACbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_ACbins_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB][numeta]; //SE/ME norm proiettato in delta Phi
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB][numeta]; //relative errors of histo projected above
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0SubFakeSB[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; 
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; 
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; //relative errors of histos above
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[nummolt+1][numzeta][numPtV0][numPtTrigger]; 
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[nummolt+1][numzeta][numPtV0][numPtTrigger]; 
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata  + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[nummolt+1][numzeta][numPtV0][numPtTrigger]; //relative error of histo above
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata  + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EtaAll[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis]; //istogramma finale su cui effettuo fit
  TH1D*	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_TrCorr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D*	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_TrCorr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D*	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_TrCorr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D*	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_TrCorr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D*	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled_TrCorr[nummolt+1][numzeta][numPtV0][numPtTrigger];


  Float_t norm_MEbins[nummolt+1][numPtV0][numSB]={0}; //valore ME isto in (0,0)
  Float_t norm_MEbins_norm[nummolt+1][numPtV0][numSB]={0}; //valore ME isto normalizzato in (0,0)
  Int_t   NTrigger[nummolt+1]={0}; //total number of trigger particles 
  Int_t   NTriggerV0[nummolt+1]={0}; //total number of trigger particles only in events with V > 0
  Float_t NSpectrum[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles
  Float_t NSpectrumError[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles
  Float_t NSpectrumV0[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles in events with V> 0
  TH1D*  fHistYvsMult[numetabis];

  TH1D*  fHistV0EfficiencyPtBins[nummolt+1]; //efficienza selezione V0 in PtV0 bins (vedi SPtV0)
  TH1D*  HistContV0PtBins[nummolt+1]; //contamination factor V0 in PtV0 bins (vedi SPtV0)

  TF1* 	  gauss[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis];
  TF1* 	  gaussint[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis];
  TF1*    pol0ZYAM[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TF1*    pol0Bulk[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TF1*    pol0BulkBis[nummolt+1][numzeta][numPtV0][numPtTrigger];

  Bool_t isProjectionPhi[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString Proj[3]={Form("|#Delta#eta| < %.1f, ",JetValue), Form("%.1f < |#Delta#eta| < %.1f, ",BulkLowValue, BulkUpValue),Form("|#Delta#eta| < %.1f, ",InclusiveUpValue)};

  TString TitleString[nummolt+1][numPtV0][2];
  TString TitleStringBis[nummolt+1][numPtV0];
  TString TitleStringTris[numPtV0];


  //********************************************************************* 
  //**************calcolo numero particelle di trigger*******************
  //********************************************************************* 
  for(Int_t m=0; m<nummolt+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //    if (m==0) continue;
    if (!isEtaEff){
      if (!(isMC && !isEfficiency)){
	fHistV0EfficiencyPtBins[m]= (TH1D*)fileinEfficiency->Get("fHistV0EfficiencyPtBins_" + Smolt[m]);
	if (!fHistV0EfficiencyPtBins[m]) {cout << "no V0EfficiencyPtBins in file " << PathInEfficiency<<endl; return;}; 
	HistContV0PtBins[m]= (TH1D*)fileinEfficiency->Get("HistContV0PtBins_" + Smolt[m]);
	if (!HistContV0PtBins[m]) {cout << "no HistContV0PtBins " << PathInEfficiency<<endl; return;}; 
      }
    }
    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.0001); j<=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.0001); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<=fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    cout << "\nn trigger in mult range (all triggers) [the one which should be used]    " << m << "  " <<  NTrigger[m] <<   endl;
  }
    
  for(Int_t sb=0; sb< numSB; sb++){
    if(sb==1 && isMC && !isEfficiency) continue;
    if(sb==1 && ishhCorr) continue;
    //    for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t m=nummolt; m>=0; m--){
      if (isHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      //if (m==0) continue;
      if (!ishhCorr){
	if(!isMC || (isMC && isEfficiency)){

	  PathInMassDef=PathInMass+ "_"+tipo[type]+Srap[israp]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG] + BkgType[isBkgParab];
	  if (isSysDef && isDefaultSel)	  PathInMassDef  += Form("_molt%i_sysT%i_sysV0Default_Sys%i_PtMin%.1f.root", m, sysTrigger,sysang,PtTrigMinFit);
	  else 	  if (isSysDef && isLoosest)	  PathInMassDef  += Form("_molt%i_sysT%i_sysV0Loosest_Sys%i_PtMin%.1f.root", m, sysTrigger,sysang,PtTrigMinFit);
	  else 	  if (isSysDef && isTightest)	  PathInMassDef  += Form("_molt%i_sysT%i_sysV0Tightest_Sys%i_PtMin%.1f.root", m, sysTrigger,sysang,PtTrigMinFit);
	  else if(isSysDef && !isDefaultSel && sysTrigger==0)	  PathInMassDef  += Form("_molt%i_sysT%i_sysV0index%i_Sys%i_PtMin%.1f.root", m, sysTrigger, indexSysV0,sysang,PtTrigMinFit);
	  else if(isSysDef && !isDefaultSel && sysTrigger==1)	  PathInMassDef  += Form("_molt%i_sysTindex%i_sysV0%i_Sys%i_PtMin%.1f.root", m, indexsysTrigger, 0,sysang,PtTrigMinFit);
	  else	{
	    PathInMassDef  += Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f", m, sysTrigger, sysV0,sysang,PtTrigMinFit);
	    if (MultBinning!=0) PathInMassDef += Form("_MultBinning%i", MultBinning);
	    PathInMassDef += ".root";
	  }

	  //	  if(type==0)	  PathInMassDef=PathInMass+"_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, sysang, PtTrigMinFit);
	  if (MasterThesisAnalysis) 	  PathInMassDef=PathInMass+"_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i_Sys%i.root", m, sysTrigger, sysV0, sysang);
	  cout << "******************************" << endl;
	  cout << "******************************" << endl;
	  cout << "\npath in mass def completo " << PathInMassDef << endl;
	  if (!IsParticleTrue){
	    filepurezza= new TFile(PathInMassDef);
	    if (!filepurezza) {cout << " il file inv mass non c'è " << endl; return;}
	    histoSSB[m]=(TH1F*)filepurezza->Get("histo_SSB");
	    histo_Bcentral[m]=(TH1F*)filepurezza->Get("histo_Bcentral");
	    histo_Bside[m]=(TH1F*)filepurezza->Get("histo_Bside");
	  }
	}
      }
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	    for (Int_t sb=0; sb<2; sb++){
	      if(sb==1 && isMC && !isEfficiency) continue;
	      TitleString[m][v][sb]= "in mult. " + Smolt[m] +"% and " + SNPtV0[v] + " GeV/c < p_{T}^{Assoc} < " + SNPtV0[v+1]+" GeV/c"+ Ssidebandtitle[sb];
	    }
	    TitleStringBis[m][v]= "in mult. " + Smolt[m] +"% and " + SNPtV0[v] + " GeV/c < p_{T}^{Assoc} < " + SNPtV0[v+1]+" GeV/c";
	    TitleStringTris[v]= "in " + SNPtV0[v] + " GeV/c < p_{T}^{Assoc} < " + SNPtV0[v+1]+" GeV/c";

	    cout << "\n\n *********************************************" << endl;
	    cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v]<< " per sideband(1) or not(0) = " << sb << endl;
	  
	    nameME[m][z][v][tr][sb]="ME_";
	    nameME[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];
	    nameSE[m][z][v][tr][sb]="SE_";
	    nameSE[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];

	    //	    if (isEtaEff && sb==0 ){
	    if (isEtaEff && (sb==0 || isSidebands)){
	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameSE[m][z][v][tr][sb]+ "_Effw");
	    }
	    else{
	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameSE[m][z][v][tr][sb] +"");
	    }
	    if (!hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]){cout << "missing histo: " << nameSE[m][z][v][tr][sb] << endl;  return;}

	    if (isEtaEff && (sb==0 || isSidebands) && !isMEFromHybrid){
	      hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]= (TH2D*)fileinME-> Get(nameME[m][z][v][tr][sb]+ "_Effw");
	    }
	    else{
	      hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]= (TH2D*)fileinME-> Get(nameME[m][z][v][tr][sb]+"");
	    }
	    if (!hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]) {cout << "missing histo: " << nameME[m][z][v][tr][sb] << endl;  return;}

	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Rebin2D(rebinx,rebiny);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Rebin2D(rebinx,rebiny);

	    hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]= (TH2D*)	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Clone(nameSE[m][z][v][tr][sb]+"_RelError");
	    hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]= (TH2D*)	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Clone(nameME[m][z][v][tr][sb]+"_RelError");

	    for (Int_t j=1; j<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsX(); j++){
	      for (Int_t l=1; l<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsY(); l++){
		Int_t Binuniv=  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBin(j,l);
		hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]->SetBinContent(Binuniv,         hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinError(Binuniv)/hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinContent(Binuniv));
		hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]->SetBinContent(Binuniv,         hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinError(Binuniv)/  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(Binuniv));
		hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]->SetBinError(Binuniv,0);
		hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]->SetBinError(Binuniv,0);
	      }
	    }

	    /////////// primo modo per normalizzare ME 
	    binwx= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->GetBinWidth(1);
	    binwy= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->GetBinWidth(1);

	    cout << "\nafter the rebin: bin x (delta phi) " << binwx <<  " bin y (delta eta) " << binwy << endl;
	    Int_t cont=0;
	    for(Int_t i=0; i < numbin ; i++ ){
	      for(Int_t j=0; j < numbin ; j++ ){
		norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001+(i-numbin/2)*binwx), hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001+(j-numbin/2)*binwy));
		cont++;
	      }
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/cont;  	    
  
	    cout << "\nvalore histo in DeltaEta=0 e DeltaPhi=0 (quadrato) " <<  norm_MEbins[m][v][sb] << endl;

	    /////////// secondo modo per normalizzare ME 
	    norm_MEbins[m][v][sb]=0;
	    cont=0;
	    for(Int_t i=0; i <hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetNbinsY()  ; i++ ){
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001),i+1);
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-0.000001),i+1);
	      cont++;
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/2/cont;
	    
	    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 (media in delta phi, 2 righe) " <<  norm_MEbins[m][v][sb] << endl;

	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_norm");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_AC");
	    if (norm_MEbins[m][v][sb]==0){
	      cout << "sb " << sb << endl;
	      cout << "normalizzazione ME non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;
	      //norm continue;  
	    }

	    //	    cout << "\n divido SE distribution per ME disatribution normalized" << endl;
	    if (norm_MEbins[m][v][sb]!=0) hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]-> Scale(1./norm_MEbins[m][v][sb]);
	    //norm	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]);
	    //****************provo a utilizzare pair acceptance della classe di molteplicità 0-100% per tutte le molt************
	    //	    cout << "I divide " << endl;
	    if (!ishhCorr){
	      if (isMEFromCorrectCentrality){
		if (isMEFromPeak && sb==1) 		hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][0]);
		else hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]);
	      }
	      else {
		if (isMEFromPeak && sb==1) 		hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]);
		else hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][sb]);
	      }
	    }
	    if (ishhCorr)	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]);
	    //	    cout << "I have divided " << endl; //error are propagated assuming non correlated histograms
	    hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_AC_RelError");

	    for (Int_t j=1; j<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsX(); j++){
	      for (Int_t l=1; l<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsY(); l++){
		Int_t Binuniv=  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBin(j,l);
		//	      cout <<	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinContent(Binuniv) << " +-  " << 	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinError(Binuniv) << endl;
		//	      cout <<	      hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]->GetBinContent(Binuniv) << " +-  " << 	      hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]->GetBinError(Binuniv) << endl;
		//	      cout <<	      hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinContent(Binuniv) << " +-  " << 	      hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinError(Binuniv) << endl;
		hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]->SetBinContent(Binuniv,         hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinError(Binuniv)/  hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinContent(Binuniv));
		hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]->SetBinError(Binuniv,0);
	      }
	    }


	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->SetTitle("Raw SE distribution "+ TitleString[m][v][sb]);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->SetTitle("ME distribution "+ TitleString[m][v][sb]);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->SetTitle("Normalized ME distribution "+ TitleString[m][v][sb]);
    	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->SetTitle("SE corrected for pair acceptance "+TitleString[m][v][sb]);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    norm_MEbins_norm[m][v][sb]= hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001), hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001));
	    cout << "valore histo  in DeltaEta=0 e DeltaPhi=0, histo normalizzato  " <<  norm_MEbins_norm[m][v][sb] <<"\n"<< endl;

	    //histograms of relative errors of SE, ME and AC (not yet here)


	    //*****************************************************************************************************************
	    //proietto in Delta Phi
	    //****************************************************************************************************************

	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][0]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaJet", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-JetValue)+1, hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(JetValue)-1, "E"));
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][2]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaAll",  hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(InclusiveLowValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(InclusiveUpValue), "E"));

	    HistBI[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaBulkI", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-BulkUpValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-BulkLowValue), "E"));
	    HistBII[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaBulkII", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(BulkLowValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(BulkUpValue), "E"));
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]= (TH1D*)	    HistBI[m][z][v][tr][sb]->Clone(nameME[m][z][v][tr][sb]+"_AC_phi_etaBI");
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->Add(	    HistBII[m][z][v][tr][sb]);

	    /*//the following gives consistent result to the Add() method used just above; the errors given by root after Add(A+B) method is the sqrt(errA^2 + errB^2)
	      for(Int_t j=1; j < HistBI[m][z][v][tr][sb]->GetNbinsX(); j++){
	      cout << "error before " <<  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->GetBinError(j) << endl;
	      //	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->SetBinContent(j, HistBI[m][z][v][tr][sb]->GetBinContent(j)+ HistBII[m][z][v][tr][sb]->GetBinContent(j));
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->SetBinError(j, sqrt(pow(HistBI[m][z][v][tr][sb]->GetBinError(j),2)+ pow(HistBII[m][z][v][tr][sb]->GetBinError(j),2)));	    
	      cout << "error after " <<  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->GetBinError(j) << endl;
	      }
	    */

	    for(Int_t eta=0; eta<numeta; eta++){
	      //	      cout << " eta " << eta << endl;
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetTitle("#Delta#phi projection in " +Proj[eta] + TitleString[m][v][sb]);
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetXaxis()->SetTitle("#Delta#phi (radians)");
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetYaxis()->SetTitle("Counts");
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetXaxis()->SetTitleOffset(1);
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetYaxis()->SetTitleOffset(1);
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetYaxis()->SetTitleSize(0.05);
	      if (sb==1) hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetLineColor(kRed);
	      Bool_t NotScale=kFALSE;
	      if (!ishhCorr){
		if (!IsParticleTrue){
		  if (sb==1){
		    if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) NotScale=kTRUE;
		    if(!NotScale){
		      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->Scale((Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001))/histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001)));
		    }
		    else {
		      cout << " I could not scale the SB distribution due to empty sidebands " << endl;
		    }
		    //		  cout << "Bcentral " << (Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001)) << " Bside " << histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001)) << " B central/Bside " << (Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001))/histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))<< endl;
		  }
		}
	      }
	      hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]= (TH1D*)	 	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->Clone(nameME[m][z][v][tr][sb]+Form("_AC_phi_eta%i_RelError",eta));
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetTitle("Relative errors of #Delta#phi projection in " +Proj[eta] + TitleString[m][v][sb]);
	      for(Int_t j=1; j < hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->GetNbinsX(); j++){
		hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->SetBinContent(j,TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetBinError(j)/hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetBinContent(j)) );
		hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->SetBinError(j,0);
	      }
	      hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->GetYaxis()->SetRangeUser(0,1);
	      hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->SetLineColor(kBlack);
	    }//chiusura ciclo eta
	  } //chiusura ciclo su pt della V0
	}//chiusura ciclo su Pt Trigger
      }//chiusura ciclo su z vertice
    }//chiusura ciclo su molteplcità
  }//chiusura ciclo su sideband or not

  cout << "********************************************" << endl;
  cout << "**                                       **" << endl;
  cout << "  **                                   **" << endl;
  cout << "    **                               **" << endl;
  cout << "      **                           **" << endl;
  cout << "        **                       **" << endl;
  cout << "          **                    **" << endl;
  cout << "            **                **" << endl;
  cout << "              **            **" << endl;
  cout << " \n\n ********** second part of the program *******************" << endl;

  for(Int_t m=0; m<nummolt+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //if (m==0) continue;
    // for(Int_t m=nummolt; m>=0; m--){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	  if (!ishhCorr){
	    cout << "\n\n *********************************************" << endl;
	    cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v] << endl;
	    if(!isMC || (isMC && isEfficiency)){
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][0]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_rap");
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][1], hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][0]);
	      Title="ME sidebands/ME peak "+TitleStringBis[m][v];
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->SetTitle(Title);
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->GetXaxis()->SetTitle("#Delta#eta");
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->GetXaxis()->SetTitleOffset(1.3);
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->GetYaxis()->SetTitle("#Delta#phi (radians)");

	    }
	  }

	  if((!isMC || (isMC && isEfficiency)) && !ishhCorr) Bool= (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0);
	  if((!isMC || (isMC && isEfficiency)) && ishhCorr) 	   Bool= (norm_MEbins[m][v][0]==0);
	  else  Bool= (norm_MEbins[m][v][0]==0);
	  if (Bool){
	    cout << "normalizzazione ME (SB o peak region) non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;
	    //norm	    continue;  
	  }
	  
	  //***************************************************************
	  //sottraggo fondo data da "finte" K0s se analizzo hV0
	  //***************************************************************
	  for(Int_t eta=0; eta< numeta; eta++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta] =  (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub", eta));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta] =  (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub_RelError", eta));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->SetName(nameME[m][z][v][tr][0]+Form("_eta%i_AC_phi_V0Sub",eta));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->SetName(nameME[m][z][v][tr][0]+Form("_eta%i_AC_phi_V0Sub_RelError",eta));

	    if (isMC && !isEfficiency) continue;
	    if (ishhCorr) continue;

	    hDeltaEtaDeltaPhi_ACbins_phi_V0SubFakeSB[m][z][v][tr][eta]=(TH1D*) 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_eta"+ RegionType[eta] +"_FakeSB");
	    if (!IsParticleTrue){
	      hDeltaEtaDeltaPhi_ACbins_phi_V0SubFakeSB[m][z][v][tr][eta]->Scale(1-histoSSB[m]->GetBinContent(v+1));
	    }

	    //error check: both methods below give the same results 
	    /*
	      for(Int_t i=1; i<hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetNbinsX(); i++ ){
	      cout << " bin " << i << ": " << 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinContent(i)<< " +- " <<     hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinError(i)<<endl;
	      cout << "bin in old way " << 	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->GetBinContent(i) -  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->GetBinContent(i) << " +- "<<  sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->GetBinError(i),2) +pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->GetBinError(i),2)) << endl;
	      }
	    */
	    cout << "\n\n\n ***+" << endl;

	    if (!IsParticleTrue)	  {
	      cout << "purity for " << m << " v " << m << " " << histoSSB[m]->GetBinContent(v+1)<< endl;
	      if (!isSidebands)  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->Scale(histoSSB[m]->GetBinContent(v+1));//this for a fast fake removal
	      else hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->Add(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta],-1); 
	    }
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetXaxis()->SetTitle("#Delta#phi (radians)");

	    for (Int_t j=1; j<= hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetNbinsX() ; j++){
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->SetBinContent(j, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinError(j)/  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinContent(j));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->SetBinError(j,0);
	    }
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->GetYaxis()->SetRangeUser(0,1);

	  }

	  cout << " 	  //sottraggo distribuzione del bulk " << endl;
	  //***************************************************************
	  //sottraggo distribuzione del bulk 
	  //***************************************************************
	  ScaleFactor[m][z][v][tr]=(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValue)))/2./ (hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)));
	  //	  cout << (hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValue))) << "   " <<(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)))<<"  " << (hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue))) << "  " << hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue));
	  cout << "\n" <<ScaleFactor[m][z][v][tr] << endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(ScaleFactor[m][z][v][tr]);

	  //**********here I rebin the two histograms before subtraction*************** 
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub_Rebin", 1));
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub_Rebin", 0));
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->Rebin(rebin);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Rebin(rebin);
	  //***************************************************************************

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->SetLineColor(kGreen+2);

	  //norebin	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSubFit_EffCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EffCorr");

	  //no rebin 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Add(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1],-1);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Add(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1],-1);

	  //fit of the bulk in the delta phi jet region to subtract bulk
	  pol0Bulk[m][z][v][tr]= new TF1("pol0",Form("pol0Bulk_m%i_v%i",m,v), -1, 1);
	  pol0Bulk[m][z][v][tr]->SetLineColor(kGreen+2);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->Fit(  pol0Bulk[m][z][v][tr], "R+");
	  for (Int_t b=1; b<= hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetNbinsX(); b++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetBinContent(b, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetBinContent(b)-  pol0Bulk[m][z][v][tr]->GetParameter(0));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetBinError(b, sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetBinError(b),2)+ pow(pol0Bulk[m][z][v][tr]->GetParError(0),2)));
	  }
	  //	  error check //!!!!!ROOT computes errors assuming uncorrelated histograms 
	  /* 
	     for(Int_t i=1; i< 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetNbinsX(); i++){
	     cout << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinContent(i) -  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinContent(i)<< endl;
	     cout << "norm factor value (should be zero) " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetNormFactor()<< endl;
	     //	    cout << "error squared unchecked " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinErrorSqUnchecked(i) << " error squared " <<  pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2)<< endl;
	     cout << "error of the first histogram of the sum " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinContent(i)<< " +- " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i)<< endl;
	     cout << "error of the second histogram of the sum " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinContent(i)<< " +- " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i)<< endl;
	     cout << "errors calculated by root: " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< " +- " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinError(i)<< " rel: " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinError(i)/hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	     cout << "errors assuming zero correlation: " << sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2) +  pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2)) << " rel error = " << sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2) +  pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2))/ 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	     cout << "errors assuming total correlation (ro =1): " << TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i)- hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i)) << " rel error = " << TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i)- hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i)) / hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	     //the following formula is not correct	    cout << "errors assuming 'partial' correlation  " <<  sqrt(TMath::Abs(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2)- pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2))) << " rel error = " <<  sqrt(TMath::Abs(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2)- pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2))) / 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	     cout << "sqrt(bin content) of BulkSub histo " << sqrt(TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)))<< " rel error " <<1./ sqrt(TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)))<<endl;
	     }
	  */

	  for(Int_t i=1; i< 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetNbinsX(); i++){
	    //	    cout << " I set errors " << endl;
	    //	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->SetBinError(i, TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetBinError(i)- hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->GetBinError(i)));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->SetBinError(i, sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetBinError(i),2)+pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->GetBinError(i),2)));
	  }
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr] = (TH1D*) 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorrNotScaled");

	  //riscalo per deltaEta width tutte e tre le proiezioni in delta phi (jet, bulk, inclusive) e clono histo. L'histo clonato sarà diviso per efficienza

	  //*****bulk*********
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(1./ScaleFactor[m][z][v][tr]);
	  ScaleFactorBulk[m][z][v][tr]=(1./2/(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)))); 
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(ScaleFactorBulk[m][z][v][tr]);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");	 
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetLineColor(kGreen+2);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetTitle("UE distribution in " +Proj[1] + TitleStringBis[m][v]);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr_RelError");	 

	  //********jet*****
	  ScaleFactorJet[m][z][v][tr]=(1./(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValueDefault)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValueDefault))));
	  cout << "	  ScaleFactorJet[m][z][v][tr] " << 	  ScaleFactorJet[m][z][v][tr]<< endl;
	  ScaleFactorJetNotDef[m][z][v][tr]=(1./(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValue)))); //this is needed fr siplaying purposes (otherwise jetnot bulksub does not correctly superimpose to bulk distribution)---no side effects are expected
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->Scale(ScaleFactorJet[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Scale(ScaleFactorJet[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Scale(ScaleFactorJetNotDef[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Scale(ScaleFactorJet[m][z][v][tr]);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]= (TH1D*)	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_RelError");
	  for (Int_t j=1; j<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->GetNbinsX() ; j++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->SetBinContent(j,hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinError(j)/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(j));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->SetBinError(j,0);
	  }
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetLineColor(kBlue);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetLineColor(kGreen+3);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetLineColor(kRed);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->SetLineColor(kRed);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetTitle("Jet distribution in " +Proj[0] + TitleStringBis[m][v]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetTitle("Jet distribution in " +Proj[0] + TitleStringBis[m][v] + "Not scaled by #Delta#eta");

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorr_RelError");

	  //********inclusive*******
	  ScaleFactorJetBulk[m][z][v][tr]=1./(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(InclusiveUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(InclusiveLowValue)));
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Scale(ScaleFactorJetBulk[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulkEffCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulkEffCorr_RelError");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->SetTitle("Jet + UE distribution in " +Proj[0] + TitleStringBis[m][v]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Scale(1./ScaleFactorJetBulk[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->SetTitle("Jet + UE distribution in " +Proj[0] + TitleStringBis[m][v]);


	  //****************************************************************************************************************
	  //divido angular correlation proietatta in deltaphi per efficienza di selezione Trigger e V0 e per contamination factors e per SSB
	  //****************************************************************************************************************

	  if (isMC && !isEfficiency) 	      isProjectionPhi[m][z][v][tr]=kTRUE;
	  else	  if (!isMC || (isMC && isEfficiency)){
	    if (!isEtaEff){
	      cout << PathInEfficiency << endl;
	      cout << "\n \n Trigger selection efficiency: " <<HistoTriggerEfficiency->GetBinContent(m+1) << endl;
	      cout << "\n V0 selection efficiency:  " << fHistV0EfficiencyPtBins[m]->GetBinContent(v+1) << endl;
	      cout << "\n Trigger contamination factor:  " << (HistContTriggerMolt->GetBinContent(m+1)) << endl;
	      cout << "\n V0 contamination factor:  " << (HistContV0PtBins[m]->GetBinContent(v+1))<< endl;

	      if (IsParticleTrue) {
		HistContV0PtBins[m]->SetBinContent(v+1, 0);
	      }
	      if((fHistV0EfficiencyPtBins[m]->GetBinContent(v+1))!=0 ){

		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));

		//setto errore all'istogramma scalato 
		Float_t Error_BulkSub_EffCorr=0;
		Float_t Error_BulkSub_EffCorrNotScaled=0;
		Float_t Error_BulkSubFit_EffCorr=0;
		Float_t Error_EffCorr=0;
		Float_t Error_Bulk_EffCorr=0;
		Float_t Error_JetBulk_EffCorr=0;
		Float_t Error_JetBulk_EffCorrNotScaled=0;

		for(Int_t i=1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetNbinsX(); i++ ){
		  Error_BulkSub_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

		  Error_BulkSubFit_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

		  Error_BulkSub_EffCorrNotScaled= 	      Error_BulkSub_EffCorr*ScaleFactorJet[m][z][v][tr];

		  Error_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

		  Error_Bulk_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

		  Error_JetBulk_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

		  Error_JetBulk_EffCorrNotScaled= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));


		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetBinError(i, Error_BulkSub_EffCorr);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetBinError(i, Error_BulkSubFit_EffCorr);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetBinError(i, Error_BulkSub_EffCorrNotScaled);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->SetBinError(i, Error_EffCorr);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetBinError(i, Error_Bulk_EffCorr);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->SetBinError(i, Error_JetBulk_EffCorr);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->SetBinError(i, Error_JetBulk_EffCorrNotScaled);

		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->SetBinContent(i, Error_BulkSub_EffCorr/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinContent(i));
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->SetBinContent(i, Error_Bulk_EffCorr/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinContent(i));
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->SetBinContent(i, Error_JetBulk_EffCorr/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinContent(i));
		  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);

		}
	    
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1);
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]   ->GetYaxis()->SetRangeUser(0,1);
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1);
		isProjectionPhi[m][z][v][tr]=kTRUE;
	      }
	      else {
		cout << "\n ****** efficienza pari a zero********* " << endl;
		isProjectionPhi[m][z][v][tr]=kFALSE;
		continue;
	      }

	    }
	    else {
	      for(Int_t i=1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetNbinsX(); i++ ){
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->SetBinContent(i,hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinError(i) / hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinContent(i));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->SetBinContent(i, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinError(i)/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinContent(i));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->SetBinContent(i, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinError(i)/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinContent(i));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);
	      }
	      isProjectionPhi[m][z][v][tr]=kTRUE;
	    }
	  } // efficiency correction

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_TrCorr[m][z][v][tr]=(TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorr_TrCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_TrCorr[m][z][v][tr]        =(TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]        ->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EffCorr_TrCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_TrCorr[m][z][v][tr]   =(TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]   ->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_Bulk_EffCorr_TrCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_TrCorr[m][z][v][tr]=(TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulk_EffCorr_TrCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled_TrCorr[m][z][v][tr]=(TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulk_EffCorrNotScaled_TrCorr");
	  //	  cout << "!!!!!!!!! bef scaling writing m " << m << " v " << v << endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_TrCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_TrCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_TrCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_TrCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled_TrCorr[m][z][v][tr]->Scale(1./NTrigger[m]);

	  

	  cout << " ZYAM " << endl;
	  
	  //fit for ZYAM*******
	  pol0ZYAM[m][z][v][tr]= new TF1(Form("pol0ZYAM_m%i_v%i",m,v), "pol0",1.1,1.9);
	  pol0ZYAM[m][z][v][tr]->SetLineColor(kRed);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Fit(pol0ZYAM[m][z][v][tr], "R+");

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]= (TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EffCorr_BulkSubZYAM");
	  for(Int_t i=1; i<=      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetNbinsX(); i++){
	    //	    cout << " I set errors " << endl;
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->SetBinContent(i, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)-pol0ZYAM[m][z][v][tr]->GetParameter(0));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->SetBinError(i,sqrt( pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinError(i),2)+pow(pol0ZYAM[m][z][v][tr]->GetParError(0),2)));
	
	  }

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->Rebin(rebin);
	  //*******************

	  //fit for BulkSubFitBis*******

	  pol0BulkBis[m][z][v][tr]= new TF1("pol0",Form("pol0BulkBis_m%i_v%i",m,v), -1,1);
	  pol0BulkBis[m][z][v][tr]->SetLineColor(kGreen+3);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->Fit(pol0BulkBis[m][z][v][tr], "R+");

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]= (TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EffCorr_BulkSubBulkFit");
	  for(Int_t i=1; i<=	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetNbinsX(); i++){
	    //	    cout << " I set errors " << endl;
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetBinContent(i, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)-pol0BulkBis[m][z][v][tr]->GetParameter(0));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetBinError(i,sqrt( pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinError(i),2)+pow(pol0BulkBis[m][z][v][tr]->GetParError(0),2)));
	
	  }

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->Rebin(rebin);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetLineColor(kGreen+3);
	  //*******************

	  /* put it above
	   */

	  //la seguente parte era utilizzata per sottrarre il contributo dell'UE nella regione di jet, utilizzando il valore medio di alcuni bin simmetricamente disposti intorno a zero. Si è scelto di cambiare metodo. Per riutilizzarlo, alcuni nomi di histo vanno cambiati. 
	  //****************************************************************************************************************
	  //bkg subtraction
	  //****************************************************************************************************************
	  /*
	    for(Int_t etabis=0; etabis<2; etabis++){
	    Float_t DeltaPhiBkg=0;
	    Float_t counter=0;
	    Float_t DeltaPhiBkgError=0;

	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->FindBin(-1.1); i++){
	    counter++;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i);
	    }
	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.84); i++){
	    counter++;
	    //	    cout <<" content " <<  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i) << endl;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i);
	    }
	    DeltaPhiBkg=DeltaPhiBkg/counter;

	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-1.1); i++){
	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i)-DeltaPhiBkg,2);
	    }
	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.84); i++){//simmetrico attorno a pi/mezzi

	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i)-DeltaPhiBkg,2);
	    }
	    DeltaPhiBkgError=sqrt(DeltaPhiBkgError/counter/(counter+1));


	    cout << "numero di bin su cui si media per sottrarre fondo: " << counter << "\n valore del fondo di coppie non correlate: " << DeltaPhiBkg << " +- " << DeltaPhiBkgError<<  endl;
	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]= (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_BkgSub"+Seta[etabis]);
	    //hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->Sumw2();

	    if(etabis==1){
	    for(Int_t i =1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->SetBinContent(i, (Float_t)(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i))-DeltaPhiBkg);
	    //the errors below are exactly the same
	    //	    cout << "\nbin n: " << i << endl;
	    //	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinError(i)<< endl;
	    //	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->GetBinError(i)<< endl;
	    }//
	    }


	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->Rebin(rebin);
	    } //etabis
	  */
	}//v
      }//tr
    }//z
  }//m


  cout << "\n\n\n Studio la dipendenza della pair correction dalla classe di molteplicità"<< endl;
  for(Int_t sb=0; sb< numSB; sb++){
    if(sb==1 && isMC && !isEfficiency && !ishhCorr) continue;
    if(sb==1 && ishhCorr) continue;
    for(Int_t m=0; m<nummolt+1; m++){
      if (isHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      //if (m==0) continue;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=PtV0Min; v<numPtV0Max; v++){

	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Clone(nameME[m][z][v][tr][sb] + "_rapMolt");
	    if (!	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]) return;
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->SetName(nameME[m][z][v][tr][sb] + "_rapMolt");
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->SetTitle("ME " +Smolt[m]+"% / ME 0-100%"+TitleStringTris[v] );
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->Divide( hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_MEbins[5][z][v][tr][sb]);
	  }
	}
      }
    }
  }

 
  auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->SetHeader("Selezioni applicate");     


  TCanvas *canvasJetFit[nummolt+1];
  for(Int_t m=0; m< nummolt+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //if (m==0) continue;
    canvasJetFit[m]= new TCanvas (Form("canvasJet_%i",m),Form("canvasJet_%i",m) , 1200, 1000);
    canvasJetFit[m]->Divide(4,2);
    for(Int_t v=PtV0Min; v< numPtV0Max; v++){
      //      cout << "v" << v << endl;
      canvasJetFit[m]->cd(v+1);

      if((!isMC || (isMC && isEfficiency)) && !ishhCorr) Bool= (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0);
      if((!isMC || (isMC && isEfficiency)) && ishhCorr) Bool= (norm_MEbins[m][v][0]==0);
      else Bool= (norm_MEbins[m][v][0]==0);
    
      if (Bool) {
	cout << " normalization not performed " << endl;	   
	//norm 	continue;
      }
      if(!isProjectionPhi[m][0][v][0]){
	cout << "projection not performed " << endl;
	continue;
      }	   
      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][0][v][0]->DrawCopy("");
    
    }
  }

  cout << "\n sto per scrivere su file " << endl;
  // salvo solo istogrammi con norm != 0 sia per sideband che per regione centrale, con efficienze diverse da zero e con integrale siiideband non =0 (ossia solo histo con significato fisico) 

  TString PathOut1;
  //  if (isMC && !isEfficiency) PathOut1=Dir+"/histo/AngularCorrelation" + file  +"_Output.root";
  //  else {
  PathOut1 = Dir+"/histo/AngularCorrelation"  +file+Path1 ;
  if (PtBinning>0)    PathOut1 +=Form("_PtBinning%i",PtBinning);
  //    if (sys==0) PathOut1 +=Form("_Jet%.2f", hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(JetValue))) ;
  if(type>=0){
    if (!ishhCorr)      PathOut1 +="_"+tipo[type];
    PathOut1 +=Srap[israp];
    PathOut1 +=SSkipAssoc[SkipAssoc];
  }
  if(isBkgParab)     PathOut1+=  "_isBkgParab";
  if (isSysDef && isDefaultSel)    PathOut1+=   Form("_SysT%i_SysV0Default_Sys%i_PtMin%.1f_Output.root", sysTrigger,  sys, PtTrigMin);
  else     if (isSysDef && isLoosest)    PathOut1+=   Form("_SysT%i_SysV0Loosest_Sys%i_PtMin%.1f_Output.root", sysTrigger,  sys, PtTrigMin);
  else     if (isSysDef && isTightest)    PathOut1+=   Form("_SysT%i_SysV0Tightest_Sys%i_PtMin%.1f_Output.root", sysTrigger,  sys, PtTrigMin);
  else if(isSysDef && !isDefaultSel &&sysTrigger==0)     PathOut1+=   Form("_SysT%i_SysV0index%i_Sys%i_PtMin%.1f_Output.root", sysTrigger, indexSysV0, sys, PtTrigMin);
  else if(isSysDef && !isDefaultSel &&sysTrigger==1)     PathOut1+=   Form("_SysTindex%i_SysV0%i_Sys%i_PtMin%.1f_Output.root", indexsysTrigger, 0, sys, PtTrigMin);
  else {
    PathOut1+=   Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f_Output", sysTrigger, sysV0, sys, PtTrigMin);
    if (yearMC == "17d20bEPOS_hK0s")    PathOut1+= "_EPOS";
    if (yearMC == "17d20bEPOS_hK0s_EtaEff")    PathOut1+= "_EPOS";
    if (yearMC == "17d20bEPOS_hXi")    PathOut1+= "_EPOS";
    if (IsPtTrigMultDep)    PathOut1+= "_IsPtTrigMultDep";
    if (IsParticleTrue) PathOut1+= "_IsParticleTrue";
    if (IsEfficiencyMassSel)    PathOut1+= "_IsEfficiencyMassSel";
    if (isSidebands)    PathOut1+= "_Sidebands";
    if (isNotSigmaTrigger) PathOut1+= "_IsNotSigmaTrigger";
    if (isMEFromHybrid) PathOut1+= "_IsMEFromHybrid";
    if (isMEFromK0s) PathOut1+= "_IsMEFromK0s";
    if (isMEFrom13TeV) PathOut1+= "_IsMEFrom13TeV";
    if (isMEFromPeak) PathOut1+= "_IsMEFromPeak";
    if (isMEFromCorrectCentrality) PathOut1+= "_IsMEFromCorrectCentrality";
    if (isEtaEff) PathOut1+= "_IsEtaEff";
    if (isAllDeltaEta) PathOut1+= "_isAllDeltaEta";
    PathOut1+= sTriggerPDGChoice[TriggerPDGChoice];
    if (MultBinning!=0) PathOut1 += Form("_MultBinning%i", MultBinning);
    //    PathOut1+="_TryisEta05";
    //    PathOut1+="_thinptbinsbis";
    //    PathOut1+= "_DCAz0.5";
    //    PathOut1+= "_isPrimaryTrigger";
    //    PathOut1 += "_Eff18g4extra";
    //    PathOut1 += "_AllEff2019";
    //    PathOut1+="_EPOS";
    PathOut1+=".root";
  }
  if (isEnlargedDeltaEtaPhi)    PathOut1 = Dir+"/histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sys, PtTrigMin)+"_DeltaEtaPhiEnlarged_Output.root";
  //    PathOut1 = Dir+"/histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sys, PtTrigMin)+"_Output.root";
  if (MasterThesisAnalysis)    PathOut1 = Dir+"/histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i_Sys%i", sysTrigger, sysV0, sys)+"_Output.root";
  //  }
  
  TFile *fileout = new TFile(PathOut1, "RECREATE");
  for(Int_t m=0; m<nummolt+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //if (m==0) continue;
    //    cout << "\n" <<m << endl;
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	  //	  cout << "v " << v << endl;
	  if(!isMC || (isMC && isEfficiency)){
	    if (!ishhCorr){
	      if (!IsParticleTrue){
		if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) {
		  cout << " B side =0 " << endl;
		  continue;
		}
	      }
	    }
	    
	    if (!ishhCorr)	     hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Write();
	  }
	  
	  if((!isMC || (isMC && isEfficiency)) && !ishhCorr) Bool= (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0);
	  if((!isMC || (isMC && isEfficiency)) && ishhCorr) Bool= (norm_MEbins[m][v][0]==0);

	  if (Bool) {
	    cout << " pt interval" << v << endl;
	    cout << " normalization not performed " << endl;	   
	    //norm 	    continue;
	   
	  }
	  if(!isProjectionPhi[m][z][v][tr]){
	    cout << "projection not performed " << endl;
	    continue;
	  }	   
	  for(Int_t sb=0; sb<numSB; sb++){
	    if(sb==1 && isMC && !isEfficiency && !ishhCorr) continue;
	    if(sb==1 && ishhCorr) continue;
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]->Write();
	    for (Int_t eta=0; eta<=2; eta++){
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->Write();
	      hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->Write();
	    }
	  }

	  //decommentare per avere proiezione in out-of-jet region sovrapposto a proiezione in jet-region
	  // hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(1./ScaleFactorBulk[m][z][v][tr]);
	  // hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(ScaleFactor[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->SetTitle("");

	  for (Int_t eta=0; eta<=2; eta++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->Write();
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->Write();
	    if (!(isMC && !isEfficiency)){
	      hDeltaEtaDeltaPhi_ACbins_phi_V0SubFakeSB[m][z][v][tr][eta]->Write();
	    }
	  }

	  //	  cout << "!!!!!!!!! bef writing m " << m << " v " << v << endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->SetTitle("Jet region projection corrected for fake "+tipoTitle[type] + ", efficiency and DeltaEta width");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->SetTitle("Jet region projection corrected for fake "+ tipoTitle[type]+", bulk subtracted and DeltaEta width");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->SetTitle("Raltive errors of Jet region projection corrected for fake "+ tipoTitle[type]+", bulk subtracted and DeltaEta width");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetTitle("Jet region projection bulk-sub, (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->Write();


	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetTitle("Jet region projection bulk-sub (bulk fit), (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->SetTitle("Jet region projection bulk-sub (ZYAM), (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetTitle("Jet region projection bulk-sub (bulk fit), (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetTitle("Jet region projection bulk-sub, (fake-"+tipoTitle[type]+" + Eff) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->Write();

	  //
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetTitle("OJ region projection, (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->Write();


	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->SetTitle("JOJ region projection, (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->SetTitle("JOJ region projection, (fake-"+tipoTitle[type]+" + Eff) corrected");


	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Write();


	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_TrCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_TrCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_TrCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_TrCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled_TrCorr[m][z][v][tr]->Write();
	  pol0ZYAM[m][z][v][tr]->Write();
	  
	}
      }
    }
  }

  cout << "Drawing canvas" << endl;
  TCanvas *canvasDraw[nummolt+1][numzeta][numPtTrigger][numPtV0][5];
  for(Int_t m=0; m<nummolt+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //if (m==0) continue;
    //    cout << "\n" <<m << endl;
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	  cout << "la normalizzazione del ME non è stato effettuata per: " << endl;
	  if (norm_MEbins[m][v][0]==0) cout << " m " << m << " PtV0 " << v << " sideband 0"<< endl;
	  if (norm_MEbins[m][v][1]==0) cout << " m " << m << " PtV0 " << v << " sideband 1"<< endl;
 
	  if (!ishhCorr){
	    if(!isMC || (isMC && isEfficiency)){
	      if (!IsParticleTrue){
		if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) {
		  cout << " B side =0 " << endl;
		  continue;
	      
		}
	      }
	    }
	    canvasDraw[m][z][v][tr][0]=new TCanvas (Form("CanvasDraw%i_%i_%i_%i_0", m,z,v,tr),Form("CanvasDraw%i_%i_%i_%i_0", m,z,v,tr), 800, 600 );
	    canvasDraw[m][z][v][tr][0]->cd();
	    if (ishhCorr) 	    hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Draw();
	      
	    //	    if (m==0 && v ==1) canvasDraw[m][z][v][tr][0]->SaveAs("MESBRatio05.pdf");
	    canvasDraw[m][z][v][tr][0]->Close();	      
	    //fileout->WriteTObject(canvasDraw[m][z][v][tr][0]);
	  }
	}

      }
    }
  }

  fHistEtaLimitsOfRegion->SetBinContent(1,hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(-JetValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(1, "JetLowValue");
  fHistEtaLimitsOfRegion->SetBinContent(2,hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(JetValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(2, "JetUpValue");
  fHistEtaLimitsOfRegion->SetBinContent(3,hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(BulkLowValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(3, "BulkLowValue");
  fHistEtaLimitsOfRegion->SetBinContent(4,hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(BulkUpValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(4, "BulkUpValue");
  fHistEtaLimitsOfRegion->SetBinContent(5,hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(InclusiveLowValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(5, "InclusiveLowValue");
  fHistEtaLimitsOfRegion->SetBinContent(6,hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(InclusiveUpValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(6, "InclusiveUpValue");

  fHistEtaLimitsOfRegion->Write();

  fileout->Close(); 
  
  if (!isMEFromCorrectCentrality)  cout << "****** sto utilizzando pair acceptance della classe molt 0-100%*****\n"<< endl;

  cout << "DeltaEta intervals used:" << endl;
  cout << "--> DeltaEta Jet interval " << hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(-JetValue)) << " - " << hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(JetValue)) << endl;

  cout << "--> DeltaEta Bulk interval " << hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(BulkLowValue)) << " - " << hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(BulkUpValue)) << endl;

  cout << "--> DeltaEta Full interval " << hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(InclusiveLowValue)) << " - " << hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[nummolt][0][1][0][0]->GetXaxis()->FindBin(InclusiveUpValue)) << endl;

  cout << "! se valori si sovrappongono, modificare intervallo DeltaEta \n" << endl;

  Int_t OJEntries=0;
  Int_t NBins=0;
  for (Int_t m=0; m<nummolt+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //cout << "\n m " << m << endl;
    for (Int_t v=PtV0Min; v<numPtV0Max; v++){
      //cout << " v" << v << endl;
      //cout << " scale factor bulk "<<       ScaleFactorBulk[m][0][v][0]<< endl;
      //cout << " scale factor jetbulk "<<       ScaleFactorJetBulk[m][0][v][0]<< endl;
      //cout << " scale factor jet "<<       ScaleFactorJet[m][0][v][0]<< endl;
      OJEntries=0;
      NBins=0;
      //      cout << SPtV0[v] << " " << hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetEntries() << endl;
      for (Int_t  bEta = hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetXaxis()->FindBin(BulkLowValue); bEta <= hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetXaxis()->FindBin(BulkUpValue); bEta++) {
	for (Int_t  bPhi = hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetYaxis()->FindBin(1); bPhi <= hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetYaxis()->FindBin(2); bPhi++) {
	  NBins++;
	  //	  cout << "bEta " << bEta << " bPhi " << bPhi << " " << hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetBin(bEta, bPhi) << endl;
	  OJEntries +=	hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetBinContent( hDeltaEtaDeltaPhi_SEbins[m][0][v][0][0]->GetBin(bEta, bPhi));
	}
      }
      //      cout << SPtV0[v] << " n. bins over which sum is performed " << NBins << "; bin sum: " << OJEntries<< endl;
    }
  }

  if (InclusiveLowValue!= -InclusiveUpValue) cout << "\n\n!!!!!!!! asymmetric inclusive interval is being used!!!!!\n\n"<<endl;

  cout << "\nRatio between ME and SE entries: " << endl;
  for(Int_t v=PtV0Min; v<numPtV0Max; v++){
    cout <<  "Pt: " << NPtV0[v] << " ME entries: " <<  hDeltaEtaDeltaPhi_MEbins[nummolt][0][v][0][0]->GetEntries()<< " SE entries: " <<  hDeltaEtaDeltaPhi_SEbins[nummolt][0][v][0][0]->GetEntries()<< " ratio ME/SE entries " << (float)hDeltaEtaDeltaPhi_MEbins[nummolt][0][v][0][0]->GetEntries()/ hDeltaEtaDeltaPhi_SEbins[nummolt][0][v][0][0]->GetEntries()  << endl;
  }
 
  for(Int_t m=0; m<nummolt+1; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    //cout << "m " << m << " NTrigger " << NTrigger[m] << endl;
  } 

  cout << "\nPartendo dai file \n" << PathIn << "\n"<< PathInME << " (for ME) \n" << PathInBis << endl;
  cout << "\nFor purity and SB " << PathInMassDef << endl;
  if (!isEtaEff)  cout << PathInEfficiency;
  cout << "\n\nho creato il file: "<<  PathOut1 << endl;
 
  if (ishhCorr) cout << " use sysV0 = 0,1,2 and sys==0 " << endl;

  cout << "Bin width of 2D angular correlation histograms " << endl;
  cout << "\nbin x (delta phi): " << binwx <<  " bin y (delta eta): " << binwy << endl;

  if (MasterThesisAnalysis) cout << " *****************************Be aware master thesis analysis has been selected!!" << endl;

}
